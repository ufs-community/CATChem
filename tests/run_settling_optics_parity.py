#!/usr/bin/env python3
"""Certify the C++ optics-table settling path against the legacy Fortran oracle.

This mirrors ``run_default_numerical_parity.py`` but isolates the
``simple_scheme: true`` (NetCDF optics-table) branch of the GOCART settling
process.  Both runners consume the same shared MET profile and the same
deterministic ``1e-12*(1+0.01*level+0.001*species)`` initial chemistry, so any
post-timestep divergence is attributable to the settling core, not to inputs.

The staged configuration rewrites the checked-in Default config in four ways:

  * schedule only the ``settling`` process in the ``test1`` phase;
  * flip ``gocart/simple_scheme`` to ``true`` so both cores read the Mie tables;
  * point ``mie/directory`` at the supplied optics directory;
  * replace the incomplete ``_3``/``_2`` optics filenames with the complete
    ``_5`` variants that carry the ``rhop``/``growth_factor`` variables the
    simple kernel requires (a missing ``rhop`` silently freezes a species).

Run it inside a container that mounts both source trees at identical absolute
paths (the legacy tree is a separate worktree), e.g.::

    docker run --rm \\
      -v /Users/barry/Documents/CATChem:/Users/barry/Documents/CATChem \\
      -v /Users/barry/Documents/catchem-legacy-ref:/Users/barry/Documents/catchem-legacy-ref \\
      -w /Users/barry/Documents/CATChem cece-dev:latest \\
      python3 tests/run_settling_optics_parity.py \\
        --candidate-runner build-parity/tests/default_parity_candidate_runner \\
        --legacy-runner /Users/barry/Documents/catchem-legacy-ref/build-legacy/tests/parity/legacy_parity_runner \\
        --workdir /tmp/parity_settling_optics \\
        --optics-dir /Users/barry/Documents/CATChem/specs/tmp
"""

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

# The complete optics tables (with rhop + growth_factor) use the trailing "_5"
# version suffix; the Default config ships the incomplete "_3"/"_2" variants.
COMPLETE_MIE_FILES = {
    "optics_SS.v3_3.nc": "optics_SS.v3_5.nc",
    "optics_DU.v15_3.nc": "optics_DU.v15_5.nc",
    "optics_BC.v1_3.nc": "optics_BC.v1_5.nc",
    "optics_OC.v1_3.nc": "optics_OC.v1_5.nc",
    "optics_SU.v1_3.nc": "optics_SU.v1_5.nc",
    "optics_NI.v2_5.nc": "optics_NI.v2_5.nc",  # already complete
}


def stage(run_dir: Path, source_root: Path, columns: int, levels: int, optics_dir: str) -> tuple[Path, Path]:
    run_dir.mkdir(parents=True, exist_ok=True)
    config_dir = source_root / "tests" / "Configs" / "Default"
    for name in ("CATChem_new_config.yml", "CATChem_species.yml"):
        shutil.copy2(config_dir / name, run_dir / name)

    config = run_dir / "CATChem_new_config.yml"
    text = config.read_text(encoding="utf-8")
    # No external emissions: settling tendencies must not depend on inventory I/O.
    text = text.replace("emission_filename: ./CATChem_emission.yml",
                        "emission_filename: ./CATChem_parity_zero_emissions.yml", 1)
    text = text.replace("diagnostics: true", "diagnostics: false")
    # Isolate the settling process in the test1 phase.
    text, replaced = re.subn(
        r"(  test1:\n    description: \"Test phase 1\"\n    processes:\n)(?:      - .+\n)+",
        r"\1      - settling\n",
        text,
        count=1,
    )
    if replaced != 1:
        raise RuntimeError("could not isolate the test1 process phase to settling")
    # Enable the optics-table (Mie) branch of the GOCART settling scheme.
    # The checked-in Default config now ships this baseline itself, so accept
    # either current value and force the simple path idempotently.
    text, replaced = re.subn(r"(simple_scheme: )\w+", r"\1true", text, count=1)
    if replaced != 1:
        raise RuntimeError("could not enable simple_scheme in the staged configuration")
    # Point mie/directory at the shared optics directory (trailing slash required
    # by both the legacy and C++ path-join rules).  Target the mie block's own
    # key by its checked-in value so an unrelated "directory:" earlier in the
    # file cannot be rewritten instead.
    directory = optics_dir if optics_dir.endswith("/") else optics_dir + "/"
    mie_directory_line = '  directory: "./ExtData/monochromatic/"'
    if mie_directory_line not in text:
        raise RuntimeError("could not locate the mie/directory line in the staged configuration")
    text = text.replace(mie_directory_line, '  directory: "' + directory + '"', 1)
    # Swap the incomplete optics filenames for the complete _5 variants.
    for old, new in COMPLETE_MIE_FILES.items():
        text = text.replace(old, new)
    config.write_text(text, encoding="utf-8")
    (run_dir / "CATChem_parity_zero_emissions.yml").write_text("categories: {}\n", encoding="utf-8")

    profile = run_dir / "settling_optics_met.json"
    subprocess.run(
        [sys.executable, str(source_root / "tests" / "build_parity_met_profile.py"),
         "--profile", str(source_root / "tests" / "MetProfiles" / "Profile_NCWCP.csv"),
         "--columns", str(columns), "--levels", str(levels), "--output", str(profile)],
        check=True,
    )
    return config, profile


def run(label: str, executable: Path, run_dir: Path, config: Path, profile: Path, steps: int, dt: float) -> Path:
    snapshot = run_dir / f"{label}.json"
    # The runners execute inside their staged run_dir, so the binary path must
    # be absolute; a relative --candidate-runner would not resolve there.
    subprocess.run(
        [str(executable.resolve()), "--config", str(config), "--met-profile", str(profile),
         "--snapshot", str(snapshot), "--steps", str(steps), "--dt", str(dt)],
        cwd=run_dir,
        check=True,
    )
    if not snapshot.is_file():
        raise RuntimeError(f"{label} runner did not write {snapshot}")
    return snapshot


def settling_species_divergence(legacy: Path, candidate: Path) -> dict:
    """Report per-species max abs/rel change from the initial state for both cores.

    Confirms the optics path is actually doing physics (dust/sea-salt move)
    rather than passing trivially on a frozen field.
    """

    def load(path: Path):
        doc = json.loads(path.read_text(encoding="utf-8"))
        field = doc["snapshots"][0]["fields"]["concentration"]
        cols, levels, nspec = field["shape"]
        flat = field["values"]
        # Canonical snapshot index: sp + nspec*(lev + levels*col).
        grid = {}
        for col in range(cols):
            for lev in range(levels):
                for sp in range(nspec):
                    grid[(sp, lev, col)] = flat[sp + nspec * (lev + levels * col)]
        return field["species"], cols, levels, nspec, grid

    lspecies, lc, ll, ln, lgrid = load(legacy)
    cspecies, cc, cl, cn, cgrid = load(candidate)
    if lspecies != cspecies or (lc, ll, ln) != (cc, cl, cn):
        raise RuntimeError("legacy and candidate snapshots describe different grids/species")
    aerosols = [name for name in lspecies if name.startswith(("dust", "seas", "bc", "oc", "so4", "msa"))]
    report = {}
    for name in aerosols:
        sp = lspecies.index(name)
        worst = 0.0
        for key in lgrid:
            if key[0] != sp:
                continue
            initial = 1.0e-12 * (1.0 + 0.01 * key[1] + 0.001 * sp)
            legacy_change = abs(lgrid[key] - initial)
            cand_change = abs(cgrid[key] - initial)
            worst = max(worst, abs(legacy_change - cand_change))
        report[name] = worst
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy-runner", type=Path, required=True)
    parser.add_argument("--candidate-runner", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--optics-dir", type=str, required=True,
                        help="absolute directory containing the complete _5 optics tables")
    parser.add_argument("--columns", type=int, default=2)
    parser.add_argument("--levels", type=int, default=20)
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument("--dt", type=float, default=3600.0)
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    parser.add_argument("--atol", type=float, default=1.0e-12)
    args = parser.parse_args()

    legacy_config, legacy_profile = stage(args.workdir / "legacy", args.source_root, args.columns, args.levels,
                                          args.optics_dir)
    candidate_config, candidate_profile = stage(args.workdir / "candidate", args.source_root, args.columns,
                                                args.levels, args.optics_dir)
    legacy_snapshot = run("legacy", args.legacy_runner, args.workdir / "legacy", legacy_config, legacy_profile,
                          args.steps, args.dt)
    candidate_snapshot = run("candidate", args.candidate_runner, args.workdir / "candidate", candidate_config,
                             candidate_profile, args.steps, args.dt)

    divergence = settling_species_divergence(legacy_snapshot, candidate_snapshot)
    print("settling optics parity: per-species |candidate-legacy| change magnitude (max over grid):")
    for name, value in sorted(divergence.items()):
        print(f"  {name:8s} {value:.3e}")

    status = subprocess.run(
        [sys.executable, str(args.source_root / "tests" / "compare_parity_snapshots.py"),
         str(legacy_snapshot), str(candidate_snapshot), "--rtol", str(args.rtol), "--atol", str(args.atol)],
        check=False,
    ).returncode

    report = args.workdir / "report.json"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        json.dumps({
            "schema_version": 1,
            "process": "settling",
            "scheme": "gocart/simple_scheme(optics)",
            "steps": args.steps,
            "dt": args.dt,
            "columns": args.columns,
            "levels": args.levels,
            "rtol": args.rtol,
            "atol": args.atol,
            "passed": status == 0,
            "per_species_max_abs_divergence": divergence,
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"parity report written: {report}")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
