#!/usr/bin/env python3
"""Run native legacy and C++ Default-config parity executables in isolation.

Each runner is a standalone executable built against exactly one core. It must
accept ``--config``, ``--met-profile``, ``--snapshot``, and ``--steps``. It
must write the normalized snapshot schema consumed by compare_parity_snapshots.
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path


def stage(run_dir: Path, source_root: Path, columns: int, levels: int, processes: tuple[str, ...]) -> tuple[Path, Path]:
    run_dir.mkdir(parents=True, exist_ok=True)
    config_dir = source_root / "tests" / "Configs" / "Default"
    for name in ("CATChem_new_config.yml", "CATChem_species.yml"):
        shutil.copy2(config_dir / name, run_dir / name)
    # The parity control case intentionally contains no externally supplied
    # emissions.  Sea-salt and dust remain enabled: their process-generated
    # tendencies are part of the physics comparison.  Pointing both runners
    # at a minimal inventory prevents external inventory parsing/file I/O from
    # becoming part of a core-numerics test.
    config = run_dir / "CATChem_new_config.yml"
    config_text = config.read_text(encoding="utf-8").replace(
            "emission_filename: ./CATChem_emission.yml",
            "emission_filename: ./CATChem_parity_zero_emissions.yml",
            1,
    )
    # Snapshot concentrations are the parity oracle.  Disable optional
    # diagnostics so unallocated legacy diagnostic buffers cannot affect an
    # otherwise independent process comparison.
    config_text = config_text.replace("diagnostics: true", "diagnostics: false")
    phase_processes = "".join(f"      - {name}\n" for name in processes)
    config_text, replacements = re.subn(
        r"(  test1:\n    description: \"Test phase 1\"\n    processes:\n)(?:      - .+\n)+",
        r"\1" + phase_processes,
        config_text,
        count=1,
    )
    if replacements != 1:
        raise RuntimeError("could not isolate the test1 process phase in the staged configuration")
    config.write_text(config_text, encoding="utf-8")
    (run_dir / "CATChem_parity_zero_emissions.yml").write_text("categories: {}\n", encoding="utf-8")
    # The Default config's settling baseline reads GOCART optics tables from
    # a path relative to the run directory, so every staged run needs them.
    # Copy exactly the filenames the staged config references out of the
    # shared optics fixture directory; a missing table aborts the C++ core
    # (fail-loud) and would silently freeze legacy settling, so fail here.
    if "settling" in processes:
        optics_source = source_root / "specs" / "tmp"
        optics_target = run_dir / "ExtData" / "monochromatic"
        optics_target.mkdir(parents=True, exist_ok=True)
        for filename in re.findall(r"(optics_\S+\.nc)", config_text):
            source = optics_source / filename
            if not source.is_file():
                raise SystemExit(f"settling optics fixture missing: {source}")
            shutil.copy2(source, optics_target / filename)
    profile = run_dir / "default_parity_met.json"
    subprocess.run(
        [sys.executable, str(source_root / "tests" / "build_parity_met_profile.py"),
         "--profile", str(source_root / "tests" / "MetProfiles" / "Profile_NCWCP.csv"),
         "--columns", str(columns), "--levels", str(levels), "--output", str(profile)],
        check=True,
    )
    return config, profile


def run(label: str, executable: Path, run_dir: Path, config: Path, profile: Path, steps: int,
        zero_initial_chemistry: bool) -> Path:
    snapshot = run_dir / f"{label}.json"
    command = [str(executable), "--config", str(config), "--met-profile", str(profile),
               "--snapshot", str(snapshot), "--steps", str(steps)]
    # The C++ candidate exposes this explicit control for the source-only
    # baseline; the legacy API driver already initializes chemistry to zero.
    if label == "candidate" and zero_initial_chemistry:
        command.extend(["--initial-chemistry", "zero"])
    subprocess.run(
        command,
        cwd=run_dir,
        check=True,
    )
    # The upstream legacy API integration driver writes this fixed-name
    # snapshot while it is being converted to the common runner interface.
    # Preserve the normalized driver contract for all other runners.
    legacy_fixed_snapshot = run_dir / "legacy_source_only.json"
    if label == "legacy" and not snapshot.is_file() and legacy_fixed_snapshot.is_file():
        shutil.copy2(legacy_fixed_snapshot, snapshot)
    if not snapshot.is_file():
        raise RuntimeError(f"{label} runner did not write {snapshot}")
    return snapshot


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy-runner", type=Path, required=True)
    parser.add_argument("--candidate-runner", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--columns", type=int, default=2,
                        help="columns in the shared profile (default: 2, ocean + dry land)")
    parser.add_argument("--levels", type=int, default=64)
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument("--zero-initial-chemistry", action="store_true",
                        help="match the legacy driver's zero-concentration source-only baseline")
    parser.add_argument("--processes", default="seasalt,dust",
                        help="comma-separated test1 process schedule (default: seasalt,dust)")
    parser.add_argument("--rtol", type=float, default=1.0e-12)
    parser.add_argument("--atol", type=float, default=1.0e-14)
    args = parser.parse_args()
    if args.columns < 1 or args.levels < 2 or args.steps < 1:
        raise SystemExit("columns and steps must be positive; levels must be at least two")
    processes = tuple(name.strip() for name in args.processes.split(",") if name.strip())
    if not processes:
        raise SystemExit("at least one process must be scheduled")

    legacy_config, legacy_profile = stage(args.workdir / "legacy", args.source_root, args.columns, args.levels, processes)
    candidate_config, candidate_profile = stage(args.workdir / "candidate", args.source_root, args.columns, args.levels,
                                                processes)
    legacy_snapshot = run("legacy", args.legacy_runner, args.workdir / "legacy", legacy_config, legacy_profile,
                          args.steps, args.zero_initial_chemistry)
    candidate_snapshot = run("candidate", args.candidate_runner, args.workdir / "candidate", candidate_config,
                             candidate_profile, args.steps, args.zero_initial_chemistry)
    return subprocess.run(
        [sys.executable, str(args.source_root / "tests" / "compare_parity_snapshots.py"),
         str(legacy_snapshot), str(candidate_snapshot), "--rtol", str(args.rtol), "--atol", str(args.atol)],
        check=False,
    ).returncode


if __name__ == "__main__":
    raise SystemExit(main())
