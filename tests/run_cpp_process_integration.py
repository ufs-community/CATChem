#!/usr/bin/env python3
"""Run one Default-config process through the C++ core integration boundary.

This is the C++ analogue of the legacy per-process integration tests.  It
stages the checked-in Default configuration, shared MET profile, species
metadata, and a single-process phase, then validates the normalized C++ state
snapshot.  Numerical legacy-vs-C++ comparison remains a separate parity test.
"""

import argparse
import json
import math
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from run_default_numerical_parity import stage


SOURCE_PROCESSES = {"seasalt", "dust"}
ALL_PROCESSES = {"seasalt", "dust", "carbchem", "settling", "drydep", "so4chem", "wetdep"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-runner", required=True, type=Path)
    parser.add_argument("--process", required=True, action="append", choices=sorted(ALL_PROCESSES))
    parser.add_argument("--seasalt-scheme", choices=("gong97", "gong03", "geos12"),
                        help="override the Default sea-salt scheme for a one-process case")
    parser.add_argument("--drop-field", default=None,
                        help="remove a met_3d field from the staged profile to exercise the "
                             "process's missing-input guard (expects nonzero exit naming the field)")
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    if args.seasalt_scheme and args.process != ["seasalt"]:
        parser.error("--seasalt-scheme requires exactly one --process seasalt")
    if args.drop_field and len(args.process) != 1:
        parser.error("--drop-field requires exactly one --process")

    for process in args.process:
        with tempfile.TemporaryDirectory(prefix=f"catchem-cpp-{process}-") as directory:
            run_dir = Path(directory)
            config, profile = stage(run_dir, args.source_root, columns=2, levels=64, processes=(process,))
            if args.seasalt_scheme:
                # Upstream's sea-salt integration test uses this synthetic
                # 10-column, 20-level ocean profile and a 3600 s step.
                profile = run_dir / "seasalt_integration_met.json"
                subprocess.run(["python3", str(args.source_root / "tests" / "build_seasalt_integration_profile.py"),
                                "--output", str(profile)], check=True)
                config_text, replacements = re.subn(
                    r"(  seasalt:\n(?:    .*\n)*?    scheme:)\s*'?(?:gong97|gong03|geos12)'?",
                    rf"\1 '{args.seasalt_scheme}'", config.read_text(encoding="utf-8"), count=1)
                if replacements != 1:
                    raise RuntimeError("could not select the requested sea-salt scheme")
                config.write_text(config_text, encoding="utf-8")
            if args.drop_field:
                # Negative case: remove a host-supplied MET field so the
                # process's missing-input guard must abort with a named error.
                fields = json.loads(profile.read_text(encoding="utf-8"))
                section = "met_3d" if args.drop_field in fields.get("met_3d", {}) else "met_interface"
                if args.drop_field not in fields.get(section, {}):
                    raise RuntimeError(f"staged profile has no {section} field {args.drop_field}")
                del fields[section][args.drop_field]
                profile.write_text(json.dumps(fields), encoding="utf-8")
            snapshot = run_dir / "candidate.json"
            command = [str(args.candidate_runner), "--config", str(config), "--met-profile", str(profile),
                       "--snapshot", str(snapshot), "--steps", "1"]
            if args.seasalt_scheme:
                command.extend(["--dt", "3600"])
            if process in SOURCE_PROCESSES:
                command.extend(["--initial-chemistry", "zero"])
            if args.drop_field:
                result = subprocess.run(command, cwd=run_dir, capture_output=True, text=True)
                combined = result.stdout + result.stderr
                if result.returncode == 0:
                    raise RuntimeError(f"--drop-field {args.drop_field} unexpectedly succeeded for {process}")
                if args.drop_field not in combined or process.lower() not in combined.lower():
                    raise RuntimeError(f"--drop-field {args.drop_field} did not abort with a named "
                                       f"{process}/{args.drop_field} error for {process}")
                print(f"PASS: C++ {process} integration aborts without {args.drop_field}")
                continue
            subprocess.run(command, cwd=run_dir, check=True)

            data = json.loads(snapshot.read_text(encoding="utf-8"))
            values = data["snapshots"][0]["fields"]["concentration"]["values"]
            if not values or not all(math.isfinite(value) for value in values):
                raise RuntimeError(f"{process} produced missing or non-finite concentrations")

            # Emission schemes must produce a signal from the source-active MET
            # columns even when all chemistry begins at zero.
            if process == "seasalt" and not any(value != 0.0 for value in values):
                raise RuntimeError(f"{process} produced no source response")
        print(f"PASS: C++ {process} integration")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
