#!/usr/bin/env python3
"""Compare normalized legacy and C++ CATChem parity snapshots.

Each JSON file must contain a top-level object with ``schema_version: 1`` and
``snapshots``.  Each snapshot has a unique ``checkpoint`` string and a
``fields`` object.  Field values are flat numeric arrays; dimensions, units,
and species order are part of each field's metadata and must match exactly.
"""

import argparse
import json
import math
import sys
from pathlib import Path


def load_snapshot(path: Path):
    with path.open(encoding="utf-8") as stream:
        result = json.load(stream)
    if result.get("schema_version") != 1 or not isinstance(result.get("snapshots"), list):
        raise ValueError(f"{path}: expected schema_version=1 and a snapshots array")
    return result


def index_snapshots(document, label):
    indexed = {}
    for snapshot in document["snapshots"]:
        checkpoint = snapshot.get("checkpoint")
        fields = snapshot.get("fields")
        if not isinstance(checkpoint, str) or not isinstance(fields, dict):
            raise ValueError(f"{label}: every snapshot needs checkpoint and fields")
        if checkpoint in indexed:
            raise ValueError(f"{label}: duplicate checkpoint {checkpoint!r}")
        indexed[checkpoint] = fields
    return indexed


def fail(message):
    print(f"PARITY FAILURE: {message}", file=sys.stderr)
    return 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("legacy", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--rtol", type=float, default=1.0e-12)
    parser.add_argument("--atol", type=float, default=1.0e-14)
    args = parser.parse_args()
    if args.rtol < 0 or args.atol < 0:
        return fail("tolerances must be non-negative")

    try:
        legacy = index_snapshots(load_snapshot(args.legacy), "legacy")
        candidate = index_snapshots(load_snapshot(args.candidate), "candidate")
    except (OSError, ValueError, json.JSONDecodeError) as error:
        return fail(str(error))

    if legacy.keys() != candidate.keys():
        return fail(f"checkpoint mismatch: legacy={sorted(legacy)}, candidate={sorted(candidate)}")

    worst = (0.0, "")
    for checkpoint, legacy_fields in legacy.items():
        candidate_fields = candidate[checkpoint]
        if legacy_fields.keys() != candidate_fields.keys():
            return fail(f"{checkpoint}: field mismatch")
        for name, old in legacy_fields.items():
            new = candidate_fields[name]
            metadata = ("units", "shape", "species")
            if any(old.get(key) != new.get(key) for key in metadata):
                return fail(f"{checkpoint}/{name}: metadata mismatch")
            old_values, new_values = old.get("values"), new.get("values")
            if not isinstance(old_values, list) or not isinstance(new_values, list):
                return fail(f"{checkpoint}/{name}: values must be arrays")
            if len(old_values) != len(new_values):
                return fail(f"{checkpoint}/{name}: value count mismatch")
            for index, (old_value, new_value) in enumerate(zip(old_values, new_values)):
                if not isinstance(old_value, (int, float)) or not isinstance(new_value, (int, float)):
                    return fail(f"{checkpoint}/{name}[{index}]: non-numeric value")
                if not math.isfinite(old_value) or not math.isfinite(new_value):
                    return fail(f"{checkpoint}/{name}[{index}]: non-finite value")
                error = abs(new_value - old_value)
                limit = args.atol + args.rtol * max(abs(old_value), abs(new_value))
                if error > limit:
                    return fail(
                        f"{checkpoint}/{name}[{index}]: legacy={old_value:.17g}, "
                        f"candidate={new_value:.17g}, abs_error={error:.3e}, limit={limit:.3e}"
                    )
                if error > worst[0]:
                    worst = (error, f"{checkpoint}/{name}[{index}]")
    print(f"PASS: numerical parity; maximum absolute error {worst[0]:.3e} at {worst[1] or 'n/a'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
