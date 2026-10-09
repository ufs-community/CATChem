#!/usr/bin/env python3
"""Compare final UFS tracer extrema against a legacy reference log."""

import argparse
import math
import re
import sys
from pathlib import Path

REQUIRED = ("so2", "so4", "dms", "msa", "bc1", "bc2", "oc1", "oc2", "dust1", "dust2", "dust3", "dust4", "dust5", "seas1", "seas2", "seas3", "seas4", "seas5", "pm25", "pm10")
LINE = re.compile(r"^\s*\d+:\s+(\w+)\s+max\s+=\s*([-+0-9.Ee]+)\s+min\s+=\s*([-+0-9.Ee]+)")


def parse_extrema(path: Path) -> dict[str, tuple[float, float]]:
    values: dict[str, tuple[float, float]] = {}
    for line in path.read_text().splitlines():
        match = LINE.match(line)
        if match:
            values[match.group(1).lower()] = (float(match.group(2)), float(match.group(3)))
    return values


def relative_difference(current: float, reference: float) -> float:
    if reference == 0.0:
        return 0.0 if current == 0.0 else math.inf
    return abs(current - reference) / abs(reference)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("legacy", type=Path)
    parser.add_argument("current", type=Path)
    parser.add_argument("--max-relative-difference", type=float, default=10.0)
    args = parser.parse_args()
    legacy, current = parse_extrema(args.legacy), parse_extrema(args.current)
    missing = [name for name in REQUIRED if name not in legacy or name not in current]
    if missing:
        print("missing required extrema: " + ", ".join(missing), file=sys.stderr)
        return 2
    failed = False
    print("species       legacy max       current max     absolute difference     relative difference")
    for name in REQUIRED:
        absolute = abs(current[name][0] - legacy[name][0])
        difference = relative_difference(current[name][0], legacy[name][0])
        print(f"{name:8s} {legacy[name][0]:16.7g} {current[name][0]:16.7g} {absolute:23.7g} {difference:20.7g}")
        failed |= difference > args.max_relative_difference
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
