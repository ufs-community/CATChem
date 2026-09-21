#!/usr/bin/env python3
"""Create the synthetic 10-column, 20-level MET case used upstream sea-salt tests."""

import argparse
import json
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    columns, levels = 10, 20
    latitudes = [-30.0] * columns
    longitudes = [-120.0 + index * 240.0 / (columns - 1) for index in range(columns)]
    wind = 8.0 + 2.0 * math.cos(math.radians(-30.0))
    delp = [10000.0 * math.exp(-(index * 20.0 / (levels - 1)) / 8.0) for index in range(levels)]
    pedge = [101325.0]
    for thickness in delp:
        pedge.append(pedge[-1] - thickness)
    data = {
        "schema_version": 1,
        "grid": {"columns": columns, "levels": levels},
        "met_2d": {
            "LAT": latitudes, "LON": longitudes, "FROCEAN": [1.0] * columns,
            "FRSEAICE": [0.0] * columns,
            "TS": [298.0 + 5.0 * math.cos(math.radians(-30.0))] * columns,
            "U10M": [-wind * 0.8] * columns, "V10M": [wind * 0.3] * columns,
            "USTAR": [0.03 * math.sqrt((wind * 0.8) ** 2 + (wind * 0.3) ** 2)] * columns,
        },
        "met_3d": {"DELP": delp, "AIRDEN": [1.2 * math.exp(-(index * 20.0 / (levels - 1)) / 8.0)
                                                    for index in range(levels)]},
        "met_interface": {"PEDGE": pedge}, "met_soil": {},
    }
    args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
