#!/usr/bin/env python3
"""Create the deterministic meteorology input shared by parity runners.

The checked-in MET profile is a 127-level sounding, whereas the Default
CATChem configuration has 64 levels.  Both the legacy and C++ runners consume
this generated JSON, which makes the remap and derived AIRDEN field identical
before either core is entered.
"""

import argparse
import json
import math
from pathlib import Path


THREE_D_FIELDS = ("BXHEIGHT", "DELP", "PMID", "ZMID", "CLDF", "QL", "QV", "SPHU", "T", "U", "V")
TWO_D_FIELDS = ("EFLUX", "FRSNO", "HFLUX", "LWGDN", "PBLH", "PS", "SNODP", "SNOMAS", "SWGDN", "T2M", "TS", "U10M", "USTAR", "V10M", "WILT", "Z0", "GWETTOP", "LWI")

# Required by the complete Default process schedule but not present in the
# one-column atmospheric sounding. Values are deterministic neutral-land
# conditions; each runner receives these exact values.
SURFACE_DEFAULTS = {
    "CLAYFRAC": 0.25, "DLUSE": 1.0, "FROCEAN": 1.0, "FRLAKE": 0.0,
    "FRSEAICE": 0.0, "GVF": 0.5, "LAI": 2.0, "LAT": 40.0, "LON": -75.0,
    "RDRAG": 0.1, "SNDFRC": 0.5, "SSM": 0.15, "USTAR_THRESHOLD": 0.2,
    "SALINITY": 35.0, "Z0H": 0.0001,
}


def source_active_surface(columns: int) -> dict[str, list[float]]:
    """Return deterministic ocean/land source conditions for parity runs.

    Even-numbered columns are ocean points for sea-salt production; odd-numbered
    columns are dry, bare-land points for windblown dust.  A two-column run
    therefore exercises both source schemes without making either process
    infer a source regime from external inventory data.
    """
    ocean = [index % 2 == 0 for index in range(columns)]
    return {
        "FROCEAN": [1.0 if value else 0.0 for value in ocean],
        "FRSEAICE": [0.0] * columns,
        "U10M": [15.0] * columns,
        "V10M": [0.0] * columns,
        "USTAR": [0.65] * columns,
        "CLAYFRAC": [0.0 if value else 0.30 for value in ocean],
        "SNDFRC": [0.0 if value else 0.60 for value in ocean],
        "SSM": [0.15 if value else 0.05 for value in ocean],
        "GWETTOP": [0.15 if value else 0.01 for value in ocean],
        "USTAR_THRESHOLD": [0.20] * columns,
        "LAI": [2.0 if value else 0.0 for value in ocean],
        "GVF": [0.5 if value else 0.0 for value in ocean],
        "LWI": [0.0 if value else 1.0 for value in ocean],
    }


def read_blocks(path: Path) -> dict[str, list[float]]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    blocks: dict[str, list[float]] = {}
    index = 0
    while index < len(lines):
        name, count = lines[index].split(",", 1)
        count_i = int(count)
        index += 1
        values: list[float] = []
        while len(values) < count_i and index < len(lines):
            values.extend(float(value) for value in lines[index].split(",") if value)
            index += 1
        if len(values) != count_i:
            raise ValueError(f"{name}: expected {count_i} values, found {len(values)}")
        blocks[name.lower()] = values
    return blocks


def remap(values: list[float], levels: int) -> list[float]:
    """Linearly interpolate cell-centre data by normalized vertical index."""
    if len(values) == levels:
        return values
    return [
        values[min(len(values) - 1, math.floor(position))] * (1.0 - (position % 1.0))
        + values[min(len(values) - 1, math.ceil(position))] * (position % 1.0)
        for position in (index * (len(values) - 1) / (levels - 1) for index in range(levels))
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", type=Path, default=Path("tests/MetProfiles/Profile_NCWCP.csv"))
    parser.add_argument("--levels", type=int, default=64)
    parser.add_argument("--columns", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.levels < 2 or args.columns < 1:
        raise SystemExit("levels must be >= 2 and columns must be >= 1")

    blocks = read_blocks(args.profile)
    required = {name.lower() for name in THREE_D_FIELDS + TWO_D_FIELDS}
    missing = sorted(required - blocks.keys())
    if missing:
        raise SystemExit(f"profile is missing required fields: {', '.join(missing)}")

    fields_3d = {name: remap(blocks[name.lower()], args.levels) for name in THREE_D_FIELDS}
    # Profile humidity is g/kg.  Use kg/kg in the moist-air ideal-gas law.
    qv = [max(0.0, value * 1.0e-3) for value in fields_3d["QV"]]
    fields_3d["QV"] = qv
    fields_3d["SPHU"] = qv[:]
    # PMID is hPa in Profile_NCWCP; AIRDEN is kg m-3.
    fields_3d["PMID"] = [pressure * 100.0 for pressure in fields_3d["PMID"]]
    fields_3d["AIRDEN"] = [pressure / (287.05 * temperature * (1.0 + 0.61 * humidity))
                            for pressure, temperature, humidity in zip(fields_3d["PMID"], fields_3d["T"], qv)]
    fields_3d["AIRDEN_DRY"] = [pressure / (287.05 * temperature)
                                for pressure, temperature in zip(fields_3d["PMID"], fields_3d["T"])]
    fields_3d["MAIRDEN"] = fields_3d["AIRDEN"][:]
    fields_3d["RH"] = [min(1.0, max(0.0, humidity / 0.02)) for humidity in qv]
    fields_3d["REEVAPLS"] = [0.0] * args.levels

    # Construct a strictly positive, descending pressure interface from the
    # remapped midpoint pressure. Profile DELP is retained as an imported
    # field, but both cores derive their execution DELP from this interface.
    pmid = fields_3d["PMID"]
    pedge = [pmid[0] * 1.02]
    pedge.extend(math.sqrt(lower * upper) for lower, upper in zip(pmid[:-1], pmid[1:]))
    pedge.append(pmid[-1] * 0.5)

    document = {
        "schema_version": 1,
        "source_profile": str(args.profile),
        "grid": {"columns": args.columns, "levels": args.levels},
        "met_3d": fields_3d,
        "met_interface": {"PEDGE": pedge, "Z": remap(blocks["z"], args.levels + 1),
                          "PFILSAN": [0.0] * (args.levels + 1), "PFLLSAN": [0.0] * (args.levels + 1)},
        # A dry profile activates Fengsha at the land column.  Sea-salt does
        # not consume soil moisture at the companion ocean column.
        # Axis-2 categorical fields are imported in the same layout as
        # legacy MetState (column, singleton, land-use category).
        "met_soil": {
            "SOILM": [0.01] * 4,
            "FRLANDUSE": [1.0] + [0.0] * 19,
            "FRLAI": [2.0] + [0.0] * 19,
            "ILAND": list(range(1, 21)),
        },
        "met_2d": {
            **{name: [blocks[name.lower()][0]] * args.columns for name in TWO_D_FIELDS},
            **{name: [value] * args.columns for name, value in SURFACE_DEFAULTS.items()},
            **source_active_surface(args.columns),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
