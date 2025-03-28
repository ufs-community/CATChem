#!/usr/bin/env python3
"""
Extract a column of input data from various sources.
"""
import re
from pathlib import Path
from textwrap import indent

import numpy as np
import pandas as pd
import xarray as xr
import yaml

HERE = Path(__file__).parent

with open(HERE / "vars.yml") as f:
    ctl = yaml.safe_load(f)
    var_info = ctl["variables"]

unique_src_ids = {d["src"].split(":")[0] for d in var_info.values()}
print("srcs:", unique_src_ids)

for d in ctl["cases"]:
    lat = d["lat"]
    lon = d["lon"]
    nz = d.get("nz")
    atm_fp = Path(d["atm"])
    sfc_fp = Path(d["sfc"])
    out_fp = Path(d["out"])
    if not out_fp.is_absolute():
        out_fp = HERE / out_fp
    print(f"lat={lat}, lon={lon}, nz={nz}, atm={atm_fp}, sfc={sfc_fp}, out={out_fp}")

    src = {}

    src["atm"] = None
    if "atm" in unique_src_ids:
        ds = (
            xr.open_dataset(atm_fp)
            .drop_vars(["lon", "lat"])
            .rename({"grid_xt": "lon", "grid_yt": "lat"})
            .sel(lat=lat, lon=lon, method="nearest")
            .squeeze()
        )

        # Switch to sfc first
        pmid, p = ds["pfull"].values, ds["phalf"].values
        assert pmid[0] < pmid[-1] and p[0] < p[-1]
        ds = ds.isel(pfull=slice(None, None, -1), phalf=slice(None, None, -1))

        # Compute height
        assert "z" not in ds.variables
        assert "zmid" not in ds.variables
        zsfc = ds["hgtsfc"]
        dz = ds["delz"]
        if (dz < 0).all():
            dz *= -1
        ztop = (zsfc + dz).cumsum("pfull")
        z = xr.DataArray(
            data=np.concatenate(([zsfc.values], ztop.values)),
            dims="phalf",
            # coords={"phalf": ds.coords["phalf"]},
        )
        zmid = ((z + z.shift(phalf=1)) / 2).isel(phalf=slice(1, None)).swap_dims(phalf="pfull")
        ds = ds.assign_coords(z=z, zmid=zmid)

        # Select levels
        if nz is not None:
            ds = ds.isel(pfull=slice(None, nz), phalf=slice(None, nz + 1))

        # Store
        src["atm"] = ds

    src["sfc"] = None
    if "sfc" in unique_src_ids:
        ds = (
            xr.open_dataset(sfc_fp)
            .drop_vars(["lon", "lat"])
            .rename({"grid_xt": "lon", "grid_yt": "lat"})
            .sel(lat=lat, lon=lon, method="nearest")
            .squeeze()
        )

        # Switch to sfc first
        pmid, p = ds["pfull"].values, ds["phalf"].values
        assert pmid[0] < pmid[-1] and p[0] < p[-1]
        ds = ds.isel(pfull=slice(None, None, -1), phalf=slice(None, None, -1))

        # Select levels (for cloud fraction)
        if nz is not None:
            ds = ds.isel(pfull=slice(None, nz), phalf=slice(None, nz + 1))

        # Compute soil moisture variables
        #
        # soilw{1..4}: volumetric (fraction) soil moisture at 10, 30, 60, 100 cm
        # Noah predicts values at layer midpoints:
        # - [0, 10)    5
        # - [10, 40)   25
        # - [40, 100)  70
        # - [100, 200) 150
        #
        # GWETTOP: 0--5 cm avg?
        # GWETROOT: 10--100 cm avg?
        ds["gwettop"] = ds["soilw1"]  # FIXME: this is value at 5 cm, not 0--5 cm avg
        ds["gwetroot"] = (ds["soilw2"] * 30 + ds["soilw3"] * 60) / 90

        # Combine soil moisture into one array
        soilws = []
        for vn in ds.data_vars:
            m = re.fullmatch(r"soilw([0-9]+)", vn)
            if m is not None:
                lev = int(m.group(1))
                soilws.append((vn, lev))
        if soilws:
            soilws.sort(key=lambda x: x[1])
            ds["soilm"] = xr.concat([ds[vn] for vn, _ in soilws], dim="soil")
        else:
            print("warning: soilw variables not identified")

        # sotyp, land, vtype are stored as float, maybe to support null mask?
        # But they seem to be all non-null, though with zeros (land 53.8%, sotyp/vtype 66.3%)
        for vn in ["dluse", "sotyp", "land", "vtype"]:
            if vn in ds.variables:
                ds[vn] = ds[vn].astype(int)

        # Store
        src["sfc"] = ds

    das = []
    for vn, d in var_info.items():
        if vn in {"z", "zmid"}:
            assert d["src"] == "atm", "diagnosed"
            src_id, src_vn = "atm", vn
        elif vn in {"gwettop", "gwetroot", "soilm"}:
            assert d["src"] == "sfc", "diagnosed"
            src_id, src_vn = "sfc", vn
        else:
            src_id, src_vn = d["src"].split(":")

        da = src[src_id][src_vn]

        # Maybe convert units
        f = d.get("f")
        if f is not None:
            da *= f

        das.append(da.rename(vn))

    # Write to text file
    with open(out_fp, "w") as f:
        for da in das:
            if pd.api.types.is_float_dtype(da.dtype):
                fmt = ".4e"
            elif pd.api.types.is_integer_dtype(da.dtype):
                fmt = "d"
            else:
                raise AssertionError
            if da.ndim == 0:
                data = f"{da.item():{fmt}}"
            elif da.ndim == 1:
                data = ",".join(f"{x:{fmt}}" for x in da.values)
            else:
                raise AssertionError
            info = ",".join([da.name, str(da.size)])
            f.write(info + "\n")
            f.write(data + "\n")

# Modify testing module

fp = HERE / "../tests/testing_mod.f90"
with open(fp) as f:
    lines = f.readlines()

read_var_tpl = """\
case ("{vn}")
  read(unum, *, iostat=ios) MetState%{upper_vn}
  if (ios /= 0) then
     print *, "Error reading {upper_vn}:", ios
     rc = 1
     return
  end if
"""
read_var_indent = 3 * 4 + 1
read_var_title = "read column data"

print_var_tpl = """\
print *, "{upper_vn}:", MetState%{upper_vn}
"""
print_var_arr_tpl = """\
print *, "{upper_vn}:", size(MetState%{upper_vn}), MetState%{upper_vn}
"""
print_var_indent = 3 * 3
print_var_title = "print column data"

read_var_blocks = []
print_var_lines = []
for da in sorted(das, key=lambda x: x.name.lower()):
    vn = da.name
    upper_vn = vn.upper()
    skip = var_info[vn].get("skip", False)
    if skip:  # not in MetState yet
        continue

    read_var_blocks.append(read_var_tpl.format(vn=vn, upper_vn=upper_vn))

    if da.ndim == 0:
        line = print_var_tpl.format(upper_vn=upper_vn)
    elif da.ndim == 1:
        line = print_var_arr_tpl.format(upper_vn=upper_vn, size=da.size)
    else:
        raise AssertionError
    print_var_lines.append(line)

read_var_str = "".join(indent(block, " " * read_var_indent) for block in read_var_blocks)
print_var_str = "".join(indent(line, " " * print_var_indent) for line in print_var_lines)

lines_strip = [line.strip() for line in lines]
a = lines_strip.index(f"! >>> {read_var_title} >>>")
b = lines_strip.index(f"! <<< {read_var_title} <<<")
assert b > a
lines = lines[:a + 1] + [read_var_str] + lines[b:]

lines_strip = [line.strip() for line in lines]
a = lines_strip.index(f"! >>> {print_var_title} >>>")
b = lines_strip.index(f"! <<< {print_var_title} <<<")
assert b > a
lines = lines[:a + 1] + [print_var_str] + lines[b:]

with open(fp, "w") as f:
    f.writelines(lines)
