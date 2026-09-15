# Standalone Driver

The **standalone driver** lets you run CATChem as a self-contained ESMF/NUOPC
application, without a host model (such as UFS/FV3) providing the grid, clock,
or meteorology. It is the easiest way to exercise CATChem end-to-end for
development, testing, and small offline experiments.

In a coupled configuration the atmospheric driver supplies the ESMF grid and
imports/exports fields. Running standalone, the driver creates its own grid and
clock and reads any required meteorology from offline files through the emission
reader.

## Component hierarchy

```text
catchem_app (main program)
     │
catchem_driver (NUOPC_Driver)     ← owns the clock, builds the grid,
     │                              adds the cap as a child component
cc_nuopc (CATChem NUOPC cap / NUOPC_Model)
```

The application is intentionally minimal: the driver (`catchem_driver`) owns the
top-level clock and adds the CATChem cap (`cc_nuopc`) as a child. Additional
components (a data/forcing component, a mediator, ...) can be added later
without changing the application.

## Building

The standalone driver is built by the CMake option
`CATCHEM_BUILD_STANDALONE`, which **requires** the NUOPC cap
(`CATCHEM_BUILD_NUOPC`):

```bash
cmake -S . -B build \
  -DCATCHEM_BUILD_NUOPC=ON \
  -DCATCHEM_BUILD_STANDALONE=ON
cmake --build build
```

This produces the executable `build/bin/catchem_app`. CMake also copies the
example driver configure file and the CATChem science YAML files (from
`tests/Configs/Default`) into the build root so a run can find them.

!!! note
    Both options default to `OFF`. `CATCHEM_BUILD_STANDALONE` on its own does
    nothing unless `CATCHEM_BUILD_NUOPC` is also enabled.

## Running

```bash
cd build
mpirun -np 1 ./bin/catchem_app catchem_standalone.configure
```

- The single command-line argument is the path to the driver configure file.
  If omitted, `catchem_standalone.configure` (in the current working directory)
  is used.
- File paths inside the configure file (the CATChem YAML, and the species /
  emission files it references) are resolved **relative to the current working
  directory**, which is why the run is launched from the build root where those
  files were copied.

## Driver configure file

The driver is controlled by an **ESMF_Config** file (default
`catchem_standalone.configure`). It sets the clock, the grid, and the path to
the CATChem *science* configuration. The science configuration itself
(processes, species, emissions, ...) stays in its own YAML file, referenced by
`catchem_config_file:`.

### Clock

| Key                | Type    | Default               | Description                                   |
| ------------------ | ------- | --------------------- | --------------------------------------------- |
| `start_time:`      | string  | `2020-01-01T00:00:00` | Start time (ISO 8601, `YYYY-MM-DDThh:mm:ss`). |
| `stop_time:`       | string  | `2020-01-01T06:00:00` | Stop time (ISO 8601).                         |
| `timestep_seconds:`| integer | `3600`                | Model time step in seconds.                   |

### Science configuration files

| Key                    | Type   | Default                  | Description                                                                 |
| ---------------------- | ------ | ------------------------ | --------------------------------------------------------------------------- |
| `catchem_config_file:` | string | `CATChem_new_config.yml` | Path to the CATChem YAML (processes, species, emissions, diagnostics).      |
| `field_mapping_file:`  | string | *(unset)*                | Only needed for coupled runs that exchange fields. Leave unset when standalone. |

### Grid

| Key               | Type    | Default  | Description                                                        |
| ----------------- | ------- | -------- | ------------------------------------------------------------------ |
| `grid_mode:`      | string  | `column` | `column` (single 1×1 column) or `gridded` (regular lat-lon).       |
| `grid_nz:`        | integer | `72`     | Number of vertical levels (an ungridded field dimension).          |
| `column_lon:`     | real    | `0.0`    | Column longitude (deg), used when `grid_mode: column`.             |
| `column_lat:`     | real    | `0.0`    | Column latitude (deg), used when `grid_mode: column`.              |
| `column_dlon:`    | real    | `0.0`    | Nominal column cell width in longitude (deg). See caveat below.    |
| `column_dlat:`    | real    | `0.0`    | Nominal column cell width in latitude (deg). See caveat below.     |
| `grid_nx:`        | integer | `1`      | Number of cells in longitude, used when `grid_mode: gridded`.      |
| `grid_ny:`        | integer | `1`      | Number of cells in latitude, used when `grid_mode: gridded`.       |
| `grid_lon_start:` | real    | `0.0`    | Western edge (deg), gridded mode.                                  |
| `grid_lon_end:`   | real    | `360.0`  | Eastern edge (deg), gridded mode.                                  |
| `grid_lat_start:` | real    | `-90.0`  | Southern edge (deg), gridded mode.                                 |
| `grid_lat_end:`   | real    | `90.0`   | Northern edge (deg), gridded mode.                                 |

!!! warning "Column corners and grid-cell area"
    A single column is a point and has no intrinsic horizontal size. Corner
    coordinates are created **only when both** `column_dlon` and `column_dlat`
    are greater than `0`. If either is `0` (the default), the column grid has
    no cell corners, so **conservative regridding** and **grid-cell area**
    (`AREA_M2`, needed for point-source emissions) are unavailable. In that case
    use `neareststod`/`bilinear` regridding for emission inputs, or set both
    widths `> 0` (e.g. to match your emission input resolution) to enable
    conservative regridding and `AREA_M2`.

### Run sequence (optional)

With a single component the default run sequence is sufficient, so no `runSeq::`
block is required. Add one when you introduce more components:

```text
runSeq::
  @3600
    CATCHEM
  @
::
```

## Meteorology and inputs

A standalone run has no host model to supply meteorology. Any required met
fields are read from offline files through the emission reader: emission mapping
targets prefixed with `MET_`/`met_` are written into the MetState, and the
remaining pressure-derived fields (`DELP`, `AIRDEN`, ...) are derived
automatically. See the [Input Files](input-files.md) and
[Configuration System](configuration.md) pages for the YAML details.

## Output

Diagnostic output is controlled by the CATChem YAML (`diagnostics/output`), not
by the driver configure file. Output is written to the configured directory
(default `./output`) as time-stamped NetCDF files. See the
[Diagnostic System](../core-concepts/diagnostics.md) for details.

## Example

A minimal single-column run:

```text
# catchem_standalone.configure
start_time:        2020-01-01T00:00:00
stop_time:         2020-01-01T06:00:00
timestep_seconds:  3600

catchem_config_file:  CATChem_new_config.yml

grid_mode:  column
grid_nz:    72
column_lon: 0.0
column_lat: 0.0
```

```bash
cd build
mpirun -np 1 ./bin/catchem_app catchem_standalone.configure
```

## See also

- [Build System](build-system.md) — full CMake option reference.
- [NUOPC Integration](../developer-guide/integration/nuopc.md) — the driver
  component architecture and how to extend it with more components.
