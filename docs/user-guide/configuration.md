# Configuration System

CATChem is configured through a small set of YAML files parsed by the C++
`ConfigManager` (`src/core/catchem_config_manager.{hpp,cpp}`). This page
describes the format the code actually reads.

!!! important
    There is **no** environment-variable substitution, no `!ENV`/`!CASE` tag
    support, and no `INHERIT`/file-merging in the current loader. Files are read
    verbatim with `YAML::LoadFile`, and an unrecognized tag such as `!ENV` makes
    the load fail. Any earlier documentation describing those features was
    aspirational and does not reflect the shipped code.

## Files and how they are loaded

A run is driven by one **main configuration file** plus two companion files
named from within it:

| File | Selected by | Loaded with |
| --- | --- | --- |
| Main config (`CATChem_new_config.yml`) | the driver / host | `load_from_file()` |
| Species / mechanism file | `simulation.species_filename` | `load_species_file()` |
| Emission mapping file | `simulation.emission_filename` | `load_emission_mapping_file()` |

`Core::initialize()` performs the sequence (see `src/core/catchem_core.cpp`):

1. `load_from_file(config_file)` — parse the main file and echo the effective
   YAML to `stdout`, so the run log records exactly what was parsed.
2. `load_species_file(...)` and `load_emission_mapping_file(...)` — the two
   companion paths are resolved **relative to the main config file's location**.
3. `validate_or_throw()` — strict schema and range validation (see
   "Validation summary"). A failure throws and aborts initialization.

The grid owned by a coupled host (column count, level count) overrides
`simulation/nx` and `simulation/ny` and `grid/number_of_levels`; for a
standalone run the YAML values win.

## Top-level keys

Strict validation accepts **only** these top-level keys; an unknown one is
reported as `"unknown top-level configuration key"`:

```
simulation  mechanism  physical_validation  grid  timesteps
diagnostics  mie  run_phases  processes  process
```

`process` is a legacy alias for `processes`; both are parsed into the same map.

## `simulation`

```yaml
simulation:
  name: test
  start_date: 20240501 0000      # stored as an opaque string
  end_date:   20240501 0100
  species_filename:  ./CATChem_species.yml     # companion file
  emission_filename: ./CATChem_emission.yml    # companion file
  verbose:
    activate: true               # -> data.simulation.verbose_enabled
    # log_level: debug           # optional; see "Logging" below
  # The keys below are optional runtime overrides
  nx: 4
  ny: 1
  nz: 64
  timestep: 3600                 # seconds; "dt" is accepted as a synonym
  nsteps: 1
```

- `name`, `start_date`, `end_date` are strings; the core does not parse the dates.
- `species_filename` / `emission_filename` point at the companion files.
- `nx`, `ny`, `nz`, `timestep`/`dt`, `nsteps` populate the runtime block used by
  standalone runs; a coupled host supplies its own grid instead.

## `grid` and `timesteps`

```yaml
grid:
  number_of_levels: 64           # also sets the runtime level count
  number_of_soil_layers: 4
timesteps:
  transport_timestep_in_s: 10
  chemistry_timestep_in_s: 60
```

## `physical_validation`

Controls how out-of-physical-range state is handled at runtime. Must be one of
three values; anything else throws during parse.

```yaml
physical_validation:
  policy: reject                 # reject | warn_and_clamp | count_and_continue
```

## `mechanism`

Optional metadata describing the chemical mechanism.

```yaml
mechanism:
  identity: "RADM2"                     # -> data.mechanism_identity
  capabilities: [photolysis, aqueous]   # -> data.mechanism_capabilities
```

## `diagnostics`

```yaml
diagnostics:
  output:
    enabled: true
    directory: "./output"
    prefix: "catchem_diag"
    frequency: 3600              # seconds between output writes
    format: "netcdf"             # only "netcdf" passes validation
    compress_lev: 2              # 0 = off, 1-9 = increasing compression
    diag_list: [so2, so4, dust1] # species to write; empty = all
  collection:
    enabled: true
    buffer_size: 1000
```

!!! warning
    `diagnostics.output.latlon_output` and the top-level `mie` block appear in
    some shipped example configs but are **not read by the current loader**.
    `mie` is accepted by validation, so it is harmless; `latlon_output` is simply
    ignored. Do not rely on either key to change behavior.

!!! deprecated
    `diagnostics.output/process_diagnostics` no longer suppresses process
    diagnostic output — it is a **deprecated no-op** retained only for
    configuration compatibility. Remove it from control files; it has no
    effect. Per-process diagnostics are now controlled entirely by each
    process's own `diagnostics:` flag (see *Process diagnostics* below).

### Process diagnostics (all schemes)

Setting `processes/<name>/diagnostics: true` registers that process's scheme
diagnostics; the optional `processes/<name>/diag_species:` list narrows them to
the named species (default = the process's own species subset, resolved at
runtime against the active mechanism). A `diag_species` name outside the
process's species set fails at initialization. Registered fields are written to
`catchem_diag*.nc` whenever runtime diagnostics are enabled — there is no
separate output switch.

```yaml
processes:
  settling:
    diagnostics: true
    diag_species: [so4, bc1, dust3, seas3]   # optional; default = process set
  dust:
    diagnostics: true                         # diag_species optional (dust bins only)
```

## `run_phases`

Defines the ordered process schedule. Each phase lists processes by name, and
the names must match keys under `processes`.

```yaml
run_phases:
  test1:
    description: "Test phase 1"
    processes:
      - seasalt
      - dust
      - carbchem
      - settling
      - drydep
      - so4chem
      - wetdep
```

## `processes` (and legacy `process`)

A map of process name to configuration block. The framework consumes these keys
directly: `activate`, `diagnostics`, `scheme`, `gas_scheme`, `aero_scheme`,
`diag_species`. Any **nested map** under a process is that scheme's option
block, and each leaf option is checked against the scheme's accepted options
during process initialization — so a typo fails loudly instead of silently
keeping a compiled default.

```yaml
processes:
  seasalt:
    activate: true
    diagnostics: true
    diag_species: []
    scheme: 'geos12'
    geos12:                      # scheme option block
      scale_factor: 1.0
      weibull_flag: false
  drydep:
    activate: true
    gas_scheme: 'wesely'         # drydep uses two scheme selectors
    aero_scheme: 'gocart'
    gocart:
      scale_factor: 1.0
      resuspension: false
  dust:
    activate: true
    scheme: 'fengsha'
    fengsha:
      alpha: 0.20
      gamma: 1.0
      drag_option: 1
```

### External emissions (`extemis`)

`extemis` groups named source categories (`anthro1`, `bio`, `fire`, `dust`,
`fengsha`, ...). The emission driver reads each category through the path-based
getters rather than the generic process parser, so category keys are not
validated against a fixed schema. Common category keys:

```yaml
processes:
  extemis:
    activate: true
    diagnostics: true
    global_factor: 1.0
    fire:
      activate: true
      scale_factor: 0.7778
      source_file: "ExtData/QFED/%y4/%m2/qfed2.emis_so2.006.%y4%m2%d2.nc4"
      format: "netcdf"           # "netcdf", or "volcano" for the point-source reader
      gridded: true
      is_2d: true
      regrid_method: conserve    # bilinear (default) | neareststod | conserve | patch | none
      time_interpolation: linear # "linear" enables time interpolation; otherwise none
      lat_name: "lat"
      lon_name: "lon"
      vertical_dist: "Ppbl"      # none | P100 | P500 | Ppbl | aviation* ...
      frequency: "daily"         # daily | monthly | hourly | static
      apply_method: "replace"    # add (default) | replace
      diagnostics: true
      diag_list: [biomass]
```

`source_file` supports the `%y4`/`%m2`/`%d2` date tokens (year, month, day),
expanded by the emission reader. This is the only templating the configuration
system supports.

## Companion files

### Species / mechanism file

A YAML **sequence**. Each entry has a `name` plus optional `__`-prefixed
attribute keys and physical metadata:

```yaml
- name: so2
  __description: Sulfur dioxide
  __is_gas: true
  __is_drydep: true
  __is_wetdep: true
  molecular weight [kg mol-1]: 64.04e-3
  __henry_k0: 1.22
  __henry_cr: 3100.0
  __dd_f0: 0.0
  __dd_hstar: 1.0e+5
```

- Species names must be unique after case normalization; duplicates are an error.
- Quote `NO` and similar values — YAML otherwise reads them as booleans.

### Emission mapping file

A map of **category → field → mapping**. Each field declares units, one scale
factor per mapped species, and the species (or `MET_`-namespace target) it feeds:

```yaml
anthro1:
  SO2:
    long_name: "Sulfur Dioxide"
    units: "kg/m2/s"
    scale: [0.97, 0.045]         # one entry per mapped species
    map: ["so2", "so4"]
```

Mapped targets that are not `MET_`-prefixed must exist in the active mechanism,
otherwise validation reports `"target is absent from active mechanism"`.

## Logging

Runtime log verbosity is resolved in priority order:

1. **YAML wins:** `simulation.verbose.log_level` sets the threshold through
   `Logger::set_level()`. Accepted values, case-insensitive: `debug`, `info`,
   `warn` (`warning` accepted), `error`. An unrecognized value throws during
   parse. When the key is absent, the current setting is kept.
2. **Environment fallback:** if no YAML level has been set, `CATCHEM_LOG_LEVEL`
   is used with the same accepted values, read once at first use.
3. **Default:** `info` when neither is provided, keeping production runs quiet.

```yaml
simulation:
  verbose:
    activate: true
    log_level: debug             # overrides CATCHEM_LOG_LEVEL
```

`Logger::debug()`, `info()`, `warn()`, and `error()` respect the threshold.
Guard expensive debug-only work with `Logger::enabled(Logger::Level::Debug)`.

## Accessing configuration at runtime

### C++ (`ConfigManager`)

Parsed values are available on `config.data` (`data.simulation`, `data.grid`,
`data.diagnostics`, `data.processes`, ...). Path queries use `/`-separated keys
and return the supplied default when a node is missing or cannot convert:

```cpp
catchem::ConfigManager config;
config.load_from_file("CATChem_new_config.yml");

bool sea_active = config.get_bool("processes/seasalt/activate", false);
std::string scheme = config.get_string("processes/seasalt/scheme", "");
int levels = config.get_int("grid/number_of_levels", 0);
auto diag_list = config.get_string_list("processes/extemis/anthro1/diag_list");
```

### Fortran / NUOPC

The `CATChem_Model` facade loads the main file during `initialize()` and exposes
typed accessors (`is_process_active`, `get_output_directory`,
`is_latlon_output_enabled`, ...). The emission and interface modules read
arbitrary YAML through the C-bound path getters
(`catchem_config_get_yaml_bool`, `_double`, `_int`, `_string`, `_list_count`,
`_list_at`).

```fortran
type(CATChem_Model) :: model
integer :: rc
call model%initialize("CATChem_new_config.yml", nx, ny, nz, rc=rc)
```

## Validation summary

`validate_or_throw()` reports an issue for any of:

- Root is not a mapping; unknown top-level key.
- `simulation/nx`, `simulation/ny`, `grid/number_of_levels`, or
  `simulation/nsteps` not positive; `simulation/timestep` outside `(0, 86400]`.
- Negative `grid/number_of_soil_layers`, or negative `timesteps/*_in_s`.
- `diagnostics/output.enabled` with a non-positive `frequency`, `compress_lev`
  outside `0..9`, or a `format` other than `netcdf`.
- `diagnostics/collection.enabled` with a non-positive `buffer_size`.
- Missing `species_filename`, an empty species list, or an empty/duplicate
  species name.
- Emission `scale`/`map` length mismatch, or a mapped non-`MET_` target absent
  from the mechanism.
- `physical_validation/policy` or `simulation/verbose/log_level` not a
  recognized value — these throw immediately during parse.

## Environment variables

Only two variables affect behavior, and neither modifies the configuration tree:

| Variable | Effect |
| --- | --- |
| `CATCHEM_LOG_LEVEL` | Logger threshold, used when no YAML `log_level` is set. |
| `NO_COLOR` | If non-empty, disables ANSI color in log output. |
