# Settling Process

**Process Type:** Deposition
**Description:** Process for computing gravitational settling of aerosol species
**Author:** Wei Li
**Generated:** 2025-12-18T14:12:33.301923

## Overview

The Settling process implements Process for computing gravitational settling of aerosol species. This process provides a modular, extensible framework for deposition calculations within the CATChem chemical transport model.

## Available Schemes

### GOCART Scheme

**Name:** `gocart`
**Description:** GOCART gravitational settling scheme
**Author:** Wei Li
**Reference:** GOCART2G process library `Chem_Settling` (metadata path) and
`Chem_SettlingSimple` (optics-table path) functions
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `scale_factor` | 1.0 |  -  | settling velocity factor |
| `simple_scheme` | False |  -  | Read wet-particle radius and density from NetCDF Mie optics tables when `true`; otherwise compute wet swelling internally from species metadata (Gerber/Fitzgerald). Both paths are supported by the C++ core and are numerically identical to `upstream/develop`. See [Optics tables](#optics-tables-simple_scheme-true) below. |
| `swelling_method` | 1 |  -  | method for calculating particle swelling: 1 Fitzgerald 1975; 2 for Gerber 1985 |
| `correction_maring` | False |  -  | correct the settling velocity following Maring et al, 2003 |

#### Required Meteorological Fields

- `T` - Meteorological field required for scheme computation
- `TSTEP` - Meteorological field required for scheme computation
- `AIRDEN` - Meteorological field required for scheme computation
- `RH` - Meteorological field required for scheme computation
- `Z` - Geometric height on the **interface** grid (n_levels + 1), required for scheme computation
- `PMID` - Layer mid-level pressure [Pa], required for scheme computation
- `DELP` - Layer pressure thickness [Pa], required for scheme computation

> The `gocart` scheme delegates to the upstream GOCART2G `Chem_Settling` /
> `Chem_SettlingSimple` kernels so results are numerically identical to
> `upstream/develop`. Missing or stale fields raise an explicit error before the
> kernel is invoked.

### Optics tables (`simple_scheme: true`)

The optics-table path reproduces the legacy GOCART2G wet-particle settling by
reading a per-aerosol-type **Mie lookup table** instead of computing swelling in
code. Configuration lives in a top-level `mie:` section, and each aerosol species
selects its table through its `__mie_name` attribute:

```yaml
mie:
  directory: "./ExtData/monochromatic/"   # joined with each file below
  files:
    DU: optics_DU.v15_5.nc    # dust   (bins 1-5)
    SS: optics_SS.v3_5.nc     # sea salt
    BC: optics_BC.v1_5.nc     # black carbon
    OC: optics_OC.v1_5.nc     # organic carbon
    SU: optics_SU.v1_5.nc     # sulfate (so4, msa)
    NI: optics_NI.v2_5.nc     # nitrate
    BRC: optics_BRC.v1_5.nc   # brown carbon

processes:
  settling:
    scheme: gocart
    gocart:
      simple_scheme: true
```

- **Type → table binding.** A species' `__mie_name` (e.g. `DU`) is matched
  (trimmed, case-sensitive) against the `mie.files` keys; the species then uses
  that table and selects its bin from the trailing digit of its short name
  (`dust4` → bin 4). This mirrors the legacy `chemstate_init_mie_data` map exactly.
- **Required table variables.** The `simple_scheme` kernel reads `rhop` (wet
  particle density) and `growth_factor` from the NetCDF table. Use the complete
  table versions that carry them (the trailing `_5` GOCART releases); a table
  missing `rhop` makes the reader fall back to its `-999.` sentinel and the
  species silently fail to settle.
- **Failure behaviour.** At `init` (before any time step) the process aborts
  loudly, naming the offending item, if: `simple_scheme: true` with an empty
  `mie.files`; the `mie.directory` is not readable; a listed file is missing; the
  NetCDF table fails to open; or a settling species has an empty/unmatched
  `__mie_name`. There is no silent fallback to the metadata path.
- **Parity.** The optics path is certified against the legacy Fortran oracle
  (`tests/run_settling_optics_parity.py`) to single-precision round-off; see
  `specs/012-settling-optics-parity/quickstart.md`.



## Process Interface

### Species

The settling process operates on the following chemical species:


### Required Inputs



### Process Diagnostics

| Diagnostic | Units | Description |
|------------|-------|-------------|
| `settling_velocity_per_species_per_level` | m/s | settling velocity per species per level |
| `settling_flux_per_species` | kg/m2/s | settling flux per species across column |

## Usage

### Basic Integration

```fortran
use SettlingProcessCreator_Mod
use SettlingCommon_Mod

! Create process instance
type(SettlingProcess_t) :: process
call create_settling_process(process, config_data)

! Use process in model time step
call process%run(state, dt)
```

### Scheme Selection

The process supports multiple schemes. Select your desired scheme:

```fortran
! Use GOCART scheme
process%scheme_name = "gocart"
```

## Implementation Details

### Pure Science Kernels

Each scheme is implemented as a pure science kernel with no infrastructure dependencies:

```fortran
! GOCART scheme
pure subroutine compute_gocart( &
   num_layers, num_species, params, &
   T, &   TSTEP, &   AIRDEN, &   RH, &   Z, &   PMID, &   DELP, &
   species_conc, emission_flux)
```

### Host Model Responsibilities

The host model (CATChem infrastructure) handles:

- Parameter initialization and validation
- Input array validation and error handling
- Memory management and array allocation
- Integration with model time stepping
- Diagnostic output management

## Configuration

### YAML Configuration Example

```yaml
processes:
  settling:
    enabled: true
    scheme: "gocart"
    parameters:
      scale_factor: 1.0
      simple_scheme: False
      swelling_method: 1
      correction_maring: False
    diagnostics:
      enabled: true
      output_frequency: "daily"
```

## Technical Specifications

- **Parallelization:** Column
- **Memory Requirements:** Low
- **Timestep Dependency:** Independent
- **Multiphase Support:** No
- **Size Bin Support:** No
- **Vectorization:** Supported

## Files Generated

### Source Code
- `src/process/settling/ProcessSettlingInterface_Mod.F90` - Main process interface
- `src/process/settling/SettlingCommon_Mod.F90` - Common types and parameters
- `src/process/settling/SettlingProcessCreator_Mod.F90` - Process factory
- `src/process/settling/schemes/SettlingScheme_GOCART_Mod.F90` - GOCART gravitational settling scheme

### Tests
- `tests/process/settling/unit/` - Unit tests
- `tests/process/settling/integration/` - Integration tests

### Documentation
- `docs/processes/settling/settling.md` - This documentation

## Contributing

When modifying or extending this process:

1. **Science Changes:** Modify the scheme modules in `schemes/`
2. **Interface Changes:** Update the main interface module
3. **New Schemes:** Add new scheme modules and update the creator
4. **Tests:** Add corresponding unit and integration tests
5. **Documentation:** Update this documentation file

## References

- GOCART: GOCART2G process library Chem_SettlingSimple function

---
*This documentation was automatically generated by the CATChem Process Generator on 2025-12-18T14:12:33.301923*
