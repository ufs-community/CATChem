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
`Chem_SettlingSimple` (optics-table path) functions; GOCART-2G described in Collow et al. [2024], <https://doi.org/10.5194/gmd-17-1443-2024>
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `scale_factor` | 1.0 |  -  | settling velocity factor |
| `simple_scheme` | False |  -  | Read wet-particle radius and density from NetCDF Mie optics tables when `true`; otherwise compute wet swelling internally from species metadata (Gerber/Fitzgerald). Both paths are supported by the C++ core and are numerically identical to the legacy Fortran core (`ufs-community/CATChem` `develop` before the C++ rework). See [Optics tables](#optics-tables-simple_scheme-true) below. |
| `swelling_method` | 1 |  -  | method for calculating particle swelling: 1 Fitzgerald 1975; 2 for Gerber 1985 |
| `correction_maring` | False |  -  | correct the settling velocity following Maring et al, 2003 |

#### Required Meteorological Fields

- `T` - Temperature [K]
- `TSTEP` - Model time step [s]
- `AIRDEN` - Wet air density [kg/m3]
- `RH` - Relative humidity [fraction, not %]
- `Z` - Geometric height on the **interface** grid [m] (n_levels + 1)
- `PMID` - Layer mid-level pressure [Pa], the arithmetic mean of the edge pressures
- `DELP` - Layer pressure thickness [Pa]

> The `gocart` scheme delegates to the upstream GOCART2G `Chem_Settling` /
> `Chem_SettlingSimple` kernels so results are numerically identical to
> the legacy Fortran core. Missing or stale fields raise an explicit error before the
> kernel is invoked.

### Optics tables (`simple_scheme: true`)

The optics-table path reproduces the legacy GOCART2G wet-particle settling by
reading a per-aerosol-type **Mie lookup table** instead of computing swelling in
code. Configuration lives in a top-level `mie:` section, and each aerosol species
selects its table through its `__mie_name` attribute:

The tables are the GOCART2G monochromatic Mie lookup files (`optics_*.nc`).
They are not part of the CATChem source tree (the full set is about 1.4 GB).
Where to get them:

- **UFS Weather Model regression tests** stage them with the GOCART input data:
  `${INPUTDATA_ROOT}/GOCART/p8/ExtData` is linked into the run directory
  (`tests/fv3_conf/control_run.IN`), so `./ExtData/monochromatic/` resolves there.
- **Other hosts** need the same `ExtData/monochromatic/` directory from a GOCART
  or GEOS data distribution (the UFS `INPUTDATA_ROOT` tarballs or a GEOS-ESM
  GOCART checkout's `ExtData`), copied or linked into the run directory.
- **CATChem's own tests** look for them under
  `CATCHEM_SETTLING_OPTICS_FIXTURE_DIR` (`tests/CMakeLists.txt`, default
  `specs/tmp`); the optics-table tests are registered only when the files are present.

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

The process is a C++ class, `catchem::SettlingProcess` (`src/process/settling/catchem_process_settling.cpp`), registered with the `ProcessRegistry` at start-up; a host model does not instantiate it directly. It is enabled and configured in the runtime YAML and runs when the core advances the configured `run_phases`:

```yaml
processes:
  settling:
    activate: true
    scheme: 'gocart'
```

Hosts drive the core through the Fortran proxy `CATChem_Model` (`src/api/CATChem_API.F90`) or `catchem::Core` in C++; see the [Process Interface API](../../api/process-interface.md). The science kernel itself stays in Fortran and is called through `SettlingScienceBridge.F90`.

### Scheme Selection

Available schemes: `gocart`. Select one in the YAML configuration; scheme-specific parameters live under a key with the scheme's name (see the parameter tables above):

```yaml
processes:
  settling:
    scheme: 'gocart'
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

## Files

### Source Code
- `src/process/settling/SettlingCommon_Mod.F90` - Scheme options, constants, and shared types used by the schemes
- `src/process/settling/SettlingScienceBridge.F90` - `bind(C)` bridge that unpacks the C++ state and calls the Fortran scheme
- `src/process/settling/catchem_process_settling.cpp` - C++ process class `catchem::SettlingProcess`: configuration, species resolution, diagnostics registration, and the per-step call into the bridge
- `src/process/settling/catchem_process_settling.hpp` - C++ process class `catchem::SettlingProcess` (declaration)
- `src/process/settling/schemes/SettlingScheme_GOCART_Mod.F90` - GOCART gravitational settling scheme

### Tests
- `tests/run_settling_optics_parity.py`
- `tests/test_catchem_settling.cpp`
- `tests/test_catchem_settling_optics.cpp`

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

- GOCART: GOCART2G process library Chem_SettlingSimple function; GOCART-2G described in Collow, A. B., et al. [2024], Geosci. Model Dev., 17, 1443–1468, <https://doi.org/10.5194/gmd-17-1443-2024>

---
*This documentation was automatically generated by the CATChem Process Generator on 2025-12-18T14:12:33.301923*
