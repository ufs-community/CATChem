# Dust Process

**Process Type:** Emission
**Description:** Process for computing windblown dust emissions
**Author:** Wei Li & Barry Baker
**Generated:** 2026-04-17T13:57:10.461600

## Overview

The Dust process implements Process for computing windblown dust emissions. This process provides a modular, extensible framework for emission calculations within the CATChem chemical transport model.

## Available Schemes

### FENGSHA Scheme

**Name:** `fengsha`
**Description:** Fengsha Dust emission scheme developed at NOAA ARL for use at NOAA NWS
**Author:** Barry Baker & Wei Li
**Reference:** Zhang et al. 2022
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `alpha` | 0.2 | 0.001 - 100.0 | linear scaling factor |
| `gamma` | 1.0 | 0.0 - 10.0 | Exponential scaling factor on source parameter |
| `drylimit_factor` | 1.0 | 0 - 10 | Dry Limit factor modifying the Fecan dry limit following Zender 2003 |
| `moist_correction_factor` | 1.0 | 0 - 10 | Moisture correction factor |
| `kvhmax` | 0.0002 | 0.0 - 1.0 | Maximum vertical to horizontal flux ratio |
| `drag_option` | 1 | 1 - 4 | Drag Partition Option: 1 - use input drag, 2 - Darmenova, 3 - Leung 2022, 4 - MB95 |
| `horizflux_option` | 1 | 1 - 3 | Horizontal flux option: 1 - White (1979), 2 - Draxler (2001), 3 - Kawamura (1964) |
| `moist_option` | 1 | 1 - 2 | Moisture parameterization: 1 - Fecan, 2 - Zhao (not implemented yet) |
| `distribution_option` | 1 | 1 - 1 | Dust Distribution option: 1 - Kok 2011, 2 - Meng 2022 (not implemented yet) |

#### Required Meteorological Fields

- `USTAR` - Meteorological field required for scheme computation
- `TSKIN` - Meteorological field required for scheme computation
- `AIRDEN` - Meteorological field required for scheme computation
- `SOILM` - Meteorological field required for scheme computation
- `Z0` - Meteorological field required for scheme computation
- `GVF` - Meteorological field required for scheme computation
- `LWI` - Meteorological field required for scheme computation
- `LAI` - Meteorological field required for scheme computation
- `FRLAKE` - Meteorological field required for scheme computation
- `FRSNO` - Meteorological field required for scheme computation
- `CLAYFRAC` - Meteorological field required for scheme computation
- `SANDFRAC` - Meteorological field required for scheme computation
- `RDRAG` - Meteorological field required for scheme computation
- `SSM` - Meteorological field required for scheme computation
- `USTAR_THRESHOLD` - Meteorological field required for scheme computation


### GINOUX Scheme

**Name:** `ginoux`
**Description:** Ginoux dust emission scheme
**Author:** Barry Baker & Wei Li
**Reference:** Ginoux et al. [2001]
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `Ch_DU` | [1.0, 1.0, 1.0, 1.0, 1.0] | 0.001 - 100.0 | Dust tuning coefficient per species bin |

#### Required Meteorological Fields

- `FRLAKE` - Meteorological field required for scheme computation
- `FRSNO` - Meteorological field required for scheme computation
- `TSKIN` - Meteorological field required for scheme computation
- `AIRDEN` - Meteorological field required for scheme computation
- `LWI` - Meteorological field required for scheme computation
- `GWETTOP` - Meteorological field required for scheme computation
- `U10M` - Meteorological field required for scheme computation
- `V10M` - Meteorological field required for scheme computation
- `SSM` - Meteorological field required for scheme computation



## Process Interface

### Species

The dust process operates on the following chemical species:


### Required Inputs

#### Meteorological Fields
- `DELP` - Required meteorological input


### Process Diagnostics

| Diagnostic | Units | Description |
|------------|-------|-------------|
| `dust_emission_total` | kg/m2/s | Total dust emissions for all bins |
| `dust_emission_per_bin` | kg/m2/s | Dust emission flux per bin |

## Usage

### Basic Integration

```fortran
use DustProcessCreator_Mod
use DustCommon_Mod

! Create process instance
type(DustProcess_t) :: process
call create_dust_process(process, config_data)

! Use process in model time step
call process%run(state, dt)
```

### Scheme Selection

The process supports multiple schemes. Select your desired scheme:

```fortran
! Use FENGSHA scheme
process%scheme_name = "fengsha"
```
```fortran
! Use GINOUX scheme
process%scheme_name = "ginoux"
```

## Implementation Details

### Pure Science Kernels

Each scheme is implemented as a pure science kernel with no infrastructure dependencies:

```fortran
! FENGSHA scheme
pure subroutine compute_fengsha( &
   num_layers, num_species, params, &
   USTAR, &   TSKIN, &   AIRDEN, &   SOILM, &   Z0, &   GVF, &   LWI, &   LAI, &   FRLAKE, &   FRSNO, &   CLAYFRAC, &   SANDFRAC, &   RDRAG, &   SSM, &   USTAR_THRESHOLD, &
   species_conc, emission_flux)
```
```fortran
! GINOUX scheme
pure subroutine compute_ginoux( &
   num_layers, num_species, params, &
   FRLAKE, &   FRSNO, &   TSKIN, &   AIRDEN, &   LWI, &   GWETTOP, &   U10M, &   V10M, &   SSM, &
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
  dust:
    enabled: true
    scheme: "fengsha"
    parameters:
      alpha: 0.2
      gamma: 1.0
      drylimit_factor: 1.0
      moist_correction_factor: 1.0
      kvhmax: 0.0002
      drag_option: 1
      horizflux_option: 1
      moist_option: 1
      distribution_option: 1
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
- `src/process/dust/ProcessDustInterface_Mod.F90` - Main process interface
- `src/process/dust/DustCommon_Mod.F90` - Common types and parameters
- `src/process/dust/DustProcessCreator_Mod.F90` - Process factory
- `src/process/dust/schemes/DustScheme_FENGSHA_Mod.F90` - Fengsha Dust emission scheme developed at NOAA ARL for use at NOAA NWS
- `src/process/dust/schemes/DustScheme_GINOUX_Mod.F90` - Ginoux dust emission scheme

### Tests
- `tests/process/dust/unit/` - Unit tests
- `tests/process/dust/integration/` - Integration tests

### Documentation
- `docs/processes/dust/dust.md` - This documentation

## Contributing

When modifying or extending this process:

1. **Science Changes:** Modify the scheme modules in `schemes/`
2. **Interface Changes:** Update the main interface module
3. **New Schemes:** Add new scheme modules and update the creator
4. **Tests:** Add corresponding unit and integration tests
5. **Documentation:** Update this documentation file

## References

- FENGSHA: Zhang et al. 2022
- GINOUX: Ginoux et al. [2001]

---
*This documentation was automatically generated by the CATChem Process Generator on 2026-04-17T13:57:10.461600*
