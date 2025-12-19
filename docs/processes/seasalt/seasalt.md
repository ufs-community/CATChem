# SeaSalt Process

**Process Type:** Emission
**Description:** Process for computing sea salt aerosol emissions over ocean surfaces
**Author:** Barry Baker & Wei Li
**Generated:** 2025-11-14T23:01:21.952377

## Overview

The SeaSalt process implements Process for computing sea salt aerosol emissions over ocean surfaces. This process provides a modular, extensible framework for emission calculations within the CATChem chemical transport model.

## Available Schemes

### GONG97 Scheme

**Name:** `gong97`
**Description:** Gong 1997 sea salt emission scheme
**Author:** Barry Baker
**Reference:** Gong et al. [1997]
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `scale_factor` | 1.0 | 0.1 - 10.0 | Emission scale factor |
| `weibull_flag` | False |  -  | Apply Weibull distribution for particle size |

#### Required Meteorological Fields

- `FROCEAN` - Meteorological field required for scheme computation
- `FRSEAICE` - Meteorological field required for scheme computation
- `SST` - Meteorological field required for scheme computation
- `U10M` - Meteorological field required for scheme computation
- `V10M` - Meteorological field required for scheme computation


### GONG03 Scheme

**Name:** `gong03`
**Description:** Gong 2003 sea salt emission scheme with improved sub- and super-micron treatment
**Author:** Barry Baker
**Reference:** Gong [2003]
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `scale_factor` | 1.0 | 0.1 - 10.0 | Emission scale factor |
| `weibull_flag` | False |  -  | Apply Weibull distribution for particle size |

#### Required Meteorological Fields

- `FROCEAN` - Meteorological field required for scheme computation
- `FRSEAICE` - Meteorological field required for scheme computation
- `SST` - Meteorological field required for scheme computation
- `U10M` - Meteorological field required for scheme computation
- `V10M` - Meteorological field required for scheme computation


### GEOS12 Scheme

**Name:** `geos12`
**Description:** GEOS-Chem 2012 sea salt emission scheme with observational constraints
**Author:** Barry Baker
**Reference:** Jaeglé et al. [2011]
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `scale_factor` | 1.0 | 0.1 - 10.0 | Emission scale factor |

#### Required Meteorological Fields

- `FROCEAN` - Meteorological field required for scheme computation
- `FRSEAICE` - Meteorological field required for scheme computation
- `SST` - Meteorological field required for scheme computation
- `USTAR` - Meteorological field required for scheme computation



## Process Interface

### Species

The seasalt process operates on the following chemical species:


### Required Inputs

#### Meteorological Fields
- `DELP` - Required meteorological input


### Process Diagnostics

| Diagnostic | Units | Description |
|------------|-------|-------------|
| `seasalt_mass_emission_total` | kg/m2/s | Sea salt mass emission flux total |
| `seasalt_number_emission_total` | kg/m2/s | Sea salt number emission flux total |

## Usage

### Basic Integration

```fortran
use SeaSaltProcessCreator_Mod
use SeaSaltCommon_Mod

! Create process instance
type(SeaSaltProcess_t) :: process
call create_seasalt_process(process, config_data)

! Use process in model time step
call process%run(state, dt)
```

### Scheme Selection

The process supports multiple schemes. Select your desired scheme:

```fortran
! Use GONG97 scheme
process%scheme_name = "gong97"
```
```fortran
! Use GONG03 scheme
process%scheme_name = "gong03"
```
```fortran
! Use GEOS12 scheme
process%scheme_name = "geos12"
```

## Implementation Details

### Pure Science Kernels

Each scheme is implemented as a pure science kernel with no infrastructure dependencies:

```fortran
! GONG97 scheme
pure subroutine compute_gong97( &
   num_layers, num_species, params, &
   FROCEAN, &   FRSEAICE, &   SST, &   U10M, &   V10M, &
   species_conc, emission_flux)
```
```fortran
! GONG03 scheme
pure subroutine compute_gong03( &
   num_layers, num_species, params, &
   FROCEAN, &   FRSEAICE, &   SST, &   U10M, &   V10M, &
   species_conc, emission_flux)
```
```fortran
! GEOS12 scheme
pure subroutine compute_geos12( &
   num_layers, num_species, params, &
   FROCEAN, &   FRSEAICE, &   SST, &   USTAR, &
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
  seasalt:
    enabled: true
    scheme: "gong97"
    parameters:
      scale_factor: 1.0
      weibull_flag: False
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
- `src/process/seasalt/ProcessSeaSaltInterface_Mod.F90` - Main process interface
- `src/process/seasalt/SeaSaltCommon_Mod.F90` - Common types and parameters
- `src/process/seasalt/SeaSaltProcessCreator_Mod.F90` - Process factory
- `src/process/seasalt/schemes/SeaSaltScheme_GONG97_Mod.F90` - Gong 1997 sea salt emission scheme
- `src/process/seasalt/schemes/SeaSaltScheme_GONG03_Mod.F90` - Gong 2003 sea salt emission scheme with improved sub- and super-micron treatment
- `src/process/seasalt/schemes/SeaSaltScheme_GEOS12_Mod.F90` - GEOS-Chem 2012 sea salt emission scheme with observational constraints

### Tests
- `tests/process/seasalt/unit/` - Unit tests
- `tests/process/seasalt/integration/` - Integration tests

### Documentation
- `docs/processes/seasalt/seasalt.md` - This documentation

## Contributing

When modifying or extending this process:

1. **Science Changes:** Modify the scheme modules in `schemes/`
2. **Interface Changes:** Update the main interface module
3. **New Schemes:** Add new scheme modules and update the creator
4. **Tests:** Add corresponding unit and integration tests
5. **Documentation:** Update this documentation file

## References

- GONG97: Gong et al. [1997]
- GONG03: Gong [2003]
- GEOS12: Jaeglé et al. [2011]

---
*This documentation was automatically generated by the CATChem Process Generator on 2025-11-14T23:01:21.952377*
