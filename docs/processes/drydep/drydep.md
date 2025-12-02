# DryDep Process

**Process Type:** Deposition
**Description:** Process for computing dry deposition of gas and aerosol species
**Author:** Wei Li
**Generated:** 2025-11-25T22:20:02.547522

## Overview

The DryDep process implements Process for computing dry deposition of gas and aerosol species. This process provides a modular, extensible framework for deposition calculations within the CATChem chemical transport model.

## Available Schemes

### WESELY Scheme

**Name:** `wesely`
**Description:** Wesely 1989 gas dry deposition scheme
**Author:** Wei Li
**Reference:** Wesely, M. L. [1989] Parameterization of surface resistances to gaseous dry deposition...
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `scale_factor` | 1.0 |  -  | DryDep velocity scale factor |
| `co2_effect` | True |  -  | Apply CO2 effect on stomatal conductance |
| `co2_level` | 600.0 |  -  | Ambient CO2 level for stomatal conductance adjustment |
| `co2_reference` | 380.0 |  -  | Reference CO2 level for stomatal conductance adjustment |

#### Required Meteorological Fields

- `USTAR` - Meteorological field required for scheme computation
- `TSTEP` - Meteorological field required for scheme computation
- `TS` - Meteorological field required for scheme computation
- `SWGDN` - Meteorological field required for scheme computation
- `SUNCOSmid` - Meteorological field required for scheme computation
- `OBK` - Meteorological field required for scheme computation
- `CLDFRC` - Meteorological field required for scheme computation
- `BXHEIGHT` - Meteorological field required for scheme computation
- `Z0` - Meteorological field required for scheme computation
- `PS` - Meteorological field required for scheme computation
- `FRLAI` - Meteorological field required for scheme computation
- `ILAND` - Meteorological field required for scheme computation
- `SALINITY` - Meteorological field required for scheme computation
- `FRLANDUSE` - Meteorological field required for scheme computation
- `TSKIN` - Meteorological field required for scheme computation
- `LON` - Meteorological field required for scheme computation
- `LAT` - Meteorological field required for scheme computation
- `LUCNAME` - Meteorological field required for scheme computation
- `IsSnow` - Meteorological field required for scheme computation
- `IsIce` - Meteorological field required for scheme computation
- `IsLand` - Meteorological field required for scheme computation


### GOCART Scheme

**Name:** `gocart`
**Description:** GOCART-2G aerosol dry deposition scheme
**Author:** Wei Li & Lacey Holland
**Reference:** Allison et al. [2024] Benchmarking GOCART-2G in GEOS
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `scale_factor` | 1.0 |  -  | Dry deposition velocity scale factor |
| `resuspension` | False |  -  | Apply resuspension for dry deposition |

#### Required Meteorological Fields

- `USTAR` - Meteorological field required for scheme computation
- `TSTEP` - Meteorological field required for scheme computation
- `T` - Meteorological field required for scheme computation
- `AIRDEN` - Meteorological field required for scheme computation
- `Z` - Meteorological field required for scheme computation
- `LWI` - Meteorological field required for scheme computation
- `PBLH` - Meteorological field required for scheme computation
- `HFLUX` - Meteorological field required for scheme computation
- `Z0H` - Meteorological field required for scheme computation
- `U10M` - Meteorological field required for scheme computation
- `V10M` - Meteorological field required for scheme computation
- `FRLAKE` - Meteorological field required for scheme computation
- `GWETTOP` - Meteorological field required for scheme computation


### ZHANG Scheme

**Name:** `zhang`
**Description:** Zhang et al. [2001] scheme with Emerson et al. [2020] updates
**Author:** Wei Li
**Reference:** Zhang et al., 2001; Emerson et al., 2020
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `scale_factor` | 1.0 |  -  | Dry deposition velocity scale factor |

#### Required Meteorological Fields

- `USTAR` - Meteorological field required for scheme computation
- `TSTEP` - Meteorological field required for scheme computation
- `TS` - Meteorological field required for scheme computation
- `OBK` - Meteorological field required for scheme computation
- `BXHEIGHT` - Meteorological field required for scheme computation
- `Z0` - Meteorological field required for scheme computation
- `RH` - Meteorological field required for scheme computation
- `PS` - Meteorological field required for scheme computation
- `U10M` - Meteorological field required for scheme computation
- `V10M` - Meteorological field required for scheme computation
- `FRLANDUSE` - Meteorological field required for scheme computation
- `ILAND` - Meteorological field required for scheme computation
- `LUCNAME` - Meteorological field required for scheme computation
- `IsSnow` - Meteorological field required for scheme computation
- `IsIce` - Meteorological field required for scheme computation



## Process Interface

### Species

The drydep process operates on the following chemical species:


### Required Inputs



### Process Diagnostics

| Diagnostic | Units | Description |
|------------|-------|-------------|
| `drydep_con_per_species` | ug/kg or ppm | Dry deposition concentration per species |
| `drydep_velocity_per_species` | m/s | Dry deposition velocity |

## Usage

### Basic Integration

```fortran
use DryDepProcessCreator_Mod
use DryDepCommon_Mod

! Create process instance
type(DryDepProcess_t) :: process
call create_drydep_process(process, config_data)

! Use process in model time step
call process%run(state, dt)
```

### Scheme Selection

The process supports multiple schemes. Select your desired scheme:

```fortran
! Use WESELY scheme
process%scheme_name = "wesely"
```
```fortran
! Use GOCART scheme
process%scheme_name = "gocart"
```
```fortran
! Use ZHANG scheme
process%scheme_name = "zhang"
```

## Implementation Details

### Pure Science Kernels

Each scheme is implemented as a pure science kernel with no infrastructure dependencies:

```fortran
! WESELY scheme
pure subroutine compute_wesely( &
   num_layers, num_species, params, &
   USTAR, &   TSTEP, &   TS, &   SWGDN, &   SUNCOSmid, &   OBK, &   CLDFRC, &   BXHEIGHT, &   Z0, &   PS, &   FRLAI, &   ILAND, &   SALINITY, &   FRLANDUSE, &   TSKIN, &   LON, &   LAT, &   LUCNAME, &   IsSnow, &   IsIce, &   IsLand, &
   species_conc, emission_flux)
```
```fortran
! GOCART scheme
pure subroutine compute_gocart( &
   num_layers, num_species, params, &
   USTAR, &   TSTEP, &   T, &   AIRDEN, &   Z, &   LWI, &   PBLH, &   HFLUX, &   Z0H, &   U10M, &   V10M, &   FRLAKE, &   GWETTOP, &
   species_conc, emission_flux)
```
```fortran
! ZHANG scheme
pure subroutine compute_zhang( &
   num_layers, num_species, params, &
   USTAR, &   TSTEP, &   TS, &   OBK, &   BXHEIGHT, &   Z0, &   RH, &   PS, &   U10M, &   V10M, &   FRLANDUSE, &   ILAND, &   LUCNAME, &   IsSnow, &   IsIce, &
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
  drydep:
    enabled: true
    scheme: "wesely"
    parameters:
      scale_factor: 1.0
      co2_effect: True
      co2_level: 600.0
      co2_reference: 380.0
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
- `src/process/drydep/ProcessDryDepInterface_Mod.F90` - Main process interface
- `src/process/drydep/DryDepCommon_Mod.F90` - Common types and parameters
- `src/process/drydep/DryDepProcessCreator_Mod.F90` - Process factory
- `src/process/drydep/schemes/DryDepScheme_WESELY_Mod.F90` - Wesely 1989 gas dry deposition scheme
- `src/process/drydep/schemes/DryDepScheme_GOCART_Mod.F90` - GOCART-2G aerosol dry deposition scheme
- `src/process/drydep/schemes/DryDepScheme_ZHANG_Mod.F90` - Zhang et al. [2001] scheme with Emerson et al. [2020] updates

### Tests
- `tests/process/drydep/unit/` - Unit tests
- `tests/process/drydep/integration/` - Integration tests

### Documentation
- `docs/processes/drydep/drydep.md` - This documentation

## Contributing

When modifying or extending this process:

1. **Science Changes:** Modify the scheme modules in `schemes/`
2. **Interface Changes:** Update the main interface module
3. **New Schemes:** Add new scheme modules and update the creator
4. **Tests:** Add corresponding unit and integration tests
5. **Documentation:** Update this documentation file

## References

- WESELY: Wesely, M. L. [1989] Parameterization of surface resistances to gaseous dry deposition...
- GOCART: Allison et al. [2024] Benchmarking GOCART-2G in GEOS
- ZHANG: Zhang et al., 2001; Emerson et al., 2020

---
*This documentation was automatically generated by the CATChem Process Generator on 2025-11-25T22:20:02.547522*
