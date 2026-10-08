# Dry Deposition Process

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
**Reference:** Wesely, M. L. [1989] Parameterization of surface resistances to gaseous dry deposition in regional-scale numerical models, Atmos. Environ., 23(6), 1293–1304, <https://doi.org/10.1016/0004-6981(89)90153-4>
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `scale_factor` | 1.0 |  -  | DryDep velocity scale factor |
| `co2_effect` | True |  -  | Apply CO2 effect on stomatal conductance |
| `co2_level` | 600.0 |  -  | Ambient CO2 level for stomatal conductance adjustment |
| `co2_reference` | 380.0 |  -  | Reference CO2 level for stomatal conductance adjustment |

#### Required Meteorological Fields

- `USTAR` - Friction velocity [m/s]
- `TSTEP` - Model time step [s]
- `TS` - Surface temperature [K]
- `SWGDN` - Incident radiation @ ground [W/m2]
- `SUNCOSmid` - COS(solar zenith angle) at midpoint of chem timestep
- `OBK` - Monin-Obhukov length [m]
- `CLDFRC` - Column cloud fraction [1]
- `BXHEIGHT` - Grid box height [m] (dry air)
- `Z0` - Surface roughness height [m]
- `PS` - Surface Pressure [Pa]
- `FRLAI` - LAI in each Fractional Land use type [m2/m2] (nx,ny,nlanduse)
- `ILAND` - Land type ID in current grid box (nx,ny,nlanduse)
- `SALINITY` - Salinity of the ocean [part per thousand]
- `FRLANDUSE` - Fractional Land Use (nx,ny,nlanduse)
- `TSKIN` - Surface skin temperature [K]
- `LON` - Longitude
- `LAT` - Latitude
- `LUCNAME` - name of land use category
- `IsSnow` - Is this a snow grid box?
- `IsIce` - Is this an ice grid box?
- `IsLand` - Is this a land grid box?


### GOCART Scheme

**Name:** `gocart`
**Description:** GOCART-2G aerosol dry deposition scheme
**Author:** Wei Li & Lacey Holland
**Reference:** Collow, A. B., et al. [2024] Benchmarking GOCART-2G in the Goddard Earth Observing System (GEOS), Geosci. Model Dev., 17, 1443–1468, https://doi.org/10.5194/gmd-17-1443-2024
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `scale_factor` | 1.0 |  -  | Dry deposition velocity scale factor |
| `resuspension` | False |  -  | Apply resuspension for dry deposition |

#### Required Meteorological Fields

- `USTAR` - Friction velocity [m/s]
- `TSTEP` - Model time step [s]
- `T` - Temperature [K]
- `AIRDEN` - Wet air density [kg/m3]
- `Z` - Geopotential Height @ level edges [m] (nx,ny,nz+1)
- `LWI` - Land water ice mask (0-sea, 1-land, 2-ice)
- `PBLH` - PBL height [m]
- `HFLUX` - Sensible heat flux [W/m2]
- `Z0H` - Surface roughness height, for heat (thermal roughness) [m]
- `U10M` - E/W wind speed @ 10m ht [m/s]
- `V10M` - N/S wind speed @ 10m ht [m/s]
- `FRLAKE` - Fraction of lake [1]
- `GWETTOP` - Top soil moisture [1]


### ZHANG Scheme

**Name:** `zhang`
**Description:** Zhang et al. [2001] scheme with Emerson et al. [2020] updates
**Author:** Wei Li
**Reference:** Zhang et al. [2001], Atmos. Environ., 35(3), 549–560, <https://doi.org/10.1016/S1352-2310(00)00326-5>; Emerson et al. [2020], PNAS, 117(42), 26076–26082, <https://doi.org/10.1073/pnas.2014761117>
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `scale_factor` | 1.0 |  -  | Dry deposition velocity scale factor |

#### Required Meteorological Fields

- `USTAR` - Friction velocity [m/s]
- `TSTEP` - Model time step [s]
- `TS` - Surface temperature [K]
- `OBK` - Monin-Obhukov length [m]
- `BXHEIGHT` - Grid box height [m] (dry air)
- `Z0` - Surface roughness height [m]
- `RH` - Relative humidity [fraction, not %]
- `PS` - Surface Pressure [Pa]
- `U10M` - E/W wind speed @ 10m ht [m/s]
- `V10M` - N/S wind speed @ 10m ht [m/s]
- `FRLANDUSE` - Fractional Land Use (nx,ny,nlanduse)
- `ILAND` - Land type ID in current grid box (nx,ny,nlanduse)
- `LUCNAME` - name of land use category
- `IsSnow` - Is this a snow grid box?
- `IsIce` - Is this an ice grid box?



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

- WESELY: Wesely, M. L. [1989] Parameterization of surface resistances to gaseous dry deposition in regional-scale numerical models, Atmos. Environ., 23(6), 1293–1304, <https://doi.org/10.1016/0004-6981(89)90153-4>
- GOCART: Collow, A. B., et al. [2024] Benchmarking GOCART-2G in the Goddard Earth Observing System (GEOS), Geosci. Model Dev., 17, 1443–1468, <https://doi.org/10.5194/gmd-17-1443-2024>
- ZHANG: Zhang, L., et al. [2001] A size-segregated particle dry deposition scheme, Atmos. Environ., 35(3), 549–560, <https://doi.org/10.1016/S1352-2310(00)00326-5>; Emerson, E. W., et al. [2020] Revisiting particle dry deposition, PNAS, 117(42), 26076–26082, <https://doi.org/10.1073/pnas.2014761117>

---
*This documentation was automatically generated by the CATChem Process Generator on 2025-11-25T22:20:02.547522*
