# Wet Deposition Process

**Process Type:** Deposition
**Description:** Process for computing wet deposition of gas and aerosol species
**Author:** Wei Li
**Generated:** 2026-09-22T13:06:55.756402

## Overview

The WetDep process implements Process for computing wet deposition of gas and aerosol species. This process provides a modular, extensible framework for deposition calculations within the CATChem chemical transport model.

## Available Schemes

### JACOB Scheme

**Name:** `jacob`
**Description:** Jacob et al. [2000] wet deposition scheme
**Author:** Wei Li
**Reference:** Jacob, D. J. et al., [2000] Harvard wet deposition scheme for GMI; peer-reviewed description in Liu et al. [2001], J. Geophys. Res., 106(D11), 12109–12128, <https://doi.org/10.1029/2000JD900839>
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `scale_factor` | 1.0 |  -  | Washout tuning factor |
| `radius_threshold` | 1.0 |  -  | Radius threshold for aerosol wet deposition (um) |
| `so4_gocart_resusp` | True |  -  | Sulfate-only GOCART-style resuspension toggle (default on) |
| `so4_washout_eff` | 1.0 |  -  | Sulfate-only below-cloud washout efficiency (1.0=unchanged); reduce to cut excess SO4 washout |

#### Required Meteorological Fields

- `T` - Temperature [K]
- `TSTEP` - Model time step [s]
- `AIRDEN_DRY` - Dry air density [kg/m3]
- `MAIRDEN` - Moist air density (same as AIRDEN to cover possible use cases) [kg/m3]
- `PFLLSAN` - Dwn flux liq prec:LS+anv [kg/m2/s] (nx,ny,nz+1)
- `PFILSAN` - Dwn flux ice prec:LS+anv [kg/m2/s] (nx,ny,nz+1)
- `PEDGE` - Air partial pressure @ level edges [Pa] (nx,ny,nz+1)
- `REEVAPLS` - Evap of precip LS+anvil [kg/kg/s] (assume per dry air)


### GOCART Scheme

**Name:** `gocart`
**Description:** GOCART2G wet removal scheme: SU_Wet_Removal for sulfate species (DMS/SO2/SO4/MSA) and WetRemovalUFS for all other species
**Author:** Wei Li
**Reference:** GOCART2G Process Library: SU_Wet_Removal (Chin et al. [2000], J. Geophys. Res., 105(D20), 24671–24687, <https://doi.org/10.1029/2000JD900384>) and WetRemovalUFS (Liu et al. [2001], J. Geophys. Res., 106(D11), 12109–12128, <https://doi.org/10.1029/2000JD900839>)
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `scale_factor` | 1.0 |  -  | Overall washout tuning factor |
| `washout_tuning` | 1.0 |  -  | WetRemovalUFS below-cloud washout tuning factor (wtune) |
| `radius_threshold` | 1.0 |  -  | Radius threshold for aerosol washout (um) (WetRemovalUFS radius_thr) |

#### Required Meteorological Fields

- `T` - Temperature [K]
- `TSTEP` - Model time step [s]
- `MAIRDEN` - Moist air density (same as AIRDEN to cover possible use cases) [kg/m3]
- `PEDGE` - Air partial pressure @ level edges [Pa] (nx,ny,nz+1)
- `PFLLSAN` - Dwn flux liq prec:LS+anv [kg/m2/s] (nx,ny,nz+1)
- `PFILSAN` - Dwn flux ice prec:LS+anv [kg/m2/s] (nx,ny,nz+1)
- `PRECCON` - Conv  precip @ ground [kg/m2/s] -> [mm/day]
- `PRECLSC` - Large-scale precip @ ground kg/m2/s] -> [mm/day]



## Process Interface

### Species

The wetdep process operates on the following chemical species:


### Required Inputs



### Process Diagnostics

| Diagnostic | Units | Description |
|------------|-------|-------------|
| `wetdep_mass_per_species_per_level` | kg/m2 | Wet deposition mass loss per species per level |
| `wetdep_flux_per_species_per_level` | kg/m2/s | Wet deposition flux per species per level |

## Usage

### Basic Integration

```fortran
use WetDepProcessCreator_Mod
use WetDepCommon_Mod

! Create process instance
type(WetDepProcess_t) :: process
call create_wetdep_process(process, config_data)

! Use process in model time step
call process%run(state, dt)
```

### Scheme Selection

The process supports multiple schemes. Select your desired scheme:

```fortran
! Use JACOB scheme
process%scheme_name = "jacob"
```
```fortran
! Use GOCART scheme
process%scheme_name = "gocart"
```

## Implementation Details

### Pure Science Kernels

Each scheme is implemented as a pure science kernel with no infrastructure dependencies:

```fortran
! JACOB scheme
pure subroutine compute_jacob( &
   num_layers, num_species, params, &
   T, &   TSTEP, &   AIRDEN_DRY, &   MAIRDEN, &   PFLLSAN, &   PFILSAN, &   PEDGE, &   REEVAPLS, &
   species_conc, emission_flux)
```
```fortran
! GOCART scheme
pure subroutine compute_gocart( &
   num_layers, num_species, params, &
   T, &   TSTEP, &   MAIRDEN, &   PEDGE, &   PFLLSAN, &   PFILSAN, &   PRECCON, &   PRECLSC, &
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
  wetdep:
    enabled: true
    scheme: "jacob"
    parameters:
      scale_factor: 1.0
      radius_threshold: 1.0
      so4_gocart_resusp: True
      so4_washout_eff: 1.0
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
- `src/process/wetdep/ProcessWetDepInterface_Mod.F90` - Main process interface
- `src/process/wetdep/WetDepCommon_Mod.F90` - Common types and parameters
- `src/process/wetdep/WetDepProcessCreator_Mod.F90` - Process factory
- `src/process/wetdep/schemes/WetDepScheme_JACOB_Mod.F90` - Jacob et al. [2000] wet deposition scheme
- `src/process/wetdep/schemes/WetDepScheme_GOCART_Mod.F90` - GOCART2G wet removal scheme: SU_Wet_Removal for sulfate species (DMS/SO2/SO4/MSA) and WetRemovalUFS for all other species

### Tests
- `tests/process/wetdep/unit/` - Unit tests
- `tests/process/wetdep/integration/` - Integration tests

### Documentation
- `docs/processes/wetdep/wetdep.md` - This documentation

## Contributing

When modifying or extending this process:

1. **Science Changes:** Modify the scheme modules in `schemes/`
2. **Interface Changes:** Update the main interface module
3. **New Schemes:** Add new scheme modules and update the creator
4. **Tests:** Add corresponding unit and integration tests
5. **Documentation:** Update this documentation file

## References

- JACOB: Jacob, D. J. et al., [2000] Harvard wet deposition scheme for GMI; peer-reviewed description in Liu, H., Jacob, D. J., Bey, I., & Yantosca, R. M. [2001], J. Geophys. Res., 106(D11), 12109–12128, <https://doi.org/10.1029/2000JD900839>

---
*This documentation was automatically generated by the CATChem Process Generator on 2026-09-22T13:06:55.756402*
