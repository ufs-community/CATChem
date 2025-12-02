# ExternalEmissionData Process

**Process Type:** Emission
**Description:** External emission data management process for CATChem. Handles loading, validation, interpolation, and access to external emission inventories from files. This process provides a standardized interface for external emission data that can be used by emission processes and the driver.
**Author:** CATChem Development Team
**Generated:** 2025-07-07T15:37:42.346554

## Overview

The ExternalEmissionData process implements External emission data management process for CATChem. Handles loading, validation, interpolation, and access to external emission inventories from files. This process provides a standardized interface for external emission data that can be used by emission processes and the driver.. This process provides a modular, extensible framework for emission calculations within the CATChem chemical transport model.

## Available Schemes

### Surface2D Scheme

**Name:** `surface_2d`
**Description:** Handler for 2D surface emission data from external inventories
**Author:** CATChem Development Team

#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `interpolation_method` |  |  -  |  |
| `time_interpolation` |  |  -  |  |
| `spatial_bounds_check` |  |  -  |  |
| `missing_value_handling` |  |  -  |  |


#### Available Diagnostics

| Diagnostic | Units | Description |
|------------|-------|-------------|
| `data_load_status` | dimensionless | Status of data loading (0=not loaded, 1=loaded, 2=error) |
| `memory_usage` | bytes | Memory usage for emission data storage |
| `field_count` | count | Number of emission fields loaded |

### Elevated3D Scheme

**Name:** `elevated_3d`
**Description:** Handler for 3D elevated emission data with vertical distribution
**Author:** CATChem Development Team

#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `interpolation_method` |  |  -  |  |
| `time_interpolation` |  |  -  |  |
| `vertical_interpolation` |  |  -  |  |
| `plume_rise_coupling` |  |  -  |  |

#### Required Meteorological Fields

- `pressure` - Meteorological field required for scheme computation
- `temperature` - Meteorological field required for scheme computation

#### Available Diagnostics

| Diagnostic | Units | Description |
|------------|-------|-------------|
| `vertical_distribution` | kg/m3/s | Vertical emission distribution |
| `plume_height` | m | Effective emission height |

### TemporalProfiles Scheme

**Name:** `temporal_profiles`
**Description:** Handler for time-varying emission profiles and scaling factors
**Author:** CATChem Development Team

#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `profile_resolution` |  |  -  |  |
| `seasonal_variation` |  |  -  |  |
| `diurnal_variation` |  |  -  |  |
| `weekly_variation` |  |  -  |  |


#### Available Diagnostics

| Diagnostic | Units | Description |
|------------|-------|-------------|
| `temporal_scaling` | dimensionless | Current temporal scaling factor |


## Process Interface

### Species

The external_emission_data process operates on the following chemical species:


### Required Inputs



### Process Diagnostics

| Diagnostic | Units | Description |
|------------|-------|-------------|
| `total_memory_usage` | MB | Total memory usage for all emission data |
| `total_categories` | count | Number of emission categories loaded |
| `total_fields` | count | Total number of emission fields across all categories |
| `data_validation_status` | dimensionless | Overall data validation status |
| `file_load_errors` | count | Number of file loading errors encountered |

## Usage

### Basic Integration

```fortran
use ExternalEmissionDataProcessCreator_Mod
use ExternalEmissionDataCommon_Mod

! Create process instance
type(ExternalEmissionDataProcess_t) :: process
call create_external_emission_data_process(process, config_data)

! Use process in model time step
call process%run(state, dt)
```

### Scheme Selection

The process supports multiple schemes. Select your desired scheme:

```fortran
! Use Surface2D scheme
process%scheme_name = "surface_2d"
```
```fortran
! Use Elevated3D scheme
process%scheme_name = "elevated_3d"
```
```fortran
! Use TemporalProfiles scheme
process%scheme_name = "temporal_profiles"
```

## Implementation Details

### Pure Science Kernels

Each scheme is implemented as a pure science kernel with no infrastructure dependencies:

```fortran
! Surface2D scheme
pure subroutine compute_surface_2d( &
   num_layers, num_species, params, &
   species_conc, emission_flux)
```
```fortran
! Elevated3D scheme
pure subroutine compute_elevated_3d( &
   num_layers, num_species, params, &
   pressure, &   temperature, &
   species_conc, emission_flux)
```
```fortran
! TemporalProfiles scheme
pure subroutine compute_temporal_profiles( &
   num_layers, num_species, params, &
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
  external_emission_data:
    enabled: true
    scheme: "surface_2d"
    parameters:
      interpolation_method:
      time_interpolation:
      spatial_bounds_check:
      missing_value_handling:
    diagnostics:
      enabled: true
      output_frequency: "daily"
```

## Technical Specifications

- **Parallelization:** Column
- **Memory Requirements:** Medium
- **Timestep Dependency:** Independent
- **Multiphase Support:** No
- **Size Bin Support:** No
- **Vectorization:** Supported

## Files Generated

### Source Code
- `src/process/external_emission_data/ProcessExternalEmissionDataInterface_Mod.F90` - Main process interface
- `src/process/external_emission_data/ExternalEmissionDataCommon_Mod.F90` - Common types and parameters
- `src/process/external_emission_data/ExternalEmissionDataProcessCreator_Mod.F90` - Process factory
- `src/process/external_emission_data/schemes/ExternalEmissionDataScheme_Surface2D_Mod.F90` - Handler for 2D surface emission data from external inventories
- `src/process/external_emission_data/schemes/ExternalEmissionDataScheme_Elevated3D_Mod.F90` - Handler for 3D elevated emission data with vertical distribution
- `src/process/external_emission_data/schemes/ExternalEmissionDataScheme_TemporalProfiles_Mod.F90` - Handler for time-varying emission profiles and scaling factors

### Tests
- `tests/process/external_emission_data/unit/` - Unit tests
- `tests/process/external_emission_data/integration/` - Integration tests

### Documentation
- `docs/processes/external_emission_data/external_emission_data.md` - This documentation

## Contributing

When modifying or extending this process:

1. **Science Changes:** Modify the scheme modules in `schemes/`
2. **Interface Changes:** Update the main interface module
3. **New Schemes:** Add new scheme modules and update the creator
4. **Tests:** Add corresponding unit and integration tests
5. **Documentation:** Update this documentation file

## References


---
*This documentation was automatically generated by the CATChem Process Generator on 2025-07-07T15:37:42.346554*
