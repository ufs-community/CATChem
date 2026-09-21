# External Emission Data Process - Integration Guide

## Overview

The External Emission Data Process is a special type of process in CATChem that handles external emission data management rather than computing emissions directly. It provides a standardized interface for loading, validating, interpolating, and accessing external emission inventories from files.

## Key Differences from Regular Emission Processes

### 1. **Purpose**
- **Regular Emission Processes**: Compute emissions based on meteorological conditions, chemical concentrations, and parameterizations
- **External Emission Data Process**: Manages pre-computed emission data from external inventories (e.g., EDGAR, CEDS, regional inventories)

### 2. **No Traditional "Schemes"**
- Instead of computational schemes, this process has **data handler schemes**:
  - `surface_2d`: For 2D surface emission data
  - `elevated_3d`: For 3D elevated emissions with vertical distribution
  - `temporal_profiles`: For time-varying emission profiles and scaling factors

### 3. **File I/O Delegation**
- File reading is delegated to the driver, not handled directly by the process
- The process provides data structures and validation for loaded data

## Architecture

### Core Components

1. **catchem_emis_data_mod.F90** (existing core module)
   - `ExtEmisDataType`: Main container for external emission data
   - `ExtEmisCategoryType`: Groups emissions by source category
   - `ExtEmisFieldType`: Individual emission fields

2. **Generated Process Files**
   - `ProcessExternalEmissionDataInterface_Mod.F90`: Main process interface
   - `ExternalEmissionDataCommon_Mod.F90`: Common data structures
   - `ExternalEmissionDataProcessCreator_Mod.F90`: Factory for creating process instances

3. **Data Handler Schemes**
   - `Surface2D`: Handles 2D surface emission data
   - `Elevated3D`: Handles 3D elevated emissions with vertical distribution
   - `TemporalProfiles`: Manages temporal scaling and profiles

## Integration with Existing catchem_emis_data_mod

The generated process works alongside the existing `catchem_emis_data_mod.F90` module:

```fortran
! The existing module provides the core data structures
USE catchem_emis_data_mod, only: ExtEmisDataType, ExtEmisCategoryType, ExtEmisFieldType

! The new process provides the process interface
USE ProcessExternalEmissionDataInterface_Mod, only: ProcessExternalEmissionDataInterface
```

## Usage Patterns

### 1. Driver Integration

```fortran
! Process validates and manages the data
call external_emission_process%initialize(config, ext_emis_data, rc)
call external_emission_process%validate_data(rc)
```

### 2. Emission Process Integration

```fortran
! Other emission processes can access external data
emission_rate = ext_emis_data%get_emission_rate('NOx', i, j, k)

! Or get column data for vectorized operations
column_ptr => ext_emis_data%get_column_ptr('anthropogenic', 'CO', i, j)
```

### 3. Temporal Interpolation

```fortran
! Update time indices for all external data
call ext_emis_data%update_time(new_time_idx, rc)

! Access temporally interpolated data
call temporal_profile_scheme%apply_scaling(base_emissions, scaled_emissions, rc)
```

## Separation of Concerns

### Current Implementation
The current approach combines external emission data management with computational emission processes, which can lead to:
- Mixed responsibilities in emission modules
- Difficulty in testing and validation
- Coupling between file I/O and emission calculations

### With External Emission Data Process
- **Data Management**: Handled by the external emission data process
- **Computational Emissions**: Handled by dedicated emission processes (biogenic, dust, etc.)
- **File I/O**: Delegated to the driver
- **Validation**: Standardized across all external data sources

## Future Extensions

### 1. Separate Processes for Specialized Emissions

#### Elevated Emissions Process
- Handle point source emissions with stack parameters
- Integrate with plume rise calculations
- Manage temporal profiles for industrial sources

#### Plume Rise Process
- Calculate effective emission heights
- Handle meteorological dependencies
- Integrate with 3D emission distribution

```yaml
# Future plume_rise_config.yaml
name: plume_rise
class_name: PlumeRise
description: Plume rise calculations for elevated point sources
process_type: emission
schemes:
  - name: briggs
    description: Briggs plume rise formulation
  - name: holland
    description: Holland plume rise formulation
```

### 2. Enhanced Data Management

#### Species Mapping
- Map inventory species to model species
- Handle speciation profiles
- Support multiple chemical mechanisms

#### Spatial Processing
- Regridding and interpolation
- Coordinate system transformations
- Boundary condition handling

## Configuration Example

```yaml
# external_emission_data_config.yaml
name: external_emission_data
schemes:
  - name: surface_2d
    parameters:
      interpolation_method: bilinear
      time_interpolation: true
      missing_value_handling: zero_fill

  - name: elevated_3d
    parameters:
      vertical_interpolation: true
      plume_rise_coupling: false

categories:
  - name: anthropogenic
    source_files:
      - EDGAR_v6_CO_2015.nc
      - EDGAR_v6_NOx_2015.nc
    scaling_factor: 1.0

  - name: biogenic
    source_files:
      - MEGAN_isoprene_2015.nc
    temporal_profiles: true
```

## Testing Strategy

### Unit Tests
- Data structure validation
- Interpolation algorithms
- Memory management
- Error handling

### Integration Tests
- Driver integration
- Process factory creation
- Multi-category data loading
- Time interpolation accuracy

### Performance Tests
- Memory usage optimization
- Large file loading
- Spatial interpolation performance
- Column access patterns

## Best Practices

1. **Modular Design**: Keep data management separate from emission calculations
2. **Validation**: Always validate external data before use
3. **Memory Management**: Monitor memory usage for large datasets
4. **Error Handling**: Provide clear error messages for file I/O issues
5. **Documentation**: Document data sources, units, and processing steps

## Migration Guide

### From Mixed Emission/Data Modules
1. Extract data management code to external emission data process
2. Keep computational emission algorithms in dedicated processes
3. Update interfaces to use standardized data access methods
4. Add proper validation and error handling

### Driver Updates
1. Implement file I/O delegation patterns
2. Add external data validation calls
3. Update time stepping to handle data interpolation
4. Enhance memory management for large datasets
