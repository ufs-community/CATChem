# Configuration API

This section covers the configuration management APIs that handle YAML-based setup and runtime configuration in CATChem.

## Overview

The Configuration system provides:

- **ConfigManager**: YAML configuration file parsing and management
- **ConfigData**: Runtime configuration data access
- **Parameter validation**: Type checking and constraint validation
- **Dynamic configuration**: Runtime parameter updates

## Core Components

### ConfigManager

The main configuration management interface:

```fortran
use ConfigManager_Mod
type(ConfigManagerType) :: config_mgr

! Load configuration from file
call config_mgr%load_config('catchem_config.yml', rc)

! Get configuration data
call config_mgr%get_config_data(config_data, rc)

! Validate configuration
call config_mgr%validate_config(rc)
```

<!-- **Auto-Generated Documentation:** [Config Manager Reference](../CATChem/namespaceconfigmanager__mod.md) -->

### ConfigData

Runtime configuration data access:

```fortran
use Config_Mod
type(ConfigDataType) :: config

! Get scalar parameters
call config%get_parameter('model.timestep', timestep, rc)
call config%get_parameter('processes.settling.enabled', enabled, rc)

! Get array parameters
call config%get_array('grid.pressure_levels', pressure_levels, rc)

! Get nested configuration
call config%get_section('processes.emission', emission_config, rc)
```

## Configuration File Format

### Basic Structure

```yaml
# CATChem Configuration File
model:
  name: "CATChem Simulation"
  timestep: 300.0  # seconds
  start_time: "2024-01-01T00:00:00"
  end_time: "2024-01-02T00:00:00"

grid:
  nx: 100
  ny: 100
  nz: 50
  pressure_levels: [1000.0, 950.0, 900.0, 850.0]  # hPa

processes:
  settling:
    enabled: true
    scheme: "stokes"
    parameters:
      slip_correction: true

  emission:
    enabled: true
    scheme: "external"
    data_file: "emissions.nc"
```

### Data Types

**Supported Types:**

- `integer` - Integer values
- `real` - Floating-point values
- `logical` - Boolean values
- `character` - String values
- `arrays` - Lists of values
- `sections` - Nested configuration objects

## Parameter Access

### Type-Safe Access

```fortran
! Integer parameters
integer :: nx, ny, nz
call config%get_parameter('grid.nx', nx, rc)
call config%get_parameter('grid.ny', ny, rc)
call config%get_parameter('grid.nz', nz, rc)

! Real parameters
real(fp) :: timestep, temperature
call config%get_parameter('model.timestep', timestep, rc)
call config%get_parameter('initial.temperature', temperature, rc)

! Logical parameters
logical :: enabled, debug_mode
call config%get_parameter('processes.settling.enabled', enabled, rc)
call config%get_parameter('debug.verbose', debug_mode, rc)

! String parameters
character(len=256) :: data_file, scheme_name
call config%get_parameter('emission.data_file', data_file, rc)
call config%get_parameter('chemistry.scheme', scheme_name, rc)
```

### Array Parameters

```fortran
! Real arrays
real(fp), allocatable :: pressure_levels(:), heights(:)
call config%get_array('grid.pressure_levels', pressure_levels, rc)
call config%get_array('grid.heights', heights, rc)

! Integer arrays
integer, allocatable :: species_indices(:)
call config%get_array('chemistry.species_indices', species_indices, rc)

! String arrays
character(len=64), allocatable :: species_names(:)
call config%get_array('chemistry.species_names', species_names, rc)
```

### Nested Configuration

```fortran
! Get configuration section
type(ConfigDataType) :: process_config
call config%get_section('processes.settling', process_config, rc)

! Access parameters within section
logical :: slip_correction
real(fp) :: density_threshold
call process_config%get_parameter('slip_correction', slip_correction, rc)
call process_config%get_parameter('density_threshold', density_threshold, rc)
```

## Validation

### Parameter Validation

```fortran
! Set validation constraints
call config%set_constraint('model.timestep', min_value=1.0_fp, &
                          max_value=3600.0_fp, rc)
call config%set_constraint('grid.nx', min_value=1, max_value=10000, rc)

! Validate configuration
call config%validate(rc)
if (rc /= CC_SUCCESS) then
    call error_mgr%report_error(ERROR_CONFIG_VALIDATION, &
                               'Configuration validation failed', rc)
endif
```

### Required Parameters

```fortran
! Mark parameters as required
call config%require_parameter('model.timestep', rc)
call config%require_parameter('grid.nx', rc)
call config%require_parameter('grid.ny', rc)

! Check required parameters
call config%check_required_parameters(rc)
```

### Custom Validation

```fortran
! Custom validation function
logical function validate_timestep(config) result(valid)
    type(ConfigDataType), intent(in) :: config
    real(fp) :: timestep

    call config%get_parameter('model.timestep', timestep, rc)
    valid = (rc == CC_SUCCESS) .and. (timestep > 0.0_fp) .and. &
            (timestep <= 3600.0_fp)
end function

! Register custom validator
call config%add_validator('model.timestep', validate_timestep, rc)
```

## Dynamic Configuration

### Runtime Updates

```fortran
! Update parameters at runtime
call config%set_parameter('model.timestep', new_timestep, rc)
call config%set_parameter('debug.verbose', .true., rc)

! Notify components of changes
call config%notify_parameter_change('model.timestep', rc)
```

### Configuration Merging

```fortran
! Load base configuration
call config_mgr%load_config('base_config.yml', rc)

! Load overlay configuration
call config_mgr%load_overlay('user_config.yml', rc)

! Merge configurations
call config_mgr%merge_configs(rc)
```

## Process Integration

### Process Configuration

```fortran
! Process-specific configuration
module MyProcess_Mod
    type :: MyProcessType
        real(fp) :: process_parameter
        logical :: enable_diagnostics
    contains
        procedure :: load_config => my_process_load_config
    end type

contains

    subroutine my_process_load_config(this, config, rc)
        class(MyProcessType), intent(inout) :: this
        type(ConfigDataType), intent(in) :: config
        integer, intent(out) :: rc

        type(ConfigDataType) :: process_config

        ! Get process-specific configuration
        call config%get_section('processes.my_process', process_config, rc)
        if (rc /= CC_SUCCESS) return

        ! Load process parameters
        call process_config%get_parameter('parameter', &
                                        this%process_parameter, rc)
        call process_config%get_parameter('enable_diagnostics', &
                                        this%enable_diagnostics, rc)
    end subroutine
end module
```

### Species Configuration

```fortran
! Chemical species configuration
type(ConfigDataType) :: species_config
character(len=64), allocatable :: species_names(:)
real(fp), allocatable :: molecular_weights(:)

call config%get_section('chemistry.species', species_config, rc)
call species_config%get_array('names', species_names, rc)
call species_config%get_array('molecular_weights', molecular_weights, rc)
```

## File Handling

### Multiple Configuration Files

```fortran
! Load from multiple files
call config_mgr%load_config('base.yml', rc)
call config_mgr%load_additional('processes.yml', rc)
call config_mgr%load_additional('species.yml', rc)
```

### Environment Variable Substitution

```yaml
# Configuration with environment variables
data:
  input_path: "${CATCHEM_DATA_DIR}/input"
  output_path: "${CATCHEM_OUTPUT_DIR}/results"

model:
  timestep: ${CATCHEM_TIMESTEP:-300.0}  # Default value
```

### Include Files

```yaml
# Main configuration
model:
  name: "CATChem Simulation"

# Include other configuration files
include:
  - "grid_config.yml"
  - "process_config.yml"
  - "species_config.yml"
```

## Error Handling

```fortran
use Error_Mod

! Configuration loading with error handling
call config_mgr%load_config('config.yml', rc)
if (rc /= CC_SUCCESS) then
    call error_mgr%report_error(ERROR_CONFIG_LOAD, &
                               'Failed to load configuration file', rc, &
                               file_name='config.yml')
    return
endif

! Parameter access with error context
call error_mgr%push_context('config_access', 'Getting model timestep')
call config%get_parameter('model.timestep', timestep, rc)
if (rc /= CC_SUCCESS) then
    call error_mgr%report_error(ERROR_PARAMETER_NOT_FOUND, &
                               'Required parameter not found', rc, &
                               parameter_name='model.timestep')
endif
call error_mgr%pop_context()
```

<!-- **Auto-Generated Documentation:** [Error Handling Reference](../CATChem/namespaceerror__mod.md) -->

## Best Practices

### Configuration Design

1. **Use hierarchical structure** for logical organization
2. **Provide sensible defaults** for optional parameters
3. **Document all parameters** with comments
4. **Group related parameters** into sections

### Parameter Access

1. **Check return codes** for all configuration operations
2. **Use type-safe access methods**
3. **Validate parameters** before use
4. **Cache frequently accessed parameters**

### Performance

1. **Load configuration once** at initialization
2. **Cache parameter values** for repeated access
3. **Use sections** to avoid repeated path parsing
4. **Minimize configuration updates** at runtime

## Configuration Examples

### Complete Model Configuration

```yaml
# Complete CATChem configuration example
model:
  name: "Regional Air Quality Simulation"
  timestep: 300.0
  start_time: "2024-01-01T00:00:00Z"
  end_time: "2024-01-03T00:00:00Z"

grid:
  projection: "lambert_conformal"
  nx: 200
  ny: 150
  nz: 40
  dx: 12000.0  # meters
  dy: 12000.0  # meters
  center_lat: 40.0
  center_lon: -95.0

chemistry:
  mechanism: "cb6r3_ae7"
  species:
    names: ["O3", "NO", "NO2", "CO", "SO2"]
    initial_concentrations: [40.0e-9, 5.0e-9, 10.0e-9, 200.0e-9, 2.0e-9]

processes:
  emission:
    enabled: true
    scheme: "external"
    data_file: "emissions/nei2017_12km.nc"

  settling:
    enabled: true
    scheme: "stokes"
    slip_correction: true

  dry_deposition:
    enabled: true
    scheme: "wesely"
    landuse_file: "landuse/nlcd_12km.nc"

output:
  frequency: 3600.0  # seconds
  variables: ["O3", "NO2", "PM25"]
  format: "netcdf"
  path: "output/"
```

## See Also

- [State Management API](state-management.md) - How configuration integrates with state
- [Process Interface API](process-interface.md) - Process configuration patterns
- [Configuration Guide](../user-guide/advanced_topics/configuration-management.md) - Advanced configuration topics
- [User Guide: Configuration](../user-guide/configuration.md) - User-level configuration help

---

<!-- **Auto-Generated Documentation:** [Complete Configuration Reference](../CATChem/namespaceconfig__manager__mod.md) -->
