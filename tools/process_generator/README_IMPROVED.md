# CATChem Process Generator - Improved Templates

## Overview

The CATChem Process Generator has been enhanced with improved templates that leverage the core CATChem infrastructure and follow modern software architecture best practices. These improvements significantly reduce boilerplate code, improve maintainability, and provide better separation of concerns.

## Key Features

### 🏗️ **Architecture Improvements**
- **Inheritance-based Design**: Templates extend core `ProcessInterface` classes
- **Utility Integration**: Leverage existing `ChemSpeciesUtils`, `UnitConversion`, and `MetFieldUtils`
- **Column Virtualization**: Automatic support for 1D column processing
- **Configuration Management**: Structured YAML-based configuration parsing

### 🚀 **Development Benefits**
- **Reduced Code Generation**: ~70% less generated boilerplate code
- **Simplified Maintenance**: Changes to infrastructure don't require template updates
- **Consistent Patterns**: All processes follow the same architectural patterns
- **Better Testing**: Clear separation enables focused unit testing

### 🔬 **Science Focus**
- **Pure Science Kernels**: Scheme modules remain infrastructure-independent
- **Algorithm Portability**: Science code can be reused in other frameworks
- **Performance Optimization**: Optimized data access through core utilities

## Template Structure

### Process Interface Template (`process_interface.F90.j2`)

```jinja
module Process{{ config.class_name }}Interface_Mod
   ! Core infrastructure imports
   use ProcessInterface_Mod, only: ProcessInterface, ColumnProcessInterface
   use ChemSpeciesUtils_Mod, only: ChemSpeciesUtilsType
   use UnitConversion_Mod, only: UnitConverterType
   use MetFieldUtils_Mod, only: MetFieldUtilsType

   ! Process extends base infrastructure
   type, extends(ProcessInterface) :: Process{{ config.class_name }}Interface
      private
      ! Only process-specific configuration
      character(len=64) :: active_scheme
      type(ChemSpeciesUtilsType) :: species_utils
      type(UnitConverterType) :: unit_converter
      type(MetFieldUtilsType) :: met_utils
   contains
      ! Minimal required implementations
      procedure :: init => process_init
      procedure :: run => process_run
      procedure :: finalize => process_finalize
   end type
```

### Scheme Module Template (`scheme_module.F90.j2`)

The scheme template remains focused on pure science:

```jinja
module {{ config.class_name }}Scheme_{{ scheme.class_name }}_Mod
   use iso_fortran_env, only: fp => real64

   ! Pure science interface
   public :: compute_{{ scheme.name }}
   public :: {{ scheme.name }}_params_t

   ! Science-only implementation
   pure subroutine compute_{{ scheme.name }}(...)
      ! Pure computational kernel
   end subroutine
```

## Configuration System

### Process Configuration Example

```yaml
processes:
  dust_emission:
    name: "dust_emission"
    description: "Mineral dust emission process"
    process_type: "emission"
    enable_column_processing: true

    schemes:
      - name: "fengsha"
        description: "Fengsha dust emission scheme"
        author: "NOAA/NCEP"
        reference: "Feng et al. (2007)"

        parameters:
          threshold_velocity:
            value: 0.2
            units: "m/s"
            description: "Threshold wind velocity"

          emission_factor:
            value: 1.0e-6
            units: "kg/m²/s"
            description: "Base emission factor"

        required_met_fields:
          - name: "surface_wind_speed"
            dimensions: "scalar"
            units: "m/s"

          - name: "soil_moisture"
            dimensions: "1d"
            units: "m³/m³"

    species: ["DUST_1", "DUST_2", "DUST_3", "DUST_4", "DUST_5"]

    diagnostics:
      - name: "dust_emission_flux"
        description: "Total dust emission flux"
        units: "kg/m²/s"

      - name: "dust_threshold_velocity"
        description: "Threshold velocity diagnostic"
        units: "m/s"
```

## Generated Code Examples

### Before (Old Template)
```fortran
! Manual data marshaling
allocate(temperature(n_levels))
allocate(wind_speed(n_levels))
allocate(species_conc(this%n_species))
allocate(species_tendencies(this%n_species))

! Manual field retrieval
call state_manager%get_met_field('temperature', i_col, temperature, error_handler)
if (error_handler%has_error()) return

call state_manager%get_met_field('wind_speed', i_col, wind_speed, error_handler)
if (error_handler%has_error()) return

! Manual species concentration retrieval
do i = 1, this%n_species
   call state_manager%get_species_concentration(this%species_indices(i), i_col, species_conc(i:i), error_handler)
   if (error_handler%has_error()) return
end do

! Manual tendency application
do i_spec = 1, this%n_species
   do i_col = 1, n_columns
      call state_manager%get_species_concentration(this%species_indices(i_spec), i_col, current_conc, error_handler)
      new_conc(1) = current_conc(1) + this%process_tendencies(i_spec, i_col) * dt
      new_conc(1) = max(0.0_fp, new_conc(1))
      call state_manager%set_species_concentration(this%species_indices(i_spec), i_col, new_conc, error_handler)
   end do
end do
```

### After (New Template)
```fortran
! Utility-based data access
call this%met_utils%get_field('temperature', virtual_column, temperature, error_manager)
call this%met_utils%get_field('wind_speed', virtual_column, wind_speed, error_manager)
call this%species_utils%get_concentrations(species_indices, state_manager, i_col, species_conc, error_manager)

! Base class handles tendency application
call this%apply_tendency(species_indices, state_manager, i_col, species_tendencies, dt, error_manager)
```

## Usage Instructions

### 1. Process Definition

Create a YAML configuration file defining your process:

```bash
# Create process configuration
cat > my_process.yaml << EOF
processes:
  my_process:
    name: "my_process"
    description: "My atmospheric process"
    # ... configuration details
EOF
```

### 2. Generate Process Interface

```bash
# Generate process interface code
cd tools/process_generator
python generate_process.py --config my_process.yaml --output ../../src/process/my_process/
```

### 3. Implement Science Schemes

The generator creates skeleton scheme modules. Implement the science algorithms:

```fortran
pure subroutine compute_my_scheme(n_levels, n_species, params, met_fields, species_conc, species_tendencies)
   ! Your science algorithm here
   ! TODO: Replace with actual physics
end subroutine
```

### 4. Build and Test

```bash
# Build with CMake
mkdir build && cd build
cmake .. -DENABLE_MY_PROCESS=ON
make

# Run tests
ctest -R my_process
```

## Advanced Features

### Column Processing Support

Enable column virtualization for better performance:

```yaml
processes:
  my_process:
    enable_column_processing: true  # Enable 1D column processing
    column_batch_size: 100         # Process columns in batches
```

Generated code automatically supports both modes:

```fortran
{% if config.enable_column_processing %}
   ! Use column virtualization
   call this%run_column_batch(state_manager, dt, error_manager)
{% else %}
   ! Use standard 3D processing
   call this%run_3d_processing(state_manager, dt, error_manager)
{% endif %}
```

### Multi-Scheme Support

Define multiple schemes within a process:

```yaml
schemes:
  - name: "scheme_a"
    description: "Fast scheme for operational use"
  - name: "scheme_b"
    description: "Detailed scheme for research"
```

Runtime scheme selection:

```fortran
select case (trim(this%active_scheme))
case ('scheme_a')
   call this%run_scheme_a_scheme(virtual_column, dt, error_manager)
case ('scheme_b')
   call this%run_scheme_b_scheme(virtual_column, dt, error_manager)
end select
```

### Diagnostic Integration

Automatic diagnostic registration and updates:

```yaml
diagnostics:
  - name: "emission_rate"
    description: "Emission rate diagnostic"
    units: "kg/m²/s"
```

Generated diagnostic handling:

```fortran
! Base class handles diagnostic storage and updates
call this%update_process_diagnostics(species_tendencies, virtual_column, error_manager)
```

## Best Practices

### 1. **Science Module Design**
- Keep scheme modules pure science kernels
- No CATChem infrastructure dependencies in schemes
- Use only basic Fortran types for portability
- Focus on algorithms, let infrastructure handle data management

### 2. **Configuration Organization**
- Use descriptive names for schemes and parameters
- Include units and descriptions for all parameters
- Group related parameters logically
- Validate parameter ranges in scheme code

### 3. **Error Handling**
- Use `ErrorManagerType` consistently
- Check for errors after each utility call
- Provide meaningful error messages
- Let base class handle common error patterns

### 4. **Performance Optimization**
- Enable column processing for performance-critical processes
- Use utility functions for data access (they're optimized)
- Minimize allocations in hot loops
- Leverage base class optimizations

### 5. **Testing Strategy**
- Test scheme modules independently (unit tests)
- Test process interfaces with mock data (integration tests)
- Use CATChem testing framework for validation
- Compare against reference implementations

## Migration from Old Templates

### Automatic Migration

For existing processes generated with old templates:

```bash
# Regenerate with new templates
python tools/process_generator/migrate_process.py --process my_old_process --output src/process/my_process/
```

### Manual Migration Steps

1. **Update Type Definition**:
   ```fortran
   ! Old
   type :: ProcessMyInterface

   ! New
   type, extends(ProcessInterface) :: ProcessMyInterface
   ```

2. **Add Utility Members**:
   ```fortran
   type(ChemSpeciesUtilsType) :: species_utils
   type(UnitConverterType) :: unit_converter
   type(MetFieldUtilsType) :: met_utils
   ```

3. **Simplify Initialization**:
   ```fortran
   ! Call base class initialization
   call this%ProcessInterface%init(config_data, state_manager, error_manager)

   ! Initialize utilities
   call this%species_utils%init(state_manager, error_manager)
   ```

4. **Replace Data Access**:
   ```fortran
   ! Old manual access
   call state_manager%get_met_field('temperature', i_col, temperature, error_handler)

   ! New utility access
   call this%met_utils%get_field('temperature', virtual_column, temperature, error_manager)
   ```

## Troubleshooting

### Common Issues

1. **Missing Core Modules**: Ensure CATChem core modules are built first
2. **Template Errors**: Check YAML configuration syntax
3. **Missing Utilities**: Verify utility modules are available in core
4. **Column Processing**: Ensure `enable_column_processing` is correctly set

### Debug Tips

1. **Enable Verbose Output**: Use `--verbose` flag with process generator
2. **Check Generated Code**: Review generated files before building
3. **Test Incrementally**: Test scheme modules independently first
4. **Use CATChem Debugging**: Enable CATChem debug modes for runtime issues

## Future Enhancements

- **Automatic Testing Generation**: Generate unit tests alongside process code
- **Performance Profiling Integration**: Built-in performance monitoring
- **External Library Integration**: Simplified integration with MICM, ESMF, etc.
- **Interactive Configuration**: Web-based configuration editor
- **Validation Framework**: Automatic validation against reference data

## Contributing

Contributions to improve the process generator templates are welcome:

1. **Template Improvements**: Enhance generated code quality
2. **New Features**: Add support for new process types or capabilities
3. **Documentation**: Improve examples and guides
4. **Testing**: Expand test coverage and validation

See `CONTRIBUTING.md` for detailed contribution guidelines.
