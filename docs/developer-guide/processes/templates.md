# Process Templates

CATChem provides sophisticated code generation templates to streamline process development. These templates ensure consistency, reduce boilerplate code, and incorporate best practices.

## Template System Overview

The template system uses Jinja2 templating to generate:

- **C++ Process Modules**: C++ process wrapper classes (`catchem_process_<name>.hpp/cpp`) inheriting from `catchem::ProcessInterface`.
- **Fortran Science Bridges**: C-interoperable Fortran ScienceBridge subroutines (`<Name>ScienceBridge.F90`) using `BIND(C)`.
- **Fortran Common & Scheme Modules**: Process configuration modules and pure science scheme implementations.
- **Build Scripts**: Process package `CMakeLists.txt` supporting `CATCHEM_ENABLE_KOKKOS`.
- **Tests**: Standalone Fortran CTest unit tests (`test_<name>_science.f90`).
- **Documentation**: Markdown API documentation templates.

## Generator Tool

The process_generator.py tool automates process creation:

```bash
cd tools/process_generator
python process_generator.py --help
```

### Basic Usage

```bash
python process_generator.py generate \
    --config configs/dust_generic_process.yaml
```

## Template Files

### C++ Process Wrapper Template

Location: `tools/process_generator/templates/catchem_process.hpp.j2` & `catchem_process.cpp.j2`

Generates the C++ `catchem::ProcessInterface` class and registers the process with `catchem::ProcessRegistry`.

### Fortran Science Bridge Template

Location: `tools/process_generator/templates/science_bridge.F90.j2`

Generates the `BIND(C)` Fortran bridge module mapping C pointers to Fortran column arrays.

```fortran
module MyProcessScienceBridge_Mod
   use iso_c_binding
   use precision_mod, only: fp
   use MyProcessCommon_Mod
   use MyProcessScheme_DEFAULT_Mod, only: compute_default

   implicit none
   private

   public :: run_myprocess_science_bridge

contains

   subroutine run_myprocess_science_bridge( &
      n_cols, n_levels, n_species, dt, &
      u10m, v10m, conc, tendency, rc &
   ) bind(C, name="run_myprocess_science_bridge")
      integer(c_int), value :: n_cols, n_levels, n_species
      real(c_double), value :: dt
      type(c_ptr), value :: u10m, v10m, conc, tendency
      integer(c_int), intent(out) :: rc
      ! Associate pointers & run column loop
   end subroutine run_myprocess_science_bridge

end module MyProcessScienceBridge_Mod
```
{%- for species in species %}
      this%species_names({{ loop.index }}) = '{{ species }}'
{%- endfor %}
{%- endif %}

{%- if size_bins > 0 %}
      ! Size bin configuration
      call this%setup_size_bins(local_rc)
      if (local_rc /= CC_SUCCESS) then
         call error_mgr%report_error("Size bin setup failed")
         rc = local_rc
         call error_mgr%pop_context()
         return
      end if
{%- endif %}

      ! Get and validate configuration
      config => container%get_config_ptr()
      call this%load_configuration(config, local_rc)
      if (local_rc /= CC_SUCCESS) then
         call error_mgr%report_error("Configuration loading failed")
         rc = local_rc
         call error_mgr%pop_context()
         return
      end if

      call this%validate_configuration(container, local_rc)
      if (local_rc /= CC_SUCCESS) then
         call error_mgr%report_error("Configuration validation failed")
         rc = local_rc
         call error_mgr%pop_context()
         return
      end if

{%- if diagnostics %}
      ! Set up diagnostics
      call this%setup_diagnostics(container, local_rc)
      if (local_rc /= CC_SUCCESS) then
         call error_mgr%report_error("Diagnostic setup failed")
         rc = local_rc
         call error_mgr%pop_context()
         return
      end if
{%- endif %}

      this%is_initialized = .true.
      this%is_active = .true.

      write(message, '(A,A,A)') '{{ description }} initialized with scheme: ', &
                                trim(this%selected_scheme)
      call error_mgr%report_info(message)

      call error_mgr%pop_context()

   end subroutine {{ name }}_process_init

   !> Run {{ name }} process
   subroutine {{ name }}_process_run(this, container, rc)
      class({{ class_name }}), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      integer :: local_rc

      rc = CC_SUCCESS

      if (.not. this%is_ready()) then
         rc = CC_FAILURE
         return
      end if

      error_mgr => container%get_error_manager()
      call error_mgr%push_context("{{ name }}_process_run")

      ! Execute scheme-specific calculations
      select case (trim(this%selected_scheme))
{%- for scheme in schemes %}
      case ('{{ scheme }}')
         call {{ scheme|lower }}_calculate(container, this, local_rc)
{%- endfor %}
      case default
         call error_mgr%report_error("Unknown scheme: " // &
                                    trim(this%selected_scheme))
         rc = CC_FAILURE
         call error_mgr%pop_context()
         return
      end select

      if (local_rc /= CC_SUCCESS) then
         call error_mgr%report_error("Scheme calculation failed")
         rc = local_rc
         call error_mgr%pop_context()
         return
      end if

{%- if diagnostics %}
      ! Update diagnostics
      call this%update_diagnostics(container, local_rc)
      if (local_rc /= CC_SUCCESS) then
         call error_mgr%report_warning("Diagnostic update failed")
      end if
{%- endif %}

      call error_mgr%pop_context()

   end subroutine {{ name }}_process_run

   !> Finalize {{ name }} process
   subroutine {{ name }}_process_finalize(this, rc)
      class({{ class_name }}), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Deallocate arrays
      if (allocated(this%species_names)) then
         deallocate(this%species_names)
      end if

{%- if size_bins > 0 %}
      if (allocated(this%size_bin_bounds)) then
         deallocate(this%size_bin_bounds)
      end if
      if (allocated(this%size_bin_centers)) then
         deallocate(this%size_bin_centers)
      end if
{%- endif %}

      this%is_initialized = .false.
      this%is_active = .false.

   end subroutine {{ name }}_process_finalize

   ! Additional implementation methods...

end module {{ module_name }}
```

### Scheme Template

Location: `util/templates/scheme_template.f90.j2`

Generates individual scheme implementations:

```fortran
{%- set module_name = scheme_name|title + "Scheme_Mod" %}
!> \file {{ scheme_name|title }}Scheme_Mod.F90
!! \brief {{ scheme_name|title }} scheme for {{ process_name }}
!! \ingroup scheme_modules
!!
!! \author {{ author }}
!! \date {{ date }}
!! \version 1.0
!!
!! {{ scheme_description }}
!!
module {{ module_name }}
   use precision_mod, only: fp
   use {{ process_name }}Common_Mod

   implicit none
   private

   public :: compute_{{ scheme_name|lower }}

{%- if scheme_parameters %}
   ! Scheme-specific parameters
{%- for param in scheme_parameters %}
   real(fp), parameter :: {{ param.name }} = {{ param.value }}
{%- endfor %}
{%- endif %}

contains

   !> Calculate using {{ scheme_name }} algorithm
   subroutine compute_{{ scheme_name|lower }}(n_levels, n_species, dt, conc, tendency)
      integer, intent(in) :: n_levels, n_species
      real(fp), intent(in) :: dt
      real(fp), intent(in) :: conc(n_levels, n_species)
      real(fp), intent(out) :: tendency(n_levels, n_species)

      integer :: k, ispec

{%- for field in required_fields %}
      {{ field.name }} => {{ field.source }}%get_field('{{ field.name }}')
{%- endfor %}

      ! {{ scheme_name }} calculations
      do k = 1, size({{ required_fields[0].name }}, 3)
         do j = 1, size({{ required_fields[0].name }}, 2)
            do i = 1, size({{ required_fields[0].name }}, 1)

               ! Scheme-specific calculation goes here
               ! TODO: Implement {{ scheme_name }} algorithm

            end do
         end do
      end do

      call error_mgr%pop_context()

   end subroutine {{ scheme_name|lower }}_calculate

end module {{ module_name }}
```

### CMake Template

Location: `util/templates/cmake_process.txt.j2`

```cmake
# Process: {{ name }}
set(PROCESS_NAME {{ name }})

# Source files
set(PROCESS_SOURCES
    {{ name }}Process_Mod.F90
    {{ name }}Common_Mod.F90
)

# Create process library
add_library(${PROCESS_NAME}_process ${PROCESS_SOURCES})

# Set module output directory
set_target_properties(${PROCESS_NAME}_process PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
)

# Link dependencies
target_link_libraries(${PROCESS_NAME}_process
    PRIVATE
        catchem_core
        catchem_process_interface
{%- for dependency in dependencies %}
        {{ dependency }}
{%- endfor %}
)

# Include directories
target_include_directories(${PROCESS_NAME}_process
    PRIVATE
        ${CMAKE_BINARY_DIR}/modules
)

# Add schemes subdirectory
add_subdirectory(schemes)

# Link scheme libraries
target_link_libraries(${PROCESS_NAME}_process
    PRIVATE
        ${PROCESS_NAME}_schemes
)

{%- if tests %}
# Add process tests
if(ENABLE_TESTING)
    add_subdirectory(tests)
endif()
{%- endif %}
```

## Configuration Templates

### Process Configuration

Location: `util/templates/process_config.yml.j2`

```yaml
# {{ name|title }} Process Configuration
{{ name }}:
  enabled: {{ enabled|default('true') }}
  scheme: {{ default_scheme }}
  timestep: {{ timestep|default('60.0') }}  # seconds

  parameters:
{%- for param in parameters %}
    {{ param.name }}: {{ param.default }}
{%- endfor %}

{%- if species %}
  species:
{%- for species in species %}
    - {{ species }}
{%- endfor %}
{%- endif %}

{%- if size_bins > 0 %}
  size_bins:
    count: {{ size_bins }}
    type: {{ size_type|default('radius') }}  # 'radius' or 'diameter'
    bounds: [{{ size_bin_bounds|join(', ') }}]  # micrometers
{%- endif %}

{%- if diagnostics %}
  diagnostics:
{%- for diag in diagnostics %}
    - name: {{ diag.name }}
      description: "{{ diag.description }}"
      units: "{{ diag.units }}"
      output_frequency: {{ diag.frequency|default('hourly') }}
{%- endfor %}
{%- endif %}

  # Scheme-specific configurations
{%- for scheme in schemes %}
  {{ scheme|lower }}_scheme:
    # {{ scheme }} specific parameters
{%- endfor %}
```

## Test Templates

### Science Unit Test Template

Location: `tools/process_generator/templates/test_science.f90.j2`

```fortran
!> \file test_{{ name }}_science.f90
!! \brief Unit test for {{ name }} science bridge and schemes
!!
program test_{{ name }}_science
   use testing_mod, only: assert
   use iso_c_binding
   use precision_mod, only: fp

   implicit none

   ! ScienceBridge BIND(C) interface test
   call run_{{ name }}_science_bridge(...)
   call assert(.true., "{{ name }} science bridge executed successfully")

end program test_{{ name }}_science
```

   subroutine test_process_init()
      call testing_start_test("Process Initialization")

      ! Create test state container
      call create_test_state_container(container)

      ! Initialize process
      call process%init(container, rc)
      call assert_equal(rc, CC_SUCCESS, "Process initialization failed")
      call assert_true(process%is_ready(), "Process not ready after init")

      call testing_end_test()
   end subroutine test_process_init

   ! Additional test subroutines...

end program test_{{ name }}
```

## Customizing Templates

### Template Variables

Common template variables:

| Variable | Description | Example |
|----------|-------------|---------|
| `name` | Process name | `settling` |
| `description` | Process description | `Particle settling` |
| `author` | Author name | `Jane Doe` |
| `date` | Creation date | `2025-01-15` |
| `version` | Version string | `1.0` |
| `schemes` | List of schemes | `['Stokes', 'Allen']` |
| `species` | List of species | `['PM25', 'DUST']` |
| `diagnostics` | Diagnostic outputs | See config template |

### Custom Templates

You can create custom templates by:

1. Copying existing templates to a new location
2. Modifying the Jinja2 syntax and content
3. Using the `--template-dir` option with the generator

```bash
python catchem_generate_process.py \
    --name myprocess \
    --template-dir ./custom_templates \
    ...
```

## Template Best Practices

### Design Principles

1. **Consistency**: Use consistent naming and structure across all templates
2. **Flexibility**: Support various process types through template variables
3. **Completeness**: Generate all necessary files for a working process
4. **Documentation**: Include comprehensive comments and documentation
5. **Testing**: Always generate unit test scaffolding

### Template Maintenance

- **Version control**: Keep templates in version control
- **Testing**: Test template changes with multiple process types
- **Documentation**: Document template variables and usage
- **Validation**: Validate generated code compiles and runs

### Advanced Features

The template system supports:

- **Conditional blocks**: Include/exclude code based on features
- **Loops**: Generate repeated structures (species, diagnostics)
- **Filters**: Transform variables (capitalize, format)
- **Macros**: Reusable code snippets
- **Inheritance**: Template inheritance for common patterns

## Examples

### Simple Process

```bash
python catchem_generate_process.py \
    --name simple \
    --description "Simple atmospheric process" \
    --schemes basic \
    --species O3
```

### Complex Aerosol Process

```bash
python catchem_generate_process.py \
    --name aerosol \
    --description "Comprehensive aerosol process" \
    --schemes modal,sectional \
    --species PM25,PM10,SO4,NH4,NO3 \
    --size-bins 10 \
    --multiphase \
    --diagnostics mass_conc,number_conc,surface_area
```

### Chemistry Process

```bash
python catchem_generate_process.py \
    --name chemistry \
    --description "Gas-phase atmospheric chemistry" \
    --schemes CB6,RACM2 \
    --species O3,NO,NO2,CO,SO2,HCHO \
    --diagnostics reaction_rates,photolysis_rates \
    --parameters kinetic_solver:LSODES,relative_tolerance:1e-3
```

The template system significantly accelerates process development while ensuring consistency and best practices across all CATChem processes.
