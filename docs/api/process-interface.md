# Process Interface API

This section covers the process interface APIs that define how atmospheric processes integrate with the CATChem framework.

## Overview

The Process Interface system provides:

- **ProcessInterface**: Base interface for all atmospheric processes
- **ColumnProcessInterface**: Column-optimized process interface
- **ProcessManager**: Process orchestration and lifecycle management
- **ProcessRegistry**: Dynamic process discovery and loading

## Core Interfaces

### ProcessInterface

The base interface that all processes must implement:

```fortran
use ProcessInterface_Mod

type, extends(ProcessInterface_t) :: MyProcessType
contains
    procedure :: initialize => my_process_init
    procedure :: run => my_process_run
    procedure :: finalize => my_process_finalize
end type
```

**Required Methods:**

- `initialize(container, rc)` - Setup process with configuration
- `run(container, rc)` - Execute process calculations
- `finalize(rc)` - Clean up process resources

**Auto-Generated Documentation:** [Process Interface Reference](../CATChem/namespaceprocessinterface__mod.md)

### ColumnProcessInterface

Optimized interface for column-based processing:

```fortran
use ColumnProcessInterface_Mod

type, extends(ColumnProcessInterface_t) :: MyColumnProcessType
contains
    procedure :: run_column => my_column_process
    procedure :: can_use_column_processing => my_column_check
end type

subroutine my_column_process(this, column, rc)
    class(MyColumnProcessType), intent(inout) :: this
    type(ColumnType), intent(inout) :: column
    integer, intent(out) :: rc

    ! Process atmospheric column
    do k = 1, column%nz
        ! Atmospheric level calculations
        call process_level(column, k, rc)
    end do
end subroutine
```

**Performance Benefits:**

- Automatic parallelization
- Cache-optimized memory access
- Linear scaling with grid size

**Auto-Generated Documentation:** [Column Interface Reference](../CATChem/namespacecolumninterface__mod.md)

## Process Implementation

### Basic Process Structure

```fortran
module MyProcess_Mod
    use ProcessInterface_Mod
    use State_Mod
    use Error_Mod
    implicit none

    type, extends(ProcessInterface_t) :: MyProcessType
        private
        ! Process-specific data
        real(fp) :: process_parameter
        logical :: is_initialized = .false.
    contains
        procedure :: initialize => my_process_init
        procedure :: run => my_process_run
        procedure :: finalize => my_process_finalize
    end type

contains

    subroutine my_process_init(this, container, rc)
        class(MyProcessType), intent(inout) :: this
        type(StateContainerType), intent(in) :: container
        integer, intent(out) :: rc

        rc = CC_SUCCESS

        ! Initialize process
        call this%load_configuration(container, rc)
        if (rc /= CC_SUCCESS) return

        ! Setup process state
        this%is_initialized = .true.
    end subroutine

    subroutine my_process_run(this, container, rc)
        class(MyProcessType), intent(inout) :: this
        type(StateContainerType), intent(inout) :: container
        integer, intent(out) :: rc

        rc = CC_SUCCESS

        if (.not. this%is_initialized) then
            rc = ERROR_NOT_INITIALIZED
            return
        endif

        ! Process implementation
        call this%execute_calculations(container, rc)
    end subroutine

end module
```

### Column-Based Process

```fortran
module MyColumnProcess_Mod
    use ColumnProcessInterface_Mod
    use ColumnInterface_Mod
    implicit none

    type, extends(ColumnProcessInterface_t) :: MyColumnProcessType
        private
        real(fp) :: column_parameter
    contains
        procedure :: run_column => my_column_run
        procedure :: can_use_column_processing => my_column_capable
    end type

contains

    subroutine my_column_run(this, column, rc)
        class(MyColumnProcessType), intent(inout) :: this
        type(ColumnType), intent(inout) :: column
        integer, intent(out) :: rc

        integer :: k

        rc = CC_SUCCESS

        ! Process each atmospheric level
        do k = 1, column%nz
            call this%process_level(column, k, rc)
            if (rc /= CC_SUCCESS) return
        end do

        ! Apply boundary conditions
        call this%apply_boundary_conditions(column, rc)
    end subroutine

    logical function my_column_capable(this)
        class(MyColumnProcessType), intent(in) :: this
        my_column_capable = .true.  ! This process supports column processing
    end function

end module
```

## Process Manager

The ProcessManager orchestrates process execution:

```fortran
use ProcessManager_Mod

type(ProcessManagerType) :: process_mgr
type(StateContainerType) :: container

! Initialize process manager
call process_mgr%init(config, rc)

! Register processes
call process_mgr%register_process('my_process', my_process_factory, rc)

! Execute all processes
call process_mgr%run_all_processes(container, rc)

! Clean up
call process_mgr%finalize(rc)
```

**Key Features:**

- Automatic process discovery
- Dependency resolution
- Parallel process execution
- Error recovery and reporting

**Auto-Generated Documentation:** [Process Manager Reference](../CATChem/namespaceprocessmanager__mod.md)

## Process Registry

Dynamic process loading and management:

```fortran
use ProcessRegistry_Mod

! Register process factory
call process_registry%register('emission_process', &
                              create_emission_process, rc)

! Create process instance
call process_registry%create_process('emission_process', &
                                   process_instance, rc)

! List available processes
call process_registry%list_processes(process_names, rc)
```

**Auto-Generated Documentation:** [Process Registry Reference](../CATChem/namespaceprocessregistry__mod.md)

## Error Handling

Comprehensive error handling for process operations:

```fortran
use Error_Mod

subroutine my_process_run(this, container, rc)
    ! ...

    ! Use error context for debugging
    call error_mgr%push_context('my_process_run', &
                               'Executing atmospheric process')

    ! Process operations with error checking
    call this%calculate_tendencies(container, rc)
    if (rc /= CC_SUCCESS) then
        call error_mgr%report_error(ERROR_CALCULATION, &
                                   'Failed to calculate tendencies', rc, &
                                   additional_info='Check input data validity')
        call error_mgr%pop_context()
        return
    endif

    call error_mgr%pop_context()
end subroutine
```

**Auto-Generated Documentation:** [Error Handling Reference](../CATChem/namespaceerror__mod.md)

## Configuration Integration

Processes integrate with the configuration system:

```fortran
! Process configuration
subroutine load_process_config(this, container, rc)
    class(MyProcessType), intent(inout) :: this
    type(StateContainerType), intent(in) :: container
    integer, intent(out) :: rc

    type(ConfigDataType) :: config

    ! Get process configuration
    call container%get_config_data(config, rc)
    if (rc /= CC_SUCCESS) return

    ! Load process-specific parameters
    call config%get_parameter('my_process.parameter1', &
                             this%parameter1, rc)
    call config%get_parameter('my_process.parameter2', &
                             this%parameter2, rc)
end subroutine
```

## Diagnostic Integration

Processes can output diagnostics:

```fortran
! Diagnostic output
subroutine process_with_diagnostics(this, container, rc)
    class(MyProcessType), intent(inout) :: this
    type(StateContainerType), intent(inout) :: container
    integer, intent(out) :: rc

    real(fp) :: diagnostic_value

    ! Calculate process diagnostics
    diagnostic_value = this%calculate_diagnostic()

    ! Output to diagnostic system
    call container%set_diagnostic('my_process_rate', &
                                 diagnostic_value, rc)
end subroutine
```

**Auto-Generated Documentation:** [Diagnostic Manager Reference](../CATChem/namespacediagnosticmanager__mod.md)

**Auto-Generated Documentation:** [Diagnostic Interface Reference](../CATChem/namespacediagnosticinterface__mod.md)

## Testing Framework

Process testing utilities:

```fortran
! Unit test example
program test_my_process
    use MyProcess_Mod
    use TestFramework_Mod
    implicit none

    type(MyProcessType) :: process
    type(StateContainerType) :: container
    integer :: rc

    ! Setup test
    call test_framework%setup('my_process_test', rc)
    call container%init_for_testing(rc)

    ! Initialize process
    call process%initialize(container, rc)
    call test_framework%assert_success(rc, 'Process initialization')

    ! Run process
    call process%run(container, rc)
    call test_framework%assert_success(rc, 'Process execution')

    ! Validate results
    call test_framework%validate_output(container, 'expected_output.nc', rc)

    ! Cleanup
    call process%finalize(rc)
    call test_framework%cleanup(rc)
end program
```

## Best Practices

### Performance

1. **Use column processing** when possible
2. **Minimize memory allocations** in run methods
3. **Cache frequently accessed data**
4. **Use efficient algorithms** for atmospheric calculations

### Code Quality

1. **Follow interface contracts** exactly
2. **Handle all error conditions** appropriately
3. **Document process physics** clearly
4. **Write comprehensive tests**

### Integration

1. **Use StateContainer** for all data access
2. **Register diagnostics** for monitoring
3. **Support configuration** parameters
4. **Handle initialization** dependencies

## Process Templates

Use the process generator for rapid development:

```bash
# Generate new process template
python util/catchem_generate_process.py \
    --name MyProcess \
    --type transport \
    --column-capable \
    --output-dir src/process/myprocess/
```

This creates a complete process template with:

- Process interface implementation
- Column processing support
- Configuration integration
- Error handling
- Unit tests
- Documentation

## See Also

- [State Management API](state-management.md) - Data access and manipulation
- [Column Interface API](column-interface.md) - Column-based processing
- [Configuration API](configuration.md) - Process configuration
- [Process Development Guide](../developer-guide/processes/creating.md) - Step-by-step process creation

---

**Auto-Generated Documentation:** [Complete Process Interface Reference](../CATChem/namespaceprocessinterface__mod.md)
