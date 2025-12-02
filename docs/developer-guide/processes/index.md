# Process Development

This section covers developing new processes and schemes for CATChem. See the **[CATChem User Guide](../../user-guide/index.md#process-documentation)** for a description of all processes available in CATChem.

## Quick Links

- **[Architecture Overview](architecture.md)** - Understanding process design patterns
- **[Process Generator](process-generator.md)** - Complete guide to using the automated process generator
- **[Creating Custom Processes](creating.md)** - Manual process development
- **[Templates and Patterns](templates.md)** - Code templates and best practices
- **[Testing Processes](testing.md)** - Testing strategies and frameworks

## Overview

CATChem processes are modular components that implement specific atmospheric transport, chemical, emission, or loss schemes. Each process follows a standardized interface and lifecycle.

## Process Architecture

All processes inherit from the base `ProcessInterface` and must implement `init`, `run`, and `finalize`. Column-based processes should extend `ColumnProcessInterface`.

```fortran
module MyProcess_Mod
  use ProcessInterface_Mod, only: ColumnProcessInterface
  use StateManager_Mod, only: StateManagerType
  use VirtualColumn_Mod, only: VirtualColumnType
  implicit none

  type, extends(ColumnProcessInterface) :: MyProcessType
    ! Private data members
  contains
    procedure :: init => MyProcess_init
    procedure :: run => MyProcess_run
    procedure :: finalize => MyProcess_finalize
    procedure :: run_column => MyProcess_run_column
  end type

contains

  subroutine MyProcess_init(this, container, rc)
    class(MyProcessType), intent(inout) :: this
    type(StateManagerType), intent(inout) :: container
    integer, intent(out) :: rc
    ! Initialize process
  end subroutine

  subroutine MyProcess_run(this, container, rc)
    class(MyProcessType), intent(inout) :: this
    type(StateManagerType), intent(inout) :: container
    integer, intent(out) :: rc
    ! 3D operations (before/after column processing)
  end subroutine

  subroutine MyProcess_finalize(this, rc)
    class(MyProcessType), intent(inout) :: this
    integer, intent(out) :: rc
    ! Clean up resources
  end subroutine

  subroutine MyProcess_run_column(this, column, container, rc)
    class(MyProcessType), intent(inout) :: this
    type(VirtualColumnType), intent(inout) :: column
    type(StateManagerType), intent(inout) :: container
    integer, intent(out) :: rc
    ! Process a single column
  end subroutine

end module
```

## Creating New Processes

To create a new process, follow the [Creating Custom Processes](creating.md) guide.

## Scheme Development

Schemes implement specific algorithms within a process. See the `seasalt` process for an example of how to structure schemes.

## Testing Processes

All new processes must include unit and integration tests. See the [Testing Processes](testing.md) guide for more information.

## Documentation

All new processes must be documented using Doxygen-style comments.

## See Also

- [Process Architecture](architecture.md)
- [Testing Processes](testing.md)
- [Process Templates](templates.md)
