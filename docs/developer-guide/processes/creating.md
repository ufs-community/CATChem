# Creating a New Column-Based Process

This guide provides a step-by-step tutorial for creating a new column-based atmospheric process in CATChem. We will use a simplified version of the `seasalt` process as an example.

## 1. Directory and File Structure

First, create a directory for your new process within `src/process`. The structure should look like this:

```
src/process/newprocess/
├── schemes/
│   ├── NewProcessScheme_Mod.F90
│   └── CMakeLists.txt
├── NewProcessInterface_Mod.F90
├── NewProcessCommon_Mod.F90
├── NewProcessCreator_Mod.F90
└── CMakeLists.txt
```

-   `NewProcessInterface_Mod.F90`: The main process module, extending `ColumnProcessInterface`.
-   `NewProcessCommon_Mod.F90`: A module for common data structures and utilities, like the process configuration.
-   `NewProcessCreator_Mod.F90`: A module responsible for creating an instance of the process.
-   `schemes/`: A directory to hold the different algorithmic implementations (schemes) for the process.

## 2. The Process Interface

The process interface is the main entry point for the process. It extends `ColumnProcessInterface` for column-based processes.

`NewProcessInterface_Mod.F90`:
```fortran
module NewProcessInterface_Mod
  use precision_mod, only: fp
  use ProcessInterface_Mod, only: ColumnProcessInterface
  use StateManager_Mod, only: StateManagerType
  use VirtualColumn_Mod, only: VirtualColumnType
  use NewProcessCommon_Mod, only: NewProcessConfigType
  use NewProcessScheme_Mod, only: compute_scheme

  implicit none
  private
  public :: NewProcessInterfaceType

  type, extends(ColumnProcessInterface) :: NewProcessInterfaceType
    private
    type(NewProcessConfigType) :: config
  contains
    procedure :: init => newprocess_init
    procedure :: run => newprocess_run
    procedure :: finalize => newprocess_finalize
    procedure :: run_column => newprocess_run_column
    procedure :: get_required_met_fields => newprocess_get_required_met_fields
  end type NewProcessInterfaceType

contains

  subroutine newprocess_init(this, container, rc)
    class(NewProcessInterfaceType), intent(inout) :: this
    type(StateManagerType), intent(inout) :: container
    integer, intent(out) :: rc
    ! Get config, load species, etc.
  end subroutine

  subroutine newprocess_run(this, container, rc)
    class(NewProcessInterfaceType), intent(inout) :: this
    type(StateManagerType), intent(inout) :: container
    integer, intent(out) :: rc
    ! No 3D operations needed for this example
    rc = 0
  end subroutine

  subroutine newprocess_finalize(this, rc)
    class(NewProcessInterfaceType), intent(inout) :: this
    integer, intent(out) :: rc
    ! Deallocate resources
    rc = 0
  end subroutine

  subroutine newprocess_run_column(this, column, container, rc)
    class(NewProcessInterfaceType), intent(inout) :: this
    type(VirtualColumnType), intent(inout) :: column
    type(StateManagerType), intent(inout) :: container
    integer, intent(out) :: rc

    ! Get met data from the virtual column
    real(fp) :: u10m, v10m, sst
    u10m = column%get_met()%U10M
    v10m = column%get_met()%V10M
    sst = column%get_met()%SST

    ! Call the scheme
    call compute_scheme(u10m, v10m, sst, ...)

    ! Update chemical state
    ! ...
  end subroutine

  function newprocess_get_required_met_fields(this) result(field_names)
    class(NewProcessInterfaceType), intent(in) :: this
    character(len=32), allocatable :: field_names(:)
    allocate(field_names(3))
    field_names(1) = 'U10M'
    field_names(2) = 'V10M'
    field_names(3) = 'SST'
  end function newprocess_get_required_met_fields

end module NewProcessInterface_Mod
```

### Key Points:

-   **`extends(ColumnProcessInterface)`**: This is crucial for column-based processes. It provides the `run_column` method and other column-related utilities.
-   **`get_required_met_fields`**: This function tells the `StateManager` which meteorological fields this process needs. The fields will then be available in the `VirtualColumnType`'s `met` object.
-   **`run_column`**: This is where the core logic for a single column is executed. You can get met data and chemical species from the `column` object, call your scheme, and then update the chemical state.

## 3. The Scheme

The scheme contains the actual scientific algorithm.

`schemes/NewProcessScheme_Mod.F90`:
```fortran
module NewProcessScheme_Mod
  use precision_mod, only: fp
  implicit none
  private
  public :: compute_scheme

contains

  subroutine compute_scheme(u10m, v10m, sst, tendency, rc)
    real(fp), intent(in) :: u10m, v10m, sst
    real(fp), intent(out) :: tendency
    integer, intent(out) :: rc

    ! Calculate tendency based on met inputs
    tendency = (u10m**2 + v10m**2) * sst * 1.0e-5
    rc = 0
  end subroutine

end module NewProcessScheme_Mod
```

## 4. Configuration and Creator

You'll need a way to configure and create your process.

`NewProcessCommon_Mod.F90`:
```fortran
module NewProcessCommon_Mod
  use precision_mod, only: fp
  implicit none
  private
  public :: NewProcessConfigType

  type :: NewProcessConfigType
    character(len=32) :: scheme = 'default'
    logical :: diagnostics = .false.
  end type NewProcessConfigType

end module
```

`NewProcessCreator_Mod.F90`:
```fortran
module NewProcessCreator_Mod
  use ProcessInterface_Mod, only: ProcessInterface
  use NewProcessInterface_Mod, only: NewProcessInterfaceType
  implicit none
  private
  public :: create_newprocess

contains

  function create_newprocess() result(process)
    type(ProcessInterface), pointer :: process
    type(NewProcessInterfaceType), pointer :: new_process
    allocate(new_process)
    process => new_process
  end function create_newprocess

end module
```

## 5. Build System (CMake)

Update the CMake files to include your new process.

`src/process/newprocess/CMakeLists.txt`:
```cmake
add_library(newprocess_interface OBJECT ${CMAKE_CURRENT_SOURCE_DIR}/NewProcessInterface_Mod.F90)
add_library(newprocess_common OBJECT ${CMAKE_CURRENT_SOURCE_DIR}/NewProcessCommon_Mod.F90)
add_library(newprocess_creator OBJECT ${CMAKE_CURRENT_SOURCE_DIR}/NewProcessCreator_Mod.F90)

add_subdirectory(schemes)

add_library(newprocess STATIC)
target_sources(newprocess PRIVATE
    $<TARGET_OBJECTS:newprocess_interface>
    $<TARGET_OBJECTS:newprocess_common>
    $<TARGET_OBJECTS:newprocess_creator>
    $<TARGET_OBJECTS:newprocess_schemes>
)

target_link_libraries(newprocess PUBLIC catchem_core)
```

`src/process/newprocess/schemes/CMakeLists.txt`:
```cmake
add_library(newprocess_schemes OBJECT ${CMAKE_CURRENT_SOURCE_DIR}/NewProcessScheme_Mod.F90)
```

Finally, add your process to `src/process/CMakeLists.txt`:
```cmake
add_subdirectory(newprocess)
list(APPEND PROCESS_LIBRARIES newprocess)
```

And register it in `src/core/ProcessRegistry_Mod.F90`.

## 6. Testing

Create a new test file in the `tests/` directory to test your process. You can use the `test_ProcessFactory.f90` as a starting point to see how to create and run a process.

This guide provides a basic skeleton. For more advanced features like diagnostics and handling multiple species, refer to the `seasalt` process implementation in `src/process/seasalt`.
