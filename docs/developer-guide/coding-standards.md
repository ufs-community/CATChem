# Coding Standards

These coding standards are the foundation of the CATChem project. They are designed to ensure that the codebase is consistent, maintainable, readable, and performant. All contributions to CATChem must adhere to these standards.

## Guiding Principles

- **Clarity over cleverness**: Write code that is easy for others to understand.
- **Consistency is key**: A consistent coding style makes the codebase easier to read and maintain.
- **Performance matters**: CATChem is a high-performance computing application, and code should be written with efficiency in mind.
- **Document your work**: Good documentation is as important as good code.

## Fortran Language Standards

- **Language Version**: CATChem targets the **Fortran 2008** standard, with some selected features from Fortran 2018.
- **Compiler Support**: The code must compile and run correctly with the latest versions of `gfortran` and `ifort`.
- **Modern Fortran**: Use modern Fortran features such as `allocatable` arrays, derived types with type-bound procedures, and modules. Avoid obsolete features like `COMMON` blocks and `GOTO` statements.

```fortran
! Good: Use modern Fortran constructs
use, intrinsic :: iso_fortran_env, only: real64
implicit none

real(real64), allocatable :: my_array(:,:)

type :: MyType
contains
  procedure :: my_procedure
end type
```

## Naming Conventions

Consistent naming is crucial for readability.

- **Modules**: `PascalCase` with a `_Mod` suffix (e.g., `SettlingProcess_Mod`).
- **Derived Types**: `PascalCase` with a `Type` suffix (e.g., `SettlingProcessType`).
- **Procedures (Subroutines and Functions)**: `snake_case` (e.g., `compute_settling_velocity`).
- **Variables**: `snake_case` (e.g., `num_species`, `settling_velocity`).
- **Constants**: `UPPER_CASE` with a `_` separator (e.g., `MAX_SPECIES`, `GRAVITATIONAL_ACCELERATION`).

## Code Structure and Formatting

- **Indentation**: Use 2 spaces for indentation. Do not use tabs.
- **Line Length**: Keep lines under 100 characters where possible.
- **Module Structure**: Organize modules with a clear separation of the public interface and the private implementation.

```fortran
module MyProcess_Mod
  use precision_mod, only: fp
  implicit none
  private

  ! Public interface
  public :: MyProcessType

  ! Type definition
  type, public :: MyProcessType
    private
    ! Private data members
  contains
    procedure :: init => my_process_init
    procedure :: run => my_process_run
  end type

contains

  subroutine my_process_init(this, ...)
    ! Implementation
  end subroutine

  subroutine my_process_run(this, ...)
    ! Implementation
  end subroutine

end module MyProcess_Mod
```

- **Procedure Structure**: All procedures should have explicit `intent` for all arguments and a clear documentation block.

## Documentation Standards

All code must be documented using Doxygen-style comments. See the [Documentation Guide](documentation.md) for more details.

- **File Headers**: Every source file should begin with a header block that describes the file's purpose.
- **Module Documentation**: Every module should have a documentation block explaining its role.
- **Procedure Documentation**: Every public procedure must have a Doxygen block that describes its purpose, parameters, and return values.
- **Inline Comments**: Use inline comments to explain complex or non-obvious parts of the code. Focus on the *why*, not the *what*.

## Error Handling

Robust error handling is critical.

- **Use the `ErrorManagerType`**: For all new code, use the modern error handling system as described in the [Error Handling Guide](core/error-handling.md).
- **Return Codes**: All subroutines that can fail must have a return code argument (`rc`).
- **Check Return Codes**: Always check the return code of any subroutine that you call.
- **Provide Context**: Use the `push_context` and `pop_context` procedures of the `ErrorManagerType` to provide a clear context for errors.

## Performance Guidelines

- **Memory Management**:
    - Use `allocatable` arrays instead of `pointer` arrays where possible.
    - Allocate memory once and reuse it if possible, rather than allocating and deallocating within a loop.
- **Loop Optimization**:
    - The inner-most loop should be over the first index of an array (column-major order).
    - Use compiler directives (e.g., `!$OMP SIMD`) to encourage vectorization for performance-critical loops.
- **Precision**:
    - Use the `fp` kind specifier from `precision_mod` for all real literals.
    - Use named constants instead of magic numbers.

## Testing Standards

All contributions must include tests.

- **Unit Tests**: Every new or modified public procedure must have a corresponding unit test.
- **Test Coverage**: Tests should cover normal operation, edge cases, and error conditions.
- **Test Structure**: Follow the structure of existing tests in the `tests/` directory.

See the [Testing Guide](testing.md) for more details.
