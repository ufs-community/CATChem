# Testing Guide

Comprehensive testing is a cornerstone of the CATChem development process. It is essential for ensuring the scientific accuracy, numerical stability, and software reliability of the model. This guide provides an overview of the testing philosophy, the different types of tests, and how to write and run tests for CATChem.

## Testing Philosophy

Our testing philosophy is based on a multi-layered approach, where each layer of testing builds upon the previous one.

- **Correctness**: Does the code produce the correct scientific results?
- **Robustness**: Does the code handle edge cases and invalid inputs gracefully?
- **Maintainability**: Are the tests easy to understand, run, and maintain?

All contributions to CATChem must include appropriate tests.

## Types of Tests

### 1. Unit Tests

- **Purpose**: To test individual components (subroutines, functions, or modules) in isolation.
- **Location**: `tests/unit/`

Unit tests are the foundation of our testing strategy. They should be fast, focused, and easy to run. We use a custom Fortran testing framework that provides a simple set of assertion procedures.

**Writing a Unit Test:**

```fortran
! In tests/unit/test_my_module.F90
program test_my_module
  use testing_mod
  use MyModule_Mod
  implicit none

  call test_suite_begin("MyModule Tests")
  call test_my_procedure()
  call test_suite_end()

contains

  subroutine test_my_procedure()
    use MyModule_Mod, only: my_procedure
    implicit none

    real :: expected_result, actual_result
    integer :: rc

    call test_begin("Test case for my_procedure")

    ! 1. Setup: Prepare the inputs
    expected_result = 42.0

    ! 2. Execute: Call the procedure to be tested
    call my_procedure(input_data, actual_result, rc)

    ! 3. Verify: Check the results
    call assert_success(rc, "Procedure should succeed")
    call assert_near(expected_result, actual_result, 1.0e-9, "Result should match expected value")

    call test_end()
  end subroutine

end program
```

### 2. Integration Tests

- **Purpose**: To test the interaction between different components of the model.
- **Location**: `tests/integration/`

Integration tests are used to verify that different parts of the model work together as expected. For example, an integration test might check the interaction between a process module and the state manager, or the full chain of processes from emissions to chemistry to deposition.

### 3. System Tests

- **Purpose**: To test the entire model in a realistic configuration.
- **Location**: `tests/system/`

System tests involve running the full CATChem model with a specific configuration and comparing the output to a set of reference results. These tests are crucial for validating the overall scientific behavior of the model.

### 4. Performance Tests

- **Purpose**: To benchmark the performance of the model and track it over time.
- **Location**: `tests/performance/`

Performance tests are used to identify performance regressions and to evaluate the impact of optimizations. See the [Performance Guide](performance.md) for more details.

## How to Write a Good Test

- **Isolate the Component**: A unit test should test one thing and one thing only. Avoid testing multiple procedures in a single test.
- **Be Descriptive**: Use descriptive names for your test suites and test cases. The name should make it clear what is being tested.
- **Write Clear Assertions**: The assertion message should clearly explain what is being checked and what the expected result is.
- **Test Edge Cases**: Don't just test the "happy path". Make sure you also test edge cases, boundary conditions, and invalid inputs.
- **Keep it Fast**: Unit tests should be fast. Slow tests are less likely to be run regularly.

## Running Tests

Tests are run using `ctest`, which is integrated with our CMake build system.

- **Run all tests**:
  ```bash
  cd build
  ctest
  ```
- **Run a specific test**:
  ```bash
  ctest -R my_test_name
  ```
- **Run tests with verbose output**:
  ```bash
  ctest --output-on-failure
  ```

## Testing and the Development Workflow

Testing is an integral part of the development workflow described in the [Contributing Guide](contributing.md).

1.  **Before you start coding**: Think about how you will test your changes.
2.  **As you code**: Write tests alongside your code. This is often referred to as Test-Driven Development (TDD).
3.  **Before you create a pull request**: Run all the tests to make sure your changes have not introduced any regressions.
4.  **In the pull request**: Your pull request will be automatically tested by our continuous integration (CI) system. All tests must pass before your pull request can be merged.

By following these guidelines, you can help us maintain the high quality and reliability of CATChem.
