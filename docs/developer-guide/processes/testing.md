# Testing Processes

Comprehensive testing is essential for reliable atmospheric process development in CATChem. This guide covers the testing framework, methodologies, and best practices for process validation.

## Testing Framework Overview

CATChem uses a multi-layered testing approach:

1. **Unit Tests**: Test individual functions and methods
2. **Process Tests**: Test complete process modules
3. **Integration Tests**: Test process interactions with the full model
4. **Validation Tests**: Compare against reference data and observations
5. **Performance Tests**: Measure computational efficiency

## Unit Testing Framework

### Testing Module

The `testing_mod` provides utilities for Fortran unit testing:

```fortran
use testing_mod

! Initialize test suite
call testing_init("Process Unit Tests")

! Start individual test
call testing_start_test("Test Description")

! Assertions
call assert_equal(actual, expected, "Values should match")
call assert_true(condition, "Condition should be true")
call assert_false(condition, "Condition should be false")
call assert_near(actual, expected, tolerance, "Values should be close")

! End test
call testing_end_test()

! Finalize test suite
call testing_finalize()
```

### Writing Unit Tests

Create unit tests for each process component:

```fortran
program test_settling_process
   use settlingProcess_Mod
   use StokesschemeScheme_Mod
   use testing_mod
   use state_mod

   implicit none

   call testing_init("Settling Process Tests")

   ! Test process initialization
   call test_process_initialization()

   ! Test Stokes scheme calculations
   call test_stokes_calculations()

   ! Test error handling
   call test_error_conditions()

   ! Test configuration validation
   call test_configuration_validation()

   call testing_finalize()

contains

   subroutine test_process_initialization()
      type(settlingProcessType) :: process
      type(StateContainerType) :: container
      integer :: rc

      call testing_start_test("Process Initialization")

      ! Create minimal test container
      call create_test_state_container(container)

      ! Test successful initialization
      call process%init(container, rc)
      call assert_equal(rc, CC_SUCCESS, "Init should succeed")
      call assert_true(process%is_ready(), "Process should be ready")
      call assert_equal(process%get_name(), 'settling', "Name should match")

      call testing_end_test()
   end subroutine test_process_initialization

   subroutine test_stokes_calculations()
      type(StateContainerType) :: container
      real(fp) :: particle_radius, air_density, temperature
      real(fp) :: settling_velocity, expected_velocity
      real(fp), parameter :: tolerance = 1.0e-6_fp
      integer :: rc

      call testing_start_test("Stokes Settling Calculations")

      ! Set up test conditions
      particle_radius = 1.0e-6_fp  ! 1 μm
      air_density = 1.2_fp         ! kg/m³
      temperature = 298.15_fp      ! K

      ! Calculate expected Stokes velocity (analytical solution)
      expected_velocity = calculate_analytical_stokes_velocity( &
         particle_radius, air_density, temperature)

      ! Test scheme calculation
      call create_test_state_container(container)
      call set_test_conditions(container, particle_radius, air_density, temperature)

      call stokes_calculate(container, rc)
      call assert_equal(rc, CC_SUCCESS, "Stokes calculation should succeed")

      ! Get calculated settling velocity
      settling_velocity = get_test_settling_velocity(container)

      ! Compare with analytical solution
      call assert_near(settling_velocity, expected_velocity, tolerance, &
                      "Settling velocity should match analytical solution")

      call testing_end_test()
   end subroutine test_stokes_calculations

   subroutine test_error_conditions()
      type(settlingProcessType) :: process
      type(StateContainerType) :: container
      integer :: rc

      call testing_start_test("Error Condition Handling")

      ! Test initialization with invalid configuration
      call create_invalid_test_container(container)
      call process%init(container, rc)
      call assert_not_equal(rc, CC_SUCCESS, "Init should fail with invalid config")

      ! Test run without initialization
      process%is_initialized = .false.
      call process%run(container, rc)
      call assert_not_equal(rc, CC_SUCCESS, "Run should fail without init")

      call testing_end_test()
   end subroutine test_error_conditions

end program test_settling_process
```

## Process Integration Testing

### Full Process Testing

Test processes within the complete model context:

```fortran
program test_settling_integration
   use catchem
   use testing_mod

   implicit none

   type(CATChemType) :: model
   character(len=256) :: config_file
   integer :: rc

   call testing_init("Settling Integration Tests")

   ! Test 1: Single timestep integration
   call test_single_timestep()

   ! Test 2: Multi-timestep stability
   call test_multi_timestep()

   ! Test 3: Mass conservation
   call test_mass_conservation()

   call testing_finalize()

contains

   subroutine test_single_timestep()
      call testing_start_test("Single Timestep Integration")

      config_file = "test_configs/settling_single_step.yml"

      ! Initialize model
      call model%init(config_file, rc)
      call assert_equal(rc, CC_SUCCESS, "Model init should succeed")

      ! Run single timestep
      call model%run(1, rc)
      call assert_equal(rc, CC_SUCCESS, "Single timestep should succeed")

      ! Check that settling occurred
      call verify_settling_occurred(model)

      call model%finalize(rc)

      call testing_end_test()
   end subroutine test_single_timestep

   subroutine test_mass_conservation()
      real(fp) :: initial_mass, final_mass
      real(fp), parameter :: conservation_tolerance = 1.0e-12_fp

      call testing_start_test("Mass Conservation")

      config_file = "test_configs/settling_conservation.yml"

      call model%init(config_file, rc)
      call assert_equal(rc, CC_SUCCESS, "Model init should succeed")

      ! Calculate initial total mass
      initial_mass = calculate_total_mass(model)

      ! Run multiple timesteps
      call model%run(100, rc)  ! 100 timesteps
      call assert_equal(rc, CC_SUCCESS, "Multi-timestep run should succeed")

      ! Calculate final total mass
      final_mass = calculate_total_mass(model)

      ! Check mass conservation (within numerical precision)
      call assert_near(final_mass, initial_mass, conservation_tolerance, &
                      "Total mass should be conserved")

      call model%finalize(rc)

      call testing_end_test()
   end subroutine test_mass_conservation

end program test_settling_integration
```

## Validation Testing

### Reference Data Comparison

Compare process outputs with established reference solutions:

```fortran
subroutine test_reference_validation()
   type(StateContainerType) :: container
   real(fp), allocatable :: reference_data(:,:,:)
   real(fp), allocatable :: computed_data(:,:,:)
   real(fp), parameter :: validation_tolerance = 0.05_fp  ! 5% tolerance
   integer :: rc

   call testing_start_test("Reference Data Validation")

   ! Load reference data
   call load_reference_data("reference/settling_test_case.nc", reference_data)

   ! Set up test case to match reference conditions
   call setup_reference_conditions(container)

   ! Run process
   call run_process_for_validation(container, rc)
   call assert_equal(rc, CC_SUCCESS, "Process should run successfully")

   ! Extract computed results
   computed_data = extract_computed_data(container)

   ! Compare with reference (statistical comparison)
   call validate_against_reference(computed_data, reference_data, &
                                  validation_tolerance)

   call testing_end_test()
end subroutine test_reference_validation
```

### Analytical Solutions

Test against known analytical solutions where available:

```fortran
subroutine test_analytical_validation()
   real(fp) :: particle_radius, particle_density, air_viscosity
   real(fp) :: analytical_velocity, computed_velocity
   real(fp), parameter :: analytical_tolerance = 1.0e-10_fp

   call testing_start_test("Analytical Solution Validation")

   ! Test conditions for which analytical solution exists
   particle_radius = 1.0e-6_fp    ! 1 μm
   particle_density = 2650.0_fp   ! kg/m³ (quartz)
   air_viscosity = 1.8e-5_fp      ! Pa·s

   ! Analytical Stokes settling velocity
   analytical_velocity = (2.0_fp * particle_radius**2 * particle_density * 9.81_fp) / &
                        (9.0_fp * air_viscosity)

   ! Compute using CATChem
   computed_velocity = compute_stokes_velocity(particle_radius, particle_density)

   ! Should match exactly for pure Stokes regime
   call assert_near(computed_velocity, analytical_velocity, analytical_tolerance, &
                   "Computed velocity should match analytical Stokes solution")

   call testing_end_test()
end subroutine test_analytical_validation
```

## Performance Testing

### Computational Benchmarking

Measure and validate computational performance:

```fortran
program benchmark_settling
   use settlingProcess_Mod
   use state_mod
   use iso_fortran_env, only : real64

   implicit none

   type(settlingProcessType) :: process
   type(StateContainerType) :: container
   integer, parameter :: n_runs = 1000
   integer :: i, rc
   real(real64) :: start_time, end_time, total_time

   ! Set up large test case
   call setup_large_test_case(container, 100, 100, 50)  ! 100x100x50 grid

   call process%init(container, rc)

   ! Warm up
   do i = 1, 10
      call process%run(container, rc)
   end do

   ! Benchmark
   call cpu_time(start_time)
   do i = 1, n_runs
      call process%run(container, rc)
   end do
   call cpu_time(end_time)

   total_time = end_time - start_time

   print *, "Settling process performance:"
   print *, "  Total time for", n_runs, "runs:", total_time, "seconds"
   print *, "  Average time per run:", total_time / n_runs, "seconds"
   print *, "  Throughput:", n_runs / total_time, "runs/second"

   call process%finalize(rc)

end program benchmark_settling
```

### Memory Usage Testing

Monitor memory usage and detect leaks:

```fortran
subroutine test_memory_usage()
   type(settlingProcessType) :: process
   type(StateContainerType) :: container
   integer :: initial_memory, final_memory, rc, i

   call testing_start_test("Memory Usage Testing")

   ! Get initial memory usage
   initial_memory = get_memory_usage()

   ! Initialize and run process multiple times
   call process%init(container, rc)
   do i = 1, 1000
      call process%run(container, rc)
   end do
   call process%finalize(rc)

   ! Get final memory usage
   final_memory = get_memory_usage()

   ! Check for memory leaks
   call assert_equal(final_memory, initial_memory, &
                    "Memory usage should return to initial level")

   call testing_end_test()
end subroutine test_memory_usage
```

## Test Configuration

### Test Data Generation

Create reproducible test datasets:

```yaml
# test_configs/settling_unit_test.yml
test_configuration:
  grid:
    nx: 10
    ny: 10
    nz: 20
    dx: 1000.0  # meters
    dy: 1000.0
    dz: 50.0

  meteorology:
    temperature: 298.15  # K
    pressure: 101325.0   # Pa
    air_density: 1.2     # kg/m³

  particles:
    - species: PM25
      radius: 1.25e-6   # m
      density: 1500.0    # kg/m³
      concentration: 10.0 # μg/m³

  settling:
    scheme: Stokes
    timestep: 60.0  # seconds
```

### Test Utilities

Common utilities for test setup:

```fortran
module test_utilities
   use precision_mod
   use state_mod

   implicit none
   private

   public :: create_test_state_container
   public :: setup_uniform_conditions
   public :: calculate_total_mass
   public :: validate_against_reference

contains

   subroutine create_test_state_container(container, nx, ny, nz)
      type(StateContainerType), intent(out) :: container
      integer, intent(in), optional :: nx, ny, nz

      integer :: local_nx, local_ny, local_nz

      local_nx = 10
      local_ny = 10
      local_nz = 20
      if (present(nx)) local_nx = nx
      if (present(ny)) local_ny = ny
      if (present(nz)) local_nz = nz

      ! Create minimal container for testing
      call container%init_for_testing(local_nx, local_ny, local_nz)

   end subroutine create_test_state_container

   subroutine setup_uniform_conditions(container, temp, press, density)
      type(StateContainerType), intent(inout) :: container
      real(fp), intent(in) :: temp, press, density

      type(MetStateType), pointer :: met_state

      met_state => container%get_met_state_ptr()
      call met_state%set_uniform('temperature', temp)
      call met_state%set_uniform('pressure', press)
      call met_state%set_uniform('air_density', density)

   end subroutine setup_uniform_conditions

end module test_utilities
```

## Test Organization

### Directory Structure

```
tests/
├── unit/                    # Unit tests
│   ├── test_settling.F90
│   ├── test_chemistry.F90
│   └── test_emissions.F90
├── integration/             # Integration tests
│   ├── test_full_model.F90
│   └── test_process_coupling.F90
├── validation/              # Validation tests
│   ├── test_reference_cases.F90
│   └── analytical_validation.F90
├── performance/             # Performance benchmarks
│   └── benchmark_processes.F90
├── data/                    # Test data
│   ├── reference/
│   └── validation/
└── configs/                 # Test configurations
    ├── unit_test.yml
    └── validation.yml
```

### Test Execution

Run tests using CMake/CTest:

```bash
# Build all tests
make build-tests

# Run all tests
make test

# Run specific test categories
ctest -L unit           # Unit tests only
ctest -L integration    # Integration tests only
ctest -L validation     # Validation tests only
ctest -L performance    # Performance tests only

# Run tests with verbose output
ctest -V

# Run specific test
ctest -R test_settling
```

## Continuous Integration

### Automated Testing

Set up CI pipelines to run tests automatically:

```yaml
# .github/workflows/tests.yml
name: CATChem Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-24.04
    strategy:
      matrix:
        compiler: ["12", "13", "14"]

    steps:
    - uses: actions/checkout@v6

    - name: Install dependencies
      run: |
        sudo apt-get install libnetcdf-dev libnetcdff-dev

    - name: Build
      run: |
        mkdir build && cd build
        cmake .. -DCMAKE_Fortran_COMPILER=gfortran-${{ matrix.compiler }}
        make -j2

    - name: Run tests
      run: |
        cd build
        ctest --output-on-failure
```

## Test Quality Assurance

### Coverage Analysis

Monitor test coverage:

```bash
# Build with coverage flags
cmake .. -DCMAKE_Fortran_FLAGS="--coverage"
make

# Run tests
make test

# Generate coverage report
gcov src/**/*.F90
lcov --capture --directory . --output-file coverage.info
genhtml coverage.info --output-directory coverage_html
```

### Test Metrics

Track important test metrics:

- **Code coverage**: Percentage of code exercised by tests
- **Test execution time**: Performance of test suite
- **Test reliability**: Frequency of false positives/negatives
- **Validation accuracy**: Agreement with reference data

## Best Practices

### Test Design

1. **Test pyramid**: Many unit tests, fewer integration tests, minimal end-to-end tests
2. **Fast execution**: Keep tests fast to enable frequent running
3. **Deterministic**: Tests should produce consistent results
4. **Independent**: Tests should not depend on each other
5. **Clear naming**: Test names should describe what is being tested

### Test Maintenance

1. **Keep tests simple**: Complex tests are hard to maintain and debug
2. **Update with code changes**: Modify tests when functionality changes
3. **Remove obsolete tests**: Delete tests for removed functionality
4. **Document test intent**: Explain what each test validates

### Error Handling

1. **Test error conditions**: Verify proper error handling
2. **Use meaningful assertions**: Make test failures informative
3. **Clean up resources**: Ensure tests don't leak memory or files
4. **Handle test failures gracefully**: Provide useful debugging information

The comprehensive testing framework ensures that CATChem processes are reliable, performant, and scientifically accurate across all supported platforms and configurations.
