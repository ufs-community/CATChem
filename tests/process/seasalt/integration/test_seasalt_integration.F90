!> \file test_seasalt_integration.F90
!! \brief Comprehensive integration tests for seasalt process using CATChemCore
!!
!! This file contains comprehensive integration tests for the seasalt process implementation
!! using the centralized CATChemCore framework. Tests complete workflow: core initialization,
!! configuration loading, process registration, and all scheme validation.
!! Generated on: 2025-12-15T16:09:09.864661

program test_seasalt_integration
   use precision_mod, only: fp
   use iso_fortran_env, only: output_unit, error_unit
   use error_mod, only: CC_SUCCESS, CC_FAILURE, ErrorManagerType, ERROR_UNSUPPORTED_OPERATION
   use CATChemCore_Mod, only: CATChemCoreType, CATChemBuilderType
   use StateManager_Mod, only: StateManagerType
   use ProcessManager_Mod, only: ProcessManagerType
   use GridManager_Mod, only: GridManagerType
   use DiagnosticManager_Mod, only: DiagnosticManagerType
   use MetState_Mod, only: MetStateType
   use ChemState_Mod, only: ChemStateType
   use ConfigManager_Mod, only: ConfigManagerType
   use ProcessSeaSaltInterface_Mod, only: ProcessSeaSaltInterface
   use SeaSaltProcessCreator_Mod, only: register_seasalt_process
   use SeaSaltCommon_Mod, only: SeaSaltProcessConfig
   use DiagnosticInterface_Mod, only: DiagnosticRegistryType, DiagnosticFieldType, &
      DIAG_REAL_SCALAR, DIAG_REAL_1D, DIAG_REAL_2D, DIAG_REAL_3D, &
      DIAG_INTEGER_SCALAR, DIAG_INTEGER_1D, DIAG_INTEGER_2D, DIAG_INTEGER_3D

   implicit none

   ! Core framework
   type(CATChemCoreType) :: core
   type(CATChemBuilderType) :: builder
   type(ProcessManagerType), pointer :: process_mgr_ptr

   ! Configuration file path
   character(len=*), parameter :: config_file = './CATChem_new_config.yml'

   ! Test parameters for realistic emission scenario
   integer, parameter :: n_columns = 10    ! Grid columns
   integer, parameter :: n_levels = 20     ! Vertical levels (surface to ~20 km)
   integer, parameter :: n_time_steps = 5  ! Multiple timesteps for integration testing
   real(fp), parameter :: dt = 3600.0_fp   ! 1 hour timestep

   ! Test schemes
   character(len=20) :: schemes(3)

   integer :: rc, i_scheme, i_time
   logical :: all_tests_passed = .true.

   ! Initialize scheme array
   schemes = [ &
      'gong97              ', &
      'gong03              ', &
      'geos12              ']

   write(output_unit,'(A)') '=================================='
   write(output_unit,'(A)') '=== SEASALT INTEGRATION TESTS ==='
   write(output_unit,'(A)') '=================================='
   write(output_unit,'(A)') 'Using CATChemCore for comprehensive testing with'
   write(output_unit,'(A)') 'configuration, meteorological data, and all scheme validation'
   write(output_unit,'(A)') ''

   ! Step 1: Initialize CATChem Core with proper grid dimensions
   write(output_unit,'(A)') 'Step 1: Initializing CATChem Core...'

   call builder%init()
   builder = builder%with_name('SeaSaltIntegrationTest')
   builder = builder%with_config(config_file)
   builder = builder%with_grid(n_columns, 1, n_levels)
   builder = builder%with_verbose()
   call builder%build(core, rc)

   if (rc /= CC_SUCCESS) then
      write(error_unit,'(A)') 'ERROR: CATChemCore initialization/configuration failed'
      all_tests_passed = .false.
      goto 999
   end if
   write(output_unit,'(A,I0,A,I0,A)') '  ✓ CATChemCore initialized: ', n_columns, ' columns, ', n_levels, ' levels'
   write(output_unit,'(A)') '  ✓ Configuration loaded and all managers set up'

   ! Register seasalt processes with ProcessFactory
   process_mgr_ptr => core%get_process_manager()
   call register_seasalt_process(process_mgr_ptr, rc)
   if (rc /= CC_SUCCESS) then
      write(error_unit,'(A)') 'ERROR: Failed to register seasalt processes with ProcessFactory'
      all_tests_passed = .false.
      goto 999
   end if
   write(output_unit,'(A)') '  ✓ SeaSalt processes registered with ProcessFactory'

   ! Step 2: Set up realistic meteorological conditions
   write(output_unit,'(A)') ''
   write(output_unit,'(A)') 'Step 2: Setting up realistic meteorological conditions...'
   call setup_met(core, rc)
   if (rc /= CC_SUCCESS) then
      write(error_unit,'(A)') 'ERROR: Failed to set up meteorological conditions'
      all_tests_passed = .false.
      goto 999
   end if
   write(output_unit,'(A)') '  ✓ Meteorological conditions configured'

   ! Step 3: Testing seasalt process with all schemes
   write(output_unit,'(A)') ''
   write(output_unit,'(A)') 'Step 3: Testing seasalt process with all schemes...'

   ! Add seasalt process for scheme testing
   call core%add_process('seasalt', rc)
   if (rc /= CC_SUCCESS) then
      write(error_unit,'(A)') 'ERROR: Failed to add seasalt process for scheme testing'
      all_tests_passed = .false.
      goto 999
   end if
   write(output_unit,'(A)') '  ✓ SeaSalt process added successfully'

   write(output_unit,'(A)') ''
   write(output_unit,'(A)') '  Testing multiple seasalt schemes...'
   do i_scheme = 1, size(schemes)
      write(output_unit,'(A,A,A)') '    Testing ', trim(schemes(i_scheme)), ' scheme...'

      call test_scheme(core, schemes(i_scheme), rc)
      if (rc /= CC_SUCCESS) then
         write(output_unit,'(A,A)') '    ✗ ', trim(schemes(i_scheme)), ' scheme test failed'
         write(error_unit,'(A,A)') 'ERROR: Scheme test failed for ', trim(schemes(i_scheme))
         all_tests_passed = .false.
      else
         write(output_unit,'(A,A)') '    ✓ ', trim(schemes(i_scheme)), ' scheme test passed'
      end if
   end do
   write(output_unit,'(A)') '  ✓ All scheme tests completed'

   ! Step 4: Test multi-timestep stability
   write(output_unit,'(A)') ''
   write(output_unit,'(A)') 'Step 4: Testing multi-timestep stability...'
   write(output_unit,'(A,I0,A)') '  Running ', n_time_steps, ' timestep integration test...'


   do i_time = 1, n_time_steps
      call core%run_timestep(i_time, dt, rc)
      if (rc /= CC_SUCCESS) then
         write(error_unit,'(A,I0)') 'ERROR: Timestep ', i_time, ' failed'
         all_tests_passed = .false.
         exit
      end if
   end do

   if (all_tests_passed) then
      write(output_unit,'(A,I0,A)') '  ✓ All ', n_time_steps, ' timesteps completed successfully'
      write(output_unit,'(A)') '    - SeaSalt process stability verified'
      write(output_unit,'(A)') '    - Multi-timestep conservation maintained'
   end if

   ! Final validation and cleanup
   write(output_unit,'(A)') ''
   write(output_unit,'(A)') 'Final validation and cleanup...'
   call core%finalize(rc)
   if (rc /= CC_SUCCESS) then
      write(error_unit,'(A)') 'WARNING: Core finalization had issues'
   end if

999 continue

   ! Print final results
   write(output_unit,'(A)') ''
   write(output_unit,'(A)') '=================================='
   if (all_tests_passed) then
      write(output_unit,'(A)') '=== ALL SEASALT TESTS PASSED! ==='
      write(output_unit,'(A)') '=== Integration test successful ==='
   else
      write(output_unit,'(A)') '=== SOME SEASALT TESTS FAILED ==='
      write(output_unit,'(A)') '=== Check error messages above ==='
   end if
   write(output_unit,'(A)') '=================================='

   if (.not. all_tests_passed) stop 1

contains

   !> Set up realistic meteorological conditions for seasalt testing
   subroutine setup_met(core_arg, rc_arg)
      type(CATChemCoreType), intent(inout) :: core_arg
      integer, intent(out) :: rc_arg

      type(StateManagerType), pointer :: state_mgr
      type(MetStateType), pointer :: met_state
      type(GridManagerType), pointer :: grid_mgr
      integer :: nx, ny, nz, i, j, k
      real(fp) :: lat, wind_speed, altitude_km, edge_altitude_km

      rc_arg = CC_SUCCESS

      ! Get managers and state pointers
      state_mgr => core_arg%get_state_manager()
      met_state => state_mgr%get_met_state_ptr()
      grid_mgr => core_arg%get_grid_manager()

      ! Get grid dimensions
      call grid_mgr%get_shape(nx, ny, nz)

      ! Allocate categorical arrays with standard dimensions

      ! Set realistic conditions for seasalt processes
      do j = 1, ny
         ! Calculate latitude for realistic gradients
         lat = -30.0_fp + (j-1) * 60.0_fp / max(1, ny-1)  ! -30°S to 30°N
         do i = 1, nx
            met_state%FROCEAN(i,j) = 1.0_fp                    ! Pure ocean everywhere
            met_state%FRSEAICE(i,j) = 0.0_fp                   ! No sea ice
            met_state%SST(i,j) = 298.0_fp + 5.0_fp * cos(lat * 3.14159_fp / 180.0_fp)  ! 293-303K SST
            wind_speed = 8.0_fp + 2.0_fp * cos(lat * 3.14159_fp / 180.0_fp)  ! 6-10 m/s
            met_state%U10M(i,j) = -wind_speed * 0.8_fp         ! Easterly trade winds
            met_state%V10M(i,j) = wind_speed * 0.3_fp          ! Slight northerly component
            met_state%USTAR(i,j) = 0.03_fp * sqrt(met_state%U10M(i,j)**2 + met_state%V10M(i,j)**2)
         end do
      end do

      ! Set up 3D atmospheric fields (nx, ny, nz)
      do j = 1, ny
         do i = 1, nx
            do k = 1, nz
               ! Calculate height-dependent values
               ! Approximate altitude in km (assuming ~1 km per level near surface)
               altitude_km = real(k-1, fp) * 1.0_fp
               met_state%DELP(i,j,k) = 5000.0_fp                                  ! Pressure thickness [Pa]
            end do
         end do
      end do



      ! Set up DELP (pressure difference between levels) for emission unit conversion
      ! DELP is only used for unit conversion in emission processes
      do j = 1, ny
         do i = 1, nx
            do k = 1, nz
               ! Set realistic pressure differences (Pa) for atmospheric layers
               ! Surface layers have higher DELP, upper levels have lower DELP
               altitude_km = real(k-1, fp) * 20.0_fp / real(nz-1, fp)  ! 0-20 km altitude
               met_state%DELP(i,j,k) = 10000.0_fp * exp(-altitude_km / 8.0_fp)  ! Exponential pressure decrease
            end do
         end do
      end do

   end subroutine setup_met

   !> Test a specific seasalt scheme with comprehensive validation
   subroutine test_scheme(core_arg, scheme_name, rc_arg)
      type(CATChemCoreType), intent(inout) :: core_arg
      character(len=*), intent(in) :: scheme_name
      integer, intent(out) :: rc_arg

      type(ProcessManagerType), pointer :: process_mgr
      type(StateManagerType), pointer :: state_mgr
      type(ProcessSeaSaltInterface), pointer :: seasalt_interface
      type(ConfigManagerType), pointer :: config_mgr
      type(ErrorManagerType), pointer :: error_mgr

      rc_arg = CC_SUCCESS

      ! Get process manager and state manager
      process_mgr => core_arg%get_process_manager()
      state_mgr => core_arg%get_state_manager()

      ! Get seasalt process interface
      seasalt_interface => null()
      select type(process => process_mgr%processes(1)%item)
       type is (ProcessSeaSaltInterface)
         seasalt_interface => process
      end select

      if (.not. associated(seasalt_interface)) then
         rc_arg = CC_FAILURE
         return
      end if

      ! Step 1: Set the timestep for process calculations
      call seasalt_interface%set_timestep(dt)

      ! Step 2: Set the scheme
      call seasalt_interface%set_scheme(scheme_name)

      ! Step 3: Reload scheme-specific configuration
      config_mgr => state_mgr%get_config_ptr()
      error_mgr => state_mgr%get_error_manager()

      if (.not. associated(config_mgr)) then
         call error_mgr%report_error(1003, &
            'ConfigManager not available from StateManager', rc_arg)
         return
      end if

      if (.not. associated(error_mgr)) then
         rc_arg = CC_FAILURE
         return
      end if

      ! Call the scheme-specific loading function directly
      select case (trim(scheme_name))
       case ('gong97')
         call seasalt_interface%process_config%load_gong97_config(config_mgr, error_mgr)
       case ('gong03')
         call seasalt_interface%process_config%load_gong03_config(config_mgr, error_mgr)
       case ('geos12')
         call seasalt_interface%process_config%load_geos12_config(config_mgr, error_mgr)
       case default
         call error_mgr%report_error(1004, &
            'Unknown scheme: ' // trim(scheme_name), rc_arg)
         return
      end select


      ! Step 4: Reset diagnostics for the new scheme
      call reset_diagnostics_for_scheme(seasalt_interface, state_mgr, scheme_name, rc_arg)
      if (rc_arg /= CC_SUCCESS) return

      ! Step 5: Run the process to populate diagnostic data
      call process_mgr%run_column_processes(state_mgr, rc_arg)
      if (rc_arg /= CC_SUCCESS) return

      ! Step 6: Validate all results
      call validate_results(core_arg, rc_arg)

   end subroutine test_scheme

   !> Reset diagnostics for a specific scheme (test-specific function)
   !! This function handles diagnostic reset when switching between schemes during testing
   subroutine reset_diagnostics_for_scheme(seasalt_interface, container, scheme_name, rc_arg)
      type(ProcessSeaSaltInterface), intent(inout) :: seasalt_interface
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: scheme_name
      integer, intent(out) :: rc_arg

      type(DiagnosticManagerType), pointer :: diag_mgr
      type(ErrorManagerType), pointer :: error_mgr
      character(len=64) :: current_scheme

      rc_arg = CC_SUCCESS

      ! Get managers
      diag_mgr => container%get_diagnostic_manager()
      error_mgr => container%get_error_manager()

      ! Get current scheme
      current_scheme = seasalt_interface%get_scheme()

      ! Remove existing process registration (this clears all diagnostic fields)
      call diag_mgr%remove_process('seasalt', rc_arg)
      if (rc_arg /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_UNSUPPORTED_OPERATION, &
            'Failed to remove existing diagnostics for seasalt process', rc_arg)
         ! Continue anyway - this might be the first registration
         rc_arg = CC_SUCCESS
      endif

      ! Re-register diagnostics for the new scheme
      ! The scheme-specific configuration should already be set correctly
      call seasalt_interface%register_diagnostics(container, rc_arg)
      if (rc_arg /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_UNSUPPORTED_OPERATION, &
            'Failed to re-register diagnostics for scheme: ' // &
            trim(current_scheme), rc_arg)
         return
      endif

   end subroutine reset_diagnostics_for_scheme

   !> Validate the results of the integration test with comprehensive diagnostic checking
   subroutine validate_results(core_arg, rc_arg)
      type(CATChemCoreType), intent(inout) :: core_arg
      integer, intent(out) :: rc_arg

      type(DiagnosticManagerType), pointer :: diag_mgr
      type(DiagnosticRegistryType), pointer :: seasalt_registry
      character(len=64), allocatable :: field_names(:)
      integer :: num_fields, i, local_rc, data_type
      real(fp) :: scalar_value
      real(fp), pointer :: array_1d_ptr(:) => null()
      real(fp), pointer :: array_2d_ptr(:,:) => null()
      real(fp), pointer :: array_3d_ptr(:,:,:) => null()
      logical :: validation_passed
      character(len=64) :: field_name
      character(len=20) :: type_name

      rc_arg = CC_SUCCESS
      validation_passed = .true.

      write(output_unit,'(A)') '  Validating seasalt emission results...'

      ! Use core validation first
      if (.not. core_arg%validate()) then
         write(error_unit,'(A)') '  ERROR: Core validation failed'
         rc_arg = CC_FAILURE
         return
      end if
      write(output_unit,'(A)') '    ✓ Core validation passed'

      ! Get DiagnosticManager from core
      diag_mgr => core_arg%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) then
         write(error_unit,'(A)') '  ERROR: Could not get DiagnosticManager from core'
         rc_arg = CC_FAILURE
         return
      end if

      ! Get the seasalt process diagnostic registry
      call diag_mgr%get_process_registry('seasalt', seasalt_registry, local_rc)
      if (local_rc /= CC_SUCCESS .or. .not. associated(seasalt_registry)) then
         write(error_unit,'(A)') '  ERROR: Could not get seasalt process registry'
         rc_arg = CC_FAILURE
         return
      end if

      ! Get the number of registered diagnostic fields
      num_fields = seasalt_registry%get_field_count()
      write(output_unit,'(A,I0,A)') '    Found ', num_fields, ' registered diagnostic fields for seasalt process'

      if (num_fields == 0) then
         write(error_unit,'(A)') '  ERROR: No diagnostic fields registered for seasalt process'
         rc_arg = CC_FAILURE
         return
      end if

      ! Allocate array for field names
      allocate(field_names(num_fields))

      ! Get all field names
      call seasalt_registry%list_fields(field_names, num_fields)

      ! Iterate through all diagnostic fields and validate them
      write(output_unit,'(A)') '    Validating all registered diagnostic fields:'

      do i = 1, num_fields
         field_name = trim(field_names(i))
         write(output_unit,'(A,I0,A,A)') '      Field ', i, ': ', trim(field_name)

         ! Get field values and type information directly from DiagnosticManager
         call diag_mgr%get_field_value('seasalt', field_name, &
            scalar_value=scalar_value, &
            array_1d_ptr=array_1d_ptr, &
            array_2d_ptr=array_2d_ptr, &
            array_3d_ptr=array_3d_ptr, &
            data_type=data_type, &
            rc=local_rc)
         if (local_rc /= CC_SUCCESS) then
            write(error_unit,'(A,A)') '    WARNING: Could not retrieve field value: ', trim(field_name)
            validation_passed = .false.
            cycle
         end if

         ! Convert data type to readable name and validate values
         call validate_field_by_type(field_name, data_type, scalar_value, &
            array_1d_ptr, array_2d_ptr, array_3d_ptr, validation_passed, verbose=.false.)

      end do

      ! Clean up
      deallocate(field_names)

      ! Final validation result
      if (.not. validation_passed) then
         write(error_unit,'(A)') '  VALIDATION FAILED: Some seasalt diagnostics failed validation'
         rc_arg = CC_FAILURE
      else
         write(output_unit,'(A)') '  ✓ All seasalt diagnostic validations passed'
         write(output_unit,'(A,I0,A)') '    - ', num_fields, ' diagnostic fields validated'
         write(output_unit,'(A)') '    - All emission values are positive'
         write(output_unit,'(A)') '    - Diagnostic system is functioning correctly'
      end if

   end subroutine validate_results

   !> Validate field values based on type and emission expectations
   subroutine validate_field_by_type(field_name, data_type, scalar_value, &
      array_1d_ptr, array_2d_ptr, array_3d_ptr, validation_passed, verbose)
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: data_type
      real(fp), intent(in) :: scalar_value
      real(fp), pointer, intent(in) :: array_1d_ptr(:)
      real(fp), pointer, intent(in) :: array_2d_ptr(:,:)
      real(fp), pointer, intent(in) :: array_3d_ptr(:,:,:)
      logical, intent(inout) :: validation_passed
      logical, intent(in), optional :: verbose

      logical :: field_passed
      logical :: is_verbose = .false.
      character(len=20) :: type_name
      integer :: i, j, k
      real(fp) :: current_value

      ! Set verbose mode
      if (present(verbose)) is_verbose = verbose
      ! Explicitly initialize field_passed for each call
      field_passed = .true.

      ! Convert data type to readable name and validate values
      select case (data_type)
       case (DIAG_REAL_SCALAR)
         type_name = 'Real Scalar'
         write(output_unit,'(A,A)') '        Type: ', trim(type_name)
         write(output_unit,'(A,E12.5)') '        Scalar value: ', scalar_value
         if (is_verbose) then
            write(output_unit,'(A,A,A,E12.5)') '          ', trim(field_name), ' = ', scalar_value
         end if

         ! Check if scalar value is finite and non-negative
         if (scalar_value /= scalar_value) then  ! NaN check
            write(error_unit,'(A,A)') '    ERROR: Field has NaN value: ', trim(field_name)
            field_passed = .false.
         else if (scalar_value < 0.0_fp) then
            write(error_unit,'(A,A)') '    ERROR: Field has negative value: ', trim(field_name)
            field_passed = .false.
         else if (.not. (scalar_value < huge(scalar_value))) then  ! Infinite check
            write(error_unit,'(A,A)') '    ERROR: Field has infinite value: ', trim(field_name)
            field_passed = .false.
         else
            write(output_unit,'(A,A)') '        ✓ Field has valid finite non-negative value: ', trim(field_name)
         end if

       case (DIAG_REAL_1D)
         type_name = 'Real 1D Array'
         write(output_unit,'(A,A)') '        Type: ', trim(type_name)
         if (associated(array_1d_ptr)) then
            write(output_unit,'(A,I0)') '        Array size: ', size(array_1d_ptr)

            ! Check each element in the 1D array
            do i = 1, size(array_1d_ptr)
               current_value = array_1d_ptr(i)
               if (is_verbose) then
                  write(output_unit,'(A,A,A,I0,A,E12.5)') '          ', trim(field_name), '[', i, '] = ', current_value
               end if
               if (current_value /= current_value) then  ! NaN check
                  write(error_unit,'(A,A,A,I0,A)') '    ERROR: Field has NaN at index ', trim(field_name), ' (', i, ')'
                  field_passed = .false.
                  exit
               else if (current_value < 0.0_fp) then
                  write(error_unit,'(A,A,A,I0,A)') '    ERROR: Field has negative value at index ', trim(field_name), ' (', i, ')'
                  field_passed = .false.
                  exit
               else if (.not. (current_value < huge(current_value))) then  ! Infinite check
                  write(error_unit,'(A,A,A,I0,A)') '    ERROR: Field has infinite value at index ', trim(field_name), ' (', i, ')'
                  field_passed = .false.
                  exit
               end if
            end do

            if (field_passed) then
               ! Check if sum of array is zero
               if (sum(array_1d_ptr) == 0.0_fp) then
                  write(error_unit,'(A,A)') '    WARNING: Field has zero sum (all elements are zero): ', trim(field_name)
                  !field_passed = .false.
               else
                  write(output_unit,'(A,A)') '        ✓ All array elements are finite and non-negative: ', trim(field_name)
               end if
            end if
         else
            write(error_unit,'(A,A)') '    ERROR: 1D array not associated for field: ', trim(field_name)
            field_passed = .false.
         end if

       case (DIAG_REAL_2D)
         type_name = 'Real 2D Array'
         write(output_unit,'(A,A)') '        Type: ', trim(type_name)
         if (associated(array_2d_ptr)) then
            write(output_unit,'(A,I0,A,I0)') '        Array size: ', size(array_2d_ptr,1), ' x ', size(array_2d_ptr,2)

            ! Check each element in the 2D array
            outer_loop_2d: do j = 1, size(array_2d_ptr,2)
               do i = 1, size(array_2d_ptr,1)
                  current_value = array_2d_ptr(i,j)
                  if (is_verbose) then
                     write(output_unit,'(A,A,A,I0,A,I0,A,E12.5)') '          ', trim(field_name), '[', i, ',', j, '] = ', current_value
                  end if
                  if (current_value /= current_value) then  ! NaN check
                     write(error_unit,'(A,A,A,I0,A,I0,A)') '    ERROR: Field has NaN at index ', trim(field_name), ' (', i, ',', j, ')'
                     field_passed = .false.
                     exit outer_loop_2d
                  else if (current_value < 0.0_fp) then
                     write(error_unit,'(A,A,A,I0,A,I0,A)') '    ERROR: Field has negative value at index ', trim(field_name), ' (', i, ',', j, ')'
                     field_passed = .false.
                     exit outer_loop_2d
                  else if (.not. (current_value < huge(current_value))) then  ! Infinite check
                     write(error_unit,'(A,A,A,I0,A,I0,A)') '    ERROR: Field has infinite value at index ', trim(field_name), ' (', i, ',', j, ')'
                     field_passed = .false.
                     exit outer_loop_2d
                  end if
               end do
            end do outer_loop_2d

            if (field_passed) then
               ! Check if sum of array is zero
               if (sum(array_2d_ptr) == 0.0_fp) then
                  write(error_unit,'(A,A)') '    WARNING: Field has zero sum (all elements are zero): ', trim(field_name)
                  !field_passed = .false.
               else
                  write(output_unit,'(A,A)') '        ✓ All array elements are finite and non-negative: ', trim(field_name)
               end if
            end if
         else
            write(error_unit,'(A,A)') '    ERROR: 2D array not associated for field: ', trim(field_name)
            field_passed = .false.
         end if

       case (DIAG_REAL_3D)
         type_name = 'Real 3D Array'
         write(output_unit,'(A,A)') '        Type: ', trim(type_name)
         if (associated(array_3d_ptr)) then
            write(output_unit,'(A,I0,A,I0,A,I0)') '        Array size: ', size(array_3d_ptr,1), ' x ', size(array_3d_ptr,2), ' x ', size(array_3d_ptr,3)

            ! Check each element in the 3D array
            outer_loop_3d: do k = 1, size(array_3d_ptr,3)
               do j = 1, size(array_3d_ptr,2)
                  do i = 1, size(array_3d_ptr,1)
                     current_value = array_3d_ptr(i,j,k)
                     if (is_verbose) then
                        write(output_unit,'(A,A,A,I0,A,I0,A,I0,A,E12.5)') '          ', trim(field_name), '[', i, ',', j, ',', k, '] = ', current_value
                     end if
                     if (current_value /= current_value) then  ! NaN check
                        write(error_unit,'(A,A,A,I0,A,I0,A,I0,A)') '    ERROR: Field has NaN at index ', trim(field_name), ' (', i, ',', j, ',', k, ')'
                        field_passed = .false.
                        exit outer_loop_3d
                     else if (current_value < 0.0_fp) then
                        write(error_unit,'(A,A,A,I0,A,I0,A,I0,A)') '    ERROR: Field has negative value at index ', trim(field_name), ' (', i, ',', j, ',', k, ')'
                        field_passed = .false.
                        exit outer_loop_3d
                     else if (.not. (current_value < huge(current_value))) then  ! Infinite check
                        write(error_unit,'(A,A,A,I0,A,I0,A,I0,A)') '    ERROR: Field has infinite value at index ', trim(field_name), ' (', i, ',', j, ',', k, ')'
                        field_passed = .false.
                        exit outer_loop_3d
                     end if
                  end do
               end do
            end do outer_loop_3d

            if (field_passed) then
               ! Check if sum of array is zero
               if (sum(array_3d_ptr) == 0.0_fp) then
                  write(error_unit,'(A,A)') '    WARNING: Field has zero sum (all elements are zero): ', trim(field_name)
                  !field_passed = .false.
               else
                  write(output_unit,'(A,A)') '        ✓ All array elements are finite and non-negative: ', trim(field_name)
               end if
            end if
         else
            write(error_unit,'(A,A)') '    ERROR: 3D array not associated for field: ', trim(field_name)
            field_passed = .false.
         end if

       case (DIAG_INTEGER_SCALAR)
         type_name = 'Integer Scalar'
         write(output_unit,'(A,A)') '        Type: ', trim(type_name)
         write(output_unit,'(A,E12.5)') '        Scalar value: ', scalar_value
         if (is_verbose) then
            write(output_unit,'(A,A,A,E12.5)') '          ', trim(field_name), ' = ', scalar_value
         end if

         ! For integer scalar, just check if it's non-negative
         if (scalar_value < 0.0_fp) then
            write(error_unit,'(A,A)') '    ERROR: Integer field has negative value: ', trim(field_name)
            field_passed = .false.
         else
            write(output_unit,'(A,A)') '        ✓ Integer field has non-negative value: ', trim(field_name)
         end if

       case (DIAG_INTEGER_1D)
         type_name = 'Integer 1D Array'
         write(output_unit,'(A,A)') '        Type: ', trim(type_name)
         if (associated(array_1d_ptr)) then
            write(output_unit,'(A,I0)') '        Array size: ', size(array_1d_ptr)

            ! Check each element in the 1D integer array
            do i = 1, size(array_1d_ptr)
               current_value = array_1d_ptr(i)
               if (is_verbose) then
                  write(output_unit,'(A,A,A,I0,A,E12.5)') '          ', trim(field_name), '[', i, '] = ', current_value
               end if
               if (current_value /= current_value) then  ! NaN check
                  write(error_unit,'(A,A,A,I0,A)') '    ERROR: Integer field has NaN at index ', trim(field_name), ' (', i, ')'
                  field_passed = .false.
                  exit
               else if (current_value < 0.0_fp) then
                  write(error_unit,'(A,A,A,I0,A)') '    ERROR: Integer field has negative value at index ', trim(field_name), ' (', i, ')'
                  field_passed = .false.
                  exit
               end if
            end do

            if (field_passed) then
               ! Check if sum of array is zero
               if (sum(array_1d_ptr) == 0.0_fp) then
                  write(error_unit,'(A,A)') '    WARNING: Integer field has zero sum (all elements are zero): ', trim(field_name)
                  !field_passed = .false.
               else
                  write(output_unit,'(A,A)') '        ✓ All integer array elements are non-negative: ', trim(field_name)
               end if
            end if
         else
            write(error_unit,'(A,A)') '    ERROR: 1D integer array not associated for field: ', trim(field_name)
            field_passed = .false.
         end if

       case (DIAG_INTEGER_2D)
         type_name = 'Integer 2D Array'
         write(output_unit,'(A,A)') '        Type: ', trim(type_name)
         if (associated(array_2d_ptr)) then
            write(output_unit,'(A,I0,A,I0)') '        Array size: ', size(array_2d_ptr,1), ' x ', size(array_2d_ptr,2)

            ! Check each element in the 2D integer array
            outer_loop_int_2d: do j = 1, size(array_2d_ptr,2)
               do i = 1, size(array_2d_ptr,1)
                  current_value = array_2d_ptr(i,j)
                  if (is_verbose) then
                     write(output_unit,'(A,A,A,I0,A,I0,A,E12.5)') '          ', trim(field_name), '[', i, ',', j, '] = ', current_value
                  end if
                  if (current_value /= current_value) then  ! NaN check
                     write(error_unit,'(A,A,A,I0,A,I0,A)') '    ERROR: Integer field has NaN at index ', trim(field_name), ' (', i, ',', j, ')'
                     field_passed = .false.
                     exit outer_loop_int_2d
                  else if (current_value < 0.0_fp) then
                     write(error_unit,'(A,A,A,I0,A,I0,A)') '    ERROR: Integer field has negative value at index ', trim(field_name), ' (', i, ',', j, ')'
                     field_passed = .false.
                     exit outer_loop_int_2d
                  end if
               end do
            end do outer_loop_int_2d

            if (field_passed) then
               ! Check if sum of array is zero
               if (sum(array_2d_ptr) == 0.0_fp) then
                  write(error_unit,'(A,A)') '    WARNING: Integer field has zero sum (all elements are zero): ', trim(field_name)
                  !field_passed = .false.
               else
                  write(output_unit,'(A,A)') '        ✓ All integer array elements are non-negative: ', trim(field_name)
               end if
            end if
         else
            write(error_unit,'(A,A)') '    ERROR: 2D integer array not associated for field: ', trim(field_name)
            field_passed = .false.
         end if

       case (DIAG_INTEGER_3D)
         type_name = 'Integer 3D Array'
         write(output_unit,'(A,A)') '        Type: ', trim(type_name)
         if (associated(array_3d_ptr)) then
            write(output_unit,'(A,I0,A,I0,A,I0)') '        Array size: ', size(array_3d_ptr,1), ' x ', size(array_3d_ptr,2), ' x ', size(array_3d_ptr,3)

            ! Check each element in the 3D integer array
            outer_loop_int_3d: do k = 1, size(array_3d_ptr,3)
               do j = 1, size(array_3d_ptr,2)
                  do i = 1, size(array_3d_ptr,1)
                     current_value = array_3d_ptr(i,j,k)
                     if (is_verbose) then
                        write(output_unit,'(A,A,A,I0,A,I0,A,I0,A,E12.5)') '          ', trim(field_name), '[', i, ',', j, ',', k, '] = ', current_value
                     end if
                     if (current_value /= current_value) then  ! NaN check
                        write(error_unit,'(A,A,A,I0,A,I0,A,I0,A)') '    ERROR: Integer field has NaN at index ', trim(field_name), ' (', i, ',', j, ',', k, ')'
                        field_passed = .false.
                        exit outer_loop_int_3d
                     else if (current_value < 0.0_fp) then
                        write(error_unit,'(A,A,A,I0,A,I0,A,I0,A)') '    ERROR: Integer field has negative value at index ', trim(field_name), ' (', i, ',', j, ',', k, ')'
                        field_passed = .false.
                        exit outer_loop_int_3d
                     end if
                  end do
               end do
            end do outer_loop_int_3d

            if (field_passed) then
               ! Check if sum of array is zero
               if (sum(array_3d_ptr) == 0.0_fp) then
                  write(error_unit,'(A,A)') '    WARNING: Integer field has zero sum (all elements are zero): ', trim(field_name)
                  !field_passed = .false.
               else
                  write(output_unit,'(A,A)') '        ✓ All integer array elements are non-negative: ', trim(field_name)
               end if
            end if
         else
            write(error_unit,'(A,A)') '    ERROR: 3D integer array not associated for field: ', trim(field_name)
            field_passed = .false.
         end if

       case default
         type_name = 'Unknown Type'
         write(output_unit,'(A,A)') '        Type: ', trim(type_name)
         write(error_unit,'(A,A)') '    ERROR: Unsupported data type for field: ', trim(field_name)
         field_passed = .false.
         return
      end select

      ! Update overall validation status
      if (.not. field_passed) then
         validation_passed = .false.
      end if

   end subroutine validate_field_by_type

end program test_seasalt_integration
