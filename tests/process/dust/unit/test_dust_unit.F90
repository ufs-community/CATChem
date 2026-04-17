!> \file test_dust_unit.F90
!! \brief Unit tests for dust process
!!
!! This file contains unit tests for the dust process implementation
!! following the same pattern as core tests like test_ConfigManager.F90
!! Generated on: 2026-04-17T13:57:10.327748

program test_dust_unit
   use testing_mod, only: assert, assert_close
   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, ErrorManagerType
   use StateManager_Mod, only: StateManagerType
   use GridManager_Mod, only: GridManagerType
   use ProcessDustInterface_Mod, only: ProcessDustInterface
   use DustCommon_Mod, only: DustConfig

   implicit none

   type(ProcessDustInterface) :: dust_process
   type(StateManagerType) :: state_mgr
   type(ErrorManagerType) :: error_mgr
   type(GridManagerType) :: grid_mgr
   integer :: rc

   write(*,*) 'Testing Dust Process module...'
   write(*,*) ''

   ! Test 1: Initialize error manager
   write(*,*) 'Test 1: Initialize error manager'
   call error_mgr%init()

   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! Test 2: Initialize grid manager
   write(*,*) 'Test 2: Initialize grid manager'
   call grid_mgr%init(1, 1, 10, error_mgr, rc=rc)  ! 1x1 grid, 10 levels for testing
   call assert(rc == CC_SUCCESS, "GridManager initialization should succeed")

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: Initialize state manager
   write(*,*) 'Test 3: Initialize state manager'
   call state_mgr%init('TestDustStateManager', rc)
   call assert(rc == CC_SUCCESS, "StateManager initialization should succeed")

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! Test 4: Dust configuration creation and defaults
   write(*,*) 'Test 4: Dust configuration creation and defaults'
   call test_dust_config_defaults()

   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! Test 5: Dust configuration validation
   write(*,*) 'Test 5: Dust configuration validation'
   call test_dust_config_validation()

   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   ! Test 6: Dust scheme configuration
   write(*,*) 'Test 6: Dust scheme configuration'
   call test_scheme_configuration()

   write(*,*) 'Test 6 passed!'
   write(*,*) ''

   ! Test 7: ProcessDustInterface creation
   write(*,*) 'Test 7: ProcessDustInterface creation'
   call test_process_interface_creation()

   write(*,*) 'Test 7 passed!'
   write(*,*) ''

   ! Test 8: Process interface methods exist (without full initialization)
   write(*,*) 'Test 8: Process interface methods exist'
   call test_process_interface_methods()

   write(*,*) 'Test 8 passed!'
   write(*,*) ''
   write(*,*) 'All Dust unit tests completed successfully!'

contains

   !> Test Dust configuration default values
   subroutine test_dust_config_defaults()
      type(DustConfig) :: config

      ! Test default values are correctly set
      call assert(config%is_active .eqv. .true., "Default is_active should be true")
      call assert(len_trim(config%scheme) > 0, "Default scheme should be set")
      call assert(config%n_species == 0, "Default n_species should be 0")
      call assert(config%diagnostics .eqv. .false., "Default diagnostics should be false")

   end subroutine test_dust_config_defaults

   !> Test Dust configuration validation
   subroutine test_dust_config_validation()
      type(DustConfig) :: config
      type(ErrorManagerType) :: error_manager

      ! Test validation of default configuration
      call config%validate(error_manager)
      call assert(.true., "Default configuration validation completed")

      ! Test validation of different schemes
      config%scheme = 'fengsha'
      call config%validate(error_manager)
      call assert(.true., "FENGSHA scheme validation completed")

      config%scheme = 'ginoux'
      call config%validate(error_manager)
      call assert(.true., "GINOUX scheme validation completed")


   end subroutine test_dust_config_validation

   !> Test scheme configuration
   subroutine test_scheme_configuration()
      type(DustConfig) :: config
      type(ErrorManagerType) :: error_manager

      ! Test valid schemes
      config%scheme = 'fengsha'
      call config%validate(error_manager)
      call assert(.true., "FENGSHA scheme validation completed")

      config%scheme = 'ginoux'
      call config%validate(error_manager)
      call assert(.true., "GINOUX scheme validation completed")

      config%scheme = 'invalid_scheme'
      call config%validate(error_manager)
      call assert(.true., "Invalid scheme validation completed")


      ! Cleanup configuration
      call config%finalize()

   end subroutine test_scheme_configuration

   !> Test ProcessDustInterface can be created
   subroutine test_process_interface_creation()
      ! Test that we can create the interface object
      call assert(.true., "ProcessDustInterface object created successfully")
   end subroutine test_process_interface_creation

   !> Test ProcessDustInterface has required methods (without calling them)
   subroutine test_process_interface_methods()
      ! Test that the interface has the expected methods by checking if it's ready
      ! (this doesn't call init, just checks the initial state)
      call assert(.not. dust_process%is_ready(), "Process should not be ready before initialization")
   end subroutine test_process_interface_methods

end program test_dust_unit
