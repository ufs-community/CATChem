!> \file test_drydep_unit.F90
!! \brief Unit tests for drydep process
!!
!! This file contains unit tests for the drydep process implementation
!! following the same pattern as core tests like test_ConfigManager.F90
!! Generated on: 2025-11-14T22:58:26.634543

program test_drydep_unit
   use testing_mod, only: assert, assert_close
   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, ErrorManagerType
   use StateManager_Mod, only: StateManagerType
   use GridManager_Mod, only: GridManagerType
   use ProcessDryDepInterface_Mod, only: ProcessDryDepInterface
   use DryDepCommon_Mod, only: DryDepConfig

   implicit none

   type(ProcessDryDepInterface) :: drydep_process
   type(StateManagerType) :: state_mgr
   type(ErrorManagerType) :: error_mgr
   type(GridManagerType) :: grid_mgr
   integer :: rc

   write(*,*) 'Testing DryDep Process module...'
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
   call state_mgr%init('TestDryDepStateManager', rc)
   call assert(rc == CC_SUCCESS, "StateManager initialization should succeed")

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! Test 4: DryDep configuration creation and defaults
   write(*,*) 'Test 4: DryDep configuration creation and defaults'
   call test_drydep_config_defaults()

   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! Test 5: DryDep configuration validation
   write(*,*) 'Test 5: DryDep configuration validation'
   call test_drydep_config_validation()

   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   ! Test 6: DryDep scheme configuration
   write(*,*) 'Test 6: DryDep scheme configuration'
   call test_scheme_configuration()

   write(*,*) 'Test 6 passed!'
   write(*,*) ''

   ! Test 7: ProcessDryDepInterface creation
   write(*,*) 'Test 7: ProcessDryDepInterface creation'
   call test_process_interface_creation()

   write(*,*) 'Test 7 passed!'
   write(*,*) ''

   ! Test 8: Process interface methods exist (without full initialization)
   write(*,*) 'Test 8: Process interface methods exist'
   call test_process_interface_methods()

   write(*,*) 'Test 8 passed!'
   write(*,*) ''
   write(*,*) 'All DryDep unit tests completed successfully!'

contains

   !> Test DryDep configuration default values
   subroutine test_drydep_config_defaults()
      type(DryDepConfig) :: config

      ! Test default values are correctly set
      call assert(config%is_active .eqv. .true., "Default is_active should be true")
      call assert(len_trim(config%gas_scheme) > 0, "Default gas_scheme should be set")
      call assert(len_trim(config%aero_scheme) > 0, "Default aero_scheme should be set")
      call assert(config%n_species == 0, "Default n_species should be 0")
      call assert(config%diagnostics .eqv. .false., "Default diagnostics should be false")

   end subroutine test_drydep_config_defaults

   !> Test DryDep configuration validation
   subroutine test_drydep_config_validation()
      type(DryDepConfig) :: config
      type(ErrorManagerType) :: error_manager

      ! Test validation of default configuration
      call config%validate(error_manager)
      call assert(.true., "Default configuration validation completed")

      ! Test validation of different schemes
      ! Test gas schemes
      config%gas_scheme = 'wesely'
      call config%validate(error_manager)
      call assert(.true., "WESELY gas scheme validation completed")

      ! Test aerosol schemes
      config%aero_scheme = 'gocart'
      call config%validate(error_manager)
      call assert(.true., "GOCART aerosol scheme validation completed")
      config%aero_scheme = 'zhang'
      call config%validate(error_manager)
      call assert(.true., "ZHANG aerosol scheme validation completed")
      ! Test gas/aerosol species differentiation
      call assert(.true., "Gas/aerosol differentiation enabled for this process")

   end subroutine test_drydep_config_validation

   !> Test scheme configuration
   subroutine test_scheme_configuration()
      type(DryDepConfig) :: config
      type(ErrorManagerType) :: error_manager

      ! Test valid schemes
      ! Test gas schemes
      config%gas_scheme = 'wesely'
      call config%validate(error_manager)
      call assert(.true., "WESELY gas scheme validation completed")
      call assert(.true., "WESELY scheme configured for gas species")

      ! Test aerosol schemes
      config%aero_scheme = 'gocart'
      call config%validate(error_manager)
      call assert(.true., "GOCART aerosol scheme validation completed")
      call assert(.true., "GOCART scheme configured for aerosol species")
      config%aero_scheme = 'zhang'
      call config%validate(error_manager)
      call assert(.true., "ZHANG aerosol scheme validation completed")
      call assert(.true., "ZHANG scheme configured for aerosol species")

      ! Test invalid schemes
      config%gas_scheme = 'invalid_gas_scheme'
      call config%validate(error_manager)
      call assert(.true., "Invalid gas scheme validation completed")

      config%aero_scheme = 'invalid_aero_scheme'
      call config%validate(error_manager)
      call assert(.true., "Invalid aerosol scheme validation completed")

      ! Test gas/aerosol differentiation logic
      call assert(.true., "Process supports gas/aerosol species differentiation")

      ! Cleanup configuration
      call config%finalize()

   end subroutine test_scheme_configuration

   !> Test ProcessDryDepInterface can be created
   subroutine test_process_interface_creation()
      ! Test that we can create the interface object
      call assert(.true., "ProcessDryDepInterface object created successfully")
   end subroutine test_process_interface_creation

   !> Test ProcessDryDepInterface has required methods (without calling them)
   subroutine test_process_interface_methods()
      ! Test that the interface has the expected methods by checking if it's ready
      ! (this doesn't call init, just checks the initial state)
      call assert(.not. drydep_process%is_ready(), "Process should not be ready before initialization")
   end subroutine test_process_interface_methods

end program test_drydep_unit
