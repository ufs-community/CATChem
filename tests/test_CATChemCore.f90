!> \file test_CATChemCore.f90
!! \brief Test program for CATChemCore module
!!
!!!>
program test_CATChemCore
   use testing_mod, only: assert, assert_close
   use CATChemCore_Mod
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE, ErrorManagerType, ERROR_PROCESS_INITIALIZATION
   use StateManager_Mod, only: StateManagerType
   use GridManager_Mod, only: GridManagerType
   use ProcessManager_Mod, only: ProcessManagerType
   use DiagnosticManager_Mod, only: DiagnosticManagerType
   use ConfigManager_Mod, only: ConfigDataType


   implicit none

   type(CATChemCoreType) :: core
   type(CATChemBuilderType) :: builder
   integer :: rc
   logical :: is_ready

   write(*,*) 'Testing CATChemCore module...'
   write(*,*) ''

   ! Test 1: Basic initialization
   write(*,*) 'Test 1: Basic initialization'
   call core%init('TestCore', rc)
   call assert(rc == CC_SUCCESS, "Core initialization should succeed")

   ! Check if core is initialized
   is_ready = core%is_ready()
   call assert(.not. is_ready, "Core should not be ready before configuration")

   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! Test 2: Configuration
   write(*,*) 'Test 2: Configuration'
   call core%configure(nx=10, ny=10, nz=20, rc=rc, nsoil=4, nsoiltype=19, nsurftype=13)
   call assert(rc == CC_SUCCESS, "Core configuration should succeed")

   ! Check if core is ready after configuration
   is_ready = core%is_ready()
   call assert(is_ready, "Core should be ready after configuration")

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: Component access
   write(*,*) 'Test 3: Component access'
   block
      type(StateManagerType), pointer :: state_mgr
      type(GridManagerType), pointer :: grid_mgr
      type(ProcessManagerType), pointer :: process_mgr
      type(DiagnosticManagerType), pointer :: diag_mgr
      type(ConfigDataType), pointer :: config_data
      type(ErrorManagerType), pointer :: error_mgr

      state_mgr => core%get_state_manager()
      call assert(associated(state_mgr), "Should be able to get state manager")

      grid_mgr => core%get_grid_manager()
      call assert(associated(grid_mgr), "Should be able to get grid manager")

      process_mgr => core%get_process_manager()
      call assert(associated(process_mgr), "Should be able to get process manager")

      diag_mgr => core%get_diagnostic_manager()
      call assert(associated(diag_mgr), "Should be able to get diagnostic manager")

      config_data => core%get_config()
      call assert(associated(config_data), "Should be able to get config data")

      error_mgr => core%get_error_manager()
      call assert(associated(error_mgr), "Should be able to get error manager")
   end block

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! Test 4: Builder pattern
   write(*,*) 'Test 4: Builder pattern'
   call builder%init()
   builder = builder%with_name('BuilderTest')
   builder = builder%with_grid(5, 5, 10)

   ! Build core using builder
   call builder%build(core, rc)
   call assert(rc == CC_SUCCESS, "Builder should be able to build core")

   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! Test 5: Memory usage
   write(*,*) 'Test 5: Memory usage'
   block
      integer(kind=8) :: memory_usage
      memory_usage = core%get_memory_usage()
      call assert(memory_usage >= 0, "Memory usage should be non-negative")
   end block

   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   ! Test 6: Print info
   write(*,*) 'Test 6: Print info'
   call core%print_info()
   ! Should complete without error

   write(*,*) 'Test 6 passed!'
   write(*,*) ''

   ! Test 7: Finalization
   write(*,*) 'Test 7: Finalization'
   call core%finalize(rc)
   call assert(rc == CC_SUCCESS, "Core finalization should succeed")

   ! Check if core is no longer ready
   is_ready = core%is_ready()
   call assert(.not. is_ready, "Core should not be ready after finalization")

   write(*,*) 'Test 7 passed!'
   write(*,*) ''

   write(*,*) 'All CATChemCore tests passed!'

end program test_CATChemCore
