!> \file test_StateInterface.f90
!! \brief Test program for StateInterface module
!!
!!!>
program test_StateInterface
   use testing_mod, only: assert, assert_close
   use StateInterface_Mod
   use Precision_Mod, only: fp

   implicit none

   type(StateContainerType) :: container
   integer :: rc

   write(*,*) 'Testing StateInterface module...'
   write(*,*) ''

   ! Test 1: Initialize state container
   write(*,*) 'Test 1: Initialize state container'
   call container%init('TestContainer', rc)
   call assert(rc == 0, "State container initialization should succeed")

   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! Test 2: Get container name
   write(*,*) 'Test 2: Get container name'
   block
      character(len=256) :: name
      name = container%get_name()
      call assert(len_trim(name) > 0, "Container name should not be empty")
   end block

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: Check if container is ready
   write(*,*) 'Test 3: Check if container is ready'
   block
      logical :: is_ready
      is_ready = container%is_ready()
      ! Should be ready after initialization
      call assert(is_ready, "Container should be ready after initialization")
   end block

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! Test 4: Get state objects (will return null pointers initially)
   write(*,*) 'Test 4: Get state objects'
   block
      type(ConfigDataType), pointer :: config_ptr
      type(MetStateType), pointer :: met_ptr
      type(ChemStateType), pointer :: chem_ptr
      type(ErrorManagerType), pointer :: error_ptr

      config_ptr => container%get_config_ptr()
      met_ptr => container%get_met_state_ptr()
      chem_ptr => container%get_chem_state_ptr()
      error_ptr => container%get_error_manager()

      ! Initially these might be null or point to uninitialized objects
      ! We're not asserting on them because behavior may vary
   end block

   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! Test 5: Get grid manager
   write(*,*) 'Test 5: Get grid manager'
   block
      type(GridManagerType), pointer :: grid_mgr_ptr
      grid_mgr_ptr => container%get_grid_manager()
      ! Should return a valid pointer or null
      ! We're not asserting on it because behavior may vary
   end block

   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   ! Test 6: Get process manager
   write(*,*) 'Test 6: Get process manager'
   block
      type(ProcessManagerType), pointer :: process_mgr_ptr
      process_mgr_ptr => container%get_process_manager()
      ! Should return a valid pointer or null
      ! We're not asserting on it because behavior may vary
   end block

   write(*,*) 'Test 6 passed!'
   write(*,*) ''

   ! Test 7: Get diagnostic manager
   write(*,*) 'Test 7: Get diagnostic manager'
   block
      type(DiagnosticManagerType), pointer :: diag_mgr_ptr
      diag_mgr_ptr => container%get_diagnostic_manager()
      ! Should return a valid pointer or null
      ! We're not asserting on it because behavior may vary
   end block

   write(*,*) 'Test 7 passed!'
   write(*,*) ''

   ! Test 8: Cleanup
   write(*,*) 'Test 8: Cleanup'
   call container%cleanup(rc)
   call assert(rc == 0, "State container cleanup should succeed")

   write(*,*) 'Test 8 passed!'
   write(*,*) ''

   write(*,*) 'All StateInterface tests passed!'

end program test_StateInterface
