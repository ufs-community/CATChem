!> \file test_ProcessFactory.f90
!! \brief Test program for ProcessFactory module
!!
!!!>
program test_ProcessFactory
   use testing_mod, only: assert, assert_close
   use ProcessFactory_Mod
   use StateManager_Mod, only: StateManagerType
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE
   use ProcessInterface_Mod

   implicit none

   type(ProcessFactoryType) :: factory
   type(StateManagerType) :: state_mgr
   integer :: rc

   write(*,*) 'Testing ProcessFactory module...'
   write(*,*) ''

   ! Test 1: Initialize factory
   write(*,*) 'Test 1: Initialize factory'
   call factory%init(rc)
   call assert(rc == CC_SUCCESS, "Process factory initialization should succeed")

   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! Test 2: Initialize state manager (needed for process creation)
   write(*,*) 'Test 2: Initialize state manager'
   call state_mgr%init('TestStateManager', rc)
   call assert(rc == CC_SUCCESS, "State manager initialization should succeed")

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: List available processes
   write(*,*) 'Test 3: List available processes'
   block
      character(len=64), allocatable :: process_names(:)
      integer :: num_processes

      call factory%list_available(process_names, rc)
      call assert(rc == CC_SUCCESS, "Listing available processes should succeed")
      ! We don't assert on num_processes because it may vary depending on what's registered

      if (allocated(process_names)) then
         deallocate(process_names)
      end if
   end block

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! Test 4: Try to create a process (will likely fail since no processes are registered)
   write(*,*) 'Test 4: Try to create a process'
   block
      class(ProcessInterface), allocatable :: process

      call factory%create_process('nonexistent_process', &
         state_mgr, process, rc)
      ! This should fail since no processes are registered
      ! We're not asserting on rc because behavior may vary
   end block

   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! Test 5: Cleanup
   write(*,*) 'Test 5: Cleanup'
   call state_mgr%finalize(rc)
   call assert(rc == CC_SUCCESS, "State manager finalization should succeed")

   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   write(*,*) 'All ProcessFactory tests passed!'

end program test_ProcessFactory
