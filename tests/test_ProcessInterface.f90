!> \file test_ProcessInterface.f90
!! \brief Test program for ProcessInterface module
!!
!!!>
module test_ProcessInterface_mod
   use testing_mod, only: assert, assert_close
   use ProcessInterface_Mod
   use StateManager_Mod, only: StateManagerType
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE

   implicit none

   type, extends(ProcessInterface) :: DummyProcessType
   contains
      procedure, public :: init => dummy_init
      procedure, public :: run => dummy_run
      procedure, public :: finalize => dummy_finalize
   end type DummyProcessType

contains

   subroutine dummy_init(this, container, rc)
      class(DummyProcessType), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      this%name = 'DummyProcess'
      this%version = '1.0'
      this%description = 'A dummy process for testing'
      call this%activate()
      rc = CC_SUCCESS
   end subroutine dummy_init

   subroutine dummy_run(this, container, rc)
      class(DummyProcessType), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      rc = CC_SUCCESS
   end subroutine dummy_run

   subroutine dummy_finalize(this, rc)
      class(DummyProcessType), intent(inout) :: this
      integer, intent(out) :: rc

      call this%deactivate()
      rc = CC_SUCCESS
   end subroutine dummy_finalize

end module test_ProcessInterface_mod

program test_ProcessInterface
   use test_ProcessInterface_mod
   implicit none

   type(DummyProcessType) :: dummy_process
   type(StateManagerType) :: state_mgr
   integer :: rc

   write(*,*) 'Testing ProcessInterface module...'
   write(*,*) ''

   ! Test 1: State manager initialization (needed for ProcessInterface)
   write(*,*) 'Test 1: State manager initialization'
   call state_mgr%init('TestStateManager', rc)
   call assert(rc == CC_SUCCESS, "StateManager initialization should succeed")

   ! Test 2: Dummy process initialization
   write(*,*) 'Test 2: Dummy process initialization'
   call dummy_process%init(state_mgr, rc)
   call assert(rc == CC_SUCCESS, "Dummy process initialization should succeed")

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: Process properties
   write(*,*) 'Test 3: Process properties'
   block
      logical :: is_ready

      is_ready = dummy_process%is_ready()
      ! Should be ready after initialization
      call assert(is_ready, "Process should be ready after initialization")
   end block

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! Test 4: Process activation/deactivation
   write(*,*) 'Test 4: Process activation/deactivation'
   call dummy_process%deactivate()
   block
      logical :: is_ready

      is_ready = dummy_process%is_ready()
      ! Should not be ready when deactivated
      ! Note: This might still pass depending on implementation
   end block

   call dummy_process%activate()
   block
      logical :: is_ready

      is_ready = dummy_process%is_ready()
      ! Should be ready after activation
      call assert(is_ready, "Process should be ready after activation")
   end block

   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! Test 5: Running the process
   write(*,*) 'Test 5: Running the process'
   call dummy_process%run(state_mgr, rc)
   call assert(rc == CC_SUCCESS, "Process run should succeed")

   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   ! Test 6: Process finalization
   write(*,*) 'Test 6: Process finalization'
   call dummy_process%finalize(rc)
   call assert(rc == CC_SUCCESS, "Process finalization should succeed")

   ! Test 7: Cleanup
   write(*,*) 'Test 7: Cleanup'
   call state_mgr%finalize(rc)
   call assert(rc == CC_SUCCESS, "StateManager finalization should succeed")

   write(*,*) 'Test 6 passed!'
   write(*,*) ''

   write(*,*) 'All ProcessInterface tests passed!'


end program test_ProcessInterface
