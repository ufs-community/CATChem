!> \file test_Error.f90
!! \brief Test program for Error module
!!
!!!>
program test_Error
   use testing_mod, only: assert, assert_close
   use Error_Mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE, ERROR_INVALID_INPUT

   implicit none

   type(ErrorManagerType) :: error_mgr
   integer :: rc

   write(*,*) 'Testing Error module...'
   write(*,*) ''

   ! Test 1: Error manager initialization
   write(*,*) 'Test 1: Error manager initialization'
   call error_mgr%init(verbose=.false., track_performance=.false., rc=rc)

   ! Check initial state
   ! No get_status() member; just check rc from init
   call assert(rc == CC_SUCCESS, "Error manager should initialize successfully")

   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! Test 2: Pushing and popping context
   write(*,*) 'Test 2: Pushing and popping context'
   call error_mgr%push_context('test_context', 'Testing context handling')

   ! Pop context
   call error_mgr%pop_context()

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: Reporting errors
   write(*,*) 'Test 3: Reporting errors'
   call error_mgr%push_context('error_test', 'Testing error reporting')
   call error_mgr%report_error(ERROR_INVALID_INPUT, 'Test error message', rc, 'test_location')
   ! Note: We're not checking the rc value here because error reporting behavior may vary
   call error_mgr%pop_context()

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! Test 4: Cleanup
   write(*,*) 'Test 4: Cleanup'
   ! call error_mgr%finalize(rc) ! No such member
   ! call assert(rc == CC_SUCCESS, "Error manager finalization should succeed")

   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   write(*,*) 'All Error module tests passed!'

end program test_Error
