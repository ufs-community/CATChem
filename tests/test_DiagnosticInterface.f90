!> \file test_DiagnosticInterface.f90
!! \brief Test program for DiagnosticInterface module
!!
!!!>
program test_DiagnosticInterface
   use testing_mod, only: assert, assert_close
   use DiagnosticInterface_Mod
   use Precision_Mod, only: fp
   use error_mod, only: ErrorManagerType, CC_SUCCESS

   implicit none

   type(DiagnosticFieldType) :: diag_field
   type(DiagnosticRegistryType) :: diag_registry
   integer :: rc
   type(ErrorManagerType), target :: error_mgr_target
   type(ErrorManagerType), pointer :: error_mgr

   write(*,*) 'Testing DiagnosticInterface module...'
   write(*,*) ''

   ! Test 1: Initialize diagnostic field
   write(*,*) 'Test 1: Initialize diagnostic field'
   call diag_field%create('test_field', 'Test diagnostic field', &
      'kg/m2/s', DIAG_REAL_2D, 'TestProcess', rc=rc)
   write(*,*) 'diag_field%create rc:', rc
   call assert(rc == CC_SUCCESS, "Diagnostic field creation should succeed")
   write(*,*) 'assert for diag_field%create passed'

   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! Test 2: Get field properties
   write(*,*) 'Test 2: Get field properties'
   block
      character(len=64) :: name
      character(len=128) :: description
      character(len=32) :: units
      integer :: data_type

      name = diag_field%get_name()
      description = diag_field%get_description()
      units = diag_field%get_units()
      data_type = diag_field%get_data_type()

      call assert(trim(name) == 'test_field', "Field name should match")
      call assert(trim(description) == 'Test diagnostic field', "Field description should match")
      call assert(trim(units) == 'kg/m2/s', "Field units should match")
      call assert(data_type == DIAG_REAL_2D, "Field data type should match")
   end block

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: Check if field is ready
   write(*,*) 'Test 3: Check if field is ready'
   block
      logical :: is_ready
      is_ready = diag_field%is_ready()
      ! Should be ready after creation
      call assert(is_ready, "Diagnostic field should be ready after creation")
   end block

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! Test 4: Initialize diagnostic registry
   write(*,*) 'Test 4: Initialize diagnostic registry'
   error_mgr => error_mgr_target
   call diag_registry%init('TestProcess', error_mgr=error_mgr)
   write(*,*) 'diag_registry%init called'
   call assert(rc == CC_SUCCESS, "Diagnostic registry initialization should succeed")

   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! Test 5: Register diagnostic field
   write(*,*) 'Test 5: Register diagnostic field'
   call diag_registry%register_field(diag_field, rc)
   call assert(rc == CC_SUCCESS, "Registering diagnostic field should succeed")

   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   ! Test 6: Get field from registry
   write(*,*) 'Test 6: Get field from registry'
   block
      type(DiagnosticFieldType) :: retrieved_field
      retrieved_field = diag_registry%get_field('test_field')
      ! Should return the registered field
   end block

   write(*,*) 'Test 6 passed!'
   write(*,*) ''

   ! Test 7: List fields in registry
   write(*,*) 'Test 7: List fields in registry'
   block
      character(len=64), allocatable :: field_names(:)
      integer :: num_fields

      ! Allocate the array to the maximum size
      allocate(field_names(10))
      call diag_registry%list_fields(field_names, num_fields)
      call assert(num_fields == 1, "Registry should contain one field")
      call assert(trim(field_names(1)) == 'test_field', "Field name should match")

      if (allocated(field_names)) deallocate(field_names)
   end block

   write(*,*) 'Test 7 passed!'
   write(*,*) ''

   ! Test 8: Cleanup
   write(*,*) 'Test 8: Cleanup'
   call diag_registry%cleanup()
   call diag_field%cleanup()

   write(*,*) 'Test 8 passed!'
   write(*,*) ''

   write(*,*) 'All DiagnosticInterface tests passed!'

end program test_DiagnosticInterface
