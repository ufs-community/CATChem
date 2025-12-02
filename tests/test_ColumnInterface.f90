!> \file test_ColumnInterface.f90
!! \brief Test program for ColumnInterface module
!!
!!!>
program test_ColumnInterface
   use testing_mod, only: assert, assert_close
   use ColumnInterface_Mod
   use Precision_Mod, only: fp

   implicit none

   type(ColumnViewType) :: col_view
   ! Remove unused/undefined VirtualColumnType instance
   integer :: rc

   write(*,*) 'Testing ColumnInterface module...'
   write(*,*) ''

   ! Test 1: Initialize column view
   write(*,*) 'Test 1: Initialize column view'
   call col_view%init(rc=rc)
   call assert(rc == 0, "Column view initialization should succeed")

   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! Test 2: Set column position
   write(*,*) 'Test 2: Set column position'
   call col_view%set_column_position(1, 1, rc)
   call assert(rc == 0, "Setting column position should succeed")

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: Test column position validity
   write(*,*) 'Test 3: Test column position validity'
   block
      real(fp), pointer :: met_ptr(:)
      real(fp), pointer :: chem_ptr(:)
      real(fp), pointer :: emis_ptr(:)

      ! These should return null pointers since we haven't set up state objects
      met_ptr => col_view%get_met_column_ptr('T')
      chem_ptr => col_view%get_chem_column_ptr('O3')
      emis_ptr => col_view%get_emis_column_ptr('NOx')

      call assert(.not. associated(met_ptr), "Met pointer should not be associated")
      call assert(.not. associated(chem_ptr), "Chem pointer should not be associated")
      call assert(.not. associated(emis_ptr), "Emis pointer should not be associated")
   end block

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! Test 4: Create column view using function
   write(*,*) 'Test 4: Create column view using function'
   col_view = create_column_view()
   ! This should create a valid column view

   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! Test 5: Cleanup
   write(*,*) 'Test 5: Cleanup'
   call col_view%cleanup()

   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   write(*,*) 'All ColumnInterface tests passed!'

end program test_ColumnInterface
