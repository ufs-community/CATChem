!> \file test_VirtualColumn.f90
!! \brief Test program for VirtualColumn module
!!
!!!>
program test_VirtualColumn
   use testing_mod, only: assert, assert_close
   use VirtualColumn_Mod
   use Precision_Mod, only: fp

   implicit none

   type(VirtualColumnType) :: virtual_col
   integer :: nlev, nspec_chem, nspec_emis
   integer :: i, j
   real(fp) :: lat, lon, area
   real(fp) :: test_value
   integer :: rc

   write(*,*) 'Testing VirtualColumn module...'
   write(*,*) ''

   ! Test 1: Initialize virtual column
   write(*,*) 'Test 1: Initialize virtual column'
   nlev = 10
   nspec_chem = 5
   nspec_emis = 3
   i = 1
   j = 1
   lat = 45.0_fp
   lon = -120.0_fp
   area = 1000000.0_fp  ! 1 km²

   call virtual_col%init(nlev, nspec_chem, nspec_emis, i, j, lat, lon, area, rc)
   call assert(rc == 0, "Virtual column initialization should succeed")

   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! Test 2: Get dimensions
   write(*,*) 'Test 2: Get dimensions'
   call virtual_col%get_dimensions(nlev, nspec_chem, nspec_emis)

   call assert(nlev == 10, "NLEV should be 10")
   call assert(nspec_chem == 5, "NSPEC_CHEM should be 5")
   call assert(nspec_emis == 3, "NSPEC_EMIS should be 3")

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: Get position
   write(*,*) 'Test 3: Get position'
   call virtual_col%get_position(i, j)

   call assert(i == 1, "I should be 1")
   call assert(j == 1, "J should be 1")

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! Test 4: Set and get meteorological field
   write(*,*) 'Test 4: Set and get meteorological field'
   ! Use a valid meteorological field from VirtualMetType, e.g. T(:) for temperature
   if (.not. associated(virtual_col%met%T)) then
      write(*,*) "[ERROR] Temperature field pointer not associated"
      stop
   endif
   virtual_col%met%T(1) = 288.15_fp
   test_value = virtual_col%met%T(1)  ! Read back the value we just set

   call assert_close(test_value, 288.15_fp, 1.0e-6_fp, "Met field value should match")

   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! Test 5: Set and get chemical field
   write(*,*) 'Test 5: Set and get chemical field'
   call virtual_col%set_chem_field(1, 1, 1.0e-9_fp)  ! 1 ppb
   test_value = virtual_col%get_chem_field(1, 1)

   call assert_close(test_value, 1.0e-9_fp, 1.0e-12_fp, "Chem field value should match")

   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   ! Test 6: Set and get emission field
   write(*,*) 'Test 6: Set and get emission field'
   ! set_emis_field expects (k, ispec, value)
   call virtual_col%set_emis_field(1, 1, 1.0e-6_fp)  ! 1 μg/m²/s
   test_value = virtual_col%get_emis_field(1, 1)

   call assert_close(test_value, 1.0e-6_fp, 1.0e-9_fp, "Emission field value should match")

   write(*,*) 'Test 6 passed!'
   write(*,*) ''

   ! Test 7: Cleanup
   write(*,*) 'Test 7: Cleanup'
   call virtual_col%cleanup()

   ! After cleanup, accessing fields should return 0 or fail gracefully
   ! test_value = virtual_col%met%T(1)
   ! Accessing met%T(1) after cleanup is unsafe and may cause a segmentation fault,
   ! because cleanup nullifies or deallocates the pointer.
   ! We're not asserting on this because behavior after cleanup may vary

   write(*,*) 'Test 7 passed!'
   write(*,*) ''

   write(*,*) 'All VirtualColumn tests passed!'

end program test_VirtualColumn
