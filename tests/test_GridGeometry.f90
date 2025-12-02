!> \file test_GridGeometry.f90
!! \brief Test program for GridGeometry module
!!
!!!>
program test_GridGeometry
   use testing_mod, only: assert, assert_close
   use GridGeometry_Mod
   use GridGeometry_Mod, only: GridGeometryType

   implicit none

   type(GridGeometryType) :: geometry
   integer :: nx, ny, nz

   write(*,*) 'Testing GridGeometry module...'
   write(*,*) ''

   ! Test 1: Default initialization
   write(*,*) 'Test 1: Default initialization'
   nx = geometry%get_nx()
   ny = geometry%get_ny()
   nz = geometry%get_nz()

   call assert(nx == 1, "Default NX should be 1")
   call assert(ny == 1, "Default NY should be 1")
   call assert(nz == 1, "Default NZ should be 1")

   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! Test 2: Setting dimensions
   write(*,*) 'Test 2: Setting dimensions'
   call geometry%set(10, 20, 30)

   nx = geometry%get_nx()
   ny = geometry%get_ny()
   nz = geometry%get_nz()

   call assert(nx == 10, "NX should be 10 after setting")
   call assert(ny == 20, "NY should be 20 after setting")
   call assert(nz == 30, "NZ should be 30 after setting")

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: Getting dimensions as array
   write(*,*) 'Test 3: Getting dimensions as array'
   call geometry%get_dimensions(nx, ny, nz)

   call assert(nx == 10, "NX should be 10")
   call assert(ny == 20, "NY should be 20")
   call assert(nz == 30, "NZ should be 30")

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   write(*,*) 'All GridGeometry tests passed!'

end program test_GridGeometry
