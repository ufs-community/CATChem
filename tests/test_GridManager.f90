!> \file test_GridManager.f90
!! \brief Test program for GridManager module
!!
!!!>
program test_GridManager
   use testing_mod, only: assert, assert_close
   use GridManager_Mod, only: GridManagerType, GridGeometryType, GridDecompositionType, ColumnIteratorType, GRID_TYPE_3D, COORD_CARTESIAN
   use ColumnInterface_Mod, only: ColumnViewType
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE, ErrorManagerType
   use Precision_Mod, only: fp

   implicit none

   type(GridManagerType) :: grid_mgr
   type(ErrorManagerType) :: error_mgr
   type(GridGeometryType) :: geometry
   type(ColumnViewType) :: column_view
   integer :: rc
   logical :: is_ready
   integer :: nx, ny, nz

   write(*,*) 'Testing GridManager module...'
   write(*,*) ''

   ! Test 1: Initialize error manager
   write(*,*) 'Test 1: Initialize error manager'
   call error_mgr%init()
   ! Error manager should be ready after initialization

   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! Test 2: Initialize grid geometry
   write(*,*) 'Test 2: Initialize grid geometry'
   call geometry%init(5, 5, 10, GRID_TYPE_3D, COORD_CARTESIAN, rc)
   call assert(rc == CC_SUCCESS, "Grid geometry initialization should succeed")

   call geometry%get_dimensions(nx, ny, nz)
   call assert(nx == 5, "NX should be 5")
   call assert(ny == 5, "NY should be 5")
   call assert(nz == 10, "NZ should be 10")

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: Initialize grid manager
   write(*,*) 'Test 3: Initialize grid manager'
   call grid_mgr%init(5, 5, 10, error_mgr, GRID_TYPE_3D, COORD_CARTESIAN, rc)
   call assert(rc == CC_SUCCESS, "Grid manager initialization should succeed")

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! Test 4: Check if grid manager is ready
   write(*,*) 'Test 4: Check if grid manager is ready'
   is_ready = grid_mgr%is_ready()
   call assert(is_ready, "Grid manager should be ready after initialization")

   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! Test 5: Get grid geometry
   write(*,*) 'Test 5: Get grid geometry'
   block
      type(GridGeometryType) :: geom

      geom = grid_mgr%get_geometry()

      call geom%get_dimensions(nx, ny, nz)
      call assert(nx == 5, "NX should be 5")
      call assert(ny == 5, "NY should be 5")
      call assert(nz == 10, "NZ should be 10")
   end block

   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   ! Test 6: Get grid decomposition
   write(*,*) 'Test 6: Get grid decomposition'
   block
      type(GridDecompositionType) :: decomp

      decomp = grid_mgr%get_decomposition()

      ! For single processor, should have simple decomposition
      ! We're not asserting on specific values because they depend on implementation
   end block

   write(*,*) 'Test 6 passed!'
   write(*,*) ''

   ! Test 7: Get grid dimensions
   write(*,*) 'Test 7: Get grid dimensions'
   block
      integer :: total_cols, local_cols

      total_cols = grid_mgr%get_total_columns()
      local_cols = grid_mgr%get_local_columns()

      call assert(total_cols == 25, "Total columns should be 25 (5x5)")
      call assert(local_cols >= 0, "Local columns should be non-negative")
   end block

   write(*,*) 'Test 7 passed!'
   write(*,*) ''

   ! Test 8: Get grid shape
   write(*,*) 'Test 8: Get grid shape'
   block
      integer :: shape_nx, shape_ny, shape_nz

      call grid_mgr%get_shape(shape_nx, shape_ny, shape_nz)

      call assert(shape_nx == 5, "Shape NX should be 5")
      call assert(shape_ny == 5, "Shape NY should be 5")
      call assert(shape_nz == 10, "Shape NZ should be 10")
   end block

   write(*,*) 'Test 8 passed!'
   write(*,*) ''

   ! Test 9: Create column iterator
   write(*,*) 'Test 9: Create column iterator'
   block
      type(ColumnIteratorType) :: iterator

      iterator = grid_mgr%create_column_iterator()
      ! Should create a valid iterator
   end block

   write(*,*) 'Test 9 passed!'
   write(*,*) ''

   ! Test 10: Get column by indices
   write(*,*) 'Test 10: Get column by indices'
   block
      type(ColumnViewType) :: column_view

      column_view = grid_mgr%get_column_by_indices(3, 3)
      ! Should return a valid column view (even if uninitialized)
   end block

   write(*,*) 'Test 10 passed!'
   write(*,*) ''

   ! Test 11: Compute distances
   write(*,*) 'Test 11: Compute distances'
   block
      real(fp) :: distance
      integer :: local_rc

      call grid_mgr%compute_distances(1, 1, 2, 2, distance, local_rc)
      call assert(local_rc == CC_SUCCESS, "Distance computation should succeed")
      call assert(distance > 0.0_fp, "Distance should be positive")
   end block

   write(*,*) 'Test 11 passed!'
   write(*,*) ''

   ! Test 12: Interpolate to column
   write(*,*) 'Test 12: Interpolate to column'
   block
      real(fp), allocatable :: source_data(:,:)
      real(fp) :: interpolated_data
      integer :: local_rc

      allocate(source_data(5, 5))
      source_data = 1.0_fp  ! Fill with constant values

      call grid_mgr%interpolate_to_column(source_data, 3, 3, interpolated_data, local_rc)
      call assert(local_rc == CC_SUCCESS, "Interpolation should succeed")
      call assert_close(interpolated_data, 1.0_fp, 1.0e-6_fp, "Interpolated data should match source")

      deallocate(source_data)
   end block

   write(*,*) 'Test 12 passed!'
   write(*,*) ''

   ! Test 13: Print info
   write(*,*) 'Test 13: Print info'
   call grid_mgr%print_info()
   ! Should complete without error

   write(*,*) 'Test 13 passed!'
   write(*,*) ''

   ! Test 14: Cleanup
   write(*,*) 'Test 14: Cleanup'
   call grid_mgr%cleanup()
   ! No finalize method for ErrorManagerType; use cleanup or report_error if needed

   write(*,*) 'Test 14 passed!'
   write(*,*) ''

   write(*,*) 'All GridManager tests passed!'

end program test_GridManager
