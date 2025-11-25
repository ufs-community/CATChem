

# File GridManager\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**GridManager\_Mod.F90**](_grid_manager___mod_8_f90.md)

[Go to the documentation of this file](_grid_manager___mod_8_f90.md)


```Fortran

module gridmanager_mod
   use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure, errormanagertype, error_process_initialization
   use columninterface_mod, only: columnviewtype
   use virtualcolumn_mod, only: virtualcolumntype

   implicit none
   private

   public :: gridmanagertype
   public :: gridgeometrytype
   public :: columniteratortype
   public :: griddecompositiontype

   ! Grid type constants
   integer, parameter, public :: GRID_TYPE_COLUMN = 1
   integer, parameter, public :: GRID_TYPE_2D     = 2
   integer, parameter, public :: GRID_TYPE_3D     = 3

   ! Coordinate system constants
   integer, parameter, public :: COORD_CARTESIAN    = 1
   integer, parameter, public :: COORD_LONLAT       = 2
   integer, parameter, public :: COORD_PROJECTED    = 3

   type :: gridgeometrytype
      private

      ! Grid dimensions
      integer :: nx = 1
      integer :: ny = 1
      integer :: nz = 1
      integer :: grid_type = grid_type_column 
      integer :: coord_system = coord_cartesian 

      ! Grid spacing and domain
      real(fp) :: dx = 1000.0_fp  
      real(fp) :: dy = 1000.0_fp  
      real(fp), allocatable :: dz(:)
      real(fp), allocatable :: z_levels(:)

      ! Domain bounds
      real(fp) :: x_min = 0.0_fp  
      real(fp) :: x_max = 1000.0_fp 
      real(fp) :: y_min = 0.0_fp  
      real(fp) :: y_max = 1000.0_fp 

      ! Geographic information (for lon/lat grids)
      real(fp), allocatable :: lon(:,:)
      real(fp), allocatable :: lat(:,:)
      real(fp), allocatable :: grid_area(:,:)

      logical :: is_initialized = .false. 

   contains
      procedure :: init => geometry_init
      procedure :: cleanup => geometry_cleanup
      procedure :: get_nx => geometry_get_nx
      procedure :: get_ny => geometry_get_ny
      procedure :: get_nz => geometry_get_nz
      procedure :: get_grid_type => geometry_get_grid_type
      procedure :: get_dimensions => geometry_get_dimensions
      procedure :: get_cell_area => geometry_get_cell_area
      procedure :: get_cell_volume => geometry_get_cell_volume
      procedure :: is_valid_position => geometry_is_valid_position
   end type gridgeometrytype

   type :: griddecompositiontype
      private

      ! Decomposition parameters
      integer :: n_procs = 1
      integer :: my_rank = 0
      integer :: i_start = 1
      integer :: i_end = 1
      integer :: j_start = 1
      integer :: j_end = 1

      ! Halo/ghost cell information
      integer :: halo_width = 0
      logical :: has_halos = .false. 

   contains
      procedure :: init => decomp_init
      procedure :: get_local_bounds => decomp_get_local_bounds
      procedure :: is_local_column => decomp_is_local_column
   end type griddecompositiontype

   type :: columniteratortype
      private

      type(GridManagerType), pointer :: grid_mgr => null()
      integer :: current_i = 1
      integer :: current_j = 1
      integer :: total_columns = 1
      integer :: current_column = 0

   contains
      procedure :: init => iterator_init
      procedure :: has_next => iterator_has_next
      procedure :: next => iterator_next
      procedure :: get_current_column => iterator_get_current_column
      procedure :: get_current_indices => iterator_get_current_indices
      procedure :: reset => iterator_reset
   end type columniteratortype

   type :: gridmanagertype
      private

      type(GridGeometryType) :: geometry
      type(GridDecompositionType) :: decomp
      type(ErrorManagerType), pointer :: error_mgr => null() 

      ! Grid state
      logical :: is_initialized = .false.       
      character(len=256) :: name = 'GridManager'

      ! Column management
      integer :: n_local_columns = 0

   contains
      ! Initialization and cleanup
      procedure :: init => grid_manager_init
      procedure :: cleanup => grid_manager_cleanup
      procedure :: validate => grid_manager_validate

      ! Grid access methods
      procedure :: get_geometry => grid_manager_get_geometry
      procedure :: get_decomposition => grid_manager_get_decomposition
      procedure :: get_total_columns => grid_manager_get_total_columns
      procedure :: get_local_columns => grid_manager_get_local_columns
      procedure :: get_shape => grid_manager_get_shape

      ! Column virtualization interface
      procedure :: create_column_iterator => grid_manager_create_column_iterator
      procedure :: get_column_by_indices => grid_manager_get_column_by_indices
      procedure :: get_column_by_location => grid_manager_get_column_by_location
      procedure :: create_virtual_column => grid_manager_create_virtual_column

      ! Grid operations
      procedure :: compute_distances => grid_manager_compute_distances
      procedure :: interpolate_to_column => grid_manager_interpolate_to_column
      procedure :: exchange_halo_data => grid_manager_exchange_halo_data

      ! Utilities
      procedure :: print_info => grid_manager_print_info
      procedure :: is_ready => grid_manager_is_ready
   end type gridmanagertype

contains

   !========================================================================
   ! GridGeometryType Implementation
   !========================================================================

   subroutine geometry_init(this, nx, ny, nz, grid_type, coord_system, rc)
      class(GridGeometryType), intent(inout) :: this
      integer, intent(in) :: nx, ny, nz
      integer, intent(in), optional :: grid_type, coord_system
      integer, intent(out) :: rc

      integer :: alloc_stat
      integer :: i

      rc = cc_success

      ! Set dimensions
      this%nx = max(1, nx)
      this%ny = max(1, ny)
      this%nz = max(1, nz)

      ! Set grid type
      if (present(grid_type)) then
         this%grid_type = grid_type
      else
         if (this%nx == 1 .and. this%ny == 1) then
            this%grid_type = grid_type_column
         else if (this%ny == 1) then
            this%grid_type = grid_type_2d
         else
            this%grid_type = grid_type_3d
         endif
      endif

      ! Set coordinate system
      if (present(coord_system)) then
         this%coord_system = coord_system
      else
         this%coord_system = coord_cartesian
      endif

      ! Allocate vertical arrays
      allocate(this%dz(this%nz), this%z_levels(this%nz+1), stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = cc_failure
         return
      endif

      ! Initialize with default values
      this%dz = 100.0_fp  ! Default 100m layers
      this%z_levels(1) = 0.0_fp
      do i = 1, this%nz
         this%z_levels(i+1) = this%z_levels(i) + this%dz(i)
      enddo

      ! Allocate geographic arrays if needed
      if (this%coord_system == coord_lonlat) then
         allocate(this%lon(this%nx, this%ny), this%lat(this%nx, this%ny), &
            this%grid_area(this%nx, this%ny), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = cc_failure
            return
         endif
      endif

      this%is_initialized = .true.

   end subroutine geometry_init

   subroutine geometry_cleanup(this)
      class(GridGeometryType), intent(inout) :: this

      if (allocated(this%dz)) deallocate(this%dz)
      if (allocated(this%z_levels)) deallocate(this%z_levels)
      if (allocated(this%lon)) deallocate(this%lon)
      if (allocated(this%lat)) deallocate(this%lat)
      if (allocated(this%grid_area)) deallocate(this%grid_area)

      this%is_initialized = .false.

   end subroutine geometry_cleanup

   function geometry_get_nx(this) result(nx)
      class(GridGeometryType), intent(in) :: this
      integer :: nx
      nx = this%nx
   end function geometry_get_nx

   function geometry_get_ny(this) result(ny)
      class(GridGeometryType), intent(in) :: this
      integer :: ny
      ny = this%ny
   end function geometry_get_ny

   function geometry_get_nz(this) result(nz)
      class(GridGeometryType), intent(in) :: this
      integer :: nz
      nz = this%nz
   end function geometry_get_nz

   function geometry_get_grid_type(this) result(grid_type)
      class(GridGeometryType), intent(in) :: this
      integer :: grid_type
      grid_type = this%grid_type
   end function geometry_get_grid_type

   subroutine geometry_get_dimensions(this, nx, ny, nz)
      class(GridGeometryType), intent(in) :: this
      integer, intent(out) :: nx, ny, nz

      nx = this%nx
      ny = this%ny
      nz = this%nz
   end subroutine geometry_get_dimensions

   function geometry_get_cell_area(this, i, j) result(area)
      class(GridGeometryType), intent(in) :: this
      integer, intent(in) :: i, j
      real(fp) :: area

      if (allocated(this%grid_area)) then
         area = this%grid_area(i, j)
      else
         area = this%dx * this%dy  ! Default rectangular area
      endif
   end function geometry_get_cell_area

   function geometry_get_cell_volume(this, i, j, k) result(volume)
      class(GridGeometryType), intent(in) :: this
      integer, intent(in) :: i, j, k
      real(fp) :: volume

      real(fp) :: area

      area = this%get_cell_area(i, j)
      volume = area * this%dz(k)
   end function geometry_get_cell_volume

   function geometry_is_valid_position(this, i, j) result(is_valid)
      class(GridGeometryType), intent(in) :: this
      integer, intent(in) :: i, j
      logical :: is_valid

      is_valid = (i >= 1 .and. i <= this%nx .and. j >= 1 .and. j <= this%ny)
   end function geometry_is_valid_position

   !========================================================================
   ! GridDecompositionType Implementation
   !========================================================================

   subroutine decomp_init(this, geometry, n_procs, my_rank, rc)
      class(GridDecompositionType), intent(inout) :: this
      type(GridGeometryType), intent(in) :: geometry
      integer, intent(in) :: n_procs, my_rank
      integer, intent(out) :: rc

      integer :: nx, ny, nz
      integer :: cols_per_proc, extra_cols

      rc = cc_success

      this%n_procs = n_procs
      this%my_rank = my_rank

      call geometry%get_dimensions(nx, ny, nz)

      ! Simple 1D decomposition along x-direction for now
      cols_per_proc = nx / n_procs
      extra_cols = mod(nx, n_procs)

      this%i_start = my_rank * cols_per_proc + 1 + min(my_rank, extra_cols)
      this%i_end = this%i_start + cols_per_proc - 1
      if (my_rank < extra_cols) this%i_end = this%i_end + 1

      this%j_start = 1
      this%j_end = ny

   end subroutine decomp_init

   subroutine decomp_get_local_bounds(this, i_start, i_end, j_start, j_end)
      class(GridDecompositionType), intent(in) :: this
      integer, intent(out) :: i_start, i_end, j_start, j_end

      i_start = this%i_start
      i_end = this%i_end
      j_start = this%j_start
      j_end = this%j_end
   end subroutine decomp_get_local_bounds

   function decomp_is_local_column(this, i, j) result(is_local)
      class(GridDecompositionType), intent(in) :: this
      integer, intent(in) :: i, j
      logical :: is_local

      is_local = (i >= this%i_start .and. i <= this%i_end .and. &
         j >= this%j_start .and. j <= this%j_end)
   end function decomp_is_local_column

   !========================================================================
   ! ColumnIteratorType Implementation
   !========================================================================

   subroutine iterator_init(this, grid_mgr, rc)
      class(ColumnIteratorType), intent(inout) :: this
      type(GridManagerType), intent(in), target :: grid_mgr
      integer, intent(out) :: rc

      rc = cc_success

      this%grid_mgr => grid_mgr
      this%total_columns = grid_mgr%get_total_columns()
      call this%reset()
   end subroutine iterator_init

   function iterator_has_next(this) result(has_next)
      class(ColumnIteratorType), intent(in) :: this
      logical :: has_next

      has_next = (this%current_column < this%total_columns)
   end function iterator_has_next

   subroutine iterator_next(this, rc)
      class(ColumnIteratorType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: nx, ny, nz

      rc = cc_success

      if (.not. this%has_next()) then
         rc = cc_failure
         return
      endif

      this%current_column = this%current_column + 1

      ! Convert linear index to 2D indices
      call this%grid_mgr%geometry%get_dimensions(nx, ny, nz)
      this%current_j = (this%current_column - 1) / nx + 1
      this%current_i = this%current_column - (this%current_j - 1) * nx
   end subroutine iterator_next

   function iterator_get_current_column(this) result(column_view)
      class(ColumnIteratorType), intent(in) :: this
      type(ColumnViewType) :: column_view
      ! Returns the column view for the current (i,j) indices. If indices are invalid, returns an invalid column_view.
      column_view = this%grid_mgr%get_column_by_indices(this%current_i, this%current_j)
   end function iterator_get_current_column

   subroutine iterator_get_current_indices(this, i, j)
      class(ColumnIteratorType), intent(in) :: this
      integer, intent(out) :: i, j

      i = this%current_i
      j = this%current_j
   end subroutine iterator_get_current_indices

   subroutine iterator_reset(this)
      class(ColumnIteratorType), intent(inout) :: this

      this%current_column = 0
      this%current_i = 1
      this%current_j = 1
   end subroutine iterator_reset

   !========================================================================
   ! GridManagerType Implementation
   !========================================================================

   subroutine grid_manager_init(this, nx, ny, nz, error_mgr, grid_type, coord_system, rc)
      class(GridManagerType), intent(inout) :: this
      integer, intent(in) :: nx, ny, nz
      type(ErrorManagerType), target, intent(inout) :: error_mgr
      integer, intent(in), optional :: grid_type, coord_system
      integer, intent(out) :: rc

      character(len=256) :: thisLoc
      integer :: local_rc

      thisloc = 'grid_manager_init (in core/GridManager_Mod.F90)'
      rc = cc_success

      this%error_mgr => error_mgr

      ! Initialize geometry
      call this%geometry%init(nx, ny, nz, grid_type, coord_system, local_rc)
      if (local_rc /= cc_success) then
         call error_mgr%report_error(error_process_initialization, &
            'Failed to initialize grid geometry', rc, thisloc)
         return
      endif

      ! Initialize decomposition (default: single processor)
      call this%decomp%init(this%geometry, 1, 0, local_rc)
      if (local_rc /= cc_success) then
         call error_mgr%report_error(error_process_initialization, &
            'Failed to initialize grid decomposition', rc, thisloc)
         return
      endif

      ! Calculate number of local columns
      this%n_local_columns = (this%decomp%i_end - this%decomp%i_start + 1) * &
         (this%decomp%j_end - this%decomp%j_start + 1)

      this%is_initialized = .true.

   end subroutine grid_manager_init

   subroutine grid_manager_cleanup(this)
      class(GridManagerType), intent(inout) :: this

      call this%geometry%cleanup()

      this%is_initialized = .false.
      this%error_mgr => null()
   end subroutine grid_manager_cleanup

   function grid_manager_validate(this) result(is_valid)
      class(GridManagerType), intent(in) :: this
      logical :: is_valid

      is_valid = this%is_initialized .and. &
         this%geometry%is_initialized .and. &
         associated(this%error_mgr)
   end function grid_manager_validate

   function grid_manager_get_geometry(this) result(geometry)
      class(GridManagerType), intent(in) :: this
      type(GridGeometryType) :: geometry

      geometry = this%geometry
   end function grid_manager_get_geometry

   function grid_manager_get_decomposition(this) result(decomp)
      class(GridManagerType), intent(in) :: this
      type(GridDecompositionType) :: decomp

      decomp = this%decomp
   end function grid_manager_get_decomposition

   function grid_manager_get_total_columns(this) result(total_columns)
      class(GridManagerType), intent(in) :: this
      integer :: total_columns

      total_columns = this%geometry%nx * this%geometry%ny
   end function grid_manager_get_total_columns

   function grid_manager_get_local_columns(this) result(local_columns)
      class(GridManagerType), intent(in) :: this
      integer :: local_columns

      local_columns = this%n_local_columns
   end function grid_manager_get_local_columns

   subroutine grid_manager_get_shape(this, nx, ny, nz)
      class(GridManagerType), intent(in) :: this
      integer, intent(out) :: nx, ny, nz
      call this%geometry%get_dimensions(nx, ny, nz)
   end subroutine grid_manager_get_shape

   function grid_manager_get_column_by_indices(this, i, j) result(column_view)
      class(GridManagerType), intent(in) :: this
      integer, intent(in) :: i, j
      type(ColumnViewType) :: column_view
      integer :: rc

      ! Check bounds for safety
      if (.not. this%geometry%is_valid_position(i, j)) then
         ! Return an uninitialized column_view to mark as invalid
         ! The column_view will have its default null() pointers
         return
      endif

      ! For now, return an uninitialized column_view
      ! TODO: Implement proper column_view initialization when StateManager integration is complete
      ! call column_view%init(container, rc)
   end function grid_manager_get_column_by_indices

   function grid_manager_get_column_by_location(this, lon, lat) result(column_view)
      class(GridManagerType), intent(in) :: this
      real(fp), intent(in) :: lon, lat
      type(ColumnViewType) :: column_view

      ! This would find the nearest grid cell to the given lon/lat
      ! and return the corresponding column view
      ! Implementation depends on coordinate system and grid projection
   end function grid_manager_get_column_by_location

   function grid_manager_create_column_iterator(this) result(iterator)
      class(GridManagerType), intent(in) :: this
      type(ColumnIteratorType) :: iterator

      integer :: rc
      call iterator%init(this, rc)
   end function grid_manager_create_column_iterator

   subroutine grid_manager_compute_distances(this, i1, j1, i2, j2, distance, rc)
      class(GridManagerType), intent(in) :: this
      integer, intent(in) :: i1, j1, i2, j2
      real(fp), intent(out) :: distance
      integer, intent(out) :: rc

      real(fp) :: dx, dy

      rc = cc_success

      ! Simple Euclidean distance for Cartesian grids
      dx = real(i2 - i1, fp) * this%geometry%dx
      dy = real(j2 - j1, fp) * this%geometry%dy
      distance = sqrt(dx*dx + dy*dy)

      ! For lon/lat grids, would use great circle distance
   end subroutine grid_manager_compute_distances

   subroutine grid_manager_interpolate_to_column(this, source_data, target_i, target_j, &
      interpolated_data, rc)
      class(GridManagerType), intent(in) :: this
      real(fp), intent(in) :: source_data(:,:)
      integer, intent(in) :: target_i, target_j
      real(fp), intent(out) :: interpolated_data
      integer, intent(out) :: rc

      rc = cc_success

      ! Simple nearest neighbor for now
      if (this%geometry%is_valid_position(target_i, target_j)) then
         interpolated_data = source_data(target_i, target_j)
      else
         rc = cc_failure
      endif
   end subroutine grid_manager_interpolate_to_column

   subroutine grid_manager_exchange_halo_data(this, data, rc)
      class(GridManagerType), intent(in) :: this
      real(fp), intent(inout) :: data(:,:,:)
      integer, intent(out) :: rc

      rc = cc_success
      ! Placeholder for MPI halo exchange implementation
   end subroutine grid_manager_exchange_halo_data

   subroutine grid_manager_print_info(this)
      class(GridManagerType), intent(in) :: this

      integer :: nx, ny, nz

      call this%geometry%get_dimensions(nx, ny, nz)

      write(*,'(A)') '==============================================='
      write(*,'(A)') 'Grid Manager Information'
      write(*,'(A)') '==============================================='
      write(*,'(A,A)') 'Name: ', trim(this%name)
      write(*,'(A,L1)') 'Initialized: ', this%is_initialized
      write(*,'(A,3I6)') 'Dimensions (nx,ny,nz): ', nx, ny, nz
      write(*,'(A,I0)') 'Grid type: ', this%geometry%grid_type
      write(*,'(A,I0)') 'Total columns: ', this%get_total_columns()
      write(*,'(A,I0)') 'Local columns: ', this%get_local_columns()
      write(*,'(A)') '==============================================='
   end subroutine grid_manager_print_info

   function grid_manager_is_ready(this) result(is_ready)
      class(GridManagerType), intent(in) :: this
      logical :: is_ready

      is_ready = this%validate()
   end function grid_manager_is_ready

   subroutine grid_manager_create_virtual_column(this, i, j, virtual_col, rc)
      class(GridManagerType), intent(in) :: this
      integer, intent(in) :: i, j
      type(VirtualColumnType), intent(out) :: virtual_col
      integer, intent(out) :: rc

      ! TODO: Implement actual virtual column construction
      rc = cc_success
   end subroutine grid_manager_create_virtual_column

end module gridmanager_mod
```


