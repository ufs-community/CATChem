

# File ColumnInterface\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**ColumnInterface\_Mod.F90**](_column_interface___mod_8_f90.md)

[Go to the documentation of this file](_column_interface___mod_8_f90.md)


```Fortran
!========================================================================
!! Virtual Column Interface Module
!!
!! Provides a column-based interface to multi-dimensional atmospheric data.
!! This allows column-based processes to work transparently with both 1D
!! (column) and multi-dimensional (2D/3D) host models.
!!
!! Key features:
!! - Zero-copy access via pointers
!! - OpenMP/OpenACC friendly
!! - Support for 1D, 2D, and 3D grids
!!
!! @author CATChem Development Team
!! @date 2025
!! @version 1.0
!!
!! This module provides a column-based interface to multi-dimensional atmospheric data,
!! allowing column-based processes to work transparently with both 1D (column) and
!! multi-dimensional (2D/3D) host models.
!!
!! @details
!! - Zero-copy access via pointers
!! - OpenMP/OpenACC friendly
!! - Support for 1D, 2D, and 3D grids
!!
module columninterface_mod
   use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure
   use chemstate_mod, only: chemstatetype
   use metstate_mod, only: metstatetype
   use gridgeometry_mod, only: gridgeometrytype
   use virtualcolumn_mod, only: virtualcolumntype
   ! use EmisState_Mod, only: EmisStateType

   implicit none
   private

   ! Forward declaration for GridManagerType to avoid circular dependency
   type :: gridmanagertype
   end type gridmanagertype

   public :: columnviewtype
   public :: create_column_view
   ! public :: VirtualColumnType
   public :: columnprocessortype

   !========================================================================
   !! Column View Type
   !!
   !! Provides column-based access to multi-dimensional state data.
   !! Supports both direct column models and virtual columns from 3D data.
   !========================================================================

   type :: columnviewtype
      private

      ! Grid information
      integer :: grid_type
      integer :: nx
      integer :: ny
      integer :: nlev
      integer :: current_i
      integer :: current_j

      ! State object references (optional - can be null for simple column access)
      type(MetStateType), pointer :: met_state => null()      
      type(ChemStateType), pointer :: chem_state => null()    

      ! Column data pointers (for current column)
      real(fp), pointer :: met_column(:) => null()            
      real(fp), pointer :: chem_column(:,:) => null()         
      real(fp), pointer :: emis_column(:,:) => null()         

   contains
      procedure :: init => column_view_init
      procedure :: set_column_position => column_view_set_position
      procedure :: get_met_column_ptr => column_view_get_met_ptr
      procedure :: get_chem_column_ptr => column_view_get_chem_ptr
      procedure :: get_emis_column_ptr => column_view_get_emis_ptr
      procedure :: apply_column_tendency => column_view_apply_tendency
      procedure :: cleanup => column_view_cleanup

   end type columnviewtype

   type :: columnprocessortype
      private

      type(VirtualColumnType), allocatable :: columns(:)
      integer :: n_columns = 0
      integer :: current_column = 1

   contains
      procedure :: init => processor_init
      procedure :: add_column => processor_add_column
      procedure :: process_all => processor_process_all
      procedure :: get_column => processor_get_column
      procedure :: get_current_column => processor_get_current_column
      procedure :: next_column => processor_next_column
      procedure :: reset => processor_reset
      procedure :: cleanup => processor_cleanup
   end type columnprocessortype

   abstract interface
      subroutine column_process_interface(column, rc)
         import :: virtualcolumntype
         type(VirtualColumnType), intent(inout) :: column
         integer, intent(out) :: rc
      end subroutine column_process_interface
   end interface

contains

   !========================================================================
   !! Initialize column view with state objects
   !========================================================================
   subroutine column_view_init(this, met_state, chem_state, rc)
      class(ColumnViewType), intent(inout) :: this
      type(MetStateType), intent(in), target, optional :: met_state
      type(ChemStateType), intent(in), target, optional :: chem_state
      integer, intent(out) :: rc

      rc = cc_success

      ! Store state references if provided
      if (present(met_state)) this%met_state => met_state
      if (present(chem_state)) this%chem_state => chem_state

      ! Get grid dimensions from MetState if available
      if (associated(this%met_state)) then
         ! For now, use placeholder dimensions
         this%nx = 1
         this%ny = 1
         this%nlev = 50
      else
         ! Use default dimensions
         this%nx = 1
         this%ny = 1
         this%nlev = 50
      endif

      if (this%nx == 1 .and. this%ny == 1) then
         this%grid_type = 1  ! Pure column model
      else if (this%ny == 1) then
         this%grid_type = 2  ! 2D (x-z) model
      else
         this%grid_type = 3  ! 3D (x-y-z) model
      endif

      ! Initialize position
      this%current_i = 1
      this%current_j = 1

   end subroutine column_view_init

   !========================================================================
   !! Set current column position (for multi-dimensional grids)
   !========================================================================
   subroutine column_view_set_position(this, i, j, rc)
      class(ColumnViewType), intent(inout) :: this
      integer, intent(in) :: i, j
      integer, intent(out) :: rc

      rc = cc_success

      ! Validate position
      if (i < 1 .or. i > this%nx .or. j < 1 .or. j > this%ny) then
         print *, 'ERROR: Invalid column position:', i, j
         rc = cc_failure
         return
      endif

      this%current_i = i
      this%current_j = j

      ! Update pointers to current column - TODO: implement if needed
      rc = cc_success

   end subroutine column_view_set_position

   !========================================================================
   !! Get meteorological column pointer
   !========================================================================
   function column_view_get_met_ptr(this, field_name) result(column_ptr)
      class(ColumnViewType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      real(fp), pointer :: column_ptr(:)

      column_ptr => null()

      if (.not. associated(this%met_state)) return

      ! Real implementation would get field from MetState
      ! For now, return null pointer
      column_ptr => null()

   end function column_view_get_met_ptr

   !========================================================================
   !! Get chemistry column pointer for a species
   !========================================================================
   function column_view_get_chem_ptr(this, species_name) result(column_ptr)
      class(ColumnViewType), intent(in) :: this
      character(len=*), intent(in) :: species_name
      real(fp), pointer :: column_ptr(:)

      column_ptr => null()

      if (.not. associated(this%chem_state)) return

      ! Real implementation would get species from ChemState
      ! For now, return null pointer
      column_ptr => null()

   end function column_view_get_chem_ptr

   !========================================================================
   !! Get emissions column pointer for a species
   !========================================================================
   function column_view_get_emis_ptr(this, species_name) result(column_ptr)
      class(ColumnViewType), intent(in) :: this
      character(len=*), intent(in) :: species_name
      real(fp), pointer :: column_ptr(:)

      ! For now, return null pointer since EmisState is commented out
      column_ptr => null()

   end function column_view_get_emis_ptr

   !========================================================================
   !! Apply tendency to current column and accumulate in emissions
   !========================================================================
   subroutine column_view_apply_tendency(this, species_name, tendency_column, &
      dt, layer_height, rc)
      class(ColumnViewType), intent(inout) :: this
      character(len=*), intent(in) :: species_name
      real(fp), intent(in) :: tendency_column(:)  ! Tendency flux (kg/m²/s)
      real(fp), intent(in) :: dt
      real(fp), intent(in), optional :: layer_height
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state
      real(fp), pointer :: chem_column(:), emis_column(:)
      real(fp) :: height
      integer :: k

      rc = cc_success

      height = 50.0_fp  ! Default surface layer height
      if (present(layer_height)) height = layer_height

      ! Get pointers to current column data
      chem_column => this%get_chem_column_ptr(species_name)
      emis_column => this%get_emis_column_ptr(species_name)

      if (.not. associated(chem_column) .or. .not. associated(emis_column)) then
         print *, 'ERROR: Could not get column pointers for ', trim(species_name)
         rc = cc_failure
         return
      endif

      ! Get chemistry state for species metadata
      ! Real implementation would use this%chem_state directly
      ! For now, skip this step

      ! Apply tendency to each vertical level
      do k = 1, this%nlev
         ! Apply tendency directly to chemistry column
         chem_column(k) = chem_column(k) + tendency_column(k) * dt

         ! Accumulate in emissions for diagnostics
         emis_column(k) = emis_column(k) + tendency_column(k)
      end do

   end subroutine column_view_apply_tendency

   !========================================================================
   !! Cleanup column view
   !========================================================================
   subroutine column_view_cleanup(this)
      class(ColumnViewType), intent(inout) :: this

      ! Nullify pointers
      this%met_state => null()
      this%chem_state => null()
      this%met_column => null()
      this%chem_column => null()
      this%emis_column => null()

   end subroutine column_view_cleanup

   !========================================================================
   !! Convenience function to create column view
   !========================================================================
   function create_column_view(met_state, chem_state) result(col_view)
      type(MetStateType), intent(in), target, optional :: met_state
      type(ChemStateType), intent(in), target, optional :: chem_state
      type(ColumnViewType) :: col_view

      integer :: rc

      call col_view%init(met_state, chem_state, rc)
      if (rc /= cc_success) then
         print *, 'WARNING: Failed to initialize column view'
      endif

   end function create_column_view

   !========================================================================
   ! ColumnProcessorType Implementation
   !========================================================================

   subroutine processor_init(this, max_columns, rc)
      class(ColumnProcessorType), intent(inout) :: this
      integer, intent(in) :: max_columns
      integer, intent(out) :: rc

      rc = cc_success

      allocate(this%columns(max_columns))
      this%n_columns = 0
      this%current_column = 1

   end subroutine processor_init

   subroutine processor_add_column(this, virtual_col, rc)
      class(ColumnProcessorType), intent(inout) :: this
      type(VirtualColumnType), intent(in) :: virtual_col
      integer, intent(out) :: rc

      rc = cc_success

      if (this%n_columns >= size(this%columns)) then
         rc = cc_failure
         return
      endif

      this%n_columns = this%n_columns + 1
      this%columns(this%n_columns) = virtual_col

   end subroutine processor_add_column

   subroutine processor_process_all(this, process_proc, rc)
      class(ColumnProcessorType), intent(inout) :: this
      procedure(column_process_interface) :: process_proc
      integer, intent(out) :: rc

      integer :: i, local_rc

      rc = cc_success

      do i = 1, this%n_columns
         call process_proc(this%columns(i), local_rc)
         if (local_rc /= cc_success) then
            rc = local_rc
            ! Continue processing other columns or exit based on error policy
         endif
      enddo

   end subroutine processor_process_all

   function processor_get_column(this, index) result(column)
      class(ColumnProcessorType), intent(in) :: this
      integer, intent(in) :: index
      type(VirtualColumnType) :: column

      if (index >= 1 .and. index <= this%n_columns) then
         column = this%columns(index)
      endif

   end function processor_get_column

   function processor_get_current_column(this) result(column)
      class(ColumnProcessorType), intent(in) :: this
      type(VirtualColumnType) :: column

      column = this%get_column(this%current_column)

   end function processor_get_current_column

   subroutine processor_next_column(this, rc)
      class(ColumnProcessorType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success

      if (this%current_column < this%n_columns) then
         this%current_column = this%current_column + 1
      else
         rc = cc_failure  ! No more columns
      endif

   end subroutine processor_next_column

   subroutine processor_reset(this)
      class(ColumnProcessorType), intent(inout) :: this

      this%current_column = 1

   end subroutine processor_reset

   subroutine processor_cleanup(this)
      class(ColumnProcessorType), intent(inout) :: this

      integer :: i

      if (allocated(this%columns)) then
         do i = 1, this%n_columns
            call this%columns(i)%cleanup()
         enddo
         deallocate(this%columns)
      endif

      this%n_columns = 0

   end subroutine processor_cleanup

end module columninterface_mod
```


