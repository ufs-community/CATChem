

# File VirtualColumn\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**VirtualColumn\_Mod.F90**](_virtual_column___mod_8_f90.md)

[Go to the documentation of this file](_virtual_column___mod_8_f90.md)


```Fortran

module virtualcolumn_mod
   use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure

   implicit none
   private

   public :: virtualcolumntype, virtualmettype

#include "virtualmet_type.inc"

   type :: virtualcolumntype
      ! Meteorological data (pointers managed by VirtualMetType)
      type(VirtualMetType) :: met

      ! Chemical and emission data (still use arrays for modification)
      real(fp), allocatable :: chem_data(:,:)
      real(fp), allocatable :: emis_data(:,:)

      ! Grid position and metadata
      integer :: grid_i = 0
      integer :: grid_j = 0
      real(fp) :: lat = 0.0_fp                               
      real(fp) :: lon = 0.0_fp                               
      real(fp) :: area = 0.0_fp                              

      ! Dimensions
      integer :: nlev = 0
      integer :: nspec_chem = 0
      integer :: nspec_emis = 0

      ! Status
      logical :: is_valid = .false.                          

   contains
      procedure :: init => virtual_column_init
      procedure :: get_met => virtual_column_get_met
      procedure :: get_chem_field => virtual_column_get_chem_field
      procedure :: get_emis_field => virtual_column_get_emis_field
      procedure :: set_chem_field => virtual_column_set_chem_field
      procedure :: set_emis_field => virtual_column_set_emis_field
      procedure :: get_position => virtual_column_get_position
      procedure :: get_metadata => virtual_column_get_metadata
      procedure :: get_dimensions => virtual_column_get_dimensions
      procedure :: is_initialized => virtual_column_is_initialized
      procedure :: cleanup => virtual_column_cleanup
   end type virtualcolumntype

contains

   !=========================================================================
   ! VirtualMetType Implementation
   !=========================================================================

   subroutine virtual_met_cleanup(this)
      class(VirtualMetType), intent(inout) :: this

      !print *, '[DEBUG] Entering virtual_met_cleanup'
      if (associated(this%T)) then
         !   print *, '[DEBUG] Cleaning up T, associated before nullify'
      else
         !   print *, '[DEBUG] T not associated'
      endif

      ! Generated cleanup code from MetState field definitions
#include "virtualmet_cleanup.inc"

      !print *, '[DEBUG] Exiting virtual_met_cleanup'
   end subroutine virtual_met_cleanup

   !=========================================================================
   ! VirtualColumnType Implementation
   !=========================================================================

   subroutine virtual_column_init(this, nlev, nspec_chem, nspec_emis, grid_i, grid_j, lat, lon, area, rc)
      class(VirtualColumnType), intent(inout) :: this
      integer, intent(in) :: nlev, nspec_chem, nspec_emis
      integer, intent(in) :: grid_i, grid_j
      real(fp), intent(in) :: lat, lon, area
      integer, intent(out) :: rc

      rc = cc_success

      ! Store dimensions
      this%nlev = nlev
      this%nspec_chem = nspec_chem
      this%nspec_emis = nspec_emis

      ! Store position/metadata
      this%grid_i = grid_i
      this%grid_j = grid_j
      this%lat = lat
      this%lon = lon
      this%area = area

      ! Allocate meteorological field arrays for testing
      if (nlev > 0) then
         if (.not. associated(this%met%T)) then
            allocate(this%met%T(nlev))
            this%met%T = 288.15_fp  ! Initialize with default temperature value
         endif
      endif

      ! Allocate chemical and emission data arrays
      if (nlev > 0 .and. nspec_chem > 0) then
         allocate(this%chem_data(nlev, nspec_chem), stat=rc)
         if (rc /= 0) then
            rc = cc_failure
            return
         endif
         this%chem_data = 0.0_fp
      endif

      if (nlev > 0 .and. nspec_emis > 0) then
         allocate(this%emis_data(nlev, nspec_emis), stat=rc)
         if (rc /= 0) then
            rc = cc_failure
            return
         endif
         this%emis_data = 0.0_fp
      endif

      this%is_valid = .true.
      rc = cc_success
   end subroutine virtual_column_init

   function virtual_column_get_met(this) result(met_ptr)
      class(VirtualColumnType), intent(in), target :: this
      type(VirtualMetType), pointer :: met_ptr

      met_ptr => this%met
   end function virtual_column_get_met

   function virtual_column_get_chem_field(this, ispec, k) result(value)
      class(VirtualColumnType), intent(in) :: this
      integer, intent(in) :: ispec, k
      real(fp) :: value

      if (allocated(this%chem_data) .and. &
         k >= 1 .and. k <= this%nlev .and. &
         ispec >= 1 .and. ispec <= this%nspec_chem) then
         value = this%chem_data(k, ispec)
      else
         value = 0.0_fp
      endif
   end function virtual_column_get_chem_field

   function virtual_column_get_emis_field(this, ispec, k) result(value)
      class(VirtualColumnType), intent(in) :: this
      integer, intent(in) :: ispec, k
      real(fp) :: value

      if (allocated(this%emis_data) .and. &
         k >= 1 .and. k <= this%nlev .and. &
         ispec >= 1 .and. ispec <= this%nspec_emis) then
         value = this%emis_data(k, ispec)
      else
         value = 0.0_fp
      endif
   end function virtual_column_get_emis_field

   subroutine virtual_column_set_chem_field(this, k, ispec, value)
      class(VirtualColumnType), intent(inout) :: this
      integer, intent(in) :: k, ispec
      real(fp), intent(in) :: value

      if (allocated(this%chem_data) .and. &
         k >= 1 .and. k <= this%nlev .and. &
         ispec >= 1 .and. ispec <= this%nspec_chem) then
         this%chem_data(k, ispec) = value
      endif
   end subroutine virtual_column_set_chem_field

   subroutine virtual_column_set_emis_field(this, k, ispec, value)
      class(VirtualColumnType), intent(inout) :: this
      integer, intent(in) :: k, ispec
      real(fp), intent(in) :: value

      if (allocated(this%emis_data) .and. &
         k >= 1 .and. k <= this%nlev .and. &
         ispec >= 1 .and. ispec <= this%nspec_emis) then
         this%emis_data(k, ispec) = value
      endif
   end subroutine virtual_column_set_emis_field

   subroutine virtual_column_get_position(this, grid_i, grid_j)
      class(VirtualColumnType), intent(in) :: this
      integer, intent(out) :: grid_i, grid_j

      grid_i = this%grid_i
      grid_j = this%grid_j
   end subroutine virtual_column_get_position

   subroutine virtual_column_get_metadata(this, lat, lon, area)
      class(VirtualColumnType), intent(in) :: this
      real(fp), intent(out) :: lat, lon, area

      lat = this%lat
      lon = this%lon
      area = this%area
   end subroutine virtual_column_get_metadata

   subroutine virtual_column_get_dimensions(this, nlev, nspec_chem, nspec_emis)
      class(VirtualColumnType), intent(in) :: this
      integer, intent(out) :: nlev, nspec_chem, nspec_emis

      nlev = this%nlev
      nspec_chem = this%nspec_chem
      nspec_emis = this%nspec_emis
   end subroutine virtual_column_get_dimensions

   function virtual_column_is_initialized(this) result(initialized)
      class(VirtualColumnType), intent(in) :: this
      logical :: initialized

      initialized = this%is_valid
   end function virtual_column_is_initialized

   subroutine virtual_column_cleanup(this)
      class(VirtualColumnType), intent(inout) :: this

      !print *, '[DEBUG] Entering virtual_column_cleanup'

      ! Clean up meteorological pointers
      call this%met%cleanup()

      ! Deallocate chemical and emission data
      if (allocated(this%chem_data)) then
         !print *, '[DEBUG] Deallocating chem_data'
         deallocate(this%chem_data)
      endif
      if (allocated(this%emis_data)) then
         !print *, '[DEBUG] Deallocating emis_data'
         deallocate(this%emis_data)
      endif

      this%nlev = 0
      this%nspec_chem = 0
      this%nspec_emis = 0
      this%is_valid = .false.

      !print *, '[DEBUG] Exiting virtual_column_cleanup'
   end subroutine virtual_column_cleanup

end module virtualcolumn_mod
```


