

# File ExtEmisData\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**ExtEmisData\_Mod.F90**](_ext_emis_data___mod_8_f90.md)

[Go to the documentation of this file](_ext_emis_data___mod_8_f90.md)


```Fortran

MODULE extemisdata_mod

   !=========================================================================
   ! Module uses
   !=========================================================================
   USE error_mod
   ! USE logging_mod, only: log_message, LOG_ERROR, LOG_WARNING, LOG_INFO, LOG_DEBUG
   USE precision_mod
   USE species_mod, only: speciestype
   USE ieee_arithmetic

   IMPLICIT NONE
   PRIVATE

   !=========================================================================
   ! Public interfaces
   !=========================================================================
   !PUBLIC :: ExtEmisDataType
   !PUBLIC :: ExtEmisFieldType
   !PUBLIC :: ExtEmisCategoryType

   TYPE, PUBLIC :: extemisfieldtype
      CHARACTER(LEN=64)             :: field_name = ''
      CHARACTER(LEN=128)            :: long_name = ''
      CHARACTER(LEN=32)             :: units = ''
      !CHARACTER(LEN=256)            :: source_file = ''    !< Source file path
      INTEGER                       :: nx = 0
      INTEGER                       :: ny = 0
      INTEGER                       :: nz = 1
      REAL(fp)                      :: factors = 1.0_fp    
      REAL(fp), ALLOCATABLE         :: lat(:)
      REAL(fp), ALLOCATABLE         :: lon(:)
      REAL(fp), ALLOCATABLE         :: stkdm(:)
      REAL(fp), ALLOCATABLE         :: stkht(:)
      REAL(fp), ALLOCATABLE         :: stktk(:)
      REAL(fp), ALLOCATABLE         :: stkve(:)
      INTEGER,  ALLOCATABLE         :: ip(:)
      INTEGER,  ALLOCATABLE         :: jp(:)
      INTEGER,  ALLOCATABLE         :: ijmap(:)
      INTEGER                       :: n_times = 0
      INTEGER                       :: current_time_idx = 1
      LOGICAL                       :: time_interpolate = .true. 
      LOGICAL                       :: diagnostic = .false. 
      REAL(fp), ALLOCATABLE         :: emission_data(:,:,:,:)
      ! REAL(fp), ALLOCATABLE         :: longitude(:,:)         !< Longitude coordinates [degrees]
      ! REAL(fp), ALLOCATABLE         :: latitude(:,:)          !< Latitude coordinates [degrees]
      ! REAL(fp), ALLOCATABLE         :: vertical(:,:)          !< Vertical coordinates (if applicable)
      ! REAL(fp), ALLOCATABLE         :: time_coords(:,:)       !< Time coordinates
      LOGICAL                       :: is_loaded = .false. 
      LOGICAL                       :: is_valid = .false.  
      CHARACTER(LEN=32)             :: interpolation_method = 'bilinear'
   CONTAINS
      PROCEDURE :: init => extemifield_init
      PROCEDURE :: cleanup => extemifield_cleanup
      PROCEDURE :: validate => extemifield_validate
      PROCEDURE :: get_emission_at_point => extemifield_get_emission_at_point
      PROCEDURE :: load_from_file => extemifield_load_from_file
      PROCEDURE :: get_column_ptr => extemifield_get_column_ptr
   END TYPE extemisfieldtype

   TYPE, PUBLIC :: extemiscategorytype
      CHARACTER(LEN=64)                         :: category_name = ''
      CHARACTER(LEN=256)                        :: description = ''
      INTEGER                                   :: n_fields = 0
      INTEGER                                   :: irec = 0
      TYPE(ExtEmisFieldType), ALLOCATABLE       :: fields(:)
      LOGICAL                                   :: is_active = .true.  
      LOGICAL                                   :: gridded = .true.    
      LOGICAL                                   :: diagnostic = .true.  
      REAL(fp)                                  :: global_scale = 1.0_fp 
      REAL(fp)                                  :: topfraction = -1.0_fp 
      CHARACTER(LEN=128)                        :: source_file = ''
      CHARACTER(LEN=128)                        :: format = ''
      CHARACTER(LEN=128)                        :: frequency = ''
      CHARACTER(LEN=128)                        :: latname = ''
      CHARACTER(LEN=128)                        :: lonname = ''
      CHARACTER(LEN=128)                        :: stkdmname = ''
      CHARACTER(LEN=128)                        :: stkhtname = ''
      CHARACTER(LEN=128)                        :: stktkname = ''
      CHARACTER(LEN=128)                        :: stkvename = ''
      CHARACTER(LEN=128)                        :: plumerise = ''


   CONTAINS
      PROCEDURE :: init => extemicat_init
      PROCEDURE :: cleanup => extemicat_cleanup
      PROCEDURE :: validate => extemicat_validate
      PROCEDURE :: add_field => extemicat_add_field
      PROCEDURE :: find_field => extemicat_find_field
   END TYPE extemiscategorytype

   TYPE, PUBLIC :: extemisdatatype
      INTEGER                                   :: n_categories = 0
      TYPE(ExtEmisCategoryType), ALLOCATABLE    :: categories(:)
      LOGICAL                                   :: is_active = .true.  
      LOGICAL                                   :: diagnostic = .true.  
      INTEGER                                   :: total_fields = 0
      CHARACTER(LEN=128)                        :: data_source = ''
      REAL(fp)                                  :: global_scale = 1.0_fp 
   CONTAINS
      PROCEDURE :: init => extemidata_init
      PROCEDURE :: cleanup => extemidata_cleanup
      PROCEDURE :: validate => extemidata_validate
      PROCEDURE :: add_category => extemidata_add_category
      PROCEDURE :: load_emission_files => extemidata_load_files
      PROCEDURE :: find_emission_field => extemidata_find_field
      PROCEDURE :: get_emission_rate => extemidata_get_emission_rate
      PROCEDURE :: update_time => extemidata_update_time
      PROCEDURE :: get_memory_usage => extemidata_get_memory_usage
      PROCEDURE :: get_column_ptr => extemidata_get_column_ptr
   END TYPE extemisdatatype

CONTAINS

   !=========================================================================
   ! ExtEmisFieldType procedures
   !=========================================================================

   subroutine extemifield_init(this, field_name, nx, ny, nz, n_times, units, rc)
      implicit none
      class(ExtEmisFieldType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: nx, ny
      integer, intent(in), optional :: nz, n_times
      character(len=*), intent(in), optional :: units
      integer, intent(out) :: rc

      rc = cc_success

      this%field_name = trim(field_name)
      this%nx = nx
      this%ny = ny

      if (present(nz)) then
         this%nz = nz
      else
         this%nz = 1
      endif

      if (present(n_times)) then
         this%n_times = n_times
      else
         this%n_times = 1
      endif

      if (present(units)) then
         this%units = trim(units)
      else
         this%units = 'kg/m2/s'
      endif

      if (allocated(this%emission_data)) deallocate(this%emission_data)
      allocate(this%emission_data(this%nx, this%ny, this%nz, this%n_times))

      ! Initialize emission data to zero to avoid garbage values
      this%emission_data = 0.0_fp

      this%current_time_idx = 1
      this%factors = 1.0_fp
      this%time_interpolate = .true.
      this%diagnostic = .false.
      this%is_loaded = .false.
      this%is_valid = .false.

   end subroutine extemifield_init

   subroutine extemifield_cleanup(this, rc)
      implicit none
      class(ExtEmisFieldType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success

      if (allocated(this%emission_data)) deallocate(this%emission_data)
      if (allocated(this%lat)) deallocate(this%lat)
      if (allocated(this%lon)) deallocate(this%lon)
      if (allocated(this%stkdm)) deallocate(this%stkdm)
      if (allocated(this%stkht)) deallocate(this%stkht)
      if (allocated(this%stktk)) deallocate(this%stktk)
      if (allocated(this%stkve)) deallocate(this%stkve)
      if (allocated(this%ip)) deallocate(this%ip)
      if (allocated(this%jp)) deallocate(this%jp)
      if (allocated(this%ijmap)) deallocate(this%ijmap)

      this%field_name = ''
      this%long_name = ''
      this%units = ''
      this%nx = 0
      this%ny = 0
      this%nz = 1
      this%factors = 1.0_fp
      this%n_times = 0
      this%current_time_idx = 1
      this%time_interpolate = .true.
      this%diagnostic = .false.
      this%is_loaded = .false.
      this%is_valid = .false.
      this%interpolation_method = 'bilinear'

   end subroutine extemifield_cleanup

   subroutine extemifield_validate(this, error_mgr, rc)
      implicit none
      class(ExtEmisFieldType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = cc_success

      if (.not. this%is_loaded) then
         call error_mgr%report_error(error_invalid_state, &
            'Field not loaded: ' // trim(this%field_name), rc)
         return
      endif

      if (.not. allocated(this%emission_data)) then
         call error_mgr%report_error(error_invalid_state, &
            'Emission data not allocated: ' // trim(this%field_name), rc)
         return
      endif

      ! Check for reasonable values (non-negative)
      if (any(this%emission_data < 0.0_fp)) then
         call error_mgr%report_error(error_invalid_state, &
            'Invalid emission values in field: ' // trim(this%field_name), rc)
         return
      endif

      this%is_valid = .true.

   end subroutine extemifield_validate

   function extemifield_get_emission_at_point(this, i, j, k, time_idx) result(emission_rate)
      implicit none
      class(ExtEmisFieldType), intent(in) :: this
      integer, intent(in) :: i, j
      integer, intent(in), optional :: k, time_idx
      real(fp) :: emission_rate

      integer :: k_use, time_use

      k_use = 1
      if (present(k)) k_use = k

      time_use = this%current_time_idx
      if (present(time_idx)) time_use = time_idx

      emission_rate = 0.0_fp

      if (this%is_loaded .and. allocated(this%emission_data)) then
         if (i >= 1 .and. i <= this%nx .and. &
            j >= 1 .and. j <= this%ny .and. &
            k_use >= 1 .and. k_use <= this%nz .and. &
            time_use >= 1 .and. time_use <= this%n_times) then
            emission_rate = this%emission_data(i, j, k_use, time_use)
         endif
      endif

   end function extemifield_get_emission_at_point

   subroutine extemifield_load_from_file(this, filename, error_mgr, rc)
      implicit none
      class(ExtEmisFieldType), intent(inout) :: this
      character(len=*), intent(in) :: filename
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = cc_success
      call error_mgr%push_context('extemifield_load_from_file', 'Loading emission field from: ' // trim(filename))

      ! Placeholder implementation - actual file I/O would be handled by the driver
      ! using NetCDF libraries or other format-specific readers

      this%is_loaded = .false.  ! Will be set to true by driver after successful load

      print *, 'WARNING: ExtEmisField file loading is handled by driver'
      call error_mgr%pop_context()

   end subroutine extemifield_load_from_file

   function extemifield_get_column_ptr(this, i, j) result(column_ptr)
      implicit none
      class(ExtEmisFieldType), intent(in), target :: this
      integer, intent(in) :: i, j
      real(fp), pointer :: column_ptr(:,:)

      column_ptr => null()

      if (this%is_loaded .and. allocated(this%emission_data)) then
         if (i >= 1 .and. i <= this%nx .and. &
            j >= 1 .and. j <= this%ny) then
            column_ptr => this%emission_data(i, j, :, :)
         endif
      endif

   end function extemifield_get_column_ptr

   !=========================================================================
   ! ExtEmisCategoryType procedures
   !=========================================================================

   subroutine extemicat_init(this, category_name, n_fields, description, rc)
      implicit none
      class(ExtEmisCategoryType), intent(inout) :: this
      character(len=*), intent(in) :: category_name
      integer, intent(in), optional :: n_fields
      character(len=*), intent(in), optional :: description
      integer, intent(out) :: rc

      integer :: alloc_stat

      rc = cc_success

      this%category_name = trim(category_name)

      if (present(description)) then
         this%description = trim(description)
      endif

      if (present(n_fields)) then
         this%n_fields = n_fields
      else
         this%n_fields = 0
      endif

      if (this%n_fields > 0) then
         if (allocated(this%fields)) deallocate(this%fields)
         allocate(this%fields(this%n_fields), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = cc_failure
            return
         endif
      endif

      this%global_scale = 1.0_fp
      this%is_active = .true.

   end subroutine extemicat_init

   subroutine extemicat_cleanup(this, rc)
      implicit none
      class(ExtEmisCategoryType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i

      rc = cc_success

      if (allocated(this%fields)) then
         do i = 1, this%n_fields
            call this%fields(i)%cleanup(rc)
            if (rc /= cc_success) return
         end do
         deallocate(this%fields)
      endif

      this%category_name = ''
      this%description = ''
      this%n_fields = 0
      this%irec = 0
      this%is_active = .true.
      this%gridded = .true.
      this%global_scale = 1.0_fp
      this%topfraction = -1.0_fp
      this%source_file = ''
      this%format = ''
      this%frequency = ''
      this%latname = ''
      this%lonname = ''
      this%stkdmname = ''
      this%stkhtname = ''
      this%stktkname = ''
      this%stkvename = ''
      this%plumerise = ''

   end subroutine extemicat_cleanup

   subroutine extemicat_validate(this, error_mgr, rc)
      implicit none
      class(ExtEmisCategoryType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      integer :: i

      rc = cc_success

      if (.not. allocated(this%fields)) then
         rc = cc_success  ! Empty category is valid
         return
      endif

      do i = 1, this%n_fields
         call this%fields(i)%validate(error_mgr, rc)
         if (rc /= cc_success) return
      end do

   end subroutine extemicat_validate

   subroutine extemicat_add_field(this, field, rc)
      implicit none
      class(ExtEmisCategoryType), intent(inout) :: this
      type(ExtEmisFieldType), intent(in) :: field
      integer, intent(out) :: rc

      type(ExtEmisFieldType), allocatable :: temp_fields(:)
      integer :: i, alloc_stat

      rc = cc_success

      ! Save existing fields to temporary array
      if (allocated(this%fields)) then
         allocate(temp_fields(this%n_fields), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = cc_failure
            return
         endif

         do i = 1, this%n_fields
            temp_fields(i) = this%fields(i)
         end do

         deallocate(this%fields)
      endif

      ! Allocate expanded array
      this%n_fields = this%n_fields + 1
      allocate(this%fields(this%n_fields), stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = cc_failure
         return
      endif

      ! Copy back existing fields
      if (allocated(temp_fields)) then
         do i = 1, this%n_fields - 1
            this%fields(i) = temp_fields(i)
         end do
         deallocate(temp_fields)
      endif

      ! Add new field
      this%fields(this%n_fields) = field

   end subroutine extemicat_add_field

   function extemicat_find_field(this, field_name) result(field_idx)
      implicit none
      class(ExtEmisCategoryType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      integer :: field_idx

      integer :: i

      field_idx = 0

      if (.not. allocated(this%fields)) return

      do i = 1, this%n_fields
         if (trim(this%fields(i)%field_name) == trim(field_name)) then
            field_idx = i
            return
         endif
      end do

   end function extemicat_find_field

   !=========================================================================
   ! ExtEmisDataType procedures
   !=========================================================================

   subroutine extemidata_init(this, n_categories, data_source, rc)
      implicit none
      class(ExtEmisDataType), intent(inout) :: this
      integer, intent(in), optional :: n_categories
      character(len=*), intent(in), optional :: data_source
      integer, intent(out) :: rc

      integer :: alloc_stat

      rc = cc_success

      if (present(data_source)) then
         this%data_source = trim(data_source)
      endif

      if (present(n_categories)) then
         this%n_categories = n_categories
      else
         this%n_categories = 0
      endif

      if (this%n_categories > 0) then
         if (allocated(this%categories)) deallocate(this%categories)
         allocate(this%categories(this%n_categories), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = cc_failure
            return
         endif
      endif

      this%is_active = .true.
      this%total_fields = 0
      this%global_scale = 1.0_fp

   end subroutine extemidata_init

   subroutine extemidata_cleanup(this, rc)
      implicit none
      class(ExtEmisDataType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i

      rc = cc_success

      if (allocated(this%categories)) then
         do i = 1, this%n_categories
            call this%categories(i)%cleanup(rc)
            if (rc /= cc_success) return
         end do
         deallocate(this%categories)
      endif

      this%n_categories = 0
      this%is_active = .true.
      this%total_fields = 0
      this%data_source = ''
      this%global_scale = 1.0_fp

   end subroutine extemidata_cleanup

   subroutine extemidata_validate(this, error_mgr, rc)
      implicit none
      class(ExtEmisDataType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      integer :: i

      rc = cc_success

      if (.not. allocated(this%categories)) then
         rc = cc_success  ! Empty container is valid
         return
      endif

      do i = 1, this%n_categories
         call this%categories(i)%validate(error_mgr, rc)
         if (rc /= cc_success) return
      end do

   end subroutine extemidata_validate

   subroutine extemidata_add_category(this, category, rc)
      implicit none
      class(ExtEmisDataType), intent(inout) :: this
      type(ExtEmisCategoryType), intent(in) :: category
      integer, intent(out) :: rc

      type(ExtEmisCategoryType), allocatable :: temp_categories(:)
      integer :: i, alloc_stat

      rc = cc_success

      ! Save existing categories to temporary array
      if (allocated(this%categories)) then
         allocate(temp_categories(this%n_categories), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = cc_failure
            return
         endif

         do i = 1, this%n_categories
            temp_categories(i) = this%categories(i)
         end do

         deallocate(this%categories)
      endif

      ! Allocate expanded array
      this%n_categories = this%n_categories + 1
      allocate(this%categories(this%n_categories), stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = cc_failure
         return
      endif

      ! Copy back existing categories
      if (allocated(temp_categories)) then
         do i = 1, this%n_categories - 1
            this%categories(i) = temp_categories(i)
         end do
         deallocate(temp_categories)
      endif

      ! Add new category
      this%categories(this%n_categories) = category
      this%total_fields = this%total_fields + category%n_fields

   end subroutine extemidata_add_category

   subroutine extemidata_load_files(this, file_list, error_mgr, rc)
      implicit none
      class(ExtEmisDataType), intent(inout) :: this
      character(len=*), intent(in) :: file_list(:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = cc_success
      call error_mgr%push_context('extemidata_load_files', 'Loading emission files')

      ! Placeholder - actual file loading is handled by the driver
      print *, 'WARNING: ExtEmisData file loading is handled by driver'

      call error_mgr%pop_context()

   end subroutine extemidata_load_files

   function extemidata_find_field(this, field_name) result(field_ptr)
      implicit none
      class(ExtEmisDataType), intent(in), target :: this
      character(len=*), intent(in) :: field_name
      type(ExtEmisFieldType), pointer :: field_ptr

      integer :: i, field_idx

      field_ptr => null()

      if (.not. allocated(this%categories)) return

      do i = 1, this%n_categories
         field_idx = this%categories(i)%find_field(field_name)
         if (field_idx > 0) then
            field_ptr => this%categories(i)%fields(field_idx)
            return
         endif
      end do

   end function extemidata_find_field

   function extemidata_get_emission_rate(this, field_name, i, j, k) result(emission_rate)
      implicit none
      class(ExtEmisDataType), intent(in), target :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: i, j
      integer, intent(in), optional :: k
      real(fp) :: emission_rate

      type(ExtEmisFieldType), pointer :: field_ptr

      emission_rate = 0.0_fp

      field_ptr => this%find_emission_field(field_name)
      if (associated(field_ptr)) then
         emission_rate = field_ptr%get_emission_at_point(i, j, k)
      endif

   end function extemidata_get_emission_rate

   subroutine extemidata_update_time(this, time_idx, rc)
      implicit none
      class(ExtEmisDataType), intent(inout) :: this
      integer, intent(in) :: time_idx
      integer, intent(out) :: rc

      integer :: i, j

      rc = cc_success

      if (.not. allocated(this%categories)) return

      do i = 1, this%n_categories
         if (.not. allocated(this%categories(i)%fields)) cycle
         do j = 1, this%categories(i)%n_fields
            this%categories(i)%fields(j)%current_time_idx = time_idx
         end do
      end do

   end subroutine extemidata_update_time

   function extemidata_get_memory_usage(this) result(memory_bytes)
      implicit none
      class(ExtEmisDataType), intent(in) :: this
      integer(kind=8) :: memory_bytes

      integer :: i, j
      real(fp) :: field_size_mb

      memory_bytes = 0

      do i = 1, this%n_categories
         do j = 1, this%categories(i)%n_fields
            if (allocated(this%categories(i)%fields(j)%emission_data)) then
               field_size_mb = real(this%categories(i)%fields(j)%nx * &
                  this%categories(i)%fields(j)%ny * &
                  this%categories(i)%fields(j)%nz * &
                  this%categories(i)%fields(j)%n_times, fp) * 8.0_fp  ! 8 bytes per real(fp)
               memory_bytes = memory_bytes + int(field_size_mb, kind=8)
            endif
         end do
      end do

   end function extemidata_get_memory_usage

   function extemidata_get_column_ptr(this, category_name, field_name, i, j) result(column_ptr)
      implicit none
      class(ExtEmisDataType), intent(in), target :: this
      character(len=*), intent(in) :: category_name, field_name
      integer, intent(in) :: i, j
      real(fp), pointer :: column_ptr(:,:)

      integer :: cat_idx, field_idx
      type(ExtEmisFieldType), pointer :: field_ptr

      column_ptr => null()

      ! If category name provided, search only that category
      if (len_trim(category_name) > 0) then
         do cat_idx = 1, this%n_categories
            if (trim(this%categories(cat_idx)%category_name) == trim(category_name)) then
               field_idx = this%categories(cat_idx)%find_field(field_name)
               if (field_idx > 0) then
                  column_ptr => this%categories(cat_idx)%fields(field_idx)%get_column_ptr(i, j)
               endif
               return
            endif
         end do
      else
         ! Search all categories
         field_ptr => this%find_emission_field(field_name)
         if (associated(field_ptr)) then
            column_ptr => field_ptr%get_column_ptr(i, j)
         endif
      endif

   end function extemidata_get_column_ptr

END MODULE extemisdata_mod
```


