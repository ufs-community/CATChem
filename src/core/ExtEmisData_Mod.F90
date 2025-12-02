!> \file ExtEmisData_Mod.F90
!! \brief Module for external emission data storage
!!
!! This module contains the ExtEmisDataType derived type for managing
!! external emission data loaded from files by the driver. This is separate
!! from the EmisState which is used for accumulating computed emissions.
!!
!! \ingroup core_modules
!! \author CATChem Development Team
!! \date 2023
!! \version 1.0
!!
MODULE ExtEmisData_Mod

   !=========================================================================
   ! Module uses
   !=========================================================================
   USE Error_Mod
   ! USE logging_mod, only: log_message, LOG_ERROR, LOG_WARNING, LOG_INFO, LOG_DEBUG
   USE Precision_Mod
   USE species_mod, only: SpeciesType
   USE IEEE_ARITHMETIC

   IMPLICIT NONE
   PRIVATE

   !=========================================================================
   ! Public interfaces
   !=========================================================================
   !PUBLIC :: ExtEmisDataType
   !PUBLIC :: ExtEmisFieldType
   !PUBLIC :: ExtEmisCategoryType

   !> \brief Derived type for individual external emission field
   !!
   !! Contains data for a single emission field from external files,
   !! including spatial and temporal information.
   !!
   !! \ingroup core_modules
   !!!>
   TYPE, PUBLIC :: ExtEmisFieldType
      CHARACTER(LEN=64)             :: field_name = ''     !< Field name (e.g., "EMIS_NO", "EMIS_SO2")
      CHARACTER(LEN=128)            :: long_name = ''      !< Descriptive name
      CHARACTER(LEN=32)             :: units = ''          !< Units (typically kg/m2/s)
      !CHARACTER(LEN=256)            :: source_file = ''    !< Source file path
      INTEGER                       :: nx = 0              !< Number of longitude points
      INTEGER                       :: ny = 0              !< Number of latitude points
      INTEGER                       :: nz = 1              !< Number of vertical levels (usually 1 for surface)
      REAL(fp)                      :: factors = 1.0_fp    !< Scaling factors for non-chemical variables
      REAL(fp), ALLOCATABLE         :: lat(:)              !< Latitude coordinates [degrees]
      REAL(fp), ALLOCATABLE         :: lon(:)              !< Longitude coordinates [degrees]
      REAL(fp), ALLOCATABLE         :: stkdm(:)            !< Stack diameter [m]
      REAL(fp), ALLOCATABLE         :: stkht(:)            !< Stack height [m]
      REAL(fp), ALLOCATABLE         :: stktk(:)            !< Stack temperature [K]
      REAL(fp), ALLOCATABLE         :: stkve(:)            !< Stack velocity [m/s]
      INTEGER,  ALLOCATABLE         :: ip(:)               !< i-indices for point sources
      INTEGER,  ALLOCATABLE         :: jp(:)               !< j-indices for point sources
      INTEGER,  ALLOCATABLE         :: ijmap(:)            !< Number of point sources within the model domain
      INTEGER                       :: n_times = 0         !< Number of time steps in file
      INTEGER                       :: current_time_idx = 1 !< Current time index
      LOGICAL                       :: time_interpolate = .true. !< Enable time interpolation
      LOGICAL                       :: diagnostic = .false. !< Enable diagnostic output of this field
      REAL(fp), ALLOCATABLE         :: emission_data(:,:,:,:) !< Emission flux [kg/m2/s] (nx,ny,nz,n_times)
      ! REAL(fp), ALLOCATABLE         :: longitude(:,:)         !< Longitude coordinates [degrees]
      ! REAL(fp), ALLOCATABLE         :: latitude(:,:)          !< Latitude coordinates [degrees]
      ! REAL(fp), ALLOCATABLE         :: vertical(:,:)          !< Vertical coordinates (if applicable)
      ! REAL(fp), ALLOCATABLE         :: time_coords(:,:)       !< Time coordinates
      LOGICAL                       :: is_loaded = .false. !< Data loading status
      LOGICAL                       :: is_valid = .false.  !< Data validation status
      CHARACTER(LEN=32)             :: interpolation_method = 'bilinear' !< Spatial interpolation method
   CONTAINS
      !> \brief Initialize emission field with metadata
      !! \copydoc extemifield_init
      PROCEDURE :: init => extemifield_init
      !> \brief Clean up resources and deallocate arrays
      !! \copydoc extemifield_cleanup
      PROCEDURE :: cleanup => extemifield_cleanup
      !> \brief Validate emission field data
      !! \copydoc extemifield_validate
      PROCEDURE :: validate => extemifield_validate
      !> \brief Get emission rate at a specific grid point
      !! \copydoc extemifield_get_emission_at_point
      PROCEDURE :: get_emission_at_point => extemifield_get_emission_at_point
      !> \brief Load emission field from file
      !! \copydoc extemifield_load_from_file
      PROCEDURE :: load_from_file => extemifield_load_from_file
      !> \brief Get a column slice of emission data
      !! \copydoc extemifield_get_column_ptr
      PROCEDURE :: get_column_ptr => extemifield_get_column_ptr
   END TYPE ExtEmisFieldType

   !> \brief Derived type for external emission categories
   !!
   !! Groups related emission fields by source category
   !! (e.g., anthropogenic, biogenic, fires, etc.)
   !!
   !! \ingroup core_modules
   !!!>
   TYPE, PUBLIC :: ExtEmisCategoryType
      CHARACTER(LEN=64)                         :: category_name = ''  !< Category name
      CHARACTER(LEN=256)                        :: description = ''    !< Category description
      INTEGER                                   :: n_fields = 0        !< Number of emission fields
      INTEGER                                   :: irec = 0            !< time slice index
      TYPE(ExtEmisFieldType), ALLOCATABLE       :: fields(:)           !< Emission fields array
      LOGICAL                                   :: is_active = .true.  !< Category enabled/disabled
      LOGICAL                                   :: gridded = .true.   !< Is this a gridded emission category
      REAL(fp)                                  :: global_scale = 1.0_fp !< Global scaling factor
      REAL(fp)                                  :: topfraction = -1.0_fp !< Top fraction for plumerise
      CHARACTER(LEN=128)                        :: source_file = ''    !< Source file path and name
      CHARACTER(LEN=128)                        :: format = ''         !< Format of file (only netcdf for now)
      CHARACTER(LEN=128)                        :: frequency = ''      !< Frequency of file (e.g., hourly, daily, weekly, monthly,static)
      CHARACTER(LEN=128)                        :: latname = ''       !< Latitude variable name in the file
      CHARACTER(LEN=128)                        :: lonname = ''       !< Longitude variable name in the file
      CHARACTER(LEN=128)                        :: stkdmname = ''      !< Stack dimension name in the file
      CHARACTER(LEN=128)                        :: stkhtname = ''      !< Stack height variable name in the file
      CHARACTER(LEN=128)                        :: stktkname = ''      !< Stack temperature variable name in the file
      CHARACTER(LEN=128)                        :: stkvename = ''      !< Stack velocity variable name in the file
      CHARACTER(LEN=128)                        :: plumerise = ''      !< plumerise scheme


   CONTAINS
      !> \brief Initialize emission category with metadata
      !! \copydoc extemicat_init
      PROCEDURE :: init => extemicat_init
      !> \brief Clean up resources and deallocate arrays
      !! \copydoc extemicat_cleanup
      PROCEDURE :: cleanup => extemicat_cleanup
      !> \brief Validate emission category data
      !! \copydoc extemicat_validate
      PROCEDURE :: validate => extemicat_validate
      !> \brief Add a new field to the category
      !! \copydoc extemicat_add_field
      PROCEDURE :: add_field => extemicat_add_field
      !> \brief Find a field by name
      !! \copydoc extemicat_find_field
      PROCEDURE :: find_field => extemicat_find_field
   END TYPE ExtEmisCategoryType

   !> \brief Main container for external emission data
   !!
   !! Manages emission data loaded from external files by the driver
   !! organized into categories and fields.
   !!
   !! \ingroup core_modules
   !!!>
   TYPE, PUBLIC :: ExtEmisDataType
      INTEGER                                   :: n_categories = 0    !< Number of emission categories
      TYPE(ExtEmisCategoryType), ALLOCATABLE    :: categories(:)       !< Emission categories array
      LOGICAL                                   :: is_active = .true.  !< All emissions enabled/disabled
      LOGICAL                                   :: diagnostic = .true.  !< Enable diagnostic output for external emissions?
      INTEGER                                   :: total_fields = 0    !< Total number of emission fields
      CHARACTER(LEN=128)                        :: data_source = ''    !< Data source information
      REAL(fp)                                  :: global_scale = 1.0_fp !< Global scaling factor for all emissions
   CONTAINS
      !> \brief Initialize emission data container
      !! \copydoc extemidata_init
      PROCEDURE :: init => extemidata_init
      !> \brief Clean up resources and deallocate arrays
      !! \copydoc extemidata_cleanup
      PROCEDURE :: cleanup => extemidata_cleanup
      !> \brief Validate all emission data
      !! \copydoc extemidata_validate
      PROCEDURE :: validate => extemidata_validate
      !> \brief Add a new category to the container
      !! \copydoc extemidata_add_category
      PROCEDURE :: add_category => extemidata_add_category
      !> \brief Load emission files
      !! \copydoc extemidata_load_files
      PROCEDURE :: load_emission_files => extemidata_load_files
      !> \brief Find emission field across all categories
      !! \copydoc extemidata_find_field
      PROCEDURE :: find_emission_field => extemidata_find_field
      !> \brief Get emission rate for a specific field and location
      !! \copydoc extemidata_get_emission_rate
      PROCEDURE :: get_emission_rate => extemidata_get_emission_rate
      !> \brief Update current time for interpolation
      !! \copydoc extemidata_update_time
      PROCEDURE :: update_time => extemidata_update_time
      !> \brief Get memory usage estimate
      !! \copydoc extemidata_get_memory_usage
      PROCEDURE :: get_memory_usage => extemidata_get_memory_usage
      !> \brief Get a pointer to the (k, t) slab for a given (i, j) column in a named field
      !!
      !! Returns a pointer to the emission_data(i, j, :, :) for the specified field and column.
      !!
      !! \param[in] this The ExtEmisDataType object
      !! \param[in] category_name Name of the emission category
      !! \param[in] field_name Name of the emission field
      !! \param[in] i Longitude index
      !! \param[in] j Latitude index
      !! \return Pointer to emission_data(i, j, :, :) (vertical x time), or null if not found
      PROCEDURE :: get_column_ptr => extemidata_get_column_ptr
   END TYPE ExtEmisDataType

CONTAINS

   !=========================================================================
   ! ExtEmisFieldType procedures
   !=========================================================================

   !> \brief Initialize an emission field with metadata
   !!
   !! Sets up field metadata including dimensions and coordinates
   !!
   !! \param[inout] this The ExtEmisFieldType object
   !! \param[in] field_name Field name
   !! \param[in] nx Number of longitude points
   !! \param[in] ny Number of latitude points
   !! \param[in] nz Number of vertical levels (optional, default=1)
   !! \param[in] n_times Number of time steps (optional, default=1)
   !! \param[in] units Units string (optional, default='kg/m2/s')
   !! \param[out] rc Return code
   subroutine extemifield_init(this, field_name, nx, ny, nz, n_times, units, rc)
      implicit none
      class(ExtEmisFieldType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: nx, ny
      integer, intent(in), optional :: nz, n_times
      character(len=*), intent(in), optional :: units
      integer, intent(out) :: rc

      rc = CC_SUCCESS

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

      this%current_time_idx = 1
      this%factors = 1.0_fp
      this%time_interpolate = .true.
      this%diagnostic = .false.
      this%is_loaded = .false.
      this%is_valid = .false.

   end subroutine extemifield_init

   !> \brief Clean up resources and deallocate arrays
   !!
   !! Deallocates all arrays and resets fields to default values
   !!
   !! \param[inout] this The ExtEmisFieldType object
   !! \param[out] rc Return code
   subroutine extemifield_cleanup(this, rc)
      implicit none
      class(ExtEmisFieldType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

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

   !> \brief Validate emission field data
   !!
   !! Checks that the emission field data is properly loaded and valid
   !!
   !! \param[inout] this The ExtEmisFieldType object
   !! \param[inout] error_mgr Error manager for context and reporting
   !! \param[out] rc Return code
   subroutine extemifield_validate(this, error_mgr, rc)
      implicit none
      class(ExtEmisFieldType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      if (.not. this%is_loaded) then
         call error_mgr%report_error(ERROR_INVALID_STATE, &
            'Field not loaded: ' // trim(this%field_name), rc)
         return
      endif

      if (.not. allocated(this%emission_data)) then
         call error_mgr%report_error(ERROR_INVALID_STATE, &
            'Emission data not allocated: ' // trim(this%field_name), rc)
         return
      endif

      ! Check for reasonable values (non-negative)
      if (any(this%emission_data < 0.0_fp)) then
         call error_mgr%report_error(ERROR_INVALID_STATE, &
            'Invalid emission values in field: ' // trim(this%field_name), rc)
         return
      endif

      this%is_valid = .true.

   end subroutine extemifield_validate

   !> \brief Get emission rate at a specific grid point
   !!
   !! Returns the emission rate for the specified grid point
   !!
   !! \param[in] this The ExtEmisFieldType object
   !! \param[in] i Longitude index
   !! \param[in] j Latitude index
   !! \param[in] k Vertical index (optional, default=1)
   !! \param[in] time_idx Time index (optional, uses current_time_idx if not specified)
   !! \return Emission rate [kg/m2/s] at the specified point
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

   !> \brief Load emission data from file (placeholder)
   !!
   !! Placeholder implementation for loading field data from a file.
   !! Actual file I/O would be handled by the driver.
   !!
   !! \param[inout] this The ExtEmisFieldType object
   !! \param[in] filename File path to load data from
   !! \param[inout] error_mgr Error manager for context and reporting
   !! \param[out] rc Return code
   subroutine extemifield_load_from_file(this, filename, error_mgr, rc)
      implicit none
      class(ExtEmisFieldType), intent(inout) :: this
      character(len=*), intent(in) :: filename
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS
      call error_mgr%push_context('extemifield_load_from_file', 'Loading emission field from: ' // trim(filename))

      ! Placeholder implementation - actual file I/O would be handled by the driver
      ! using NetCDF libraries or other format-specific readers

      this%is_loaded = .false.  ! Will be set to true by driver after successful load

      print *, 'WARNING: ExtEmisField file loading is handled by driver'
      call error_mgr%pop_context()

   end subroutine extemifield_load_from_file

   !> \brief Get a pointer to the (k, t) slab for a given (i, j) column in emission_data
   !!
   !! Returns a pointer to emission_data(i, j, :, :) for the specified (i, j) column.
   !!
   !! \param[in] this The ExtEmisFieldType object
   !! \param[in] i Longitude index
   !! \param[in] j Latitude index
   !! \return Pointer to emission_data(i, j, :, :) (vertical x time), or null if not found
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

   !> \brief Initialize an emission category
   !!
   !! Sets up category metadata and allocates field arrays
   !!
   !! \param[inout] this The ExtEmisCategoryType object
   !! \param[in] category_name Category name
   !! \param[in] n_fields Number of emission fields (optional, default=0)
   !! \param[in] description Category description (optional)
   !! \param[out] rc Return code
   subroutine extemicat_init(this, category_name, n_fields, description, rc)
      implicit none
      class(ExtEmisCategoryType), intent(inout) :: this
      character(len=*), intent(in) :: category_name
      integer, intent(in), optional :: n_fields
      character(len=*), intent(in), optional :: description
      integer, intent(out) :: rc

      integer :: alloc_stat

      rc = CC_SUCCESS

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
            rc = CC_FAILURE
            return
         endif
      endif

      this%global_scale = 1.0_fp
      this%is_active = .true.

   end subroutine extemicat_init

   !> \brief Clean up resources and deallocate arrays
   !!
   !! Deallocates all fields and resets category to default values
   !!
   !! \param[inout] this The ExtEmisCategoryType object
   !! \param[out] rc Return code
   subroutine extemicat_cleanup(this, rc)
      implicit none
      class(ExtEmisCategoryType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      if (allocated(this%fields)) then
         do i = 1, this%n_fields
            call this%fields(i)%cleanup(rc)
            if (rc /= CC_SUCCESS) return
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

   !> \brief Validate all fields in the category
   !!
   !! Checks that all emission fields are properly loaded and valid
   !!
   !! \param[inout] this The ExtEmisCategoryType object
   !! \param[inout] error_mgr Error manager for context and reporting
   !! \param[out] rc Return code
   subroutine extemicat_validate(this, error_mgr, rc)
      implicit none
      class(ExtEmisCategoryType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      if (.not. allocated(this%fields)) then
         rc = CC_SUCCESS  ! Empty category is valid
         return
      endif

      do i = 1, this%n_fields
         call this%fields(i)%validate(error_mgr, rc)
         if (rc /= CC_SUCCESS) return
      end do

   end subroutine extemicat_validate

   !> \brief Add a new field to the category
   !!
   !! Expands the fields array and adds a new emission field
   !!
   !! \param[inout] this The ExtEmisCategoryType object
   !! \param[in] field ExtEmisFieldType to add
   !! \param[out] rc Return code
   subroutine extemicat_add_field(this, field, rc)
      implicit none
      class(ExtEmisCategoryType), intent(inout) :: this
      type(ExtEmisFieldType), intent(in) :: field
      integer, intent(out) :: rc

      type(ExtEmisFieldType), allocatable :: temp_fields(:)
      integer :: i, alloc_stat

      rc = CC_SUCCESS

      ! Save existing fields to temporary array
      if (allocated(this%fields)) then
         allocate(temp_fields(this%n_fields), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = CC_FAILURE
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
         rc = CC_FAILURE
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

   !> \brief Find a field by name
   !!
   !! Searches for a field with the given name in this category.
   !! Returns the index if found, or 0 if not found.
   !!
   !! \param[in] this The ExtEmisCategoryType object
   !! \param[in] field_name Name of the emission field to find
   !! \return Field index (1-based) or 0 if not found
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

   !> \brief Initialize emission data container
   !!
   !! Sets up container metadata and allocates category arrays
   !!
   !! \param[inout] this The ExtEmisDataType object
   !! \param[in] n_categories Number of emission categories (optional, default=0)
   !! \param[in] data_source Source information (optional)
   !! \param[out] rc Return code
   subroutine extemidata_init(this, n_categories, data_source, rc)
      implicit none
      class(ExtEmisDataType), intent(inout) :: this
      integer, intent(in), optional :: n_categories
      character(len=*), intent(in), optional :: data_source
      integer, intent(out) :: rc

      integer :: alloc_stat

      rc = CC_SUCCESS

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
            rc = CC_FAILURE
            return
         endif
      endif

      this%is_active = .true.
      this%total_fields = 0
      this%global_scale = 1.0_fp

   end subroutine extemidata_init

   !> \brief Clean up resources and deallocate arrays
   !!
   !! Deallocates all categories and fields and resets container to default values
   !!
   !! \param[inout] this The ExtEmisDataType object
   !! \param[out] rc Return code
   subroutine extemidata_cleanup(this, rc)
      implicit none
      class(ExtEmisDataType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      if (allocated(this%categories)) then
         do i = 1, this%n_categories
            call this%categories(i)%cleanup(rc)
            if (rc /= CC_SUCCESS) return
         end do
         deallocate(this%categories)
      endif

      this%n_categories = 0
      this%is_active = .true.
      this%total_fields = 0
      this%data_source = ''
      this%global_scale = 1.0_fp

   end subroutine extemidata_cleanup

   !> \brief Validate all emission data
   !!
   !! Checks that all categories and fields are properly loaded and valid
   !!
   !! \param[inout] this The ExtEmisDataType object
   !! \param[inout] error_mgr Error manager for context and reporting
   !! \param[out] rc Return code
   subroutine extemidata_validate(this, error_mgr, rc)
      implicit none
      class(ExtEmisDataType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      if (.not. allocated(this%categories)) then
         rc = CC_SUCCESS  ! Empty container is valid
         return
      endif

      do i = 1, this%n_categories
         call this%categories(i)%validate(error_mgr, rc)
         if (rc /= CC_SUCCESS) return
      end do

   end subroutine extemidata_validate

   !> \brief Add a new category to the container
   !!
   !! Expands the categories array and adds a new emission category
   !!
   !! \param[inout] this The ExtEmisDataType object
   !! \param[in] category ExtEmisCategoryType to add
   !! \param[out] rc Return code
   subroutine extemidata_add_category(this, category, rc)
      implicit none
      class(ExtEmisDataType), intent(inout) :: this
      type(ExtEmisCategoryType), intent(in) :: category
      integer, intent(out) :: rc

      type(ExtEmisCategoryType), allocatable :: temp_categories(:)
      integer :: i, alloc_stat

      rc = CC_SUCCESS

      ! Save existing categories to temporary array
      if (allocated(this%categories)) then
         allocate(temp_categories(this%n_categories), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = CC_FAILURE
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
         rc = CC_FAILURE
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

   !> \brief Load emission files (placeholder)
   !!
   !! Placeholder implementation for loading all emission data from files.
   !! Actual file I/O would be handled by the driver.
   !!
   !! \param[inout] this The ExtEmisDataType object
   !! \param[in] file_list List of files to load
   !! \param[inout] error_mgr Error manager for context and reporting
   !! \param[out] rc Return code
   subroutine extemidata_load_files(this, file_list, error_mgr, rc)
      implicit none
      class(ExtEmisDataType), intent(inout) :: this
      character(len=*), intent(in) :: file_list(:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS
      call error_mgr%push_context('extemidata_load_files', 'Loading emission files')

      ! Placeholder - actual file loading is handled by the driver
      print *, 'WARNING: ExtEmisData file loading is handled by driver'

      call error_mgr%pop_context()

   end subroutine extemidata_load_files

   !> \brief Find emission field across all categories
   !!
   !! Returns a pointer to the emission field with the given name, searching all categories.
   !!
   !! \param[in] this The ExtEmisDataType object
   !! \param[in] field_name Name of the emission field
   !! \return Pointer to the emission field (or null)
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

   !> \brief Get emission rate for a specific field and location
   !!
   !! Returns the emission rate for a given field at the specified grid point.
   !!
   !! \param[in] this The ExtEmisDataType object
   !! \param[in] field_name Name of the emission field
   !! \param[in] i Longitude index
   !! \param[in] j Latitude index
   !! \param[in] k Vertical index (optional)
   !! \return Emission rate [kg/m2/s] at the specified location
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

   !> \brief Update current time for interpolation
   !!
   !! Updates the current time index for all fields for time interpolation.
   !!
   !! \param[inout] this The ExtEmisDataType object
   !! \param[in] time_idx New time index
   !! \param[out] rc Return code
   subroutine extemidata_update_time(this, time_idx, rc)
      implicit none
      class(ExtEmisDataType), intent(inout) :: this
      integer, intent(in) :: time_idx
      integer, intent(out) :: rc

      integer :: i, j

      rc = CC_SUCCESS

      if (.not. allocated(this%categories)) return

      do i = 1, this%n_categories
         if (.not. allocated(this%categories(i)%fields)) cycle
         do j = 1, this%categories(i)%n_fields
            this%categories(i)%fields(j)%current_time_idx = time_idx
         end do
      end do

   end subroutine extemidata_update_time

   !> \brief Get memory usage estimate
   !!
   !! Returns the estimated memory usage in bytes for all emission data.
   !!
   !! \param[in] this The ExtEmisDataType object
   !! \return Memory usage in bytes
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

   !> \brief Get a pointer to the (k, t) slab for a given (i, j) column in a named field
   !!
   !! Returns a pointer to the emission_data(i, j, :, :) for the specified field and column.
   !!
   !! \param[in] this The ExtEmisDataType object
   !! \param[in] category_name Name of the emission category
   !! \param[in] field_name Name of the emission field
   !! \param[in] i Longitude index
   !! \param[in] j Latitude index
   !! \return Pointer to emission_data(i, j, :, :) (vertical x time), or null if not found
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

END MODULE ExtEmisData_Mod
