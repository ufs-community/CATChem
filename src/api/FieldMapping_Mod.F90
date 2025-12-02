!> \file FieldMapping_Mod.F90
!! \brief Field mapping system for CATChem high-level API
!! \ingroup catchem_api
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module provides a flexible field mapping system that allows host models
!! to map their field names to CATChem field names. This enables seamless
!! integration without requiring host models to use CATChem's internal naming
!! conventions.
!!
!! Example usage:
!! \code{.f90}
!!   type(FieldMappingType) :: mapper
!!   call mapper%init()
!!
!!   ! Map host model fields to CATChem fields
!!   call mapper%add_mapping('host_temp', 'temperature', 'meteo', rc)
!!   call mapper%add_mapping('host_pres', 'pressure', 'meteo', rc)
!!   call mapper%add_mapping('host_o3', 'O3', 'chemistry', rc)
!!
!!   ! Use mapping to set data
!!   call catchem%set_field_by_mapping(mapper, 'host_temp', temp_data, rc)
!! \endcode
!!
module FieldMapping_Mod
   use precision_mod
   use error_mod

   implicit none
   private

   public :: FieldMappingType
   public :: FieldMappingEntryType
   public :: MAPPING_SUCCESS, MAPPING_FAILURE
   public :: MAPPING_FIELD_NOT_FOUND, MAPPING_INVALID_CATEGORY

   ! Return codes
   integer, parameter :: MAPPING_SUCCESS = 0
   integer, parameter :: MAPPING_FAILURE = -1
   integer, parameter :: MAPPING_FIELD_NOT_FOUND = -2
   integer, parameter :: MAPPING_INVALID_CATEGORY = -3

   ! Field categories
   character(len=*), parameter :: CATEGORY_METEO = 'meteo'
   character(len=*), parameter :: CATEGORY_CHEMISTRY = 'chemistry'
   character(len=*), parameter :: CATEGORY_EMISSIONS = 'emissions'
   character(len=*), parameter :: CATEGORY_DIAGNOSTICS = 'diagnostics'

   ! Field dimensionality types
   integer, parameter :: FIELD_1D = 1
   integer, parameter :: FIELD_2D = 2
   integer, parameter :: FIELD_3D = 3

   ! Mapping direction types
   integer, parameter :: MAPPING_HOST_TO_CATCHEM = 1
   integer, parameter :: MAPPING_CATCHEM_TO_HOST = 2
   integer, parameter :: MAPPING_BIDIRECTIONAL = 3

   !> Individual field mapping entry
   type :: FieldMappingEntryType
      character(len=64) :: host_name = ''      !< Host model field name
      character(len=64) :: catchem_name = ''   !< CATChem field name
      character(len=32) :: category = ''       !< Field category (meteo, chemistry, etc.)
      character(len=32) :: units = ''          !< Expected units
      character(len=128) :: description = ''   !< Field description
      logical :: is_required = .false.         !< Whether field is required
      logical :: is_mapped = .false.           !< Whether mapping is active
      integer :: mapping_direction = MAPPING_HOST_TO_CATCHEM  !< Mapping direction

      ! Field dimensionality and shape information
      integer :: host_dims = 0                 !< Host field dimensionality (1D, 2D, 3D)
      integer :: catchem_dims = 0              !< CATChem field dimensionality
      integer :: host_shape(3) = 0             !< Host field shape [dim1, dim2, dim3]
      integer :: catchem_shape(3) = 0          !< CATChem field shape [nx, ny, nz]
      logical :: needs_reshape = .false.       !< Whether reshaping is needed

      ! Dimension mapping for flexible reshaping
      ! host_dim_map(i) indicates which CATChem dimension corresponds to host dimension i
      integer :: host_dim_map(3) = [1, 2, 3]   !< Default: 1D->nx, 2D->nx,ny, 3D->nx,ny,nz

      ! For specialized mappings (e.g., nHoriz -> nx*ny, nLev -> nz)
      logical :: is_column_data = .false.      !< True if host field is column-based (nHoriz, nLev)
      logical :: expand_horizontal = .false.   !< True if nHoriz should expand to nx*ny

      ! Data validation parameters
      real(fp) :: min_value = -huge(1.0_fp)    !< Minimum valid value
      real(fp) :: max_value = huge(1.0_fp)     !< Maximum valid value
      logical :: check_bounds = .false.        !< Whether to check bounds
   end type FieldMappingEntryType

   !> Main field mapping manager
   type :: FieldMappingType
      private

      integer :: num_mappings = 0
      integer, parameter :: max_mappings = 1000
      type(FieldMappingEntryType) :: mappings(max_mappings)
      logical :: is_initialized = .false.

   contains
      ! Core mapping operations
      procedure :: init => mapping_init
      procedure :: cleanup => mapping_cleanup
      procedure :: add_mapping => mapping_add
      procedure :: remove_mapping => mapping_remove
      procedure :: clear_all => mapping_clear_all

      ! Query operations
      procedure :: find_mapping => mapping_find
      procedure :: get_catchem_name => mapping_get_catchem_name
      procedure :: get_host_name => mapping_get_host_name
      procedure :: get_category => mapping_get_category
      procedure :: is_mapped => mapping_is_mapped
      procedure :: is_required => mapping_is_required

      ! Bulk operations
      procedure :: get_all_mappings => mapping_get_all
      procedure :: get_mappings_by_category => mapping_get_by_category
      procedure :: load_from_file => mapping_load_from_file
      procedure :: save_to_file => mapping_save_to_file

      ! Validation
      procedure :: validate_mapping => mapping_validate
      procedure :: validate_all => mapping_validate_all
      procedure :: check_required_fields => mapping_check_required

      ! Utility
      procedure :: print_summary => mapping_print_summary
      procedure :: get_stats => mapping_get_stats

      ! Enhanced field mapping operations
      procedure :: set_field_dimensions => mapping_set_dimensions
      procedure :: setup_reshaping => mapping_setup_reshaping
      procedure :: reshape_data_1d_to_3d => mapping_reshape_1d_to_3d
      procedure :: reshape_data_2d_to_3d => mapping_reshape_2d_to_3d
      procedure :: reshape_data_3d_to_1d => mapping_reshape_3d_to_1d
      procedure :: reshape_data_3d_to_2d => mapping_reshape_3d_to_2d

      ! Data transformation utilities
      procedure :: validate_shapes => mapping_validate_shapes
      procedure :: compute_reshape_params => mapping_compute_reshape_params
   end type FieldMappingType

contains

   !> Initialize the field mapping system
   subroutine mapping_init(this)
      class(FieldMappingType), intent(inout) :: this

      this%num_mappings = 0
      this%is_initialized = .true.

      ! Initialize all mappings
      this%mappings(:)%host_name = ''
      this%mappings(:)%catchem_name = ''
      this%mappings(:)%category = ''
      this%mappings(:)%is_mapped = .false.
      this%mappings(:)%is_required = .false.
   end subroutine mapping_init

   !> Clean up the field mapping system
   subroutine mapping_cleanup(this)
      class(FieldMappingType), intent(inout) :: this

      call this%clear_all()
      this%is_initialized = .false.
   end subroutine mapping_cleanup

   !> Add a new field mapping
   subroutine mapping_add(this, host_name, catchem_name, category, rc, &
      units, description, is_required, min_value, max_value)
      class(FieldMappingType), intent(inout) :: this
      character(len=*), intent(in) :: host_name
      character(len=*), intent(in) :: catchem_name
      character(len=*), intent(in) :: category
      integer, intent(out) :: rc
      character(len=*), intent(in), optional :: units
      character(len=*), intent(in), optional :: description
      logical, intent(in), optional :: is_required
      real(fp), intent(in), optional :: min_value, max_value

      integer :: idx

      rc = MAPPING_SUCCESS

      ! Check if initialized
      if (.not. this%is_initialized) then
         rc = MAPPING_FAILURE
         return
      endif

      ! Check if we have space
      if (this%num_mappings >= this%max_mappings) then
         rc = MAPPING_FAILURE
         return
      endif

      ! Validate category
      if (.not. is_valid_category(category)) then
         rc = MAPPING_INVALID_CATEGORY
         return
      endif

      ! Check if mapping already exists
      idx = this%find_mapping(host_name)
      if (idx > 0) then
         ! Update existing mapping
         this%mappings(idx)%catchem_name = trim(catchem_name)
         this%mappings(idx)%category = trim(category)
      else
         ! Add new mapping
         this%num_mappings = this%num_mappings + 1
         idx = this%num_mappings
         this%mappings(idx)%host_name = trim(host_name)
         this%mappings(idx)%catchem_name = trim(catchem_name)
         this%mappings(idx)%category = trim(category)
      endif

      ! Set optional parameters
      this%mappings(idx)%is_mapped = .true.

      if (present(units)) this%mappings(idx)%units = trim(units)
      if (present(description)) this%mappings(idx)%description = trim(description)
      if (present(is_required)) this%mappings(idx)%is_required = is_required

      if (present(min_value) .and. present(max_value)) then
         this%mappings(idx)%min_value = min_value
         this%mappings(idx)%max_value = max_value
         this%mappings(idx)%check_bounds = .true.
      endif
   end subroutine mapping_add

   !> Remove a field mapping
   subroutine mapping_remove(this, host_name, rc)
      class(FieldMappingType), intent(inout) :: this
      character(len=*), intent(in) :: host_name
      integer, intent(out) :: rc

      integer :: idx, i

      rc = MAPPING_SUCCESS

      idx = this%find_mapping(host_name)
      if (idx <= 0) then
         rc = MAPPING_FIELD_NOT_FOUND
         return
      endif

      ! Shift remaining mappings down
      do i = idx, this%num_mappings - 1
         this%mappings(i) = this%mappings(i + 1)
      end do

      ! Clear the last mapping
      this%mappings(this%num_mappings)%host_name = ''
      this%mappings(this%num_mappings)%catchem_name = ''
      this%mappings(this%num_mappings)%category = ''
      this%mappings(this%num_mappings)%is_mapped = .false.

      this%num_mappings = this%num_mappings - 1
   end subroutine mapping_remove

   !> Clear all mappings
   subroutine mapping_clear_all(this)
      class(FieldMappingType), intent(inout) :: this

      this%num_mappings = 0
      this%mappings(:)%host_name = ''
      this%mappings(:)%catchem_name = ''
      this%mappings(:)%category = ''
      this%mappings(:)%is_mapped = .false.
   end subroutine mapping_clear_all

   !> Find mapping index for a host field name
   function mapping_find(this, host_name) result(idx)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: host_name
      integer :: idx

      integer :: i

      idx = 0
      do i = 1, this%num_mappings
         if (trim(this%mappings(i)%host_name) == trim(host_name)) then
            idx = i
            return
         endif
      end do
   end function mapping_find

   !> Get CATChem field name for a host field name
   function mapping_get_catchem_name(this, host_name, rc) result(catchem_name)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: host_name
      integer, intent(out) :: rc
      character(len=64) :: catchem_name

      integer :: idx

      idx = this%find_mapping(host_name)
      if (idx > 0) then
         catchem_name = this%mappings(idx)%catchem_name
         rc = MAPPING_SUCCESS
      else
         catchem_name = ''
         rc = MAPPING_FIELD_NOT_FOUND
      endif
   end function mapping_get_catchem_name

   !> Get host field name for a CATChem field name
   function mapping_get_host_name(this, catchem_name, rc) result(host_name)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: catchem_name
      integer, intent(out) :: rc
      character(len=64) :: host_name

      integer :: i

      host_name = ''
      rc = MAPPING_FIELD_NOT_FOUND

      do i = 1, this%num_mappings
         if (trim(this%mappings(i)%catchem_name) == trim(catchem_name)) then
            host_name = this%mappings(i)%host_name
            rc = MAPPING_SUCCESS
            return
         endif
      end do
   end function mapping_get_host_name

   !> Get category for a host field name
   function mapping_get_category(this, host_name, rc) result(category)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: host_name
      integer, intent(out) :: rc
      character(len=32) :: category

      integer :: idx

      idx = this%find_mapping(host_name)
      if (idx > 0) then
         category = this%mappings(idx)%category
         rc = MAPPING_SUCCESS
      else
         category = ''
         rc = MAPPING_FIELD_NOT_FOUND
      endif
   end function mapping_get_category

   !> Check if a host field is mapped
   function mapping_is_mapped(this, host_name) result(is_mapped)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: host_name
      logical :: is_mapped

      integer :: idx

      idx = this%find_mapping(host_name)
      is_mapped = (idx > 0 .and. this%mappings(idx)%is_mapped)
   end function mapping_is_mapped

   !> Check if a host field is required
   function mapping_is_required(this, host_name) result(is_required)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: host_name
      logical :: is_required

      integer :: idx

      idx = this%find_mapping(host_name)
      is_required = (idx > 0 .and. this%mappings(idx)%is_required)
   end function mapping_is_required

   !> Get all mappings
   subroutine mapping_get_all(this, mappings, count)
      class(FieldMappingType), intent(in) :: this
      type(FieldMappingEntryType), allocatable, intent(out) :: mappings(:)
      integer, intent(out) :: count

      count = this%num_mappings
      if (count > 0) then
         allocate(mappings(count))
         mappings(1:count) = this%mappings(1:count)
      endif
   end subroutine mapping_get_all

   !> Get mappings by category
   subroutine mapping_get_by_category(this, category, mappings, count)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: category
      type(FieldMappingEntryType), allocatable, intent(out) :: mappings(:)
      integer, intent(out) :: count

      integer :: i, j

      ! Count mappings in category
      count = 0
      do i = 1, this%num_mappings
         if (trim(this%mappings(i)%category) == trim(category)) then
            count = count + 1
         endif
      end do

      ! Extract mappings
      if (count > 0) then
         allocate(mappings(count))
         j = 0
         do i = 1, this%num_mappings
            if (trim(this%mappings(i)%category) == trim(category)) then
               j = j + 1
               mappings(j) = this%mappings(i)
            endif
         end do
      endif
   end subroutine mapping_get_by_category

   !> Load mappings from a YAML configuration file
   subroutine mapping_load_from_file(this, filename, rc)
      class(FieldMappingType), intent(inout) :: this
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      ! TODO: Implement YAML parsing for field mappings
      ! This would read a configuration file like:
      ! field_mappings:
      !   meteorology:
      !     - host_name: "host_temp"
      !       catchem_name: "temperature"
      !       units: "K"
      !       required: true
      !   chemistry:
      !     - host_name: "host_o3"
      !       catchem_name: "O3"
      !       units: "mol/mol"

      rc = MAPPING_SUCCESS
      ! Placeholder implementation
   end subroutine mapping_load_from_file

   !> Save mappings to a YAML configuration file
   subroutine mapping_save_to_file(this, filename, rc)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      rc = MAPPING_SUCCESS
      ! Placeholder implementation
   end subroutine mapping_save_to_file

   !> Validate a specific mapping
   function mapping_validate(this, host_name, data, rc) result(is_valid)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: host_name
      real(fp), intent(in) :: data(:,:,:)  ! Generic 3D data array
      integer, intent(out) :: rc
      logical :: is_valid

      integer :: idx

      rc = MAPPING_SUCCESS
      is_valid = .true.

      idx = this%find_mapping(host_name)
      if (idx <= 0) then
         rc = MAPPING_FIELD_NOT_FOUND
         is_valid = .false.
         return
      endif

      ! Check bounds if enabled
      if (this%mappings(idx)%check_bounds) then
         if (any(data < this%mappings(idx)%min_value) .or. &
            any(data > this%mappings(idx)%max_value)) then
            is_valid = .false.
         endif
      endif
   end function mapping_validate

   !> Validate all mappings
   subroutine mapping_validate_all(this, rc, error_message)
      class(FieldMappingType), intent(in) :: this
      integer, intent(out) :: rc
      character(len=*), intent(out), optional :: error_message

      integer :: i
      logical :: all_valid

      rc = MAPPING_SUCCESS
      all_valid = .true.

      do i = 1, this%num_mappings
         if (.not. is_valid_category(this%mappings(i)%category)) then
            all_valid = .false.
            rc = MAPPING_INVALID_CATEGORY
            if (present(error_message)) then
               error_message = 'Invalid category: ' // trim(this%mappings(i)%category)
            endif
            return
         endif

         if (len_trim(this%mappings(i)%host_name) == 0 .or. &
            len_trim(this%mappings(i)%catchem_name) == 0) then
            all_valid = .false.
            rc = MAPPING_FAILURE
            if (present(error_message)) then
               error_message = 'Empty field names not allowed'
            endif
            return
         endif
      end do
   end subroutine mapping_validate_all

   !> Check that all required fields are mapped
   subroutine mapping_check_required(this, missing_fields, count, rc)
      class(FieldMappingType), intent(in) :: this
      character(len=64), allocatable, intent(out) :: missing_fields(:)
      integer, intent(out) :: count, rc

      integer :: i, j

      rc = MAPPING_SUCCESS

      ! Count missing required fields
      count = 0
      do i = 1, this%num_mappings
         if (this%mappings(i)%is_required .and. .not. this%mappings(i)%is_mapped) then
            count = count + 1
         endif
      end do

      ! Extract missing field names
      if (count > 0) then
         allocate(missing_fields(count))
         j = 0
         do i = 1, this%num_mappings
            if (this%mappings(i)%is_required .and. .not. this%mappings(i)%is_mapped) then
               j = j + 1
               missing_fields(j) = this%mappings(i)%host_name
            endif
         end do
         rc = MAPPING_FAILURE
      endif
   end subroutine mapping_check_required

   !> Print summary of all mappings
   subroutine mapping_print_summary(this)
      class(FieldMappingType), intent(in) :: this

      integer :: i

      print *, 'Field Mapping Summary:'
      print *, '====================='
      print *, 'Total mappings: ', this%num_mappings
      print *, ''

      do i = 1, this%num_mappings
         print *, 'Mapping ', i, ':'
         print *, '  Host name:    ', trim(this%mappings(i)%host_name)
         print *, '  CATChem name: ', trim(this%mappings(i)%catchem_name)
         print *, '  Category:     ', trim(this%mappings(i)%category)
         print *, '  Required:     ', this%mappings(i)%is_required
         print *, '  Units:        ', trim(this%mappings(i)%units)
         print *, ''
      end do
   end subroutine mapping_print_summary

   !> Get mapping statistics
   subroutine mapping_get_stats(this, total, by_category, categories)
      class(FieldMappingType), intent(in) :: this
      integer, intent(out) :: total
      integer, allocatable, intent(out) :: by_category(:)
      character(len=32), allocatable, intent(out) :: categories(:)

      character(len=32) :: unique_cats(this%max_mappings)
      integer :: cat_counts(this%max_mappings)
      integer :: num_cats, i, j, cat_idx

      total = this%num_mappings
      num_cats = 0
      unique_cats = ''
      cat_counts = 0

      ! Count by category
      do i = 1, this%num_mappings
         cat_idx = 0
         do j = 1, num_cats
            if (trim(unique_cats(j)) == trim(this%mappings(i)%category)) then
               cat_idx = j
               exit
            endif
         end do

         if (cat_idx == 0) then
            num_cats = num_cats + 1
            unique_cats(num_cats) = this%mappings(i)%category
            cat_counts(num_cats) = 1
         else
            cat_counts(cat_idx) = cat_counts(cat_idx) + 1
         endif
      end do

      if (num_cats > 0) then
         allocate(by_category(num_cats))
         allocate(categories(num_cats))
         by_category(1:num_cats) = cat_counts(1:num_cats)
         categories(1:num_cats) = unique_cats(1:num_cats)
      endif
   end subroutine mapping_get_stats

   !> Check if a category is valid
   function is_valid_category(category) result(is_valid)
      character(len=*), intent(in) :: category
      logical :: is_valid

      is_valid = (trim(category) == CATEGORY_METEO .or. &
         trim(category) == CATEGORY_CHEMISTRY .or. &
         trim(category) == CATEGORY_EMISSIONS .or. &
         trim(category) == CATEGORY_DIAGNOSTICS)
   end function is_valid_category

end module FieldMapping_Mod
