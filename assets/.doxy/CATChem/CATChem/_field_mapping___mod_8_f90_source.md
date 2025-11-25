

# File FieldMapping\_Mod.F90

[**File List**](files.md) **>** [**api**](dir_da61e3e9a357748887e3ca8d7c5a0c16.md) **>** [**FieldMapping\_Mod.F90**](_field_mapping___mod_8_f90.md)

[Go to the documentation of this file](_field_mapping___mod_8_f90.md)


```Fortran

module fieldmapping_mod
   use precision_mod
   use error_mod

   implicit none
   private

   public :: fieldmappingtype
   public :: fieldmappingentrytype
   public :: mapping_success, mapping_failure
   public :: mapping_field_not_found, mapping_invalid_category

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

   type :: fieldmappingentrytype
      character(len=64) :: host_name = ''
      character(len=64) :: catchem_name = ''
      character(len=32) :: category = ''
      character(len=32) :: units = ''
      character(len=128) :: description = ''
      logical :: is_required = .false.         
      logical :: is_mapped = .false.           
      integer :: mapping_direction = mapping_host_to_catchem  

      ! Field dimensionality and shape information
      integer :: host_dims = 0
      integer :: catchem_dims = 0
      integer :: host_shape(3) = 0
      integer :: catchem_shape(3) = 0
      logical :: needs_reshape = .false.       

      ! Dimension mapping for flexible reshaping
      ! host_dim_map(i) indicates which CATChem dimension corresponds to host dimension i
      integer :: host_dim_map(3) = [1, 2, 3]

      ! For specialized mappings (e.g., nHoriz -> nx*ny, nLev -> nz)
      logical :: is_column_data = .false.      
      logical :: expand_horizontal = .false.   

      ! Data validation parameters
      real(fp) :: min_value = -huge(1.0_fp)    
      real(fp) :: max_value = huge(1.0_fp)     
      logical :: check_bounds = .false.        
   end type fieldmappingentrytype

   type :: fieldmappingtype
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
   end type fieldmappingtype

contains

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

   subroutine mapping_cleanup(this)
      class(FieldMappingType), intent(inout) :: this

      call this%clear_all()
      this%is_initialized = .false.
   end subroutine mapping_cleanup

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

      rc = mapping_success

      ! Check if initialized
      if (.not. this%is_initialized) then
         rc = mapping_failure
         return
      endif

      ! Check if we have space
      if (this%num_mappings >= this%max_mappings) then
         rc = mapping_failure
         return
      endif

      ! Validate category
      if (.not. is_valid_category(category)) then
         rc = mapping_invalid_category
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

   subroutine mapping_remove(this, host_name, rc)
      class(FieldMappingType), intent(inout) :: this
      character(len=*), intent(in) :: host_name
      integer, intent(out) :: rc

      integer :: idx, i

      rc = mapping_success

      idx = this%find_mapping(host_name)
      if (idx <= 0) then
         rc = mapping_field_not_found
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

   subroutine mapping_clear_all(this)
      class(FieldMappingType), intent(inout) :: this

      this%num_mappings = 0
      this%mappings(:)%host_name = ''
      this%mappings(:)%catchem_name = ''
      this%mappings(:)%category = ''
      this%mappings(:)%is_mapped = .false.
   end subroutine mapping_clear_all

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

   function mapping_get_catchem_name(this, host_name, rc) result(catchem_name)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: host_name
      integer, intent(out) :: rc
      character(len=64) :: catchem_name

      integer :: idx

      idx = this%find_mapping(host_name)
      if (idx > 0) then
         catchem_name = this%mappings(idx)%catchem_name
         rc = mapping_success
      else
         catchem_name = ''
         rc = mapping_field_not_found
      endif
   end function mapping_get_catchem_name

   function mapping_get_host_name(this, catchem_name, rc) result(host_name)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: catchem_name
      integer, intent(out) :: rc
      character(len=64) :: host_name

      integer :: i

      host_name = ''
      rc = mapping_field_not_found

      do i = 1, this%num_mappings
         if (trim(this%mappings(i)%catchem_name) == trim(catchem_name)) then
            host_name = this%mappings(i)%host_name
            rc = mapping_success
            return
         endif
      end do
   end function mapping_get_host_name

   function mapping_get_category(this, host_name, rc) result(category)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: host_name
      integer, intent(out) :: rc
      character(len=32) :: category

      integer :: idx

      idx = this%find_mapping(host_name)
      if (idx > 0) then
         category = this%mappings(idx)%category
         rc = mapping_success
      else
         category = ''
         rc = mapping_field_not_found
      endif
   end function mapping_get_category

   function mapping_is_mapped(this, host_name) result(is_mapped)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: host_name
      logical :: is_mapped

      integer :: idx

      idx = this%find_mapping(host_name)
      is_mapped = (idx > 0 .and. this%mappings(idx)%is_mapped)
   end function mapping_is_mapped

   function mapping_is_required(this, host_name) result(is_required)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: host_name
      logical :: is_required

      integer :: idx

      idx = this%find_mapping(host_name)
      is_required = (idx > 0 .and. this%mappings(idx)%is_required)
   end function mapping_is_required

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

      rc = mapping_success
      ! Placeholder implementation
   end subroutine mapping_load_from_file

   subroutine mapping_save_to_file(this, filename, rc)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      rc = mapping_success
      ! Placeholder implementation
   end subroutine mapping_save_to_file

   function mapping_validate(this, host_name, data, rc) result(is_valid)
      class(FieldMappingType), intent(in) :: this
      character(len=*), intent(in) :: host_name
      real(fp), intent(in) :: data(:,:,:)  ! Generic 3D data array
      integer, intent(out) :: rc
      logical :: is_valid

      integer :: idx

      rc = mapping_success
      is_valid = .true.

      idx = this%find_mapping(host_name)
      if (idx <= 0) then
         rc = mapping_field_not_found
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

   subroutine mapping_validate_all(this, rc, error_message)
      class(FieldMappingType), intent(in) :: this
      integer, intent(out) :: rc
      character(len=*), intent(out), optional :: error_message

      integer :: i
      logical :: all_valid

      rc = mapping_success
      all_valid = .true.

      do i = 1, this%num_mappings
         if (.not. is_valid_category(this%mappings(i)%category)) then
            all_valid = .false.
            rc = mapping_invalid_category
            if (present(error_message)) then
               error_message = 'Invalid category: ' // trim(this%mappings(i)%category)
            endif
            return
         endif

         if (len_trim(this%mappings(i)%host_name) == 0 .or. &
            len_trim(this%mappings(i)%catchem_name) == 0) then
            all_valid = .false.
            rc = mapping_failure
            if (present(error_message)) then
               error_message = 'Empty field names not allowed'
            endif
            return
         endif
      end do
   end subroutine mapping_validate_all

   subroutine mapping_check_required(this, missing_fields, count, rc)
      class(FieldMappingType), intent(in) :: this
      character(len=64), allocatable, intent(out) :: missing_fields(:)
      integer, intent(out) :: count, rc

      integer :: i, j

      rc = mapping_success

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
         rc = mapping_failure
      endif
   end subroutine mapping_check_required

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

   function is_valid_category(category) result(is_valid)
      character(len=*), intent(in) :: category
      logical :: is_valid

      is_valid = (trim(category) == category_meteo .or. &
         trim(category) == category_chemistry .or. &
         trim(category) == category_emissions .or. &
         trim(category) == category_diagnostics)
   end function is_valid_category

end module fieldmapping_mod
```


