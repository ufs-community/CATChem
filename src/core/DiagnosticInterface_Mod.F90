!> \file DiagnosticInterface_Mod.F90
!! \brief Dynamic diagnostic system interfaces and types
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module provides interfaces and types for the dynamic diagnostic system,
!! allowing processes to register and manage their own diagnostic outputs at runtime.
!!
!! \details
!! The diagnostic system supports:
!! - Multiple data types (scalar, 1D, 2D, 3D arrays)
!! - Flexible metadata (units, description, output frequency)
!! - Process-specific diagnostic registration
!! - Runtime diagnostic query and collection
!! - Optional diagnostic output control
!!
!! \section diag_usage Usage Example
!! \code{.f90}
!! use DiagnosticInterface_Mod
!! type(DiagnosticFieldType) :: diag_field
!! integer :: rc
!!
!! call diag_field%create('dust_flux', 'Total dust emission flux', &
!!                        'kg m-2 s-1', DIAG_REAL_2D, rc)
!! call diag_mgr%register_diagnostic('dust_process', diag_field, rc)
!! \endcode
!!
module DiagnosticInterface_Mod
   use precision_mod, only: fp
   use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE, &
      ERROR_MEMORY_ALLOCATION, ERROR_INVALID_INPUT, &
      ERROR_DUPLICATE_ENTRY, ERROR_NOT_FOUND

   implicit none
   private

   integer, parameter :: max_fields = 100      !< Maximum fields per process

   public :: DiagnosticFieldType, DiagnosticRegistryType
   public :: DiagnosticDataType
   public :: DIAG_REAL_SCALAR, DIAG_REAL_1D, DIAG_REAL_2D, DIAG_REAL_3D
   public :: DIAG_INTEGER_SCALAR, DIAG_INTEGER_1D, DIAG_INTEGER_2D, DIAG_INTEGER_3D
   public :: DIAG_LOGICAL_SCALAR, DIAG_LOGICAL_1D, DIAG_LOGICAL_2D, DIAG_LOGICAL_3D
   public :: DIAG_FREQ_NEVER, DIAG_FREQ_TIMESTEP, DIAG_FREQ_HOURLY, DIAG_FREQ_DAILY, DIAG_FREQ_CUSTOM

   ! Diagnostic data type enumerations
   integer, parameter :: DIAG_REAL_SCALAR = 1
   integer, parameter :: DIAG_REAL_1D = 2
   integer, parameter :: DIAG_REAL_2D = 3
   integer, parameter :: DIAG_REAL_3D = 4
   integer, parameter :: DIAG_INTEGER_SCALAR = 11
   integer, parameter :: DIAG_INTEGER_1D = 12
   integer, parameter :: DIAG_INTEGER_2D = 13
   integer, parameter :: DIAG_INTEGER_3D = 14
   integer, parameter :: DIAG_LOGICAL_SCALAR = 21
   integer, parameter :: DIAG_LOGICAL_1D = 22
   integer, parameter :: DIAG_LOGICAL_2D = 23
   integer, parameter :: DIAG_LOGICAL_3D = 24

   ! Diagnostic frequency enumerations
   integer, parameter :: DIAG_FREQ_NEVER = 0      !< Never output
   integer, parameter :: DIAG_FREQ_TIMESTEP = 1   !< Every timestep
   integer, parameter :: DIAG_FREQ_HOURLY = 2     !< Hourly output
   integer, parameter :: DIAG_FREQ_DAILY = 3      !< Daily output
   integer, parameter :: DIAG_FREQ_CUSTOM = 99    !< Custom frequency

   !> \brief Type for diagnostic data storage
   !!
   !! Union-like type that can hold different data types and dimensions.
   !! Only one component should be allocated at a time.
   type :: DiagnosticDataType
      private
      integer :: data_type = 0                     !< Data type identifier
      logical :: is_allocated = .false.           !< Allocation status

      ! Real data storage
      real(fp) :: real_scalar = 0.0_fp
      real(fp), allocatable :: real_1d(:)
      real(fp), allocatable :: real_2d(:,:)
      real(fp), allocatable :: real_3d(:,:,:)

      ! Integer data storage
      integer :: int_scalar = 0
      integer, allocatable :: int_1d(:)
      integer, allocatable :: int_2d(:,:)
      integer, allocatable :: int_3d(:,:,:)

      ! Logical data storage
      logical :: logical_scalar = .false.
      logical, allocatable :: logical_1d(:)
      logical, allocatable :: logical_2d(:,:)
      logical, allocatable :: logical_3d(:,:,:)

   contains
      procedure :: allocate_data => diag_data_allocate
      procedure :: deallocate_data => diag_data_deallocate
      procedure :: get_data_type => diag_data_get_type
      procedure :: is_data_allocated => diag_data_is_allocated
      procedure :: set_real_scalar => diag_data_set_real_scalar
      procedure :: get_real_scalar => diag_data_get_real_scalar
      procedure :: set_real_1d => diag_data_set_real_1d
      procedure :: get_real_1d_ptr => diag_data_get_real_1d_ptr
      procedure :: set_real_2d => diag_data_set_real_2d
      procedure :: get_real_2d_ptr => diag_data_get_real_2d_ptr
      procedure :: set_real_3d => diag_data_set_real_3d
      procedure :: get_real_3d_ptr => diag_data_get_real_3d_ptr
      final :: diag_data_finalize
   end type DiagnosticDataType

   !> \brief Individual diagnostic field specification
   !!
   !! Contains metadata and data storage for a single diagnostic field.
   type :: DiagnosticFieldType
      private
      character(len=64) :: field_name = ''         !< Diagnostic field name
      character(len=128) :: description = ''       !< Field description
      character(len=32) :: units = ''              !< Physical units
      character(len=64) :: process_name = ''       !< Process that owns this diagnostic
      integer :: data_type = 0                     !< Data type identifier
      integer :: output_frequency = DIAG_FREQ_NEVER !< Output frequency
      real(fp) :: custom_frequency = 0.0_fp        !< Custom output frequency (seconds)
      logical :: is_enabled = .true.               !< Whether diagnostic is enabled
      logical :: is_initialized = .false.          !< Initialization status
      type(DiagnosticDataType) :: data             !< Actual data storage

      ! Diagnostic species filtering support
      character(len=32), allocatable :: diagnostic_species(:)  !< User-defined species for diagnostics
      integer, allocatable :: diagnostic_species_id(:)         !< Indices mapping diagnostic_species to species names

   contains
      procedure :: create => diag_field_create
      procedure :: initialize_data => diag_field_init_data
      procedure :: cleanup => diag_field_cleanup
      procedure :: is_ready => diag_field_is_ready
      procedure :: get_name => diag_field_get_name
      procedure :: get_description => diag_field_get_description
      procedure :: get_units => diag_field_get_units
      procedure :: get_process_name => diag_field_get_process_name
      procedure :: get_data_type => diag_field_get_data_type
      procedure :: get_data_ptr => diag_field_get_data_ptr
      procedure :: set_enabled => diag_field_set_enabled
      procedure :: get_is_enabled => diag_field_is_enabled
      procedure :: should_output => diag_field_should_output
      procedure :: update_data => diag_field_update_data
      procedure :: reset_data => diag_field_reset_data
      procedure :: validate_field => diag_field_validate_field
      procedure :: get_diagnostic_species => diag_field_get_diagnostic_species
      procedure :: get_diagnostic_species_id => diag_field_get_diagnostic_species_id
      final :: diag_field_finalize
   end type DiagnosticFieldType

   !> \brief Registry of diagnostic fields for a process
   !!
   !! Manages a collection of diagnostic fields for a single process.
   type :: DiagnosticRegistryType
      private
      character(len=64) :: process_name = ''       !< Process name
      integer :: n_fields = 0                     !< Number of registered fields
      type(DiagnosticFieldType) :: fields(max_fields)    !< Diagnostic fields array
      logical :: is_initialized = .false.         !< Initialization status

   contains
      procedure :: init => diag_registry_init
      procedure :: cleanup => diag_registry_cleanup
      procedure :: finalize => diag_registry_finalize
      procedure :: register_field => diag_registry_register
      procedure :: get_field => diag_registry_get_field
      procedure :: get_field_ptr => diag_registry_get_field_ptr
      procedure :: list_fields => diag_registry_list_fields
      procedure :: get_field_count => diag_registry_get_count
      procedure :: get_num_diagnostics => diag_registry_get_num_diagnostics
      procedure :: field_exists => diag_registry_field_exists
      procedure :: enable_field => diag_registry_enable_field
      procedure :: disable_field => diag_registry_disable_field
      procedure :: set_output_frequency => diag_registry_set_frequency
      procedure :: reset => diag_registry_reset
      procedure :: validate => diag_registry_validate
   end type DiagnosticRegistryType

contains

   !> \brief Allocate diagnostic data storage
   !!
   !! \param[inout] this DiagnosticDataType instance
   !! \param[in] data_type Type of data to allocate
   !! \param[in] dims Dimensions for multi-dimensional arrays
   !! \param[out] rc Return code
   subroutine diag_data_allocate(this, data_type, dims, rc)
      class(DiagnosticDataType), intent(inout) :: this
      integer, intent(in) :: data_type
      integer, intent(in), optional :: dims(:)
      integer, intent(out) :: rc
      integer :: alloc_stat

      rc = CC_SUCCESS
      alloc_stat = 0

      ! Deallocate any existing data first
      call this%deallocate_data()

      this%data_type = data_type

      select case (data_type)
       case (DIAG_REAL_SCALAR)
         this%real_scalar = 0.0_fp
         this%is_allocated = .true.

       case (DIAG_REAL_1D)
         if (.not. present(dims) .or. size(dims) < 1) then
            rc = ERROR_INVALID_INPUT
            write(*,*) "Error: Invalid dimensions for DIAG_REAL_1D"
            return
         end if
         allocate(this%real_1d(dims(1)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = ERROR_MEMORY_ALLOCATION
            write(*,*) "Error: Memory allocation failed for DIAG_REAL_1D"
            return
         end if
         this%real_1d = 0.0_fp
         this%is_allocated = .true.

       case (DIAG_REAL_2D)
         if (.not. present(dims) .or. size(dims) < 2) then
            rc = ERROR_INVALID_INPUT
            write(*,*) "Error: Invalid dimensions for DIAG_REAL_2D"
            return
         end if
         allocate(this%real_2d(dims(1), dims(2)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = ERROR_MEMORY_ALLOCATION
            write(*,*) "Error: Memory allocation failed for DIAG_REAL_2D"
            return
         end if
         this%real_2d = 0.0_fp
         this%is_allocated = .true.

       case (DIAG_REAL_3D)
         if (.not. present(dims) .or. size(dims) < 3) then
            rc = ERROR_INVALID_INPUT
            write(*,*) "Error: Invalid dimensions for DIAG_REAL_3D"
            return
         end if
         allocate(this%real_3d(dims(1), dims(2), dims(3)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = ERROR_MEMORY_ALLOCATION
            write(*,*) "Error: Memory allocation failed for DIAG_REAL_3D"
            return
         end if
         this%real_3d = 0.0_fp
         this%is_allocated = .true.

       case (DIAG_INTEGER_SCALAR)
         this%int_scalar = 0
         this%is_allocated = .true.

       case (DIAG_INTEGER_1D)
         if (.not. present(dims) .or. size(dims) < 1) then
            rc = ERROR_INVALID_INPUT
            write(*,*) "Error: Invalid dimensions for DIAG_INTEGER_1D"
            return
         end if
         allocate(this%int_1d(dims(1)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = ERROR_MEMORY_ALLOCATION
            write(*,*) "Error: Memory allocation failed for DIAG_INTEGER_1D"
            return
         end if
         this%int_1d = 0
         this%is_allocated = .true.

       case (DIAG_INTEGER_2D)
         if (.not. present(dims) .or. size(dims) < 2) then
            rc = ERROR_INVALID_INPUT
            write(*,*) "Error: Invalid dimensions for DIAG_INTEGER_2D"
            return
         end if
         allocate(this%int_2d(dims(1), dims(2)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = ERROR_MEMORY_ALLOCATION
            write(*,*) "Error: Memory allocation failed for DIAG_INTEGER_2D"
            return
         end if
         this%int_2d = 0
         this%is_allocated = .true.

       case (DIAG_INTEGER_3D)
         if (.not. present(dims) .or. size(dims) < 3) then
            rc = ERROR_INVALID_INPUT
            write(*,*) "Error: Invalid dimensions for DIAG_INTEGER_3D"
            return
         end if
         allocate(this%int_3d(dims(1), dims(2), dims(3)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = ERROR_MEMORY_ALLOCATION
            write(*,*) "Error: Memory allocation failed for DIAG_INTEGER_3D"
            return
         end if
         this%int_3d = 0
         this%is_allocated = .true.

       case (DIAG_LOGICAL_SCALAR)
         this%logical_scalar = .false.
         this%is_allocated = .true.

       case (DIAG_LOGICAL_1D)
         if (.not. present(dims) .or. size(dims) < 1) then
            rc = ERROR_INVALID_INPUT
            write(*,*) "Error: Invalid dimensions for DIAG_LOGICAL_1D"
            return
         end if
         allocate(this%logical_1d(dims(1)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = ERROR_MEMORY_ALLOCATION
            write(*,*) "Error: Memory allocation failed for DIAG_LOGICAL_1D"
            return
         end if
         this%logical_1d = .false.
         this%is_allocated = .true.

       case (DIAG_LOGICAL_2D)
         if (.not. present(dims) .or. size(dims) < 2) then
            rc = ERROR_INVALID_INPUT
            write(*,*) "Error: Invalid dimensions for DIAG_LOGICAL_2D"
            return
         end if
         allocate(this%logical_2d(dims(1), dims(2)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = ERROR_MEMORY_ALLOCATION
            write(*,*) "Error: Memory allocation failed for DIAG_LOGICAL_2D"
            return
         end if
         this%logical_2d = .false.
         this%is_allocated = .true.

       case (DIAG_LOGICAL_3D)
         if (.not. present(dims) .or. size(dims) < 3) then
            rc = ERROR_INVALID_INPUT
            write(*,*) "Error: Invalid dimensions for DIAG_LOGICAL_3D"
            return
         end if
         allocate(this%logical_3d(dims(1), dims(2), dims(3)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = ERROR_MEMORY_ALLOCATION
            write(*,*) "Error: Memory allocation failed for DIAG_LOGICAL_3D"
            return
         end if
         this%logical_3d = .false.
         this%is_allocated = .true.

       case default
         rc = ERROR_INVALID_INPUT
         write(*,*) "Error: Unsupported data type"
      end select

   end subroutine diag_data_allocate

   !> \brief Deallocate diagnostic data storage
   subroutine diag_data_deallocate(this)
      class(DiagnosticDataType), intent(inout) :: this

      if (allocated(this%real_1d)) deallocate(this%real_1d)
      if (allocated(this%real_2d)) deallocate(this%real_2d)
      if (allocated(this%real_3d)) deallocate(this%real_3d)
      if (allocated(this%int_1d)) deallocate(this%int_1d)
      if (allocated(this%int_2d)) deallocate(this%int_2d)
      if (allocated(this%int_3d)) deallocate(this%int_3d)
      if (allocated(this%logical_1d)) deallocate(this%logical_1d)
      if (allocated(this%logical_2d)) deallocate(this%logical_2d)
      if (allocated(this%logical_3d)) deallocate(this%logical_3d)

      this%is_allocated = .false.
      this%data_type = 0

   end subroutine diag_data_deallocate

   !> \brief Get data type identifier
   function diag_data_get_type(this) result(data_type)
      class(DiagnosticDataType), intent(in) :: this
      integer :: data_type
      data_type = this%data_type
   end function diag_data_get_type

   !> \brief Check if data is allocated
   function diag_data_is_allocated(this) result(is_allocated)
      class(DiagnosticDataType), intent(in) :: this
      logical :: is_allocated
      is_allocated = this%is_allocated
   end function diag_data_is_allocated

   !> \brief Set real scalar value
   subroutine diag_data_set_real_scalar(this, value)
      class(DiagnosticDataType), intent(inout) :: this
      real(fp), intent(in) :: value
      if (this%data_type == DIAG_REAL_SCALAR) then
         this%real_scalar = value
      end if
   end subroutine diag_data_set_real_scalar

   !> \brief Get real scalar value
   function diag_data_get_real_scalar(this) result(value)
      class(DiagnosticDataType), intent(in) :: this
      real(fp) :: value
      if (this%data_type == DIAG_REAL_SCALAR) then
         value = this%real_scalar
      else
         value = 0.0_fp
      end if
   end function diag_data_get_real_scalar

   !> \brief Set real 1D array values
   subroutine diag_data_set_real_1d(this, values)
      class(DiagnosticDataType), intent(inout) :: this
      real(fp), intent(in) :: values(:)
      if (this%data_type == DIAG_REAL_1D .and. allocated(this%real_1d)) then
         if (size(values) == size(this%real_1d)) then
            this%real_1d = values
         end if
      end if
   end subroutine diag_data_set_real_1d

   !> \brief Get pointer to real 1D array
   function diag_data_get_real_1d_ptr(this) result(ptr)
      class(DiagnosticDataType), intent(in), target :: this
      real(fp), pointer :: ptr(:)
      if (this%data_type == DIAG_REAL_1D .and. allocated(this%real_1d)) then
         ptr => this%real_1d
      else
         nullify(ptr)
      end if
   end function diag_data_get_real_1d_ptr

   !> \brief Set real 2D array values
   subroutine diag_data_set_real_2d(this, values)
      class(DiagnosticDataType), intent(inout) :: this
      real(fp), intent(in) :: values(:,:)
      if (this%data_type == DIAG_REAL_2D .and. allocated(this%real_2d)) then
         if (size(values,1) == size(this%real_2d,1) .and. size(values,2) == size(this%real_2d,2)) then
            this%real_2d = values
         end if
      end if
   end subroutine diag_data_set_real_2d

   !> \brief Get pointer to real 2D array
   function diag_data_get_real_2d_ptr(this) result(ptr)
      class(DiagnosticDataType), intent(in), target :: this
      real(fp), pointer :: ptr(:,:)
      if (this%data_type == DIAG_REAL_2D .and. allocated(this%real_2d)) then
         ptr => this%real_2d
      else
         nullify(ptr)
      end if
   end function diag_data_get_real_2d_ptr

   !> \brief Set real 3D array values
   subroutine diag_data_set_real_3d(this, values)
      class(DiagnosticDataType), intent(inout) :: this
      real(fp), intent(in) :: values(:,:,:)
      if (this%data_type == DIAG_REAL_3D .and. allocated(this%real_3d)) then
         if (size(values,1) == size(this%real_3d,1) .and. &
            size(values,2) == size(this%real_3d,2) .and. &
            size(values,3) == size(this%real_3d,3)) then
            this%real_3d = values
         end if
      end if
   end subroutine diag_data_set_real_3d

   !> \brief Get pointer to real 3D array
   function diag_data_get_real_3d_ptr(this) result(ptr)
      class(DiagnosticDataType), intent(in), target :: this
      real(fp), pointer :: ptr(:,:,:)
      if (this%data_type == DIAG_REAL_3D .and. allocated(this%real_3d)) then
         ptr => this%real_3d
      else
         nullify(ptr)
      end if
   end function diag_data_get_real_3d_ptr

   !> \brief Finalizer for DiagnosticDataType
   subroutine diag_data_finalize(this)
      type(DiagnosticDataType), intent(inout) :: this
      call this%deallocate_data()
   end subroutine diag_data_finalize

   !> \brief Create and configure a diagnostic field
   !!
   !! \param[inout] this DiagnosticFieldType instance
   !! \param[in] field_name Name of the diagnostic field
   !! \param[in] description Human-readable description
   !! \param[in] units Physical units
   !! \param[in] data_type Type of data (scalar, 1D, 2D, 3D)
   !! \param[in] process_name Name of the process owning this diagnostic
   !! \param[in] diagnostic_species Optional array of species names for diagnostics
   !! \param[in] diagnostic_species_id Optional array of species indices for diagnostics
   !! \param[out] rc Return code
   subroutine diag_field_create(this, field_name, description, units, data_type, process_name, &
      diagnostic_species, diagnostic_species_id, rc)
      class(DiagnosticFieldType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      character(len=*), intent(in) :: description
      character(len=*), intent(in) :: units
      integer, intent(in) :: data_type
      character(len=*), intent(in), optional :: process_name
      character(len=*), intent(in), optional :: diagnostic_species(:)
      integer, intent(in), optional :: diagnostic_species_id(:)
      integer, intent(out) :: rc

      ! No local variables needed

      rc = CC_SUCCESS

      ! Validate inputs
      if (len_trim(field_name) == 0) then
         rc = ERROR_INVALID_INPUT
         return
      end if

      ! Set field metadata
      this%field_name = trim(field_name)
      this%description = trim(description)
      this%units = trim(units)
      this%data_type = data_type
      if (present(process_name)) then
         this%process_name = trim(process_name)
      end if

      ! Store diagnostic species arrays if provided
      if (present(diagnostic_species)) then
         allocate(this%diagnostic_species(size(diagnostic_species)))
         this%diagnostic_species = diagnostic_species
      end if

      if (present(diagnostic_species_id)) then
         allocate(this%diagnostic_species_id(size(diagnostic_species_id)))
         this%diagnostic_species_id = diagnostic_species_id
      end if

      ! Set defaults
      this%output_frequency = DIAG_FREQ_TIMESTEP
      this%is_enabled = .true.
      this%is_initialized = .true.

      ! Initialize data storage for all types with default dimensions for arrays
      select case (data_type)
       case (DIAG_REAL_SCALAR, DIAG_INTEGER_SCALAR, DIAG_LOGICAL_SCALAR)
         ! For scalar types, we can initialize without dimensions
         call this%data%allocate_data(this%data_type, rc=rc)
         if (rc /= CC_SUCCESS) then
            write(*,*) "Error: Could not allocate data for diagnostic field"
            this%is_initialized = .false.
            return
         end if
       case (DIAG_REAL_1D, DIAG_INTEGER_1D, DIAG_LOGICAL_1D)
         ! For 1D arrays, initialize with default size of 1
         call this%data%allocate_data(this%data_type, [1], rc=rc)
         if (rc /= CC_SUCCESS) then
            write(*,*) "Error: Could not allocate data for diagnostic field"
            this%is_initialized = .false.
            return
         end if
       case (DIAG_REAL_2D, DIAG_INTEGER_2D, DIAG_LOGICAL_2D)
         ! For 2D arrays, initialize with default size of 1x1
         call this%data%allocate_data(this%data_type, [1, 1], rc=rc)
         if (rc /= CC_SUCCESS) then
            write(*,*) "Error: Could not allocate data for diagnostic field"
            this%is_initialized = .false.
            return
         end if
       case (DIAG_REAL_3D, DIAG_INTEGER_3D, DIAG_LOGICAL_3D)
         ! For 3D arrays, initialize with default size of 1x1x1
         call this%data%allocate_data(this%data_type, [1, 1, 1], rc=rc)
         if (rc /= CC_SUCCESS) then
            write(*,*) "Error: Could not allocate data for diagnostic field"
            this%is_initialized = .false.
            return
         end if
       case default
         rc = ERROR_INVALID_INPUT
         this%is_initialized = .false.
         return
      end select

      ! Verify that data was allocated successfully
      if (.not. this%data%is_data_allocated()) then
         rc = ERROR_MEMORY_ALLOCATION
         write(*,*) "Error: Data allocation failed during create()"
         this%is_initialized = .false.
         return
      end if

   end subroutine diag_field_create

   !> \brief Initialize data storage for diagnostic field
   !!
   !! \param[inout] this DiagnosticFieldType instance
   !! \param[in] dims Dimensions for array data types
   !! \param[out] rc Return code
   subroutine diag_field_init_data(this, dims, rc)
      class(DiagnosticFieldType), intent(inout) :: this
      integer, intent(in), optional :: dims(:)
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      if (.not. this%is_initialized) then
         rc = ERROR_INVALID_INPUT
         return
      end if

      call this%data%allocate_data(this%data_type, dims, rc)

   end subroutine diag_field_init_data

   !> \brief Clean up diagnostic field
   subroutine diag_field_cleanup(this)
      class(DiagnosticFieldType), intent(inout) :: this
      call this%data%deallocate_data()
      if (allocated(this%diagnostic_species)) deallocate(this%diagnostic_species)
      if (allocated(this%diagnostic_species_id)) deallocate(this%diagnostic_species_id)
      this%is_initialized = .false.
   end subroutine diag_field_cleanup

   !> \brief Check if diagnostic field is ready for use
   function diag_field_is_ready(this) result(is_ready)
      class(DiagnosticFieldType), intent(in) :: this
      logical :: is_ready
      is_ready = this%is_initialized .and. this%data%is_data_allocated()
   end function diag_field_is_ready

   !> \brief Get field name
   function diag_field_get_name(this) result(name)
      class(DiagnosticFieldType), intent(in) :: this
      character(len=64) :: name
      name = this%field_name
   end function diag_field_get_name

   !> \brief Get field description
   function diag_field_get_description(this) result(description)
      class(DiagnosticFieldType), intent(in) :: this
      character(len=128) :: description
      description = this%description
   end function diag_field_get_description

   !> \brief Get field units
   function diag_field_get_units(this) result(units)
      class(DiagnosticFieldType), intent(in) :: this
      character(len=32) :: units
      units = this%units
   end function diag_field_get_units

   !> \brief Get process name
   function diag_field_get_process_name(this) result(process_name)
      class(DiagnosticFieldType), intent(in) :: this
      character(len=64) :: process_name
      process_name = this%process_name
   end function diag_field_get_process_name

   !> \brief Get data type
   function diag_field_get_data_type(this) result(data_type)
      class(DiagnosticFieldType), intent(in) :: this
      integer :: data_type
      data_type = this%data_type
   end function diag_field_get_data_type

   !> \brief Get pointer to data storage
   function diag_field_get_data_ptr(this) result(data_ptr)
      class(DiagnosticFieldType), intent(in), target :: this
      type(DiagnosticDataType), pointer :: data_ptr
      data_ptr => this%data
   end function diag_field_get_data_ptr

   !> \brief Get diagnostic species array
   function diag_field_get_diagnostic_species(this) result(diagnostic_species)
      class(DiagnosticFieldType), intent(in) :: this
      character(len=32), allocatable :: diagnostic_species(:)
      if (allocated(this%diagnostic_species)) then
         allocate(diagnostic_species(size(this%diagnostic_species)))
         diagnostic_species = this%diagnostic_species
      else
         allocate(diagnostic_species(0))
      end if
   end function diag_field_get_diagnostic_species

   !> \brief Get diagnostic species ID array
   function diag_field_get_diagnostic_species_id(this) result(diagnostic_species_id)
      class(DiagnosticFieldType), intent(in) :: this
      integer, allocatable :: diagnostic_species_id(:)
      if (allocated(this%diagnostic_species_id)) then
         allocate(diagnostic_species_id(size(this%diagnostic_species_id)))
         diagnostic_species_id = this%diagnostic_species_id
      else
         allocate(diagnostic_species_id(0))
      end if
   end function diag_field_get_diagnostic_species_id

   !> \brief Enable/disable diagnostic field
   subroutine diag_field_set_enabled(this, enabled)
      class(DiagnosticFieldType), intent(inout) :: this
      logical, intent(in) :: enabled
      this%is_enabled = enabled
   end subroutine diag_field_set_enabled

   !> \brief Check if field is enabled
   function diag_field_is_enabled(this) result(is_enabled)
      class(DiagnosticFieldType), intent(in) :: this
      logical :: is_enabled
      is_enabled = this%is_enabled
   end function diag_field_is_enabled

   !> \brief Check if field should be output at current time
   function diag_field_should_output(this, current_time, dt) result(should_output)
      class(DiagnosticFieldType), intent(in) :: this
      real(fp), intent(in) :: current_time
      real(fp), intent(in) :: dt
      logical :: should_output

      should_output = .false.

      if (.not. this%is_enabled) return

      select case (this%output_frequency)
       case (DIAG_FREQ_NEVER)
         should_output = .false.
       case (DIAG_FREQ_TIMESTEP)
         should_output = .true.
       case (DIAG_FREQ_HOURLY)
         should_output = (mod(current_time, 3600.0_fp) < dt)
       case (DIAG_FREQ_DAILY)
         should_output = (mod(current_time, 86400.0_fp) < dt)
       case (DIAG_FREQ_CUSTOM)
         if (this%custom_frequency > 0.0_fp) then
            should_output = (mod(current_time, this%custom_frequency) < dt)
         end if
      end select

   end function diag_field_should_output

   !> \brief Update diagnostic field data
   !!
   !! This is a generic interface that processes can use to update their diagnostics
   subroutine diag_field_update_data(this, scalar_val, array_1d, array_2d, array_3d)
      class(DiagnosticFieldType), intent(inout) :: this
      real(fp), intent(in), optional :: scalar_val
      real(fp), intent(in), optional :: array_1d(:)
      real(fp), intent(in), optional :: array_2d(:,:)
      real(fp), intent(in), optional :: array_3d(:,:,:)

      if (.not. this%is_ready()) return

      if (present(scalar_val)) then
         call this%data%set_real_scalar(scalar_val)
      else if (present(array_1d)) then
         call this%data%set_real_1d(array_1d)
      else if (present(array_2d)) then
         call this%data%set_real_2d(array_2d)
      else if (present(array_3d)) then
         call this%data%set_real_3d(array_3d)
      end if

   end subroutine diag_field_update_data

   !> \brief Finalizer for DiagnosticFieldType
   subroutine diag_field_finalize(this)
      type(DiagnosticFieldType), intent(inout) :: this
      call this%cleanup()
   end subroutine diag_field_finalize

   !> \brief Reset diagnostic field data to zero/default values
   !!
   !! \param[inout] this DiagnosticFieldType instance
   !! \param[out] rc Return code
   subroutine diag_field_reset_data(this, rc)
      class(DiagnosticFieldType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      if (.not. this%is_ready()) then
         rc = CC_FAILURE
         return
      endif

      ! Reset data based on type
      select case (this%data%get_data_type())
       case (DIAG_REAL_SCALAR)
         call this%data%set_real_scalar(0.0_fp)
       case (DIAG_REAL_1D)
         ! Would set 1D array to zero
       case (DIAG_REAL_2D)
         ! Would set 2D array to zero
       case (DIAG_REAL_3D)
         ! Would set 3D array to zero
      end select

   end subroutine diag_field_reset_data

   !> \brief Validate diagnostic field
   !!
   !! \param[in] this DiagnosticFieldType instance
   !! \param[inout] error_mgr Error manager for reporting
   !! \param[out] rc Return code
   subroutine diag_field_validate_field(this, error_mgr, rc)
      class(DiagnosticFieldType), intent(in) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Validate field name
      if (len_trim(this%field_name) == 0) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, 'Diagnostic field name is empty', rc)
         return
      endif

      ! Validate data type
      if (this%data_type <= 0) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, 'Invalid diagnostic data type', rc)
         return
      endif

      ! Validate initialization
      if (.not. this%is_initialized) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, 'Diagnostic field not initialized', rc)
         return
      endif

   end subroutine diag_field_validate_field

   !-------------------
   ! DiagnosticRegistryType procedure implementations
   !-------------------

   subroutine diag_registry_init(this, process_name, error_mgr, rc)
      class(DiagnosticRegistryType), intent(inout) :: this
      character(len=*), intent(in), optional :: process_name
      type(ErrorManagerType), pointer, intent(inout), optional :: error_mgr
      integer, intent(out), optional :: rc

      if (present(rc)) rc = 0
      this%process_name = ''
      if (present(process_name)) this%process_name = trim(process_name)
      this%n_fields = 0
      call this%cleanup()  ! Clean up any previous state
      this%is_initialized = .true.
   end subroutine diag_registry_init

   subroutine diag_registry_cleanup(this)
      class(DiagnosticRegistryType), intent(inout) :: this
      integer :: i
      do i = 1, this%n_fields
         call this%fields(i)%cleanup()
      end do
      this%n_fields = 0
      this%is_initialized = .false.
   end subroutine diag_registry_cleanup

   !> \brief Finalize diagnostic registry (alias for cleanup)
   subroutine diag_registry_finalize(this, rc)
      class(DiagnosticRegistryType), intent(inout) :: this
      integer, intent(out) :: rc
      rc = 0
      call this%cleanup()
   end subroutine diag_registry_finalize

   subroutine diag_registry_register(this, field, rc)
      class(DiagnosticRegistryType), intent(inout) :: this
      type(DiagnosticFieldType), intent(in) :: field
      integer, intent(out) :: rc
      integer :: i
      rc = CC_SUCCESS
      if (.not. this%is_initialized) then
         rc = ERROR_INVALID_INPUT
         return
      end if
      ! Check if field is valid
      if (.not. field%is_ready()) then
         rc = ERROR_INVALID_INPUT
         return
      end if
      ! Check for duplicate
      do i = 1, this%n_fields
         if (trim(this%fields(i)%field_name) == trim(field%field_name)) then
            rc = ERROR_DUPLICATE_ENTRY
            return
         end if
      end do
      if (this%n_fields >= max_fields) then
         rc = ERROR_MEMORY_ALLOCATION
         return
      end if
      this%n_fields = this%n_fields + 1
      this%fields(this%n_fields) = field
   end subroutine diag_registry_register

   function diag_registry_get_field(this, name) result(field)
      class(DiagnosticRegistryType), intent(in) :: this
      character(len=*), intent(in) :: name
      type(DiagnosticFieldType) :: field
      integer :: i
      field = DiagnosticFieldType()
      do i = 1, this%n_fields
         if (trim(this%fields(i)%field_name) == trim(name)) then
            field = this%fields(i)
            return
         end if
      end do
   end function diag_registry_get_field

   function diag_registry_get_field_ptr(this, name) result(field_ptr)
      class(DiagnosticRegistryType), intent(in), target :: this
      character(len=*), intent(in) :: name
      type(DiagnosticFieldType), pointer :: field_ptr
      integer :: i
      nullify(field_ptr)
      do i = 1, this%n_fields
         if (trim(this%fields(i)%field_name) == trim(name)) then
            field_ptr => this%fields(i)
            return
         end if
      end do
   end function diag_registry_get_field_ptr

   subroutine diag_registry_list_fields(this, names, n)
      class(DiagnosticRegistryType), intent(in) :: this
      character(len=*), intent(out) :: names(:)
      integer, intent(out) :: n
      integer :: i

      ! Ensure we don't exceed array bounds
      n = min(this%n_fields, size(names))

      do i = 1, n
         names(i) = this%fields(i)%field_name
      end do
   end subroutine diag_registry_list_fields

   function diag_registry_get_count(this) result(count)
      class(DiagnosticRegistryType), intent(in) :: this
      integer :: count
      count = this%n_fields
   end function diag_registry_get_count

   function diag_registry_get_num_diagnostics(this) result(num)
      class(DiagnosticRegistryType), intent(in) :: this
      integer :: num
      num = this%n_fields
   end function diag_registry_get_num_diagnostics

   function diag_registry_field_exists(this, name) result(exists)
      class(DiagnosticRegistryType), intent(in) :: this
      character(len=*), intent(in) :: name
      logical :: exists
      integer :: i
      exists = .false.
      do i = 1, this%n_fields
         if (trim(this%fields(i)%field_name) == trim(name)) then
            exists = .true.
            return
         end if
      end do
   end function diag_registry_field_exists

   subroutine diag_registry_enable_field(this, name)
      class(DiagnosticRegistryType), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer :: i
      do i = 1, this%n_fields
         if (trim(this%fields(i)%field_name) == trim(name)) then
            call this%fields(i)%set_enabled(.true.)
            return
         end if
      end do
   end subroutine diag_registry_enable_field

   subroutine diag_registry_disable_field(this, name)
      class(DiagnosticRegistryType), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer :: i
      do i = 1, this%n_fields
         if (trim(this%fields(i)%field_name) == trim(name)) then
            call this%fields(i)%set_enabled(.false.)
            return
         end if
      end do
   end subroutine diag_registry_disable_field

   subroutine diag_registry_set_frequency(this, name, freq)
      class(DiagnosticRegistryType), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer, intent(in) :: freq
      integer :: i
      do i = 1, this%n_fields
         if (trim(this%fields(i)%field_name) == trim(name)) then
            this%fields(i)%output_frequency = freq
            return
         end if
      end do
   end subroutine diag_registry_set_frequency

   subroutine diag_registry_reset(this, rc)
      class(DiagnosticRegistryType), intent(inout) :: this
      integer, intent(out), optional :: rc
      integer :: i, local_rc
      if (present(rc)) rc = 0
      do i = 1, this%n_fields
         call this%fields(i)%reset_data(local_rc)
         if (present(rc) .and. local_rc /= 0) rc = local_rc
      end do
   end subroutine diag_registry_reset

   subroutine diag_registry_validate(this, error_mgr, rc)
      class(DiagnosticRegistryType), intent(in) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc
      integer :: i
      rc = CC_SUCCESS
      do i = 1, this%n_fields
         call this%fields(i)%validate_field(error_mgr, rc)
         if (rc /= CC_SUCCESS) return
      end do
   end subroutine diag_registry_validate

   !-------------------
   ! End DiagnosticRegistryType implementations
   !-------------------

end module DiagnosticInterface_Mod
