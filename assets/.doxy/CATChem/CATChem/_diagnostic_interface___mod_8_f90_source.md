

# File DiagnosticInterface\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**DiagnosticInterface\_Mod.F90**](_diagnostic_interface___mod_8_f90.md)

[Go to the documentation of this file](_diagnostic_interface___mod_8_f90.md)


```Fortran

module diagnosticinterface_mod
   use precision_mod, only: fp
   use error_mod, only: errormanagertype, cc_success, cc_failure, &
      error_memory_allocation, error_invalid_input, &
      error_duplicate_entry, error_not_found

   implicit none
   private

   integer, parameter :: max_fields = 100

   public :: diagnosticfieldtype, diagnosticregistrytype
   public :: diagnosticdatatype
   public :: diag_real_scalar, diag_real_1d, diag_real_2d, diag_real_3d
   public :: diag_integer_scalar, diag_integer_1d, diag_integer_2d, diag_integer_3d
   public :: diag_logical_scalar, diag_logical_1d, diag_logical_2d, diag_logical_3d
   public :: diag_freq_never, diag_freq_timestep, diag_freq_hourly, diag_freq_daily, diag_freq_custom

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
   integer, parameter :: DIAG_FREQ_NEVER = 0
   integer, parameter :: DIAG_FREQ_TIMESTEP = 1
   integer, parameter :: DIAG_FREQ_HOURLY = 2
   integer, parameter :: DIAG_FREQ_DAILY = 3
   integer, parameter :: DIAG_FREQ_CUSTOM = 99

   type :: diagnosticdatatype
      private
      integer :: data_type = 0
      logical :: is_allocated = .false.           

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
   end type diagnosticdatatype

   type :: diagnosticfieldtype
      private
      character(len=64) :: field_name = ''
      character(len=128) :: description = ''
      character(len=32) :: units = ''
      character(len=64) :: process_name = ''
      integer :: data_type = 0
      integer :: output_frequency = diag_freq_never 
      real(fp) :: custom_frequency = 0.0_fp        
      logical :: is_enabled = .true.               
      logical :: is_initialized = .false.          
      type(DiagnosticDataType) :: data

      ! Diagnostic species filtering support
      character(len=32), allocatable :: diagnostic_species(:)
      integer, allocatable :: diagnostic_species_id(:)

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
   end type diagnosticfieldtype

   type :: diagnosticregistrytype
      private
      character(len=64) :: process_name = ''
      integer :: n_fields = 0
      type(DiagnosticFieldType) :: fields(max_fields)
      logical :: is_initialized = .false.         

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
   end type diagnosticregistrytype

contains

   subroutine diag_data_allocate(this, data_type, dims, rc)
      class(DiagnosticDataType), intent(inout) :: this
      integer, intent(in) :: data_type
      integer, intent(in), optional :: dims(:)
      integer, intent(out) :: rc
      integer :: alloc_stat

      rc = cc_success
      alloc_stat = 0

      ! Deallocate any existing data first
      call this%deallocate_data()

      this%data_type = data_type

      select case (data_type)
       case (diag_real_scalar)
         this%real_scalar = 0.0_fp
         this%is_allocated = .true.

       case (diag_real_1d)
         if (.not. present(dims) .or. size(dims) < 1) then
            rc = error_invalid_input
            write(*,*) "Error: Invalid dimensions for DIAG_REAL_1D"
            return
         end if
         allocate(this%real_1d(dims(1)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = error_memory_allocation
            write(*,*) "Error: Memory allocation failed for DIAG_REAL_1D"
            return
         end if
         this%real_1d = 0.0_fp
         this%is_allocated = .true.

       case (diag_real_2d)
         if (.not. present(dims) .or. size(dims) < 2) then
            rc = error_invalid_input
            write(*,*) "Error: Invalid dimensions for DIAG_REAL_2D"
            return
         end if
         allocate(this%real_2d(dims(1), dims(2)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = error_memory_allocation
            write(*,*) "Error: Memory allocation failed for DIAG_REAL_2D"
            return
         end if
         this%real_2d = 0.0_fp
         this%is_allocated = .true.

       case (diag_real_3d)
         if (.not. present(dims) .or. size(dims) < 3) then
            rc = error_invalid_input
            write(*,*) "Error: Invalid dimensions for DIAG_REAL_3D"
            return
         end if
         allocate(this%real_3d(dims(1), dims(2), dims(3)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = error_memory_allocation
            write(*,*) "Error: Memory allocation failed for DIAG_REAL_3D"
            return
         end if
         this%real_3d = 0.0_fp
         this%is_allocated = .true.

       case (diag_integer_scalar)
         this%int_scalar = 0
         this%is_allocated = .true.

       case (diag_integer_1d)
         if (.not. present(dims) .or. size(dims) < 1) then
            rc = error_invalid_input
            write(*,*) "Error: Invalid dimensions for DIAG_INTEGER_1D"
            return
         end if
         allocate(this%int_1d(dims(1)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = error_memory_allocation
            write(*,*) "Error: Memory allocation failed for DIAG_INTEGER_1D"
            return
         end if
         this%int_1d = 0
         this%is_allocated = .true.

       case (diag_integer_2d)
         if (.not. present(dims) .or. size(dims) < 2) then
            rc = error_invalid_input
            write(*,*) "Error: Invalid dimensions for DIAG_INTEGER_2D"
            return
         end if
         allocate(this%int_2d(dims(1), dims(2)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = error_memory_allocation
            write(*,*) "Error: Memory allocation failed for DIAG_INTEGER_2D"
            return
         end if
         this%int_2d = 0
         this%is_allocated = .true.

       case (diag_integer_3d)
         if (.not. present(dims) .or. size(dims) < 3) then
            rc = error_invalid_input
            write(*,*) "Error: Invalid dimensions for DIAG_INTEGER_3D"
            return
         end if
         allocate(this%int_3d(dims(1), dims(2), dims(3)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = error_memory_allocation
            write(*,*) "Error: Memory allocation failed for DIAG_INTEGER_3D"
            return
         end if
         this%int_3d = 0
         this%is_allocated = .true.

       case (diag_logical_scalar)
         this%logical_scalar = .false.
         this%is_allocated = .true.

       case (diag_logical_1d)
         if (.not. present(dims) .or. size(dims) < 1) then
            rc = error_invalid_input
            write(*,*) "Error: Invalid dimensions for DIAG_LOGICAL_1D"
            return
         end if
         allocate(this%logical_1d(dims(1)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = error_memory_allocation
            write(*,*) "Error: Memory allocation failed for DIAG_LOGICAL_1D"
            return
         end if
         this%logical_1d = .false.
         this%is_allocated = .true.

       case (diag_logical_2d)
         if (.not. present(dims) .or. size(dims) < 2) then
            rc = error_invalid_input
            write(*,*) "Error: Invalid dimensions for DIAG_LOGICAL_2D"
            return
         end if
         allocate(this%logical_2d(dims(1), dims(2)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = error_memory_allocation
            write(*,*) "Error: Memory allocation failed for DIAG_LOGICAL_2D"
            return
         end if
         this%logical_2d = .false.
         this%is_allocated = .true.

       case (diag_logical_3d)
         if (.not. present(dims) .or. size(dims) < 3) then
            rc = error_invalid_input
            write(*,*) "Error: Invalid dimensions for DIAG_LOGICAL_3D"
            return
         end if
         allocate(this%logical_3d(dims(1), dims(2), dims(3)), stat=alloc_stat)
         if (alloc_stat /= 0) then
            rc = error_memory_allocation
            write(*,*) "Error: Memory allocation failed for DIAG_LOGICAL_3D"
            return
         end if
         this%logical_3d = .false.
         this%is_allocated = .true.

       case default
         rc = error_invalid_input
         write(*,*) "Error: Unsupported data type"
      end select

   end subroutine diag_data_allocate

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

   function diag_data_get_type(this) result(data_type)
      class(DiagnosticDataType), intent(in) :: this
      integer :: data_type
      data_type = this%data_type
   end function diag_data_get_type

   function diag_data_is_allocated(this) result(is_allocated)
      class(DiagnosticDataType), intent(in) :: this
      logical :: is_allocated
      is_allocated = this%is_allocated
   end function diag_data_is_allocated

   subroutine diag_data_set_real_scalar(this, value)
      class(DiagnosticDataType), intent(inout) :: this
      real(fp), intent(in) :: value
      if (this%data_type == diag_real_scalar) then
         this%real_scalar = value
      end if
   end subroutine diag_data_set_real_scalar

   function diag_data_get_real_scalar(this) result(value)
      class(DiagnosticDataType), intent(in) :: this
      real(fp) :: value
      if (this%data_type == diag_real_scalar) then
         value = this%real_scalar
      else
         value = 0.0_fp
      end if
   end function diag_data_get_real_scalar

   subroutine diag_data_set_real_1d(this, values)
      class(DiagnosticDataType), intent(inout) :: this
      real(fp), intent(in) :: values(:)
      if (this%data_type == diag_real_1d .and. allocated(this%real_1d)) then
         if (size(values) == size(this%real_1d)) then
            this%real_1d = values
         end if
      end if
   end subroutine diag_data_set_real_1d

   function diag_data_get_real_1d_ptr(this) result(ptr)
      class(DiagnosticDataType), intent(in), target :: this
      real(fp), pointer :: ptr(:)
      if (this%data_type == diag_real_1d .and. allocated(this%real_1d)) then
         ptr => this%real_1d
      else
         nullify(ptr)
      end if
   end function diag_data_get_real_1d_ptr

   subroutine diag_data_set_real_2d(this, values)
      class(DiagnosticDataType), intent(inout) :: this
      real(fp), intent(in) :: values(:,:)
      if (this%data_type == diag_real_2d .and. allocated(this%real_2d)) then
         if (size(values,1) == size(this%real_2d,1) .and. size(values,2) == size(this%real_2d,2)) then
            this%real_2d = values
         end if
      end if
   end subroutine diag_data_set_real_2d

   function diag_data_get_real_2d_ptr(this) result(ptr)
      class(DiagnosticDataType), intent(in), target :: this
      real(fp), pointer :: ptr(:,:)
      if (this%data_type == diag_real_2d .and. allocated(this%real_2d)) then
         ptr => this%real_2d
      else
         nullify(ptr)
      end if
   end function diag_data_get_real_2d_ptr

   subroutine diag_data_set_real_3d(this, values)
      class(DiagnosticDataType), intent(inout) :: this
      real(fp), intent(in) :: values(:,:,:)
      if (this%data_type == diag_real_3d .and. allocated(this%real_3d)) then
         if (size(values,1) == size(this%real_3d,1) .and. &
            size(values,2) == size(this%real_3d,2) .and. &
            size(values,3) == size(this%real_3d,3)) then
            this%real_3d = values
         end if
      end if
   end subroutine diag_data_set_real_3d

   function diag_data_get_real_3d_ptr(this) result(ptr)
      class(DiagnosticDataType), intent(in), target :: this
      real(fp), pointer :: ptr(:,:,:)
      if (this%data_type == diag_real_3d .and. allocated(this%real_3d)) then
         ptr => this%real_3d
      else
         nullify(ptr)
      end if
   end function diag_data_get_real_3d_ptr

   subroutine diag_data_finalize(this)
      type(DiagnosticDataType), intent(inout) :: this
      call this%deallocate_data()
   end subroutine diag_data_finalize

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

      rc = cc_success

      ! Validate inputs
      if (len_trim(field_name) == 0) then
         rc = error_invalid_input
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
      this%output_frequency = diag_freq_timestep
      this%is_enabled = .true.
      this%is_initialized = .true.

      ! Initialize data storage for all types with default dimensions for arrays
      select case (data_type)
       case (diag_real_scalar, diag_integer_scalar, diag_logical_scalar)
         ! For scalar types, we can initialize without dimensions
         call this%data%allocate_data(this%data_type, rc=rc)
         if (rc /= cc_success) then
            write(*,*) "Error: Could not allocate data for diagnostic field"
            this%is_initialized = .false.
            return
         end if
       case (diag_real_1d, diag_integer_1d, diag_logical_1d)
         ! For 1D arrays, initialize with default size of 1
         call this%data%allocate_data(this%data_type, [1], rc=rc)
         if (rc /= cc_success) then
            write(*,*) "Error: Could not allocate data for diagnostic field"
            this%is_initialized = .false.
            return
         end if
       case (diag_real_2d, diag_integer_2d, diag_logical_2d)
         ! For 2D arrays, initialize with default size of 1x1
         call this%data%allocate_data(this%data_type, [1, 1], rc=rc)
         if (rc /= cc_success) then
            write(*,*) "Error: Could not allocate data for diagnostic field"
            this%is_initialized = .false.
            return
         end if
       case (diag_real_3d, diag_integer_3d, diag_logical_3d)
         ! For 3D arrays, initialize with default size of 1x1x1
         call this%data%allocate_data(this%data_type, [1, 1, 1], rc=rc)
         if (rc /= cc_success) then
            write(*,*) "Error: Could not allocate data for diagnostic field"
            this%is_initialized = .false.
            return
         end if
       case default
         rc = error_invalid_input
         this%is_initialized = .false.
         return
      end select

      ! Verify that data was allocated successfully
      if (.not. this%data%is_data_allocated()) then
         rc = error_memory_allocation
         write(*,*) "Error: Data allocation failed during create()"
         this%is_initialized = .false.
         return
      end if

   end subroutine diag_field_create

   subroutine diag_field_init_data(this, dims, rc)
      class(DiagnosticFieldType), intent(inout) :: this
      integer, intent(in), optional :: dims(:)
      integer, intent(out) :: rc

      rc = cc_success

      if (.not. this%is_initialized) then
         rc = error_invalid_input
         return
      end if

      call this%data%allocate_data(this%data_type, dims, rc)

   end subroutine diag_field_init_data

   subroutine diag_field_cleanup(this)
      class(DiagnosticFieldType), intent(inout) :: this
      call this%data%deallocate_data()
      if (allocated(this%diagnostic_species)) deallocate(this%diagnostic_species)
      if (allocated(this%diagnostic_species_id)) deallocate(this%diagnostic_species_id)
      this%is_initialized = .false.
   end subroutine diag_field_cleanup

   function diag_field_is_ready(this) result(is_ready)
      class(DiagnosticFieldType), intent(in) :: this
      logical :: is_ready
      is_ready = this%is_initialized .and. this%data%is_data_allocated()
   end function diag_field_is_ready

   function diag_field_get_name(this) result(name)
      class(DiagnosticFieldType), intent(in) :: this
      character(len=64) :: name
      name = this%field_name
   end function diag_field_get_name

   function diag_field_get_description(this) result(description)
      class(DiagnosticFieldType), intent(in) :: this
      character(len=128) :: description
      description = this%description
   end function diag_field_get_description

   function diag_field_get_units(this) result(units)
      class(DiagnosticFieldType), intent(in) :: this
      character(len=32) :: units
      units = this%units
   end function diag_field_get_units

   function diag_field_get_process_name(this) result(process_name)
      class(DiagnosticFieldType), intent(in) :: this
      character(len=64) :: process_name
      process_name = this%process_name
   end function diag_field_get_process_name

   function diag_field_get_data_type(this) result(data_type)
      class(DiagnosticFieldType), intent(in) :: this
      integer :: data_type
      data_type = this%data_type
   end function diag_field_get_data_type

   function diag_field_get_data_ptr(this) result(data_ptr)
      class(DiagnosticFieldType), intent(in), target :: this
      type(DiagnosticDataType), pointer :: data_ptr
      data_ptr => this%data
   end function diag_field_get_data_ptr

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

   subroutine diag_field_set_enabled(this, enabled)
      class(DiagnosticFieldType), intent(inout) :: this
      logical, intent(in) :: enabled
      this%is_enabled = enabled
   end subroutine diag_field_set_enabled

   function diag_field_is_enabled(this) result(is_enabled)
      class(DiagnosticFieldType), intent(in) :: this
      logical :: is_enabled
      is_enabled = this%is_enabled
   end function diag_field_is_enabled

   function diag_field_should_output(this, current_time, dt) result(should_output)
      class(DiagnosticFieldType), intent(in) :: this
      real(fp), intent(in) :: current_time
      real(fp), intent(in) :: dt
      logical :: should_output

      should_output = .false.

      if (.not. this%is_enabled) return

      select case (this%output_frequency)
       case (diag_freq_never)
         should_output = .false.
       case (diag_freq_timestep)
         should_output = .true.
       case (diag_freq_hourly)
         should_output = (mod(current_time, 3600.0_fp) < dt)
       case (diag_freq_daily)
         should_output = (mod(current_time, 86400.0_fp) < dt)
       case (diag_freq_custom)
         if (this%custom_frequency > 0.0_fp) then
            should_output = (mod(current_time, this%custom_frequency) < dt)
         end if
      end select

   end function diag_field_should_output

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

   subroutine diag_field_finalize(this)
      type(DiagnosticFieldType), intent(inout) :: this
      call this%cleanup()
   end subroutine diag_field_finalize

   subroutine diag_field_reset_data(this, rc)
      class(DiagnosticFieldType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success

      if (.not. this%is_ready()) then
         rc = cc_failure
         return
      endif

      ! Reset data based on type
      select case (this%data%get_data_type())
       case (diag_real_scalar)
         call this%data%set_real_scalar(0.0_fp)
       case (diag_real_1d)
         ! Would set 1D array to zero
       case (diag_real_2d)
         ! Would set 2D array to zero
       case (diag_real_3d)
         ! Would set 3D array to zero
      end select

   end subroutine diag_field_reset_data

   subroutine diag_field_validate_field(this, error_mgr, rc)
      class(DiagnosticFieldType), intent(in) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = cc_success

      ! Validate field name
      if (len_trim(this%field_name) == 0) then
         call error_mgr%report_error(error_invalid_input, 'Diagnostic field name is empty', rc)
         return
      endif

      ! Validate data type
      if (this%data_type <= 0) then
         call error_mgr%report_error(error_invalid_input, 'Invalid diagnostic data type', rc)
         return
      endif

      ! Validate initialization
      if (.not. this%is_initialized) then
         call error_mgr%report_error(error_invalid_input, 'Diagnostic field not initialized', rc)
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
      rc = cc_success
      if (.not. this%is_initialized) then
         rc = error_invalid_input
         return
      end if
      ! Check if field is valid
      if (.not. field%is_ready()) then
         rc = error_invalid_input
         return
      end if
      ! Check for duplicate
      do i = 1, this%n_fields
         if (trim(this%fields(i)%field_name) == trim(field%field_name)) then
            rc = error_duplicate_entry
            return
         end if
      end do
      if (this%n_fields >= max_fields) then
         rc = error_memory_allocation
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
      field = diagnosticfieldtype()
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
      rc = cc_success
      do i = 1, this%n_fields
         call this%fields(i)%validate_field(error_mgr, rc)
         if (rc /= cc_success) return
      end do
   end subroutine diag_registry_validate

   !-------------------
   ! End DiagnosticRegistryType implementations
   !-------------------

end module diagnosticinterface_mod
```


