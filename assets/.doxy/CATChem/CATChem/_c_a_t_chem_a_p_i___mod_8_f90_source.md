

# File CATChemAPI\_Mod.F90

[**File List**](files.md) **>** [**api**](dir_da61e3e9a357748887e3ca8d7c5a0c16.md) **>** [**CATChemAPI\_Mod.F90**](_c_a_t_chem_a_p_i___mod_8_f90.md)

[Go to the documentation of this file](_c_a_t_chem_a_p_i___mod_8_f90.md)


```Fortran

module catchemapi_mod
   use precision_mod
   use error_mod
   use state_mod, only: statecontainertype, statebuildertype
   use configmanager_mod, only: configmanagertype
   use processmanager_mod, only: processmanagertype
   use gridmanager_mod, only: gridmanagertype
   use columninterface_mod, only: virtualcolumntype
   use chemstate_mod, only: chemstatetype
   use metstate_mod, only: metstatetype
   use emisstate_mod, only: emisstatetype
   use diagstate_mod, only: diagstatetype
   use diagnosticinterface_mod, only: diagnosticmanagertype, diagnosticregistrytype, &
      diagnosticfieldtype

   implicit none
   private

   public :: catcheminstancetype
   public :: catchemconfigtype
   public :: catchemdatatype
   public :: catchemdiagnostictype
   public :: catchem_success, catchem_failure

   ! Return codes
   integer, parameter :: CATCHEM_SUCCESS = 0
   integer, parameter :: CATCHEM_FAILURE = -1

   type :: catchemconfigtype
      ! Grid dimensions
      integer :: nx = 0
      integer :: ny = 0
      integer :: nz = 0
      ! NOTE: nspecies is now read from configuration files, not set directly

      ! Time parameters
      real(fp) :: dt = 3600.0_fp           
      integer :: nsteps = 1

      ! Process configuration
      logical :: enable_dust = .false.      
      logical :: enable_seasalt = .false.   
      logical :: enable_drydep = .false.    
      logical :: enable_chemistry = .false. 
      logical :: enable_external_emis = .false. 

      ! Run phase configuration
      logical :: enable_run_phases = .false.   
      character(len=64), allocatable :: run_phase_names(:)

      ! Performance options
      logical :: use_column_processing = .true.  
      logical :: enable_diagnostics = .true.     

      ! Configuration files - these define the initialization order
      character(len=256) :: config_file = ''
      character(len=256) :: species_file = ''
      character(len=256) :: emission_file = ''
   end type catchemconfigtype

   type :: fieldmappingtype
      character(len=64) :: host_field_name = ''
      character(len=64) :: catchem_field_name = ''
      character(len=32) :: field_type = ''
      character(len=32) :: units = ''
      real(fp) :: scale_factor = 1.0_fp            
      real(fp) :: offset = 0.0_fp                  
      logical :: is_required = .true.              
      logical :: is_3d = .true.                    
   end type fieldmappingtype

   type :: fieldmappingregistrytype
      type(FieldMappingType), allocatable :: met_mappings(:)
      type(FieldMappingType), allocatable :: chem_mappings(:)
      type(FieldMappingType), allocatable :: emis_mappings(:)
      type(FieldMappingType), allocatable :: diag_mappings(:)
      integer :: num_met_mappings = 0
      integer :: num_chem_mappings = 0
      integer :: num_emis_mappings = 0
      integer :: num_diag_mappings = 0
   contains
      procedure :: initialize => registry_initialize
      procedure :: add_met_mapping => registry_add_met_mapping
      procedure :: add_chem_mapping => registry_add_chem_mapping
      procedure :: add_emis_mapping => registry_add_emis_mapping
      procedure :: add_diag_mapping => registry_add_diag_mapping
   end type fieldmappingregistrytype

   type :: flexibledataexchangetype
      ! Field registry
      type(FieldMappingRegistryType) :: field_mapping_registry

      ! Dynamic data storage
      real(fp), allocatable :: met_data_3d(:,:,:,:)
      real(fp), allocatable :: met_data_2d(:,:,:)
      real(fp), allocatable :: chem_data(:,:,:,:)
      real(fp), allocatable :: emis_data_3d(:,:,:,:)
      real(fp), allocatable :: emis_data_2d(:,:,:)
      real(fp), allocatable :: diag_data_3d(:,:,:,:)
      real(fp), allocatable :: diag_data_2d(:,:,:)

      ! Field name arrays for dynamic mapping
      character(len=64), allocatable :: met_field_names(:)
      character(len=64), allocatable :: chem_field_names(:)
      character(len=64), allocatable :: emis_field_names(:)
      character(len=64), allocatable :: diag_field_names(:)

      ! Dimensions
      integer :: nx, ny, nz
      integer :: n_met_3d = 0, n_met_2d = 0
      integer :: n_chem = 0, n_emis_3d = 0, n_emis_2d = 0
      integer :: n_diag_3d = 0, n_diag_2d = 0

   contains
      procedure :: setup_flexible_exchange => flex_setup_exchange
      procedure :: map_host_field => flex_map_host_field
      procedure :: get_mapped_field => flex_get_mapped_field
      procedure :: set_mapped_field => flex_set_mapped_field
      procedure :: validate_mappings => flex_validate_mappings
   end type flexibledataexchangetype

   type :: catchemdatatype
      ! Meteorological inputs
      real(fp), allocatable :: temperature(:,:,:)
      real(fp), allocatable :: pressure(:,:,:)
      real(fp), allocatable :: humidity(:,:,:)
      real(fp), allocatable :: wind_u(:,:,:)
      real(fp), allocatable :: wind_v(:,:,:)
      real(fp), allocatable :: surface_pressure(:,:)
      real(fp), allocatable :: roughness_length(:,:)

      ! Chemical species concentrations
      real(fp), allocatable :: concentrations(:,:,:,:)

      ! Emission data (optional)
      real(fp), allocatable :: emission_rates(:,:,:,:)

      ! Diagnostic outputs
      real(fp), allocatable :: dust_emissions(:,:,:)
      real(fp), allocatable :: seasalt_emissions(:,:,:)
      real(fp), allocatable :: drydep_velocity(:,:,:)
   end type catchemdatatype

   type :: catchemdiagnostictype
      character(len=64) :: field_name = ''
      character(len=128) :: description = ''
      character(len=32) :: units = ''
      character(len=32) :: process_name = ''
      integer :: dimensions = 0
      logical :: is_available = .false.           
      real(fp), allocatable :: data_2d(:,:)
      real(fp), allocatable :: data_3d(:,:,:)
      real(fp) :: scalar_data = 0.0_fp            
   end type catchemdiagnostictype

   type :: catcheminstancetype
      private

      ! Core components
      type(StateContainerType) :: container
      type(ProcessManagerType) :: process_manager
      type(ConfigManagerType)  :: config_manager
      type(ErrorManagerType), pointer :: error_mgr => null()

      ! Configuration
      type(CATChemConfigType) :: config
      logical :: is_initialized = .false.
      logical :: is_ready = .false.

      ! Grid and species information (read from config files)
      integer :: nx, ny, nz, nspecies

      ! Run phase support
      logical :: multi_phase_mode = .false.
      character(len=64), allocatable :: active_phases(:)
      integer :: current_phase = 0

   contains
      ! Simplified interface methods
      procedure :: init => instance_init
      procedure :: init_from_config_files => instance_init_from_config_files
      procedure :: setup_grid => instance_setup_grid
      procedure :: add_process => instance_add_process
      procedure :: run_timestep => instance_run_timestep
      procedure :: run_multiple_steps => instance_run_multiple_steps
      procedure :: get_concentrations => instance_get_concentrations
      procedure :: set_concentrations => instance_set_concentrations
      procedure :: set_meteorology => instance_set_meteorology
      procedure :: set_emissions => instance_set_emissions
      procedure :: get_diagnostics => instance_get_diagnostics
      procedure :: finalize => instance_finalize

      ! Run phase management
      procedure :: setup_run_phases => instance_setup_run_phases
      procedure :: run_phase => instance_run_phase
      procedure :: run_all_phases => instance_run_all_phases
      procedure :: get_phase_names => instance_get_phase_names
      procedure :: get_current_phase => instance_get_current_phase

      ! Process manager configuration
      procedure :: configure_process_manager => instance_configure_process_manager
      procedure :: set_process_manager => instance_set_process_manager
      ! Status and utility methods
      procedure :: is_ready_to_run => instance_is_ready
      procedure :: get_species_names => instance_get_species_names
      procedure :: get_error_message => instance_get_error_message
      procedure :: validate_data => instance_validate_data

      ! Process-level diagnostic methods
      procedure :: get_available_diagnostics => instance_get_available_diagnostics
      procedure :: get_process_diagnostic => instance_get_process_diagnostic
      procedure :: get_all_process_diagnostics => instance_get_all_process_diagnostics
      procedure :: list_process_diagnostics => instance_list_process_diagnostics

      ! High-level field mapping methods
      procedure :: fill_met_state_from_host => instance_fill_met_state_from_host
      procedure :: extract_met_state_to_host => instance_extract_met_state_to_host
      procedure :: fill_chem_state_from_host => instance_fill_chem_state_from_host
      procedure :: extract_chem_state_to_host => instance_extract_chem_state_to_host
      procedure :: fill_emis_state_from_host => instance_fill_emis_state_from_host
      procedure :: extract_emis_state_to_host => instance_extract_emis_state_to_host
      procedure :: extract_diagnostics_to_host => instance_extract_diagnostics_to_host

      ! Field mapping registry management
      procedure :: setup_field_mapping => instance_setup_field_mapping
      procedure :: add_field_mapping => instance_add_field_mapping
      procedure :: validate_field_mappings => instance_validate_field_mappings

      ! Advanced access (for power users)
      procedure :: get_state_container => instance_get_container
      procedure :: get_process_manager => instance_get_process_manager
      procedure :: get_diagnostic_manager => instance_get_diagnostic_manager
   end type catcheminstancetype

contains

   subroutine instance_init(this, config, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(CATChemConfigType), intent(in) :: config
      integer, intent(out) :: rc

      type(StateBuilderType) :: builder
      character(len=256) :: error_msg

      rc = catchem_success
      this%config = config

      ! Validate basic configuration
      if (config%nx <= 0 .or. config%ny <= 0 .or. config%nz <= 0) then
         rc = catchem_failure
         return
      endif

      ! Store grid dimensions
      this%nx = config%nx
      this%ny = config%ny
      this%nz = config%nz

      ! Initialize configuration manager
      call this%config_manager%init(rc)
      if (rc /= catchem_success) return

      ! Load configuration file if provided
      if (len_trim(config%config_file) > 0) then
         call this%config_manager%load_config(config%config_file, rc)
         if (rc /= catchem_success) return
      endif

      ! Build state container with sensible defaults
      call builder%init()
      call builder%set_grid_dimensions(config%nx, config%ny, config%nz)
      call builder%set_num_species(config%nspecies)
      call builder%enable_column_processing(config%use_column_processing)
      call builder%enable_diagnostics(config%enable_diagnostics)

      ! Build the container
      call builder%build(this%container, rc)
      if (rc /= catchem_success) return

      ! Get error manager for future use
      this%error_mgr => this%container%get_error_manager()

      ! Initialize process manager
      call this%process_manager%init(this%container, rc)
      if (rc /= catchem_success) return

      this%is_initialized = .true.

   end subroutine instance_init

   subroutine instance_init_from_config_files(this, config, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(CATChemConfigType), intent(in) :: config
      integer, intent(out) :: rc

      type(StateBuilderType) :: builder
      integer :: temp_nspecies

      rc = catchem_success
      this%config = config

      ! Validate that all required config files are provided
      if (len_trim(config%config_file) == 0) then
         rc = catchem_failure
         return
      endif

      ! Initialize configuration manager
      call this%config_manager%init(rc)
      if (rc /= catchem_success) return

      ! Step 1: Read main configuration file (sets grid, processes, etc.)
      call this%config_manager%load_config(config%config_file, rc)
      if (rc /= catchem_success) return

      ! Extract grid dimensions from main config
      call this%config_manager%get_grid_dimensions(this%nx, this%ny, this%nz, rc)
      if (rc /= catchem_success) return

      ! Step 2: Read species configuration file (determines nspecies)
      if (len_trim(config%species_file) > 0) then
         call this%config_manager%load_species_data(config%species_file, rc)
         if (rc /= catchem_success) return
      else
         ! Try to get species file from main config
         call this%config_manager%get_species_filename(config%species_file, rc)
         if (rc == catchem_success .and. len_trim(config%species_file) > 0) then
            call this%config_manager%load_species_data(config%species_file, rc)
            if (rc /= catchem_success) return
         endif
      endif

      ! Get number of species from loaded configuration
      call this%config_manager%get_num_species(this%nspecies, rc)
      if (rc /= catchem_success) return

      ! Step 3: Read emission configuration file (if provided)
      if (len_trim(config%emission_file) > 0) then
         call this%config_manager%load_emission_config(config%emission_file, rc)
         if (rc /= catchem_success) return
      endif

      ! Build state container with configuration-driven parameters
      call builder%init()
      call builder%set_grid_dimensions(this%nx, this%ny, this%nz)
      call builder%set_num_species(this%nspecies)
      call builder%enable_column_processing(config%use_column_processing)
      call builder%enable_diagnostics(config%enable_diagnostics)

      ! Build the container
      call builder%build(this%container, rc)
      if (rc /= catchem_success) return

      ! Get error manager for future use
      this%error_mgr => this%container%get_error_manager()

      ! Initialize process manager
      call this%process_manager%init(this%container, rc)
      if (rc /= catchem_success) return

      ! Setup run phases if enabled
      if (config%enable_run_phases) then
         call this%setup_run_phases(config%run_phase_names, rc)
         if (rc /= catchem_success) return
      endif

      this%is_initialized = .true.

   end subroutine instance_init_from_config_files

   subroutine instance_setup_grid(this, lats, lons, levels, rc)
      class(CATChemInstanceType), intent(inout) :: this
      real(fp), intent(in) :: lats(:,:)
      real(fp), intent(in) :: lons(:,:)
      real(fp), intent(in) :: levels(:)
      integer, intent(out) :: rc

      type(GridManagerType), pointer :: grid_mgr

      rc = catchem_success

      if (.not. this%is_initialized) then
         rc = catchem_failure
         return
      endif

      ! Get grid manager and setup coordinates
      grid_mgr => this%container%get_grid_manager_ptr()
      if (.not. associated(grid_mgr)) then
         rc = catchem_failure
         return
      endif

      ! Setup coordinate system
      call grid_mgr%set_coordinates(lats, lons, levels, rc)
      if (rc /= catchem_success) return

      ! Enable column virtualization if requested
      if (this%config%use_column_processing) then
         call grid_mgr%enable_virtualization(rc)
         if (rc /= catchem_success) return
      endif

   end subroutine instance_setup_grid

   subroutine instance_add_process(this, process_name, process_config, rc)
      class(CATChemInstanceType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      character(len=*), intent(in), optional :: process_config
      integer, intent(out) :: rc

      rc = catchem_success

      if (.not. this%is_initialized) then
         rc = catchem_failure
         return
      endif

      ! Add process to manager
      if (present(process_config)) then
         call this%process_manager%add_process(process_name, process_config, rc)
      else
         call this%process_manager%add_process(process_name, rc=rc)
      endif

      if (rc == catchem_success) then
         this%is_ready = .true.
      endif

   end subroutine instance_add_process

   subroutine instance_setup_run_phases(this, phase_names, rc)
      class(CATChemInstanceType), intent(inout) :: this
      character(len=64), intent(in) :: phase_names(:)
      integer, intent(out) :: rc

      integer :: i, num_phases

      rc = catchem_success

      if (.not. this%is_initialized) then
         rc = catchem_failure
         return
      endif

      num_phases = size(phase_names)
      if (num_phases <= 0) then
         rc = catchem_failure
         return
      endif

      ! Store phase names
      if (allocated(this%active_phases)) deallocate(this%active_phases)
      allocate(this%active_phases(num_phases))
      this%active_phases = phase_names

      ! Configure process manager for multi-phase execution
      call this%process_manager%configure_run_phases(phase_names, rc)
      if (rc /= catchem_success) return

      this%multi_phase_mode = .true.
      this%current_phase = 0

   end subroutine instance_setup_run_phases

   subroutine instance_run_phase(this, phase_name, rc)
      class(CATChemInstanceType), intent(inout) :: this
      character(len=*), intent(in) :: phase_name
      integer, intent(out) :: rc

      integer :: phase_index

      rc = catchem_success

      if (.not. this%is_ready_to_run()) then
         rc = catchem_failure
         return
      endif

      if (.not. this%multi_phase_mode) then
         rc = catchem_failure
         return
      endif

      ! Find phase index
      phase_index = 0
      do i = 1, size(this%active_phases)
         if (trim(this%active_phases(i)) == trim(phase_name)) then
            phase_index = i
            exit
         endif
      enddo

      if (phase_index == 0) then
         rc = catchem_failure
         return
      endif

      ! Run the specific phase
      call this%process_manager%run_phase(phase_name, this%container, rc)
      if (rc /= catchem_success) return

      this%current_phase = phase_index

   end subroutine instance_run_phase

   subroutine instance_run_all_phases(this, rc)
      class(CATChemInstanceType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i

      rc = catchem_success

      if (.not. this%multi_phase_mode) then
         ! Fall back to regular run_all if phases not configured
         call this%process_manager%run_all_processes(this%container, rc)
         return
      endif

      ! Run each phase in sequence
      do i = 1, size(this%active_phases)
         call this%instance_run_phase(this%active_phases(i), rc)
         if (rc /= catchem_success) return
      enddo

   end subroutine instance_run_all_phases

   subroutine instance_get_phase_names(this, phase_names, rc)
      class(CATChemInstanceType), intent(in) :: this
      character(len=64), allocatable, intent(out) :: phase_names(:)
      integer, intent(out) :: rc

      rc = catchem_success

      if (.not. this%multi_phase_mode) then
         rc = catchem_failure
         return
      endif

      if (allocated(this%active_phases)) then
         allocate(phase_names(size(this%active_phases)))
         phase_names = this%active_phases
      else
         rc = catchem_failure
      endif

   end subroutine instance_get_phase_names

   function instance_get_current_phase(this) result(phase_index)
      class(CATChemInstanceType), intent(in) :: this
      integer :: phase_index

      phase_index = this%current_phase

   end function instance_get_current_phase

   subroutine instance_configure_process_manager(this, max_processes, enable_column_batching, rc)
      class(CATChemInstanceType), intent(inout) :: this
      integer, intent(in), optional :: max_processes
      logical, intent(in), optional :: enable_column_batching
      integer, intent(out) :: rc

      rc = catchem_success

      if (.not. this%is_initialized) then
         rc = catchem_failure
         return
      endif

      if (present(max_processes)) then
         call this%process_manager%set_max_processes(max_processes, rc)
         if (rc /= catchem_success) return
      endif

      if (present(enable_column_batching)) then
         call this%process_manager%enable_column_batching(enable_column_batching, rc)
         if (rc /= catchem_success) return
      endif

   end subroutine instance_configure_process_manager

   subroutine instance_set_process_manager(this, process_manager, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(ProcessManagerType), intent(in) :: process_manager
      integer, intent(out) :: rc

      rc = catchem_success

      if (.not. this%is_initialized) then
         rc = catchem_failure
         return
      endif

      this%process_manager = process_manager

   end subroutine instance_set_process_manager

   subroutine instance_run_timestep(this, input_data, output_data, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(CATChemDataType), intent(in) :: input_data
      type(CATChemDataType), intent(out) :: output_data
      integer, intent(out) :: rc

      rc = catchem_success

      if (.not. this%is_ready_to_run()) then
         rc = catchem_failure
         return
      endif

      ! Set input data
      call this%set_meteorology(input_data, rc)
      if (rc /= catchem_success) return

      call this%set_concentrations(input_data, rc)
      if (rc /= catchem_success) return

      if (allocated(input_data%emission_rates)) then
         call this%set_emissions(input_data, rc)
         if (rc /= catchem_success) return
      endif

      ! Run processes - use phase-based execution if configured
      if (this%multi_phase_mode) then
         call this%run_all_phases(rc)
      else
         call this%process_manager%run_all_processes(this%container, rc)
      endif
      if (rc /= catchem_success) return

      ! Get output data
      call this%get_concentrations(output_data, rc)
      if (rc /= catchem_success) return

      call this%get_diagnostics(output_data, rc)
      if (rc /= catchem_success) return

   end subroutine instance_run_timestep

   subroutine instance_run_multiple_steps(this, nsteps, input_data, output_data, rc)
      class(CATChemInstanceType), intent(inout) :: this
      integer, intent(in) :: nsteps
      type(CATChemDataType), intent(in) :: input_data(:)
      type(CATChemDataType), intent(out) :: output_data(:)
      integer, intent(out) :: rc

      integer :: i

      rc = catchem_success

      if (size(input_data) /= nsteps .or. size(output_data) /= nsteps) then
         rc = catchem_failure
         return
      endif

      do i = 1, nsteps
         call this%run_timestep(input_data(i), output_data(i), rc)
         if (rc /= catchem_success) return
      enddo

   end subroutine instance_run_multiple_steps

   function instance_is_ready(this) result(is_ready)
      class(CATChemInstanceType), intent(in) :: this
      logical :: is_ready

      is_ready = this%is_initialized .and. this%is_ready

   end function instance_is_ready

   subroutine instance_set_meteorology(this, data, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(CATChemDataType), intent(in) :: data
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state

      rc = catchem_success

      met_state => this%container%get_met_state_ptr()
      if (.not. associated(met_state)) then
         rc = catchem_failure
         return
      endif

      ! Copy meteorological data to internal state
      if (allocated(data%temperature)) then
         met_state%temperature = data%temperature
      endif

      if (allocated(data%pressure)) then
         met_state%pressure = data%pressure
      endif

      if (allocated(data%humidity)) then
         met_state%humidity = data%humidity
      endif

      if (allocated(data%wind_u) .and. allocated(data%wind_v)) then
         met_state%wind_u = data%wind_u
         met_state%wind_v = data%wind_v
      endif

   end subroutine instance_set_meteorology

   subroutine instance_set_concentrations(this, data, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(CATChemDataType), intent(in) :: data
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state

      rc = catchem_success

      chem_state => this%container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         rc = catchem_failure
         return
      endif

      if (allocated(data%concentrations)) then
         chem_state%concentrations = data%concentrations
      endif

   end subroutine instance_set_concentrations

   subroutine instance_set_emissions(this, data, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(CATChemDataType), intent(in) :: data
      integer, intent(out) :: rc

      type(EmisStateType), pointer :: emis_state

      rc = catchem_success

      emis_state => this%container%get_emis_state_ptr()
      if (.not. associated(emis_state)) then
         rc = catchem_failure
         return
      endif

      if (allocated(data%emission_rates)) then
         emis_state%emission_rates = data%emission_rates
      endif

   end subroutine instance_set_emissions

   subroutine instance_get_concentrations(this, data, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(CATChemDataType), intent(out) :: data
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state

      rc = catchem_success

      chem_state => this%container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         rc = catchem_failure
         return
      endif

      if (allocated(chem_state%concentrations)) then
         if (.not. allocated(data%concentrations)) then
            allocate(data%concentrations(size(chem_state%concentrations,1), &
               size(chem_state%concentrations,2), &
               size(chem_state%concentrations,3), &
               size(chem_state%concentrations,4)))
         endif
         data%concentrations = chem_state%concentrations
      endif

   end subroutine instance_get_concentrations

   subroutine instance_get_diagnostics(this, data, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(CATChemDataType), intent(out) :: data
      integer, intent(out) :: rc

      type(DiagStateType), pointer :: diag_state

      rc = catchem_success

      diag_state => this%container%get_diag_state_ptr()
      if (.not. associated(diag_state)) then
         rc = catchem_failure
         return
      endif

      ! Extract diagnostic fields if available
      ! TODO: Add specific diagnostic field extraction based on available processes

   end subroutine instance_get_diagnostics

   subroutine instance_get_error_message(this, error_msg)
      class(CATChemInstanceType), intent(in) :: this
      character(len=*), intent(out) :: error_msg

      if (associated(this%error_mgr)) then
         call this%error_mgr%get_error_message(error_msg)
      else
         error_msg = 'No error manager available'
      endif

   end subroutine instance_get_error_message

   subroutine instance_validate_data(this, data, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(CATChemDataType), intent(in) :: data
      integer, intent(out) :: rc

      rc = catchem_success

      ! Basic validation - check array dimensions match expected grid
      if (allocated(data%temperature)) then
         if (size(data%temperature,1) /= this%nx .or. &
            size(data%temperature,2) /= this%ny .or. &
            size(data%temperature,3) /= this%nz) then
            rc = catchem_failure
            return
         endif
      endif

      if (allocated(data%concentrations)) then
         if (size(data%concentrations,1) /= this%nx .or. &
            size(data%concentrations,2) /= this%ny .or. &
            size(data%concentrations,3) /= this%nz .or. &
            size(data%concentrations,4) /= this%nspecies) then
            rc = catchem_failure
            return
         endif
      endif

   end subroutine instance_validate_data

   subroutine instance_get_species_names(this, species_names, rc)
      class(CATChemInstanceType), intent(in) :: this
      character(len=64), allocatable, intent(out) :: species_names(:)
      integer, intent(out) :: rc

      rc = catchem_success

      ! TODO: Extract species names from configuration
      allocate(species_names(this%nspecies))
      ! For now, use generic names
      ! In real implementation, this would come from the species configuration file

   end subroutine instance_get_species_names

   subroutine instance_finalize(this, rc)
      class(CATChemInstanceType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = catchem_success

      call this%process_manager%finalize(rc)
      call this%container%finalize(rc)

      this%is_initialized = .false.
      this%is_ready = .false.

   end subroutine instance_finalize

   ! Placeholder implementations for missing methods
   ! These would be implemented based on the actual StateContainer and ProcessManager interfaces

   subroutine instance_get_available_diagnostics(this, diagnostics, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(CATChemDiagnosticType), allocatable, intent(out) :: diagnostics(:)
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Implement diagnostic retrieval
   end subroutine instance_get_available_diagnostics

   subroutine instance_get_process_diagnostic(this, process_name, diag_name, diagnostic, rc)
      class(CATChemInstanceType), intent(in) :: this
      character(len=*), intent(in) :: process_name, diag_name
      type(CATChemDiagnosticType), intent(out) :: diagnostic
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Implement specific diagnostic retrieval
   end subroutine instance_get_process_diagnostic

   subroutine instance_get_all_process_diagnostics(this, diagnostics, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(CATChemDiagnosticType), allocatable, intent(out) :: diagnostics(:)
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Implement all diagnostics retrieval
   end subroutine instance_get_all_process_diagnostics

   subroutine instance_list_process_diagnostics(this, process_names, diagnostic_info, rc)
      class(CATChemInstanceType), intent(in) :: this
      character(len=64), allocatable, intent(out) :: process_names(:)
      character(len=256), allocatable, intent(out) :: diagnostic_info(:)
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Implement diagnostic listing
   end subroutine instance_list_process_diagnostics

   ! Field mapping methods (placeholder implementations)
   subroutine instance_fill_met_state_from_host(this, data_exchange, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(FlexibleDataExchangeType), intent(in) :: data_exchange
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Implement field mapping from host
   end subroutine instance_fill_met_state_from_host

   subroutine instance_extract_met_state_to_host(this, data_exchange, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(FlexibleDataExchangeType), intent(out) :: data_exchange
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Implement field mapping to host
   end subroutine instance_extract_met_state_to_host

   subroutine instance_fill_chem_state_from_host(this, data_exchange, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(FlexibleDataExchangeType), intent(in) :: data_exchange
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Implement chemical state mapping from host
   end subroutine instance_fill_chem_state_from_host

   subroutine instance_extract_chem_state_to_host(this, data_exchange, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(FlexibleDataExchangeType), intent(out) :: data_exchange
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Implement chemical state mapping to host
   end subroutine instance_extract_chem_state_to_host

   subroutine instance_fill_emis_state_from_host(this, data_exchange, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(FlexibleDataExchangeType), intent(in) :: data_exchange
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Implement emission state mapping from host
   end subroutine instance_fill_emis_state_from_host

   subroutine instance_extract_emis_state_to_host(this, data_exchange, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(FlexibleDataExchangeType), intent(out) :: data_exchange
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Implement emission state mapping to host
   end subroutine instance_extract_emis_state_to_host

   subroutine instance_extract_diagnostics_to_host(this, data_exchange, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(FlexibleDataExchangeType), intent(out) :: data_exchange
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Implement diagnostic mapping to host
   end subroutine instance_extract_diagnostics_to_host

   subroutine instance_setup_field_mapping(this, mappings, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(FieldMappingType), intent(in) :: mappings(:)
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Implement field mapping setup
   end subroutine instance_setup_field_mapping

   subroutine instance_add_field_mapping(this, field_type, catchem_name, host_name, data_type, rc)
      class(CATChemInstanceType), intent(inout) :: this
      character(len=*), intent(in) :: field_type, catchem_name, host_name, data_type
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Implement individual field mapping addition
   end subroutine instance_add_field_mapping

   subroutine instance_validate_field_mappings(this, rc)
      class(CATChemInstanceType), intent(in) :: this
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Implement field mapping validation
   end subroutine instance_validate_field_mappings

   subroutine instance_get_container(this, container)
      class(CATChemInstanceType), intent(in) :: this
      type(StateContainerType), intent(out) :: container
      container = this%container
   end subroutine instance_get_container

   subroutine instance_get_process_manager(this, process_manager)
      class(CATChemInstanceType), intent(in) :: this
      type(ProcessManagerType), intent(out) :: process_manager
      process_manager = this%process_manager
   end subroutine instance_get_process_manager

   subroutine instance_get_diagnostic_manager(this, diag_manager)
      class(CATChemInstanceType), intent(in) :: this
      type(DiagnosticManagerType), intent(out) :: diag_manager
      ! TODO: Get diagnostic manager from container
   end subroutine instance_get_diagnostic_manager

   ! Field mapping registry methods (placeholder implementations)

   subroutine registry_initialize(this, rc)
      class(FieldMappingRegistryType), intent(inout) :: this
      integer, intent(out) :: rc
      rc = catchem_success
      ! Initialize empty registry
   end subroutine registry_initialize

   subroutine registry_add_met_mapping(this, mapping, rc)
      class(FieldMappingRegistryType), intent(inout) :: this
      type(FieldMappingType), intent(in) :: mapping
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Add met mapping to registry
   end subroutine registry_add_met_mapping

   subroutine registry_add_chem_mapping(this, mapping, rc)
      class(FieldMappingRegistryType), intent(inout) :: this
      type(FieldMappingType), intent(in) :: mapping
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Add chem mapping to registry
   end subroutine registry_add_chem_mapping

   subroutine registry_add_emis_mapping(this, mapping, rc)
      class(FieldMappingRegistryType), intent(inout) :: this
      type(FieldMappingType), intent(in) :: mapping
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Add emis mapping to registry
   end subroutine registry_add_emis_mapping

   subroutine registry_add_diag_mapping(this, mapping, rc)
      class(FieldMappingRegistryType), intent(inout) :: this
      type(FieldMappingType), intent(in) :: mapping
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Add diag mapping to registry
   end subroutine registry_add_diag_mapping

   ! FlexibleDataExchangeType methods (placeholder implementations)

   subroutine flex_setup_exchange(this, nx, ny, nz, rc)
      class(FlexibleDataExchangeType), intent(inout) :: this
      integer, intent(in) :: nx, ny, nz
      integer, intent(out) :: rc
      rc = catchem_success
      this%nx = nx; this%ny = ny; this%nz = nz
   end subroutine flex_setup_exchange

   subroutine flex_map_host_field(this, field_name, host_data, rc)
      class(FlexibleDataExchangeType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: host_data(:,:,:)
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Map host field to internal storage
   end subroutine flex_map_host_field

   subroutine flex_get_mapped_field(this, field_name, data, rc)
      class(FlexibleDataExchangeType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      real(fp), allocatable, intent(out) :: data(:,:,:)
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Retrieve mapped field data
   end subroutine flex_get_mapped_field

   subroutine flex_set_mapped_field(this, field_name, data, rc)
      class(FlexibleDataExchangeType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: data(:,:,:)
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Set mapped field data
   end subroutine flex_set_mapped_field

   subroutine flex_validate_mappings(this, rc)
      class(FlexibleDataExchangeType), intent(in) :: this
      integer, intent(out) :: rc
      rc = catchem_success
      ! TODO: Validate all field mappings
   end subroutine flex_validate_mappings
end module catchemapi_mod
```


