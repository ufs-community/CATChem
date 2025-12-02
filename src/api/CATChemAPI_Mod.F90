!> \file CATChemAPI_Mod.F90
!! \brief High-level CATChem API for easy integration
!! \ingroup catchem_api
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module provides a simplified, high-level API for integrating CATChem
!! into different modeling architectures. It abstracts away the complexity
!! of the StateContainer architecture while providing clean, simple interfaces
!! for the most common use cases.
!!
!! The API supports:
!! - Simple initialization with sensible defaults
!! - Easy data exchange with host models
!! - Automatic process management
!! - Column virtualization for performance
!! - Error handling with clear messages
!!
!! Advanced users can still access the detailed interfaces through the
!! underlying modules for fine-grained control.
!!
module CATChemAPI_Mod
   use precision_mod
   use error_mod
   use state_mod, only: StateContainerType, StateBuilderType
   use ConfigManager_Mod, only: ConfigManagerType
   use ProcessManager_Mod, only: ProcessManagerType
   use GridManager_Mod, only: GridManagerType
   use ColumnInterface_Mod, only: VirtualColumnType
   use ChemState_Mod, only: ChemStateType
   use MetState_Mod, only: MetStateType
   use EmisState_Mod, only: EmisStateType
   use DiagState_Mod, only: DiagStateType
   use DiagnosticInterface_Mod, only: DiagnosticManagerType, DiagnosticRegistryType, &
      DiagnosticFieldType

   implicit none
   private

   public :: CATChemInstanceType
   public :: CATChemConfigType
   public :: CATChemDataType
   public :: CATChemDiagnosticType
   public :: CATCHEM_SUCCESS, CATCHEM_FAILURE

   ! Return codes
   integer, parameter :: CATCHEM_SUCCESS = 0
   integer, parameter :: CATCHEM_FAILURE = -1

   !> Simple configuration type for basic CATChem setup
   type :: CATChemConfigType
      ! Grid dimensions
      integer :: nx = 0                    !< Number of grid points in x
      integer :: ny = 0                    !< Number of grid points in y
      integer :: nz = 0                    !< Number of grid points in z
      ! NOTE: nspecies is now read from configuration files, not set directly

      ! Time parameters
      real(fp) :: dt = 3600.0_fp           !< Time step [s]
      integer :: nsteps = 1                !< Number of time steps

      ! Process configuration
      logical :: enable_dust = .false.      !< Enable dust emissions
      logical :: enable_seasalt = .false.   !< Enable sea salt emissions
      logical :: enable_drydep = .false.    !< Enable dry deposition
      logical :: enable_chemistry = .false. !< Enable chemistry (future)
      logical :: enable_external_emis = .false. !< Enable external emissions

      ! Run phase configuration
      logical :: enable_run_phases = .false.   !< Enable multi-phase execution
      character(len=64), allocatable :: run_phase_names(:)  !< Names of run phases to execute

      ! Performance options
      logical :: use_column_processing = .true.  !< Use column virtualization
      logical :: enable_diagnostics = .true.     !< Enable diagnostic output

      ! Configuration files - these define the initialization order
      character(len=256) :: config_file = ''     !< Main config file path (read first)
      character(len=256) :: species_file = ''    !< Species definition file (read second)
      character(len=256) :: emission_file = ''   !< Emission config file (read third)
   end type CATChemConfigType

   !> Field mapping type for flexible data exchange
   type :: FieldMappingType
      character(len=64) :: host_field_name = ''    !< Host model field name
      character(len=64) :: catchem_field_name = '' !< CATChem internal field name
      character(len=32) :: field_type = ''         !< 'met', 'chem', 'emis', 'diag'
      character(len=32) :: units = ''              !< Field units
      real(fp) :: scale_factor = 1.0_fp            !< Unit conversion factor
      real(fp) :: offset = 0.0_fp                  !< Unit conversion offset
      logical :: is_required = .true.              !< Whether field is required
      logical :: is_3d = .true.                    !< Whether field is 3D (false = 2D)
   end type FieldMappingType

   !> Field mapping registry for managing host-CATChem field relationships
   type :: FieldMappingRegistryType
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
   end type FieldMappingRegistryType

   !> Enhanced data exchange type with flexible field mapping
   type :: FlexibleDataExchangeType
      ! Field registry
      type(FieldMappingRegistryType) :: field_mapping_registry

      ! Dynamic data storage
      real(fp), allocatable :: met_data_3d(:,:,:,:)    !< 3D meteorological data
      real(fp), allocatable :: met_data_2d(:,:,:)      !< 2D meteorological data
      real(fp), allocatable :: chem_data(:,:,:,:)      !< Chemical species data
      real(fp), allocatable :: emis_data_3d(:,:,:,:)   !< 3D emission data
      real(fp), allocatable :: emis_data_2d(:,:,:)     !< 2D emission data
      real(fp), allocatable :: diag_data_3d(:,:,:,:)   !< 3D diagnostic data
      real(fp), allocatable :: diag_data_2d(:,:,:)     !< 2D diagnostic data

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
   end type FlexibleDataExchangeType

   !> Data exchange type for simple data passing
   type :: CATChemDataType
      ! Meteorological inputs
      real(fp), allocatable :: temperature(:,:,:)     !< Temperature [K]
      real(fp), allocatable :: pressure(:,:,:)        !< Pressure [Pa]
      real(fp), allocatable :: humidity(:,:,:)        !< Specific humidity [kg/kg]
      real(fp), allocatable :: wind_u(:,:,:)          !< U wind component [m/s]
      real(fp), allocatable :: wind_v(:,:,:)          !< V wind component [m/s]
      real(fp), allocatable :: surface_pressure(:,:)  !< Surface pressure [Pa]
      real(fp), allocatable :: roughness_length(:,:)  !< Surface roughness [m]

      ! Chemical species concentrations
      real(fp), allocatable :: concentrations(:,:,:,:) !< Species concentrations [mol/mol]

      ! Emission data (optional)
      real(fp), allocatable :: emission_rates(:,:,:,:) !< Emission rates [mol/m2/s or mol/m3/s]

      ! Diagnostic outputs
      real(fp), allocatable :: dust_emissions(:,:,:)   !< Dust emission flux [kg/m2/s]
      real(fp), allocatable :: seasalt_emissions(:,:,:) !< Sea salt emission flux [kg/m2/s]
      real(fp), allocatable :: drydep_velocity(:,:,:)  !< Dry deposition velocity [m/s]
   end type CATChemDataType

   !> Diagnostic field information type
   type :: CATChemDiagnosticType
      character(len=64) :: field_name = ''        !< Diagnostic field name
      character(len=128) :: description = ''      !< Field description
      character(len=32) :: units = ''             !< Field units
      character(len=32) :: process_name = ''      !< Source process name
      integer :: dimensions = 0                   !< Number of dimensions (2D or 3D)
      logical :: is_available = .false.           !< Whether field is available
      real(fp), allocatable :: data_2d(:,:)       !< 2D diagnostic data
      real(fp), allocatable :: data_3d(:,:,:)     !< 3D diagnostic data
      real(fp) :: scalar_data = 0.0_fp            !< Scalar diagnostic data
   end type CATChemDiagnosticType

   !> Main CATChem instance type - simplified interface
   type :: CATChemInstanceType
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
   end type CATChemInstanceType

contains

   !> Initialize CATChem with basic configuration
   subroutine instance_init(this, config, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(CATChemConfigType), intent(in) :: config
      integer, intent(out) :: rc

      type(StateBuilderType) :: builder
      character(len=256) :: error_msg

      rc = CATCHEM_SUCCESS
      this%config = config

      ! Validate basic configuration
      if (config%nx <= 0 .or. config%ny <= 0 .or. config%nz <= 0) then
         rc = CATCHEM_FAILURE
         return
      endif

      ! Store grid dimensions
      this%nx = config%nx
      this%ny = config%ny
      this%nz = config%nz

      ! Initialize configuration manager
      call this%config_manager%init(rc)
      if (rc /= CATCHEM_SUCCESS) return

      ! Load configuration file if provided
      if (len_trim(config%config_file) > 0) then
         call this%config_manager%load_config(config%config_file, rc)
         if (rc /= CATCHEM_SUCCESS) return
      endif

      ! Build state container with sensible defaults
      call builder%init()
      call builder%set_grid_dimensions(config%nx, config%ny, config%nz)
      call builder%set_num_species(config%nspecies)
      call builder%enable_column_processing(config%use_column_processing)
      call builder%enable_diagnostics(config%enable_diagnostics)

      ! Build the container
      call builder%build(this%container, rc)
      if (rc /= CATCHEM_SUCCESS) return

      ! Get error manager for future use
      this%error_mgr => this%container%get_error_manager()

      ! Initialize process manager
      call this%process_manager%init(this%container, rc)
      if (rc /= CATCHEM_SUCCESS) return

      this%is_initialized = .true.

   end subroutine instance_init

   !> Modern initialization from configuration files (recommended)
   !! Follows proper initialization order: main config → species YAML → emission YAML
   subroutine instance_init_from_config_files(this, config, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(CATChemConfigType), intent(in) :: config
      integer, intent(out) :: rc

      type(StateBuilderType) :: builder
      integer :: temp_nspecies

      rc = CATCHEM_SUCCESS
      this%config = config

      ! Validate that all required config files are provided
      if (len_trim(config%config_file) == 0) then
         rc = CATCHEM_FAILURE
         return
      endif

      ! Initialize configuration manager
      call this%config_manager%init(rc)
      if (rc /= CATCHEM_SUCCESS) return

      ! Step 1: Read main configuration file (sets grid, processes, etc.)
      call this%config_manager%load_config(config%config_file, rc)
      if (rc /= CATCHEM_SUCCESS) return

      ! Extract grid dimensions from main config
      call this%config_manager%get_grid_dimensions(this%nx, this%ny, this%nz, rc)
      if (rc /= CATCHEM_SUCCESS) return

      ! Step 2: Read species configuration file (determines nspecies)
      if (len_trim(config%species_file) > 0) then
         call this%config_manager%load_species_data(config%species_file, rc)
         if (rc /= CATCHEM_SUCCESS) return
      else
         ! Try to get species file from main config
         call this%config_manager%get_species_filename(config%species_file, rc)
         if (rc == CATCHEM_SUCCESS .and. len_trim(config%species_file) > 0) then
            call this%config_manager%load_species_data(config%species_file, rc)
            if (rc /= CATCHEM_SUCCESS) return
         endif
      endif

      ! Get number of species from loaded configuration
      call this%config_manager%get_num_species(this%nspecies, rc)
      if (rc /= CATCHEM_SUCCESS) return

      ! Step 3: Read emission configuration file (if provided)
      if (len_trim(config%emission_file) > 0) then
         call this%config_manager%load_emission_config(config%emission_file, rc)
         if (rc /= CATCHEM_SUCCESS) return
      endif

      ! Build state container with configuration-driven parameters
      call builder%init()
      call builder%set_grid_dimensions(this%nx, this%ny, this%nz)
      call builder%set_num_species(this%nspecies)
      call builder%enable_column_processing(config%use_column_processing)
      call builder%enable_diagnostics(config%enable_diagnostics)

      ! Build the container
      call builder%build(this%container, rc)
      if (rc /= CATCHEM_SUCCESS) return

      ! Get error manager for future use
      this%error_mgr => this%container%get_error_manager()

      ! Initialize process manager
      call this%process_manager%init(this%container, rc)
      if (rc /= CATCHEM_SUCCESS) return

      ! Setup run phases if enabled
      if (config%enable_run_phases) then
         call this%setup_run_phases(config%run_phase_names, rc)
         if (rc /= CATCHEM_SUCCESS) return
      endif

      this%is_initialized = .true.

   end subroutine instance_init_from_config_files

   !> Setup grid geometry and coordinate system
   subroutine instance_setup_grid(this, lats, lons, levels, rc)
      class(CATChemInstanceType), intent(inout) :: this
      real(fp), intent(in) :: lats(:,:)      !< Latitude coordinates [degrees]
      real(fp), intent(in) :: lons(:,:)      !< Longitude coordinates [degrees]
      real(fp), intent(in) :: levels(:)      !< Vertical levels [Pa or m]
      integer, intent(out) :: rc

      type(GridManagerType), pointer :: grid_mgr

      rc = CATCHEM_SUCCESS

      if (.not. this%is_initialized) then
         rc = CATCHEM_FAILURE
         return
      endif

      ! Get grid manager and setup coordinates
      grid_mgr => this%container%get_grid_manager_ptr()
      if (.not. associated(grid_mgr)) then
         rc = CATCHEM_FAILURE
         return
      endif

      ! Setup coordinate system
      call grid_mgr%set_coordinates(lats, lons, levels, rc)
      if (rc /= CATCHEM_SUCCESS) return

      ! Enable column virtualization if requested
      if (this%config%use_column_processing) then
         call grid_mgr%enable_virtualization(rc)
         if (rc /= CATCHEM_SUCCESS) return
      endif

   end subroutine instance_setup_grid

   !> Add a process to the simulation
   subroutine instance_add_process(this, process_name, process_config, rc)
      class(CATChemInstanceType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      character(len=*), intent(in), optional :: process_config
      integer, intent(out) :: rc

      rc = CATCHEM_SUCCESS

      if (.not. this%is_initialized) then
         rc = CATCHEM_FAILURE
         return
      endif

      ! Add process to manager
      if (present(process_config)) then
         call this%process_manager%add_process(process_name, process_config, rc)
      else
         call this%process_manager%add_process(process_name, rc=rc)
      endif

      if (rc == CATCHEM_SUCCESS) then
         this%is_ready = .true.
      endif

   end subroutine instance_add_process

   !> Setup run phases for multi-phase execution
   subroutine instance_setup_run_phases(this, phase_names, rc)
      class(CATChemInstanceType), intent(inout) :: this
      character(len=64), intent(in) :: phase_names(:)
      integer, intent(out) :: rc

      integer :: i, num_phases

      rc = CATCHEM_SUCCESS

      if (.not. this%is_initialized) then
         rc = CATCHEM_FAILURE
         return
      endif

      num_phases = size(phase_names)
      if (num_phases <= 0) then
         rc = CATCHEM_FAILURE
         return
      endif

      ! Store phase names
      if (allocated(this%active_phases)) deallocate(this%active_phases)
      allocate(this%active_phases(num_phases))
      this%active_phases = phase_names

      ! Configure process manager for multi-phase execution
      call this%process_manager%configure_run_phases(phase_names, rc)
      if (rc /= CATCHEM_SUCCESS) return

      this%multi_phase_mode = .true.
      this%current_phase = 0

   end subroutine instance_setup_run_phases

   !> Run a specific phase
   subroutine instance_run_phase(this, phase_name, rc)
      class(CATChemInstanceType), intent(inout) :: this
      character(len=*), intent(in) :: phase_name
      integer, intent(out) :: rc

      integer :: phase_index

      rc = CATCHEM_SUCCESS

      if (.not. this%is_ready_to_run()) then
         rc = CATCHEM_FAILURE
         return
      endif

      if (.not. this%multi_phase_mode) then
         rc = CATCHEM_FAILURE
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
         rc = CATCHEM_FAILURE
         return
      endif

      ! Run the specific phase
      call this%process_manager%run_phase(phase_name, this%container, rc)
      if (rc /= CATCHEM_SUCCESS) return

      this%current_phase = phase_index

   end subroutine instance_run_phase

   !> Run all configured phases in sequence
   subroutine instance_run_all_phases(this, rc)
      class(CATChemInstanceType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i

      rc = CATCHEM_SUCCESS

      if (.not. this%multi_phase_mode) then
         ! Fall back to regular run_all if phases not configured
         call this%process_manager%run_all_processes(this%container, rc)
         return
      endif

      ! Run each phase in sequence
      do i = 1, size(this%active_phases)
         call this%instance_run_phase(this%active_phases(i), rc)
         if (rc /= CATCHEM_SUCCESS) return
      enddo

   end subroutine instance_run_all_phases

   !> Get names of configured run phases
   subroutine instance_get_phase_names(this, phase_names, rc)
      class(CATChemInstanceType), intent(in) :: this
      character(len=64), allocatable, intent(out) :: phase_names(:)
      integer, intent(out) :: rc

      rc = CATCHEM_SUCCESS

      if (.not. this%multi_phase_mode) then
         rc = CATCHEM_FAILURE
         return
      endif

      if (allocated(this%active_phases)) then
         allocate(phase_names(size(this%active_phases)))
         phase_names = this%active_phases
      else
         rc = CATCHEM_FAILURE
      endif

   end subroutine instance_get_phase_names

   !> Get current phase index
   function instance_get_current_phase(this) result(phase_index)
      class(CATChemInstanceType), intent(in) :: this
      integer :: phase_index

      phase_index = this%current_phase

   end function instance_get_current_phase

   !> Configure process manager with custom settings
   subroutine instance_configure_process_manager(this, max_processes, enable_column_batching, rc)
      class(CATChemInstanceType), intent(inout) :: this
      integer, intent(in), optional :: max_processes
      logical, intent(in), optional :: enable_column_batching
      integer, intent(out) :: rc

      rc = CATCHEM_SUCCESS

      if (.not. this%is_initialized) then
         rc = CATCHEM_FAILURE
         return
      endif

      if (present(max_processes)) then
         call this%process_manager%set_max_processes(max_processes, rc)
         if (rc /= CATCHEM_SUCCESS) return
      endif

      if (present(enable_column_batching)) then
         call this%process_manager%enable_column_batching(enable_column_batching, rc)
         if (rc /= CATCHEM_SUCCESS) return
      endif

   end subroutine instance_configure_process_manager

   !> Set a custom process manager (advanced usage)
   subroutine instance_set_process_manager(this, process_manager, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(ProcessManagerType), intent(in) :: process_manager
      integer, intent(out) :: rc

      rc = CATCHEM_SUCCESS

      if (.not. this%is_initialized) then
         rc = CATCHEM_FAILURE
         return
      endif

      this%process_manager = process_manager

   end subroutine instance_set_process_manager

   !> Run a single time step
   subroutine instance_run_timestep(this, input_data, output_data, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(CATChemDataType), intent(in) :: input_data
      type(CATChemDataType), intent(out) :: output_data
      integer, intent(out) :: rc

      rc = CATCHEM_SUCCESS

      if (.not. this%is_ready_to_run()) then
         rc = CATCHEM_FAILURE
         return
      endif

      ! Set input data
      call this%set_meteorology(input_data, rc)
      if (rc /= CATCHEM_SUCCESS) return

      call this%set_concentrations(input_data, rc)
      if (rc /= CATCHEM_SUCCESS) return

      if (allocated(input_data%emission_rates)) then
         call this%set_emissions(input_data, rc)
         if (rc /= CATCHEM_SUCCESS) return
      endif

      ! Run processes - use phase-based execution if configured
      if (this%multi_phase_mode) then
         call this%run_all_phases(rc)
      else
         call this%process_manager%run_all_processes(this%container, rc)
      endif
      if (rc /= CATCHEM_SUCCESS) return

      ! Get output data
      call this%get_concentrations(output_data, rc)
      if (rc /= CATCHEM_SUCCESS) return

      call this%get_diagnostics(output_data, rc)
      if (rc /= CATCHEM_SUCCESS) return

   end subroutine instance_run_timestep

   !> Run multiple time steps
   subroutine instance_run_multiple_steps(this, nsteps, input_data, output_data, rc)
      class(CATChemInstanceType), intent(inout) :: this
      integer, intent(in) :: nsteps
      type(CATChemDataType), intent(in) :: input_data(:)
      type(CATChemDataType), intent(out) :: output_data(:)
      integer, intent(out) :: rc

      integer :: i

      rc = CATCHEM_SUCCESS

      if (size(input_data) /= nsteps .or. size(output_data) /= nsteps) then
         rc = CATCHEM_FAILURE
         return
      endif

      do i = 1, nsteps
         call this%run_timestep(input_data(i), output_data(i), rc)
         if (rc /= CATCHEM_SUCCESS) return
      enddo

   end subroutine instance_run_multiple_steps

   !> Check if ready to run
   function instance_is_ready(this) result(is_ready)
      class(CATChemInstanceType), intent(in) :: this
      logical :: is_ready

      is_ready = this%is_initialized .and. this%is_ready

   end function instance_is_ready

   !> Set meteorological data
   subroutine instance_set_meteorology(this, data, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(CATChemDataType), intent(in) :: data
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state

      rc = CATCHEM_SUCCESS

      met_state => this%container%get_met_state_ptr()
      if (.not. associated(met_state)) then
         rc = CATCHEM_FAILURE
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

   !> Set chemical concentrations
   subroutine instance_set_concentrations(this, data, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(CATChemDataType), intent(in) :: data
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state

      rc = CATCHEM_SUCCESS

      chem_state => this%container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         rc = CATCHEM_FAILURE
         return
      endif

      if (allocated(data%concentrations)) then
         chem_state%concentrations = data%concentrations
      endif

   end subroutine instance_set_concentrations

   !> Set emission data
   subroutine instance_set_emissions(this, data, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(CATChemDataType), intent(in) :: data
      integer, intent(out) :: rc

      type(EmisStateType), pointer :: emis_state

      rc = CATCHEM_SUCCESS

      emis_state => this%container%get_emis_state_ptr()
      if (.not. associated(emis_state)) then
         rc = CATCHEM_FAILURE
         return
      endif

      if (allocated(data%emission_rates)) then
         emis_state%emission_rates = data%emission_rates
      endif

   end subroutine instance_set_emissions

   !> Get concentrations from CATChem
   subroutine instance_get_concentrations(this, data, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(CATChemDataType), intent(out) :: data
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state

      rc = CATCHEM_SUCCESS

      chem_state => this%container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         rc = CATCHEM_FAILURE
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

   !> Get diagnostic data
   subroutine instance_get_diagnostics(this, data, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(CATChemDataType), intent(out) :: data
      integer, intent(out) :: rc

      type(DiagStateType), pointer :: diag_state

      rc = CATCHEM_SUCCESS

      diag_state => this%container%get_diag_state_ptr()
      if (.not. associated(diag_state)) then
         rc = CATCHEM_FAILURE
         return
      endif

      ! Extract diagnostic fields if available
      ! TODO: Add specific diagnostic field extraction based on available processes

   end subroutine instance_get_diagnostics

   !> Get error message
   subroutine instance_get_error_message(this, error_msg)
      class(CATChemInstanceType), intent(in) :: this
      character(len=*), intent(out) :: error_msg

      if (associated(this%error_mgr)) then
         call this%error_mgr%get_error_message(error_msg)
      else
         error_msg = 'No error manager available'
      endif

   end subroutine instance_get_error_message

   !> Validate input data
   subroutine instance_validate_data(this, data, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(CATChemDataType), intent(in) :: data
      integer, intent(out) :: rc

      rc = CATCHEM_SUCCESS

      ! Basic validation - check array dimensions match expected grid
      if (allocated(data%temperature)) then
         if (size(data%temperature,1) /= this%nx .or. &
            size(data%temperature,2) /= this%ny .or. &
            size(data%temperature,3) /= this%nz) then
            rc = CATCHEM_FAILURE
            return
         endif
      endif

      if (allocated(data%concentrations)) then
         if (size(data%concentrations,1) /= this%nx .or. &
            size(data%concentrations,2) /= this%ny .or. &
            size(data%concentrations,3) /= this%nz .or. &
            size(data%concentrations,4) /= this%nspecies) then
            rc = CATCHEM_FAILURE
            return
         endif
      endif

   end subroutine instance_validate_data

   !> Get species names
   subroutine instance_get_species_names(this, species_names, rc)
      class(CATChemInstanceType), intent(in) :: this
      character(len=64), allocatable, intent(out) :: species_names(:)
      integer, intent(out) :: rc

      rc = CATCHEM_SUCCESS

      ! TODO: Extract species names from configuration
      allocate(species_names(this%nspecies))
      ! For now, use generic names
      ! In real implementation, this would come from the species configuration file

   end subroutine instance_get_species_names

   !> Finalize CATChem instance
   subroutine instance_finalize(this, rc)
      class(CATChemInstanceType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CATCHEM_SUCCESS

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
      rc = CATCHEM_SUCCESS
      ! TODO: Implement diagnostic retrieval
   end subroutine instance_get_available_diagnostics

   subroutine instance_get_process_diagnostic(this, process_name, diag_name, diagnostic, rc)
      class(CATChemInstanceType), intent(in) :: this
      character(len=*), intent(in) :: process_name, diag_name
      type(CATChemDiagnosticType), intent(out) :: diagnostic
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Implement specific diagnostic retrieval
   end subroutine instance_get_process_diagnostic

   subroutine instance_get_all_process_diagnostics(this, diagnostics, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(CATChemDiagnosticType), allocatable, intent(out) :: diagnostics(:)
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Implement all diagnostics retrieval
   end subroutine instance_get_all_process_diagnostics

   subroutine instance_list_process_diagnostics(this, process_names, diagnostic_info, rc)
      class(CATChemInstanceType), intent(in) :: this
      character(len=64), allocatable, intent(out) :: process_names(:)
      character(len=256), allocatable, intent(out) :: diagnostic_info(:)
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Implement diagnostic listing
   end subroutine instance_list_process_diagnostics

   ! Field mapping methods (placeholder implementations)
   subroutine instance_fill_met_state_from_host(this, data_exchange, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(FlexibleDataExchangeType), intent(in) :: data_exchange
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Implement field mapping from host
   end subroutine instance_fill_met_state_from_host

   subroutine instance_extract_met_state_to_host(this, data_exchange, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(FlexibleDataExchangeType), intent(out) :: data_exchange
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Implement field mapping to host
   end subroutine instance_extract_met_state_to_host

   subroutine instance_fill_chem_state_from_host(this, data_exchange, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(FlexibleDataExchangeType), intent(in) :: data_exchange
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Implement chemical state mapping from host
   end subroutine instance_fill_chem_state_from_host

   subroutine instance_extract_chem_state_to_host(this, data_exchange, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(FlexibleDataExchangeType), intent(out) :: data_exchange
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Implement chemical state mapping to host
   end subroutine instance_extract_chem_state_to_host

   subroutine instance_fill_emis_state_from_host(this, data_exchange, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(FlexibleDataExchangeType), intent(in) :: data_exchange
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Implement emission state mapping from host
   end subroutine instance_fill_emis_state_from_host

   subroutine instance_extract_emis_state_to_host(this, data_exchange, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(FlexibleDataExchangeType), intent(out) :: data_exchange
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Implement emission state mapping to host
   end subroutine instance_extract_emis_state_to_host

   subroutine instance_extract_diagnostics_to_host(this, data_exchange, rc)
      class(CATChemInstanceType), intent(in) :: this
      type(FlexibleDataExchangeType), intent(out) :: data_exchange
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Implement diagnostic mapping to host
   end subroutine instance_extract_diagnostics_to_host

   subroutine instance_setup_field_mapping(this, mappings, rc)
      class(CATChemInstanceType), intent(inout) :: this
      type(FieldMappingType), intent(in) :: mappings(:)
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Implement field mapping setup
   end subroutine instance_setup_field_mapping

   subroutine instance_add_field_mapping(this, field_type, catchem_name, host_name, data_type, rc)
      class(CATChemInstanceType), intent(inout) :: this
      character(len=*), intent(in) :: field_type, catchem_name, host_name, data_type
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Implement individual field mapping addition
   end subroutine instance_add_field_mapping

   subroutine instance_validate_field_mappings(this, rc)
      class(CATChemInstanceType), intent(in) :: this
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
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
      rc = CATCHEM_SUCCESS
      ! Initialize empty registry
   end subroutine registry_initialize

   subroutine registry_add_met_mapping(this, mapping, rc)
      class(FieldMappingRegistryType), intent(inout) :: this
      type(FieldMappingType), intent(in) :: mapping
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Add met mapping to registry
   end subroutine registry_add_met_mapping

   subroutine registry_add_chem_mapping(this, mapping, rc)
      class(FieldMappingRegistryType), intent(inout) :: this
      type(FieldMappingType), intent(in) :: mapping
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Add chem mapping to registry
   end subroutine registry_add_chem_mapping

   subroutine registry_add_emis_mapping(this, mapping, rc)
      class(FieldMappingRegistryType), intent(inout) :: this
      type(FieldMappingType), intent(in) :: mapping
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Add emis mapping to registry
   end subroutine registry_add_emis_mapping

   subroutine registry_add_diag_mapping(this, mapping, rc)
      class(FieldMappingRegistryType), intent(inout) :: this
      type(FieldMappingType), intent(in) :: mapping
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Add diag mapping to registry
   end subroutine registry_add_diag_mapping

   ! FlexibleDataExchangeType methods (placeholder implementations)

   subroutine flex_setup_exchange(this, nx, ny, nz, rc)
      class(FlexibleDataExchangeType), intent(inout) :: this
      integer, intent(in) :: nx, ny, nz
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      this%nx = nx; this%ny = ny; this%nz = nz
   end subroutine flex_setup_exchange

   subroutine flex_map_host_field(this, field_name, host_data, rc)
      class(FlexibleDataExchangeType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: host_data(:,:,:)
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Map host field to internal storage
   end subroutine flex_map_host_field

   subroutine flex_get_mapped_field(this, field_name, data, rc)
      class(FlexibleDataExchangeType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      real(fp), allocatable, intent(out) :: data(:,:,:)
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Retrieve mapped field data
   end subroutine flex_get_mapped_field

   subroutine flex_set_mapped_field(this, field_name, data, rc)
      class(FlexibleDataExchangeType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: data(:,:,:)
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Set mapped field data
   end subroutine flex_set_mapped_field

   subroutine flex_validate_mappings(this, rc)
      class(FlexibleDataExchangeType), intent(in) :: this
      integer, intent(out) :: rc
      rc = CATCHEM_SUCCESS
      ! TODO: Validate all field mappings
   end subroutine flex_validate_mappings
end module CATChemAPI_Mod
