

# File CATChem\_HighLevel\_API.F90

[**File List**](files.md) **>** [**api**](dir_da61e3e9a357748887e3ca8d7c5a0c16.md) **>** [**CATChem\_HighLevel\_API.F90**](_c_a_t_chem___high_level___a_p_i_8_f90.md)

[Go to the documentation of this file](_c_a_t_chem___high_level___a_p_i_8_f90.md)


```Fortran

module catchem_highlevel_api
   use precision_mod
   use error_mod
   use state_mod, only: statecontainertype, statebuildertype
   use processmanager_mod, only: processmanagertype
   use processfactory_mod, only: processfactorytype
   use gridmanager_mod, only: gridmanagertype
   use chemstate_mod, only: chemstatetype
   use metstate_mod, only: metstatetype
   use emisstate_mod, only: emisstatetype
   use diagstate_mod, only: diagstatetype
   use config_opt_mod, only: configtype
   use configmanager_mod, only: configmanagertype

   implicit none
   private

   public :: catchemmodeltype
   public :: catchemconfigtype
   public :: catchemdatatype
   public :: catchemdiagnosticstype

   ! High-level API procedures
   public :: catchem_init
   public :: catchem_run_timestep
   public :: catchem_finalize
   public :: catchem_add_process
   public :: catchem_set_emissions
   public :: catchem_set_meteorology
   public :: catchem_get_concentrations
   public :: catchem_get_diagnostics

   type :: catchemconfigtype
      ! Grid configuration
      integer :: nx = 0
      integer :: ny = 0
      integer :: nz = 0
      integer :: n_species = 0

      ! File paths
      character(len=256) :: config_file = ''
      character(len=256) :: species_file = ''
      character(len=256) :: mechanism_file = ''

      ! Processing options
      logical :: use_column_processing = .true.     
      logical :: enable_diagnostics = .true.        
      logical :: validate_inputs = .true.           
      logical :: verbose_output = .false.           

      ! Timing
      real(fp) :: timestep = 3600.0_fp             

      ! Optional advanced settings (can be ignored for simple use cases)
      character(len=64) :: chemistry_solver = 'default'
      character(len=64) :: transport_scheme = 'default'
      real(fp) :: solver_tolerance = 1.0e-6_fp          
   end type catchemconfigtype

   type :: catchemdatatype
      ! Meteorological data
      real(fp), allocatable :: temperature(:,:,:)
      real(fp), allocatable :: pressure(:,:,:)
      real(fp), allocatable :: humidity(:,:,:)
      real(fp), allocatable :: wind_u(:,:,:)
      real(fp), allocatable :: wind_v(:,:,:)
      real(fp), allocatable :: wind_w(:,:,:)
      real(fp), allocatable :: surface_pressure(:,:)
      real(fp), allocatable :: surface_temperature(:,:)

      ! Chemical concentrations
      real(fp), allocatable :: concentrations(:,:,:,:)
      character(len=32), allocatable :: species_names(:)

      ! Emission data
      real(fp), allocatable :: surface_emissions(:,:,:)
      real(fp), allocatable :: volume_emissions(:,:,:,:)

      ! Grid information
      real(fp), allocatable :: latitude(:,:)
      real(fp), allocatable :: longitude(:,:)
      real(fp), allocatable :: grid_area(:,:)
      real(fp), allocatable :: level_heights(:)

      logical :: is_allocated = .false.               
   end type catchemdatatype

   type :: catchemdiagnosticstype
      ! Process diagnostics
      real(fp), allocatable :: dust_emissions(:,:)
      real(fp), allocatable :: seasalt_emissions(:,:)
      real(fp), allocatable :: dry_deposition(:,:,:)
      real(fp), allocatable :: wet_deposition(:,:,:)

      ! Chemistry diagnostics
      real(fp), allocatable :: reaction_rates(:,:,:,:)
      real(fp), allocatable :: photolysis_rates(:,:,:,:)

      ! Performance diagnostics
      real(fp) :: total_runtime = 0.0_fp              
      real(fp) :: chemistry_time = 0.0_fp             
      real(fp) :: emission_time = 0.0_fp              
      real(fp) :: transport_time = 0.0_fp             

      logical :: is_allocated = .false.               
   end type catchemdiagnosticstype

   type :: catchemmodeltype
      private

      ! Core CATChem components (hidden from user)
      type(StateContainerType) :: container
      type(ProcessManagerType) :: process_manager
      type(ProcessFactoryType) :: process_factory
      type(ConfigManagerType) :: config_manager

      ! Configuration
      type(CATChemConfigType) :: config

      ! Status
      logical :: is_initialized = .false.
      logical :: has_processes = .false.
      character(len=256) :: last_error = ''

      ! Performance tracking
      real(fp) :: total_runtime = 0.0_fp
      integer :: n_timesteps = 0

   contains
      ! Main interface methods
      procedure :: init => model_init
      procedure :: add_process => model_add_process
      procedure :: run_timestep => model_run_timestep
      procedure :: finalize => model_finalize

      ! Data interface methods
      procedure :: set_meteorology => model_set_meteorology
      procedure :: set_emissions => model_set_emissions
      procedure :: set_concentrations => model_set_concentrations
      procedure :: get_concentrations => model_get_concentrations
      procedure :: get_diagnostics => model_get_diagnostics

      ! Utility methods
      procedure :: is_ready => model_is_ready
      procedure :: get_status => model_get_status
      procedure :: get_performance_stats => model_get_performance_stats
      procedure :: validate_inputs => model_validate_inputs

      ! Advanced interface (access to detailed components)
      procedure :: get_state_container => model_get_state_container
      procedure :: get_process_manager => model_get_process_manager
      procedure :: get_config_manager => model_get_config_manager
   end type catchemmodeltype

contains

   !========================================================================
   ! High-Level API Functions (Simplified Interface)
   !========================================================================

   subroutine catchem_init(model, config, rc)
      type(CATChemModelType), intent(out) :: model
      type(CATChemConfigType), intent(in) :: config
      integer, intent(out) :: rc

      call model%init(config, rc)
   end subroutine catchem_init

   subroutine catchem_run_timestep(model, met_data, emis_data, rc)
      type(CATChemModelType), intent(inout) :: model
      type(CATChemDataType), intent(in) :: met_data
      type(CATChemDataType), intent(in) :: emis_data
      integer, intent(out) :: rc

      ! Set input data
      call model%set_meteorology(met_data, rc)
      if (rc /= cc_success) return

      call model%set_emissions(emis_data, rc)
      if (rc /= cc_success) return

      ! Run timestep
      call model%run_timestep(rc)
   end subroutine catchem_run_timestep

   subroutine catchem_finalize(model, rc)
      type(CATChemModelType), intent(inout) :: model
      integer, intent(out) :: rc

      call model%finalize(rc)
   end subroutine catchem_finalize

   subroutine catchem_add_process(model, process_name, process_config, rc)
      type(CATChemModelType), intent(inout) :: model
      character(len=*), intent(in) :: process_name
      character(len=*), intent(in), optional :: process_config
      integer, intent(out) :: rc

      call model%add_process(process_name, process_config, rc)
   end subroutine catchem_add_process

   subroutine catchem_set_emissions(model, emis_data, rc)
      type(CATChemModelType), intent(inout) :: model
      type(CATChemDataType), intent(in) :: emis_data
      integer, intent(out) :: rc

      call model%set_emissions(emis_data, rc)
   end subroutine catchem_set_emissions

   subroutine catchem_set_meteorology(model, met_data, rc)
      type(CATChemModelType), intent(inout) :: model
      type(CATChemDataType), intent(in) :: met_data
      integer, intent(out) :: rc

      call model%set_meteorology(met_data, rc)
   end subroutine catchem_set_meteorology

   subroutine catchem_get_concentrations(model, concentrations, species_names, rc)
      type(CATChemModelType), intent(in) :: model
      real(fp), allocatable, intent(out) :: concentrations(:,:,:,:)
      character(len=32), allocatable, intent(out) :: species_names(:)
      integer, intent(out) :: rc

      call model%get_concentrations(concentrations, species_names, rc)
   end subroutine catchem_get_concentrations

   subroutine catchem_get_diagnostics(model, diagnostics, rc)
      type(CATChemModelType), intent(in) :: model
      type(CATChemDiagnosticsType), intent(out) :: diagnostics
      integer, intent(out) :: rc

      call model%get_diagnostics(diagnostics, rc)
   end subroutine catchem_get_diagnostics

   !========================================================================
   ! CATChemModelType Implementation
   !========================================================================

   subroutine model_init(this, config, rc)
      class(CATChemModelType), intent(out) :: this
      type(CATChemConfigType), intent(in) :: config
      integer, intent(out) :: rc

      type(StateBuilderType) :: builder
      type(ErrorManagerType), pointer :: error_mgr
      character(len=256) :: message

      rc = cc_success
      this%config = config

      ! Validate configuration
      if (config%nx <= 0 .or. config%ny <= 0 .or. config%nz <= 0) then
         this%last_error = 'Invalid grid dimensions'
         rc = cc_failure
         return
      endif

      if (len_trim(config%config_file) == 0) then
         this%last_error = 'Configuration file not specified'
         rc = cc_failure
         return
      endif

      ! Initialize configuration manager
      call this%config_manager%init(rc)
      if (rc /= cc_success) then
         this%last_error = 'Failed to initialize configuration manager'
         return
      endif

      call this%config_manager%load_config(config%config_file, rc)
      if (rc /= cc_success) then
         this%last_error = 'Failed to load configuration file: ' // trim(config%config_file)
         return
      endif

      ! Build state container
      call builder%init()
      call builder%set_grid_dimensions(config%nx, config%ny, config%nz)
      call builder%set_species_count(config%n_species)
      call builder%enable_column_processing(config%use_column_processing)
      call builder%enable_diagnostics(config%enable_diagnostics)
      call builder%set_verbose(config%verbose_output)

      call builder%build(this%container, rc)
      if (rc /= cc_success) then
         this%last_error = 'Failed to build state container'
         return
      endif

      ! Initialize process factory and manager
      call this%process_factory%init(rc)
      if (rc /= cc_success) then
         this%last_error = 'Failed to initialize process factory'
         return
      endif

      call this%process_manager%init(this%container, rc)
      if (rc /= cc_success) then
         this%last_error = 'Failed to initialize process manager'
         return
      endif

      this%is_initialized = .true.

      ! Log successful initialization
      error_mgr => this%container%get_error_manager()
      write(message, '(A,I0,A,I0,A,I0,A)') 'CATChem initialized: grid (', &
         config%nx, 'x', config%ny, 'x', config%nz, '), column processing enabled'
      call error_mgr%report_info(message)

   end subroutine model_init

   subroutine model_add_process(this, process_name, process_config, rc)
      class(CATChemModelType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      character(len=*), intent(in), optional :: process_config
      integer, intent(out) :: rc

      class(ProcessInterface), pointer :: process
      character(len=256) :: config_file

      rc = cc_success

      if (.not. this%is_initialized) then
         this%last_error = 'Model not initialized'
         rc = cc_failure
         return
      endif

      ! Use provided config or default
      if (present(process_config)) then
         config_file = process_config
      else
         config_file = ''  ! Use default configuration
      endif

      ! Create process using factory
      call this%process_factory%create_process(process_name, config_file, process, rc)
      if (rc /= cc_success) then
         this%last_error = 'Failed to create process: ' // trim(process_name)
         return
      endif

      ! Register process with manager
      call this%process_manager%register_process(process, rc)
      if (rc /= cc_success) then
         this%last_error = 'Failed to register process: ' // trim(process_name)
         return
      endif

      this%has_processes = .true.

   end subroutine model_add_process

   subroutine model_run_timestep(this, rc)
      class(CATChemModelType), intent(inout) :: this
      integer, intent(out) :: rc

      real(fp) :: start_time, end_time
      type(ErrorManagerType), pointer :: error_mgr

      rc = cc_success

      if (.not. this%is_ready()) then
         this%last_error = 'Model not ready for execution'
         rc = cc_failure
         return
      endif

      ! Start timing
      call cpu_time(start_time)

      ! Run all processes
      call this%process_manager%run_all_processes(this%container, rc)
      if (rc /= cc_success) then
         this%last_error = 'Process execution failed'
         return
      endif

      ! Update timing
      call cpu_time(end_time)
      this%total_runtime = this%total_runtime + (end_time - start_time)
      this%n_timesteps = this%n_timesteps + 1

      ! Log progress (every 100 timesteps)
      if (mod(this%n_timesteps, 100) == 0) then
         error_mgr => this%container%get_error_manager()
         call error_mgr%report_info('Completed 100 timesteps')
      endif

   end subroutine model_run_timestep

   subroutine model_set_meteorology(this, met_data, rc)
      class(CATChemModelType), intent(inout) :: this
      type(CATChemDataType), intent(in) :: met_data
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state

      rc = cc_success

      if (.not. this%is_initialized) then
         this%last_error = 'Model not initialized'
         rc = cc_failure
         return
      endif

      ! Validate input data
      call this%validate_inputs(met_data, rc)
      if (rc /= cc_success) return

      ! Get meteorological state
      met_state => this%container%get_met_state_ptr()
      if (.not. associated(met_state)) then
         this%last_error = 'Meteorological state not available'
         rc = cc_failure
         return
      endif

      ! Copy data to state (simplified - in reality would handle unit conversions, etc.)
      if (allocated(met_data%temperature)) then
         met_state%temperature = met_data%temperature
      endif

      if (allocated(met_data%pressure)) then
         met_state%pressure = met_data%pressure
      endif

      if (allocated(met_data%humidity)) then
         met_state%humidity = met_data%humidity
      endif

      ! Set additional fields as needed...

   end subroutine model_set_meteorology

   subroutine model_set_emissions(this, emis_data, rc)
      class(CATChemModelType), intent(inout) :: this
      type(CATChemDataType), intent(in) :: emis_data
      integer, intent(out) :: rc

      type(EmisStateType), pointer :: emis_state

      rc = cc_success

      if (.not. this%is_initialized) then
         this%last_error = 'Model not initialized'
         rc = cc_failure
         return
      endif

      ! Get emission state
      emis_state => this%container%get_emis_state_ptr()
      if (.not. associated(emis_state)) then
         this%last_error = 'Emission state not available'
         rc = cc_failure
         return
      endif

      ! Copy emission data (simplified - would handle species mapping, unit conversions, etc.)
      ! Implementation would depend on EmisState structure

   end subroutine model_set_emissions

   subroutine model_get_concentrations(this, concentrations, species_names, rc)
      class(CATChemModelType), intent(in) :: this
      real(fp), allocatable, intent(out) :: concentrations(:,:,:,:)
      character(len=32), allocatable, intent(out) :: species_names(:)
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state

      rc = cc_success

      if (.not. this%is_initialized) then
         this%last_error = 'Model not initialized'
         rc = cc_failure
         return
      endif

      ! Get chemical state
      chem_state => this%container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         this%last_error = 'Chemical state not available'
         rc = cc_failure
         return
      endif

      ! Copy concentration data (simplified)
      ! Implementation would depend on ChemState structure
      allocate(concentrations(this%config%nx, this%config%ny, this%config%nz, this%config%n_species))
      allocate(species_names(this%config%n_species))

      ! concentrations = chem_state%concentrations  ! Simplified
      ! species_names = chem_state%species_names    ! Simplified

   end subroutine model_get_concentrations

   subroutine model_get_diagnostics(this, diagnostics, rc)
      class(CATChemModelType), intent(in) :: this
      type(CATChemDiagnosticsType), intent(out) :: diagnostics
      integer, intent(out) :: rc

      rc = cc_success

      if (.not. this%is_initialized) then
         this%last_error = 'Model not initialized'
         rc = cc_failure
         return
      endif

      ! Collect diagnostic data from various sources
      diagnostics%total_runtime = this%total_runtime

      ! Get process-specific diagnostics
      ! Implementation would collect from diagnostic manager

      diagnostics%is_allocated = .true.

   end subroutine model_get_diagnostics

   function model_is_ready(this) result(ready)
      class(CATChemModelType), intent(in) :: this
      logical :: ready

      ready = this%is_initialized .and. this%has_processes
   end function model_is_ready

   function model_get_status(this) result(status)
      class(CATChemModelType), intent(in) :: this
      character(len=256) :: status

      if (.not. this%is_initialized) then
         status = 'Not initialized'
      else if (.not. this%has_processes) then
         status = 'Initialized, no processes registered'
      else
         write(status, '(A,I0,A)') 'Ready, ', this%n_timesteps, ' timesteps completed'
      endif
   end function model_get_status

   subroutine model_get_performance_stats(this, total_time, avg_time_per_step, n_steps)
      class(CATChemModelType), intent(in) :: this
      real(fp), intent(out) :: total_time
      real(fp), intent(out) :: avg_time_per_step
      integer, intent(out) :: n_steps

      total_time = this%total_runtime
      n_steps = this%n_timesteps

      if (n_steps > 0) then
         avg_time_per_step = total_time / real(n_steps, fp)
      else
         avg_time_per_step = 0.0_fp
      endif
   end subroutine model_get_performance_stats

   subroutine model_validate_inputs(this, data, rc)
      class(CATChemModelType), intent(in) :: this
      type(CATChemDataType), intent(in) :: data
      integer, intent(out) :: rc

      rc = cc_success

      if (.not. this%config%validate_inputs) return

      ! Check grid dimensions match
      if (allocated(data%temperature)) then
         if (size(data%temperature,1) /= this%config%nx .or. &
            size(data%temperature,2) /= this%config%ny .or. &
            size(data%temperature,3) /= this%config%nz) then
            rc = cc_failure
            return
         endif
      endif

      ! Check for valid values
      if (allocated(data%temperature)) then
         if (any(data%temperature < 100.0_fp) .or. any(data%temperature > 400.0_fp)) then
            rc = cc_failure
            return
         endif
      endif

      ! Additional validation as needed...

   end subroutine model_validate_inputs

   subroutine model_finalize(this, rc)
      class(CATChemModelType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success

      ! Finalize components
      call this%process_manager%finalize(rc)
      call this%process_factory%finalize(rc)
      call this%container%finalize(rc)
      call this%config_manager%finalize(rc)

      this%is_initialized = .false.
      this%has_processes = .false.

   end subroutine model_finalize

   !========================================================================
   ! Advanced Interface Methods (Access to Detailed Components)
   !========================================================================

   function model_get_state_container(this) result(container)
      class(CATChemModelType), intent(in) :: this
      type(StateContainerType), pointer :: container

      container => null()
      if (this%is_initialized) then
         ! Return pointer to container (would need to be implemented carefully)
         ! container => this%container
      endif
   end function model_get_state_container

   function model_get_process_manager(this) result(manager)
      class(CATChemModelType), intent(in) :: this
      type(ProcessManagerType), pointer :: manager

      manager => null()
      if (this%is_initialized) then
         ! Return pointer to manager (would need to be implemented carefully)
         ! manager => this%process_manager
      endif
   end function model_get_process_manager

   function model_get_config_manager(this) result(manager)
      class(CATChemModelType), intent(in) :: this
      type(ConfigManagerType), pointer :: manager

      manager => null()
      if (this%is_initialized) then
         ! Return pointer to manager (would need to be implemented carefully)
         ! manager => this%config_manager
      endif
   end function model_get_config_manager

end module catchem_highlevel_api
```


