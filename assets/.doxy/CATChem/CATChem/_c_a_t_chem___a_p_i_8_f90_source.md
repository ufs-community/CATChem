

# File CATChem\_API.F90

[**File List**](files.md) **>** [**api**](dir_da61e3e9a357748887e3ca8d7c5a0c16.md) **>** [**CATChem\_API.F90**](_c_a_t_chem___a_p_i_8_f90.md)

[Go to the documentation of this file](_c_a_t_chem___a_p_i_8_f90.md)


```Fortran

module catchem_api
   use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure, errormanagertype
   use catchemcore_mod, only: catchemcoretype, catchembuildertype
   use statemanager_mod, only: statemanagertype
   use processmanager_mod, only: processmanagertype
   use gridmanager_mod, only: gridmanagertype
   use diagnosticmanager_mod, only: diagnosticmanagertype
   use metstate_mod, only: metstatetype
   use chemstate_mod, only: chemstatetype
   use configmanager_mod, only: configdatatype
   use diagnosticinterface_mod, only: diagnosticregistrytype, diagnosticfieldtype, &
      diag_real_scalar, diag_real_1d, diag_real_2d, diag_real_3d, &
      diag_integer_scalar, diag_integer_1d, diag_integer_2d, diag_integer_3d
   use processinterface_mod, only: processinterface
   ! Import process registration functions
   use seasaltprocesscreator_mod, only: register_seasalt_process
   use drydepprocesscreator_mod, only: register_drydep_process
   use wetdepprocesscreator_mod, only: register_wetdep_process
   use settlingprocesscreator_mod, only: register_settling_process

   implicit none
   private

   ! Public interface types and constants
   public :: catchem_model

   type :: catchem_model
      private
      ! Core engine - uses existing CATChemCore infrastructure
      type(CATChemCoreType) :: core

      ! Configuration and status tracking
      logical :: initialized = .false.
      logical :: grid_setup = .false.
      logical :: enable_run_phase = .false.
      character(len=512) :: config_file = ''
      character(len=64), allocatable, public :: required_fields(:)
      type(ErrorManagerType) :: error_manager

      ! Grid information
      integer :: nx = 0, ny = 0, nz = 0
      integer :: nsoil = 4, nsoiltype = 19, nsurftype = 20

   contains
      ! Basic lifecycle methods
      procedure :: initialize => model_initialize
      procedure :: finalize => model_finalize

      ! Grid access (no separate setup needed now)
      procedure :: get_grid_dimensions => model_get_grid_dimensions

      ! Process management
      procedure :: add_process => model_add_process
      procedure :: get_process_names => model_get_process_names
      procedure :: get_num_processes => model_get_num_processes
      procedure, private :: model_register_process

      ! Run execution
      procedure :: run_timestep => model_run_timestep
      procedure :: run_phase => model_run_phase
      procedure :: run_all_phases => model_run_all_phases
      procedure :: get_phase_names => model_get_phase_names

      ! Data exchange methods
      procedure :: set_chemistry => model_set_chemistry
      procedure :: get_chemistry => model_get_chemistry

      ! Diagnostic methods
      procedure :: get_diagnostic_names => model_get_diagnostic_names
      procedure :: get_diagnostic => model_get_diagnostic
      procedure :: get_all_diagnostics => model_get_all_diagnostics

      ! Utility methods
      procedure :: is_ready => model_is_ready
      procedure :: is_initialized => model_is_initialized
      procedure :: get_required_met_index  => model_get_required_met_index
      procedure :: get_diag_index_from_field  => model_get_diag_index_from_field

      ! Core access methods (for advanced users)
      procedure :: get_state_manager => model_get_state_manager
      procedure :: get_process_manager => model_get_process_manager
      procedure :: get_error_manager => model_get_error_manager
      procedure :: get_grid_manager => model_get_grid_manager
      procedure :: get_diagnostic_manager => model_get_diagnostic_manager
   end type catchem_model

contains

   subroutine model_initialize(this, config_file, nx, ny, nz, nsoil, nsoiltype, nsurftype, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: config_file
      integer, intent(in) :: nx, ny, nz
      integer, intent(in), optional :: nsoil, nsoiltype, nsurftype
      integer, intent(out) :: rc

      type(CATChemBuilderType) :: builder
      type(ConfigDataType), pointer :: config_data => null()

      rc = cc_success
      call this%error_manager%init()

      ! Validate inputs
      if (len_trim(config_file) == 0) then
         call this%error_manager%push_context('model_initialize', 'validating configuration file')
         call this%error_manager%report_error(1001, 'Configuration file path is empty', rc)
         call this%error_manager%pop_context()
         return
      endif

      if (nx <= 0 .or. ny <= 0 .or. nz <= 0) then
         call this%error_manager%push_context('model_initialize', 'validating grid dimensions')
         call this%error_manager%report_error(1001, 'Grid dimensions must be positive', rc)
         call this%error_manager%pop_context()
         return
      endif

      ! Store config file path and grid dimensions
      this%config_file = trim(config_file)
      this%nx = nx
      this%ny = ny
      this%nz = nz

      ! Store soil and surface parameters (use provided values or defaults)
      if (present(nsoil)) this%nsoil = nsoil
      if (present(nsoiltype)) this%nsoiltype = nsoiltype
      if (present(nsurftype)) this%nsurftype = nsurftype

      ! Initialize core using builder pattern with grid information
      call builder%init()
      builder = builder%with_name('CATChem_API_Instance')
      builder = builder%with_config(config_file)
      builder = builder%with_grid(nx, ny, nz, this%nsoil, this%nsoiltype, this%nsurftype)
      call builder%build(this%core, rc)

      if (rc /= cc_success) then
         call this%error_manager%push_context('model_initialize', 'building CATChem core')
         call this%error_manager%report_error(1014, 'Failed to initialize CATChem core with config: ' // trim(config_file), rc)
         call this%error_manager%pop_context()
         return
      endif

      this%initialized = .true.
      this%grid_setup = .true.  ! Grid is now set up during initialization

      ! Get configuration data from core
      config_data => this%core%get_config()
      if ( .not. associated(config_data)) then
         call this%error_manager%push_context('model_initialize', 'accessing configuration data')
         call this%error_manager%report_error(1002, 'Required managers or config data not available', rc)
         call this%error_manager%pop_context()
         return
      endif
      this%enable_run_phase = config_data%run_phases_enabled

   end subroutine model_initialize

   subroutine model_finalize(this, rc)
      class(CATChem_Model), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success

      if (.not. this%initialized) then
         rc = cc_success  ! Already finalized
         return
      endif

      ! Finalize core
      call this%core%finalize(rc)
      if (rc /= cc_success) then
         call this%error_manager%push_context('model_finalize', 'finalizing CATChem core')
         call this%error_manager%report_error(1015, 'Core finalization had issues', rc)
         call this%error_manager%pop_context()
         ! Don't return failure for finalization warnings
         rc = cc_success
      endif

      ! Reset state
      this%initialized = .false.
      this%grid_setup = .false.
      this%enable_run_phase = .false.
      this%nx = 0
      this%ny = 0
      this%nz = 0
      this%nsoil = 4
      this%nsoiltype = 19
      this%nsurftype = 20
      this%config_file = ''
      if (allocated(this%required_fields)) deallocate(this%required_fields)
   end subroutine model_finalize

   subroutine model_get_grid_dimensions(this, nx, ny, nz, nsoil, nsoiltype, nsurftype)
      class(CATChem_Model), intent(in) :: this
      integer, intent(out) :: nx, ny, nz
      integer, intent(out), optional :: nsoil, nsoiltype, nsurftype

      nx = this%nx
      ny = this%ny
      nz = this%nz
      if (present(nsoil)) nsoil = this%nsoil
      if (present(nsoiltype)) nsoiltype = this%nsoiltype
      if (present(nsurftype)) nsurftype = this%nsurftype
   end subroutine model_get_grid_dimensions

   subroutine model_add_process(this, rc)
      use configmanager_mod, only: configdatatype
      class(CATChem_Model), intent(inout) :: this
      integer, intent(out) :: rc

      type(ConfigDataType), pointer :: config_data => null()
      type(ProcessManagerType), pointer :: process_mgr => null()
      integer :: i, reg_rc, add_rc, num_fields

      rc = cc_success

      if (.not. this%initialized) then
         call this%error_manager%push_context('model_add_process', 'checking initialization status')
         call this%error_manager%report_error(1003, 'Model must be initialized first', rc)
         call this%error_manager%pop_context()
         return
      endif

      ! Get configuration data from core
      config_data => this%core%get_config()
      if (.not. associated(config_data)) then
         call this%error_manager%push_context('model_add_process', 'accessing configuration data')
         call this%error_manager%report_error(1002, 'Configuration data not available', rc)
         call this%error_manager%pop_context()
         return
      endif

      ! Get process manager
      process_mgr => this%core%get_process_manager()
      if (.not. associated(process_mgr)) then
         call this%error_manager%push_context('model_add_process', 'accessing process manager')
         call this%error_manager%report_error(1014, 'ProcessManager not available from core', rc)
         call this%error_manager%pop_context()
         return
      endif

      ! Check if processes are available (either run phase or direct processes)
      if (.not. allocated(config_data%run_phase_processes)) then
         call this%error_manager%push_context('model_add_process', 'checking process configuration')
         call this%error_manager%report_error(1002, 'No processes configured', rc)
         call this%error_manager%pop_context()
         return
      endif

      ! Loop through all processes in configuration and add enabled ones
      do i = 1, size(config_data%run_phase_processes)
         if (config_data%run_phase_processes(i)%enabled) then

            ! Automatically register the process based on its name
            call this%model_register_process(config_data%run_phase_processes(i)%name, process_mgr, reg_rc)
            if (reg_rc /= cc_success) then
               write(*,'(A,A,A)') 'Warning: Failed to register enabled process: ', &
                  trim(config_data%run_phase_processes(i)%name), '. Skipping this process.'
               cycle  ! Skip this process and continue with others
            endif

            ! Add process to core
            call this%core%add_process(config_data%run_phase_processes(i)%name, add_rc)
            if (add_rc /= cc_success) then
               write(*,'(A,A,A)') 'Warning: Failed to add enabled process to core: ', &
                  trim(config_data%run_phase_processes(i)%name), '. Skipping this process.'
               cycle  ! Skip this process and continue with others
            endif

            write(*,'(A,A)') 'Successfully added enabled process: ', &
               trim(config_data%run_phase_processes(i)%name)
         else
            write(*,'(A,A)') 'Skipping disabled process: ', &
               trim(config_data%run_phase_processes(i)%name)
         endif
      end do

      ! Get required met fields from process_mgr
      if (allocated(process_mgr%required_met_fields)) then
         num_fields = size(process_mgr%required_met_fields)
         allocate(this%required_fields(num_fields))
         this%required_fields = process_mgr%required_met_fields
      else
         call this%error_manager%push_context('model_add_process', 'checking required met fields')
         call this%error_manager%report_error(1014, 'No met fields found', rc)
         call this%error_manager%pop_context()
      endif

   end subroutine model_add_process

   subroutine model_register_process(this, process_name, process_mgr, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = cc_success

      ! Dispatch to the appropriate registration function based on process name
      select case (trim(process_name))
       case ('seasalt')
         call register_seasalt_process(process_mgr, rc)
         if (rc /= cc_success) then
            call this%error_manager%push_context('model_register_process', 'registering seasalt process')
            call this%error_manager%report_error(1014, 'Failed to register seasalt process', rc)
            call this%error_manager%pop_context()
         endif

         ! Add more processes here as they become available
         ! case ('dust')
         !    call register_dust_process(process_mgr, rc)
       case ('drydep')
         call register_drydep_process(process_mgr, rc)
         if (rc /= cc_success) then
            call this%error_manager%push_context('model_register_process', 'registering drydep process')
            call this%error_manager%report_error(1014, 'Failed to register drydep process', rc)
            call this%error_manager%pop_context()
         endif
       case ('wetdep')
         call register_wetdep_process(process_mgr, rc)
         if (rc /= cc_success) then
            call this%error_manager%push_context('model_register_process', 'registering wetdep process')
            call this%error_manager%report_error(1014, 'Failed to register wetdep process', rc)
            call this%error_manager%pop_context()
         endif
       case ('settling')
         call register_settling_process(process_mgr, rc)
         if (rc /= cc_success) then
            call this%error_manager%push_context('model_register_process', 'registering settling process')
            call this%error_manager%report_error(1014, 'Failed to register settling process', rc)
            call this%error_manager%pop_context()
         endif
         ! case ('chemistry')
         !    call register_chemistry_process(process_mgr, rc)

       case default
         call this%error_manager%push_context('model_register_process', 'validating process type')
         call this%error_manager%report_error(1016, 'Unknown process type: ' // trim(process_name) // &
            '. Supported processes: seasalt', rc)
         call this%error_manager%pop_context()
      end select

   end subroutine model_register_process

   subroutine model_get_process_names(this, process_names, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), allocatable, intent(out) :: process_names(:)
      integer, intent(out) :: rc

      type(ProcessManagerType), pointer :: process_mgr => null()
      character(len=64) :: temp_names(50)  ! Temporary array with max size
      integer :: count, i

      rc = cc_success

      if (.not. this%initialized) then
         allocate(process_names(0))
         rc = cc_failure
         return
      endif

      ! Get process manager from core
      process_mgr => this%core%get_process_manager()
      if (.not. associated(process_mgr)) then
         allocate(process_names(0))
         rc = cc_failure
         return
      endif

      ! Get process list from ProcessManager
      call process_mgr%list_processes(temp_names, count)

      ! Allocate output array with actual count
      allocate(process_names(count))

      ! Copy the actual process names
      do i = 1, count
         process_names(i) = temp_names(i)
      end do

   end subroutine model_get_process_names

   function model_get_num_processes(this) result(num_processes)
      class(CATChem_Model), intent(inout) :: this
      integer :: num_processes

      type(ProcessManagerType), pointer :: process_mgr => null()
      character(len=64) :: temp_names(50)  ! Temporary array with max size

      num_processes = 0

      if (.not. this%initialized) return

      ! Get process manager from core
      process_mgr => this%core%get_process_manager()
      if (.not. associated(process_mgr)) return

      ! Get process count from ProcessManager
      call process_mgr%list_processes(temp_names, num_processes)

   end function model_get_num_processes


   subroutine model_run_timestep(this, timestep, dt, rc)
      class(CATChem_Model), intent(inout) :: this
      integer, intent(in) :: timestep         ! Current timestep number
      real(fp), intent(in) :: dt              ! Timestep size [s]
      integer, intent(out) :: rc

      rc = cc_success

      if (.not. this%is_ready()) then
         call this%error_manager%push_context('model_run_timestep', 'checking model readiness')
         call this%error_manager%report_error(1003, 'Model is not ready to run timestep', rc)
         call this%error_manager%pop_context()
         return
      endif

      if (dt <= 0.0_fp) then
         call this%error_manager%push_context('model_run_timestep', 'validating timestep size')
         call this%error_manager%report_error(1001, 'Timestep size must be positive', rc)
         call this%error_manager%pop_context()
         return
      endif

      if (this%enable_run_phase) then
         call this%run_all_phases(rc)
         if (rc /= cc_success) then
            call this%error_manager%push_context('model_run_timestep', 'running all phases')
            call this%error_manager%report_error(1015, 'Failed to run all phases during timestep', rc)
            call this%error_manager%pop_context()
            return
         endif
      else
         ! Run the core timestep
         call this%core%run_timestep(timestep, dt, rc)
         if (rc /= cc_success) then
            call this%error_manager%push_context('model_run_timestep', 'running core timestep')
            call this%error_manager%report_error(1015, 'Failed to run all processes during timestep', rc)
            call this%error_manager%pop_context()
            return
         endif
      endif
   end subroutine model_run_timestep

   subroutine model_run_phase(this, phase_name, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: phase_name
      integer, intent(out) :: rc

      type(ProcessManagerType), pointer :: process_mgr => null()
      type(StateManagerType), pointer :: state_mgr => null()
      type(ConfigDataType), pointer :: config_data => null()

      rc = cc_success

      if (.not. this%is_ready()) then
         call this%error_manager%push_context('model_run_phase', 'checking model readiness')
         call this%error_manager%report_error(1003, 'Model is not ready for phase execution', rc)
         call this%error_manager%pop_context()
         return
      endif

      ! Get managers and config data
      process_mgr => this%core%get_process_manager()
      state_mgr => this%core%get_state_manager()
      config_data => this%core%get_config()

      if (.not. associated(process_mgr) .or. .not. associated(state_mgr) .or. .not. associated(config_data)) then
         call this%error_manager%push_context('model_run_phase', 'accessing required managers')
         call this%error_manager%report_error(1014, 'Required managers or config data not available', rc)
         call this%error_manager%pop_context()
         return
      endif

      ! Run the specific phase via ProcessManager
      call process_mgr%run_phase(phase_name, config_data, state_mgr, rc)
      if (rc /= cc_success) then
         call this%error_manager%push_context('model_run_phase', 'executing phase: ' // trim(phase_name))
         call this%error_manager%report_error(1015, 'Failed to run phase: ' // trim(phase_name), rc)
         call this%error_manager%pop_context()
      endif

   end subroutine model_run_phase

   subroutine model_run_all_phases(this, rc)
      class(CATChem_Model), intent(inout) :: this
      integer, intent(out) :: rc

      type(ProcessManagerType), pointer :: process_mgr => null()
      type(StateManagerType), pointer :: state_mgr => null()
      type(ConfigDataType), pointer :: config_data => null()

      rc = cc_success

      if (.not. this%is_ready()) then
         call this%error_manager%push_context('model_run_all_phases', 'checking model readiness')
         call this%error_manager%report_error(1003, 'Model is not ready for phase execution', rc)
         call this%error_manager%pop_context()
         return
      endif

      ! Get managers and config data
      process_mgr => this%core%get_process_manager()
      state_mgr => this%core%get_state_manager()
      config_data => this%core%get_config()

      if (.not. associated(process_mgr) .or. .not. associated(state_mgr) .or. .not. associated(config_data)) then
         call this%error_manager%push_context('model_run_all_phases', 'accessing required managers')
         call this%error_manager%report_error(1014, 'Required managers or config data not available', rc)
         call this%error_manager%pop_context()
         return
      endif

      ! Run all phases in sequence using ConfigManager data
      call process_mgr%run_all_phases(config_data, state_mgr, rc)
      if (rc /= cc_success) then
         call this%error_manager%push_context('model_run_all_phases', 'executing all phases')
         call this%error_manager%report_error(1015, 'Failed to run all phases', rc)
         call this%error_manager%pop_context()
      endif

   end subroutine model_run_all_phases

   subroutine model_get_phase_names(this, phase_names, rc)
      use configmanager_mod, only: configdatatype
      class(CATChem_Model), intent(inout) :: this
      character(len=*), allocatable, intent(out) :: phase_names(:)
      integer, intent(out) :: rc

      type(ConfigDataType), pointer :: config_data => null()
      integer :: i

      rc = cc_success

      if (.not. this%initialized) then
         allocate(phase_names(0))
         rc = cc_failure
         return
      endif

      ! Get configuration data from core
      config_data => this%core%get_config()
      if (.not. associated(config_data)) then
         allocate(phase_names(0))
         rc = cc_failure
         return
      endif

      ! Check if run phases are available
      if (.not. allocated(config_data%run_phases)) then
         allocate(phase_names(0))
         return
      endif

      allocate(phase_names(size(config_data%run_phases)))
      do i = 1, size(config_data%run_phases)
         phase_names(i) = config_data%run_phases(i)%name
      end do

   end subroutine model_get_phase_names


   subroutine model_set_chemistry(this, species_names, concentrations, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: species_names(:)
      real(fp), intent(in) :: concentrations(:,:,:,:)  ! [species, nx, ny, nz]
      integer, intent(out) :: rc

      type(StateManagerType), pointer :: state_mgr => null()
      type(ChemStateType), pointer :: chem_state => null()
      integer :: num_species

      rc = cc_success
      ! Initialize error handling

      if (.not. this%is_ready()) then
         ! Error: 'Model is not ready for chemistry data'
         rc = cc_failure
         return
      endif

      num_species = size(species_names)
      if (num_species /= size(concentrations, 1)) then
         ! Error: 'Number of species names does not match concentration dimensions'
         rc = cc_failure
         return
      endif

      ! Get state manager and chemistry state
      state_mgr => this%core%get_state_manager()
      if (.not. associated(state_mgr)) then
         ! Error: 'State manager not available'
         rc = cc_failure
         return
      endif

      chem_state => state_mgr%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         ! Error: 'Chemistry state not available'
         rc = cc_failure
         return
      endif

      ! TODO: Implement chemistry data transfer
      ! This would involve mapping species names to indices and copying data
      ! Error: 'Chemistry data transfer not yet implemented'
      rc = cc_failure

   end subroutine model_set_chemistry

   subroutine model_get_chemistry(this, species_names, concentrations, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), allocatable, intent(out) :: species_names(:)
      real(fp), allocatable, intent(out) :: concentrations(:,:,:,:)  ! [species, nx, ny, nz]
      integer, intent(out) :: rc

      type(StateManagerType), pointer :: state_mgr => null()
      type(ChemStateType), pointer :: chem_state => null()

      rc = cc_success

      if (.not. this%is_ready()) then
         rc = cc_failure
         return
      endif

      ! Get state manager and chemistry state
      state_mgr => this%core%get_state_manager()
      if (.not. associated(state_mgr)) then
         rc = cc_failure
         return
      endif

      chem_state => state_mgr%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         rc = cc_failure
         return
      endif

      ! TODO: Implement chemistry data retrieval
      ! This would involve getting species names and copying concentration data
      allocate(species_names(0))
      allocate(concentrations(0, this%nx, this%ny, this%nz))

   end subroutine model_get_chemistry


   subroutine model_get_diagnostic_names(this, diagnostic_names, diagnostic_fields, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), allocatable, intent(out) :: diagnostic_names(:)
      character(len=*), allocatable, optional, intent(out) :: diagnostic_fields(:)
      integer, intent(out) :: rc

      type(DiagnosticManagerType), pointer :: diag_mgr => null()
      type(DiagnosticRegistryType), pointer :: registry => null()
      character(len=64), allocatable :: process_list(:), field_names(:)
      integer :: num_processes, i, j, field_count, total_fields, name_idx
      integer :: local_rc

      rc = cc_success
      ! Initialize error handling

      if (.not. this%initialized) then
         allocate(diagnostic_names(0))
         ! Error: 'Model not initialized'
         rc = cc_failure
         return
      endif

      diag_mgr => this%core%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) then
         allocate(diagnostic_names(0))
         ! Error: 'Diagnostic manager not available'
         rc = cc_failure
         return
      endif

      ! Get list of all processes with diagnostics
      call diag_mgr%list_processes(process_list, num_processes, local_rc)
      if (local_rc /= cc_success .or. num_processes == 0) then
         allocate(diagnostic_names(0))
         return  ! No processes or failed to get list - return empty array
      endif

      ! Count total diagnostic fields across all processes
      total_fields = 0
      do i = 1, num_processes
         call diag_mgr%get_process_registry(process_list(i), registry, local_rc)
         if (local_rc == cc_success .and. associated(registry)) then
            total_fields = total_fields + registry%get_field_count()
         endif
      end do

      ! Allocate output array
      allocate(diagnostic_names(total_fields))
      if (present(diagnostic_fields))  allocate( diagnostic_fields(total_fields) )

      ! Collect all diagnostic field names with process prefix
      name_idx = 0
      do i = 1, num_processes
         call diag_mgr%get_process_registry(process_list(i), registry, local_rc)
         if (local_rc == cc_success .and. associated(registry)) then
            field_count = registry%get_field_count()
            if (field_count > 0) then
               allocate(field_names(field_count))
               call registry%list_fields(field_names, field_count)

               do j = 1, field_count
                  name_idx = name_idx + 1
                  ! Create qualified name: process_name.field_name
                  diagnostic_names(name_idx) = trim(process_list(i)) // '.' // trim(field_names(j))
                  if (present(diagnostic_fields)) diagnostic_fields(name_idx) = trim(field_names(j))
               end do

               deallocate(field_names)
            endif
         endif
      end do

      ! Clean up
      if (allocated(process_list)) deallocate(process_list)

   end subroutine model_get_diagnostic_names

   subroutine model_get_diagnostic(this, diagnostic_name, diagnostic_data, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: diagnostic_name
      real(fp), allocatable, intent(out) :: diagnostic_data(:,:,:)
      integer, intent(out) :: rc

      type(DiagnosticManagerType), pointer :: diag_mgr => null()
      type(DiagnosticRegistryType), pointer :: registry => null()
      character(len=64), allocatable :: process_list(:), field_names(:)
      character(len=64) :: process_name, field_name
      integer :: num_processes, i, j, field_count, dot_pos, data_type
      integer :: local_rc
      real(fp) :: scalar_value
      real(fp), pointer :: array_1d_ptr(:) => null()
      real(fp), pointer :: array_2d_ptr(:,:) => null()
      real(fp), pointer :: array_3d_ptr(:,:,:) => null()
      logical :: found = .false.

      rc = cc_success
      ! Initialize error handling

      if (.not. this%initialized) then
         ! Error: 'Model not initialized'
         rc = cc_failure
         return
      endif

      diag_mgr => this%core%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) then
         ! Error: 'Diagnostic manager not available'
         rc = cc_failure
         return
      endif

      ! Parse diagnostic name (format: process_name.field_name)
      dot_pos = index(diagnostic_name, '.')
      if (dot_pos <= 1) then
         ! Error: 'Invalid diagnostic name format. Expected: process_name.field_name'
         rc = cc_failure
         return
      endif

      process_name = diagnostic_name(1:dot_pos-1)
      field_name = diagnostic_name(dot_pos+1:)

      ! Get field value using DiagnosticManager
      call diag_mgr%get_field_value(process_name, field_name, &
         scalar_value, array_1d_ptr, array_2d_ptr, array_3d_ptr, &
         data_type, rc = local_rc)

      if (local_rc /= cc_success) then
         ! Error: 'Failed to retrieve diagnostic field: ' // trim(diagnostic_name)
         rc = cc_failure
         return
      endif

      ! Allocate and populate output array based on data type
      allocate(diagnostic_data(this%nx, this%ny, this%nz))
      diagnostic_data = 0.0_fp

      select case (data_type)
       case (diag_real_scalar)
         diagnostic_data = scalar_value
         found = .true.

       case (diag_real_2d)
         if (associated(array_2d_ptr)) then
            ! Copy 2D data to first level of 3D array
            diagnostic_data(:,:,1) = array_2d_ptr
            found = .true.
         endif

       case (diag_real_3d)
         if (associated(array_3d_ptr)) then
            diagnostic_data = array_3d_ptr
            found = .true.
         endif

       case default
         ! Error: 'Unsupported diagnostic data type for: ' // trim(diagnostic_name)
         rc = cc_failure
         return
      end select

      if (.not. found) then
         ! Error: 'No data available for diagnostic: ' // trim(diagnostic_name)
         rc = cc_failure
      endif

   end subroutine model_get_diagnostic

   subroutine model_get_all_diagnostics(this, diagnostic_names, diagnostic_data, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), allocatable, intent(out) :: diagnostic_names(:)
      real(fp), allocatable, intent(out) :: diagnostic_data(:,:,:,:)  ! [diag, nx, ny, nz]
      integer, intent(out) :: rc

      type(DiagnosticManagerType), pointer :: diag_mgr => null()
      integer :: local_rc

      rc = cc_success
      ! Initialize error handling

      if (.not. this%initialized) then
         ! Error: 'Model not initialized'
         allocate(diagnostic_names(0))
         allocate(diagnostic_data(0, this%nx, this%ny, this%nz))
         rc = cc_failure
         return
      endif

      diag_mgr => this%core%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) then
         ! Error: 'Diagnostic manager not available'
         allocate(diagnostic_names(0))
         allocate(diagnostic_data(0, this%nx, this%ny, this%nz))
         rc = cc_failure
         return
      endif

      ! First get all diagnostic field names
      call this%get_diagnostic_names(diagnostic_names, rc =local_rc)
      if (local_rc /= cc_success) then
         ! Error: 'Failed to get diagnostic names'
         allocate(diagnostic_data(0, this%nx, this%ny, this%nz))
         rc = cc_failure
         return
      endif

      ! Collect all diagnostic data using DiagnosticManager
      call diag_mgr%collect_all_diagnostics(local_rc)
      if (local_rc /= cc_success) then
         ! Error: 'Failed to collect diagnostic data'
         allocate(diagnostic_data(size(diagnostic_names), this%nx, this%ny, this%nz))
         rc = cc_failure
         return
      endif

      ! For now, allocate the data array but don't populate it
      ! TODO: Extract individual diagnostic values and populate diagnostic_data array
      ! This would require iterating through diagnostic_names and calling get_field_value for each
      allocate(diagnostic_data(size(diagnostic_names), this%nx, this%ny, this%nz))
      diagnostic_data = 0.0_fp

      ! Note: The actual data collection has been performed and stored in each field's
      ! DiagnosticDataType structure. Individual fields can be accessed using get_diagnostic()

   end subroutine model_get_all_diagnostics

   function model_get_diag_index_from_field(this, field_name) result(found_index)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      character(len=128), allocatable :: diagnostic_names(:)
      character(len=128), allocatable :: diagnostic_fields(:)
      integer :: found_index, rc, i

      found_index = 0
      call this%get_diagnostic_names(diagnostic_names, diagnostic_fields, rc)
      if (allocated(diagnostic_fields)) then
         do i = 1, size(diagnostic_fields)
            if (trim(field_name) == trim(diagnostic_fields(i))) then
               found_index = i
               exit !exit the loop once found; TODO: this assumes we do not have duplicated field names
            end if
         end do
      else
         write(*,'(A)') ' Warning: No diagnostic fields are registered!!'
      end if
   end function model_get_diag_index_from_field

   function model_is_ready(this) result(is_ready)
      class(CATChem_Model), intent(inout) :: this
      logical :: is_ready

      is_ready = this%initialized .and. this%grid_setup .and. (this%get_num_processes() > 0)
   end function model_is_ready

   function model_is_initialized(this) result(is_initialized)
      class(CATChem_Model), intent(in) :: this
      logical :: is_initialized

      is_initialized = this%initialized
   end function model_is_initialized

   function model_get_required_met_index(this, var_name) result(found_index)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer :: found_index
      integer :: i
      type(ProcessManagerType), pointer :: process_mgr

      found_index = 0
      if (allocated(this%required_fields)) then
         do i = 1, size(this%required_fields)
            if (trim(var_name) == trim(this%required_fields(i))) then
               found_index = i
               exit !exit the loop once found
            end if
         end do
      else
         write(*,'(A)') ' Warning: No met fields are registered!!'
      end if
   end function model_get_required_met_index

   function model_get_error_manager(this) result(error_mgr_ptr)
      class(CATChem_Model), intent(inout), target :: this
      type(ErrorManagerType), pointer :: error_mgr_ptr

      error_mgr_ptr => this%error_manager
   end function model_get_error_manager

   ! Core access methods for advanced users

   function model_get_state_manager(this) result(state_mgr_ptr)
      class(CATChem_Model), intent(inout) :: this
      type(StateManagerType), pointer :: state_mgr_ptr

      state_mgr_ptr => this%core%get_state_manager()
   end function model_get_state_manager

   function model_get_process_manager(this) result(process_mgr_ptr)
      class(CATChem_Model), intent(inout) :: this
      type(ProcessManagerType), pointer :: process_mgr_ptr

      process_mgr_ptr => this%core%get_process_manager()
   end function model_get_process_manager

   function model_get_grid_manager(this) result(grid_mgr_ptr)
      class(CATChem_Model), intent(inout) :: this
      type(GridManagerType), pointer :: grid_mgr_ptr

      grid_mgr_ptr => this%core%get_grid_manager()
   end function model_get_grid_manager

   function model_get_diagnostic_manager(this) result(diag_mgr_ptr)
      class(CATChem_Model), intent(inout) :: this
      type(DiagnosticManagerType), pointer :: diag_mgr_ptr

      diag_mgr_ptr => this%core%get_diagnostic_manager()
   end function model_get_diagnostic_manager

end module catchem_api
```


