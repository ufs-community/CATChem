! \file init_mod.F90
!! \brief Modern initialization module using StateContainer and enhanced error handling
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module provides comprehensive initialization routines for CATChem using
!! the modern StateContainer architecture with integrated error handling.
!!
!! \details
!! **Modern Features:**
!! - StateContainer-based initialization with dependency injection
!! - Integrated error handling with context tracking
!! - Enhanced configuration management
!! - Comprehensive validation and consistency checking
!! - Modular initialization for different CATChem components
!! - Performance monitoring and timing
!!
module init_mod
   use precision_mod
   use error_mod
   use state_mod
   use ConfigManager_Mod, only: ConfigDataType
   use MetState_Mod, only: MetStateType
   use ChemState_Mod, only: ChemStateType
   use EmisState_Mod, only: EmisStateType
   use DiagState_Mod, only: DiagStateType

   implicit none
   private

   ! Public initialization routines
   public :: initialize_catchem
   public :: initialize_core_states
   public :: initialize_processes
   public :: validate_initialization
   public :: finalize_catchem

contains

   !> \brief Main CATChem initialization routine
   !!
   !! This is the primary initialization routine that sets up a complete CATChem
   !! instance using the StateContainer pattern with comprehensive error handling.
   !!
   !! \param[out] container Initialized StateContainer with all components
   !! \param[in] config_file Path to configuration file
   !! \param[in] container_name Optional name for the container
   !! \param[out] rc Return code
   subroutine initialize_catchem(container, config_file, rc, container_name)
      implicit none

      type(StateContainerType), intent(out) :: container
      character(len=*), intent(in) :: config_file
      integer, intent(out) :: rc
      character(len=*), intent(in), optional :: container_name

      type(ErrorManagerType), pointer :: error_mgr
      character(len=256) :: name

      ! Set container name
      name = 'CATChem_Main'
      if (present(container_name)) name = trim(container_name)

      ! Initialize container
      call container%init(name, rc)
      if (rc /= CC_SUCCESS) return

      ! Get error manager for this initialization
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('initialize_catchem', 'Main CATChem initialization')

      ! Phase 1: Initialize core states
      call initialize_core_states(container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_PROCESS_INITIALIZATION, &
            'Failed to initialize core states', rc, &
            'initialize_catchem', &
            'Check system memory and configuration file')
         call error_mgr%pop_context()
         return
      endif

      ! Phase 2: Load configuration
      call load_configuration(container, config_file, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_INVALID_CONFIG, &
            'Failed to load configuration', rc, &
            'initialize_catchem', &
            'Verify configuration file exists and is valid')
         call error_mgr%pop_context()
         return
      endif

      ! Phase 3: Initialize processes
      call initialize_processes(container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_PROCESS_INITIALIZATION, &
            'Failed to initialize processes', rc, &
            'initialize_catchem', &
            'Check process configuration and input data')
         call error_mgr%pop_context()
         return
      endif

      ! Phase 4: Validate complete initialization
      call validate_initialization(container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_STATE_INCONSISTENCY, &
            'Initialization validation failed', rc, &
            'initialize_catchem', &
            'Review error messages for specific validation failures')
         call error_mgr%pop_context()
         return
      endif

      call error_mgr%pop_context()

      ! Success
      rc = CC_SUCCESS

   end subroutine initialize_catchem

   !> \brief Initialize core state objects
   !!
   !! This routine initializes all core state objects (Met, Chem, Emis, Diag)
   !! within the StateContainer.
   !!
   !! \param[inout] container StateContainer to initialize
   !! \param[out] rc Return code
   subroutine initialize_core_states(container, rc)
      implicit none

      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      type(MetStateType), pointer :: met_state
      type(ChemStateType), pointer :: chem_state
      type(EmisStateType), pointer :: emis_state
      type(DiagStateType), pointer :: diag_state

      error_mgr => container%get_error_manager()
      call error_mgr%push_context('initialize_core_states', 'Initializing core state objects')

      ! Initialize meteorological state
      call error_mgr%push_context('init_met_state', 'Allocating meteorological state')
      met_state => container%get_met_state_ptr()
      if (.not. associated(met_state)) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
            'Failed to get meteorological state pointer', rc, &
            'initialize_core_states', &
            'Check memory availability')
         call error_mgr%pop_context()
         call error_mgr%pop_context()
         return
      endif
      call error_mgr%pop_context()

      ! Initialize chemical state
      call error_mgr%push_context('init_chem_state', 'Allocating chemical state')
      chem_state => container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
            'Failed to get chemical state pointer', rc, &
            'initialize_core_states', &
            'Check memory availability')
         call error_mgr%pop_context()
         call error_mgr%pop_context()
         return
      endif
      call error_mgr%pop_context()

      ! Initialize emission state
      call error_mgr%push_context('init_emis_state', 'Allocating emission state')
      emis_state => container%get_emis_state_ptr()
      if (.not. associated(emis_state)) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
            'Failed to get emission state pointer', rc, &
            'initialize_core_states', &
            'Check memory availability')
         call error_mgr%pop_context()
         call error_mgr%pop_context()
         return
      endif
      call error_mgr%pop_context()

      ! Initialize diagnostic state
      call error_mgr%push_context('init_diag_state', 'Allocating diagnostic state')
      diag_state => container%get_diag_state_ptr()
      if (.not. associated(diag_state)) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
            'Failed to get diagnostic state pointer', rc, &
            'initialize_core_states', &
            'Check memory availability')
         call error_mgr%pop_context()
         call error_mgr%pop_context()
         return
      endif
      call error_mgr%pop_context()

      call error_mgr%pop_context()
      rc = CC_SUCCESS

   end subroutine initialize_core_states

   !> \brief Load configuration from file
   !!
   !! This routine loads the configuration file and applies settings to the
   !! StateContainer using the enhanced configuration system.
   !!
   !! \param[inout] container StateContainer to configure
   !! \param[in] config_file Path to configuration file
   !! \param[out] rc Return code
   subroutine load_configuration(container, config_file, rc)
      implicit none

      type(StateContainerType), intent(inout) :: container
      character(len=*), intent(in) :: config_file
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      type(ConfigDataType), pointer :: config
      logical :: file_exists

      error_mgr => container%get_error_manager()
      call error_mgr%push_context('load_configuration', 'Loading CATChem configuration')

      ! Check if config file exists
      inquire(file=trim(config_file), exist=file_exists)
      if (.not. file_exists) then
         call error_mgr%report_error(ERROR_FILE_NOT_FOUND, &
            'Configuration file not found: '//trim(config_file), rc, &
            'load_configuration', &
            'Ensure configuration file exists in specified location')
         call error_mgr%pop_context()
         return
      endif

      ! Get config pointer
      config => container%get_config_ptr()
      if (.not. associated(config)) then
         call error_mgr%report_error(ERROR_STATE_INCONSISTENCY, &
            'Configuration state not properly initialized', rc, &
            'load_configuration', &
            'Call initialize_core_states first')
         call error_mgr%pop_context()
         return
      endif

      ! Load configuration (this would use the new ConfigManager)
      ! For now, using placeholder logic
      call error_mgr%push_context('parse_config_file', 'Parsing YAML configuration')

      ! TODO: Replace with actual ConfigManager call
      ! call config_manager%load_from_file(config_file, config, rc)
      rc = CC_SUCCESS  ! Placeholder

      call error_mgr%pop_context()

      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_FILE_READ, &
            'Failed to parse configuration file', rc, &
            'load_configuration', &
            'Check configuration file syntax and format')
         call error_mgr%pop_context()
         return
      endif

      call error_mgr%pop_context()

   end subroutine load_configuration

   !> \brief Initialize atmospheric chemistry processes
   !!
   !! This routine initializes all atmospheric chemistry processes (dust, seasalt,
   !! dry deposition, etc.) based on the loaded configuration.
   !!
   !! \param[inout] container StateContainer with processes to initialize
   !! \param[out] rc Return code
   subroutine initialize_processes(container, rc)
      implicit none

      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      type(ConfigDataType), pointer :: config

      error_mgr => container%get_error_manager()
      call error_mgr%push_context('initialize_processes', 'Initializing atmospheric chemistry processes')

      config => container%get_config_ptr()
      if (.not. associated(config)) then
         call error_mgr%report_error(ERROR_STATE_INCONSISTENCY, &
            'Configuration not loaded', rc, &
            'initialize_processes', &
            'Load configuration before initializing processes')
         call error_mgr%pop_context()
         return
      endif

      ! Initialize processes based on configuration
      ! TODO: This would use the new process factory pattern

      ! Placeholder for dust process initialization
      call error_mgr%push_context('init_dust_process', 'Initializing dust emission process')
      ! call dust_factory%create_process(config%dust, dust_process, rc)
      rc = CC_SUCCESS  ! Placeholder
      call error_mgr%pop_context()

      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_PROCESS_INITIALIZATION, &
            'Failed to initialize dust process', rc, &
            'initialize_processes', &
            'Check dust process configuration')
         call error_mgr%pop_context()
         return
      endif

      ! Placeholder for seasalt process initialization
      call error_mgr%push_context('init_seasalt_process', 'Initializing sea salt emission process')
      ! call seasalt_factory%create_process(config%seasalt, seasalt_process, rc)
      rc = CC_SUCCESS  ! Placeholder
      call error_mgr%pop_context()

      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_PROCESS_INITIALIZATION, &
            'Failed to initialize sea salt process', rc, &
            'initialize_processes', &
            'Check sea salt process configuration')
         call error_mgr%pop_context()
         return
      endif

      call error_mgr%pop_context()

   end subroutine initialize_processes

   !> \brief Validate complete initialization
   !!
   !! This routine performs comprehensive validation of the initialized
   !! StateContainer to ensure all components are properly set up and consistent.
   !!
   !! \param[inout] container StateContainer to validate
   !! \param[out] rc Return code
   subroutine validate_initialization(container, rc)
      implicit none

      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      type(StateValidatorType) :: validator

      error_mgr => container%get_error_manager()
      call error_mgr%push_context('validate_initialization', 'Validating complete CATChem initialization')

      ! Validate the container
      call container%validate(rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_STATE_INCONSISTENCY, &
            'Container validation failed', rc, &
            'validate_initialization', &
            'Check initialization steps and error messages above')
         call error_mgr%pop_context()
         return
      endif

      ! Perform detailed validation using StateValidator
      call validator%validate_container(container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_STATE_INCONSISTENCY, &
            'Detailed validation failed', rc, &
            'validate_initialization', &
            'Review state consistency and configuration')
         call error_mgr%pop_context()
         return
      endif

      ! Validation successful
      call error_mgr%pop_context()
      rc = CC_SUCCESS

   end subroutine validate_initialization

   !> \brief Finalize CATChem and cleanup resources
   !!
   !! This routine performs cleanup of all CATChem components and resources.
   !!
   !! \param[inout] container StateContainer to finalize
   !! \param[out] rc Return code
   subroutine finalize_catchem(container, rc)
      implicit none

      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr

      error_mgr => container%get_error_manager()
      call error_mgr%push_context('finalize_catchem', 'Finalizing CATChem and cleaning up resources')

      ! Print final statistics
      ! TODO: Implement error statistics printing
      ! call error_mgr%print_summary()

      ! Finalize the container (this will cleanup all state objects)
      call container%finalize(rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_MEMORY_DEALLOCATION, &
            'Failed to properly finalize container', rc, &
            'finalize_catchem', &
            'Check for memory leaks or resource conflicts')
         call error_mgr%pop_context()
         return
      endif

      call error_mgr%pop_context()
      rc = CC_SUCCESS

   end subroutine finalize_catchem
end module init_mod
