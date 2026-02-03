!> \file catchem_nuopc_interface.F90
!! \brief NUOPC interface for CATChem atmospheric chemistry model
!!
!! \details
!! This module provides the core interface routines for the CATChem NUOPC cap,
!! handling data transformations between ESMF/NUOPC fields and CATChem states.
!! It follows similar patterns as the CCPP interface but is adapted for NUOPC
!! framework requirements including ESMF grids, fields, and parallel operations.
!!
!! Key functionalities include:
!! - CATChem model initialization and finalization within NUOPC framework
!! - Data transformation between ESMF fields and CATChem state objects
!! - Field mapping configuration and management
!! - Integration with CF-compliant input and NetCDF output systems
!! - Support for flexible grid configurations and parallel decomposition
!! - Error handling and logging for NUOPC applications
!!
!! The interface supports both sequential and parallel execution modes
!! and provides standardized data exchange capabilities for coupling
!! with other Earth system model components.
!!
!! \author Barry Baker, NOAA/OAR/ARL
!! \date November 2024
!! \ingroup catchem_nuopc_group

module catchem_nuopc_interface

   use ESMF
   use NUOPC
   use MPI
   use CATChem_API, only: CATChem_Model
   ! use catchem_nuopc_cf_input
   ! use catchem_nuopc_netcdf_out
   ! use machine, only: kind_phys
   use precision_mod, only: fp
   use Constants, only: g0, Rd
   use Error_Mod, only : CC_SUCCESS, CC_FAILURE
   use StateManager_Mod, only: StateManagerType
   use ProcessManager_Mod, only: ProcessManagerType
   use ConfigManager_Mod, only: ConfigManagerType
   use error_mod, only: ErrorManagerType
   use MetState_Mod, only: MetStateType
   use ChemState_Mod, only: ChemStateType
   use TimeState_Mod, only: TimeStateType
   use ExtEmisData_Mod, only: ExtEmisDataType  ! External emissions data type
   use DiagnosticManager_Mod, only: DiagnosticManagerType
   use DiagnosticInterface_Mod, only: DiagnosticRegistryType, DIAG_REAL_SCALAR, DIAG_REAL_1D, DIAG_REAL_2D, DIAG_REAL_3D
   use aqmio, only: AQMIO_Create, AQMIO_Write, AQMIO_Close, AQMIO_Write1D, AQMIO_FMT_NETCDF
   use catchem_emis_mod

   implicit none

   private

   public :: catchem_nuopc_init
   public :: catchem_nuopc_run
   public :: catchem_nuopc_finalize
   public :: transform_nuopc_to_catchem
   public :: transform_catchem_to_nuopc
   public :: load_field_config
   public :: catchem_diagnostics_write
   !public :: get_cc_wrap  ! Accessor for process-local wrapper
   public :: get_n_import_fields, get_import_field_info  ! Safe field_config access
   public :: get_n_export_fields, get_export_field_info  ! Safe field_config access

   !> \brief Field mapping configuration structure
   !!
   !! Defines the mapping between NUOPC standard names and CATChem variables,
   !! including metadata for proper data transformation and validation.
   !! \{
   type :: field_mapping_type
      character(len=128) :: standard_name !< NUOPC/CF standard field name
      character(len=128) :: catchem_var   !< Corresponding CATChem variable path
      integer :: dimensions               !< Number of spatial dimensions (2D/3D)
      character(len=64) :: units          !< Physical units for conversion
      logical :: optional = .false.       !< Whether field is required or optional
   end type field_mapping_type
   !! \}

   !> \brief Field configuration structure for NUOPC interface
   !!
   !! Contains the complete field mapping configuration including both
   !! import and export field definitions with associated metadata.
   !! \{
   type :: field_config_type
      integer :: n_import_fields = 0                         !< Number of import fields
      integer :: n_export_fields = 0                         !< Number of export fields
      type(field_mapping_type), allocatable :: import_fields(:) !< Import field mapping array
      type(field_mapping_type), allocatable :: export_fields(:) !< Export field mapping array
   end type field_config_type

   !> field_config can be shared among MPIs
   type(field_config_type), public :: field_config

   !> \brief tracer mapping between CATChem and NUOPC
   !!
   !! Defines the tracer mapping index between NUOPC and chem_state if CATChem
   !! \{
   type :: tracer_index_map
      integer, allocatable :: nuopc_to_cc(:)  !< mapping index from NUOPC to CATChem
      character(len=128), allocatable :: names(:) !< NUOPC tracer name
      character(len=128), allocatable :: units(:) !< NUOPC tracer unit
   end type tracer_index_map
   !! \}

   !> Container for process-private CATChem state to avoid MPI sharing
   type :: cc_wrap_type
      type(CATChem_Model) :: catchem_model
      type(ExtEmisDataType) :: ext_emis ! External emissions data object
      type(field_config_type) :: field_config  ! Moved from module level for MPI safety
      type(tracer_index_map) :: tracer_map
      type(ESMF_Grid) :: grid
      logical :: initialized = .false.
      ! Diagnostic output variables (moved from module level for MPI safety)
      type(ESMF_Time) :: last_output_time
      type(ESMF_Time) :: startTime
      type(ESMF_Time) :: endTime
      type(ESMF_TimeInterval) :: output_interval
      type(ESMF_TimeInterval) :: timeStep
      logical :: output_timing_initialized = .false.
      character(len=256) :: output_directory = './output'
      character(len=64) :: output_prefix = 'catchem_diag'
      integer :: output_frequency = 3600  ! Default: 1 hour in seconds
      type(ESMF_GridComp) :: iocomp
      ! Time slice tracking for NetCDF output
      integer :: current_time_slice = 0
   end type cc_wrap_type

   type CATChem_InternalState
      type(cc_wrap_type), pointer :: wrap => null()
   end type CATChem_InternalState
   public :: cc_wrap_type, CATChem_InternalState


contains

   !> Initialize CATChem model for NUOPC interface
   !!
   !! This routine performs comprehensive initialization of the CATChem model
   !! within the NUOPC framework, including configuration loading, memory
   !! allocation, and setup of I/O systems for external data and diagnostics.
   !!
   !! @param config CATChem configuration object to initialize
   !! @param catchem_states Main container for all CATChem state variables
   !! @param dustState Dust emission and transport state object
   !! @param seaSaltState Sea salt emission and transport state object
   !! @param dryDepState Dry deposition process state object
   !! @param im Number of horizontal grid points
   !! @param config_file Path to CATChem configuration file
   !! @param grid ESMF grid for regridding and I/O operations
   !! @param errflg Error flag (CC_SUCCESS on success)
   !! @param errmsg Error message string if errflg indicates failure
   !!
   !! This routine performs the following initialization steps:
   !! - Reads and validates the CATChem configuration file
   !! - Initializes all CATChem state objects and process modules
   !! - Allocates memory for MetState, EmisState, ChemState, and DiagState arrays
   !! - Sets up field mapping configuration for NUOPC data exchange
   !! - Initializes CF-compliant input system for external data
   !! - Initializes NetCDF output system for diagnostic data
   !! - Validates grid compatibility and spatial dimensions
   !!
   !! The routine is similar to catchem_init in the CCPP interface but includes
   !! additional NUOPC-specific functionality such as ESMF grid integration
   !! and YAML-based field configuration management.
   !!
   !! @note This routine must be called before any CATChem calculations
   !!       and requires valid ESMF grid and configuration files
   !!
   !! @warning Proper error checking should be performed on errflg after calling
   subroutine catchem_nuopc_init(model, config_file, lat, lon, nlev, tracerinfo, input_grid, startTime,stopTime, timeStep, clock, nsoil, nsoiltype, nsurftype, rc)
      use ChemSpeciesUtils_Mod, only : create_species_mapping

      type(ESMF_GridComp)  :: model
      character(len=*), intent(in) :: config_file
      real(ESMF_KIND_R8), dimension(:,:), intent(in) :: lat
      real(ESMF_KIND_R8), dimension(:,:), intent(in) :: lon
      integer, intent(in) :: nlev
      type(ESMF_Info), intent(in) :: tracerinfo
      type(ESMF_Grid), intent(in) :: input_grid
      type(ESMF_Time), intent(in), optional :: startTime,stopTime
      type(ESMF_TimeInterval), intent(in), optional :: timeStep
      type(ESMF_Clock), intent(in), optional :: clock
      integer, intent(in), optional :: nsoil, nsoiltype, nsurftype
      integer, intent(out) :: rc

      ! Local variables
      type(StateManagerType), pointer :: state_mgr
      type(ConfigManagerType), pointer :: config_manager
      type(MetStateType), pointer :: met_state
      type(ChemStateType), pointer :: chem_state
      integer :: nx, ny, num_processes, stat
      integer(ESMF_KIND_I8) :: tstep_seconds
      character(len=128), allocatable :: tracer_names(:) !< NUOPC tracer name
      character(len=128), allocatable :: tracer_units(:) !< NUOPC tracer unit
      type(CATChem_InternalState) :: is
      type(cc_wrap_type), pointer:: cc_wrap

      ! Initialize
      rc = CC_SUCCESS

      ! -- allocate memory for the internal state and store it into component
      allocate(is%wrap, stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return !bail out

      cc_wrap => is%wrap

      !get nx, ny
      nx = size(lat, 1)
      ny = size(lat, 2)

      ! Initialize catchem using process-local variable
      if (present(nsoil) .and. present(nsoiltype) .and. present(nsurftype)) then
         call cc_wrap%catchem_model%initialize(config_file, nx, ny, nlev, nsoil, nsoiltype, nsurftype, rc)
      else
         call cc_wrap%catchem_model%initialize(config_file, nx, ny, nlev, rc = rc)
      end if

      if (rc /= CC_SUCCESS) then
         call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
            msg="CATChem initialization failed", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return  ! bail out
      end if

      !assign lat and lon to metstate
      state_mgr => cc_wrap%catchem_model%get_state_manager()
      met_state => state_mgr%get_met_state_ptr()
      met_state%lat = lat
      met_state%lon = lon
      ! Convert longitude from 0–360 to -180–180
      where (met_state%lon > 180.0_fp)
         met_state%lon = met_state%lon - 360.0_fp
      end where

      !initialize extemission data here
      config_manager => state_mgr%get_config_ptr()
      call catchem_emis_init(cc_wrap%ext_emis, config_manager, nx, ny, nlev, clock, rc)

      !populate tracer mapping using process-local tracer_map
      call TracerInfoGet(tracerinfo, 'tracerNames', tracer_names, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__,  file=__FILE__)) return  ! bail out

      if (.not.allocated(tracer_names)) then
         call ESMF_LogWrite("Unable to retrieve imported tracer list", &
            ESMF_LOGMSG_WARNING, line=__LINE__, file=__FILE__, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         return
      end if

      ! - import tracer units if available
      call TracerInfoGet(tracerinfo, 'tracerUnits', tracer_units, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__,  file=__FILE__)) return  ! bail out

      if (.not.allocated(tracer_units)) then
         allocate(tracer_units(size(tracer_names)), stat=stat)
         if (ESMF_LogFoundAllocError(statusToCheck=stat, &
            msg="Unable to allocate internal workspace", &
            line=__LINE__,  file=__FILE__)) return  ! bail out
         tracer_units = 'n/a'
      end if

      !copy to cc_wrap
      allocate(cc_wrap%tracer_map%names, stat=stat, source=tracer_names)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
         msg="Unable to allocate tracers name", &
         line=__LINE__,  file=__FILE__, rcToReturn=rc)) return  ! bail out

      allocate(cc_wrap%tracer_map%units, stat=stat, source=tracer_units)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
         msg="Unable to allocate tracers unit", &
         line=__LINE__,  file=__FILE__, rcToReturn=rc)) return  ! bail out

      allocate(cc_wrap%tracer_map%nuopc_to_cc(size(tracer_names)), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
         msg="Unable to allocate nuopc_to_cc mapping", &
         line=__LINE__,  file=__FILE__, rcToReturn=rc)) return  ! bail out

      ! assign mapping index
      call create_species_mapping(state_mgr, cc_wrap%tracer_map%names, cc_wrap%tracer_map%nuopc_to_cc, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__,  file=__FILE__)) return  ! bail out

      !copy fields to cc_wrap
      cc_wrap%field_config = field_config
      ! Set the process-local grid variable
      cc_wrap%grid = input_grid
      ! Set time information if provided
      if (present(stopTime)) then
         cc_wrap%endTime = stopTime
      end if
      if (present(startTime)) then
         cc_wrap%startTime = startTime
      end if

      if (present(timeStep)) then
         cc_wrap%timeStep = timeStep
         call ESMF_TimeIntervalGet(timeStep, s_i8=tstep_seconds, rc=rc)
         state_mgr%tstep = real(tstep_seconds, fp)
      end if

      ! Add all enabled processes from configuration
      call cc_wrap%catchem_model%add_process(rc)
      num_processes = cc_wrap%catchem_model%get_num_processes()
      if (rc /= CC_SUCCESS .or. num_processes <= 0) then
         call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
            msg="CATChem initialization failed", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return  ! bail out
      end if

      ! Mark this process as initialized
      cc_wrap%initialized = .true.

      !store cc_wrap into component
      call ESMF_GridCompSetInternalState(model, is, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__,  file=__FILE__)) return  ! bail out

      !deallocate variable
      if (allocated(tracer_names)) deallocate(tracer_names)
      if (allocated(tracer_units)) deallocate(tracer_units)
      if (allocated(field_config%import_fields)) deallocate(field_config%import_fields)
      if (allocated(field_config%export_fields)) deallocate(field_config%export_fields)

      ! ! Initialize CF input system
      ! call cf_input_init('catchem_input_config.yml', grid, errflg)
      ! if (errflg /= ESMF_SUCCESS) then
      !   errmsg = 'Error initializing CF input system'
      !   errflg = CC_FAILURE
      !   return
      ! end if

      ! ! Initialize NetCDF output system
      ! call output_diagnostics_init('catchem_output_config.yml', grid, errflg)
      ! if (errflg /= ESMF_SUCCESS) then
      !   errmsg = 'Error initializing NetCDF output system'
      !   errflg = CC_FAILURE
      !   return
      ! end if

   end subroutine catchem_nuopc_init

   !> Get process-local CATChem wrapper (guaranteed thread/process safe)
   !!
   !! This function provides access to the process-local CATChem state
   !! using static local variables which are guaranteed to be process-private
   !! in MPI environments.
   !!
   !! \return Pointer to process-local CATChem wrapper
   ! function get_cc_wrap() result(wrap_ptr)
   !   type(ESMF_GridComp)         :: model
   !   type(cc_wrap_type), pointer :: wrap_ptr

   !   type(CATChem_InternalState), target :: is
   !   integer                    :: verbosity, localrc
   !   character(len=128) :: name

   !   ! -- get component's information
   !   call NUOPC_CompGet(model, rc=localrc)
   !   if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
   !     line=__LINE__,  file=__FILE__)) return  ! bail out

   !   ! -- get component's internal state
   !   call ESMF_GridCompGetInternalState(model, is, localrc)
   !   if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
   !     line=__LINE__,  file=__FILE__))  return  ! bail out

   !   wrap_ptr => is%wrap
   ! end function get_cc_wrap

   ! Run CATChem processes for NUOPC interface
   !!
   !! \param config         CATChem configuration
   !! \param catchem_states CATChem container
   !! \param dustState      Dust process state
   !! \param seaSaltState   Sea salt process state
   !! \param dryDepState    Dry deposition process state
   !! \param    dt             Time step (seconds)
   !! \param    current_time   Current model time
   !! \param   errflg         Error flag
   !! \param   errmsg         Error message
   !!
   subroutine catchem_nuopc_run( cc_wrap, dt, current_time, errmsg, rc)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      real(ESMF_KIND_R8), intent(in) :: dt
      type(ESMF_Time), intent(in) :: current_time
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: rc

      ! Get process-local state
      type(StateManagerType), pointer :: state_mgr => null()
      integer, save :: timestep = 0

      !cc_wrap => get_cc_wrap()

      rc = CC_SUCCESS
      errmsg = ''

      ! Update extemission data first
      state_mgr => cc_wrap%catchem_model%get_state_manager()
      call catchem_emis_update(cc_wrap%ext_emis, current_time, state_mgr, &
         cc_wrap%iocomp, cc_wrap%grid, real(dt, fp), rc)

      !Run CATChem processes
      timestep = timestep + 1
      call cc_wrap%catchem_model%run_timestep(timestep, real(dt, fp), rc)
      if (rc /= CC_SUCCESS) then
         write(errmsg, '(A,I0)') 'Error in run_timestep at timestep = ', timestep
         return
      end if

      ! Write NetCDF output diagnostics if needed
      call catchem_diagnostics_write(cc_wrap, current_time, rc)
      if (rc /= ESMF_SUCCESS) then
         errmsg = 'Error writing NetCDF output diagnostics'
         return
      end if

   end subroutine catchem_nuopc_run

   ! Finalize CATChem for NUOPC interface
   !!
   !! \param catchem_states CATChem container
   !! \param dustState      Dust process state
   !! \param seaSaltState   Sea salt process state
   !! \param dryDepState    Dry deposition process state
   !! \param   errflg         Error flag
   !! \param   errmsg         Error message
   !!
   subroutine catchem_nuopc_finalize(cc_wrap, rc, errmsg)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      integer, intent(out) :: rc
      character(len=*), intent(out) :: errmsg

      ! Get process-local state
      !type(cc_wrap_type), pointer :: cc_wrap
      !cc_wrap => get_cc_wrap()

      rc = CC_SUCCESS
      errmsg = ''

      ! Finalize CF input system
      !call cf_input_finalize()

      ! Finalize NetCDF output system
      !call output_diagnostics_finalize()

      ! Finalize CATChem model
      if (cc_wrap%initialized) then
         !finalize extemission data
         call catchem_emis_finalize(cc_wrap%ext_emis, rc)
         !finalize catchem model
         call cc_wrap%catchem_model%finalize(rc)
         if (rc /= CC_SUCCESS) then
            errmsg = 'Error in calling cc_wrap%catchem_model%finalize!'
            return
         end if
         cc_wrap%initialized = .false.
      end if

      ! Deallocate field mappings
      if (allocated(cc_wrap%field_config%import_fields)) deallocate(cc_wrap%field_config%import_fields)
      if (allocated(cc_wrap%field_config%export_fields)) deallocate(cc_wrap%field_config%export_fields)

   end subroutine catchem_nuopc_finalize

   ! Transform NUOPC import fields to CATChem states
   !!
   !! \param    importState    NUOPC import state
   !! \param catchem_states CATChem container
   !! \param    im             Horizontal dimension
   !! \param    kme            Vertical dimension
   !! \param   rc             ESMF return code
   !!
   subroutine transform_nuopc_to_catchem(cc_wrap, importState, currTime, rc)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_State), intent(in) :: importState
      type(ESMF_Time), intent(in) :: currTime
      integer, intent(out) :: rc

      type(ESMF_Field) :: field
      type(StateManagerType), pointer :: state_mgr
      type(ErrorManagerType), pointer :: error_mgr
      type(TimeStateType), pointer :: time_state
      type(MetStateType), pointer :: met_state
      logical, allocatable :: set_required_met(:)
      integer(ESMF_KIND_I8) :: timestep_seconds
      integer :: year, month, day, hour, minute, second
      integer :: i, n, n_met
      !type(cc_wrap_type), pointer :: cc_wrap

      rc = ESMF_SUCCESS

      ! assign time to catchem model's time state
      state_mgr => cc_wrap%catchem_model%get_state_manager()
      error_mgr => state_mgr%get_error_manager()
      time_state => state_mgr%get_time_state_ptr()
      met_state => state_mgr%get_met_state_ptr()

      call ESMF_TimeGet(currTime, yy=year, mm=month, dd=day, &
         h=hour, m=minute, s=second, rc=rc)
      call ESMF_TimeIntervalGet(cc_wrap%timeStep, s_i8=timestep_seconds, rc=rc)
      call time_state%init(year, month, day, hour, minute, second, real(timestep_seconds), error_mgr, rc)
      if (rc /= CC_SUCCESS) then
         return !maybe add an error message
      end if

      ! This is to check if all required met fields in CATChem are set
      if (allocated(cc_wrap%catchem_model%required_fields)) then
         n_met = size(cc_wrap%catchem_model%required_fields)
         allocate(set_required_met(n_met))
         set_required_met = .false.
      end if

      ! Loop through all import fields and transform to CATChem states
      do n = 1, cc_wrap%field_config%n_import_fields

         ! Try to get field from import state (will fail if not present)
         call ESMF_StateGet(importState, trim(cc_wrap%field_config%import_fields(n)%standard_name), field, rc=rc)

         if (rc /= ESMF_SUCCESS) then
            if (.not. cc_wrap%field_config%import_fields(n)%optional) then
               call ESMF_LogWrite("Required field not found: "// &
                  trim(cc_wrap%field_config%import_fields(n)%standard_name), ESMF_LOGMSG_ERROR, rc=rc)
               rc = ESMF_FAILURE
               return
            else
               cycle  ! Skip optional fields that are not present
            end if
         end if

         ! Transform based on field type and dimensions
         call transform_field_to_catchem(cc_wrap, field, cc_wrap%field_config%import_fields(n), &
            cc_wrap%field_config%import_fields(n)%optional, set_required_met, rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

      end do

      !derive some met fields if required after reading from NUOPC
      if (allocated(cc_wrap%catchem_model%required_fields)) then
         do i = 1, n_met
            if (.not. set_required_met(i)) then
               call met_state%derive_field(trim(cc_wrap%catchem_model%required_fields(i)), error_mgr, time_state, rc)
               if (rc /= CC_SUCCESS) then
                  call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
                     msg="Error deriving required met field: "// trim(cc_wrap%catchem_model%required_fields(i)), &
                     line=__LINE__, file=__FILE__, rcToReturn=rc)
                  return  ! bail out
               else
                  set_required_met(i) = .true.
               end if
            end if
         end do
      end if

      !check if all require met fields are set
      if (allocated(cc_wrap%catchem_model%required_fields)) then
         do i = 1, n_met
            if (.not. set_required_met(i)) then
               !write(*,*) 'Wait. A required field is not set: ' // trim(cc_wrap%catchem_model%required_fields(i))
               call ESMF_LogWrite("Required met field not set yet: "// &
                  trim(cc_wrap%catchem_model%required_fields(i)), ESMF_LOGMSG_ERROR, rc=rc)
               rc = ESMF_FAILURE
               return
            end if
         end do
         !deallocate array
         deallocate(set_required_met)
      end if

   end subroutine transform_nuopc_to_catchem

   ! Transform CATChem states to NUOPC export fields
   !!
   !! \param exportState    NUOPC export state
   !! \param    catchem_states CATChem container
   !! \param    im             Horizontal dimension
   !! \param    kme            Vertical dimension
   !! \param   rc             ESMF return code
   !!
   subroutine transform_catchem_to_nuopc(cc_wrap, exportState, rc)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_State), intent(inout) :: exportState
      integer, intent(out) :: rc

      type(ESMF_Field) :: field
      integer :: n
      !type(cc_wrap_type), pointer :: cc_wrap

      rc = ESMF_SUCCESS

      ! Get process-local state
      !cc_wrap => get_cc_wrap()

      ! Loop through all export fields and transform from CATChem states
      do n = 1, cc_wrap%field_config%n_export_fields

         ! Try to get field from export state (will fail if not present)
         call ESMF_StateGet(exportState, trim(cc_wrap%field_config%export_fields(n)%standard_name), field, rc=rc)

         if (rc /= ESMF_SUCCESS) then
            if (.not. cc_wrap%field_config%export_fields(n)%optional) then
               call ESMF_LogWrite("Required export field not found: "// &
                  trim(cc_wrap%field_config%export_fields(n)%standard_name), ESMF_LOGMSG_ERROR, rc=rc)
               rc = ESMF_FAILURE
               return
            else
               cycle  ! Skip optional fields that are not present
            end if
         end if

         ! Transform from CATChem to field
         call transform_catchem_to_field(cc_wrap, field, cc_wrap%field_config%export_fields(n), rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

      end do

   end subroutine transform_catchem_to_nuopc

   ! Transform individual field to CATChem state
   !!
   !! \param    field          ESMF field
   !! \param    field_map      Field mapping information
   !! \param catchem_states CATChem container
   !! \param    im             Horizontal dimension
   !! \param    kme            Vertical dimension
   !! \param   rc             ESMF return code
   !!
   subroutine transform_field_to_catchem(cc_wrap, field, field_map, required, is_met_set, rc)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_Field), intent(in) :: field
      type(field_mapping_type), intent(in) :: field_map
      logical, dimension(:), intent(inout) :: is_met_set
      logical, intent(in) :: required
      integer, intent(out) :: rc

      !local vars
      !type(ProcessManagerType), pointer :: process_mgr
      type(StateManagerType), pointer :: state_mgr
      type(ErrorManagerType), pointer :: error_mgr
      type(MetStateType), pointer :: met_state
      type(ChemStateType), pointer :: chem_state
      !type(cc_wrap_type), pointer :: cc_wrap
      real(ESMF_KIND_R8), pointer :: fptr4d(:,:,:,:), fptr3d(:,:,:), fptr2d(:,:)
      real(ESMF_KIND_R8), pointer :: fptr4d_rev(:,:,:,:), fptr3d_rev(:,:,:)
      real(fp), pointer :: column_ptr(:) !catchem met column pointer to get vertical dimension for nz+1 variables
      real(ESMF_KIND_R8) :: unit_conv
      integer :: i, j, k, v, ni, nj, nk, nk1, nv, kk, v_cc, met_index

      rc = ESMF_SUCCESS

      ! Get process-local state
      !process_mgr => cc_wrap%catchem_model%get_process_manager()
      state_mgr => cc_wrap%catchem_model%get_state_manager()
      error_mgr => state_mgr%get_error_manager()
      met_state => state_mgr%get_met_state_ptr()

      if ( .not. associated(met_state) ) then
         call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
            msg="met_state is not associated in CATChem before transformation from NUOPC", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return  ! bail out
      end if

      ! Transform based on field mapping
      !select case (trim(field_map%catchem_var))
      select case (field_map%dimensions)

         ! 2D meteorological fields
       case (2)
         nullify(fptr2d)
         call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         !set to met_state in CATChem
         if (trim(field_map%catchem_var) == 'DLUSE' .or. trim(field_map%catchem_var) == 'DSOILTYPE' .or. &
            trim(field_map%catchem_var) == 'LWI') then
            !convert to integer
            call met_state%set_field(trim(field_map%catchem_var), int(fptr2d), error_mgr, rc)
         else if (trim(field_map%catchem_var) == 'Z0') then ! roughness length in cm in NUOPC but m in CATChem
            call met_state%set_field(trim(field_map%catchem_var), real(fptr2d, fp)*0.01_fp, error_mgr, rc)
         else
            call met_state%set_field(trim(field_map%catchem_var), real(fptr2d, fp), error_mgr, rc)
         end if

         if (rc == CC_SUCCESS) then
            if (allocated(cc_wrap%catchem_model%required_fields)) then
               met_index = cc_wrap%catchem_model%get_required_met_index( trim(field_map%catchem_var) )
               if (met_index >0 ) then
                  is_met_set(met_index) = .true.
               end if
            end if
            !TODO: met%set_met will stop the model run if field not matching. We may fix it later.
         else if (.not. required) then
            ! If the field is not required, we can skip the transformation
            call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
               msg="Met field is not set and its optional: " // trim(field_map%catchem_var), &
               line=__LINE__, file=__FILE__, rcToReturn=rc)
         else
            call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
               msg="Met field is not set successfully for: " // trim(field_map%catchem_var), &
               line=__LINE__, file=__FILE__, rcToReturn=rc)
            return  ! bail out
         end if

         ! 3D meteorological fields
       case (3)
         nullify(fptr3d, fptr3d_rev)
         call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         ni = size(fptr3d, 1)
         nj = size(fptr3d, 2)
         nk = size(fptr3d, 3)

         if (trim(field_map%catchem_var) .ne. 'SOILM' .and.  trim(field_map%catchem_var) .ne. 'SOILT') then
            !get catchem receriver vertical dimension for nz+1 variables while NUOPC has nz levels
            !Currently only PFILSAN and PFLLSAN are in this case following GOCART and in most cases,
            ! nk == nk1
            call met_state%get_field_ptr(trim(field_map%catchem_var), i=1, j=1, col_ptr=column_ptr, rc=rc)
            if (rc /= CC_SUCCESS) then
               call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
                  msg="Error getting met field pointer for: " // trim(field_map%catchem_var), &
                  line=__LINE__, file=__FILE__, rcToReturn=rc)
               return  ! bail out
            end if
            nk1 = size(column_ptr)
         else
            !SOILM and SOILT have not been allocated because we do not have soil layers yet and the get_field_ptr will fail
            nk1 = nk
         end if

         ! Allocate fptr3d_rev with the same dimensions as fptr3d
         allocate(fptr3d_rev(ni, nj, nk1))

         ! -- map provider field levels to receiver field levels in the same (not reverse) order
         ! -- NOTE: if provider field from NUOPC has fewer vertical levels than the receiver field in CATChem,
         ! -- the remaining receiver field levels are filled by replicating values from
         ! -- the closest available level in the provider field.
         kk = 1
         do k = 1, nk1
            !kk = nk - k + 1 !no need to reverse
            !kk = k
            do j = 1, nj
               do i = 1, ni
                  if (trim(field_map%catchem_var) == 'Z' .or. trim(field_map%catchem_var) == 'ZMID') then
                     fptr3d_rev(i,j,k) = fptr3d(i,j,kk) / g0
                  else
                     fptr3d_rev(i,j,k) = fptr3d(i,j,kk)
                  end if
               end do
            end do
            kk = min(nk, kk + 1)
         end do

         !set to met_state in CATChem
         call met_state%set_field(trim(field_map%catchem_var), real(fptr3d_rev,fp), error_mgr, rc)

         if (rc == CC_SUCCESS) then
            if (allocated(cc_wrap%catchem_model%required_fields)) then
               met_index = cc_wrap%catchem_model%get_required_met_index( trim(field_map%catchem_var) )
               if (met_index >0 ) then
                  is_met_set(met_index) = .true.
               end if
            end if
         else if (.not. required) then
            ! If the field is not required, we can skip the transformation
            call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
               msg="Met field is not set and its optional: " // trim(field_map%catchem_var), &
               line=__LINE__, file=__FILE__, rcToReturn=rc)
         else
            call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
               msg="Met field is not set successfully for: " // trim(field_map%catchem_var), &
               line=__LINE__, file=__FILE__, rcToReturn=rc)
            deallocate(fptr3d_rev)  ! Clean up before returning
            return  ! bail out
         end if

         ! Clean up allocated memory
         deallocate(fptr3d_rev)

         ! 4D tracer concentrations
       case (4)
         nullify(fptr4d, fptr4d_rev)
         call ESMF_FieldGet(field, farrayPtr=fptr4d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         chem_state => state_mgr%get_chem_state_ptr()
         if ( .not. associated(chem_state) ) then
            call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
               msg="chem_state is not associated in CATChem before transformation from NUOPC", &
               line=__LINE__, file=__FILE__, rcToReturn=rc)
            return  ! bail out
         end if

         ni = size(fptr4d, 1)
         nj = size(fptr4d, 2)
         nk = size(fptr4d, 3)
         nv = size(fptr4d, 4)

         ! Allocate fptr4d_rev with the same dimensions as fptr4d
         allocate(fptr4d_rev(ni, nj, nk, size(chem_state%ChemSpecies)))
         fptr4d_rev = 0.0_fp  ! Initialize to zero

         ! Reverse vertical layers
         do v = 1, nv
            !read in specific humidity from tracer array
            if (trim(cc_wrap%tracer_map%names(v)) == 'sphum') then
               call met_state%set_field('QV', real(fptr4d(:,:, :,v), fp), error_mgr, rc)
               if (rc == CC_SUCCESS) then
                  if (allocated(cc_wrap%catchem_model%required_fields)) then
                     met_index = cc_wrap%catchem_model%get_required_met_index( 'QV' )
                     if (met_index >0 ) then
                        is_met_set(met_index) = .true.
                     end if
                  end if
               else if (.not. required) then
                  ! If the field is not required, we can skip the transformation
                  call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
                     msg="Met field is not set and its optional: QV", &
                     line=__LINE__, file=__FILE__, rcToReturn=rc)
               else
                  call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
                     msg="Met field is not set successfully for: QV", &
                     line=__LINE__, file=__FILE__, rcToReturn=rc)
                  deallocate(fptr4d_rev)  ! Clean up before returning
                  return  ! bail out
               end if
            end if

            !map NUOPC tracer index to CATChem species index
            v_cc = cc_wrap%tracer_map%nuopc_to_cc(v)
            if (v_cc <= 0) cycle !if not a species in CATChem, go to next cycle
            !unit conversion
            if (chem_state%ChemSpecies(v_cc)%is_gas) then
               !unit_conv = 28.9644  / chem_state%ChemSpecies(v_cc)%mw_g * 1.0e-3  ! convert from ug/kg to ppm for gases
               unit_conv = 1.00  !keep it in ppmV
            else
               unit_conv = 1.00  ! convert from ug/kg to ug/kg for aerosols
            end if

            do k = 1, nk
               !kk = nk - k + 1 !no need to reverse
               kk = k
               do j = 1, nj
                  do i = 1, ni
                     fptr4d_rev(i,j,kk,v_cc) = max(fptr4d(i,j,k,v), 0.0_fp) * unit_conv
                  end do
               end do
            end do
         end do

         !set to concentrations in CATChem
         call chem_state%set_all_concentrations(real(fptr4d_rev, fp), rc)
         if (rc /= CC_SUCCESS) then
            call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
               msg="Tracer array is not retrieved successfully for: " // trim(field_map%catchem_var), &
               line=__LINE__, file=__FILE__, rcToReturn=rc)
            deallocate(fptr4d_rev)  ! Clean up before returning
            return  ! bail out
         end if

         ! Clean up allocated memory
         deallocate(fptr4d_rev)

       case default
         call ESMF_LogWrite("Unknown field mapping dimension for: " // trim(field_map%catchem_var), &
            ESMF_LOGMSG_WARNING, rc=rc)

      end select

   end subroutine transform_field_to_catchem

   ! Transform CATChem state to ESMF field
   !!
   !! \param    catchem_states CATChem container
   !! \param    field_map      Field mapping information
   !! \param field          ESMF field
   !! \param    im             Horizontal dimension
   !! \param    kme            Vertical dimension
   !! \param   rc             ESMF return code
   !!
   subroutine transform_catchem_to_field(cc_wrap, field, field_map, rc)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(field_mapping_type), intent(in) :: field_map
      type(ESMF_Field), intent(inout) :: field
      integer, intent(out) :: rc

      type(StateManagerType), pointer :: state_mgr
      type(ChemStateType), pointer :: chem_state
      !type(cc_wrap_type), pointer :: cc_wrap
      real(ESMF_KIND_R8), pointer :: fptr4d(:,:,:,:), fptr3d(:,:,:), fptr2d(:,:)
      real(fp), allocatable :: cc_diag_data(:,:,:)
      character(len=128), allocatable :: diagnostic_names(:)
      real(ESMF_KIND_R8) :: unit_conv
      integer :: i, j, k, v, ni, nj, nk, kk, nv, v_cc, found_index

      rc = ESMF_SUCCESS

      ! Get process-local state
      !cc_wrap => get_cc_wrap()

      !TODO: we assume all the export fields are from DiagManager
      call cc_wrap%catchem_model%get_diagnostic_names(diagnostic_names, rc = rc)

      state_mgr => cc_wrap%catchem_model%get_state_manager()

      ! Transform based on field mapping
      select case (field_map%dimensions)

         ! 2D  fields
       case (2)
         nullify(fptr2d)
         if(allocated(cc_diag_data)) deallocate(cc_diag_data)
         ! get field pointer
         call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         !find index of diagnostic name in the format of process_name.field_name
         found_index = cc_wrap%catchem_model%get_diag_index_from_field(field_map%catchem_var)
         if (found_index > 0) then
            call cc_wrap%catchem_model%get_diagnostic(diagnostic_names(found_index), cc_diag_data, rc)
            if (rc /= ESMF_SUCCESS) then
               call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
                  msg="Failed to get diagnostic data for: " // trim(diagnostic_names(found_index)), &
                  line=__LINE__, file=__FILE__, rcToReturn=rc)
               return
            end if

            !assign data back to NUOPC
            fptr2d = cc_diag_data(:,:,1)

         else
            fptr2d = 0.0_ESMF_KIND_R8  ! Species not found
         end if

         ! 3D  fields
       case (3)
         nullify(fptr3d)
         if(allocated(cc_diag_data)) deallocate(cc_diag_data)
         ! get field pointer
         call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         !find index of diagnostic name in the format of process_name.field_name
         found_index = cc_wrap%catchem_model%get_diag_index_from_field(field_map%catchem_var)
         if (found_index > 0) then
            call cc_wrap%catchem_model%get_diagnostic(diagnostic_names(found_index), cc_diag_data, rc)
            if (rc /= ESMF_SUCCESS) then
               call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
                  msg="Failed to get diagnostic data for: " // trim(diagnostic_names(found_index)), &
                  line=__LINE__, file=__FILE__, rcToReturn=rc)
               return
            end if

            !assign data back to NUOPC
            ni = size(fptr3d, 1)
            nj = size(fptr3d, 2)
            nk = size(fptr3d, 3)
            !revserse vertical layers
            do k = 1, nk
               !kk = nk - k + 1 !no need to reverse
               kk = k
               do j = 1, nj
                  do i = 1, ni
                     fptr3d(i,j,kk) = cc_diag_data(i,j,k)
                  end do
               end do
            end do

         else
            fptr3d = 0.0_ESMF_KIND_R8  ! Species not found
         end if

         ! 4D  fields for chemical species
       case (4)
         nullify(fptr4d)
         if(allocated(cc_diag_data)) deallocate(cc_diag_data)
         ! get field pointer
         call ESMF_FieldGet(field, farrayPtr=fptr4d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         chem_state => state_mgr%get_chem_state_ptr()
         if ( .not. associated(chem_state) ) then
            call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
               msg="chem_state is not associated in CATChem before transformation to NUOPC", &
               line=__LINE__, file=__FILE__, rcToReturn=rc)
            return  ! bail out
         end if

         ni = size(fptr4d, 1)
         nj = size(fptr4d, 2)
         nk = size(fptr4d, 3)
         nv = size(fptr4d, 4)
         ! Reverse vertical layers
         do v = 1, nv
            v_cc = cc_wrap%tracer_map%nuopc_to_cc(v)
            if (v_cc > 0) then
               cc_diag_data = chem_state%ChemSpecies(v_cc)%conc
               if (chem_state%ChemSpecies(v_cc)%is_gas) then
                  !unit_conv = 1.0e3 * chem_state%ChemSpecies(v_cc)%mw_g /28.9644  ! convert from ppm to ug/kg for gases
                  unit_conv = 1.00  !keep it in ppmV
               else
                  unit_conv = 1.00  ! convert from ug/kg to ug/kg for aerosols
               end if

               do k = 1, nk
                  !kk = nk - k + 1 !no need to reverse
                  kk = k
                  do j = 1, nj
                     do i = 1, ni
                        fptr4d(i,j,kk,v) = cc_diag_data(i,j,k) * unit_conv
                     end do
                  end do
               end do
            else
               !fptr4d(:,:,:,v) = 0.0_ESMF_KIND_R8  ! Species not found; do nothing
            end if
         end do   !nv

       case default
         call ESMF_LogWrite("Unknown export field dimension for: "//trim(field_map%catchem_var), &
            ESMF_LOGMSG_WARNING, rc=rc)

      end select

   end subroutine transform_catchem_to_field

   !> \brief Write diagnostic output using AQMIO
   !!
   !! This subroutine writes CATChem diagnostic fields to NetCDF files using the
   !! working AQMIO module. It retrieves field data, metadata (description, units),
   !! and configuration from the DiagnosticManager and creates properly documented
   !! NetCDF output files.
   !!
   !! \param current_time Current simulation time
   !! \param rc Return code
   subroutine catchem_diagnostics_write(cc_wrap, current_time, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_Time), intent(in) :: current_time
      integer, intent(out) :: rc

      !type(cc_wrap_type), pointer :: cc_wrap
      type(DiagnosticManagerType), pointer :: diag_mgr => null()
      type(ESMF_Time) :: time_on_file
      character(len=64), allocatable :: process_list(:)
      character(len=64), allocatable :: field_names(:)
      integer :: num_processes, num_fields, i, j
      logical :: time_to_write
      character(len=256) :: filename
      character(len=*), parameter :: routine = 'catchem_diagnostics_write'

      rc = CC_SUCCESS

      ! Get process-local state
      !cc_wrap => get_cc_wrap()

      ! Initialize output timing if not done
      if (.not. cc_wrap%output_timing_initialized) then
         call initialize_output_timing(cc_wrap, current_time, rc)
         if (rc /= CC_SUCCESS) return
      end if

      ! Check if it's time to write output
      call check_diagnostic_output_time(cc_wrap, current_time, time_to_write, time_on_file, rc)
      if (rc /= CC_SUCCESS) return
      if (.not. time_to_write) return

      ! Get diagnostic manager
      diag_mgr => cc_wrap%catchem_model%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) then
         rc = CC_FAILURE
         return
      end if

      ! Get list of processes with diagnostics
      call diag_mgr%list_processes(process_list, num_processes, rc)
      if (rc /= CC_SUCCESS .or. num_processes == 0) return

      ! Use grid (must be set during initialization)
      if (.not. ESMF_GridIsCreated(cc_wrap%grid)) then
         rc = CC_FAILURE
         write(*,'(A)') 'Error: grid not initialized.'
         return
      end if

      ! Initialize AQMIO component if not done
      if (.not. ESMF_GridCompIsCreated(cc_wrap%iocomp)) then
         cc_wrap%iocomp = AQMIO_Create(cc_wrap%grid, rc =rc)
         if (rc /= CC_SUCCESS) return
      end if

      ! Generate filename for current time
      call generate_diagnostic_filename(cc_wrap, time_on_file, filename, rc)
      if (rc /= CC_SUCCESS) return

      ! Update time variable in NetCDF file and get the time slice to use
      call update_time_variable(cc_wrap, filename, time_on_file, cc_wrap%current_time_slice, rc)
      if (rc /= CC_SUCCESS) return

      ! Write diagnostics for each process
      do i = 1, num_processes
         call write_process_diagnostics(cc_wrap, trim(process_list(i)), filename, rc)
         if (rc /= CC_SUCCESS) then
            ! Log error and return
            write(*,'(A,A)') 'Error: Failed to write diagnostics for process: ', trim(process_list(i))
            return
         end if
      end do

      !write extemission fields if needed
      call catchem_emis_write_diagnostics(cc_wrap%ext_emis, cc_wrap%current_time_slice, cc_wrap%iocomp, cc_wrap%grid, filename, rc)
      if (rc /= CC_SUCCESS) then
         write(*,'(A)') 'Error: Failed to write external emission diagnostics.'
         return
      end if

      ! Update last output time
      cc_wrap%last_output_time = current_time

      write(*,'(A,A)') 'CATChem: Wrote diagnostic output to ', trim(filename)

   end subroutine catchem_diagnostics_write

   !> \brief Write diagnostics for a specific process
   !!
   !! \param process_name Name of the process
   !! \param filename Output filename
   !! \param rc Return code
   subroutine write_process_diagnostics(cc_wrap, process_name, filename, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      character(len=*), intent(in) :: process_name
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      type(DiagnosticManagerType), pointer :: diag_mgr => null()
      type(DiagnosticRegistryType), pointer :: registry => null()
      !type(cc_wrap_type), pointer :: cc_wrap
      character(len=64), allocatable :: field_names(:)
      integer :: num_fields, i, data_type
      real(fp) :: scalar_value
      real(fp), pointer :: array_1d_ptr(:) => null()
      real(fp), pointer :: array_2d_ptr(:,:) => null()
      real(fp), pointer :: array_3d_ptr(:,:,:) => null()
      character(len=128) :: description
      character(len=32) :: units
      character(len=64) :: field_name

      rc = CC_SUCCESS

      ! Get process-local state
      !cc_wrap => get_cc_wrap()

      ! Get diagnostic manager and process registry
      diag_mgr => cc_wrap%catchem_model%get_diagnostic_manager()
      call diag_mgr%get_process_registry(process_name, registry, rc)
      if (rc /= CC_SUCCESS .or. .not. associated(registry)) return

      ! Get field count
      num_fields = registry%get_field_count()
      if (num_fields == 0) return

      ! Get field names
      allocate(field_names(num_fields))
      call registry%list_fields(field_names, num_fields)

      ! Write each field
      do i = 1, num_fields
         field_name = trim(field_names(i))

         ! Get field value with metadata
         call diag_mgr%get_field_value(process_name, field_name, &
            scalar_value, array_1d_ptr, array_2d_ptr, array_3d_ptr, &
            data_type, description, units, rc)

         if (rc /= CC_SUCCESS) then
            write(*,'(A,A,A,A)') 'Warning: Failed to get field: ', trim(process_name), '.', trim(field_name)
            cycle
         end if

         ! Write field to NetCDF using AQMIO
         call write_diagnostic_field(cc_wrap, field_name, data_type, scalar_value, &
            array_1d_ptr, array_2d_ptr, array_3d_ptr, &
            description, units, filename, rc)

         if (rc /= CC_SUCCESS) then
            write(*,'(A,A,A,A)') 'Error: Failed to write field: ', trim(process_name), '.', trim(field_name)
            return
         end if

         ! Clean up pointers
         if (associated(array_1d_ptr)) nullify(array_1d_ptr)
         if (associated(array_2d_ptr)) nullify(array_2d_ptr)
         if (associated(array_3d_ptr)) nullify(array_3d_ptr)
      end do

      deallocate(field_names)

   end subroutine write_process_diagnostics

   !> \brief Write individual diagnostic field to NetCDF
   !!
   !! \param field_name Field name for NetCDF variable
   !! \param data_type Type of diagnostic data
   !! \param scalar_value Scalar value (if applicable)
   !! \param array_1d_ptr 1D array pointer (if applicable)
   !! \param array_2d_ptr 2D array pointer (if applicable)
   !! \param array_3d_ptr 3D array pointer (if applicable)
   !! \param description Field description for metadata
   !! \param units Field units for metadata
   !! \param filename Output filename
   !! \param rc Return code
   subroutine write_diagnostic_field(cc_wrap, field_name, data_type, scalar_value, &
      array_1d_ptr, array_2d_ptr, array_3d_ptr, &
      description, units, filename, rc)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: data_type
      real(fp), intent(in) :: scalar_value
      real(fp), pointer, intent(in) :: array_1d_ptr(:)
      real(fp), pointer, intent(in) :: array_2d_ptr(:,:)
      real(fp), pointer, intent(in) :: array_3d_ptr(:,:,:)
      character(len=*), intent(in) :: description
      character(len=*), intent(in) :: units
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      type(ESMF_Field) :: esmf_field
      type(ESMF_Info) :: info
      real(ESMF_KIND_R4), pointer :: field_data_2d(:,:) => null()
      real(ESMF_KIND_R4), pointer :: field_data_3d(:,:,:) => null()
      integer :: i, j, k, time_slice

      rc = CC_SUCCESS

      ! Get current time slice - this will be the same for all fields in this diagnostic write
      time_slice = cc_wrap%current_time_slice

      ! Create appropriate ESMF field based on data type
      select case (data_type)
       case (DIAG_REAL_2D)
         if (.not. associated(array_2d_ptr)) then
            rc = CC_FAILURE
            return
         end if
         esmf_field = ESMF_FieldCreate(cc_wrap%grid, &
            name=trim(field_name), &
            typekind=ESMF_TYPEKIND_R4, &
            rc=rc)
         if (rc /= ESMF_SUCCESS) return

         !set some info for the field
         call ESMF_InfoGetFromHost(esmf_field, info, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return  ! bail out

         call ESMF_InfoSet(info, "units", trim(units), rc=rc)
         call ESMF_InfoSet(info, "description", trim(description), rc=rc)

         !set values
         call ESMF_FieldGet(esmf_field, farrayPtr=field_data_2d, rc=rc)
         if (rc /= ESMF_SUCCESS) return
         do j = 1, size(array_2d_ptr, 2)
            do i = 1, size(array_2d_ptr, 1)
               field_data_2d(i, j) = real(array_2d_ptr(i, j), ESMF_KIND_R4)
            end do
         end do
         call AQMIO_Write(cc_wrap%iocomp, (/esmf_field/), timeSlice=time_slice, fileName=trim(filename), &
            iofmt=AQMIO_FMT_NETCDF, rc=rc)

       case (DIAG_REAL_3D)
         if (.not. associated(array_3d_ptr)) then
            rc = CC_FAILURE
            return
         end if
         esmf_field = ESMF_FieldCreate(cc_wrap%grid, &
            name=trim(field_name), &
            typekind=ESMF_TYPEKIND_R4, &
            ungriddedLBound=(/1/), &
            ungriddedUBound=(/size(array_3d_ptr, 3)/), &
            rc=rc)
         if (rc /= ESMF_SUCCESS) return

         !set some info for the field
         call ESMF_InfoGetFromHost(esmf_field, info, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return  ! bail out

         call ESMF_InfoSet(info, "units", trim(units), rc=rc)
         call ESMF_InfoSet(info, "description", trim(description), rc=rc)

         !set values
         call ESMF_FieldGet(esmf_field, farrayPtr=field_data_3d, rc=rc)
         if (rc /= ESMF_SUCCESS) return
         do k = 1, size(array_3d_ptr, 3)
            do j = 1, size(array_3d_ptr, 2)
               do i = 1, size(array_3d_ptr, 1)
                  field_data_3d(i, j, k) = real(array_3d_ptr(i, j, k), ESMF_KIND_R4)
               end do
            end do
         end do
         call AQMIO_Write(cc_wrap%iocomp, (/esmf_field/), timeSlice=time_slice, fileName=trim(filename), &
            iofmt=AQMIO_FMT_NETCDF, rc=rc)

       case default
         rc = CC_FAILURE
         return
      end select

      ! TODO: Add NetCDF attributes for description and units
      ! This would require extending AQMIO or using NetCDF directly
      ! For now, we rely on the working AQMIO functionality

      ! Clean up
      if (ESMF_FieldIsCreated(esmf_field)) then
         call ESMF_FieldDestroy(esmf_field, rc=rc)
      end if

   end subroutine write_diagnostic_field

   !> \brief Update time variable in NetCDF file
   !!
   !! This function handles proper time series management by either:
   !! - Reading existing time values and appending a new time step, or
   !! - Creating the time variable for the first time
   !! Uses 1D LocStream-based fields and new AQMIO_Write1D function for proper
   !! 1D field handling without MPI operations.
   !!
   !! \param cc_wrap CATChem wrapper containing ESMF components
   !! \param filename NetCDF filename
   !! \param current_time Current simulation time to append
   !! \param time_slice Returns the time slice number that was written
   !! \param rc Return code
   subroutine update_time_variable(cc_wrap, filename, current_time, time_slice, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      character(len=*), intent(in) :: filename
      type(ESMF_Time), intent(in) :: current_time
      integer, intent(out) :: time_slice
      integer, intent(out) :: rc

      ! Local variables for direct data I/O approach
      integer(ESMF_KIND_I4) :: new_time_data(1)
      type(ESMF_Time) :: reference_time
      type(ESMF_TimeInterval) :: time_diff
      integer(ESMF_KIND_I8) :: time_seconds
      type(ESMF_VM) :: vm
      integer :: ibuf(1)  ! Buffer for MPI broadcast

      rc = CC_SUCCESS

      ! Create reference time (Unix epoch: 1970-01-01 00:00:00)
      call ESMF_TimeSet(reference_time, yy=1970, mm=1, dd=1, h=0, m=0, s=0, rc=rc)
      if (rc /= ESMF_SUCCESS) return

      ! Convert current time to epoch seconds
      time_diff = current_time - reference_time
      call ESMF_TimeIntervalGet(time_diff, s_i8=time_seconds, rc=rc)
      if (rc /= ESMF_SUCCESS) return

      new_time_data(1) = int(time_seconds, ESMF_KIND_I4)

      ! Use the new direct write function with append=true
      ! This automatically handles reading existing data and appending the new time
      call AQMIO_Write1D(filename, "time", append=.true., del_old_file=.true., rc=rc, &
         data_i4=new_time_data, current_size=time_slice, &
         iocomp=cc_wrap%iocomp)
      if (rc /= ESMF_SUCCESS) return

      ! Broadcast time_slice from I/O PET to all other PETs so they have the correct value
      ! Get VM from the IOComp for broadcasting
      call ESMF_GridCompGet(cc_wrap%iocomp, vm=vm, rc=rc)
      if (rc == ESMF_SUCCESS) then
         ! Use buffer for broadcast (ESMF_VMBroadcast expects arrays)
         ibuf(1) = time_slice
         call ESMF_VMBroadcast(vm, ibuf, 1, 0, rc=rc)
         time_slice = ibuf(1)
      end if

   end subroutine update_time_variable

   !> \brief Initialize output timing
   !!
   !! \param start_time Simulation start time
   !! \param rc Return code
   subroutine initialize_output_timing(cc_wrap, start_time, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_Time), intent(in) :: start_time
      integer, intent(out) :: rc

      !type(cc_wrap_type), pointer :: cc_wrap

      rc = CC_SUCCESS

      ! Get process-local state
      !cc_wrap => get_cc_wrap()

      ! Create output interval from frequency (in seconds)
      call ESMF_TimeIntervalSet(cc_wrap%output_interval, s=cc_wrap%output_frequency, rc=rc)
      if (rc /= ESMF_SUCCESS) return

      ! Set initial output time (start time minus interval so first check will trigger output)
      cc_wrap%last_output_time = start_time - cc_wrap%output_interval

      cc_wrap%output_timing_initialized = .true.

   end subroutine initialize_output_timing

   !> \brief Check if it's time to write diagnostic output
   !!
   !! \param current_time Current simulation time
   !! \param time_to_write True if it's time to write
   !! \param rc Return code
   subroutine check_diagnostic_output_time(cc_wrap, current_time, time_to_write, time_on_file, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_Time), intent(in) :: current_time
      logical, intent(out) :: time_to_write
      type(ESMF_Time), intent(out) :: time_on_file
      integer, intent(out) :: rc

      !type(cc_wrap_type), pointer :: cc_wrap
      type(ESMF_Time) :: next_output_time
      type(ESMF_Time) :: next_time

      rc = CC_SUCCESS
      time_to_write = .false.

      ! Get process-local state
      !cc_wrap => get_cc_wrap()

      ! Calculate next output time
      next_output_time = cc_wrap%last_output_time + cc_wrap%output_interval

      ! Check if current time is at or past next output time
      if (current_time >= next_output_time) then
         time_to_write = .true.
         time_on_file = current_time
      end if

      !Note that the last run cycle is one time step before the end time, so we do not have the last hour saved out.
      !To ensure the last time step is written, we can force output if we are within one time step of the end time.
      ! next_time = current_time + cc_wrap%timeStep
      ! if ( (cc_wrap%endTime == next_time) ) then
      !   time_to_write = .true.
      !   time_on_file = cc_wrap%endTime
      ! end if


      ! !!!!!log to debug
      ! block
      !   integer :: year, month, day, hour, minute, second
      !   character(len=64) :: time_str_current, time_str_next, time_str_end

      !   call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, &
      !                    h=hour, m=minute, s=second, rc=rc)
      !   write(time_str_current, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
      !         year, month, day, hour, minute, second

      !   call ESMF_TimeGet(cc_wrap%endTime, yy=year, mm=month, dd=day, &
      !                    h=hour, m=minute, s=second, rc=rc)
      !   write(time_str_next, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
      !         year, month, day, hour, minute, second

      !   write(*,'(A,A,A)') 'Debug: Current time: ', trim(time_str_current), ', End time: ', trim(time_str_next)

      ! end block

   end subroutine check_diagnostic_output_time

   !> \brief Generate filename for diagnostic output
   !!
   !! \param current_time Current simulation time
   !! \param filename Generated filename
   !! \param rc Return code
   subroutine generate_diagnostic_filename(cc_wrap, time_on_file, filename, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_Time), intent(in) :: time_on_file
      character(len=*), intent(out) :: filename
      integer, intent(out) :: rc

      !type(cc_wrap_type), pointer :: cc_wrap
      type(ESMF_VM) :: vm
      integer :: year, month, day, hour, minute, second, localPet
      character(len=256) :: time_string
      logical :: dir_exists

      rc = CC_SUCCESS

      ! Only have PET 0 create directories to avoid race conditions
      call ESMF_GridCompGet(cc_wrap%iocomp, vm=vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      call ESMF_VMGet(vm, localPet=localPet, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      if (localPet == 0) then
         ! check if output directory exists, create if not
         inquire(file=trim(cc_wrap%output_directory), exist=dir_exists)
         if (.not. dir_exists) then
            ! Create directory
            call system('mkdir -p ' // trim(cc_wrap%output_directory))
         end if
      end if

      ! Get time components
      call ESMF_TimeGet(time_on_file, yy=year, mm=month, dd=day, &
         h=hour, m=minute, s=second, rc=rc)
      if (rc /= ESMF_SUCCESS) return

      ! Create filename: output_directory/prefix_YYYYMMDD_HHMMSS.nc
      write(time_string, '(I4.4,I2.2,I2.2,A,I2.2,I2.2,I2.2)') &
         year, month, day, '_', hour, minute, second

      filename = trim(cc_wrap%output_directory) // '/' // trim(cc_wrap%output_prefix) // '_' // trim(time_string) // '.nc'

      !This will save all hours to the same file
      !filename = trim(cc_wrap%output_directory) // '/' // trim(cc_wrap%output_prefix) // '.nc'

   end subroutine generate_diagnostic_filename


   ! Load field configuration from YAML file
   !!
   !! \param info  tracerinfo from NUOPC
   !! \param key   infor intended to get
   !! \param values  infor values returned
   !! \param rc  return status
   !!
   subroutine TracerInfoGet(info, key, values, rc)
      ! -- interface variables
      type(ESMF_Info),               intent(in)  :: info
      character(len=*),              intent(in)  :: key
      character(len=*), allocatable, intent(out) :: values(:)
      integer,          optional,    intent(out) :: rc

      ! -- local variables
      integer :: localrc
      logical :: isKeyFound

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

      ! -- check if metadata key is present in ESMF_Info object
      isKeyFound = ESMF_InfoIsPresent(info, trim(key), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return

      if (isKeyFound) then
         isKeyFound = ESMF_InfoIsSet(info, trim(key), rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return
      end if

      if (isKeyFound) then
         call ESMF_InfoGetAlloc(info, trim(key), values, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return
      end if

   end subroutine TracerInfoGet

   ! Load field configuration from YAML file
   !!
   !! \param  config_file Configuration file path
   !! \param errflg      Error flag
   !! \param errmsg      Error message
   !!
   subroutine load_field_config(config_file, errflg, errmsg)

      character(len=*), intent(in) :: config_file
      integer, intent(out) :: errflg
      character(len=*), intent(out) :: errmsg

      !type(cc_wrap_type), pointer :: cc_wrap

      errflg = CC_SUCCESS
      errmsg = ''

      ! Parse import fields
      call parse_field_section(config_file, 'import_fields', field_config%import_fields, &
         field_config%n_import_fields, errflg, errmsg)
      if (errflg /= CC_SUCCESS) return

      ! Parse export fields
      call parse_field_section(config_file, 'export_fields', field_config%export_fields, &
         field_config%n_export_fields, errflg, errmsg)
      if (errflg /= CC_SUCCESS) return

   end subroutine load_field_config

   !> Parse a field section (import_fields or export_fields) from YAML file
   !!
   !! \param filename YAML configuration file
   !! \param section_name Section name ('import_fields' or 'export_fields')
   !! \param fields Array to store parsed fields
   !! \param n_fields Number of fields found
   !! \param errflg Error flag
   !! \param errmsg Error message
   !!
   subroutine parse_field_section(filename, section_name, fields, n_fields, errflg, errmsg)
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: section_name
      type(field_mapping_type), allocatable, intent(out) :: fields(:)
      integer, intent(out) :: n_fields
      integer, intent(out) :: errflg
      character(len=*), intent(out) :: errmsg

      integer :: unit_num, io_stat, indent_level, section_indent
      character(len=256) :: line, trimmed_line, field_name, field_value
      logical :: in_section, found_section
      integer :: line_number, field_idx, colon_pos
      type(field_mapping_type), allocatable :: temp_fields(:)
      type(field_mapping_type) :: current_field
      logical :: in_field_item, field_already_saved

      errflg = CC_SUCCESS
      errmsg = ''
      n_fields = 0
      in_section = .false.
      found_section = .false.
      section_indent = -1
      line_number = 0
      field_idx = 0
      in_field_item = .false.
      field_already_saved = .false.

      ! Initialize current field
      current_field%standard_name = ''
      current_field%catchem_var = ''
      current_field%dimensions = 0
      current_field%units = ''
      current_field%optional = .false.

      ! Open file for reading
      open(newunit=unit_num, file=trim(filename), status='old', action='read', iostat=io_stat)
      if (io_stat /= 0) then
         write(errmsg, '(A,A)') 'Cannot open configuration file: ', trim(filename)
         errflg = CC_FAILURE
         return
      endif

      ! Allocate temporary storage for up to 50 fields
      allocate(temp_fields(50))

      ! Read file line by line
      do
         read(unit_num, '(A)', iostat=io_stat) line
         if (io_stat /= 0) exit  ! End of file or error

         line_number = line_number + 1

         ! Remove inline comments - everything after '#' character
         if (index(line, '#') > 0) then
            line = line(1:index(line, '#')-1)
         endif

         trimmed_line = trim(adjustl(line))

         ! Skip empty lines (comments have already been stripped)
         if (len_trim(trimmed_line) == 0) cycle

         ! Calculate indentation level
         do indent_level = 1, len_trim(line)
            if (line(indent_level:indent_level) /= ' ') exit
         end do
         indent_level = indent_level - 1

         ! Look for section header
         if (index(trimmed_line, ':') > 0) then
            colon_pos = index(trimmed_line, ':')
            field_name = trimmed_line(1:colon_pos-1)
            field_name = trim(adjustl(field_name))

            ! Check if we found our target section
            if (trim(field_name) == trim(section_name) .and. indent_level == 0) then
               in_section = .true.
               found_section = .true.
               section_indent = indent_level
               cycle
            endif

            ! If we're already in a section and encounter another top-level section, exit
            if (in_section .and. indent_level == 0 .and. trim(field_name) /= trim(section_name)) then
               ! Save the last field if we're still processing one
               if (in_field_item .and. current_field%standard_name /= '') then
                  n_fields = n_fields + 1
                  if (n_fields <= size(temp_fields)) then
                     temp_fields(n_fields) = current_field
                  endif
                  field_already_saved = .true.  ! Mark that we've saved the field
               endif
               exit
            endif

            ! Process items within the section
            if (in_section .and. indent_level > section_indent) then

               ! Look for array items (lines starting with "- ")
               if (index(trimmed_line, '- ') == 1) then
                  ! Save previous field if we have one
                  if (in_field_item .and. current_field%standard_name /= '') then
                     n_fields = n_fields + 1
                     if (n_fields <= size(temp_fields)) then
                        temp_fields(n_fields) = current_field
                     endif
                  endif

                  ! Start new field item
                  in_field_item = .true.
                  current_field%standard_name = ''
                  current_field%catchem_var = ''
                  current_field%dimensions = 0
                  current_field%units = ''
                  current_field%optional = .false.

                  ! Parse the first property if it's on the same line as the dash
                  if (len_trim(trimmed_line) > 2) then
                     trimmed_line = trim(adjustl(trimmed_line(3:)))  ! Remove "- "
                     if (index(trimmed_line, ':') > 0) then
                        colon_pos = index(trimmed_line, ':')
                        field_name = trim(adjustl(trimmed_line(1:colon_pos-1)))
                        field_value = trim(adjustl(trimmed_line(colon_pos+1:)))
                        call parse_field_property(field_name, field_value, current_field)
                     endif
                  endif

               elseif (in_field_item .and. indent_level > section_indent + 2) then
                  ! Parse field properties
                  if (index(trimmed_line, ':') > 0) then
                     colon_pos = index(trimmed_line, ':')
                     field_name = trim(adjustl(trimmed_line(1:colon_pos-1)))
                     field_value = trim(adjustl(trimmed_line(colon_pos+1:)))
                     call parse_field_property(field_name, field_value, current_field)
                  endif
               endif

            elseif (in_section .and. indent_level <= section_indent) then
               ! We've left our section
               if (in_field_item .and. current_field%standard_name /= '') then
                  n_fields = n_fields + 1
                  if (n_fields <= size(temp_fields)) then
                     temp_fields(n_fields) = current_field
                  endif
                  field_already_saved = .true.  ! Mark that we've saved the field
               endif
               exit
            endif
         endif
      end do

      ! Save the last field if we're still processing one AND it hasn't been saved yet
      if (in_field_item .and. current_field%standard_name /= '' .and. .not. field_already_saved) then
         n_fields = n_fields + 1
         if (n_fields <= size(temp_fields)) then
            temp_fields(n_fields) = current_field
         endif
      endif

      close(unit_num)

      ! Check if we found the section
      if (.not. found_section) then
         write(errmsg, '(A,A,A)') 'Section "', trim(section_name), '" not found in configuration file'
         errflg = CC_FAILURE
         deallocate(temp_fields)
         return
      endif

      ! Allocate final array and copy data
      if (n_fields > 0) then
         allocate(fields(n_fields))
         fields(1:n_fields) = temp_fields(1:n_fields)
      endif

      deallocate(temp_fields)

   end subroutine parse_field_section

   !> Parse a field property and set it in the field structure
   !!
   !! \param property_name Name of the property
   !! \param property_value Value of the property
   !! \param field Field structure to update
   !!
   subroutine parse_field_property(property_name, property_value, field)
      character(len=*), intent(in) :: property_name
      character(len=*), intent(in) :: property_value
      type(field_mapping_type), intent(inout) :: field

      character(len=256) :: clean_value
      integer :: read_stat

      ! Remove quotes from string values
      clean_value = property_value
      if (len_trim(clean_value) >= 2) then
         if ((clean_value(1:1) == '"' .and. clean_value(len_trim(clean_value):len_trim(clean_value)) == '"') .or. &
            (clean_value(1:1) == "'" .and. clean_value(len_trim(clean_value):len_trim(clean_value)) == "'")) then
            clean_value = clean_value(2:len_trim(clean_value)-1)
         endif
      endif

      select case (trim(property_name))
       case ('standard_name')
         field%standard_name = trim(clean_value)
       case ('catchem_var')
         field%catchem_var = trim(clean_value)
       case ('dimensions')
         read(clean_value, *, iostat=read_stat) field%dimensions
         if (read_stat /= 0) field%dimensions = 0
       case ('units')
         field%units = trim(clean_value)
       case ('optional')
         select case (trim(clean_value))
          case ('true', 'True', 'TRUE', '.true.')
            field%optional = .true.
          case ('false', 'False', 'FALSE', '.false.')
            field%optional = .false.
          case default
            field%optional = .false.
         end select
      end select

   end subroutine parse_field_property

   !> \brief Get number of import fields (MPI-safe accessor)
   !!
   !! \return Number of import fields configured
   function get_n_import_fields(cc_wrap) result(n_fields)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      integer :: n_fields
      !type(cc_wrap_type), pointer :: cc_wrap

      !cc_wrap => get_cc_wrap()
      n_fields = cc_wrap%field_config%n_import_fields
   end function get_n_import_fields

   !> \brief Get import field information (MPI-safe accessor)
   !!
   !! \param field_index Index of the field (1-based)
   !! \param standard_name NUOPC standard name
   !! \param optional Whether field is optional
   function get_import_field_info(cc_wrap, field_index, standard_name, optional) result(success)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      integer, intent(in) :: field_index
      character(len=*), intent(out) :: standard_name
      logical, intent(out) :: optional
      logical :: success
      !type(cc_wrap_type), pointer :: cc_wrap

      success = .false.
      !cc_wrap => get_cc_wrap()

      if (field_index > 0 .and. field_index <= cc_wrap%field_config%n_import_fields) then
         standard_name = cc_wrap%field_config%import_fields(field_index)%standard_name
         optional = cc_wrap%field_config%import_fields(field_index)%optional
         success = .true.
      end if
   end function get_import_field_info

   !> \brief Get number of export fields (MPI-safe accessor)
   !!
   !! \return Number of export fields configured
   function get_n_export_fields(cc_wrap) result(n_fields)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      integer :: n_fields
      !type(cc_wrap_type), pointer :: cc_wrap

      !cc_wrap => get_cc_wrap()
      n_fields = cc_wrap%field_config%n_export_fields
   end function get_n_export_fields

   !> \brief Get export field information (MPI-safe accessor)
   !!
   !! \param field_index Index of the field (1-based)
   !! \param standard_name NUOPC standard name
   !! \param optional Whether field is optional
   function get_export_field_info(cc_wrap, field_index, standard_name, optional) result(success)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      integer, intent(in) :: field_index
      character(len=*), intent(out) :: standard_name
      logical, intent(out) :: optional
      logical :: success
      !type(cc_wrap_type), pointer :: cc_wrap

      success = .false.
      !cc_wrap => get_cc_wrap()

      if (field_index > 0 .and. field_index <= cc_wrap%field_config%n_export_fields) then
         standard_name = cc_wrap%field_config%export_fields(field_index)%standard_name
         optional = cc_wrap%field_config%export_fields(field_index)%optional
         success = .true.
      end if
   end function get_export_field_info

end module catchem_nuopc_interface
