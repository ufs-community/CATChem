!> \file catchem_driver.F90
!! \brief NUOPC_Driver component for standalone CATChem execution
!!
!! \details
!! This module specializes \c NUOPC_Driver to run CATChem as a standalone
!! application. It follows the same pattern as the UFS Driver
!! (UFSDriver.F90) and the CECE standalone driver: the driver owns the
!! top-level clock and adds one or more child model components, wiring them
!! together through a run sequence.
!!
!! For the first version a single child component is added -- the CATChem
!! NUOPC cap (\c cc_nuopc). The structure is deliberately kept extensible so
!! additional components (e.g. a data/forcing component, a second chemistry
!! component, or a mediator) can be added later by:
!!   1. \c use-ing the new component's SetServices, and
!!   2. adding a \c NUOPC_DriverAddComp call in \c SetModelServices, and
!!   3. extending the \c runSeq:: block in the driver configure file.
!!
!! Configuration is read from an ESMF configure file (default
!! \c catchem_standalone.configure) that controls the clock (start/stop
!! time, time step) and the grid (column vs gridded, dimensions, extents).
!! The CATChem science configuration remains a separate YAML file consumed
!! by the cap.
!!
!! \par Component hierarchy
!! \code
!!   CATChemApp (main program)
!!        |
!!   catchem_driver (NUOPC_Driver)   <-- this module
!!        |
!!   cc_nuopc (CATChem NUOPC cap / NUOPC_Model)
!! \endcode
!!
!! \author CATChem standalone driver
!! \ingroup catchem_nuopc_group
module catchem_driver

   use ESMF
   use NUOPC
   use NUOPC_Driver, &
      driver_routine_SS             => SetServices, &
      driver_label_SetModelServices => label_SetModelServices, &
      driver_label_SetRunSequence   => label_SetRunSequence

   use cc_nuopc, only: &
      catchem_cap_SS => SetServices, &
      catchem_cap_SV => SetVM, &
      cc_nuopc_set_config_file, &
      cc_nuopc_set_field_mapping_file

   use catchem_standalone_grid_mod, only: CATChemGridConfig, create_standalone_grid

   implicit none
   private

   public :: SetServices
   public :: set_driver_config_file

   character(len=*), parameter :: u_FILE_u = __FILE__

   !> Path to the driver ESMF configure file (set by the main application)
   character(len=512), save :: g_driver_cfg_file = "catchem_standalone.configure"

contains

   !> \brief Set the driver configure file path (called by the application)
   subroutine set_driver_config_file(config_file)
      character(len=*), intent(in) :: config_file
      g_driver_cfg_file = config_file
   end subroutine set_driver_config_file

   !> \brief Register the driver's specialized methods
   subroutine SetServices(driver, rc)
      type(ESMF_GridComp)  :: driver
      integer, intent(out) :: rc

      rc = ESMF_SUCCESS

      ! Derive from NUOPC_Driver
      call NUOPC_CompDerive(driver, driver_routine_SS, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      ! Specialize: define child components and the clock
      call NUOPC_CompSpecialize(driver, specLabel=driver_label_SetModelServices, &
         specRoutine=SetModelServices, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      ! Specialize: define the run sequence
      call NUOPC_CompSpecialize(driver, specLabel=driver_label_SetRunSequence, &
         specRoutine=SetRunSequence, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

   end subroutine SetServices

   !> \brief Add child components and create the top-level clock
   subroutine SetModelServices(driver, rc)
      type(ESMF_GridComp)  :: driver
      integer, intent(out) :: rc

      type(ESMF_Config)       :: config
      type(ESMF_GridComp)     :: child
      type(ESMF_Grid)         :: grid
      type(ESMF_Time)         :: startTime, stopTime
      type(ESMF_TimeInterval) :: timeStep
      type(ESMF_Clock)        :: internalClock
      type(CATChemGridConfig) :: gridCfg

      character(len=64)  :: start_str, stop_str
      character(len=512) :: catchem_config_file, field_mapping_file
      integer :: timestep_sec

      rc = ESMF_SUCCESS

      ! ----------------------------------------------------------------
      ! Load the driver configuration
      ! ----------------------------------------------------------------
      config = ESMF_ConfigCreate(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_ConfigLoadFile(config, trim(g_driver_cfg_file), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call read_driver_config(config, start_str, stop_str, timestep_sec, &
         catchem_config_file, field_mapping_file, gridCfg, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      ! Make the driver configure file available to the standard NUOPC
      ! attribute machinery as well.
      call ESMF_GridCompSet(driver, config=config, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      ! ----------------------------------------------------------------
      ! Hand the CATChem cap its configuration file paths
      ! ----------------------------------------------------------------
      call cc_nuopc_set_config_file(trim(catchem_config_file))
      if (len_trim(field_mapping_file) > 0) then
         call cc_nuopc_set_field_mapping_file(trim(field_mapping_file))
      end if

      ! ----------------------------------------------------------------
      ! Build the standalone grid (column or gridded) and attach it to the
      ! child component so the cap's standalone initialization path can
      ! retrieve it via ESMF_GridCompGet(model, grid=grid, ...).
      ! ----------------------------------------------------------------
      call create_standalone_grid(gridCfg, grid, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      ! ----------------------------------------------------------------
      ! Add the CATChem cap as a child component
      ! ----------------------------------------------------------------
      call NUOPC_DriverAddComp(driver, "CATCHEM", catchem_cap_SS, &
         catchem_cap_SV, comp=child, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_GridCompSet(child, grid=grid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call NUOPC_CompAttributeSet(child, name="Verbosity", value="high", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      ! ----------------------------------------------------------------
      ! Create and set the top-level clock
      ! ----------------------------------------------------------------
      call ESMF_TimeSet(startTime, timeString=trim(start_str), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_TimeSet(stopTime, timeString=trim(stop_str), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_TimeIntervalSet(timeStep, s=timestep_sec, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      internalClock = ESMF_ClockCreate(name="catchem_driver_clock", &
         timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_GridCompSet(driver, clock=internalClock, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_LogWrite("catchem_driver: SetModelServices complete", &
         ESMF_LOGMSG_INFO, rc=rc)

   end subroutine SetModelServices

   !> \brief Define the run sequence for the driver
   !!
   !! Ingests a free-format \c runSeq:: block from the driver configure file
   !! if present. With a single component this is optional -- NUOPC will build
   !! a trivial run sequence automatically -- but ingesting it here keeps the
   !! driver ready for multi-component coupling.
   subroutine SetRunSequence(driver, rc)
      type(ESMF_GridComp)  :: driver
      integer, intent(out) :: rc

      type(ESMF_Config)      :: config
      type(NUOPC_FreeFormat) :: runSeqFF
      logical :: isPresent

      rc = ESMF_SUCCESS

      call ESMF_GridCompGet(driver, config=config, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      ! Only ingest a run sequence if the configure file actually defines one.
      isPresent = config_has_label(config, "runSeq::")
      if (.not. isPresent) then
         call ESMF_LogWrite("catchem_driver: no runSeq:: in config, " // &
            "using default single-component run sequence", ESMF_LOGMSG_INFO, rc=rc)
         return
      end if

      runSeqFF = NUOPC_FreeFormatCreate(config, label="runSeq::", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call NUOPC_DriverIngestRunSequence(driver, runSeqFF, &
         autoAddConnectors=.true., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call NUOPC_FreeFormatDestroy(runSeqFF, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

   end subroutine SetRunSequence

   !> \brief Read driver settings from an ESMF configure object
   subroutine read_driver_config(config, start_str, stop_str, timestep_sec, &
      catchem_config_file, field_mapping_file, gridCfg, rc)
      type(ESMF_Config),       intent(inout) :: config
      character(len=*),        intent(out)   :: start_str, stop_str
      integer,                 intent(out)   :: timestep_sec
      character(len=*),        intent(out)   :: catchem_config_file
      character(len=*),        intent(out)   :: field_mapping_file
      type(CATChemGridConfig), intent(out)   :: gridCfg
      integer,                 intent(out)   :: rc

      rc = ESMF_SUCCESS

      ! -- Clock control
      call ESMF_ConfigGetAttribute(config, start_str, &
         label="start_time:", default="2020-01-01T00:00:00", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_ConfigGetAttribute(config, stop_str, &
         label="stop_time:", default="2020-01-01T06:00:00", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_ConfigGetAttribute(config, timestep_sec, &
         label="timestep_seconds:", default=3600, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      ! -- CATChem science configuration files
      call ESMF_ConfigGetAttribute(config, catchem_config_file, &
         label="catchem_config_file:", default="CATChem_new_config.yml", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_ConfigGetAttribute(config, field_mapping_file, &
         label="field_mapping_file:", default="", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      ! -- Grid configuration
      call ESMF_ConfigGetAttribute(config, gridCfg%mode, &
         label="grid_mode:", default="column", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_ConfigGetAttribute(config, gridCfg%nx, &
         label="grid_nx:", default=1, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_ConfigGetAttribute(config, gridCfg%ny, &
         label="grid_ny:", default=1, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_ConfigGetAttribute(config, gridCfg%nz, &
         label="grid_nz:", default=72, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_ConfigGetAttribute(config, gridCfg%lon_start, &
         label="grid_lon_start:", default=0.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_ConfigGetAttribute(config, gridCfg%lon_end, &
         label="grid_lon_end:", default=360.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_ConfigGetAttribute(config, gridCfg%lat_start, &
         label="grid_lat_start:", default=-90.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_ConfigGetAttribute(config, gridCfg%lat_end, &
         label="grid_lat_end:", default=90.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_ConfigGetAttribute(config, gridCfg%column_lon, &
         label="column_lon:", default=0.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

      call ESMF_ConfigGetAttribute(config, gridCfg%column_lat, &
         label="column_lat:", default=0.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=u_FILE_u)) return

   end subroutine read_driver_config

   !> \brief Return .true. if the configure object contains the given label
   logical function config_has_label(config, label) result(found)
      type(ESMF_Config), intent(inout) :: config
      character(len=*),  intent(in)    :: label
      integer :: rc

      found = .false.
      call ESMF_ConfigFindLabel(config, label=label, isPresent=found, rc=rc)
      if (rc /= ESMF_SUCCESS) found = .false.
   end function config_has_label

end module catchem_driver
