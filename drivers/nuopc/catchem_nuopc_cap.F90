!> \file catchem_nuopc_cap.F90
!! \brief NUOPC cap for CATChem atmospheric chemistry model
!!
!! \defgroup catchem_nuopc_group CATChem NUOPC Interface
!! \brief NUOPC interface drivers and utilities for CATChem
!! \ingroup catchem
!!
!! This group contains all NUOPC-compliant interface modules and utilities
!! for integrating CATChem with the NUOPC framework, including the main
!! cap module, data transformation utilities, and I/O capabilities.
!!
!! \details
!! This module provides the NUOPC (National Unified Operational Prediction
!! Capability) cap for the CATChem (Configurable ATmospheric Chemistry) model.
!! The cap enables CATChem to run within the NUOPC/ESMF framework as a
!! component in coupled Earth system models.
!!
!! The cap implements the standard NUOPC phases:
!! - Initialize Phase 1: Advertise import and export fields
!! - Initialize Phase 2: Realize fields and initialize the CATChem model
!! - Run: Execute chemistry calculations and data exchange
!! - Finalize: Clean up resources and finalize CATChem processes
!!
!! Key features:
!! - Standard NUOPC/ESMF compliance for easy integration
!! - Flexible field mapping and data exchange
!! - Support for various grid configurations
!! - Configurable chemistry processes and diagnostics
!! - Parallel execution capabilities
!! - Error handling and logging
!!
!! \note This cap follows NUOPC conventions and requires ESMF/NUOPC libraries
!!
!! \author Barry Baker & Wei Li, NOAA/OAR/ARL
!! \date November 2024
!! \ingroup catchem_nuopc_group

module cc_nuopc
! Renamed from catchem_nuopc_cap to aqm for UFS Driver compatibility
! UFS expects: use aqm, only: AQM_SS => SetServices (after FRONT_AQM=aqm substitution)

   use ESMF
   use NUOPC
   use NUOPC_Model, &
      modelSS        => SetServices, &
      model_label_Advertise       => label_Advertise,      &
      model_label_DataInitialize  => label_DataInitialize, &
      model_label_Advance => label_Advance, &
      model_label_CheckImport => label_CheckImport, &
      model_label_SetRunClock => label_SetRunClock, &
      model_label_Finalize        => label_Finalize

   use catchem_nuopc_interface

   implicit none

   private

   public :: SetServices

   !> \brief Component configuration parameters
   !! \{
   character(len=256), save :: config_file = 'CATChem_new_config.yml' !< Configuration file path
   character(len=256), save :: field_mapping_file = 'CATChem_field_mapping.yml' !< Field mapping file path
   !logical, save :: do_chemistry = .true.                         !< Enable chemistry calculations
   !! \}

contains

   !> Set services for the CATChem NUOPC cap
   !!
   !! This is the main entry point called by NUOPC to set up the CATChem component.
   !! It registers the initialize, run, and finalize phase entry points and sets
   !! component metadata attributes.
   !!
   !! @param model NUOPC model component to configure
   !! @param rc ESMF return code (ESMF_SUCCESS on success)
   !!
   !! This routine performs the following setup operations:
   !! - Derives the component from the NUOPC model template
   !! - Sets component metadata (name, version)
   !! - Registers entry points for all required NUOPC phases:
   !!   - IPDv00p1: Initialize Phase 1 (advertise fields)
   !!   - IPDv00p2: Initialize Phase 2 (realize fields and initialize model)
   !!   - RunPhase1: Model advance (execute chemistry)
   !!   - FinalizePhase1: Model cleanup
   !!
   !! @note This routine must be called before any other component operations
   subroutine SetServices(model, rc)
      type(ESMF_GridComp)  :: model
      integer, intent(out) :: rc

      rc = ESMF_SUCCESS

      ! Set the model services
      call NUOPC_CompDerive(model, modelSS, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! ! Set component metadata
      ! call ESMF_AttributeSet(model, name="model_name", value="CATChem", rc=rc)
      ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      !   line=__LINE__, file=__FILE__)) return

      ! call ESMF_AttributeSet(model, name="model_version", value="1.0", rc=rc)
      ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      !   line=__LINE__, file=__FILE__)) return

      !Note NUOPC_CompSetEntryPoint is deprecated for newer version of ESMF
      call NUOPC_CompSpecialize(model, specLabel=model_label_Advertise, &
         specRoutine=InitializeP1, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      call NUOPC_CompSpecialize(model, specLabel=model_label_DataInitialize, &
         specRoutine=InitializeP2, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      call NUOPC_CompSpecialize(model, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

   end subroutine SetServices

   !> \brief Initialize Phase 1 - Advertise import and export fields
   !!
   !! In this phase, the component advertises what fields it can import
   !! (receive from other components) and export (provide to other components).
   !! This establishes the interface contract without creating actual field objects.
   !!
   !! \param[inout] model NUOPC model component
   !! \param[out] rc ESMF return code (ESMF_SUCCESS on success)
   !!
   !! \details
   !! This routine performs the following operations:
   !! - Retrieves the component's import and export states
   !! - Calls the field advertisement routine to declare all required fields
   !! - Sets up the field interface for meteorological inputs and chemistry outputs
   !! - Enables field matching and connection with other components
   !!
   !! Fields advertised typically include:
   !! - Import: Temperature, humidity, pressure, winds, surface properties
   !! - Export: Chemical species concentrations, deposition fluxes, emissions
   !!
   !! \note This is a standard NUOPC initialization phase that must complete
   !!       successfully before Phase 2 can proceed
   !!
   !! \ingroup catchem_nuopc_group
   subroutine InitializeP1(model, rc)
      type(ESMF_GridComp)  :: model
      integer, intent(out) :: rc

      type(ESMF_State) :: importState, exportState
      !type(ESMF_Field), pointer :: fieldList(:)
      character(len=*), parameter :: routine = 'InitializeP1'
      integer :: i
      character(len=218) :: errmsg

      rc = ESMF_SUCCESS

      ! Get import and export states
      call NUOPC_ModelGet(model, importState=importState, exportState=exportState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Load field configuration
      call load_field_config(field_mapping_file, rc, errmsg)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      !retrieve member list from import state, if any
      !nullify(fieldList)
      !call NUOPC_GetStateMemberLists(importState, fieldList=fieldList, nestedFlag=.true., rc=rc)
      !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      !  line=__LINE__,  file=__FILE__)) return

      !call ESMF_LogWrite("Import fields number: "//real_to_string(real(size(fieldList),ESMF_KIND_R8)), ESMF_LOGMSG_INFO, rc=rc)


      ! Advertise import fields only when it has nothing
      !if (size(fieldList) == 0) then
      ! Advertise import fields using MPI-safe accessor functions
      do i = 1, size(field_config%import_fields)
         !   block
         !     character(len=128) :: standard_name
         !     logical :: optional
         !     if (get_import_field_info(i, standard_name, optional)) then
         call NUOPC_Advertise(importState, &
            StandardName=trim(field_config%import_fields(i)%standard_name), &
            TransferOfferGeomObject="cannot provide", &
            SharePolicyField="share", rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         !     end if
         !   end block
      end do
      !end if

      ! retrieve member list from export state, if any
      !nullify(fieldList)
      !call NUOPC_GetStateMemberLists(exportState, fieldList=fieldList, nestedFlag=.true., rc=rc)
      !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      !  line=__LINE__,  file=__FILE__)) return

      !call ESMF_LogWrite("Export fields number: "//real_to_string(real(size(fieldList),ESMF_KIND_R8)), ESMF_LOGMSG_INFO, rc=rc)

      ! Advertise export fields only when it has nothing
      !if (size(fieldList) == 0) then
      ! Advertise export fields using MPI-safe accessor functions
      do i = 1, size(field_config%export_fields)
         !   block
         !     character(len=128) :: standard_name
         !     logical :: optional
         !     if (get_export_field_info(i, standard_name, optional)) then
         call NUOPC_Advertise(exportState, &
            StandardName=trim(field_config%export_fields(i)%standard_name), &
            TransferOfferGeomObject="cannot provide", &
            SharePolicyField="share", rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         !     end if
         !   end block
      end do
      !end if

      ! Log successful completion
      call ESMF_LogWrite("CATChem: Completed "//routine, ESMF_LOGMSG_INFO, rc=rc)

   end subroutine InitializeP1

   !> \brief Initialize Phase 2 - Realize fields and initialize CATChem model
   !!
   !! In this phase, the component creates the actual ESMF fields for the
   !! advertised imports and exports, and performs complete initialization
   !! of the CATChem model including memory allocation and configuration.
   !!
   !! \param[inout] model NUOPC model component
   !! \param[out] rc ESMF return code (ESMF_SUCCESS on success)
   !!
   !! \details
   !! This routine performs comprehensive model initialization:
   !! - Retrieves component information (local PET, PET count)
   !! - Gets the computational grid from the driver
   !! - Determines grid dimensions for memory allocation
   !! - Creates actual ESMF field objects on the grid
   !! - Reads CATChem configuration from file
   !! - Initializes all CATChem state containers and processes
   !! - Sets up chemistry, meteorology, emissions, and diagnostic systems
   !! - Prepares the model for time stepping
   !!
   !! The grid is typically provided by the parent driver or mediator
   !! component and defines the spatial discretization for all field
   !! operations and data exchange.
   !!
   !! \note This phase must complete successfully before any model
   !!       advance operations can be performed
   !!
   !! \ingroup catchem_nuopc_group
   subroutine InitializeP2(model, rc)
      type(ESMF_GridComp)  :: model
      integer, intent(out) :: rc

      type(ESMF_State) :: importState, exportState
      type(ESMF_Grid) :: grid
      type(ESMF_Array) :: array
      type(ESMF_Info) :: tracerInfo
      type(ESMF_Field), pointer :: fieldList(:)
      type(ESMF_Clock)          :: clock
      type(ESMF_Time)           :: currTime, startTime, stopTime
      type(ESMF_TimeInterval)   :: timeStep
      real(ESMF_KIND_R8), dimension(:,:), pointer :: coord
      real(ESMF_KIND_R8), dimension(:,:), allocatable :: lon
      real(ESMF_KIND_R8), dimension(:,:), allocatable :: lat
      type(ESMF_CoordSys_Flag)   :: coordSys
      character(len=*), parameter :: routine = 'InitializeP2'
      real(ESMF_KIND_R8), parameter :: rad_to_deg = 180._ESMF_KIND_R8 / 3.14159265358979323846_ESMF_KIND_R8
      real(ESMF_KIND_R8) :: convet_unit
      character(len=512) :: errmsg
      integer :: localPet, petCount
      integer :: im, jm  ! Grid dimensions
      integer :: item, coord_item, rank, localDeCount, numLevels, localDe, localrc, stat
      integer, dimension(2) :: lb, ub
      logical :: has_tracer_array
      type(CATChem_InternalState) :: is

      rc = ESMF_SUCCESS
      has_tracer_array = .false.

      call ESMF_LogWrite("CATChem: Enter InitializeP2", ESMF_LOGMSG_INFO, rc=rc)

      ! Get component information
      call ESMF_GridCompGet(model, localPet=localPet, petCount=petCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Get import and export states
      call NUOPC_ModelGet(model, importState=importState, exportState=exportState, modelClock=clock, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! -- get clock information
      call ESMF_ClockGet(clock, startTime=startTime, stopTime=stopTime, timeStep=timeStep, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! retrieve member list from import state, if any
      nullify(fieldList)
      call NUOPC_GetStateMemberLists(importState, fieldList=fieldList, nestedFlag=.true., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__,  file=__FILE__)) return

      ! retrieve number of vertical levels from imported fields
      if (associated(fieldList)) then
         do item = 1, size(fieldList)

            call ESMF_FieldGet(fieldList(item), rank=rank, localDeCount=localDeCount, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
            ! -- validate field data decomposition
            if (localDeCount /= 1) then
               call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="localDeCount must be 1", &
                  line=__LINE__, file=__FILE__, rcToReturn=rc)
            end if

            if (rank == 4) then !use tracer array to get domain
               has_tracer_array = .true.
               call ESMF_FieldGet(fieldList(item), array=array, grid=grid, &
                  ungriddedLBound=lb, ungriddedUBound=ub, rc=rc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=__FILE__)) return  ! bail out
               ! -- populate remaining output arguments
               numLevels = ub(1) - lb(1) + 1
               call ESMF_InfoGetFromHost(array, tracerInfo, rc=localrc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=__FILE__)) return  ! bail out

               do localDe = 0, localDeCount-1
                  ! -- get local coordinate arrays
                  call ESMF_GridGet(grid, coordSys=coordSys, rc=rc)
                  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                     line=__LINE__, file=__FILE__)) return  ! bail out

                  do coord_item = 1, 2
                     call ESMF_GridGetCoord(grid, coordDim=coord_item, staggerloc=ESMF_STAGGERLOC_CENTER, &
                        localDE=localDe, farrayPtr=coord, rc=rc)
                     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                        line=__LINE__, file=__FILE__)) return  ! bail out

                     if (coordSys == ESMF_COORDSYS_SPH_DEG) then
                        !coordinates are in degrees already
                        convet_unit = 1._ESMF_KIND_R8
                     else if (coordSys == ESMF_COORDSYS_SPH_RAD) then
                        !convert radians to degrees
                        convet_unit = rad_to_deg
                     else
                        call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
                           msg="Unsupported coordinate system - Failed to set coordinates for air quality model", &
                           line=__LINE__, file=__FILE__, rcToReturn=rc)
                        return  ! bail out
                     end if

                     select case (coord_item)
                      case(1)
                        lon = coord * convet_unit
                      case(2)
                        lat = coord * convet_unit
                      case default
                        !do nothing
                     end select
                  end do ! loop over coordinate dimensions
               end do ! loop over local DEs

            end if !rank = 4
         end do
         if (.not. has_tracer_array) then
            call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="tracer array is needed!", &
               line=__LINE__, file=__FILE__, rcToReturn=rc)
         end if

         deallocate(fieldList, stat=stat)
         if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg="Unable to deallocate internal memory", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out
         nullify(fieldList)

      end if

      ! Initialize CATChem using the interface (TODO: not provide nsoil, nsoiltype and nsurftype)
      call catchem_nuopc_init(model, config_file, lat, lon, numLevels, tracerInfo, grid, &
         startTime=startTime, stopTime=stopTime, timeStep=timeStep, clock=clock, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return  ! bail out

      ! -- indicate that data initialization is complete (breaking out of init-loop)
      call NUOPC_CompAttributeSet(model, &
         name="InitializeDataComplete", value="true", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return  ! bail out

   end subroutine InitializeP2

   !> \brief Model advance routine - Execute one time step of chemistry calculations
   !!
   !! This routine is called at each time step to advance the chemistry model.
   !! It handles the complete workflow of importing meteorological data,
   !! running chemistry calculations, and exporting the computed results.
   !!
   !! \param[inout] model NUOPC model component
   !! \param[out] rc ESMF return code (ESMF_SUCCESS on success)
   !!
   !! \details
   !! This routine performs the following operations for each time step:
   !! - Retrieves the component clock and current simulation time
   !! - Calculates the time step duration for chemistry integration
   !! - Imports meteorological fields from other model components
   !! - Transforms NUOPC field data into CATChem state format
   !! - Executes all enabled chemistry processes (if do_chemistry=true):
   !!   - Dust emission and transport
   !!   - Sea salt emission and transport
   !!   - Gas-phase and aerosol chemistry
   !!   - Dry and wet deposition processes
   !! - Transforms computed results back to NUOPC field format
   !! - Exports chemistry fields to other model components
   !! - Updates diagnostic outputs and logging
   !!
   !! The routine handles both sequential and parallel execution,
   !! with appropriate logging and error checking throughout.
   !!
   !! \note This routine is called repeatedly during model integration
   !!       and must maintain consistent state between calls
   !!
   !! \ingroup catchem_nuopc_group
   subroutine ModelAdvance(model, rc)
      type(ESMF_GridComp)  :: model
      integer, intent(out) :: rc

      type(ESMF_State) :: importState, exportState
      type(ESMF_Clock) :: clock
      type(ESMF_Time) :: currTime
      type(ESMF_TimeInterval) :: timeStep
      type(CATChem_InternalState) :: is
      character(len=*), parameter :: routine = 'ModelAdvance'
      character(len=512) :: errmsg
      integer :: localPet
      real(ESMF_KIND_R8) :: dt_seconds

      rc = ESMF_SUCCESS

      ! Get component information
      call ESMF_GridCompGet(model, localPet=localPet, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Get states and clock
      call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
         exportState=exportState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Get current time and time step
      call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      call ESMF_TimeIntervalGet(timeStep, s_r8=dt_seconds, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      if (localPet == 0) then
         call ESMF_LogWrite("CATChem: Running CATChem for dt = " // &
            trim(adjustl(real_to_string(dt_seconds))) // " seconds", &
            ESMF_LOGMSG_INFO, rc=rc)
      end if

      ! -- get component's internal state
      call ESMF_GridCompGetInternalState(model, is, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__,  file=__FILE__))  return  ! bail out

      ! Import meteorological data from other components
      call transform_nuopc_to_catchem(is%wrap, importState, currTime, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Run CATChem processes with current time
      call catchem_nuopc_run(is%wrap, dt_seconds, currTime, errmsg, rc)
      if (rc /= ESMF_SUCCESS) then
         call ESMF_LogWrite("CATChem: Failed to run CATChem - " // trim(errmsg), &
            ESMF_LOGMSG_ERROR, rc=rc)
         rc = ESMF_FAILURE
         return
      end if

      ! Export results to other components
      call transform_catchem_to_nuopc(is%wrap, exportState, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Log successful completion
      if (localPet == 0) then
         call ESMF_LogWrite("CATChem: Completed "//routine, ESMF_LOGMSG_INFO, rc=rc)
      end if

   end subroutine ModelAdvance

   !> \brief Finalize the CATChem model component
   !!
   !! This routine performs cleanup operations and finalizes all CATChem
   !! processes, deallocating memory and closing any open resources.
   !!
   !! \param[inout] model NUOPC model component to finalize
   !! \param[out] rc ESMF return code (ESMF_SUCCESS on success)
   !!
   !! \details
   !! This routine performs the following cleanup operations:
   !! - Finalizes all CATChem chemistry and emission processes
   !! - Deallocates memory for all state containers
   !! - Closes any open files or external resources
   !! - Performs diagnostic output if enabled
   !! - Logs completion status and any warnings
   !!
   !! This is a standard NUOPC finalization phase that should be called
   !! when the component is no longer needed or at the end of simulation.
   !!
   !! \note Failure to call this routine may result in memory leaks
   !!       or incomplete output files
   !!
   !! \ingroup catchem_nuopc_group
   subroutine ModelFinalize(model, rc)
      type(ESMF_GridComp)  :: model
      integer, intent(out) :: rc

      type(CATChem_InternalState) :: is
      character(len=*), parameter :: routine = 'ModelFinalize'
      character(len=512) :: errmsg
      integer :: localPet

      rc = ESMF_SUCCESS

      ! Get component information
      call ESMF_GridCompGet(model, localPet=localPet, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! -- get component's internal state
      call ESMF_GridCompGetInternalState(model, is, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__,  file=__FILE__))  return  ! bail out

      ! Finalize CATChem using the interface
      call catchem_nuopc_finalize(is%wrap, rc, errmsg)
      if (rc /= ESMF_SUCCESS) then
         call ESMF_LogWrite("CATChem: Warning - " // trim(errmsg), &
            ESMF_LOGMSG_WARNING, rc=rc)
      end if

      ! Log successful completion
      if (localPet == 0) then
         call ESMF_LogWrite("CATChem: Completed "//routine, ESMF_LOGMSG_INFO, rc=rc)
      end if

      rc = ESMF_SUCCESS

   end subroutine ModelFinalize

   !> \brief Convert real number to string for logging purposes
   !!
   !! This utility function converts a real number to a formatted string
   !! representation suitable for logging and diagnostic messages.
   !!
   !! \param[in] val Real value to convert
   !! \return str String representation of the value
   !!
   !! \details
   !! The function formats the real number with appropriate precision
   !! and removes leading/trailing whitespace for clean output in log
   !! messages and diagnostic information.
   !!
   !! \ingroup catchem_nuopc_group
   function real_to_string(val) result(str)
      real(ESMF_KIND_R8), intent(in) :: val
      character(len=32) :: str

      write(str, '(f0.2)') val
      str = adjustl(str)
   end function real_to_string

end module cc_nuopc
