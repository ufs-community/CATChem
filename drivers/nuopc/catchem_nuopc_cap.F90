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

   use ESMF
   use NUOPC
   use NUOPC_Model, only: &
      NUOPC_ModelGet, &
      SetVM, &
      modelSS        => SetServices, &
      model_label_Advertise       => label_Advertise,      &
      model_label_DataInitialize  => label_DataInitialize, &
      model_label_Advance => label_Advance, &
      model_label_CheckImport => label_CheckImport, &
      model_label_SetRunClock => label_SetRunClock, &
      model_label_Finalize        => label_Finalize

   use catchem_nuopc_interface
   use, intrinsic :: iso_c_binding, only : c_int, c_size_t

   implicit none

   private

   public :: SetServices
   public :: SetVM !for GCAFS only

   !> \brief Component configuration parameters
   !! \{
   character(len=256), save :: config_file = 'CATChem_new_config.yml' !< Configuration file path
   character(len=256), save :: field_mapping_file = 'CATChem_field_mapping.yml' !< Field mapping file path
   !logical, save :: do_chemistry = .true.                         !< Enable chemistry calculations
   !! \}

   ! --- glibc allocator hooks for memory-leak diagnosis / mitigation (Linux/glibc) ---
   ! malloc_trim(0): reclaim freed arena top back to the OS (candidate FIX).
   ! mallinfo():    query live/free heap bytes (the decisive real-leak vs fragmentation
   !                signal). Classic mallinfo links on ALL glibc versions; its int fields
   !                are widened as unsigned (0..4GB) by cc_u32 below. If a given
   !                compiler/glibc rejects the struct-return binding, rebuild with
   !                -DCATCHEM_DISABLE_MALLINFO to keep malloc_trim + rss logging only.
#ifndef CATCHEM_DISABLE_MALLINFO
   type, bind(C) :: cc_mallinfo_t
      integer(c_int) :: arena     ! non-mmapped space allocated from system (bytes)
      integer(c_int) :: ordblks   ! number of free chunks
      integer(c_int) :: smblks    ! number of fastbin blocks
      integer(c_int) :: hblks     ! number of mmapped regions
      integer(c_int) :: hblkhd    ! space in mmapped regions (bytes)
      integer(c_int) :: usmblks   ! always 0, unused
      integer(c_int) :: fsmblks   ! space in freed fastbin blocks (bytes)
      integer(c_int) :: uordblks  ! total allocated (in-use) space (bytes) -- LIVE HEAP
      integer(c_int) :: fordblks  ! total free space (bytes)
      integer(c_int) :: keepcost  ! top-most, releasable (via malloc_trim) space (bytes)
   end type cc_mallinfo_t
#endif

   interface
      function cc_c_malloc_trim(pad) bind(C, name="malloc_trim") result(res)
         import :: c_int, c_size_t
         integer(c_size_t), value :: pad
         integer(c_int) :: res
      end function cc_c_malloc_trim
#ifndef CATCHEM_DISABLE_MALLINFO
      function cc_c_mallinfo() bind(C, name="mallinfo") result(mi)
         import :: cc_mallinfo_t
         type(cc_mallinfo_t) :: mi
      end function cc_c_mallinfo
#endif
   end interface

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
      type(ESMF_Time)           :: startTime, stopTime
      type(ESMF_TimeInterval)   :: timeStep
      real(ESMF_KIND_R8), dimension(:,:), pointer :: coord
      real(ESMF_KIND_R8), dimension(:,:), allocatable :: lon
      real(ESMF_KIND_R8), dimension(:,:), allocatable :: lat
      type(ESMF_CoordSys_Flag)   :: coordSys
      real(ESMF_KIND_R8), parameter :: rad_to_deg = 180._ESMF_KIND_R8 / 3.14159265358979323846_ESMF_KIND_R8
      real(ESMF_KIND_R8) :: convet_unit
      integer :: localPet, petCount
      integer :: item, coord_item, rank, localDeCount, numLevels, localDe, localrc, stat
      integer, dimension(2) :: lb, ub
      logical :: has_tracer_array

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

      ! CAPMEM leak instrumentation (raw absolute VmRSS at cap boundaries).
      ! Gated by env CATCHEM_MEM_CAP=N (print every N ModelAdvance calls; unset/<=0 = off).
      integer, save :: cap_stride = -2
      integer, save :: cap_ncall = 0
      character(len=32) :: cap_env
      integer :: cap_len, cap_stat, cap_ios
      integer :: cap_r0, cap_r1, cap_r2, cap_r3
      logical :: cap_on

      ! malloc_trim / mallinfo instrumentation (independent of CAPMEM).
      ! CATCHEM_MALLOC_TRIM=N  -> call malloc_trim(0) every N ModelAdvance calls (N=1 => every step).
      ! CATCHEM_MALLOC_STATS=N -> print rss + mallinfo line every N calls.
      integer, save :: mt_stride = -2
      integer, save :: ms_stride = -2
      integer, save :: mm_ncall = 0
      character(len=32) :: mm_env
      integer :: mm_len, mm_stat, mm_ios
      integer :: mm_rss_pre, mm_rss_post, mm_trimret
      logical :: mt_on, ms_on
#ifndef CATCHEM_DISABLE_MALLINFO
      type(cc_mallinfo_t) :: mi
#endif

      ! /proc/self/maps + VmData/VmSize instrumentation (locate mmap growth; glibc heap ruled out).
      ! CATCHEM_MAPS_STATS=N -> print maps region-count + anon/file virtual size + VmData/VmSize every N calls.
      integer, save :: mp_stride = -2
      integer, save :: mp_ncall = 0
      character(len=32) :: mp_env
      integer :: mp_len, mp_stat, mp_ios
      integer :: mp_rss, mp_data, mp_size
      integer(8) :: mp_nmaps, mp_anon_kb, mp_file_kb
      logical :: mp_on

      ! /proc/self/smaps_rollup + top-Rss-region instrumentation (locate the progressively-
      ! faulting mapping and whether it is anon/dirty (unreclaimable leak) vs clean file cache).
      ! CATCHEM_SMAPS_STATS=N -> print resident breakdown + largest-Rss region name every N calls.
      integer, save :: sp_stride = -2
      integer, save :: sp_ncall = 0
      character(len=32) :: sp_env
      integer :: sp_len, sp_stat, sp_ios
      integer(8) :: sp_rss, sp_anon, sp_pdirty, sp_pclean, sp_sclean, sp_swap, sp_toprss
      character(len=256) :: sp_topname
      logical :: sp_on

      rc = ESMF_SUCCESS

      if (cap_stride == -2) then
         call get_environment_variable('CATCHEM_MEM_CAP', cap_env, cap_len, cap_stat)
         if (cap_stat == 0 .and. cap_len > 0) then
            read(cap_env, *, iostat=cap_ios) cap_stride
            if (cap_ios /= 0) cap_stride = -1
         else
            cap_stride = -1
         end if
      end if
      cap_on = (cap_stride > 0)
      cap_r0 = -1; cap_r1 = -1; cap_r2 = -1; cap_r3 = -1
      if (cap_on) cap_r0 = cc_cap_read_vmrss_kb()   ! cap entry (after host/mediator ran)

      if (mt_stride == -2) then
         call get_environment_variable('CATCHEM_MALLOC_TRIM', mm_env, mm_len, mm_stat)
         if (mm_stat == 0 .and. mm_len > 0) then
            read(mm_env, *, iostat=mm_ios) mt_stride
            if (mm_ios /= 0) mt_stride = -1
         else
            mt_stride = -1
         end if
      end if
      if (ms_stride == -2) then
         call get_environment_variable('CATCHEM_MALLOC_STATS', mm_env, mm_len, mm_stat)
         if (mm_stat == 0 .and. mm_len > 0) then
            read(mm_env, *, iostat=mm_ios) ms_stride
            if (mm_ios /= 0) ms_stride = -1
         else
            ms_stride = -1
         end if
      end if
      mt_on = (mt_stride > 0)
      ms_on = (ms_stride > 0)

      if (mp_stride == -2) then
         call get_environment_variable('CATCHEM_MAPS_STATS', mp_env, mp_len, mp_stat)
         if (mp_stat == 0 .and. mp_len > 0) then
            read(mp_env, *, iostat=mp_ios) mp_stride
            if (mp_ios /= 0) mp_stride = -1
         else
            mp_stride = -1
         end if
      end if
      mp_on = (mp_stride > 0)

      if (sp_stride == -2) then
         call get_environment_variable('CATCHEM_SMAPS_STATS', sp_env, sp_len, sp_stat)
         if (sp_stat == 0 .and. sp_len > 0) then
            read(sp_env, *, iostat=sp_ios) sp_stride
            if (sp_ios /= 0) sp_stride = -1
         else
            sp_stride = -1
         end if
      end if
      sp_on = (sp_stride > 0)

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
#ifdef CATCHEM_TRACE_NUOPC
      call ESMF_TraceRegionEnter("transform_nuopc_to_catchem", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#endif
      call transform_nuopc_to_catchem(is%wrap, importState, currTime, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#ifdef CATCHEM_TRACE_NUOPC
      call ESMF_TraceRegionExit("transform_nuopc_to_catchem", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#endif

      if (cap_on) cap_r1 = cc_cap_read_vmrss_kb()   ! after import transform

      ! Run CATChem processes with current time
#ifdef CATCHEM_TRACE_NUOPC
      call ESMF_TraceRegionEnter("catchem_nuopc_run", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#endif
      call catchem_nuopc_run(is%wrap, dt_seconds, currTime, errmsg, rc)
      if (rc /= ESMF_SUCCESS) then
         call ESMF_LogWrite("CATChem: Failed to run CATChem - " // trim(errmsg), &
            ESMF_LOGMSG_ERROR, rc=rc)
         rc = ESMF_FAILURE
         return
      end if
#ifdef CATCHEM_TRACE_NUOPC
      call ESMF_TraceRegionExit("catchem_nuopc_run", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#endif

      if (cap_on) cap_r2 = cc_cap_read_vmrss_kb()   ! after catchem_nuopc_run

      ! Export results to other components
#ifdef CATCHEM_TRACE_NUOPC
      call ESMF_TraceRegionEnter("transform_catchem_to_nuopc", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#endif
      call transform_catchem_to_nuopc(is%wrap, exportState, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#ifdef CATCHEM_TRACE_NUOPC
      call ESMF_TraceRegionExit("transform_catchem_to_nuopc", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#endif

      if (cap_on) then
         cap_r3 = cc_cap_read_vmrss_kb()   ! after export transform (cap exit)
         cap_ncall = cap_ncall + 1
         if (mod(cap_ncall, cap_stride) == 0) then
            write(*,'(A,I0,4(A,I0))') '[CATChem CAPMEM] step=', cap_ncall, &
               '  entry=', cap_r0, '  imp=', cap_r1, '  run=', cap_r2, '  exp=', cap_r3
            flush(6)
         end if
      end if

      ! --- malloc_trim mitigation + mallinfo live-heap diagnosis (end of step) ---
      ! Placed AFTER export so it reclaims/measures all of this step's CATChem heap churn.
      ! Interpretation across steps:
      !   inuse_kB (uordblks) grows ~linearly  => GENUINE never-freed leak (trim/buffers won't help).
      !   inuse_kB flat & rss_post flat after trim => reclaimable arena => ship CATCHEM_MALLOC_TRIM=1.
      !   inuse_kB flat but rss_post still grows  => fragmentation (trapped free blocks) => reduce churn.
      if (mt_on .or. ms_on) then
         mm_ncall = mm_ncall + 1
         mm_rss_pre  = cc_cap_read_vmrss_kb()
         mm_trimret  = -1
         if (mt_on .and. mod(mm_ncall, mt_stride) == 0) then
            mm_trimret = int(cc_c_malloc_trim(0_c_size_t))
         end if
         mm_rss_post = cc_cap_read_vmrss_kb()
         if (ms_on .and. mod(mm_ncall, ms_stride) == 0) then
#ifndef CATCHEM_DISABLE_MALLINFO
            mi = cc_c_mallinfo()
            write(*,'(A,I0,8(A,I0))') '[CATChem MALLOC] step=', mm_ncall, &
               ' rss_pre_kB=', mm_rss_pre, ' rss_post_kB=', mm_rss_post, &
               ' trim_released=', mm_trimret, &
               ' arena_kB=', int(cc_u32(mi%arena)/1024_8), &
               ' mmap_kB=', int(cc_u32(mi%hblkhd)/1024_8), &
               ' inuse_kB=', int(cc_u32(mi%uordblks)/1024_8), &
               ' free_kB=', int(cc_u32(mi%fordblks)/1024_8), &
               ' keepcost_kB=', int(cc_u32(mi%keepcost)/1024_8)
#else
            write(*,'(A,I0,3(A,I0))') '[CATChem MALLOC] step=', mm_ncall, &
               ' rss_pre_kB=', mm_rss_pre, ' rss_post_kB=', mm_rss_post, &
               ' trim_released=', mm_trimret
#endif
            flush(6)
         end if
      end if

      ! --- /proc/self/maps + VmData/VmSize diagnosis (glibc heap ruled out => locate mmap growth) ---
      ! nmaps rising          => new mmap regions each step (un-freed mappings; find the caller).
      ! anon_kB rising        => anonymous virtual (heap / secondary glibc arena / direct anon mmap).
      ! file_kB rising        => file-backed mapping growth (NetCDF/HDF5/library-mapped file).
      ! vmsize flat, vmrss up => progressive faulting of an already-mapped (pre-reserved) region.
      if (mp_on) then
         mp_ncall = mp_ncall + 1
         if (mod(mp_ncall, mp_stride) == 0) then
            call cc_cap_read_status(mp_rss, mp_data, mp_size)
            call cc_cap_read_maps(mp_nmaps, mp_anon_kb, mp_file_kb)
            write(*,'(A,I0,6(A,I0))') '[CATChem MAPS] step=', mp_ncall, &
               ' vmrss_kB=', mp_rss, ' vmsize_kB=', mp_size, ' vmdata_kB=', mp_data, &
               ' nmaps=', mp_nmaps, ' anon_kB=', mp_anon_kb, ' file_kB=', mp_file_kb
            flush(6)
         end if
      end if

      ! --- /proc/self/smaps_rollup + largest-Rss region (identify the faulting mapping) ---
      ! anon_kB or pdirty_kB rising ~ vmrss => UNRECLAIMABLE (anonymous/dirty) leak -> real OOM cause.
      ! pclean_kB/sclean_kB rising, anon flat => reclaimable file-backed page cache (node pressure).
      ! top=<name> names the single mapping holding the most resident memory (the prime suspect).
      if (sp_on) then
         sp_ncall = sp_ncall + 1
         if (mod(sp_ncall, sp_stride) == 0) then
            call cc_cap_read_smaps_rollup(sp_rss, sp_anon, sp_pdirty, sp_pclean, sp_sclean, sp_swap)
            call cc_cap_top_rss_region(sp_toprss, sp_topname)
            write(*,'(8(A,I0),2A)') '[CATChem SMAPS] step=', sp_ncall, &
               ' rss_kB=', sp_rss, ' anon_kB=', sp_anon, ' pdirty_kB=', sp_pdirty, &
               ' pclean_kB=', sp_pclean, ' sclean_kB=', sp_sclean, ' swap_kB=', sp_swap, &
               ' topRss_kB=', sp_toprss, ' top=', trim(sp_topname)
            flush(6)
         end if
      end if

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

      ! Deallocate internal state wrapper
      if (associated(is%wrap)) then
         deallocate(is%wrap)
         nullify(is%wrap)
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

   !> \brief Read current resident set size (VmRSS, kB) from /proc/self/status.
   !! Returns -1 if unavailable. Used only by the CAPMEM leak instrumentation.
   integer function cc_cap_read_vmrss_kb() result(kb)
      integer :: u, ios
      character(len=256) :: line
      kb = -1
      open(newunit=u, file='/proc/self/status', status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         if (line(1:6) == 'VmRSS:') then
            read(line(7:), *, iostat=ios) kb
            if (ios /= 0) kb = -1
            exit
         end if
      end do
      close(u)
   end function cc_cap_read_vmrss_kb

#ifndef CATCHEM_DISABLE_MALLINFO
   !> \brief Widen a (possibly negative) 32-bit mallinfo field to an unsigned 0..4GB value.
   integer(8) function cc_u32(v) result(r)
      integer(c_int), intent(in) :: v
      r = int(v, 8)
      if (r < 0) r = r + 4294967296_8
   end function cc_u32
#endif

   !> \brief Read VmRSS, VmData, VmSize (kB) from /proc/self/status in one pass. -1 if missing.
   subroutine cc_cap_read_status(rss_kb, data_kb, size_kb)
      integer, intent(out) :: rss_kb, data_kb, size_kb
      integer :: u, ios
      character(len=256) :: line
      rss_kb = -1; data_kb = -1; size_kb = -1
      open(newunit=u, file='/proc/self/status', status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         if (line(1:6) == 'VmRSS:') then
            read(line(7:), *, iostat=ios) rss_kb;  if (ios /= 0) rss_kb = -1
         else if (line(1:7) == 'VmData:') then
            read(line(8:), *, iostat=ios) data_kb; if (ios /= 0) data_kb = -1
         else if (line(1:7) == 'VmSize:') then
            read(line(8:), *, iostat=ios) size_kb; if (ios /= 0) size_kb = -1
         end if
      end do
      close(u)
   end subroutine cc_cap_read_status

   !> \brief Parse /proc/self/maps: count regions (nmaps) and sum virtual size (kB) into
   !! anonymous vs file-backed buckets. A region is file-backed if its line contains a '/'
   !! pathname; anonymous otherwise (covers plain anon, [heap], [stack], [anon:*]).
   subroutine cc_cap_read_maps(nmaps, anon_kb, file_kb)
      integer(8), intent(out) :: nmaps, anon_kb, file_kb
      integer :: u, ios, dashp, sp
      character(len=512) :: line
      integer(8) :: a0, a1, sz
      nmaps = 0; anon_kb = 0; file_kb = 0
      open(newunit=u, file='/proc/self/maps', status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         nmaps = nmaps + 1
         dashp = index(line, '-')
         if (dashp <= 1) cycle
         a0 = cc_hex2i(line(1:dashp-1))
         sp = index(line(dashp+1:), ' ')
         if (sp <= 1) cycle
         a1 = cc_hex2i(line(dashp+1:dashp+sp-1))
         sz = (a1 - a0) / 1024_8
         if (index(line, '/') > 0) then
            file_kb = file_kb + sz
         else
            anon_kb = anon_kb + sz
         end if
      end do
      close(u)
   end subroutine cc_cap_read_maps

   !> \brief Convert a hex string (no 0x prefix) to integer(8). Non-hex chars are skipped.
   integer(8) function cc_hex2i(s) result(v)
      character(len=*), intent(in) :: s
      integer :: i, d
      character :: c
      v = 0_8
      do i = 1, len_trim(s)
         c = s(i:i)
         select case (c)
         case ('0':'9'); d = ichar(c) - ichar('0')
         case ('a':'f'); d = ichar(c) - ichar('a') + 10
         case ('A':'F'); d = ichar(c) - ichar('A') + 10
         case default;   cycle
         end select
         v = v * 16_8 + int(d, 8)
      end do
   end function cc_hex2i

   !> \brief Read aggregate resident breakdown (kB) from /proc/self/smaps_rollup. -1 if missing.
   !! Distinguishes anonymous/dirty (unreclaimable) from clean file-backed (reclaimable) RSS.
   subroutine cc_cap_read_smaps_rollup(rss, anon, pdirty, pclean, sclean, swap)
      integer(8), intent(out) :: rss, anon, pdirty, pclean, sclean, swap
      integer :: u, ios
      character(len=256) :: line
      rss = -1; anon = -1; pdirty = -1; pclean = -1; sclean = -1; swap = -1
      open(newunit=u, file='/proc/self/smaps_rollup', status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         if (line(1:4) == 'Rss:') then
            read(line(5:), *, iostat=ios) rss
         else if (line(1:10) == 'Anonymous:') then
            read(line(11:), *, iostat=ios) anon
         else if (line(1:14) == 'Private_Dirty:') then
            read(line(15:), *, iostat=ios) pdirty
         else if (line(1:14) == 'Private_Clean:') then
            read(line(15:), *, iostat=ios) pclean
         else if (line(1:13) == 'Shared_Clean:') then
            read(line(14:), *, iostat=ios) sclean
         else if (line(1:5) == 'Swap:') then
            read(line(6:), *, iostat=ios) swap
         end if
      end do
      close(u)
   end subroutine cc_cap_read_smaps_rollup

   !> \brief Scan /proc/self/smaps for the single mapping with the largest Rss; return its
   !! Rss (kB) and a short name ('/'-path, [bracket] region, or [anon]). -1 if unavailable.
   subroutine cc_cap_top_rss_region(top_rss_kb, top_name)
      integer(8), intent(out) :: top_rss_kb
      character(len=*), intent(out) :: top_name
      integer :: u, ios, sl, br
      character(len=512) :: line
      character(len=256) :: curname
      character :: c
      integer(8) :: rss
      top_rss_kb = -1; top_name = ''
      curname = '[anon]'
      open(newunit=u, file='/proc/self/smaps', status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         c = line(1:1)
         if ((c >= '0' .and. c <= '9') .or. (c >= 'a' .and. c <= 'f')) then
            ! region header line: derive a short name for this mapping
            sl = index(line, '/')
            br = index(line, '[')
            if (sl > 0) then
               curname = adjustl(line(sl:))
            else if (br > 0) then
               curname = adjustl(line(br:))
            else
               curname = '[anon]'
            end if
         else if (line(1:4) == 'Rss:') then
            read(line(5:), *, iostat=ios) rss
            if (ios == 0 .and. rss > top_rss_kb) then
               top_rss_kb = rss
               top_name = trim(curname)
            end if
         end if
      end do
      close(u)
   end subroutine cc_cap_top_rss_region

end module cc_nuopc
