!> \file nuopc_contract_harness.F90
!! \brief Local DataATM-style NUOPC contract harness for CATChem.
!!
!! Builds only when CATCHEM_BUILD_NUOPC=ON. Reproduces the ursa
!! run-phase mechanism in-process: load the real field mapping
!! (CATChem_field_mapping.yml), create an ESMF field for every
!! required import entry, and drive transform_nuopc_to_catchem —
!! asserting success and the Z0 cm->m conversion landing in the met
!! state (the 2026-08-11 first-timestep abort regression).
module dataatm_services
   use ESMF
   use NUOPC
   use NUOPC_Model, only: modelSS => SetServices, &
      model_label_Advertise => label_Advertise, &
      model_label_RealizeProvided => label_RealizeProvided, &
      model_label_Advance => label_Advance, &
      model_label_CheckImport => label_CheckImport, &
      NUOPC_ModelGet
   use catchem_nuopc_interface, only: load_field_config, field_config
   implicit none
contains
   subroutine SetServices(model, rc)
      type(ESMF_GridComp) :: model
      integer, intent(out) :: rc
      rc = ESMF_SUCCESS
      call NUOPC_CompDerive(model, modelSS, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call NUOPC_CompSpecialize(model, specLabel=model_label_Advertise, &
         specRoutine=Advertise, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call NUOPC_CompSpecialize(model, specLabel=model_label_RealizeProvided, &
         specRoutine=Realize, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
         specRoutine=Advance, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      ! The mock consumes output passively.  Its Advance specialization does
      ! not use CATChem fields as an input to a forecast calculation, so the
      ! generic template's pre-advance timestamp gate is not applicable.
      call NUOPC_CompSpecialize(model, specLabel=model_label_CheckImport, &
         specRoutine=CheckImport, rc=rc)
   end subroutine SetServices

   subroutine Advertise(model, rc)
      type(ESMF_GridComp) :: model
      integer, intent(out) :: rc
      type(ESMF_State) :: importState, exportState
      integer :: i
      character(len=256) :: errmsg
      rc = ESMF_SUCCESS
      call NUOPC_ModelGet(model, importState=importState, exportState=exportState, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call load_field_config('CATChem_field_mapping.yml', rc, errmsg)
      if (rc /= ESMF_SUCCESS) return
      do i = 1, field_config%n_import_fields
         if (field_config%import_fields(i)%optional .and. &
            .not. field_config%import_fields(i)%advertise) cycle
         call NUOPC_Advertise(exportState, &
            StandardName=trim(field_config%import_fields(i)%standard_name), &
            TransferOfferGeomObject='will provide', SharePolicyField='share', rc=rc)
         if (rc /= ESMF_SUCCESS) return
      end do
      ! Consume CATChem's advertised diagnostics and tracer output as a
      ! DataATM analogue.  This makes the CATChem export state real during
      ! RunPhase rather than leaving its fields as uncreated offers.
      do i = 1, field_config%n_export_fields
         call NUOPC_Advertise(importState, &
            StandardName=trim(field_config%export_fields(i)%standard_name), &
            TransferOfferGeomObject='will provide', SharePolicyField='share', rc=rc)
         if (rc /= ESMF_SUCCESS) return
      end do
   end subroutine Advertise

   subroutine CheckImport(model, rc)
      type(ESMF_GridComp) :: model
      integer, intent(out) :: rc
      ! ESMF specialization signature; this no-op never touches the component.
      associate(unused_model => model); end associate
      rc = ESMF_SUCCESS
   end subroutine CheckImport

   subroutine Advance(model, rc)
      type(ESMF_GridComp) :: model
      integer, intent(out) :: rc
      type(ESMF_State) :: exportState
      type(ESMF_Clock) :: clock
      type(ESMF_Time) :: curr_time, next_time
      type(ESMF_TimeInterval) :: time_step
      type(ESMF_Field), pointer :: fields(:)
      type(ESMF_Field) :: field
      character(len=ESMF_MAXSTR) :: field_name
      integer :: i, rank
      rc = ESMF_SUCCESS
      call NUOPC_ModelGet(model, exportState=exportState, modelClock=clock, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call ESMF_ClockGet(clock, currTime=curr_time, timeStep=time_step, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      next_time = curr_time + time_step
      ! The driver validates consumer imports at its current coupling time.
      ! Stamp before inspecting a member list: NUOPC may represent
      ! an empty or nested export state without returning a local list.
      call NUOPC_SetTimestamp(exportState, curr_time, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      nullify(fields)
      call NUOPC_GetStateMemberLists(exportState, fieldList=fields, nestedFlag=.true., rc=rc)
      if (rc /= ESMF_SUCCESS) return
      if (.not. associated(fields)) return
      do i = 1, size(fields)
         field = fields(i)
         call ESMF_FieldGet(field, name=field_name, rank=rank, rc=rc)
         if (rc /= ESMF_SUCCESS) return
         call populate_mock_field(field, rc)
         if (rc /= ESMF_SUCCESS) return
         call NUOPC_SetTimestamp(field, curr_time, rc=rc)
         if (rc /= ESMF_SUCCESS) return
      end do
      ! Preserve the driver's current coupling timestamp after updating fields.
      call NUOPC_SetTimestamp(exportState, curr_time, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      deallocate(fields)
      nullify(fields)
   end subroutine Advance

   subroutine populate_mock_field(field, rc)
      type(ESMF_Field), intent(inout) :: field
      integer, intent(out) :: rc
      real(ESMF_KIND_R8), pointer :: f2(:,:), f3(:,:,:), f4(:,:,:,:)
      character(len=ESMF_MAXSTR) :: field_name
      integer :: k, rank

      rc = ESMF_SUCCESS
      call ESMF_FieldGet(field, name=field_name, rank=rank, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      select case (rank)
       case (2)
         nullify(f2); call ESMF_FieldGet(field, farrayPtr=f2, rc=rc)
         if (rc /= ESMF_SUCCESS) return
         f2 = 1.0_ESMF_KIND_R8
       case (3)
         nullify(f3); call ESMF_FieldGet(field, farrayPtr=f3, rc=rc)
         if (rc /= ESMF_SUCCESS) return
         select case (trim(field_name))
          case ('inst_pres_interface')
            do k = 1, size(f3, 3)
               f3(:, :, k) = 100000.0_ESMF_KIND_R8 - 10000.0_ESMF_KIND_R8 * real(k - 1, ESMF_KIND_R8)
            end do
          case ('inst_pres_levels')
            do k = 1, size(f3, 3)
               f3(:, :, k) = 95000.0_ESMF_KIND_R8 - 10000.0_ESMF_KIND_R8 * real(k - 1, ESMF_KIND_R8)
            end do
          case ('inst_geop_interface')
            do k = 1, size(f3, 3)
               f3(:, :, k) = 1000.0_ESMF_KIND_R8 * real(k - 1, ESMF_KIND_R8)
            end do
          case ('inst_geop_levels')
            do k = 1, size(f3, 3)
               f3(:, :, k) = 500.0_ESMF_KIND_R8 + 1000.0_ESMF_KIND_R8 * real(k - 1, ESMF_KIND_R8)
            end do
          case ('inst_temp_levels')
            f3 = 300.0_ESMF_KIND_R8
          case default
            f3 = 1.0_ESMF_KIND_R8
         end select
       case (4)
         nullify(f4); call ESMF_FieldGet(field, farrayPtr=f4, rc=rc)
         if (rc /= ESMF_SUCCESS) return
         ! Keep host-owned and pseudo-tracer slots distinct from transported
         ! chemistry values so both exchange directions can be checked.
         f4 = 0.0_ESMF_KIND_R8
         f4(:,:,:,1) = 1.0E-2_ESMF_KIND_R8
         f4(:,:,:,2) = 1.0E-6_ESMF_KIND_R8
         f4(:,:,:,3) = 2.0E-6_ESMF_KIND_R8
         f4(:,:,:,4) = 7.0E-6_ESMF_KIND_R8
         f4(:,:,:,5) = 8.0E-6_ESMF_KIND_R8
      end select
   end subroutine populate_mock_field

   subroutine Realize(model, rc)
      type(ESMF_GridComp) :: model
      integer, intent(out) :: rc
      type(ESMF_State) :: importState, exportState
      type(ESMF_Grid) :: grid
      type(ESMF_Field) :: field
      type(ESMF_Array) :: array
      type(ESMF_Info) :: tracer_info
      integer :: i, nzf, ntr
      character(len=ESMF_MAXSTR) :: tracer_names(5), tracer_units(5)
      rc = ESMF_SUCCESS
      call NUOPC_ModelGet(model, importState=importState, exportState=exportState, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call ESMF_GridCompGet(model, grid=grid, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      ntr = 5
      tracer_names = [character(len=ESMF_MAXSTR) :: 'sphum', 'chemical_gas', 'chemical_aerosol', 'PM25', 'PM10']
      tracer_units = [character(len=ESMF_MAXSTR) :: '1', 'ppm', 'kg kg^-1', 'ug m-3', 'ug m-3']
      do i = 1, field_config%n_import_fields
         if (field_config%import_fields(i)%optional .and. &
            .not. field_config%import_fields(i)%advertise) cycle
         select case (field_config%import_fields(i)%dimensions)
          case (2)
            field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
               name=trim(field_config%import_fields(i)%standard_name), rc=rc)
          case (3)
            nzf = 5
            if (trim(field_config%import_fields(i)%vertical_axis) == 'interface') nzf = 6
            field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
               ungriddedLBound=(/1/), ungriddedUBound=(/nzf/), &
               name=trim(field_config%import_fields(i)%standard_name), rc=rc)
          case (4)
            field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
               ungriddedLBound=(/1,1/), ungriddedUBound=(/5,ntr/), &
               name=trim(field_config%import_fields(i)%standard_name), rc=rc)
          case default
            cycle
         end select
         if (rc /= ESMF_SUCCESS) return
         if (field_config%import_fields(i)%dimensions == 4) then
            call ESMF_FieldGet(field, array=array, rc=rc)
            if (rc /= ESMF_SUCCESS) return
            call ESMF_InfoGetFromHost(array, tracer_info, rc=rc)
            if (rc /= ESMF_SUCCESS) return
            call ESMF_InfoSet(tracer_info, 'tracerNames', tracer_names, rc=rc)
            if (rc /= ESMF_SUCCESS) return
            call ESMF_InfoSet(tracer_info, 'tracerUnits', tracer_units, rc=rc)
            if (rc /= ESMF_SUCCESS) return
         end if
         ! NUOPC can call a consumer's first RunPhase before the producer's
         ! Advance specialization.  Realized exports must therefore already
         ! contain a valid initial profile.
         call populate_mock_field(field, rc)
         if (rc /= ESMF_SUCCESS) return
         call NUOPC_Realize(exportState, field, rc=rc)
         if (rc /= ESMF_SUCCESS) return
      end do
      ! Realize the receiving half of CATChem's exports on the same grid.
      ! This is what the UFS atmospheric component supplies for chemistry
      ! outputs; without it a shared NUOPC export remains an uncreated offer.
      do i = 1, field_config%n_export_fields
         select case (field_config%export_fields(i)%dimensions)
          case (2)
            field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
               name=trim(field_config%export_fields(i)%standard_name), rc=rc)
          case (4)
            field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
               ungriddedLBound=(/1,1/), ungriddedUBound=(/5,ntr/), &
               name=trim(field_config%export_fields(i)%standard_name), rc=rc)
          case default
            cycle
         end select
         if (rc /= ESMF_SUCCESS) return
         if (field_config%export_fields(i)%dimensions == 4) then
            call ESMF_FieldGet(field, array=array, rc=rc)
            if (rc /= ESMF_SUCCESS) return
            call ESMF_InfoGetFromHost(array, tracer_info, rc=rc)
            if (rc /= ESMF_SUCCESS) return
            call ESMF_InfoSet(tracer_info, 'tracerNames', tracer_names, rc=rc)
            if (rc /= ESMF_SUCCESS) return
            call ESMF_InfoSet(tracer_info, 'tracerUnits', tracer_units, rc=rc)
            if (rc /= ESMF_SUCCESS) return
         end if
         call NUOPC_Realize(importState, field, rc=rc)
         if (rc /= ESMF_SUCCESS) return
      end do
   end subroutine Realize
end module dataatm_services

module contract_driver_services
   use ESMF
   use NUOPC
   use NUOPC_Driver, only: driverSS => SetServices, &
      label_SetModelServices, NUOPC_DriverAddComp
   use NUOPC_Connector, only: connectorSS => SetServices
   use cc_nuopc, only: CATChem_SetServices => SetServices
   use dataatm_services, only: DataATM_SetServices => SetServices
   implicit none
   logical, save :: model_services_called = .false.
contains
   logical function driver_services_seen()
      driver_services_seen = model_services_called
   end function driver_services_seen
   subroutine SetServices(driver, rc)
      type(ESMF_GridComp) :: driver
      integer, intent(out) :: rc
      rc = ESMF_SUCCESS
      call NUOPC_CompDerive(driver, driverSS, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call NUOPC_CompSpecialize(driver, specLabel=label_SetModelServices, &
         specRoutine=SetModelServices, rc=rc)
      if (rc /= ESMF_SUCCESS) return
   end subroutine SetServices

   subroutine SetModelServices(driver, rc)
      type(ESMF_GridComp) :: driver
      integer, intent(out) :: rc
      type(ESMF_GridComp) :: catchem
      type(ESMF_GridComp) :: dataatm
      type(ESMF_CplComp) :: connector
      type(ESMF_Grid) :: data_grid
      type(ESMF_Time) :: start_time, stop_time
      type(ESMF_TimeInterval) :: time_step
      type(ESMF_Clock) :: clock
      rc = ESMF_SUCCESS
      model_services_called = .true.
      call NUOPC_DriverAddComp(driver, "CATChem", &
         compSetServicesRoutine=CATChem_SetServices, comp=catchem, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call NUOPC_DriverAddComp(driver, "DataATM", &
         compSetServicesRoutine=DataATM_SetServices, comp=dataatm, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      data_grid = ESMF_GridCreateNoPeriDim(maxIndex=(/4, 2/), rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call ESMF_GridAddCoord(data_grid, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call ESMF_GridCompSet(dataatm, grid=data_grid, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call ESMF_GridCompSet(catchem, grid=data_grid, rc=rc)
      if (rc /= ESMF_SUCCESS) return

      call NUOPC_DriverAddComp(driver, srcCompLabel="DataATM", dstCompLabel="CATChem", &
         compSetServicesRoutine=connectorSS, comp=connector, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call NUOPC_DriverAddComp(driver, srcCompLabel="CATChem", dstCompLabel="DataATM", &
         compSetServicesRoutine=connectorSS, comp=connector, rc=rc)
      if (rc /= ESMF_SUCCESS) return

      call ESMF_TimeIntervalSet(time_step, s=600, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call ESMF_TimeSet(start_time, yy=2021, mm=3, dd=22, h=6, &
         calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call ESMF_TimeSet(stop_time, yy=2021, mm=3, dd=22, h=6, m=20, &
         calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      clock = ESMF_ClockCreate(name="DataATM-CATChem lifecycle clock", &
         startTime=start_time, stopTime=stop_time, timeStep=time_step, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call ESMF_GridCompSet(driver, clock=clock, rc=rc)
   end subroutine SetModelServices
end module contract_driver_services

program nuopc_contract_harness
   use iso_c_binding, only: c_ptr, c_char, c_double, c_int, c_associated, c_f_pointer, c_null_char
   use ESMF
   use contract_driver_services, only: ContractDriver_SetServices => SetServices, &
      driver_services_seen
   use catchem_nuopc_interface, only: load_field_config, transform_nuopc_to_catchem, &
      transform_catchem_to_nuopc, field_config, cc_wrap_type, update_pm_diagnostics
   use catchem_bridge_error, only: CC_SUCCESS
   use catchem_bridge_precision, only: fp, exact_equal
   use CATChem_API, only: CATChem_Model

   implicit none

   interface
      type(c_ptr) function catchem_state_get_pointer_2d(state_ptr, name) bind(C, name="catchem_state_get_pointer_2d")
         import :: c_ptr, c_char
         type(c_ptr), value :: state_ptr
         character(kind=c_char), intent(in) :: name(*)
      end function
      integer(c_int) function catchem_state_get_species_count_checked(state, count) &
         bind(C, name="catchem_state_get_species_count_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state
         integer(c_int), intent(out) :: count
      end function
      integer(c_int) function catchem_state_get_species_name_at_checked(state, index, name, length) &
         bind(C, name="catchem_state_get_species_name_at_checked")
         import :: c_ptr, c_int, c_char
         type(c_ptr), value :: state
         integer(c_int), value :: index, length
         character(c_char), intent(out) :: name(*)
      end function
      integer(c_int) function catchem_state_is_species_gas_checked(state, index, value) &
         bind(C, name="catchem_state_is_species_gas_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state
         integer(c_int), value :: index
         integer(c_int), intent(out) :: value
      end function
      integer(c_int) function catchem_state_is_species_aerosol_checked(state, index, value) &
         bind(C, name="catchem_state_is_species_aerosol_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state
         integer(c_int), value :: index
         integer(c_int), intent(out) :: value
      end function
      integer(c_int) function catchem_state_get_species_mw_checked(state, index, mw) &
         bind(C, name="catchem_state_get_species_mw_checked")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state
         integer(c_int), value :: index
         real(c_double), intent(out) :: mw
      end function
   end interface

   integer, parameter :: nx = 4, ny = 2, nz = 5, ntr = 5
   type(cc_wrap_type) :: cc_wrap
   type(CATChem_Model) :: parity_model
   type(ESMF_Grid) :: grid
   type(ESMF_State) :: importState, exportState
   type(ESMF_GridComp) :: driver_comp
   type(ESMF_Field) :: field
   type(ESMF_Array) :: array
   type(ESMF_Info) :: tracer_info
   type(ESMF_Time) :: currTime
   real(ESMF_KIND_R8), pointer :: fptr2(:, :), fptr3(:, :, :), fptr4(:, :, :, :)
   real(c_double), pointer :: z0_ptr(:,:) => null()
   type(c_ptr) :: raw_z0_ptr
   character(len=256) :: errmsg
   integer :: rc, urc, i, nzf, parity_issues, gas_slot, aerosol_slot
   integer(c_int) :: gas_idx, aerosol_idx, gas_flag, aerosol_flag, species_count
   real(c_double) :: gas_mw
   integer(c_int) :: parity_count, parity_status
   character(kind=c_char) :: parity_c_name(64)
   character(len=64) :: parity_name
   character(len=512) :: parity_report
   ! Same length as cc_wrap_type%tracer_map%names/units (character(len=128)).
   character(len=128) :: exchange_tracer_names(ntr), exchange_tracer_units(ntr)
   character(len=64), parameter :: parity_expected(3) = [character(len=64) :: &
      'unfamiliar_alpha', 'unfamiliar_beta', 'unfamiliar_gamma']

   call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, &
      defaultlogfilename="test_nuopc_transform.log", rc=rc)
   call check(rc, "ESMF_Initialize")

   ! CATChem is created exactly once by NUOPC_DriverAddComp in the driver's
   ! SetModelServices specialization.  Creating a detached component here
   ! would give initialization and RunPhase1 different component identities.
   driver_comp = ESMF_GridCompCreate(name="contract_driver", rc=rc)
   call check(rc, "GridCompCreate driver")
   call ESMF_GridCompSetServices(driver_comp, ContractDriver_SetServices, rc=rc)
   call check(rc, "GridCompSetServices driver")
   call ESMF_GridCompInitialize(driver_comp, userRc=urc, rc=rc)
   call check(rc, "GridCompInitialize managed driver")
   call check(urc, "Managed driver user return code")
   call ESMF_GridCompRun(driver_comp, userRc=urc, rc=rc)
   call check(rc, "GridCompRun managed driver")
   call check(urc, "Managed driver run user return code")
   call ESMF_GridCompFinalize(driver_comp, userRc=urc, rc=rc)
   call check(rc, "GridCompFinalize managed driver")
   call check(urc, "Managed driver finalize user return code")
   if (.not. driver_services_seen()) error stop 'Managed driver services were not invoked'
   print *, 'PASS: DataATM -> CATChem managed lifecycle completed'
   ! Keep ESMF initialized while the remainder of this executable drives the
   ! production transforms directly with synthetic exchange states.

   ! 1. Load the real production field mapping
   call load_field_config('CATChem_field_mapping.yml', rc, errmsg)
   if (rc /= 0) then
      print *, 'FAIL: load_field_config: ', trim(errmsg)
      error stop 1
   end if

   ! 2. Initialize the model (facade + host-local grid dims)
   call cc_wrap%catchem_model%initialize('CATChem_new_config.yml', nx, ny, nz, rc=rc)
   if (rc /= 0) then
      print *, 'FAIL: model initialize rc=', rc
      error stop 1
   end if
   cc_wrap%field_config = field_config

   ! Select representative transported species from the active state rather
   ! than baking a particular mechanism (such as GOCART) into this harness.
   parity_status = catchem_state_get_species_count_checked(cc_wrap%catchem_model%state_mgr_ptr, species_count)
   if (parity_status /= 0_c_int) error stop 'Active species count lookup failed'
   gas_slot = 0; aerosol_slot = 0
   gas_idx = 0_c_int; aerosol_idx = 0_c_int
   do i = 1, species_count
      gas_flag = 0_c_int; aerosol_flag = 0_c_int
      if (catchem_state_is_species_gas_checked(cc_wrap%catchem_model%state_mgr_ptr, int(i,c_int), gas_flag) /= 0_c_int) cycle
      if (catchem_state_is_species_aerosol_checked(cc_wrap%catchem_model%state_mgr_ptr, int(i,c_int), aerosol_flag) /= 0_c_int) cycle
      if (gas_flag /= 0_c_int .and. gas_slot == 0) then
         gas_idx = int(i,c_int); gas_slot = 2
      else if (aerosol_flag /= 0_c_int .and. aerosol_slot == 0) then
         aerosol_idx = int(i,c_int); aerosol_slot = 3
      end if
      if (gas_slot > 0 .and. aerosol_slot > 0) exit
   end do
   if (gas_slot == 0 .or. aerosol_slot == 0) error stop 'Active mechanism lacks gas/aerosol representatives'
   parity_status = catchem_state_get_species_name_at_checked(cc_wrap%catchem_model%state_mgr_ptr, gas_idx, parity_c_name, 64_c_int)
   call parity_from_c(parity_c_name, parity_name)
   exchange_tracer_names(2) = trim(parity_name)
   parity_status = catchem_state_get_species_name_at_checked(cc_wrap%catchem_model%state_mgr_ptr, aerosol_idx, parity_c_name, 64_c_int)
   call parity_from_c(parity_c_name, parity_name)
   exchange_tracer_names(3) = trim(parity_name)
   exchange_tracer_names(1) = 'sphum'; exchange_tracer_names(4) = 'PM25'; exchange_tracer_names(5) = 'PM10'
   exchange_tracer_units = 'kg kg-1'; exchange_tracer_units(1) = '1'
   ! UFS advertises gas tracers in ppmv and aerosol tracers in ug/kg.
   exchange_tracer_units(2) = 'ppm'; exchange_tracer_units(3) = 'ug/kg'
   ! PM2.5/PM10 are diagnostic mass concentrations, not transported tracers.
   exchange_tracer_units(4:5) = 'ug m-3'

   ! The direct transform fixture does not pass through the ESMF cap's
   ! catchem_nuopc_init wrapper, so install the same validated descriptor
   ! shape here. Species identity and gas molecular weight remain discovered
   ! from the active state; no mechanism-specific names or ordering are used.
   allocate(cc_wrap%tracer_map%names(ntr), cc_wrap%tracer_map%units(ntr))
   allocate(cc_wrap%tracer_map%nuopc_to_cc(ntr), cc_wrap%tracer_map%entry_kind(ntr))
   allocate(cc_wrap%tracer_map%host_to_catchem(ntr), cc_wrap%tracer_map%catchem_to_host(ntr))
   cc_wrap%tracer_map%names = exchange_tracer_names
   cc_wrap%tracer_map%units = exchange_tracer_units
   cc_wrap%tracer_map%nuopc_to_cc = [0, int(gas_idx), int(aerosol_idx), 0, 0]
   cc_wrap%tracer_map%entry_kind = [0, 1, 1, 2, 2]
   if (catchem_state_get_species_mw_checked(cc_wrap%catchem_model%state_mgr_ptr, gas_idx, gas_mw) /= 0_c_int) &
      error stop 'Representative gas molecular weight lookup failed'
   ! UFS/GOCART representative fields are already in their advertised native
   ! units, so ppm and ug/kg use identity transforms at this boundary.
   cc_wrap%tracer_map%host_to_catchem = 1.0_c_double
   cc_wrap%tracer_map%host_to_catchem(2) = 1.0_c_double
   cc_wrap%tracer_map%host_to_catchem(3) = 1.0_c_double
   cc_wrap%tracer_map%catchem_to_host = 1.0_c_double / cc_wrap%tracer_map%host_to_catchem

   ! Shared host-conformance fixture: mechanism order, failure category, and report shape.
   call parity_model%initialize('host_conformance/CATChem_config.yml', 2, 1, 3, rc=rc)
   if (rc /= 0) error stop 'NUOPC parity fixture initialization failed'
   parity_status = catchem_state_get_species_count_checked(parity_model%state_mgr_ptr, parity_count)
   if (parity_status /= 0_c_int .or. parity_count /= 3_c_int) error stop 'NUOPC mechanism count parity failed'
   do i = 1, 3
      parity_status = catchem_state_get_species_name_at_checked( &
         parity_model%state_mgr_ptr, int(i, c_int), parity_c_name, 64_c_int)
      if (parity_status /= 0_c_int) error stop 'NUOPC mechanism name query failed'
      call parity_from_c(parity_c_name, parity_name)
      if (trim(parity_name) /= trim(parity_expected(i))) error stop 'NUOPC mechanism ordering parity failed'
   end do
   call parity_model%set_physical_validation_policy(99, rc)
   if (rc /= 8) error stop 'NUOPC invalid-policy category parity failed'
   call parity_model%set_physical_validation_policy(0, rc)
   call parity_model%get_physical_validation_report(parity_issues, parity_report, rc)
   if (rc /= 0 .or. parity_issues /= 0) error stop 'NUOPC physical-report parity failed'
   call parity_model%finalize(rc)
   if (rc /= 0) error stop 'NUOPC parity fixture finalization failed'
   call ESMF_TimeIntervalSet(cc_wrap%timeStep, s=600, rc=rc)
   call check(rc, "TimeIntervalSet")

   ! 3. Build an import state holding every REQUIRED import field from
   !    the mapping, sized per its declared dimensionality. Z0
   !    (inst_surface_roughness) carries 150 cm; everything else 1.0.
   grid = ESMF_GridCreateNoPeriDim(maxIndex=(/nx, ny/), rc=rc)
   call check(rc, "GridCreate")
   importState = ESMF_StateCreate(name="import", rc=rc)
   call check(rc, "StateCreate")

   do i = 1, cc_wrap%field_config%n_import_fields
      select case (cc_wrap%field_config%import_fields(i)%dimensions)
       case (2)
         field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
            name=trim(cc_wrap%field_config%import_fields(i)%standard_name), rc=rc)
         call check(rc, "FieldCreate 2d")
         call ESMF_FieldGet(field, farrayPtr=fptr2, rc=rc)
         call check(rc, "FieldGet 2d")
         if (trim(cc_wrap%field_config%import_fields(i)%catchem_var) == 'Z0') then
            fptr2 = 150.0_ESMF_KIND_R8 ! cm; the transform converts to 1.5 m
         else
            block
               integer :: ix, iy
               do iy = 1, ny
                  do ix = 1, nx
                     fptr2(ix, iy) = real(ix + 10 * iy, ESMF_KIND_R8)
                  end do
               end do
            end block
         end if
       case (3)
         nzf = nz
         if (trim(cc_wrap%field_config%import_fields(i)%vertical_axis) == 'interface') nzf = nz + 1
         field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
            ungriddedLBound=(/1/), ungriddedUBound=(/nzf/), &
            name=trim(cc_wrap%field_config%import_fields(i)%standard_name), rc=rc)
         call check(rc, "FieldCreate 3d")
         call ESMF_FieldGet(field, farrayPtr=fptr3, rc=rc)
         call check(rc, "FieldGet 3d")
         block
            integer :: ix, iy, iz
            do iz = 1, nzf
               do iy = 1, ny
                  do ix = 1, nx
                     select case (trim(cc_wrap%field_config%import_fields(i)%catchem_var))
                      case ('T')
                        fptr3(ix, iy, iz) = 280.0_ESMF_KIND_R8
                      case ('QV')
                        fptr3(ix, iy, iz) = 0.01_ESMF_KIND_R8
                      case ('RH')
                        fptr3(ix, iy, iz) = 0.5_ESMF_KIND_R8
                      case ('PEDGE')
                        fptr3(ix, iy, iz) = 100000.0_ESMF_KIND_R8 - 10000.0_ESMF_KIND_R8 * real(iz - 1, ESMF_KIND_R8)
                      case ('PMID')
                        fptr3(ix, iy, iz) = 95000.0_ESMF_KIND_R8 - 10000.0_ESMF_KIND_R8 * real(iz - 1, ESMF_KIND_R8)
                      case ('CLDF')
                        fptr3(ix, iy, iz) = 0.2_ESMF_KIND_R8
                      case default
                        fptr3(ix, iy, iz) = real(ix + 10 * iy + 100 * iz, ESMF_KIND_R8)
                     end select
                  end do
               end do
            end do
         end block
       case (4)
         field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
            ungriddedLBound=(/1, 1/), ungriddedUBound=(/nz, ntr/), &
            name=trim(cc_wrap%field_config%import_fields(i)%standard_name), rc=rc)
         call check(rc, "FieldCreate 4d")
         call ESMF_FieldGet(field, farrayPtr=fptr4, rc=rc)
         call check(rc, "FieldGet 4d")
         call ESMF_FieldGet(field, array=array, rc=rc)
         call check(rc, "FieldGet 4d array")
         call ESMF_InfoGetFromHost(array, tracer_info, rc=rc)
         call check(rc, "Tracer metadata lookup")
         call ESMF_InfoSet(tracer_info, 'tracerNames', exchange_tracer_names, rc=rc)
         call check(rc, "Tracer names metadata")
         call ESMF_InfoSet(tracer_info, 'tracerUnits', exchange_tracer_units, rc=rc)
         call check(rc, "Tracer units metadata")
         block
            integer :: ix, iy, iz, itr
            do itr = 1, ntr
               do iz = 1, nz
                  do iy = 1, ny
                     do ix = 1, nx
                        select case (itr)
                         case (1)
                           fptr4(ix, iy, iz, itr) = 1.0E-2_ESMF_KIND_R8
                         case (2)
                           fptr4(ix, iy, iz, itr) = 1.0E-6_ESMF_KIND_R8
                         case (3)
                           fptr4(ix, iy, iz, itr) = 2.0E-6_ESMF_KIND_R8
                         case (4)
                           fptr4(ix, iy, iz, itr) = 7.0E-6_ESMF_KIND_R8
                         case (5)
                           fptr4(ix, iy, iz, itr) = 8.0E-6_ESMF_KIND_R8
                        end select
                     end do
                  end do
               end do
            end do
         end block
       case default
         cycle
      end select

      call ESMF_StateAdd(importState, (/field/), rc=rc)
      call check(rc, "StateAdd")
   end do

   ! 4. Drive the production transform — the exact ursa run-phase path
   call ESMF_TimeSet(currTime, yy=2021, mm=3, dd=22, h=6, rc=rc)
   call check(rc, "TimeSet")
   call transform_nuopc_to_catchem(cc_wrap, importState, currTime, rc)
   if (rc /= ESMF_SUCCESS) then
      print *, 'FAIL: transform_nuopc_to_catchem rc=', rc, &
         ' (see test_nuopc_transform.log)'
      error stop 1
   end if
   print *, 'PASS: transform_nuopc_to_catchem over the full required mapping'

   ! UFS supplies PF*SAN on nlev levels.  Preserve upstream CATChem's
   ! nlev+1 interface reconstruction by repeating the final host value.
   do i = 1, cc_wrap%field_config%n_import_fields
      if (trim(cc_wrap%field_config%import_fields(i)%vertical_axis) /= 'level_to_interface') cycle
      if (.not. allocated(cc_wrap%met_buf_3d(i)%data) .or. &
         size(cc_wrap%met_buf_3d(i)%data, 3) /= nz + 1) then
         error stop 'PF*SAN level_to_interface buffer extent is invalid'
      end if
      if (.not. exact_equal(cc_wrap%met_buf_3d(i)%data(1,1,1), 111.0_c_double) .or. &
         .not. exact_equal(cc_wrap%met_buf_3d(i)%data(1,1,nz+1), cc_wrap%met_buf_3d(i)%data(1,1,nz))) then
         error stop 'PF*SAN level_to_interface buffer was not reconstructed as upstream'
      end if
   end do

   ! Confirm gas and aerosol slots reach their native CATChem units while
   ! the host-owned rank-4 boundary remains otherwise untouched.
   if (abs(cc_wrap%chem_buf_4d(1,1,1,gas_idx) - 1.0E-6_c_double * &
      cc_wrap%tracer_map%host_to_catchem(2)) > 1.0E-12_c_double) error stop 'Gas import conversion failed'
   if (abs(cc_wrap%chem_buf_4d(1,1,1,aerosol_idx) - 2.0E-6_c_double * &
      cc_wrap%tracer_map%host_to_catchem(3)) > 1.0E-12_c_double) error stop 'Aerosol import conversion failed'
   if (abs(cc_wrap%tracer_map%host_to_catchem(2) * cc_wrap%tracer_map%catchem_to_host(2) - 1.0_c_double) > 1.0E-12_c_double .or. &
      abs(cc_wrap%tracer_map%host_to_catchem(3) * cc_wrap%tracer_map%catchem_to_host(3) - 1.0_c_double) > 1.0E-12_c_double) &
      error stop 'Tracer conversion factors are not reciprocal'
   block
      logical :: sphum_preserved
      sphum_preserved = .false.
      do i = 1, cc_wrap%field_config%n_import_fields
         if (trim(cc_wrap%field_config%import_fields(i)%host_tracer_var) /= 'QV') cycle
         if (allocated(cc_wrap%met_buf_3d(i)%data)) then
            sphum_preserved = abs(cc_wrap%met_buf_3d(i)%data(1,1,1) - 1.0E-2_c_double) <= 1.0E-12_c_double
         end if
      end do
      if (.not. sphum_preserved) error stop 'Host-owned sphum was not preserved'
   end block
   print *, 'PASS: representative gas/aerosol conversion and sphum ownership'

   ! 5. Z0 must have landed in the C++ met state
   raw_z0_ptr = catchem_state_get_pointer_2d(cc_wrap%catchem_model%state_mgr_ptr, "Z0" // c_null_char)
   if (.not. c_associated(raw_z0_ptr)) then
      print *, 'FAIL: Z0 not bound in C++ met state after transform'
      error stop 1
   end if
   call c_f_pointer(raw_z0_ptr, z0_ptr, [nx, ny])
   if (abs(z0_ptr(1, 1) - 1.5_fp) > 1.0e-12_fp) then
      print *, 'FAIL: Z0 cm->m conversion not applied:', z0_ptr(1, 1)
      error stop 1
   end if
   print *, 'PASS: Z0 converted (150 cm -> 1.5 m) and stored'

   ! 6. PM2.5/PM10 diagnostics update — the 2026-08-11 17:07 run-phase
   !    abort regression. Must succeed and register the fields in the
   !    C++ DiagnosticManager (where NUOPC export reads from).
   call update_pm_diagnostics(cc_wrap, rc)
   if (rc /= CC_SUCCESS) then
      print *, 'FAIL: update_pm_diagnostics rc=', rc
      error stop 1
   end if
   block
      character(len=64), allocatable :: diag_names(:)
      call cc_wrap%catchem_model%get_diagnostic_names(diag_names, rc=rc)
      if (rc /= 0 .or. .not. allocated(diag_names)) then
         print *, 'FAIL: get_diagnostic_names rc=', rc
         error stop 1
      end if
      if (.not. any(diag_names == 'pm25') .or. .not. any(diag_names == 'pm10')) then
         print *, 'FAIL: pm25/pm10 not registered in the C++ diagnostic manager'
         error stop 1
      end if
   end block
   block
      real(fp), allocatable :: pm_diag(:,:,:)

      call cc_wrap%catchem_model%get_diagnostic('pm25', pm_diag, rc)
      if (rc /= 0 .or. abs(pm_diag(1,1,1)) > 1.0e-3_fp) &
         error stop 'PM25 diagnostic has an incorrect C++ buffer precision'
      deallocate(pm_diag)
      call cc_wrap%catchem_model%get_diagnostic('pm10', pm_diag, rc)
      if (rc /= 0 .or. abs(pm_diag(1,1,1)) > 1.0e-3_fp) &
         error stop 'PM10 diagnostic has an incorrect C++ buffer precision'
      deallocate(pm_diag)
   end block
   print *, 'PASS: PM2.5/PM10 diagnostics updated and registered'

   ! 7. Test transform_catchem_to_nuopc (full C++ Core -> C API -> Fortran export state)
   exportState = ESMF_StateCreate(name="export", rc=rc)
   call check(rc, "StateCreate export")

   do i = 1, cc_wrap%field_config%n_export_fields
      select case (cc_wrap%field_config%export_fields(i)%dimensions)
       case (2)
         field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
            name=trim(cc_wrap%field_config%export_fields(i)%standard_name), rc=rc)
         call check(rc, "Export FieldCreate 2d")
       case (3)
         field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
            ungriddedLBound=(/1/), ungriddedUBound=(/nz/), &
            name=trim(cc_wrap%field_config%export_fields(i)%standard_name), rc=rc)
         call check(rc, "Export FieldCreate 3d")
       case (4)
         field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
            ungriddedLBound=(/1, 1/), ungriddedUBound=(/nz, ntr/), &
            name=trim(cc_wrap%field_config%export_fields(i)%standard_name), rc=rc)
         call check(rc, "Export FieldCreate 4d")
       case default
         cycle
      end select
      call ESMF_StateAdd(exportState, (/field/), rc=rc)
      call check(rc, "Export StateAdd")
   end do

   call transform_catchem_to_nuopc(cc_wrap, exportState, rc)
   if (rc /= ESMF_SUCCESS) then
      print *, 'FAIL: transform_catchem_to_nuopc rc=', rc
      error stop 1
   end if
   do i = 1, cc_wrap%field_config%n_export_fields
      if (cc_wrap%field_config%export_fields(i)%dimensions /= 4) cycle
      call ESMF_StateGet(exportState, trim(cc_wrap%field_config%export_fields(i)%standard_name), field, rc=rc)
      call check(rc, 'Export StateGet tracer field')
      call ESMF_FieldGet(field, farrayPtr=fptr4, rc=rc)
      call check(rc, 'Export FieldGet tracer pointer')
      if (abs(fptr4(1,1,1,1) - 1.0E-2_ESMF_KIND_R8) > 1.0E-12_ESMF_KIND_R8) &
         error stop 'Host-owned sphum changed during export'
      if (abs(fptr4(1,1,1,2) - 1.0E-6_ESMF_KIND_R8) > 1.0E-12_ESMF_KIND_R8 .or. &
         abs(fptr4(1,1,1,3) - 2.0E-6_ESMF_KIND_R8) > 1.0E-12_ESMF_KIND_R8) &
         error stop 'Chemical export conversion failed'
      block
         real(fp), allocatable :: pm_diag(:,:,:)
         call cc_wrap%catchem_model%get_diagnostic('pm25', pm_diag, rc)
         if (rc /= 0 .or. .not. allocated(pm_diag)) then
            error stop 'PM25 diagnostic was not exported in ug m-3'
         end if
         if (abs(fptr4(1,1,1,4) - pm_diag(1,1,1)) > 1.0e-12_ESMF_KIND_R8) &
            error stop 'PM25 diagnostic was not exported in ug m-3'
         deallocate(pm_diag)
         call cc_wrap%catchem_model%get_diagnostic('pm10', pm_diag, rc)
         if (rc /= 0 .or. .not. allocated(pm_diag)) then
            error stop 'PM10 diagnostic was not exported in ug m-3'
         end if
         if (abs(fptr4(1,1,1,5) - pm_diag(1,1,1)) > 1.0e-12_ESMF_KIND_R8) &
            error stop 'PM10 diagnostic was not exported in ug m-3'
         deallocate(pm_diag)
      end block
   end do
   print *, 'PASS: transform_catchem_to_nuopc over export state'

   call ESMF_Finalize(rc=rc)
   print *, 'All NUOPC contract harness tests passed!'

contains

   subroutine parity_from_c(c_value, value)
      character(kind=c_char), intent(in) :: c_value(*)
      character(len=*), intent(out) :: value
      integer :: j
      value = ''
      do j = 1, len(value)
         if (c_value(j) == c_null_char) exit
         value(j:j) = c_value(j)
      end do
   end subroutine parity_from_c

   subroutine check(rc_in, what)
      integer, intent(in) :: rc_in
      character(len=*), intent(in) :: what
      if (rc_in /= ESMF_SUCCESS) then
         print *, 'FAIL (setup): ', what, ' rc=', rc_in
         error stop 1
      end if
   end subroutine check

end program nuopc_contract_harness
