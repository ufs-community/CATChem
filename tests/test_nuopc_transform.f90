!> \file test_nuopc_transform.f90
!! \brief ESMF-backed test of the NUOPC import transform.
!!
!! Builds only when CATCHEM_BUILD_NUOPC=ON. Reproduces the ursa
!! run-phase mechanism in-process: load the real field mapping
!! (CATChem_field_mapping.yml), create an ESMF field for every
!! required import entry, and drive transform_nuopc_to_catchem —
!! asserting success and the Z0 cm->m conversion landing in the met
!! state (the 2026-08-11 first-timestep abort regression).
program test_nuopc_transform
   use iso_c_binding, only: c_ptr, c_char, c_double, c_int, c_associated, c_f_pointer, c_null_char
   use ESMF
   use catchem_nuopc_interface, only: load_field_config, transform_nuopc_to_catchem, &
      transform_catchem_to_nuopc, field_config, cc_wrap_type, update_pm_diagnostics, TRACER_CHEMICAL
   use catchem_bridge_error, only: CC_SUCCESS
   use catchem_bridge_precision, only: fp
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
   end interface

   integer, parameter :: nx = 4, ny = 2, nz = 5, ntr = 2
   type(cc_wrap_type) :: cc_wrap
   type(CATChem_Model) :: parity_model
   type(ESMF_Grid) :: grid
   type(ESMF_State) :: importState, exportState
   type(ESMF_Field) :: field
   type(ESMF_Time) :: currTime
   real(ESMF_KIND_R8), pointer :: fptr2(:, :), fptr3(:, :, :), fptr4(:, :, :, :)
   real(c_double), pointer :: z0_ptr(:,:) => null()
   type(c_ptr) :: raw_z0_ptr
   character(len=256) :: errmsg
   integer :: rc, i, nzf, parity_issues
   integer(c_int) :: parity_count, parity_status
   character(kind=c_char) :: parity_c_name(64)
   character(len=64) :: parity_name
   character(len=512) :: parity_report
   character(len=64), parameter :: parity_expected(3) = [character(len=64) :: &
      'unfamiliar_alpha', 'unfamiliar_beta', 'unfamiliar_gamma']

   call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, &
      defaultlogfilename="test_nuopc_transform.log", rc=rc)
   call check(rc, "ESMF_Initialize")

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
         block
            integer :: ix, iy, iz, itr
            do itr = 1, ntr
               do iz = 1, nz
                  do iy = 1, ny
                     do ix = 1, nx
                        fptr4(ix, iy, iz, itr) = real(ix + 10 * iy + 100 * iz + 1000 * itr, ESMF_KIND_R8)
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

   ! 3b. Production builds cc_wrap%tracer_map inside catchem_nuopc_init from
   !     tracerinfo metadata; this harness drives the transform directly, so
   !     construct a minimal valid map mirroring init: two host chemical
   !     tracers mapped 1:1 onto CATChem species with unit conversion factors.
   !     Without it the rank-4 chemistry import has no validated mapping and
   !     transform_nuopc_to_catchem fails by design.
   allocate(cc_wrap%tracer_map%names(ntr))
   allocate(cc_wrap%tracer_map%units(ntr))
   allocate(cc_wrap%tracer_map%nuopc_to_cc(ntr))
   allocate(cc_wrap%tracer_map%entry_kind(ntr))
   allocate(cc_wrap%tracer_map%host_to_catchem(ntr))
   allocate(cc_wrap%tracer_map%catchem_to_host(ntr))
   cc_wrap%tracer_map%names = [character(len=128) :: 'tracer_a', 'tracer_b']
   cc_wrap%tracer_map%units = 'kg kg-1'
   cc_wrap%tracer_map%nuopc_to_cc = [1, 2]
   cc_wrap%tracer_map%entry_kind = TRACER_CHEMICAL
   cc_wrap%tracer_map%host_to_catchem = 1.0_c_double
   cc_wrap%tracer_map%catchem_to_host = 1.0_c_double

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
   print *, 'PASS: transform_catchem_to_nuopc over export state'

   call ESMF_Finalize(rc=rc)
   print *, 'All NUOPC transform tests passed!'

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

end program test_nuopc_transform
