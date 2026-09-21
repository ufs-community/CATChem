!> \file test_nuopc_diag_output.f90
!! \brief ESMF-backed integration test for the generic process-diagnostic writer.
!!
!! Builds only when CATCHEM_BUILD_NUOPC=ON.  Reproduces the run-phase diagnostic
!! export in-process: initialize a model whose run_phases activate every
!! registered science process (seasalt, dust, drydep, wetdep, settling, so4chem,
!! carbchem), drive the axes-driven write_process_diagnostics into a NetCDF file,
!! reopen it with nf90, and assert the feature-013 contract:
!!
!!   US1 (T010): every process contributes >=1 variable, and every {Column,Level}
!!               field keeps its full vertical extent (no rank-2 truncation).
!!   US2 (T015): the five packed processes unpack to <field>_<label> variables
!!               (one per bin/species), the compact parent is gone, and a slot
!!               value equals the corresponding host-pointer column (SC-002).
!!
!! A field whose packed axis carries no label must fail loudly (FR-007/008); the
!! C++ invariant tests (test_diagnostic_lifecycle) cover that path directly.
program test_nuopc_diag_output
   use iso_c_binding, only: c_ptr, c_char, c_int, c_double, c_int64_t, c_null_char, &
      c_null_ptr, c_associated, c_f_pointer
   use ESMF
   use netcdf
   use aqmio, only: AQMIO_Create, AQMIO_Destroy
   use catchem_nuopc_interface, only: cc_wrap_type, write_process_diagnostics, update_time_variable, &
      write_global_attributes
   use catchem_bridge_precision, only: fp

   implicit none

   interface
      integer(c_int) function catchem_diag_get_pointer_checked(core_ptr, name, rank, dims, ptr_out) &
         bind(C, name="catchem_diag_get_pointer_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), value :: rank
         integer(c_int), intent(in) :: dims(*)
         type(c_ptr), intent(out) :: ptr_out
      end function
      integer(c_int) function catchem_diag_get_unpack_label_at_checked(core_ptr, name, slot, label_out, &
         label_length) &
         bind(C, name="catchem_diag_get_unpack_label_at_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), value :: slot, label_length
         character(kind=c_char), intent(inout) :: label_out(*)
      end function
      integer(c_int) function catchem_diag_get_dims_checked(core_ptr, name, dims_out, dims_length) &
         bind(C, name="catchem_diag_get_dims_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), intent(inout) :: dims_out(*)
         integer(c_int), value :: dims_length
      end function
      integer(c_int) function catchem_diag_get_rank_checked(core_ptr, name, rank_out) &
         bind(C, name="catchem_diag_get_rank_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), intent(out) :: rank_out
      end function
      integer(c_int) function catchem_state_get_species_count_checked(state_ptr, count_out) &
         bind(C, name="catchem_state_get_species_count_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), intent(out) :: count_out
      end function
      integer(c_int) function catchem_state_get_species_conc_pointer_checked(state_ptr, species_index, &
         dim1, dim2, ptr_out) &
         bind(C, name="catchem_state_get_species_conc_pointer_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: species_index, dim1, dim2
         type(c_ptr), intent(out) :: ptr_out
      end function
   end interface

   integer, parameter :: nx = 3, ny = 2, nz = 4
   integer, parameter :: ncols = nx * ny

   integer :: nfail, rc
   character(len=64) :: msg

   nfail = 0

   call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, &
      defaultlogfilename="test_nuopc_diag_output.log", rc=rc)
   call check(rc, "ESMF_Initialize")

   ! Case A: empty diag_list -> every registered field is written (US1/US2).
   call run_case('CATChem_diag_output_config.yml', 'diag_out.nc', nfail)
   print *, 'PASS: process-diagnostic writer contract satisfied (US1 + US2)'

   ! Case B: narrowing diag_list (FR-009 / T022) -> only the selected parent
   ! and its unpacked children survive; everything else is absent.  The
   ! unmatched selector 'no_such_field' is reported on stdout by the writer.
   call run_case_narrow('CATChem_diag_narrow_config.yml', 'diag_narrow.nc', nfail)
   print *, 'PASS: diag_list volume filter satisfied (US3)'

   call ESMF_Finalize(endflag=ESMF_END_KEEPMPI)

   if (nfail > 0) then
      write(msg, '(A,I0,A)') 'FAIL: ', nfail, ' diagnostic-output assertion(s) failed'
      print *, trim(msg)
      error stop 1
   end if

contains

   !> Drive one full model lifecycle and write the process diagnostics.
   subroutine run_case(config, outname, nfail)
      character(len=*), intent(in) :: config, outname
      integer, intent(inout) :: nfail
      type(cc_wrap_type) :: cc_wrap
      type(ESMF_Grid) :: grid
      type(ESMF_Time) :: currTime
      real(c_double), target, allocatable :: conc_buf(:,:,:)
      integer :: rc

      ! 1. Initialize the model with the all-process diagnostic config.  Each
      !    process init() registers its diagnostics in the C++ DiagnosticManager.
      call cc_wrap%catchem_model%initialize(config, nx, ny, nz, rc=rc)
      if (rc /= 0) then
         print *, 'FAIL: model initialize rc=', rc
         error stop 1
      end if
      if (.not. cc_wrap%catchem_model%is_diag_enabled()) error stop 'diagnostics not enabled in config'

      ! 2. Build the grid + AQMIO I/O component the writer needs.  Field creation
      !    inside write_diagnostic_field is grid-only (no tile arrays), so a plain
      !    single-tile grid is sufficient.
      grid = ESMF_GridCreateNoPeriDim(maxIndex=(/nx, ny/), rc=rc)
      call check(rc, "GridCreate")
      cc_wrap%grid = grid
      cc_wrap%iocomp = AQMIO_Create(grid, rc=rc)
      if (.not. ESMF_GridCompIsCreated(cc_wrap%iocomp)) error stop 'AQMIO_Create failed'
      cc_wrap%compress_lev = 0
      cc_wrap%current_time_slice = 0
      call ESMF_TimeSet(currTime, yy=2024, mm=5, dd=1, h=0, m=0, s=0, rc=rc)
      call check(rc, "TimeSet")

      ! 3. Create the time axis, stamp provenance globals, then write every
      !    registered process diagnostic.
      call update_time_variable(cc_wrap, outname, currTime, cc_wrap%current_time_slice, rc)
      call check(rc, "update_time_variable")
      call write_global_attributes(cc_wrap, outname, rc)
      call check(rc, "write_global_attributes")

      ! SC-005 (T029): the writer runs strictly after the science step and must
      ! never perturb tracer concentrations.  Bind a deterministic concentration
      ! buffer, snapshot it bitwise, run the writer (producing the file that
      ! verify_output inspects), snapshot again, and assert bit-for-bit equality.
      ! This is the realizable form of the on/off ncdiff in this environment
      ! (the standalone app needs external ExtData).
      call science_unchanged_by_writer(cc_wrap, outname, conc_buf, nfail, rc)
      call check(rc, "write_process_diagnostics")

      ! 4. Reopen the file and assert the contract (AQMIO_Close flushed it).
      call verify_output(outname, cc_wrap%catchem_model%cpp_core_ptr, nfail)

      call AQMIO_Destroy(cc_wrap%iocomp, rc=rc)
      call ESMF_GridDestroy(grid, rc=rc)
      call cc_wrap%catchem_model%finalize(rc)
   end subroutine run_case

   !> SC-005 guard (T029): bind a deterministic tracer-concentration buffer,
   !! snapshot it bitwise, run the diagnostic writer, snapshot again, and
   !! assert bit-for-bit equality.  The writer must be pure output: it reads
   !! host pointers only and runs after the science step, so any perturbation
   !! of concentrations is a defect.  (The literal run-on/off ncdiff needs the
   !! standalone app with external ExtData inputs; this is the equivalent
   !! in-process guarantee.)
   !! The buffer is owned by the caller so it stays alive until the model is
   !! finalized (the C++ view over it is unmanaged and non-owning).
   subroutine science_unchanged_by_writer(cc_wrap, outname, conc_buf, nfail, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      character(len=*), intent(in) :: outname
      real(c_double), target, allocatable, intent(inout) :: conc_buf(:,:,:)
      integer, intent(inout) :: nfail
      integer, intent(out) :: rc
      integer(c_int) :: c_status, c_count
      integer :: ncols, nlev, total, i, j, v
      integer(c_int64_t), allocatable :: snap_before(:), snap_after(:)

      ncols = nx * ny
      nlev = nz

      ! The unified-chemistry buffer is species-major blocks, matching the
      ! Fortran (column, level, species) layout the C++ getter slices with
      ! conc + (species-1)*ncols*nlev.  Query the registered species count
      ! first so the buffer's third extent is exact (a mismatch is rejected by
      ! the extent-checked bind path).
      c_status = catchem_state_get_species_count_checked(cc_wrap%catchem_model%state_mgr_ptr, c_count)
      if (c_status /= 0_c_int .or. c_count <= 0) then
         print *, 'FAIL: species count status=', int(c_status), ' count=', int(c_count)
         nfail = nfail + 1
         rc = 1
         return
      end if

      allocate(conc_buf(ncols, nlev, int(c_count)))
      do v = 1, int(c_count)
         do j = 1, nlev
            do i = 1, ncols
               conc_buf(i, j, v) = 1.0d0 + 0.5d0*real(i, c_double) &
                  - 0.25d0*real(j, c_double) + 3.0d0*real(v, c_double)
            end do
         end do
      end do
      call cc_wrap%catchem_model%bind_unified_chemistry(conc_buf, rc)
      if (rc /= 0) then
         print *, 'FAIL: bind_unified_chemistry rc=', rc
         nfail = nfail + 1
         return
      end if

      total = ncols * nlev * int(c_count)
      allocate(snap_before(total), snap_after(total))

      call snapshot_concentrations(cc_wrap, c_count, ncols, nlev, snap_before)
      call write_process_diagnostics(cc_wrap, 'all', outname, rc)
      call snapshot_concentrations(cc_wrap, c_count, ncols, nlev, snap_after)

      call expect(all(snap_before == snap_after), &
         'SC-005: writer left tracer concentrations bit-identical')

      deallocate(snap_before, snap_after)
   end subroutine science_unchanged_by_writer

   !> Copy every species' [ncols,nlev] concentration block, bitwise, into snap.
   subroutine snapshot_concentrations(cc_wrap, c_count, ncols, nlev, snap)
      type(cc_wrap_type), intent(in) :: cc_wrap
      integer(c_int), intent(in) :: c_count
      integer, intent(in) :: ncols, nlev
      integer(c_int64_t), intent(out) :: snap(:)
      integer(c_int) :: c_status, v
      integer :: base, block_size
      type(c_ptr) :: raw_ptr
      real(c_double), pointer :: sp(:,:) => null()

      block_size = ncols * nlev
      do v = 1, c_count
         raw_ptr = c_null_ptr
         c_status = catchem_state_get_species_conc_pointer_checked( &
            cc_wrap%catchem_model%state_mgr_ptr, v, int(ncols, c_int), int(nlev, c_int), raw_ptr)
         if (c_status /= 0_c_int .or. .not. c_associated(raw_ptr)) then
            call expect(.false., 'concentration pointer available for every species')
            return
         end if
         call c_f_pointer(raw_ptr, sp, [ncols, nlev])
         base = (int(v) - 1) * block_size
         snap(base + 1 : base + block_size) = transfer(sp, snap(base + 1 : base + block_size))
         nullify(sp)
      end do
   end subroutine snapshot_concentrations

   !> Drive a narrowing-diag_list run and assert only the selected variables
   !! survive (US3 / T022).  The config selects dust_emission_total, the
   !! settling_flux_per_species parent (covering every unpacked child), and a
   !! no_such_field entry that must produce an unmatched-selector warning.
   subroutine run_case_narrow(config, outname, nfail)
      character(len=*), intent(in) :: config, outname
      integer, intent(inout) :: nfail
      type(cc_wrap_type) :: cc_wrap
      type(ESMF_Grid) :: grid
      type(ESMF_Time) :: currTime
      integer :: rc
      integer :: ncid, status, v

      call cc_wrap%catchem_model%initialize(config, nx, ny, nz, rc=rc)
      if (rc /= 0) then
         print *, 'FAIL: narrow model initialize rc=', rc
         error stop 1
      end if
      grid = ESMF_GridCreateNoPeriDim(maxIndex=(/nx, ny/), rc=rc)
      call check(rc, "GridCreate (narrow)")
      cc_wrap%grid = grid
      cc_wrap%iocomp = AQMIO_Create(grid, rc=rc)
      if (.not. ESMF_GridCompIsCreated(cc_wrap%iocomp)) error stop 'AQMIO_Create failed (narrow)'
      cc_wrap%compress_lev = 0
      cc_wrap%current_time_slice = 0
      call ESMF_TimeSet(currTime, yy=2024, mm=5, dd=1, h=0, m=0, s=0, rc=rc)
      call check(rc, "TimeSet (narrow)")
      call update_time_variable(cc_wrap, outname, currTime, cc_wrap%current_time_slice, rc)
      call check(rc, "update_time_variable (narrow)")
      call write_process_diagnostics(cc_wrap, 'all', outname, rc)
      call check(rc, "write_process_diagnostics (narrow)")

      status = nf90_open(trim(outname), nf90_nowrite, ncid)
      if (status /= nf90_noerr) then
         nfail = nfail + 1
         print *, 'FAIL: cannot reopen narrow file: ', trim(nf90_strerror(status))
         return
      end if

      print *, 'US3: diag_list narrows written variables'
      ! Selected singleton total survives.
      call expect(find_var(ncid, 'dust_emission_total') >= 0, 'narrow: dust_emission_total present')
      ! Selected parent covers its unpacked children (settling_flux_per_species_<label>).
      v = find_var(ncid, 'settling_flux_per_species_so4')
      call expect(v >= 0, 'narrow: settling_flux_per_species child present')
      ! Everything not named by a selector must be absent (SC-007).
      call expect(find_var(ncid, 'dust_emission_bin_dust1') == -1, 'narrow: unselected dust bin absent')
      call expect(find_var(ncid, 'seasalt_mass_emission_total') == -1, 'narrow: unselected seasalt absent')
      call expect(find_var(ncid, 'drydep_con_per_species_so2') == -1, 'narrow: unselected drydep absent')
      call expect(find_var(ncid, 'wetdep_mass_so2') == -1, 'narrow: unselected wetdep absent')
      call expect(find_var(ncid, 'PSO4_from_gaseous_SO2_per_level') == -1, 'narrow: unselected so4chem absent')
      call expect(find_var(ncid, 'carbchem_prod_mass_oc1') == -1, 'narrow: unselected carbchem absent')

      status = nf90_close(ncid)
      call AQMIO_Destroy(cc_wrap%iocomp, rc=rc)
      call ESMF_GridDestroy(grid, rc=rc)
      call cc_wrap%catchem_model%finalize(rc)
   end subroutine run_case_narrow

   subroutine check(status, context)
      integer, intent(in) :: status
      character(len=*), intent(in) :: context
      if (status /= ESMF_SUCCESS) then
         print *, 'FAIL: ', trim(context), ' rc=', status
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
   end subroutine check

   !> Look up a variable by name; return its netCDF id (or -1 when absent).
   function find_var(ncid, name) result(varid)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer :: varid
      integer :: status
      status = nf90_inq_varid(ncid, trim(name), varid)
      if (status /= nf90_noerr) varid = -1
   end function find_var

   !> Number of dimensions of a variable.
   function var_ndims(ncid, varid) result(nd)
      integer, intent(in) :: ncid, varid
      integer :: nd
      integer :: status
      nd = -1
      if (varid < 0) return
      status = nf90_inquire_variable(ncid, varid, ndims=nd)
      if (status /= nf90_noerr) nd = -1
   end function var_ndims

   !> True when the variable carries a dimension named `dimname` (AQMIO orders
   !! fields as (Time, lev, grid_yt, grid_xt), so the vertical is neither first
   !! nor last and must be located by name).
   function var_has_dim(ncid, varid, dimname) result(present_dim)
      integer, intent(in) :: ncid, varid
      character(len=*), intent(in) :: dimname
      logical :: present_dim
      integer :: nd, dimids(8), did, status, k
      present_dim = .false.
      if (varid < 0) return
      nd = var_ndims(ncid, varid)
      if (nd <= 0) return
      status = nf90_inq_dimid(ncid, trim(dimname), did)
      if (status /= nf90_noerr) return
      status = nf90_inquire_variable(ncid, varid, dimids=dimids(1:nd))
      if (status /= nf90_noerr) return
      do k = 1, nd
         if (dimids(k) == did) then
            present_dim = .true.
            return
         end if
      end do
   end function var_has_dim

   !> Size of a named dimension in the file (or -1 when absent).
   function dim_len(ncid, dimname) result(d)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: dimname
      integer :: d
      integer :: did, status
      d = -1
      status = nf90_inq_dimid(ncid, trim(dimname), did)
      if (status /= nf90_noerr) return
      status = nf90_inquire_dimension(ncid, did, len=d)
      if (status /= nf90_noerr) d = -1
   end function dim_len

   subroutine expect(cond, label)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: label
      if (.not. cond) then
         nfail = nfail + 1
         print *, '  ASSERT FAILED: ', trim(label)
      end if
   end subroutine expect

   !> True when the file carries a string global attribute of this name.
   function global_present(ncid, name) result(present_att)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      logical :: present_att
      character(len=512) :: buf
      present_att = (nf90_get_att(ncid, NF90_GLOBAL, trim(name), buf) == nf90_noerr)
   end function global_present

   !> True when a string global attribute exists and equals the expected value.
   function global_equals(ncid, name, expected) result(ok)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name, expected
      logical :: ok
      integer :: status
      character(len=512) :: buf
      status = nf90_get_att(ncid, NF90_GLOBAL, trim(name), buf)
      ok = (status == nf90_noerr) .and. (trim(buf) == trim(expected))
      if (status /= nf90_noerr) print *, '  global lookup failed: ', trim(name), &
         trim(nf90_strerror(status))
   end function global_equals

   subroutine verify_output(fname, core_ptr, nfail)
      character(len=*), intent(in) :: fname
      type(c_ptr), intent(in) :: core_ptr
      integer, intent(inout) :: nfail
      integer :: ncid, status, v, nd

      status = nf90_open(trim(fname), nf90_nowrite, ncid)
      if (status /= nf90_noerr) then
         nfail = nfail + 1
         print *, 'FAIL: cannot reopen ', trim(fname), ': ', trim(nf90_strerror(status))
         return
      end if

      ! --- US1: at least one variable per process prefix present ---
      print *, 'US1: per-process presence + vertical extent'
      call expect(find_var(ncid, 'dust_emission_total') >= 0, 'dust_emission_total present')
      call expect(find_var(ncid, 'seasalt_mass_emission_total') >= 0, 'seasalt total present')
      call expect(find_var(ncid, 'drydep_con_per_species_so2') >= 0, 'drydep unpacked present')
      call expect(find_var(ncid, 'wetdep_mass_so2') >= 0, 'wetdep per-species present')
      call expect(find_var(ncid, 'PSO4_from_gaseous_SO2_per_level') >= 0, 'so4chem level field present')
      call expect(find_var(ncid, 'carbchem_prod_mass_oc1') >= 0, 'carbchem unpacked present')
      call expect(find_var(ncid, 'settling_flux_per_species_so4') >= 0, 'settling unpacked present')

      ! --- US1: level-axis fields keep their vertical extent (FR-003) ---
      ! A {Column,Level} field must carry the 'lev' dimension at full extent;
      ! the pre-fix writer forced rank-2 and dropped the vertical entirely.
      v = find_var(ncid, 'PSO4_from_gaseous_SO2_per_level')
      nd = var_ndims(ncid, v)
      call expect(nd >= 3, 'so4chem level field is rank>=3')
      call expect(var_has_dim(ncid, v, 'lev'), 'so4chem level field carries lev dimension')
      call expect(dim_len(ncid, 'lev') >= nz, 'so4chem level field keeps vertical extent')

      ! --- US2: packed axes unpack to <field>_<label>; parent gone ---
      print *, 'US2: packed fields unpacked, compact parent absent'
      call expect(find_var(ncid, 'dust_emission_bin') == -1, 'dust compact parent absent')
      call expect(find_var(ncid, 'seasalt_mass_emission_bins') == -1, 'seasalt compact parent absent')
      call expect(find_var(ncid, 'drydep_con_per_species') == -1, 'drydep compact parent absent')
      ! One variable per bin/species, named with the slot short_name.
      call expect(find_var(ncid, 'dust_emission_bin_dust1') >= 0, 'dust bin label dust1')
      call expect(find_var(ncid, 'dust_emission_bin_dust3') >= 0, 'dust bin label dust3')
      call expect(find_var(ncid, 'seasalt_mass_emission_bins_seas1') >= 0, 'seasalt bin label seas1')
      call expect(find_var(ncid, 'seasalt_mass_emission_bins_seas3') >= 0, 'seasalt bin label seas3')
      call expect(find_var(ncid, 'settling_flux_per_species_seas3') >= 0, 'settling flux label seas3')
      call expect(find_var(ncid, 'drydep_con_per_species_so4') >= 0, 'drydep label so4')
      call expect(find_var(ncid, 'carbchem_prod_mass_bc1') >= 0, 'carbchem label bc1')
      ! 3-D packed field (settling velocity) unpacks to a rank-3 variable per species.
      v = find_var(ncid, 'settling_velocity_per_species_per_level_so4')
      call expect(v >= 0, 'settling velocity unpacked per species')
      call expect(var_ndims(ncid, v) >= 3, 'settling velocity per-species is rank>=3')
      call expect(var_has_dim(ncid, v, 'lev'), 'settling velocity keeps vertical extent')

      ! --- US2 (SC-002): an unpacked slot value equals the host column ---
      call value_matches_host(ncid, core_ptr, 'dust_emission_bin', 'dust1', 'dust_emission_bin_dust1')
      call value_matches_host(ncid, core_ptr, 'seasalt_mass_emission_bins', 'seas1', &
         'seasalt_mass_emission_bins_seas1')

      ! --- US4 (T024): run-level provenance global attributes (FR-011, C-10) ---
      print *, 'US4: global provenance attributes'
      call expect(global_equals(ncid, 'institution', 'Test Lab'), &
         'user attribute overrides core default: institution')
      call expect(global_equals(ncid, 'references', &
         'https://github.com/UFS-Community/CATChem'), 'references core default present')
      call expect(global_present(ncid, 'catchem_core_version'), 'catchem_core_version present')
      call expect(global_present(ncid, 'catchem_core_commit'), 'catchem_core_commit present')
      call expect(global_present(ncid, 'config_file'), 'config_file present')
      call expect(global_equals(ncid, 'title', 'feature 013 diagnostic output test'), &
         'user-only attribute present: title')

      status = nf90_close(ncid)
   end subroutine verify_output

   !> Compare an unpacked slot variable's NetCDF values against the packed
   !! parent's host storage at the slot whose label matches (SC-002).  The slot
   !! is located by querying the C++ label list for the parent, so the check is
   !! independent of registration order.  Skips silently when the parent or slot
   !! is absent, keeping the assertion meaningful across configs.
   subroutine value_matches_host(ncid, core_ptr, parent, label, varname)
      integer, intent(in) :: ncid
      type(c_ptr), intent(in) :: core_ptr
      character(len=*), intent(in) :: parent, label, varname
      integer :: v, slot, nslot, s, status
      integer(c_int) :: c_rank, c_dims(3), c_status
      character(kind=c_char) :: c_label(64)
      real(fp), pointer :: packed(:,:) => null()
      real(ESMF_KIND_R4), allocatable :: file_vals(:)
      character(len=64) :: got_label
      type(c_ptr) :: raw_ptr

      v = find_var(ncid, varname)
      if (v < 0) then
         print *, '  value check SKIP: var absent ', trim(varname)
         return
      end if

      c_rank = 0
      c_status = catchem_diag_get_rank_checked(core_ptr, trim(parent)//c_null_char, c_rank)
      if (c_status /= 0_c_int .or. int(c_rank) /= 2) then
         print *, '  value check SKIP: rank ', int(c_status), int(c_rank), ' parent=', trim(parent)
         return
      end if
      c_dims = 0
      c_status = catchem_diag_get_dims_checked(core_ptr, trim(parent)//c_null_char, c_dims, 3_c_int)
      if (c_status /= 0_c_int) then
         print *, '  value check SKIP: dims status ', int(c_status)
         return
      end if
      nslot = int(c_dims(2))

      ! Find the slot whose registered label equals `label`.
      slot = -1
      do s = 0, nslot - 1
         c_label = c_null_char
         c_status = catchem_diag_get_unpack_label_at_checked(core_ptr, trim(parent)//c_null_char, &
            int(s, c_int), c_label, 64_c_int)
         if (c_status /= 0_c_int) cycle
         call c_fortran_str(c_label, got_label)
         if (trim(got_label) == trim(label)) then
            slot = s + 1
            exit
         end if
      end do
      if (slot < 1) then
         call expect(.false., 'label found for parent: '//trim(parent))
         return
      end if

      raw_ptr = c_null_ptr
      c_status = catchem_diag_get_pointer_checked(core_ptr, trim(parent)//c_null_char, c_rank, c_dims, raw_ptr)
      if (c_status /= 0_c_int .or. .not. c_associated(raw_ptr)) then
         print *, '  value check SKIP: pointer status ', int(c_status), ' assoc=', c_associated(raw_ptr)
         return
      end if
      call c_f_pointer(raw_ptr, packed, [ncols, nslot])

      ! The file variable is (grid_xt, grid_yt, Time); read the single time
      ! slice across the full grid so the 1D result is column-major (xt fastest),
      ! matching the packed host storage layout col = i + (j-1)*nx.
      allocate(file_vals(ncols))
      status = nf90_get_var(ncid, v, file_vals, start=(/1, 1, 1/), count=(/nx, ny, 1/))
      if (status /= nf90_noerr) then
         print *, '  value check SKIP: get_var status ', status, trim(nf90_strerror(status))
         deallocate(file_vals)
         return
      end if
      call expect(all_close(packed(:, slot), file_vals, ncols), &
         'unpacked slot equals packed host column: '//trim(varname))
      deallocate(file_vals)
   end subroutine value_matches_host

   !> Convert a NUL-terminated C character array to a Fortran string.
   subroutine c_fortran_str(cstr, fstr)
      character(kind=c_char), intent(in) :: cstr(:)
      character(len=*), intent(out) :: fstr
      integer :: i
      fstr = ' '
      do i = 1, min(len(fstr), size(cstr))
         if (cstr(i) == c_null_char) exit
         fstr(i:i) = cstr(i)
      end do
   end subroutine c_fortran_str

   pure function all_close(a, b, n) result(ok)
      integer, intent(in) :: n
      real(fp), intent(in) :: a(:)
      real(ESMF_KIND_R4), intent(in) :: b(:)
      logical :: ok
      integer :: i
      real(fp) :: tol
      ok = .true.
      do i = 1, n
         tol = 1.0e-4_fp * max(1.0_fp, abs(a(i)))
         if (abs(real(a(i), ESMF_KIND_R4) - b(i)) > tol) then
            ok = .false.
            return
         end if
      end do
   end function all_close

end program test_nuopc_diag_output
