! C ABI adapter for the legacy GOCART2G settling science path.
!
! This bridge reproduces the upstream/develop ProcessSettlingInterface
! execution exactly: one `compute_gocart` call per column covering all
! settling species, using the metadata (non-Mie) branch of the scheme.
! Units and layout follow specs/011-restore-numerical-parity/contracts/
! settling-science-bridge.md:
!   - every array is column-major with the flattened column fastest;
!   - the vertical order is bottom-to-top (compute_gocart reverses internally);
!   - concentrations are µg/kg on both sides of the boundary (kg/kg conversion
!     happens inside the scheme);
!   - radius is passed in µm (µm -> m conversion happens inside the scheme).
module SettlingScienceBridge_Mod
   use iso_c_binding, only: c_int, c_double, c_char, c_ptr, c_f_pointer
   use catchem_bridge_precision, only: fp
   use catchem_bridge_error, only: CC_SUCCESS
   use GOCART2G_MieMod, only: GOCART2G_Mie
   use SettlingCommon_Mod, only: SettlingSchemeGOCARTConfig
   use SettlingScheme_GOCART_Mod, only: compute_gocart
   implicit none

   ! Aerosol optics (Mie) tables backing the legacy simple_scheme path.
   !
   ! GOCART2G_Mie is a Fortran derived type with deferred-length allocatable and
   ! pointer components, so it cannot cross the C ABI; the store therefore lives
   ! on this side and C++ only ever passes fixed-width names and paths.  The
   ! tables are loaded once during process initialization and are read-only
   ! afterwards, which mirrors the upstream ChemStateType%MieData ownership and
   ! keeps every MPI rank's copy identical (same config, same files).
   type(GOCART2G_Mie), allocatable, target, save :: mie_store(:)
   character(len=32), allocatable, target, save :: mie_names_store(:)
   ! Zero-sized target used as the scheme's Mie actual argument on the metadata
   ! path, preserving today's "no tables" behavior.
   type(GOCART2G_Mie), target, save :: empty_mie(0)
contains
   subroutine run_settling_mie_init(n_files, type_names, file_paths, init_rc) &
      bind(C, name='run_settling_mie_init')
      integer(c_int), value :: n_files
      character(kind=c_char), intent(in) :: type_names(32,n_files)
      character(kind=c_char), intent(in) :: file_paths(512,n_files)
      integer(c_int), intent(out) :: init_rc

      integer :: idx, local_rc

      ! Idempotent re-initialization: drop any previous store so repeated
      ! process init (tests, restarts) never leaves stale tables behind.
      if (allocated(mie_store)) deallocate(mie_store)
      if (allocated(mie_names_store)) deallocate(mie_names_store)
      init_rc = 0_c_int
      if (n_files <= 0) return

      allocate(mie_store(n_files), mie_names_store(n_files))
      do idx = 1, n_files
         mie_names_store(idx) = c_name_to_fortran(type_names(:,idx))
         mie_store(idx) = GOCART2G_Mie(c_path_to_fortran(file_paths(:,idx)), rc=local_rc)
         if (local_rc /= CC_SUCCESS) then
            ! Report the 1-based index of the failing table; the GOCART reader
            ! already wrote the offending path to stderr.
            init_rc = int(idx, c_int)
            return
         end if
      end do
   contains
      function c_path_to_fortran(c_text) result(text)
         character(kind=c_char), intent(in) :: c_text(512)
         character(len=512) :: text
         integer :: i
         text = ''
         do i = 1, 512
            text(i:i) = c_text(i)
         end do
         text = trim(adjustl(text))
      end function c_path_to_fortran
      function c_name_to_fortran(c_text) result(text)
         character(kind=c_char), intent(in) :: c_text(32)
         character(len=32) :: text
         integer :: i
         text = ''
         do i = 1, 32
            text(i:i) = c_text(i)
         end do
         text = trim(adjustl(text))
      end function c_name_to_fortran
   end subroutine run_settling_mie_init

   subroutine run_settling_science_bridge(n_columns, n_levels, n_aerosols, n_total_species, &
      dt, scale_factor, swelling_rh_max, correction_maring, maring_dust_only, &
      airden, delp, pmid, rh, temperature, z_edge, &
      aerosol_species_names, species_names, species_is_dust, species_is_hydrophilic, radius, density, &
      concentration, simple_scheme, aerosol_mie_names, &
      diag_velocity, diag_flux, diagnostic_species_id, n_diag_species, bridge_rc) &
      bind(C, name='run_settling_science_bridge')
      integer(c_int), value :: n_columns, n_levels, n_aerosols, n_total_species
      integer(c_int), value :: correction_maring, maring_dust_only
      real(c_double), value :: dt, scale_factor, swelling_rh_max
      real(c_double), intent(in) :: airden(n_columns,n_levels), delp(n_columns,n_levels)
      real(c_double), intent(in) :: pmid(n_columns,n_levels), rh(n_columns,n_levels)
      real(c_double), intent(in) :: temperature(n_columns,n_levels)
      real(c_double), intent(in) :: z_edge(n_columns,n_levels+1)
      character(kind=c_char), intent(in) :: aerosol_species_names(32,n_aerosols)
      character(kind=c_char), intent(in) :: species_names(32,n_total_species)
      integer(c_int), intent(in) :: species_is_dust(n_aerosols)
      integer(c_int), intent(in) :: species_is_hydrophilic(n_aerosols)
      real(c_double), intent(in) :: radius(n_aerosols), density(n_aerosols)
      real(c_double), intent(inout) :: concentration(n_columns,n_levels,n_total_species)
      integer(c_int), value :: simple_scheme
      character(kind=c_char), intent(in) :: aerosol_mie_names(32,n_aerosols)
      ! Per-process diagnostics.  diag_velocity/diag_flux are C pointers into the
      ! DiagnosticManager field storage; they are null (and n_diag_species is 0)
      ! when diagnostics are disabled, so they are only c_f_pointer'd when
      ! nonzero to avoid a zero-extent assumed-size dummy.  diagnostic_species_id
      ! holds 1-based LOCAL positions within the aerosol subset (the species_idx
      ! space compute_gocart loops over); the size-1 dummy is never dereferenced
      ! when n_diag_species == 0.
      type(c_ptr), value :: diag_velocity, diag_flux
      integer(c_int), value :: n_diag_species
      integer(c_int), intent(in) :: diagnostic_species_id(max(n_diag_species,1))
      integer(c_int), intent(out) :: bridge_rc

      type(SettlingSchemeGOCARTConfig) :: params
      type(GOCART2G_Mie), pointer :: mie_actual(:)
      character(len=32) :: aerosol_names(n_aerosols)
      character(len=32) :: aerosol_mies(n_aerosols)
      integer :: target_species(n_aerosols)
      integer :: species_mie_map(n_aerosols)
      logical :: is_dust(n_aerosols)
      logical :: is_hydrophilic(n_aerosols)
      real(fp) :: species_radius(n_aerosols), species_density(n_aerosols)
      real(fp) :: airden_1d(n_levels), delp_1d(n_levels), pmid_1d(n_levels)
      real(fp) :: rh_1d(n_levels), t_1d(n_levels), z_1d(n_levels+1)
      real(fp) :: conc_2d(n_levels,n_aerosols), tend_2d(n_levels,n_aerosols)
      real(c_double), pointer :: f_diag_velocity(:,:,:), f_diag_flux(:,:)
      real(fp), allocatable :: col_velocity(:,:), col_flux(:)
      integer :: column, species, k

      bridge_rc = 0_c_int
      if (n_aerosols <= 0) return

      allocate(col_velocity(n_levels, n_diag_species))
      allocate(col_flux(n_diag_species))
      f_diag_velocity => null(); f_diag_flux => null()
      if (n_diag_species > 0) then
         call c_f_pointer(diag_velocity, f_diag_velocity, [n_columns, n_levels, n_diag_species])
         call c_f_pointer(diag_flux, f_diag_flux, [n_columns, n_diag_species])
      end if

      ! Resolve settling species against the full chemistry list by name
      ! (no index crossing the boundary); mirror upstream trimmed comparison.
      do species = 1, n_aerosols
         aerosol_names(species) = c_name_to_fortran(aerosol_species_names(:,species))
         target_species(species) = 0
         do k = 1, n_total_species
            if (trim(c_name_to_fortran(species_names(:,k))) == trim(aerosol_names(species))) then
               target_species(species) = k
               exit
            end if
         end do
         if (target_species(species) == 0) then
            bridge_rc = 1_c_int
            return
         end if
      end do

      ! Scheme parameters.  scale_factor is retained for configuration
      ! compatibility but, exactly like upstream, compute_gocart does not consume
      ! it on either path.  swelling_rh_max applies to the metadata path only
      ! (the scheme gates the clamp on .not. simple_scheme).
      params%scheme_name = 'gocart'
      params%scale_factor = real(scale_factor, fp)
      params%simple_scheme = (simple_scheme /= 0)
      params%swelling_rh_max = real(swelling_rh_max, fp)
      params%correction_maring = (correction_maring /= 0)
      params%maring_dust_only = (maring_dust_only /= 0)

      is_dust = (species_is_dust /= 0)
      is_hydrophilic = (species_is_hydrophilic /= 0)
      do species = 1, n_aerosols
         ! Radii stay in µm: the scheme performs the µm -> m conversion.
         species_radius(species) = real(radius(species), fp)
         species_density(species) = real(density(species), fp)
      end do

      if (params%simple_scheme) then
         ! Optics-table path: resolve each settling species' __mie_name against
         ! the loaded table names, mirroring upstream chemstate_init_mie_data's
         ! SpcMieMap construction (trimmed comparison, 1-based, 0 = unresolved).
         if (.not. allocated(mie_store)) then
            ! C++ validates this at init; defensive backstop so a mis-wired
            ! caller can never reach Chem_SettlingSimple without tables.
            bridge_rc = 2_c_int
            return
         end if
         mie_actual => mie_store
         do species = 1, n_aerosols
            aerosol_mies(species) = c_name_to_fortran(aerosol_mie_names(:,species))
            species_mie_map(species) = 0
            if (len_trim(aerosol_mies(species)) == 0) cycle
            do k = 1, size(mie_names_store)
               if (trim(mie_names_store(k)) == trim(aerosol_mies(species))) then
                  species_mie_map(species) = k
                  exit
               end if
            end do
         end do
         ! Any settling species without a table aborts before advancing columns
         ! (the scheme's own "Invalid Mie data mapping" guard is the last resort).
         if (any(species_mie_map == 0)) then
            bridge_rc = 2_c_int
            return
         end if
      else
         ! Metadata path: no Mie tables (mirrors upstream simple_scheme=false).
         mie_actual => empty_mie
         species_mie_map = 0
      end if

      do column = 1, n_columns
         do k = 1, n_levels
            airden_1d(k) = real(airden(column,k), fp)
            delp_1d(k) = real(delp(column,k), fp)
            pmid_1d(k) = real(pmid(column,k), fp)
            rh_1d(k) = real(rh(column,k), fp)
            t_1d(k) = real(temperature(column,k), fp)
         end do
         do k = 1, n_levels + 1
            z_1d(k) = real(z_edge(column,k), fp)
         end do
         do species = 1, n_aerosols
            do k = 1, n_levels
               ! µg/kg both sides; kg/kg conversion happens inside the scheme.
               conc_2d(k,species) = real(concentration(column,k,target_species(species)), fp)
               tend_2d(k,species) = 0.0_fp
            end do
         end do
         col_velocity = 0.0_fp
         col_flux = 0.0_fp

         ! One call per column for all settling species, exactly like the
         ! upstream run_gocart_scheme_column.  Scheme-internal failures report
         ! through CC_Error (stderr banner) and return without aborting, which
         ! is the legacy behavior; replacement tendencies are written back for
         ! whatever the scheme produced.  The optional diagnostic arguments are
         ! only passed when diagnostics are enabled (n_diag_species > 0); the
         ! no-diagnostic call path must stay bit-identical to the legacy one.
         if (n_diag_species > 0) then
            call compute_gocart(n_levels, n_aerosols, params, &
               airden_1d, delp_1d, pmid_1d, rh_1d, t_1d, real(dt, fp), z_1d, &
               aerosol_names, mie_actual, species_mie_map, species_radius, species_density, &
               is_dust, is_hydrophilic, conc_2d, tend_2d, &
               settling_velocity_per_species_per_level=col_velocity, &
               settling_flux_per_species=col_flux, &
               diagnostic_species_id=diagnostic_species_id)
         else
            call compute_gocart(n_levels, n_aerosols, params, &
               airden_1d, delp_1d, pmid_1d, rh_1d, t_1d, real(dt, fp), z_1d, &
               aerosol_names, mie_actual, species_mie_map, species_radius, species_density, &
               is_dust, is_hydrophilic, conc_2d, tend_2d)
         end if

         do species = 1, n_aerosols
            do k = 1, n_levels
               ! Replacement tendencies (max(0, qa) applied inside the scheme).
               ! Clamp at zero to flush denormalized FP negatives that survive
               ! the fp -> c_double cast (matches the emission bridge convention).
               concentration(column,k,target_species(species)) = &
                  max(0.0_c_double, real(tend_2d(k,species), c_double))
            end do
         end do

         if (n_diag_species > 0) then
            ! The scheme already writes velocity with the level flip
            ! (SD(1,1,num_layers:1:-1)), so surface=1 matches every other
            ! process's manager layout; no extra flipping here.
            f_diag_velocity(column,:,:) = real(col_velocity, c_double)
            f_diag_flux(column,:)        = real(col_flux, c_double)
         end if
      end do

      deallocate(col_velocity, col_flux)
   contains
      function c_name_to_fortran(c_name) result(name)
         character(kind=c_char), intent(in) :: c_name(32)
         character(len=32) :: name
         integer :: i
         name = ''
         do i = 1, size(c_name)
            name(i:i) = c_name(i)
         end do
         name = trim(adjustl(name))
      end function c_name_to_fortran
   end subroutine run_settling_science_bridge
end module SettlingScienceBridge_Mod
