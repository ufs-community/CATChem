module DustScienceBridge_Mod
   use iso_c_binding
   use catchem_bridge_precision, only: fp
   use catchem_bridge_constants, only: g0
   use DustScheme_FENGSHA_Mod, only: compute_fengsha
   use DustScheme_GINOUX_Mod, only: compute_ginoux
   use DustCommon_Mod, only: DustSchemeFENGSHAConfig, DustSchemeGINOUXConfig

   implicit none
   private

contains

   subroutine run_dust_science_bridge( &
      n_cols, n_levels, n_species, n_total_species, n_soil, dt, &
      active_scheme, diagnostics, &
      fengsha_alpha, fengsha_gamma, fengsha_drylimit_factor, fengsha_moisture_factor, fengsha_kvhmax, &
      fengsha_drag_option, fengsha_horizflux_option, fengsha_moist_option, fengsha_distribution_option, &
      ginoux_ch_du, n_ginoux_ch_du, &
      airden, delp, clayfrac, frlake, frsno, gvf, lai, lwi, rdrag, sandfrac, &
      soilm, gwettop, ssm, tskin, u10m, v10m, ustar, ustar_threshold, z0, &
      species_density, species_radius, species_lower_radius, species_upper_radius, &
      bin_species_names, species_names, &
      conc, tendency, &
      diag_emission_total, diag_emission_bin, diag_horizontal_flux, diag_moisture_correction, diag_effective_threshold, diag_utar_threshold, &
      diagnostic_species_id, n_diag_species &
      ) bind(C, name="run_dust_science_bridge")

      ! C-interoperable metadata
      ! n_species     : number of dust bin species (subset the scheme operates on)
      ! n_total_species: number of species in the full unified conc array
      integer(c_int), value :: n_cols, n_levels, n_species, n_total_species, n_soil
      real(c_double), value :: dt
      character(kind=c_char), intent(in) :: active_scheme(*)
      integer(c_int), value :: diagnostics
      ! Dust bin names (subset) and the full chemistry catalog names.  The
      ! bridge resolves each bin to its slot in the full conc array by name,
      ! so no species index crosses the C boundary (mirrors settling).
      character(kind=c_char), intent(in) :: bin_species_names(32, n_species)
      character(kind=c_char), intent(in) :: species_names(32, n_total_species)

      ! Scheme tuning options staged by DustProcess::init from the runtime
      ! YAML.  The C++ layer owns parsing; the bridge only applies them onto
      ! the scheme configuration types so Fengsha/Ginoux no longer run on
      ! compiled defaults.  Ch_DU carries one multiplier per dust size bin.
      real(c_double), value :: fengsha_alpha, fengsha_gamma, fengsha_drylimit_factor
      real(c_double), value :: fengsha_moisture_factor, fengsha_kvhmax
      integer(c_int), value :: fengsha_drag_option, fengsha_horizflux_option
      integer(c_int), value :: fengsha_moist_option, fengsha_distribution_option
      integer(c_int), value :: n_ginoux_ch_du
      real(c_double), intent(in) :: ginoux_ch_du(n_ginoux_ch_du)

      ! 2D/3D Pointers from C++ Kokkos state
      type(c_ptr), value :: airden
      type(c_ptr), value :: delp
      type(c_ptr), value :: clayfrac
      type(c_ptr), value :: frlake
      type(c_ptr), value :: frsno
      type(c_ptr), value :: gvf
      type(c_ptr), value :: lai
      type(c_ptr), value :: lwi
      type(c_ptr), value :: rdrag
      type(c_ptr), value :: sandfrac
      type(c_ptr), value :: soilm
      type(c_ptr), value :: gwettop
      type(c_ptr), value :: ssm
      type(c_ptr), value :: tskin
      type(c_ptr), value :: u10m
      type(c_ptr), value :: v10m
      type(c_ptr), value :: ustar
      type(c_ptr), value :: ustar_threshold
      type(c_ptr), value :: z0

      ! Species properties
      type(c_ptr), value :: species_density
      type(c_ptr), value :: species_radius
      type(c_ptr), value :: species_lower_radius
      type(c_ptr), value :: species_upper_radius

      ! Multi-dimensional Views
      type(c_ptr), value :: conc
      type(c_ptr), value :: tendency

      ! Diagnostics
      type(c_ptr), value :: diag_emission_total
      type(c_ptr), value :: diag_emission_bin
      type(c_ptr), value :: diag_horizontal_flux
      type(c_ptr), value :: diag_moisture_correction
      type(c_ptr), value :: diag_effective_threshold
      type(c_ptr), value :: diag_utar_threshold
      ! Per-process diagnostics.  The diag_* pointers reference
      ! DiagnosticManager field storage; they are null (and n_diag_species is
      ! 0) when diagnostics are disabled, so the per-bin pointers are only
      ! c_f_pointer'd when n_diag_species > 0.  diagnostic_species_id holds
      ! 1-based LOCAL bin positions within the canonical dust-bin subset the
      ! scheme iterates (species_idx space); the size-1 dummy is never
      ! dereferenced when n_diag_species == 0 (mirrors the settling bridge).
      integer(c_int), value :: n_diag_species
      integer(c_int), intent(in) :: diagnostic_species_id(max(n_diag_species,1))

      ! Local Fortran Pointers
      real(c_double), pointer :: f_airden(:,:), f_delp(:,:)
      real(c_double), pointer :: f_clayfrac(:), f_frlake(:), f_frsno(:), f_gvf(:), f_lai(:)
      integer(c_int), pointer :: f_lwi(:)
      real(c_double), pointer :: f_rdrag(:), f_sandfrac(:), f_soilm(:,:), f_gwettop(:), f_ssm(:), f_tskin(:)
      real(c_double), pointer :: f_u10m(:), f_v10m(:), f_ustar(:), f_ustar_threshold(:), f_z0(:)
      real(c_double), pointer :: f_conc(:,:,:), f_tendency(:,:,:)

      real(c_double), pointer :: f_density(:), f_radius(:), f_lower_radius(:), f_upper_radius(:)

      real(c_double), pointer :: f_diag_emission_total(:), f_diag_emission_bin(:,:)
      real(c_double), pointer :: f_diag_horizontal_flux(:), f_diag_moisture_correction(:)
      real(c_double), pointer :: f_diag_effective_threshold(:), f_diag_utar_threshold(:,:)

      ! Local Slices for computation
      real(fp) :: col_airden(n_levels)
      real(fp) :: col_soilm(n_soil)
      real(fp) :: col_conc(n_levels, n_species)
      real(fp) :: col_tendency(n_levels, n_species)
      ! Surface-emission schemes receive a one-layer physical surface column.
      ! The bridge, not the scheme, owns CATChem's bottom-to-top indexing.
      real(fp) :: surface_airden(1)
      real(fp) :: surface_conc(1, n_species)
      real(fp) :: surface_tendency(1, n_species)

      real(fp) :: f_species_density(n_species)
      real(fp) :: f_species_radius(n_species)
      real(fp) :: f_species_lower_radius(n_species)
      real(fp) :: f_species_upper_radius(n_species)

      real(fp) :: col_emission_total
      ! Per-bin diagnostic scratch is sized by the diag_species subset count,
      ! not the bin count: the scheme fills slot diag_idx (1..n_diag_species).
      real(fp) :: col_emission_bin(max(n_diag_species,1))
      real(fp) :: col_horizontal_flux
      real(fp) :: col_moisture_correction
      real(fp) :: col_effective_threshold
      real(fp) :: col_utar_threshold(max(n_diag_species,1))

      type(DustSchemeFENGSHAConfig) :: fengsha_config
      type(DustSchemeGINOUXConfig)  :: ginoux_config
      character(len=32) :: local_scheme
      integer :: icol, i, ispec, k
      ! Map from each dust bin (1..n_species) to its column in the full conc
      ! array (1..n_total_species), resolved by trimmed name comparison.
      integer :: target_species(n_species)
      character(len=32) :: bin_names(n_species)

      ! Map Scheme Name
      local_scheme = ""
      do i = 1, 32
         if (active_scheme(i) == c_null_char) exit
         local_scheme(i:i) = active_scheme(i)
      end do
      local_scheme = trim(local_scheme)

      ! Apply the YAML tuning options staged by the C++ process layer onto
      ! the scheme configuration types.  Ch_DU must match the dust size-bin
      ! count declared by the scheme type.
      fengsha_config%alpha = real(fengsha_alpha, fp)
      fengsha_config%gamma = real(fengsha_gamma, fp)
      fengsha_config%drylimit_factor = real(fengsha_drylimit_factor, fp)
      fengsha_config%moist_correction_factor = real(fengsha_moisture_factor, fp)
      fengsha_config%kvhmax = real(fengsha_kvhmax, fp)
      fengsha_config%drag_option = int(fengsha_drag_option)
      fengsha_config%horizflux_option = int(fengsha_horizflux_option)
      fengsha_config%moist_option = int(fengsha_moist_option)
      fengsha_config%distribution_option = int(fengsha_distribution_option)
      if (n_ginoux_ch_du /= size(ginoux_config%Ch_DU)) then
         write(*,'(A,I0,A,I0)') 'FATAL ERROR: DustScienceBridge ginoux Ch_DU length ', n_ginoux_ch_du, &
            ' does not match the scheme bin count ', size(ginoux_config%Ch_DU)
         call flush(6)
         error stop "FATAL ERROR: DustScienceBridge ginoux Ch_DU length mismatch"
      end if
      ginoux_config%Ch_DU = real(ginoux_ch_du, fp)

      ! Pointer Associations
      call c_f_pointer(airden, f_airden, [n_cols, n_levels])
      call c_f_pointer(delp, f_delp, [n_cols, n_levels])
      call c_f_pointer(clayfrac, f_clayfrac, [n_cols])
      call c_f_pointer(frlake, f_frlake, [n_cols])
      call c_f_pointer(frsno, f_frsno, [n_cols])
      call c_f_pointer(gvf, f_gvf, [n_cols])
      call c_f_pointer(lai, f_lai, [n_cols])
      call c_f_pointer(lwi, f_lwi, [n_cols])
      call c_f_pointer(rdrag, f_rdrag, [n_cols])
      call c_f_pointer(sandfrac, f_sandfrac, [n_cols])
      call c_f_pointer(soilm, f_soilm, [n_cols, n_soil])
      call c_f_pointer(gwettop, f_gwettop, [n_cols])
      call c_f_pointer(ssm, f_ssm, [n_cols])
      call c_f_pointer(tskin, f_tskin, [n_cols])
      call c_f_pointer(u10m, f_u10m, [n_cols])
      call c_f_pointer(v10m, f_v10m, [n_cols])
      call c_f_pointer(ustar, f_ustar, [n_cols])
      call c_f_pointer(ustar_threshold, f_ustar_threshold, [n_cols])
      call c_f_pointer(z0, f_z0, [n_cols])

      call c_f_pointer(conc, f_conc, [n_cols, n_levels, n_total_species])
      call c_f_pointer(tendency, f_tendency, [n_cols, n_levels, n_total_species])

      ! Resolve each dust bin to its slot in the full chemistry array by name.
      do i = 1, n_species
         bin_names(i) = c_name_to_fortran(bin_species_names(:, i))
         target_species(i) = 0
         do k = 1, n_total_species
            if (trim(c_name_to_fortran(species_names(:, k))) == trim(bin_names(i))) then
               target_species(i) = k
               exit
            end if
         end do
         if (target_species(i) == 0) then
            write(*,'(A,A)') 'FATAL ERROR: DustScienceBridge cannot resolve dust bin ', trim(bin_names(i))
            call flush(6)
            error stop "FATAL ERROR: DustScienceBridge unresolved dust bin species"
         end if
      end do

      call c_f_pointer(species_density, f_density, [n_species])
      call c_f_pointer(species_radius, f_radius, [n_species])
      call c_f_pointer(species_lower_radius, f_lower_radius, [n_species])
      call c_f_pointer(species_upper_radius, f_upper_radius, [n_species])

      if (diagnostics /= 0 .and. n_diag_species > 0) then
         call c_f_pointer(diag_emission_total, f_diag_emission_total, [n_cols])
         call c_f_pointer(diag_emission_bin, f_diag_emission_bin, [n_cols, n_diag_species])
         call c_f_pointer(diag_horizontal_flux, f_diag_horizontal_flux, [n_cols])
         call c_f_pointer(diag_moisture_correction, f_diag_moisture_correction, [n_cols])
         call c_f_pointer(diag_effective_threshold, f_diag_effective_threshold, [n_cols])
         call c_f_pointer(diag_utar_threshold, f_diag_utar_threshold, [n_cols, n_diag_species])
      end if

      ! Cast species properties to fp precision
      f_species_density = real(f_density, fp)
      f_species_radius  = real(f_radius, fp)
      f_species_lower_radius = real(f_lower_radius, fp)
      f_species_upper_radius = real(f_upper_radius, fp)

      ! Process Columns
      do icol = 1, n_cols

         col_airden(:) = real(f_airden(icol, :), fp)
         col_soilm(:)  = real(f_soilm(icol, :), fp)
         ! Gather the dust-bin concentrations out of the full unified array
         ! using the name-resolved map (mirrors settling's conc_2d gather).
         do ispec = 1, n_species
            do k = 1, n_levels
               col_conc(k, ispec) = real(f_conc(icol, k, target_species(ispec)), fp)
            end do
         end do
         col_tendency(:,:) = 0.0_fp

         col_emission_total = 0.0_fp
         col_emission_bin(:) = 0.0_fp
         col_horizontal_flux = 0.0_fp
         col_moisture_correction = 0.0_fp
         col_effective_threshold = 0.0_fp
         col_utar_threshold(:) = 0.0_fp

         if (local_scheme == "fengsha") then
            surface_airden(1) = col_airden(1)
            surface_conc(1, :) = col_conc(1, :)
            surface_tendency = 0.0_fp
            call compute_fengsha( &
               1, n_species, fengsha_config, g0, &
               surface_airden, real(f_clayfrac(icol), fp), real(f_frlake(icol), fp), real(f_frsno(icol), fp), &
               real(f_gvf(icol), fp), real(f_lai(icol), fp), int(f_lwi(icol)), real(f_rdrag(icol), fp), &
               real(f_sandfrac(icol), fp), col_soilm, real(f_ssm(icol), fp), real(f_tskin(icol), fp), &
               real(f_ustar(icol), fp), real(f_ustar_threshold(icol), fp), real(f_z0(icol), fp), &
               f_species_radius, f_species_lower_radius, f_species_upper_radius, &
               surface_conc, surface_tendency, &
               dust_emission_total=col_emission_total, &
               dust_emission_per_bin=col_emission_bin, &
               dust_horizontal_flux=col_horizontal_flux, &
               dust_moisture_correction=col_moisture_correction, &
               dust_effective_threshold=col_effective_threshold, &
               diagnostic_species_id=diagnostic_species_id)
            col_tendency(1, :) = surface_tendency(1, :)
         else if (local_scheme == "ginoux") then
            call compute_ginoux( &
               n_levels, n_species, ginoux_config, g0, &
               col_airden, real(f_frlake(icol), fp), real(f_frsno(icol), fp), real(f_gwettop(icol), fp), &
               int(f_lwi(icol)), real(f_ssm(icol), fp), real(f_tskin(icol), fp), &
               real(f_u10m(icol), fp), real(f_v10m(icol), fp), &
               f_species_density, f_species_radius, &
               col_conc, col_tendency, &
               dust_emission_total=col_emission_total, &
               dust_emission_per_bin=col_emission_bin, &
               utar_threshold_per_bin=col_utar_threshold, &
               diagnostic_species_id=diagnostic_species_id)
         end if

         ! Write tendencies (Convert kg/m2/s at surface layer 1 to ug/kg concentration change)
         ! dqa = flux * dt * g0 / delp(1) * 1.0e9.  Scatter back into the full
         ! unified conc array via the name-resolved map.
         do ispec = 1, n_species
            if (abs(col_tendency(1, ispec)) > 1.0e-32_fp) then
               col_tendency(1, ispec) = (col_tendency(1, ispec) * real(dt, fp) * g0 / real(f_delp(icol, 1), fp)) * 1.0e9_fp

               f_tendency(icol, 1, target_species(ispec)) = &
                  f_tendency(icol, 1, target_species(ispec)) + real(col_tendency(1, ispec) / real(dt, fp), c_double)
               f_conc(icol, 1, target_species(ispec)) = &
                  max(0.0_c_double, f_conc(icol, 1, target_species(ispec)) + real(col_tendency(1, ispec), c_double))
            end if
         end do

         ! Write diagnostics
         if (diagnostics /= 0 .and. n_diag_species > 0) then
            f_diag_emission_total(icol) = real(col_emission_total, c_double)
            f_diag_emission_bin(icol, :) = real(col_emission_bin(1:n_diag_species), c_double)
            f_diag_horizontal_flux(icol) = real(col_horizontal_flux, c_double)
            f_diag_moisture_correction(icol) = real(col_moisture_correction, c_double)
            f_diag_effective_threshold(icol) = real(col_effective_threshold, c_double)
            f_diag_utar_threshold(icol, :) = real(col_utar_threshold(1:n_diag_species), c_double)
         end if

      end do

   contains
      function c_name_to_fortran(c_name) result(name)
         character(kind=c_char), intent(in) :: c_name(32)
         character(len=32) :: name
         integer :: ic
         name = ''
         do ic = 1, 32
            if (c_name(ic) == c_null_char) exit
            name(ic:ic) = c_name(ic)
         end do
         name = trim(adjustl(name))
      end function c_name_to_fortran
   end subroutine run_dust_science_bridge

end module DustScienceBridge_Mod
