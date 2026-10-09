module SeaSaltScienceBridge_Mod
   use iso_c_binding, only: c_ptr, c_f_pointer, c_double, c_char, c_associated, c_null_char, c_bool, c_int
   use catchem_bridge_precision, only: fp
   use catchem_bridge_constants, only: g0, AIRMW, PI
   use SeaSaltCommon_Mod, only: SeaSaltSchemeGONG97Config, SeaSaltSchemeGONG03Config, SeaSaltSchemeGEOS12Config
   use SeaSaltScheme_GONG97_Mod, only: compute_gong97
   use SeaSaltScheme_GONG03_Mod, only: compute_gong03
   use SeaSaltScheme_GEOS12_Mod, only: compute_geos12
   implicit none
contains

   subroutine run_seasalt_science_bridge( &
      n_cols, n_levels, n_species, n_total_species, dt, &
      active_scheme, diagnostics, &
      gong97_scale_factor, gong97_weibull_flag, &
      gong03_scale_factor, gong03_weibull_flag, &
      geos12_scale_factor, geos12_weibull_flag, &
   ! Met Pointers
      c_frocean, c_frseaice, c_lat, c_lon, c_sst, c_u10m, c_v10m, c_ustar, c_delp, &
   ! Species Metadata
      species_density, species_radius, species_lower_radius, species_upper_radius, is_gas_arr, species_mw_g, &
   ! Species name maps (bin subset + full chemistry catalog)
      bin_species_names, species_names, &
   ! Concentrations and Tendency
      c_conc, c_tendency, &
   ! Diagnostics
      c_diag_mass_total, c_diag_num_total, c_diag_mass_bin, c_diag_num_bin, &
      diagnostic_species_id, n_diag_species &
      ) bind(C, name="run_seasalt_science_bridge")

      ! n_species      : number of sea-salt bin species (subset the scheme runs on)
      ! n_total_species: number of species in the full unified conc array
      integer(c_int), value :: n_cols, n_levels, n_species, n_total_species
      real(c_double), value :: dt
      character(kind=c_char), intent(in) :: active_scheme(*)
      integer(c_int), value :: diagnostics
      ! Sea-salt bin names (subset) and the full chemistry catalog names.  The
      ! bridge resolves each bin to its slot in the full conc array by name,
      ! so no species index crosses the C boundary (mirrors settling).
      character(kind=c_char), intent(in) :: bin_species_names(32, n_species)
      character(kind=c_char), intent(in) :: species_names(32, n_total_species)

      ! Scheme tuning options staged by SeaSaltProcess::init from the runtime
      ! YAML.  All three schemes are carried so the bridge dispatches on the
      ! active scheme without needing the C++ layer to know the defaults.
      real(c_double), value :: gong97_scale_factor, gong03_scale_factor, geos12_scale_factor
      integer(c_int), value :: gong97_weibull_flag, gong03_weibull_flag, geos12_weibull_flag

      type(c_ptr), value :: c_frocean, c_frseaice, c_lat, c_lon, c_sst, c_u10m, c_v10m, c_ustar, c_delp
      type(c_ptr), value :: c_conc, c_tendency
      type(c_ptr), value :: c_diag_mass_total, c_diag_num_total, c_diag_mass_bin, c_diag_num_bin

      ! Per-process diagnostics.  The diag_* pointers reference
      ! DiagnosticManager field storage; they are null (and n_diag_species is
      ! 0) when diagnostics are disabled or no bins are selected, so the
      ! per-bin pointers are only c_f_pointer'd when n_diag_species > 0.
      ! diagnostic_species_id holds 1-based LOCAL bin positions within the
      ! canonical SEAS1..SEAS5 subset the scheme iterates (species_idx space);
      ! the size-1 dummy is never dereferenced when n_diag_species == 0
      ! (mirrors the settling/dust bridges).  The total fields stay [n_cols].
      integer(c_int), value :: n_diag_species
      real(c_double), intent(in) :: species_density(n_species)
      real(c_double), intent(in) :: species_radius(n_species)
      real(c_double), intent(in) :: species_lower_radius(n_species)
      real(c_double), intent(in) :: species_upper_radius(n_species)
      logical(c_bool), intent(in) :: is_gas_arr(n_species)
      real(c_double), intent(in) :: species_mw_g(n_species)

      integer(c_int), intent(in) :: diagnostic_species_id(max(n_diag_species,1))

      ! Slicing array pointers mapping directly to C++ 8-byte double views
      real(c_double), pointer :: frocean(:), frseaice(:), lat(:), lon(:), sst(:), u10m(:), v10m(:), ustar(:), delp(:,:)
      real(c_double), pointer :: conc(:,:,:), tendency(:,:,:)
      real(c_double), pointer :: diag_mass_total(:), diag_num_total(:)
      real(c_double), pointer :: diag_mass_bin(:,:), diag_num_bin(:,:)

      ! Loop variables
      integer :: icol, i, k
      real(fp) :: dqa, converter
      character(len=64) :: local_scheme
      ! Map from each sea-salt bin (1..n_species) to its column in the full
      ! conc array (1..n_total_species), resolved by trimmed name comparison.
      integer :: target_species(n_species)
      character(len=32) :: bin_names(n_species)

      ! Local physical variables in native precision (fp) to avoid double-float mismatches inside solvers
      real(fp) :: f_frocean, f_frseaice, f_lat, f_lon, f_sst, f_ustar, f_u10m, f_v10m
      real(fp) :: f_delp_layer1
      real(fp) :: f_density(n_species)
      real(fp) :: f_radius(n_species)
      real(fp) :: f_lower_radius(n_species)
      real(fp) :: f_upper_radius(n_species)
      real(fp) :: f_mw_g(n_species)
      real(fp) :: f_conc(1, n_species)
      real(fp) :: f_tendency(1, n_species)

      ! Local diagnostic buffers in fp precision
      real(fp) :: col_mass_total, col_num_total
      real(fp) :: col_mass_bin(max(n_diag_species,1)), col_num_bin(max(n_diag_species,1))

      type(SeaSaltSchemeGONG97Config) :: gong97_config
      type(SeaSaltSchemeGONG03Config) :: gong03_config
      type(SeaSaltSchemeGEOS12Config) :: geos12_config

      ! Convert C string to Fortran
      local_scheme = ""
      icol = 1
      do while (icol < 64)
         if (active_scheme(icol) == c_null_char) exit
         local_scheme(icol:icol) = active_scheme(icol)
         icol = icol + 1
      end do
      local_scheme = trim(adjustl(local_scheme))

      ! Apply the YAML tuning options staged by the C++ process layer onto
      ! the scheme configuration types so the active scheme no longer runs
      ! on compiled defaults.
      gong97_config%scale_factor = real(gong97_scale_factor, fp)
      gong97_config%weibull_flag = (gong97_weibull_flag /= 0)
      gong03_config%scale_factor = real(gong03_scale_factor, fp)
      gong03_config%weibull_flag = (gong03_weibull_flag /= 0)
      geos12_config%scale_factor = real(geos12_scale_factor, fp)
      geos12_config%weibull_flag = (geos12_weibull_flag /= 0)

      ! Associate pointers
      call c_f_pointer(c_frocean,  frocean,  [n_cols])
      call c_f_pointer(c_frseaice, frseaice, [n_cols])
      call c_f_pointer(c_lat,      lat,      [n_cols])
      call c_f_pointer(c_lon,      lon,      [n_cols])
      call c_f_pointer(c_sst,      sst,      [n_cols])
      call c_f_pointer(c_u10m,     u10m,     [n_cols])
      call c_f_pointer(c_v10m,     v10m,     [n_cols])
      call c_f_pointer(c_ustar,    ustar,    [n_cols])
      call c_f_pointer(c_delp,     delp,     [n_cols, n_levels])

      call c_f_pointer(c_conc,     conc,     [n_cols, n_levels, n_total_species])
      call c_f_pointer(c_tendency, tendency, [n_cols, n_levels, n_total_species])

      ! Resolve each sea-salt bin to its slot in the full chemistry array by name.
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
            write(*,'(A,A)') 'FATAL ERROR: SeaSaltScienceBridge cannot resolve sea-salt bin ', trim(bin_names(i))
            call flush(6)
            error stop "FATAL ERROR: SeaSaltScienceBridge unresolved sea-salt bin species"
         end if
      end do

      if (diagnostics /= 0) then
         call c_f_pointer(c_diag_mass_total, diag_mass_total, [n_cols])
         call c_f_pointer(c_diag_num_total,  diag_num_total,  [n_cols])
         if (n_diag_species > 0) then
            call c_f_pointer(c_diag_mass_bin, diag_mass_bin, [n_cols, n_diag_species])
            call c_f_pointer(c_diag_num_bin,  diag_num_bin,  [n_cols, n_diag_species])
         end if
      end if

      ! Cast metadata properties once
      f_density      = real(species_density, fp)
      f_radius       = real(species_radius, fp)
      f_lower_radius = real(species_lower_radius, fp)
      f_upper_radius = real(species_upper_radius, fp)
      f_mw_g         = real(species_mw_g, fp)

      ! Iterate columns
      do icol = 1, n_cols
         f_frocean  = real(frocean(icol), fp)
         f_frseaice = real(frseaice(icol), fp)
         f_lat      = real(lat(icol), fp)
         f_lon      = real(lon(icol), fp)
         f_sst      = real(sst(icol), fp)
         f_u10m     = real(u10m(icol), fp)
         f_v10m     = real(v10m(icol), fp)
         f_ustar    = real(ustar(icol), fp)
         f_delp_layer1 = real(delp(icol, 1), fp)

         ! Gather the sea-salt bin concentrations out of the full unified
         ! array using the name-resolved map (mirrors settling's gather).
         do i = 1, n_species
            f_conc(1, i) = real(conc(icol, 1, target_species(i)), fp)
         end do
         f_tendency(1, :) = 0.0_fp
         col_mass_total   = 0.0_fp
         col_num_total    = 0.0_fp
         col_mass_bin     = 0.0_fp
         col_num_bin      = 0.0_fp

         if (local_scheme == "gong97") then
            call compute_gong97( &
               1, n_species, gong97_config, PI, &
               f_frocean, f_frseaice, f_lat, f_lon, f_sst, f_u10m, f_v10m, &
               f_density, f_radius, f_lower_radius, f_upper_radius, &
               f_conc, f_tendency, &
               seasalt_mass_emission_total=col_mass_total, &
               seasalt_number_emission_total=col_num_total, &
               seasalt_mass_emission_per_bin=col_mass_bin, &
               seasalt_number_emission_per_bin=col_num_bin, &
               diagnostic_species_id=diagnostic_species_id)
         else if (local_scheme == "gong03") then
            call compute_gong03( &
               1, n_species, gong03_config, PI, &
               f_frocean, f_frseaice, f_lat, f_lon, f_sst, f_u10m, f_v10m, &
               f_density, f_radius, f_lower_radius, f_upper_radius, &
               f_conc, f_tendency, &
               seasalt_mass_emission_total=col_mass_total, &
               seasalt_number_emission_total=col_num_total, &
               seasalt_mass_emission_per_bin=col_mass_bin, &
               seasalt_number_emission_per_bin=col_num_bin, &
               diagnostic_species_id=diagnostic_species_id)
         else if (local_scheme == "geos12") then
            call compute_geos12( &
               1, n_species, geos12_config, PI, &
               f_frocean, f_frseaice, f_lat, f_lon, f_sst, f_u10m, f_ustar, f_v10m, &
               f_density, f_radius, f_lower_radius, f_upper_radius, &
               f_conc, f_tendency, &
               seasalt_mass_emission_total=col_mass_total, &
               seasalt_number_emission_total=col_num_total, &
               seasalt_mass_emission_per_bin=col_mass_bin, &
               seasalt_number_emission_per_bin=col_num_bin, &
               diagnostic_species_id=diagnostic_species_id)
         end if

         ! Apply calculated emission fluxes to surface concentration (layer 1),
         ! scattering back into the full unified array via the name map.
         do i = 1, n_species
            tendency(icol, 1, target_species(i)) = f_tendency(1, i)

            ! dqa = species_tendencies(1, i) * timestep * g0 / met%DELP(1)
            dqa = f_tendency(1, i) * real(dt, fp) * g0 / f_delp_layer1

            ! converter
            if (is_gas_arr(i)) then
               converter = AIRMW / f_mw_g(i) * 1.0e6_fp
            else
               converter = 1.0e9_fp
            end if
            dqa = dqa * converter

            conc(icol, 1, target_species(i)) = &
               max(0.0_c_double, conc(icol, 1, target_species(i)) + real(dqa, c_double))
         end do

         ! Write diagnostics back to C++ pointers with double casting.  The
         ! per-bin arrays are n_diag_species wide (the selected subset); the
         ! scheme already scattered each bin into its diag_idx slot, so the
         ! copy is a straight 1:1 column write.
         if (diagnostics /= 0) then
            diag_mass_total(icol) = real(col_mass_total, c_double)
            diag_num_total(icol)  = real(col_num_total, c_double)
            if (n_diag_species > 0) then
               diag_mass_bin(icol, 1:n_diag_species) = real(col_mass_bin(1:n_diag_species), c_double)
               diag_num_bin(icol, 1:n_diag_species)  = real(col_num_bin(1:n_diag_species), c_double)
            end if
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
   end subroutine run_seasalt_science_bridge

end module SeaSaltScienceBridge_Mod
