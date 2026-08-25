module SeaSaltScienceBridge_Mod
   use iso_c_binding, only: c_ptr, c_f_pointer, c_double, c_char, c_associated, c_null_char, c_bool, c_int
   use Precision_Mod, only: fp
   use Constants, only: g0, AIRMW, PI
   use SeaSaltCommon_Mod, only: SeaSaltSchemeGONG97Config, SeaSaltSchemeGONG03Config, SeaSaltSchemeGEOS12Config
   use SeaSaltScheme_GONG97_Mod, only: compute_gong97
   use SeaSaltScheme_GONG03_Mod, only: compute_gong03
   use SeaSaltScheme_GEOS12_Mod, only: compute_geos12
   implicit none
contains

   subroutine run_seasalt_science_bridge( &
      n_cols, n_levels, n_species, dt, &
      active_scheme, diagnostics, &
   ! Met Pointers
      c_frocean, c_frseaice, c_lat, c_lon, c_sst, c_u10m, c_v10m, c_ustar, c_delp, &
   ! Species Metadata
      species_density, species_radius, species_lower_radius, species_upper_radius, is_gas_arr, species_mw_g, &
   ! Concentrations and Tendency
      c_conc, c_tendency, &
   ! Diagnostics
      c_diag_mass_total, c_diag_num_total, c_diag_mass_bin, c_diag_num_bin, &
      diagnostic_species_id, n_diag_species &
      ) bind(C, name="run_seasalt_science_bridge")

      integer(c_int), value :: n_cols, n_levels, n_species
      real(c_double), value :: dt
      character(kind=c_char), intent(in) :: active_scheme(*)
      integer(c_int), value :: diagnostics

      type(c_ptr), value :: c_frocean, c_frseaice, c_lat, c_lon, c_sst, c_u10m, c_v10m, c_ustar, c_delp
      type(c_ptr), value :: c_conc, c_tendency
      type(c_ptr), value :: c_diag_mass_total, c_diag_num_total, c_diag_mass_bin, c_diag_num_bin

      integer(c_int), value :: n_diag_species
      real(c_double), intent(in) :: species_density(n_species)
      real(c_double), intent(in) :: species_radius(n_species)
      real(c_double), intent(in) :: species_lower_radius(n_species)
      real(c_double), intent(in) :: species_upper_radius(n_species)
      logical(c_bool), intent(in) :: is_gas_arr(n_species)
      real(c_double), intent(in) :: species_mw_g(n_species)

      integer(c_int), intent(in) :: diagnostic_species_id(n_diag_species)

      ! Slicing array pointers mapping directly to C++ 8-byte double views
      real(c_double), pointer :: frocean(:), frseaice(:), lat(:), lon(:), sst(:), u10m(:), v10m(:), ustar(:), delp(:,:)
      real(c_double), pointer :: conc(:,:,:), tendency(:,:,:)
      real(c_double), pointer :: diag_mass_total(:), diag_num_total(:), diag_mass_bin(:,:), diag_num_bin(:,:)

      ! Loop variables
      integer :: icol, i
      real(fp) :: dqa, converter
      character(len=64) :: local_scheme

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
      real(fp) :: col_mass_bin(n_species), col_num_bin(n_species)

      type(SeaSaltSchemeGONG97Config) :: gong97_config
      type(SeaSaltSchemeGONG03Config) :: gong03_config
      type(SeaSaltSchemeGEOS12Config) :: geos12_config

      ! Convert C string to Fortran
      local_scheme = ""
      icol = 1
      do while (active_scheme(icol) /= c_null_char .and. icol < 64)
         local_scheme(icol:icol) = active_scheme(icol)
         icol = icol + 1
      end do
      local_scheme = trim(adjustl(local_scheme))

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

      call c_f_pointer(c_conc,     conc,     [n_cols, n_levels, n_species])
      call c_f_pointer(c_tendency, tendency, [n_cols, n_levels, n_species])

      if (diagnostics /= 0) then
         call c_f_pointer(c_diag_mass_total, diag_mass_total, [n_cols])
         call c_f_pointer(c_diag_num_total,  diag_num_total,  [n_cols])
         call c_f_pointer(c_diag_mass_bin,   diag_mass_bin,   [n_cols, n_species])
         call c_f_pointer(c_diag_num_bin,    diag_num_bin,    [n_cols, n_species])
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

         ! Cast concentrations
         f_conc(1, :)     = real(conc(icol, 1, :), fp)
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

         ! Apply calculated emission fluxes to surface concentration (layer 1)
         do i = 1, n_species
            tendency(icol, 1, i) = f_tendency(1, i)

            ! dqa = species_tendencies(1, i) * timestep * g0 / met%DELP(1)
            dqa = f_tendency(1, i) * real(dt, fp) * g0 / f_delp_layer1

            ! converter
            if (is_gas_arr(i)) then
               converter = AIRMW / f_mw_g(i) * 1.0e6_fp
            else
               converter = 1.0e9_fp
            end if
            dqa = dqa * converter

            conc(icol, 1, i) = conc(icol, 1, i) + dqa
         end do

         ! Write diagnostics back to C++ pointers with double casting
         if (diagnostics /= 0) then
            diag_mass_total(icol) = real(col_mass_total, c_double)
            diag_num_total(icol)  = real(col_num_total, c_double)
            diag_mass_bin(icol, :) = real(col_mass_bin, c_double)
            diag_num_bin(icol, :)  = real(col_num_bin, c_double)
         end if
      end do

   end subroutine run_seasalt_science_bridge

end module SeaSaltScienceBridge_Mod
