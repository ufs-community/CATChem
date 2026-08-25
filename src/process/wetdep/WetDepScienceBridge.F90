module WetDepScienceBridge_Mod
   use iso_c_binding, only: c_ptr, c_f_pointer, c_double, c_char, c_associated, c_bool, c_int, c_null_char
   use Precision_Mod, only: fp
   use Constants, only: g0, AIRMW
   use WetDepCommon_Mod, only: WetDepSchemeJACOBConfig
   use WetDepScheme_JACOB_Mod, only: compute_jacob
   implicit none
contains

   subroutine run_wetdep_science_bridge( &
      n_cols, n_levels, n_species, dt, &
      diagnostics, &
   ! 3D Met Pointers
      c_airden_dry, c_mairden, c_pedge, c_pfilsan, c_pfllsan, c_reevapls, c_t_air, &
   ! Metadata
      species_is_aerosol, species_henry_cr, species_henry_k0, species_henry_pKa, &
      species_wd_retfactor, species_wd_LiqAndGas, species_wd_convfacI2G, species_wd_rainouteff, species_wd_reevap_frac, &
      species_radius, species_mw_g, species_names, &
   ! Concentrations, Tendencies & Diagnostics
      c_conc, c_tendency, c_diag_mass, c_diag_flux, &
      diagnostic_species_id, n_diag_species &
      ) bind(C, name="run_wetdep_science_bridge")

      integer(c_int), value :: n_cols, n_levels, n_species
      real(c_double), value :: dt
      integer(c_int), value :: diagnostics

      ! C pointers
      type(c_ptr), value :: c_airden_dry, c_mairden, c_pedge, c_pfilsan, c_pfllsan, c_reevapls, c_t_air
      type(c_ptr), value :: c_conc, c_tendency, c_diag_mass, c_diag_flux

      ! Metadata dummy arrays in double precision to match C++
      logical(c_bool), intent(in) :: species_is_aerosol(n_species)
      real(c_double), intent(in) :: species_henry_cr(n_species)
      real(c_double), intent(in) :: species_henry_k0(n_species)
      real(c_double), intent(in) :: species_henry_pKa(n_species)
      real(c_double), intent(in) :: species_wd_retfactor(n_species)
      logical(c_bool), intent(in) :: species_wd_LiqAndGas(n_species)
      real(c_double), intent(in) :: species_wd_convfacI2G(n_species)
      real(c_double), intent(in) :: species_wd_rainouteff(n_species, 3) ! Exactly 3 elements matching GOCART/Jacob specifications
      real(c_double), intent(in) :: species_wd_reevap_frac(n_species)
      real(c_double), intent(in) :: species_radius(n_species)
      real(c_double), intent(in) :: species_mw_g(n_species)
      character(kind=c_char), intent(in) :: species_names(32, n_species)

      integer(c_int), value :: n_diag_species
      integer(c_int), intent(in) :: diagnostic_species_id(n_diag_species)

      ! Slicing array pointers pointing directly to double precision (c_double) C++ views
      real(c_double), pointer :: airden_dry(:,:), mairden(:,:), pedge(:,:), pfilsan(:,:), pfllsan(:,:), reevapls(:,:), t_air(:,:)
      real(c_double), pointer :: conc(:,:,:), tendency(:,:,:)
      real(c_double), pointer :: diag_mass(:,:,:), diag_flux(:,:,:)

      ! Loop variables
      integer :: icol, i, j, ispec
      character(len=32) :: dummy_sp_names(n_species)

      ! Local arrays in native solver precision (fp) to avoid double-float mismatches
      real(fp) :: f_airden_dry(n_levels)
      real(fp) :: f_mairden(n_levels)
      real(fp) :: f_pedge(n_levels+1)
      real(fp) :: f_pfilsan(n_levels+1)
      real(fp) :: f_pfllsan(n_levels+1)
      real(fp) :: f_reevapls(n_levels)
      real(fp) :: f_t_air(n_levels)

      ! Casted metadata properties
      logical :: f_is_aerosol(n_species)
      real(fp) :: f_henry_cr(n_species)
      real(fp) :: f_henry_k0(n_species)
      real(fp) :: f_henry_pKa(n_species)
      real(fp) :: f_wd_retfactor(n_species)
      logical :: f_wd_LiqAndGas(n_species)
      real(fp) :: f_wd_convfacI2G(n_species)
      real(fp) :: f_wd_rainouteff(n_species, 3) ! Shaped with species as dimension 1, and 3 efficiency factors
      real(fp) :: col_wd_reevap_frac(n_species)
      real(fp) :: f_radius(n_species)
      real(fp) :: f_mw_g(n_species)

      ! Sliced concentration and tendencies in solver precision
      real(fp) :: f_conc(n_levels, n_species)
      real(fp) :: col_tendencies(n_levels, n_species)
      real(fp) :: col_diag_mass(n_levels, n_species)
      real(fp) :: col_diag_flux(n_levels, n_species)

      type(WetDepSchemeJACOBConfig) :: jacob_config

      ! Map pointers
      call c_f_pointer(c_airden_dry, airden_dry, [n_cols, n_levels])
      call c_f_pointer(c_mairden,    mairden,    [n_cols, n_levels])
      call c_f_pointer(c_pedge,      pedge,      [n_cols, n_levels+1])
      call c_f_pointer(c_pfilsan,    pfilsan,    [n_cols, n_levels+1])
      call c_f_pointer(c_pfllsan,    pfllsan,    [n_cols, n_levels+1])
      call c_f_pointer(c_reevapls,   reevapls,   [n_cols, n_levels])
      call c_f_pointer(c_t_air,      t_air,      [n_cols, n_levels])

      call c_f_pointer(c_conc,     conc,     [n_cols, n_levels, n_species])
      call c_f_pointer(c_tendency, tendency, [n_cols, n_levels, n_species])

      if (diagnostics /= 0) then
         call c_f_pointer(c_diag_mass, diag_mass, [n_cols, n_levels, n_species])
         call c_f_pointer(c_diag_flux, diag_flux, [n_cols, n_levels, n_species])
      end if

      ! Extract real species names from flat char array passed via BIND(C)
      do i = 1, n_species
         dummy_sp_names(i) = ""
         do j = 1, 32
            if (species_names(j, i) == c_null_char) exit
            dummy_sp_names(i)(j:j) = species_names(j, i)
         end do
         dummy_sp_names(i) = trim(adjustl(dummy_sp_names(i)))
      end do

      ! Copy to standard logical arrays & cast doubles once
      do i = 1, n_species
         f_is_aerosol(i) = species_is_aerosol(i)
         f_wd_LiqAndGas(i) = species_wd_LiqAndGas(i)
      end do
      f_henry_cr        = real(species_henry_cr, fp)
      f_henry_k0        = real(species_henry_k0, fp)
      f_henry_pKa       = real(species_henry_pKa, fp)
      f_wd_retfactor    = real(species_wd_retfactor, fp)
      f_wd_convfacI2G   = real(species_wd_convfacI2G, fp)
      f_wd_rainouteff   = real(species_wd_rainouteff, fp)
      col_wd_reevap_frac = real(species_wd_reevap_frac, fp)
      f_radius          = real(species_radius, fp)
      f_mw_g            = real(species_mw_g, fp)

      ! Iterate columns
      do icol = 1, n_cols
         f_airden_dry   = real(airden_dry(icol, :), fp)
         f_mairden      = real(mairden(icol, :), fp)
         f_pedge        = real(pedge(icol, :), fp)
         f_pfilsan      = real(pfilsan(icol, :), fp)
         f_pfllsan      = real(pfllsan(icol, :), fp)
         f_reevapls     = real(reevapls(icol, :), fp)
         f_t_air        = real(t_air(icol, :), fp)

         f_conc         = real(conc(icol, :, :), fp)
         col_tendencies = 0.0_fp
         col_diag_mass  = 0.0_fp
         col_diag_flux  = 0.0_fp

         ! Execute JACOB scheme
         call compute_jacob( &
            n_levels, n_species, jacob_config, &
            f_airden_dry, f_mairden, f_pedge, f_pfilsan, f_pfllsan, &
            f_reevapls, f_t_air, real(dt, fp), &
            f_is_aerosol, dummy_sp_names, f_henry_cr, f_henry_k0, f_henry_pKa, &
            f_wd_retfactor, f_wd_LiqAndGas, f_wd_convfacI2G, f_wd_rainouteff, &
            col_wd_reevap_frac, f_radius, f_mw_g, &
            f_conc, col_tendencies, &
            wetdep_mass_per_species_per_level=col_diag_mass, &
            wetdep_flux_per_species_per_level=col_diag_flux, &
            diagnostic_species_id=diagnostic_species_id(1:n_diag_species))

         ! Convert tendencies from process-specific units (ug/kg/s or ppm/s) to kg/kg/s
         do ispec = 1, n_species
            if ( f_is_aerosol(ispec) ) then
               col_tendencies(:, ispec) = col_tendencies(:, ispec) * 1.0e-9_fp
            else
               col_tendencies(:, ispec) = col_tendencies(:, ispec) * 1.0e-6_fp * (f_mw_g(ispec) / AIRMW)
            end if
         end do

         ! Write tendencies and concentrations back in-place (casting to c_double)
         tendency(icol, :, :) = real(col_tendencies, c_double)
         conc(icol, :, :) = conc(icol, :, :) + real(dt * col_tendencies, c_double)

         if (diagnostics /= 0) then
            diag_mass(icol, :, :) = real(col_diag_mass, c_double)
            diag_flux(icol, :, :) = real(col_diag_flux, c_double)
         end if
      end do

   end subroutine run_wetdep_science_bridge

end module WetDepScienceBridge_Mod
