module WetDepScienceBridge_Mod
   use iso_c_binding, only: c_ptr, c_f_pointer, c_double, c_char, c_associated, c_bool, c_int
   use catchem_bridge_precision, only: fp
   use catchem_bridge_constants, only: g0, AIRMW
   use WetDepCommon_Mod, only: WetDepSchemeJACOBConfig
   use WetDepScheme_JACOB_Mod, only: compute_jacob
   implicit none
contains

   subroutine run_wetdep_science_bridge( &
      n_cols, n_levels, n_species, dt, &
      diagnostics, &
      jacob_scale_factor, jacob_radius_threshold, jacob_so4_gocart_resusp, jacob_so4_washout_eff, &
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

      ! Scheme tuning options staged by WetDepProcess::init from the runtime
      ! YAML.  The C++ layer owns parsing and validation; the bridge only
      ! applies them onto the Jacob configuration type.
      real(c_double), value :: jacob_scale_factor, jacob_radius_threshold, jacob_so4_washout_eff
      integer(c_int), value :: jacob_so4_gocart_resusp

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
      ! Precipitation fluxes live on vertical interfaces (the host's 0:nlev
      ! index range maps to this Fortran 1:nlev+1 range).
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

      ! Apply the YAML tuning options staged by the C++ process layer onto
      ! the Jacob configuration so wet deposition no longer runs on compiled
      ! defaults.
      jacob_config%scale_factor = real(jacob_scale_factor, fp)
      jacob_config%radius_threshold = real(jacob_radius_threshold, fp)
      jacob_config%so4_gocart_resusp = (jacob_so4_gocart_resusp /= 0)
      jacob_config%so4_washout_eff = real(jacob_so4_washout_eff, fp)

      ! Map pointers
      nullify(airden_dry, mairden, pedge, pfilsan, pfllsan, reevapls, t_air)
      nullify(conc, tendency, diag_mass, diag_flux)

      if (.not. c_associated(c_airden_dry) .and. .not. c_associated(c_mairden)) then
         write(*,'(A)') 'FATAL ERROR: WetDepScienceBridge missing required field AIRDEN / AIRDEN_DRY'
         call flush(6)
         error stop "FATAL ERROR: WetDepScienceBridge missing required field AIRDEN"
      end if
      if (.not. c_associated(c_pedge)) then
         write(*,'(A)') 'FATAL ERROR: WetDepScienceBridge missing required field PEDGE'
         call flush(6)
         error stop "FATAL ERROR: WetDepScienceBridge missing required field PEDGE"
      end if
      if (.not. c_associated(c_t_air)) then
         write(*,'(A)') 'FATAL ERROR: WetDepScienceBridge missing required field T'
         call flush(6)
         error stop "FATAL ERROR: WetDepScienceBridge missing required field T"
      end if
      if (.not. c_associated(c_conc) .or. .not. c_associated(c_tendency)) then
         write(*,'(A)') 'FATAL ERROR: WetDepScienceBridge missing required concentration or tendency pointers'
         call flush(6)
         error stop "FATAL ERROR: WetDepScienceBridge missing required concentration or tendency pointers"
      end if

      if (c_associated(c_airden_dry)) call c_f_pointer(c_airden_dry, airden_dry, [n_cols, n_levels])
      if (c_associated(c_mairden))    call c_f_pointer(c_mairden,    mairden,    [n_cols, n_levels])
      if (c_associated(c_pedge))      call c_f_pointer(c_pedge,      pedge,      [n_cols, n_levels+1])
      if (c_associated(c_pfilsan))    call c_f_pointer(c_pfilsan,    pfilsan,    [n_cols, n_levels+1])
      if (c_associated(c_pfllsan))    call c_f_pointer(c_pfllsan,    pfllsan,    [n_cols, n_levels+1])
      if (c_associated(c_reevapls))   call c_f_pointer(c_reevapls,   reevapls,   [n_cols, n_levels])
      if (c_associated(c_t_air))      call c_f_pointer(c_t_air,      t_air,      [n_cols, n_levels])

      call c_f_pointer(c_conc,     conc,     [n_cols, n_levels, n_species])
      call c_f_pointer(c_tendency, tendency, [n_cols, n_levels, n_species])

      if (diagnostics /= 0) then
         if (c_associated(c_diag_mass)) call c_f_pointer(c_diag_mass, diag_mass, [n_cols, n_levels, n_species])
         if (c_associated(c_diag_flux)) call c_f_pointer(c_diag_flux, diag_flux, [n_cols, n_levels, n_species])
      end if

      ! Extract real species names from flat char array passed via BIND(C)
      do i = 1, n_species
         dummy_sp_names(i) = ""
         do j = 1, 32
            dummy_sp_names(i)(j:j) = species_names(j, i)
         end do
         dummy_sp_names(i) = trim(adjustl(dummy_sp_names(i)))
      end do

      ! Copy to standard logical arrays & cast doubles once
      f_is_aerosol      = species_is_aerosol
      f_henry_cr        = real(species_henry_cr, fp)
      f_henry_k0        = real(species_henry_k0, fp)
      f_henry_pKa       = real(species_henry_pKa, fp)
      f_wd_retfactor    = real(species_wd_retfactor, fp)
      f_wd_LiqAndGas    = species_wd_LiqAndGas
      f_wd_convfacI2G   = real(species_wd_convfacI2G, fp)
      f_wd_rainouteff   = real(species_wd_rainouteff, fp)
      col_wd_reevap_frac = real(species_wd_reevap_frac, fp)
      f_radius          = real(species_radius, fp)
      f_mw_g            = real(species_mw_g, fp)

      ! Iterate columns
      do icol = 1, n_cols
         if (associated(airden_dry)) then
            f_airden_dry = real(airden_dry(icol, :), fp)
         else
            f_airden_dry = real(mairden(icol, :), fp)
         end if

         if (associated(mairden)) then
            f_mairden = real(mairden(icol, :), fp)
         else
            f_mairden = f_airden_dry
         end if

         f_pedge = real(pedge(icol, :), fp)
         f_t_air = real(t_air(icol, :), fp)

         if (associated(pfilsan))    then; f_pfilsan  = real(pfilsan(icol, :), fp);  else; f_pfilsan  = 0.0_fp; end if
         if (associated(pfllsan))    then; f_pfllsan  = real(pfllsan(icol, :), fp);  else; f_pfllsan  = 0.0_fp; end if
         if (associated(reevapls))   then; f_reevapls = real(reevapls(icol, :), fp); else; f_reevapls = 0.0_fp; end if

         ! Extract input concentrations (already in ug/kg for aerosols, ppmv for gases)
         do ispec = 1, n_species
            f_conc(:, ispec) = real(conc(icol, :, ispec), fp)
         end do

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
            diagnostic_species_id=diagnostic_species_id)

         ! compute_jacob preserves the legacy CATChem replacement contract:
         ! species_tendencies contains the post-wetdep concentration in the
         ! native aerosol (ug/kg) or gas (ppmv) units, rather than a rate.
         ! The upstream ProcessWetDep interface assigns this value directly.
         do ispec = 1, n_species
            tendency(icol, :, ispec) = real(col_tendencies(:, ispec), c_double)
            conc(icol, :, ispec) = real(col_tendencies(:, ispec), c_double)
         end do

         if (diagnostics /= 0) then
            diag_mass(icol, :, :) = real(col_diag_mass, c_double)
            diag_flux(icol, :, :) = real(col_diag_flux, c_double)
         end if
      end do

   end subroutine run_wetdep_science_bridge

end module WetDepScienceBridge_Mod
