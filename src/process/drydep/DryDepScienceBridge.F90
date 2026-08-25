module DryDepScienceBridge_Mod
   use iso_c_binding, only: c_ptr, c_f_pointer, c_double, c_char, c_associated, c_null_char, c_bool, c_int
   use Precision_Mod, only: fp
   use DryDepCommon_Mod, only: DryDepSchemeWESELYConfig, DryDepSchemeGOCARTConfig, DryDepSchemeZHANGConfig
   use DryDepScheme_WESELY_Mod, only: compute_wesely
   use DryDepScheme_GOCART_Mod, only: compute_gocart
   use DryDepScheme_ZHANG_Mod, only: compute_zhang
   implicit none
contains

   subroutine run_drydep_science_bridge( &
      n_cols, n_levels, n_species, dt, &
      gas_scheme, aero_scheme, diagnostics, &
   ! 3D Met Pointers
      c_bxheight, c_airden, c_t_air, c_z_edges, c_rh, &
   ! 2D/1D Met Pointers
      c_cldfrc, c_frlai, c_frlanduse, c_iland, c_is_ice, c_is_land, c_is_snow, &
      c_lat, c_lon, c_obk, c_ps, c_salinity, c_suncosmid, c_swgdn, c_ts, c_tskin, &
      c_ustar, c_z0, c_frlake, c_gwettop, c_hflux, c_lwi, c_pblh, c_u10m, c_v10m, c_z0h, &
   ! Species Metadata Dummy Arrays
      species_mw_g, species_dd_f0, species_dd_hstar, species_dd_DvzAerSnow, &
      species_dd_DvzMinVal_snow, species_dd_DvzMinVal_land, species_density, &
      species_radius, species_is_seasalt, species_is_dust, species_lower_radius, &
      species_upper_radius, is_gas_arr, &
   ! Concentrations, Tendencies & Diagnostics
      c_conc, c_tendency, c_diag_con, c_diag_vel, &
      diagnostic_species_id, n_diag_species &
      ) bind(C, name="run_drydep_science_bridge")

      integer(c_int), value :: n_cols, n_levels, n_species
      real(c_double), value :: dt
      character(kind=c_char), intent(in) :: gas_scheme(*)
      character(kind=c_char), intent(in) :: aero_scheme(*)
      integer(c_int), value :: diagnostics

      ! C pointers
      type(c_ptr), value :: c_bxheight, c_airden, c_t_air, c_z_edges, c_rh
      type(c_ptr), value :: c_cldfrc, c_frlai, c_frlanduse, c_iland, c_is_ice, c_is_land, c_is_snow
      type(c_ptr), value :: c_lat, c_lon, c_obk, c_ps, c_salinity, c_suncosmid, c_swgdn, c_ts, c_tskin
      type(c_ptr), value :: c_ustar, c_z0, c_frlake, c_gwettop, c_hflux, c_lwi, c_pblh, c_u10m, c_v10m, c_z0h
      type(c_ptr), value :: c_conc, c_tendency, c_diag_con, c_diag_vel

      ! Metadata dummy arrays in double precision to match C++ doubles
      real(c_double), intent(in) :: species_mw_g(n_species)
      real(c_double), intent(in) :: species_dd_f0(n_species)
      real(c_double), intent(in) :: species_dd_hstar(n_species)
      real(c_double), intent(in) :: species_dd_DvzAerSnow(n_species)
      real(c_double), intent(in) :: species_dd_DvzMinVal_snow(n_species)
      real(c_double), intent(in) :: species_dd_DvzMinVal_land(n_species)
      real(c_double), intent(in) :: species_density(n_species)
      real(c_double), intent(in) :: species_radius(n_species)
      logical(c_bool), intent(in) :: species_is_seasalt(n_species)
      logical(c_bool), intent(in) :: species_is_dust(n_species)
      real(c_double), intent(in) :: species_lower_radius(n_species)
      real(c_double), intent(in) :: species_upper_radius(n_species)
      logical(c_bool), intent(in) :: is_gas_arr(n_species)
      integer(c_int), value :: n_diag_species
      integer(c_int), intent(in) :: diagnostic_species_id(n_diag_species)

      ! Slicing array pointers pointing directly to double precision (c_double) C++ views
      real(c_double), pointer :: bxheight(:,:), airden(:,:), t_air(:,:), z_edges(:,:), rh(:,:)
      real(c_double), pointer :: cldfrc(:), frlai(:,:,:), frlanduse(:,:,:), lat(:), lon(:)
      integer, pointer :: iland(:,:,:)
      logical(c_bool), pointer :: is_ice(:), is_land(:), is_snow(:)
      real(c_double), pointer :: obk(:), ps(:), salinity(:), suncosmid(:), swgdn(:), ts(:), tskin(:)
      real(c_double), pointer :: ustar(:), z0(:), frlake(:), gwettop(:), hflux(:)
      integer, pointer :: lwi(:)
      real(c_double), pointer :: pblh(:), u10m(:), v10m(:), z0h(:)

      real(c_double), pointer :: conc(:,:,:), tendency(:,:,:), diag_con(:,:), diag_vel(:,:)

      ! Loop variables
      integer :: icol
      character(len=64) :: local_gas, local_aero
      character(len=255) :: local_lucname = "NOAH"
      character(len=30) :: dummy_sp_names(n_species)

      ! Local arrays in native solver precision (fp) to avoid double-float mismatches
      real(fp) :: f_bxheight(1), f_airden(1), f_t_air(1), f_z_edges(2), f_rh(1)
      real(fp) :: f_cldfrc, f_frlai(20), f_frlanduse(20), f_lat, f_lon
      integer :: f_iland(20)
      logical :: f_is_ice, f_is_land, f_is_snow
      real(fp) :: f_obk, f_ps, f_salinity, f_suncosmid, f_swgdn, f_ts, f_tskin
      real(fp) :: f_ustar, f_z0, f_frlake, f_gwettop, f_hflux, f_pblh, f_u10m, f_v10m, f_z0h
      integer :: f_lwi

      ! Casted metadata properties
      real(fp) :: f_mw_g(n_species)
      real(fp) :: f_dd_f0(n_species)
      real(fp) :: f_dd_hstar(n_species)
      real(fp) :: f_dd_DvzAerSnow(n_species)
      real(fp) :: f_dd_DvzMinVal_snow(n_species)
      real(fp) :: f_dd_DvzMinVal_land(n_species)
      real(fp) :: f_density(n_species)
      real(fp) :: f_radius(n_species)
      real(fp) :: f_lower_radius(n_species)
      real(fp) :: f_upper_radius(n_species)
      logical  :: f_is_seasalt(n_species)
      logical  :: f_is_dust(n_species)
      logical  :: f_is_gas_arr(n_species)

      ! Sliced concentration and tendencies in solver precision
      real(fp) :: f_conc(1, n_species)
      real(fp) :: col_tendencies(1, n_species)
      real(fp) :: col_diag_con(n_species)
      real(fp) :: col_diag_vel(n_species)

      type(DryDepSchemeWESELYConfig) :: wesely_config
      type(DryDepSchemeGOCARTConfig) :: gocart_config
      type(DryDepSchemeZHANGConfig) :: zhang_config

      ! Convert C strings to Fortran strings
      icol = 1
      do while (gas_scheme(icol) /= c_null_char .and. icol < 64)
         local_gas(icol:icol) = gas_scheme(icol)
         icol = icol + 1
      end do
      local_gas = trim(adjustl(local_gas))

      icol = 1
      do while (aero_scheme(icol) /= c_null_char .and. icol < 64)
         local_aero(icol:icol) = aero_scheme(icol)
         icol = icol + 1
      end do
      local_aero = trim(adjustl(local_aero))

      ! Associate pointers
      call c_f_pointer(c_bxheight, bxheight, [n_cols, n_levels])
      call c_f_pointer(c_airden,   airden,   [n_cols, n_levels])
      call c_f_pointer(c_t_air,    t_air,    [n_cols, n_levels])
      call c_f_pointer(c_z_edges,  z_edges,  [n_cols, n_levels+1])
      call c_f_pointer(c_rh,       rh,       [n_cols, n_levels])

      call c_f_pointer(c_cldfrc,    cldfrc,    [n_cols])
      call c_f_pointer(c_frlai,     frlai,     [n_cols, 1, 20]) ! n_landuse = 20
      call c_f_pointer(c_frlanduse, frlanduse, [n_cols, 1, 20])
      call c_f_pointer(c_iland,     iland,     [n_cols, 1, 20])
      call c_f_pointer(c_is_ice,    is_ice,    [n_cols])
      call c_f_pointer(c_is_land,   is_land,   [n_cols])
      call c_f_pointer(c_is_snow,   is_snow,   [n_cols])
      call c_f_pointer(c_lat,       lat,       [n_cols])
      call c_f_pointer(c_lon,       lon,       [n_cols])
      call c_f_pointer(c_obk,       obk,       [n_cols])
      call c_f_pointer(c_ps,        ps,        [n_cols])
      call c_f_pointer(c_salinity,  salinity,  [n_cols])
      call c_f_pointer(c_suncosmid, suncosmid, [n_cols])
      call c_f_pointer(c_swgdn,     swgdn,     [n_cols])
      call c_f_pointer(c_ts,        ts,        [n_cols])
      call c_f_pointer(c_tskin,     tskin,     [n_cols])
      call c_f_pointer(c_ustar,     ustar,     [n_cols])
      call c_f_pointer(c_z0,        z0,        [n_cols])
      call c_f_pointer(c_frlake,    frlake,    [n_cols])
      call c_f_pointer(c_gwettop,   gwettop,   [n_cols])
      call c_f_pointer(c_hflux,     hflux,     [n_cols])
      call c_f_pointer(c_lwi,       lwi,       [n_cols])
      call c_f_pointer(c_pblh,      pblh,      [n_cols])
      call c_f_pointer(c_u10m,      u10m,      [n_cols])
      call c_f_pointer(c_v10m,      v10m,      [n_cols])
      call c_f_pointer(c_z0h,       z0h,       [n_cols])

      call c_f_pointer(c_conc,     conc,     [n_cols, n_levels, n_species])
      call c_f_pointer(c_tendency, tendency, [n_cols, n_levels, n_species])

      if (diagnostics /= 0) then
         call c_f_pointer(c_diag_con, diag_con, [n_cols, n_species])
         call c_f_pointer(c_diag_vel, diag_vel, [n_cols, n_species])
      endif

      dummy_sp_names = "UNKNOWN"

      ! Copy metadata once
      f_mw_g = real(species_mw_g, fp)
      f_dd_f0 = real(species_dd_f0, fp)
      f_dd_hstar = real(species_dd_hstar, fp)
      f_dd_DvzAerSnow = real(species_dd_DvzAerSnow, fp)
      f_dd_DvzMinVal_snow = real(species_dd_DvzMinVal_snow, fp)
      f_dd_DvzMinVal_land = real(species_dd_DvzMinVal_land, fp)
      f_density = real(species_density, fp)
      f_radius = real(species_radius, fp)
      f_lower_radius = real(species_lower_radius, fp)
      f_upper_radius = real(species_upper_radius, fp)
      f_is_seasalt = species_is_seasalt
      f_is_dust = species_is_dust
      f_is_gas_arr = is_gas_arr

      ! Iterate columns and slice
      do icol = 1, n_cols
         ! Cast scalars and slices
         f_bxheight(1) = real(bxheight(icol, 1), fp)
         f_airden(1)   = real(airden(icol, 1), fp)
         f_t_air(1)    = real(t_air(icol, 1), fp)
         f_z_edges(1:2) = real(z_edges(icol, 1:2), fp)
         f_rh(1)       = real(rh(icol, 1), fp)

         f_cldfrc     = real(cldfrc(icol), fp)
         f_frlai      = real(frlai(icol, 1, :), fp)
         f_frlanduse  = real(frlanduse(icol, 1, :), fp)
         f_iland      = iland(icol, 1, :)
         f_is_ice     = is_ice(icol)
         f_is_land    = is_land(icol)
         f_is_snow    = is_snow(icol)
         f_lat        = real(lat(icol), fp)
         f_lon        = real(lon(icol), fp)
         f_obk        = real(obk(icol), fp)
         f_ps         = real(ps(icol), fp)
         f_salinity   = real(salinity(icol), fp)
         f_suncosmid  = real(suncosmid(icol), fp)
         f_swgdn      = real(swgdn(icol), fp)
         f_ts         = real(ts(icol), fp)
         f_tskin      = real(tskin(icol), fp)
         f_ustar      = real(ustar(icol), fp)
         f_z0         = real(z0(icol), fp)
         f_frlake     = real(frlake(icol), fp)
         f_gwettop    = real(gwettop(icol), fp)
         f_hflux      = real(hflux(icol), fp)
         f_lwi        = lwi(icol)
         f_pblh       = real(pblh(icol), fp)
         f_u10m       = real(u10m(icol), fp)
         f_v10m       = real(v10m(icol), fp)
         f_z0h        = real(z0h(icol), fp)

         f_conc(1, :)   = real(conc(icol, 1, :), fp)
         col_tendencies = 0.0_fp
         col_diag_con   = 0.0_fp
         col_diag_vel   = 0.0_fp

         ! Execute GAS schemes
         if (trim(local_gas) == "wesely") then
            call compute_wesely( &
               1, n_species, wesely_config, &
               f_bxheight, f_cldfrc, f_frlai, f_frlanduse, &
               f_iland, f_is_ice, f_is_land, f_is_snow, &
               f_lat, f_lon, local_lucname, f_obk, f_ps, f_salinity, &
               f_suncosmid, f_swgdn, f_ts, f_tskin, &
               real(dt, fp), f_ustar, f_z0, &
               f_mw_g, f_dd_f0, dummy_sp_names, f_dd_hstar, &
               f_dd_DvzAerSnow, f_dd_DvzMinVal_snow, f_dd_DvzMinVal_land, &
               f_conc, col_tendencies, f_is_gas_arr, col_diag_con, col_diag_vel, &
               diagnostic_species_id)
         endif

         ! Execute AEROSOL schemes
         if (trim(local_aero) == "gocart") then
            call compute_gocart( &
               1, n_species, gocart_config, &
               f_airden, f_frlake, f_gwettop, f_hflux, &
               f_lwi, f_pblh, f_t_air, real(dt, fp), &
               f_u10m, f_ustar, f_v10m, f_z_edges, f_z0h, &
               f_density, f_radius, f_is_dust, f_is_seasalt, &
               f_conc, col_tendencies, f_is_gas_arr, col_diag_con, col_diag_vel, &
               diagnostic_species_id)
         else if (trim(local_aero) == "zhang") then
            call compute_zhang( &
               1, n_species, zhang_config, &
               f_bxheight, f_frlanduse, f_iland, &
               f_is_ice, f_is_snow, local_lucname, f_obk, f_ps, f_rh, &
               f_ts, real(dt, fp), f_u10m, f_ustar, f_v10m, f_z0, &
               f_mw_g, f_radius, f_density, dummy_sp_names, &
               f_dd_hstar, f_dd_DvzAerSnow, f_dd_DvzMinVal_snow, &
               f_dd_DvzMinVal_land, f_lower_radius, f_upper_radius, &
               f_is_dust, f_is_seasalt, f_conc, col_tendencies, &
               f_is_gas_arr, col_diag_con, col_diag_vel, diagnostic_species_id)
         endif

         ! Write tendencies and concentrations back in-place (casting back to c_double)
         tendency(icol, 1, :) = real(col_tendencies(1, :), c_double)
         conc(icol, 1, :) = conc(icol, 1, :) + real(dt * col_tendencies(1, :), c_double)

         if (diagnostics /= 0) then
            diag_con(icol, :) = real(col_diag_con, c_double)
            diag_vel(icol, :) = real(col_diag_vel, c_double)
         endif
      end do

   end subroutine run_drydep_science_bridge

end module DryDepScienceBridge_Mod
