module SO4chemScienceBridge_Mod
   use iso_c_binding, only: c_ptr, c_f_pointer, c_double, c_char, c_associated, c_bool, c_int
   use Precision_Mod, only: fp
   use Constants, only: g0, Cpd, AVO, VON_KARMAN, AIRMW, PI
   use SO4chemCommon_Mod, only: SO4chemSchemeGOCARTConfig
   use SO4chemScheme_GOCART_Mod, only: compute_gocart
   implicit none
contains

   subroutine run_so4chem_science_bridge( &
      n_cols, n_levels, n_species, dt, &
      diagnostics, &
   ! Date & Time
      year, month, day, hour, minute, second, &
   ! 3D Met Pointers
      c_airden, c_cldf, c_delp, c_pmid, c_t_air, c_z_edges, &
   ! 2D Met Pointers
      c_hflux, c_lat, c_lon, c_lwi, c_pblh, c_u10m, c_ustar, c_v10m, c_z0h, &
   ! Metadata
      species_mw_g, species_names, &
   ! Chem and Tendency
      c_conc, c_tendency, &
   ! Persistent states
      c_firsttime, c_nymd_last, c_nhms_last_recycle, c_xh2o2_init, &
      c_pso4_so2, c_pso4_g_so2, c_pso4_aq_so2, c_pso2_dms, c_dms_flux, &
      diagnostic_species_id, n_diag_species &
      ) bind(C, name="run_so4chem_science_bridge")

      integer(c_int), value :: n_cols, n_levels, n_species
      real(c_double), value :: dt
      integer(c_int), value :: diagnostics

      integer(c_int), value :: year, month, day, hour, minute, second

      ! C pointers
      type(c_ptr), value :: c_airden, c_cldf, c_delp, c_pmid, c_t_air, c_z_edges
      type(c_ptr), value :: c_hflux, c_lat, c_lon, c_lwi, c_pblh, c_u10m, c_ustar, c_v10m, c_z0h
      type(c_ptr), value :: c_conc, c_tendency
      type(c_ptr), value :: c_firsttime, c_nymd_last, c_nhms_last_recycle, c_xh2o2_init
      type(c_ptr), value :: c_pso4_so2, c_pso4_g_so2, c_pso4_aq_so2, c_pso2_dms, c_dms_flux

      ! Metadata
      real(c_double), intent(in) :: species_mw_g(n_species)
      character(kind=c_char), intent(in) :: species_names(32, n_species)

      integer(c_int), value :: n_diag_species
      integer(c_int), intent(in) :: diagnostic_species_id(n_diag_species)

      ! Slicing array pointers pointing directly to double precision (c_double) C++ views
      real(c_double), pointer :: airden(:,:), cldf(:,:), delp(:,:), pmid(:,:), t_air(:,:), z_edges(:,:)
      real(c_double), pointer :: hflux(:), lat(:), lon(:), pblh(:), u10m(:), ustar(:), v10m(:), z0h(:)
      integer, pointer :: lwi(:)
      real(c_double), pointer :: conc(:,:,:), tendency(:,:,:)

      ! Persistent pointers pointing to double precision C++ views
      logical(c_bool), pointer :: firsttime(:)
      integer, pointer :: nymd_last(:), nhms_last_recycle(:)
      real(c_double), pointer :: xh2o2_init(:,:), pso4_so2(:,:), pso4_g_so2(:,:), pso4_aq_so2(:,:), pso2_dms(:,:), dms_flux(:)

      ! Loop variables
      integer :: icol, i, j, ispec
      character(len=32) :: dummy_sp_names(n_species)
      logical :: f_firsttime, is_aero

      ! Local arrays in native solver precision (fp) to avoid double-float mismatches
      real(fp) :: f_airden(n_levels)
      real(fp) :: f_cldf(n_levels)
      real(fp) :: f_delp(n_levels)
      real(fp) :: f_pmid(n_levels)
      real(fp) :: f_t_air(n_levels)
      real(fp) :: f_z_edges(n_levels+1)

      real(fp) :: f_hflux, f_lat, f_lon, f_pblh, f_u10m, f_ustar, f_v10m, f_z0h

      ! Casted metadata properties
      real(fp) :: f_mw_g(n_species)

      ! Sliced concentration and tendencies in solver precision
      real(fp) :: f_conc(n_levels, n_species)
      real(fp) :: col_tendencies(n_levels, n_species)
      real(fp) :: col_conc_new(n_levels)

      ! Local arrays to resolve Fortran BIND(C) allocatable constraints & rank matches
      real(fp), allocatable :: local_xh2o2_init(:)
      real(fp) :: col_prod_rate(n_levels, n_species)
      real(fp) :: col_pso4_g(n_levels)
      real(fp) :: col_pso4_aq(n_levels)
      real(fp) :: col_dms_flux

      type(SO4chemSchemeGOCARTConfig) :: gocart_config

      ! Associate pointers
      call c_f_pointer(c_airden,   airden,   [n_cols, n_levels])
      call c_f_pointer(c_cldf,     cldf,     [n_cols, n_levels])
      call c_f_pointer(c_delp,     delp,     [n_cols, n_levels])
      call c_f_pointer(c_pmid,     pmid,     [n_cols, n_levels])
      call c_f_pointer(c_t_air,    t_air,    [n_cols, n_levels])
      call c_f_pointer(c_z_edges,  z_edges,  [n_cols, n_levels+1])

      call c_f_pointer(c_hflux,    hflux,    [n_cols])
      call c_f_pointer(c_lat,      lat,      [n_cols])
      call c_f_pointer(c_lon,      lon,      [n_cols])
      call c_f_pointer(c_lwi,      lwi,      [n_cols])
      call c_f_pointer(c_pblh,     pblh,     [n_cols])
      call c_f_pointer(c_u10m,     u10m,     [n_cols])
      call c_f_pointer(c_ustar,    ustar,    [n_cols])
      call c_f_pointer(c_v10m,     v10m,     [n_cols])
      call c_f_pointer(c_z0h,      z0h,      [n_cols])

      call c_f_pointer(c_conc,       conc,       [n_cols, n_levels, n_species])
      call c_f_pointer(c_tendency,   tendency,   [n_cols, n_levels, n_species])

      ! Associate persistent pointers
      call c_f_pointer(c_firsttime,         firsttime,         [n_cols])
      call c_f_pointer(c_nymd_last,         nymd_last,         [n_cols])
      call c_f_pointer(c_nhms_last_recycle, nhms_last_recycle, [n_cols])
      call c_f_pointer(c_xh2o2_init,         xh2o2_init,         [n_cols, n_levels])
      call c_f_pointer(c_pso4_so2,          pso4_so2,          [n_cols, n_levels])
      call c_f_pointer(c_pso4_g_so2,        pso4_g_so2,        [n_cols, n_levels])
      call c_f_pointer(c_pso4_aq_so2,       pso4_aq_so2,       [n_cols, n_levels])
      call c_f_pointer(c_pso2_dms,          pso2_dms,          [n_cols, n_levels])
      call c_f_pointer(c_dms_flux,          dms_flux,          [n_cols])

      ! Extract real species names from flat char array passed via BIND(C)
      do i = 1, n_species
         do j = 1, 32
            dummy_sp_names(i)(j:j) = species_names(j, i)
         end do
         dummy_sp_names(i) = trim(adjustl(dummy_sp_names(i)))
      end do

      ! Copy metadata properties once
      f_mw_g = real(species_mw_g, fp)

      ! Iterate columns
      do icol = 1, n_cols
         ! Cast scalars and slices
         f_airden       = real(airden(icol, :), fp)
         f_cldf         = real(cldf(icol, :), fp)
         f_delp         = real(delp(icol, :), fp)
         f_pmid         = real(pmid(icol, :), fp)
         f_t_air        = real(t_air(icol, :), fp)
         f_z_edges      = real(z_edges(icol, :), fp)

         f_hflux        = real(hflux(icol), fp)
         f_lat          = real(lat(icol), fp)
         f_lon          = real(lon(icol), fp)
         f_pblh         = real(pblh(icol), fp)
         f_u10m         = real(u10m(icol), fp)
         f_ustar        = real(ustar(icol), fp)
         f_v10m         = real(v10m(icol), fp)
         f_z0h          = real(z0h(icol), fp)

         f_conc         = real(conc(icol, :, :), fp)
         col_tendencies = 0.0_fp
         col_prod_rate  = 0.0_fp
         col_pso4_g     = 0.0_fp
         col_pso4_aq    = 0.0_fp
         col_dms_flux   = 0.0_fp

         ! Copy to local allocatable H2O2 initialization buffer
         allocate(local_xh2o2_init(n_levels))
         local_xh2o2_init = real(xh2o2_init(icol, :), fp)

         ! Copy logical value
         f_firsttime = firsttime(icol)

         ! Execute GOCART sulfur chemistry solver
         call compute_gocart( &
            n_levels, n_species, gocart_config, &
            g0, Cpd, AVO, VON_KARMAN, AIRMW, PI, &
            year, month, day, hour, minute, second, &
            f_airden, f_cldf, f_delp, &
            f_hflux, f_lat, f_lon, lwi(icol), f_pblh, f_pmid, &
            f_t_air, real(dt, fp), f_u10m, f_ustar, f_v10m, &
            f_z_edges, f_z0h, &
            f_mw_g, dummy_sp_names, &
            f_conc, col_tendencies, &
            f_firsttime, nymd_last(icol), nhms_last_recycle(icol), local_xh2o2_init, &
            col_prod_rate, col_pso4_g, col_pso4_aq, col_dms_flux, &
            diagnostic_species_id=diagnostic_species_id)

         ! Convert updated GOCART concentrations (ppm or ug/kg) to actual tendencies in kg/kg/s
         do ispec = 1, n_species
            is_aero = .false.
            if ( dummy_sp_names(ispec) == 'SO4' .or. dummy_sp_names(ispec) == 'so4' .or. &
               dummy_sp_names(ispec) == 'MSA' .or. dummy_sp_names(ispec) == 'msa' .or. &
               dummy_sp_names(ispec) == 'ASO4J' .or. dummy_sp_names(ispec) == 'aso4j' ) then
               is_aero = .true.
            end if

            if ( is_aero ) then
               col_conc_new(:) = col_tendencies(:, ispec) * 1.0e-9_fp
            else
               col_conc_new(:) = col_tendencies(:, ispec) * 1.0e-6_fp * (f_mw_g(ispec) / AIRMW)
            end if

            col_tendencies(:, ispec) = (col_conc_new(:) - real(conc(icol, :, ispec), fp)) / dt
         end do

         ! Write tendencies and updated concentrations back in-place (casting to c_double)
         tendency(icol, :, :) = real(col_tendencies, c_double)
         conc(icol, :, :) = conc(icol, :, :) + real(dt * col_tendencies, c_double)

         ! Copy persistent changes and diagnostics back to C++ buffers (casting to c_double)
         firsttime(icol)     = f_firsttime
         xh2o2_init(icol, :) = real(local_xh2o2_init, c_double)
         deallocate(local_xh2o2_init)

         pso4_so2(icol, :)    = real(col_prod_rate(:, 1), c_double) ! Maps Production_rate for first diagnostic species (e.g. SO2)
         pso4_g_so2(icol, :)  = real(col_pso4_g, c_double)
         pso4_aq_so2(icol, :) = real(col_pso4_aq, c_double)
         dms_flux(icol)       = real(col_dms_flux, c_double)
      end do

   end subroutine run_so4chem_science_bridge

end module SO4chemScienceBridge_Mod
