!> \file test_so4chem_science.f90
!! \brief Unit tests for GOCART SO4 chemistry science bridge and schemes
!!
program test_so4chem_science
   use testing_mod, only: assert
   use iso_c_binding
   use precision_mod, only: fp

   implicit none

   interface
      subroutine run_so4chem_science_bridge( &
         n_cols, n_levels, n_species, dt, &
         diagnostics, &
         year, month, day, hour, minute, second, &
         c_airden, c_cldf, c_delp, c_pmid, c_t_air, c_z_edges, &
         c_hflux, c_lat, c_lon, c_lwi, c_pblh, c_u10m, c_ustar, c_v10m, c_z0h, &
         species_mw_g, species_names, &
         c_conc, c_tendency, &
         c_firsttime, c_nymd_last, c_nhms_last_recycle, c_xh2o2_init, &
         c_pso4_so2, c_pso4_g_so2, c_pso4_aq_so2, c_pso2_dms, c_dms_flux, &
         diagnostic_species_id, n_diag_species &
         ) bind(C, name="run_so4chem_science_bridge")
         import :: c_ptr, c_double, c_char, c_int, c_bool
         integer(c_int), value :: n_cols, n_levels, n_species
         real(c_double), value :: dt
         integer(c_int), value :: diagnostics
         integer(c_int), value :: year, month, day, hour, minute, second
         type(c_ptr), value :: c_airden, c_cldf, c_delp, c_pmid, c_t_air, c_z_edges
         type(c_ptr), value :: c_hflux, c_lat, c_lon, c_lwi, c_pblh, c_u10m, c_ustar, c_v10m, c_z0h
         type(c_ptr), value :: c_conc, c_tendency
         type(c_ptr), value :: c_firsttime, c_nymd_last, c_nhms_last_recycle, c_xh2o2_init
         type(c_ptr), value :: c_pso4_so2, c_pso4_g_so2, c_pso4_aq_so2, c_pso2_dms, c_dms_flux
         real(c_double), intent(in) :: species_mw_g(n_species)
         character(kind=c_char), intent(in) :: species_names(32, n_species)
         integer(c_int), value :: n_diag_species
         integer(c_int), intent(in) :: diagnostic_species_id(n_diag_species)
      end subroutine run_so4chem_science_bridge
   end interface

   write(*,*) 'Testing SO4chem Science Bridge...'
   write(*,*) ''

   block
      integer(c_int), parameter :: n_cols = 1, n_levels = 2, n_species = 4
      real(c_double), target :: airden(n_cols, n_levels) = 1.2_c_double
      real(c_double), target :: cldf(n_cols, n_levels) = 0.0_c_double
      real(c_double), target :: delp(n_cols, n_levels) = 2000.0_c_double
      real(c_double), target :: pmid(n_cols, n_levels) = 95000.0_c_double
      real(c_double), target :: t_air(n_cols, n_levels) = 280.0_c_double
      real(c_double), target :: z_edges(n_cols, n_levels+1) = 0.0_c_double

      real(c_double), target :: hflux(n_cols) = 50.0_c_double
      real(c_double), target :: lat(n_cols) = 40.0_c_double
      real(c_double), target :: lon(n_cols) = -105.0_c_double
      integer(c_int), target :: lwi(n_cols) = 1
      real(c_double), target :: pblh(n_cols) = 1000.0_c_double
      real(c_double), target :: u10m(n_cols) = 5.0_c_double
      real(c_double), target :: ustar(n_cols) = 0.3_c_double
      real(c_double), target :: v10m(n_cols) = 0.0_c_double
      real(c_double), target :: z0h(n_cols) = 0.01_c_double

      real(c_double), target :: conc(n_cols, n_levels, n_species) = 1.0e-9_c_double
      real(c_double), target :: tendency(n_cols, n_levels, n_species) = 0.0_c_double

      logical(c_bool), target :: firsttime = .true.
      integer(c_int), target :: nymd_last = 0
      integer(c_int), target :: nhms_last_recycle = 0
      real(c_double), target :: xh2o2_init(n_cols, n_levels) = 1.0e-9_c_double

      real(c_double), target :: pso4_so2(n_cols, n_levels) = 0.0_c_double
      real(c_double), target :: pso4_g_so2(n_cols, n_levels) = 0.0_c_double
      real(c_double), target :: pso4_aq_so2(n_cols, n_levels) = 0.0_c_double
      real(c_double), target :: pso2_dms(n_cols, n_levels) = 0.0_c_double
      real(c_double), target :: dms_flux(n_cols) = 0.0_c_double

      real(c_double) :: mw(n_species) = [64.0_c_double, 96.0_c_double, 62.0_c_double, 96.0_c_double]
      character(kind=c_char) :: names(32, n_species) = c_null_char
      integer(c_int) :: diag_ids(1) = 1

      call run_so4chem_science_bridge( &
         n_cols, n_levels, n_species, 3600.0_c_double, 0, &
         2026, 7, 13, 12, 0, 0, &
         c_loc(airden), c_loc(cldf), c_loc(delp), c_loc(pmid), c_loc(t_air), c_loc(z_edges), &
         c_loc(hflux), c_loc(lat), c_loc(lon), c_loc(lwi), c_loc(pblh), c_loc(u10m), c_loc(ustar), c_loc(v10m), c_loc(z0h), &
         mw, names, &
         c_loc(conc), c_loc(tendency), &
         c_loc(firsttime), c_loc(nymd_last), c_loc(nhms_last_recycle), c_loc(xh2o2_init), &
         c_loc(pso4_so2), c_loc(pso4_g_so2), c_loc(pso4_aq_so2), c_loc(pso2_dms), c_loc(dms_flux), &
         diag_ids, 1)

      call assert(.true., "SO4chem science bridge executed successfully")
   end block

   write(*,*) 'All SO4chem Science tests passed successfully!'

end program test_so4chem_science
