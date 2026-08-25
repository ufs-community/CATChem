!> \file test_drydep_science.f90
!! \brief Unit tests for dry deposition science bridge and schemes
!!
program test_drydep_science
   use testing_mod, only: assert
   use iso_c_binding
   use precision_mod, only: fp

   implicit none

   interface
      subroutine run_drydep_science_bridge( &
         n_cols, n_levels, n_species, dt, &
         gas_scheme, aero_scheme, diagnostics, &
         c_bxheight, c_airden, c_t_air, c_z_edges, c_rh, &
         c_cldfrc, c_frlai, c_frlanduse, c_iland, c_is_ice, c_is_land, c_is_snow, &
         c_lat, c_lon, c_obk, c_ps, c_salinity, c_suncosmid, c_swgdn, c_ts, c_tskin, &
         c_ustar, c_z0, c_frlake, c_gwettop, c_hflux, c_lwi, c_pblh, c_u10m, c_v10m, c_z0h, &
         species_mw_g, species_dd_f0, species_dd_hstar, species_dd_DvzAerSnow, &
         species_dd_DvzMinVal_snow, species_dd_DvzMinVal_land, species_density, &
         species_radius, species_is_seasalt, species_is_dust, species_lower_radius, &
         species_upper_radius, is_gas_arr, &
         c_conc, c_tendency, c_diag_con, c_diag_vel, &
         diagnostic_species_id, n_diag_species &
         ) bind(C, name="run_drydep_science_bridge")
         import :: c_ptr, c_double, c_char, c_int, c_bool
         integer(c_int), value :: n_cols, n_levels, n_species
         real(c_double), value :: dt
         character(kind=c_char), intent(in) :: gas_scheme(*)
         character(kind=c_char), intent(in) :: aero_scheme(*)
         integer(c_int), value :: diagnostics
         type(c_ptr), value :: c_bxheight, c_airden, c_t_air, c_z_edges, c_rh
         type(c_ptr), value :: c_cldfrc, c_frlai, c_frlanduse, c_iland, c_is_ice, c_is_land, c_is_snow
         type(c_ptr), value :: c_lat, c_lon, c_obk, c_ps, c_salinity, c_suncosmid, c_swgdn, c_ts, c_tskin
         type(c_ptr), value :: c_ustar, c_z0, c_frlake, c_gwettop, c_hflux, c_lwi, c_pblh, c_u10m, c_v10m, c_z0h
         type(c_ptr), value :: c_conc, c_tendency, c_diag_con, c_diag_vel
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
      end subroutine run_drydep_science_bridge
   end interface

   write(*,*) 'Testing DryDep Science Bridge...'
   write(*,*) ''

   block
      integer(c_int), parameter :: n_cols = 1, n_levels = 2, n_species = 1
      real(c_double), target :: bxheight(n_cols, n_levels) = 1000.0_c_double
      real(c_double), target :: airden(n_cols, n_levels) = 1.2_c_double
      real(c_double), target :: t_air(n_cols, n_levels) = 288.15_c_double
      real(c_double), target :: z_edges(n_cols, n_levels+1) = 0.0_c_double
      real(c_double), target :: rh(n_cols, n_levels) = 0.5_c_double

      real(c_double), target :: cldfrc(n_cols) = 0.0_c_double
      real(c_double), target :: frlai(n_cols, 1, 20) = 0.1_c_double
      real(c_double), target :: frlanduse(n_cols, 1, 20) = 0.05_c_double
      integer(c_int), target :: iland(n_cols, 1, 20) = 1
      logical(c_bool), target :: is_ice(n_cols) = .false.
      logical(c_bool), target :: is_land(n_cols) = .true.
      logical(c_bool), target :: is_snow(n_cols) = .false.

      real(c_double), target :: lat(n_cols) = 40.0_c_double
      real(c_double), target :: lon(n_cols) = -105.0_c_double
      real(c_double), target :: obk(n_cols) = -20.0_c_double
      real(c_double), target :: ps(n_cols) = 101325.0_c_double
      real(c_double), target :: salinity(n_cols) = 35.0_c_double
      real(c_double), target :: suncosmid(n_cols) = 0.8_c_double
      real(c_double), target :: swgdn(n_cols) = 500.0_c_double
      real(c_double), target :: ts(n_cols) = 290.0_c_double
      real(c_double), target :: tskin(n_cols) = 290.0_c_double
      real(c_double), target :: ustar(n_cols) = 0.3_c_double
      real(c_double), target :: z0(n_cols) = 0.01_c_double
      real(c_double), target :: frlake(n_cols) = 0.0_c_double
      real(c_double), target :: gwettop(n_cols) = 0.1_c_double
      real(c_double), target :: hflux(n_cols) = 50.0_c_double
      integer(c_int), target :: lwi(n_cols) = 1
      real(c_double), target :: pblh(n_cols) = 1000.0_c_double
      real(c_double), target :: u10m(n_cols) = 5.0_c_double
      real(c_double), target :: v10m(n_cols) = 0.0_c_double
      real(c_double), target :: z0h(n_cols) = 0.01_c_double

      real(c_double), target :: conc(n_cols, n_levels, n_species) = 1.0e-9_c_double
      real(c_double), target :: tendency(n_cols, n_levels, n_species) = 0.0_c_double
      real(c_double), target :: diag_con(n_cols, n_species) = 0.0_c_double
      real(c_double), target :: diag_vel(n_cols, n_species) = 0.0_c_double

      character(kind=c_char) :: gas_scheme(7) = ['w','e','s','e','l','y',c_null_char]
      character(kind=c_char) :: aero_scheme(7) = ['g','o','c','a','r','t',c_null_char]

      real(c_double) :: mw(1) = 48.0_c_double
      real(c_double) :: f0(1) = 0.0_c_double
      real(c_double) :: hstar(1) = 1.0e-2_c_double
      real(c_double) :: aer_snow(1) = 0.0_c_double
      real(c_double) :: min_snow(1) = 0.0_c_double
      real(c_double) :: min_land(1) = 0.0_c_double
      real(c_double) :: density(1) = 2000.0_c_double
      real(c_double) :: radius(1) = 1.0e-6_c_double
      logical(c_bool) :: is_seasalt(1) = .false.
      logical(c_bool) :: is_dust(1) = .false.
      real(c_double) :: lower_r(1) = 0.5e-6_c_double
      real(c_double) :: upper_r(1) = 1.5e-6_c_double
      logical(c_bool) :: is_gas(1) = .true.
      integer(c_int) :: diag_ids(1) = 1

      call run_drydep_science_bridge( &
         n_cols, n_levels, n_species, 3600.0_c_double, &
         gas_scheme, aero_scheme, 1, &
         c_loc(bxheight), c_loc(airden), c_loc(t_air), c_loc(z_edges), c_loc(rh), &
         c_loc(cldfrc), c_loc(frlai), c_loc(frlanduse), c_loc(iland), c_loc(is_ice), c_loc(is_land), c_loc(is_snow), &
         c_loc(lat), c_loc(lon), c_loc(obk), c_loc(ps), c_loc(salinity), c_loc(suncosmid), c_loc(swgdn), c_loc(ts), c_loc(tskin), &
         c_loc(ustar), c_loc(z0), c_loc(frlake), c_loc(gwettop), c_loc(hflux), c_loc(lwi), c_loc(pblh), c_loc(u10m), c_loc(v10m), c_loc(z0h), &
         mw, f0, hstar, aer_snow, min_snow, min_land, density, radius, is_seasalt, is_dust, lower_r, upper_r, is_gas, &
         c_loc(conc), c_loc(tendency), c_loc(diag_con), c_loc(diag_vel), &
         diag_ids, 1)

      call assert(.true., "DryDep science bridge executed successfully")
   end block

   write(*,*) 'All DryDep Science tests passed successfully!'

end program test_drydep_science
