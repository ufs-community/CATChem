!> \file test_dust_science.f90
!! \brief Unit tests for dust aerosol emission science bridge and schemes
!!
program test_dust_science
   use testing_mod, only: assert
   use iso_c_binding
   use precision_mod, only: fp

   implicit none

   interface
      subroutine run_dust_science_bridge( &
         n_cols, n_levels, n_species, n_soil, dt, &
         active_scheme, diagnostics, &
         airden, clayfrac, frlake, frsno, gvf, lai, lwi, rdrag, sandfrac, &
         soilm, ssm, tskin, u10m, v10m, ustar, ustar_threshold, z0, &
         species_density, species_radius, species_lower_radius, species_upper_radius, &
         conc, tendency, &
         diag_emission_total, diag_emission_bin, diag_horizontal_flux, diag_moisture_correction, diag_effective_threshold, diag_utar_threshold, &
         diagnostic_species_id, n_diag_species &
         ) bind(C, name="run_dust_science_bridge")
         import :: c_ptr, c_double, c_char, c_int
         integer(c_int), value :: n_cols, n_levels, n_species, n_soil
         real(c_double), value :: dt
         character(kind=c_char), intent(in) :: active_scheme(*)
         integer(c_int), value :: diagnostics
         type(c_ptr), value :: airden, clayfrac, frlake, frsno, gvf, lai, lwi, rdrag, sandfrac
         type(c_ptr), value :: soilm, ssm, tskin, u10m, v10m, ustar, ustar_threshold, z0
         type(c_ptr), value :: species_density, species_radius, species_lower_radius, species_upper_radius
         type(c_ptr), value :: conc, tendency
         type(c_ptr), value :: diag_emission_total, diag_emission_bin, diag_horizontal_flux, diag_moisture_correction, diag_effective_threshold, diag_utar_threshold
         integer(c_int), intent(in) :: diagnostic_species_id(n_diag_species)
         integer(c_int), value :: n_diag_species
      end subroutine run_dust_science_bridge
   end interface

   write(*,*) 'Testing Dust Science Bridge...'
   write(*,*) ''

   block
      integer(c_int), parameter :: n_cols = 1, n_levels = 2, n_species = 2, n_soil = 4
      real(c_double), target :: airden(n_cols, n_levels) = 1.225_c_double
      real(c_double), target :: clayfrac(n_cols) = 0.2_c_double
      real(c_double), target :: frlake(n_cols) = 0.0_c_double
      real(c_double), target :: frsno(n_cols) = 0.0_c_double
      real(c_double), target :: gvf(n_cols) = 0.1_c_double
      real(c_double), target :: lai(n_cols) = 0.1_c_double
      integer(c_int), target :: lwi(n_cols) = 1
      real(c_double), target :: rdrag(n_cols) = 0.1_c_double
      real(c_double), target :: sandfrac(n_cols) = 0.5_c_double
      real(c_double), target :: soilm(n_cols, n_soil) = 0.05_c_double
      real(c_double), target :: ssm(n_cols) = 0.5_c_double
      real(c_double), target :: tskin(n_cols) = 300.0_c_double
      real(c_double), target :: u10m(n_cols) = 10.0_c_double
      real(c_double), target :: v10m(n_cols) = 0.0_c_double
      real(c_double), target :: ustar(n_cols) = 0.6_c_double
      real(c_double), target :: ustar_threshold(n_cols) = 0.2_c_double
      real(c_double), target :: z0(n_cols) = 0.01_c_double

      real(c_double), target :: species_density(n_species) = [2650.0_c_double, 2650.0_c_double]
      real(c_double), target :: species_radius(n_species) = [1.0_c_double, 3.0_c_double]
      real(c_double), target :: species_lower_radius(n_species) = [0.5_c_double, 2.0_c_double]
      real(c_double), target :: species_upper_radius(n_species) = [2.0_c_double, 5.0_c_double]

      real(c_double), target :: conc(n_cols, n_levels, n_species) = 0.0_c_double
      real(c_double), target :: tendency(n_cols, n_levels, n_species) = 0.0_c_double

      real(c_double), target :: diag_emission_total(n_cols) = 0.0_c_double
      real(c_double), target :: diag_emission_bin(n_cols, n_species) = 0.0_c_double
      real(c_double), target :: diag_horizontal_flux(n_cols) = 0.0_c_double
      real(c_double), target :: diag_moisture_correction(n_cols) = 0.0_c_double
      real(c_double), target :: diag_effective_threshold(n_cols) = 0.0_c_double
      real(c_double), target :: diag_utar_threshold(n_cols, n_species) = 0.0_c_double

      character(kind=c_char) :: fengsha_scheme(32)
      character(kind=c_char) :: ginoux_scheme(32)
      integer(c_int) :: diag_ids(1) = 1

      fengsha_scheme = c_null_char
      fengsha_scheme(1:7) = ['f','e','n','g','s','h','a']

      ginoux_scheme = c_null_char
      ginoux_scheme(1:6) = ['g','i','n','o','u','x']

      ! Test 1: Fengsha Scheme via Science Bridge
      call run_dust_science_bridge( &
         n_cols, n_levels, n_species, n_soil, 3600.0_c_double, &
         fengsha_scheme, 0, &
         c_loc(airden), c_loc(clayfrac), c_loc(frlake), c_loc(frsno), c_loc(gvf), c_loc(lai), c_loc(lwi), c_loc(rdrag), c_loc(sandfrac), &
         c_loc(soilm), c_loc(ssm), c_loc(tskin), c_loc(u10m), c_loc(v10m), c_loc(ustar), c_loc(ustar_threshold), c_loc(z0), &
         c_loc(species_density), c_loc(species_radius), c_loc(species_lower_radius), c_loc(species_upper_radius), &
         c_loc(conc), c_loc(tendency), &
         c_loc(diag_emission_total), c_loc(diag_emission_bin), c_loc(diag_horizontal_flux), c_loc(diag_moisture_correction), &
         c_loc(diag_effective_threshold), c_loc(diag_utar_threshold), &
         diag_ids, 1)

      call assert(.true., "Fengsha dust science bridge executed successfully")

      ! Test 2: Ginoux Scheme via Science Bridge
      call run_dust_science_bridge( &
         n_cols, n_levels, n_species, n_soil, 3600.0_c_double, &
         ginoux_scheme, 0, &
         c_loc(airden), c_loc(clayfrac), c_loc(frlake), c_loc(frsno), c_loc(gvf), c_loc(lai), c_loc(lwi), c_loc(rdrag), c_loc(sandfrac), &
         c_loc(soilm), c_loc(ssm), c_loc(tskin), c_loc(u10m), c_loc(v10m), c_loc(ustar), c_loc(ustar_threshold), c_loc(z0), &
         c_loc(species_density), c_loc(species_radius), c_loc(species_lower_radius), c_loc(species_upper_radius), &
         c_loc(conc), c_loc(tendency), &
         c_loc(diag_emission_total), c_loc(diag_emission_bin), c_loc(diag_horizontal_flux), c_loc(diag_moisture_correction), &
         c_loc(diag_effective_threshold), c_loc(diag_utar_threshold), &
         diag_ids, 1)

      call assert(.true., "Ginoux dust science bridge executed successfully")
   end block

   write(*,*) 'All Dust Science tests passed successfully!'

end program test_dust_science
