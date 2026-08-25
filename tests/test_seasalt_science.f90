!> \file test_seasalt_science.f90
!! \brief Unit tests for sea salt aerosol emission schemes (Gong 97, Gong 03, GEOS12)
!!
program test_seasalt_science
   use testing_mod, only: assert, assert_close
   use precision_mod, only: fp
   use SeaSaltCommon_Mod, only: SeaSaltSchemeGONG97Config, SeaSaltSchemeGONG03Config, SeaSaltSchemeGEOS12Config
   use SeaSaltScheme_GONG97_Mod, only: compute_gong97
   use SeaSaltScheme_GONG03_Mod, only: compute_gong03
   use SeaSaltScheme_GEOS12_Mod, only: compute_geos12

   implicit none

   write(*,*) 'Testing SeaSalt Science Schemes...'
   write(*,*) ''

   ! Test 1: Gong 1997 scheme
   write(*,*) 'Test 1: Gong 1997 Sea Salt Emission Scheme'
   block
      integer, parameter :: num_species = 2
      type(SeaSaltSchemeGONG97Config) :: config
      real(fp) :: frocean, frseaice, lat, lon, sst, u10m, v10m
      real(fp) :: species_density(num_species)
      real(fp) :: species_radius(num_species)
      real(fp) :: species_lower_radius(num_species)
      real(fp) :: species_upper_radius(num_species)
      real(fp) :: species_mw_g(num_species)
      real(fp) :: emission_flux(1, num_species)

      frocean = 1.0_fp   ! 100% ocean
      frseaice = 0.0_fp  ! No sea ice
      lat = 20.0_fp
      lon = -150.0_fp
      sst = 295.0_fp
      u10m = 8.0_fp      ! 8 m/s wind
      v10m = 0.0_fp

      species_density = 2200.0_fp
      species_radius = [0.5e-6_fp, 2.0e-6_fp]
      species_lower_radius = [0.1e-6_fp, 1.0e-6_fp]
      species_upper_radius = [1.0e-6_fp, 5.0e-6_fp]
      species_mw_g = 58.44_fp

      call compute_gong97(1, num_species, config, 3.14159265_fp, &
         frocean, frseaice, lat, lon, sst, u10m, v10m, &
         species_density, species_radius, species_lower_radius, species_upper_radius, &
         species_mw_g, emission_flux)

      call assert(emission_flux(1, 1) >= 0.0_fp, "Gong97 emission flux species 1 must be non-negative")
      call assert(emission_flux(1, 2) >= 0.0_fp, "Gong97 emission flux species 2 must be non-negative")
   end block
   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! Test 2: Gong 2003 scheme
   write(*,*) 'Test 2: Gong 2003 Sea Salt Emission Scheme'
   block
      integer, parameter :: num_species = 1
      type(SeaSaltSchemeGONG03Config) :: config
      real(fp) :: frocean, frseaice, lat, lon, sst, u10m, v10m
      real(fp) :: species_density(1), species_radius(1), species_lower_radius(1), species_upper_radius(1), species_mw_g(1)
      real(fp) :: emission_flux(1, 1)

      frocean = 1.0_fp
      frseaice = 0.0_fp
      u10m = 10.0_fp
      v10m = 0.0_fp
      species_density = 2200.0_fp
      species_radius = 1.0e-6_fp
      species_lower_radius = 0.5e-6_fp
      species_upper_radius = 1.5e-6_fp
      species_mw_g = 58.44_fp

      call compute_gong03(1, 1, config, 3.14159265_fp, &
         frocean, frseaice, lat, lon, sst, u10m, v10m, &
         species_density, species_radius, species_lower_radius, species_upper_radius, &
         species_mw_g, emission_flux)

      call assert(emission_flux(1, 1) >= 0.0_fp, "Gong03 emission flux must be non-negative")
   end block
   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   write(*,*) 'All SeaSalt Science tests passed successfully!'

end program test_seasalt_science
