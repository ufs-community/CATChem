!> \file test_wetdep_science.f90
!! \brief Unit tests for Jacob wet deposition scheme
!!
program test_wetdep_science
   use testing_mod, only: assert
   use precision_mod, only: fp
   use WetDepCommon_Mod, only: WetDepSchemeJACOBConfig
   use WetDepScheme_JACOB_Mod, only: compute_jacob

   implicit none

   write(*,*) 'Testing WetDep Science Schemes...'
   write(*,*) ''

   ! Test 1: Jacob Wet Deposition Scheme
   block
      integer, parameter :: num_layers = 2
      integer, parameter :: num_species = 3
      type(WetDepSchemeJACOBConfig) :: config
      real(fp) :: dt
      real(fp) :: airden_dry(num_layers)
      real(fp) :: mairden(num_layers)
      real(fp) :: pedge(num_layers + 1)
      real(fp) :: pfilsan(num_layers + 1)
      real(fp) :: pfllsan(num_layers + 1)
      real(fp) :: reevapls(num_layers)
      real(fp) :: t(num_layers)

      logical :: species_is_aerosol(num_species)
      character(len=32) :: species_short_name(num_species)
      real(fp) :: species_henry_cr(num_species)
      real(fp) :: species_henry_k0(num_species)
      real(fp) :: species_henry_pKa(num_species)
      real(fp) :: species_wd_retfactor(num_species)
      logical :: species_wd_LiqAndGas(num_species)
      real(fp) :: species_wd_convfacI2G(num_species)
      real(fp) :: species_wd_rainouteff(num_species, 3)
      real(fp) :: species_wd_reevap_frac(num_species)
      real(fp) :: species_radius(num_species)
      real(fp) :: species_mw_g(num_species)

      real(fp) :: species_conc(num_layers, num_species)
      real(fp) :: species_tendencies(num_layers, num_species)

      dt = 3600.0_fp
      airden_dry = 1.2_fp
      mairden = 1.2_fp
      pedge = [101325.0_fp, 90000.0_fp, 80000.0_fp]
      pfilsan = 0.0_fp
      pfllsan = 0.0_fp
      reevapls = 0.0_fp
      t = [280.0_fp, 275.0_fp]

      species_is_aerosol = [.false., .true., .false.]
      species_short_name = ['SO2 ', 'SO4 ', 'H2O2']
      species_henry_cr = 0.0_fp
      species_henry_k0 = 1.0e2_fp
      species_henry_pKa = 0.0_fp
      species_wd_retfactor = 1.0_fp
      species_wd_LiqAndGas = .false.
      species_wd_convfacI2G = 0.0_fp
      species_wd_rainouteff = 1.0_fp
      species_wd_reevap_frac = 0.5_fp
      species_radius = 1.0e-6_fp
      species_mw_g = [64.0_fp, 96.0_fp, 34.0_fp]

      species_conc = 1.0e-9_fp
      species_tendencies = 0.0_fp

      call compute_jacob( &
         num_layers, num_species, config, &
         airden_dry, mairden, pedge, pfilsan, pfllsan, &
         reevapls, t, dt, &
         species_is_aerosol, species_short_name, &
         species_henry_cr, species_henry_k0, species_henry_pKa, &
         species_wd_retfactor, species_wd_LiqAndGas, &
         species_wd_convfacI2G, species_wd_rainouteff, &
         species_wd_reevap_frac, species_radius, species_mw_g, &
         species_conc, species_tendencies)

      call assert(.true., "Jacob wet deposition compute executed successfully")
   end block

   write(*,*) 'All WetDep Science tests passed successfully!'

end program test_wetdep_science
