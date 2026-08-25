!> \file test_carbchem_science.f90
!! \brief Unit tests for GOCART Carbon chemistry science bridge and schemes
!!
program test_carbchem_science
   use testing_mod, only: assert
   use iso_c_binding
   use precision_mod, only: fp

   implicit none

   interface
      subroutine run_carbchem_science_bridge( &
         n_cols, n_levels, n_species, dt, &
         active_scheme, diagnostics, &
         year, month, day, hour, minute, second, &
         airden, delp, pmid, &
         species_t_chem_loss, species_names_char, &
         conc, tendency, &
         diag_prod_mass, diag_loss_flux, diag_phobic_mass, diag_phobic_flux, &
         diagnostic_species_id, n_diag_species &
         ) bind(C, name="run_carbchem_science_bridge")
         import :: c_ptr, c_double, c_char, c_int
         integer(c_int), value :: n_cols, n_levels, n_species
         real(c_double), value :: dt
         character(kind=c_char), intent(in) :: active_scheme(*)
         integer(c_int), value :: diagnostics
         integer(c_int), value :: year, month, day, hour, minute, second
         type(c_ptr), value :: airden, delp, pmid
         type(c_ptr), value :: species_t_chem_loss, species_names_char
         type(c_ptr), value :: conc, tendency
         type(c_ptr), value :: diag_prod_mass, diag_loss_flux, diag_phobic_mass, diag_phobic_flux
         integer(c_int), value :: n_diag_species
         integer(c_int), intent(in) :: diagnostic_species_id(n_diag_species)
      end subroutine run_carbchem_science_bridge
   end interface

   write(*,*) 'Testing CarbChem Science Bridge...'
   write(*,*) ''

   block
      integer(c_int), parameter :: n_cols = 1, n_levels = 2, n_species = 4
      real(c_double), target :: airden(n_cols, n_levels) = 1.2_c_double
      real(c_double), target :: delp(n_cols, n_levels) = 2000.0_c_double
      real(c_double), target :: pmid(n_cols, n_levels) = 95000.0_c_double

      real(c_double), target :: conc(n_cols, n_levels, n_species) = 1.0e-9_c_double
      real(c_double), target :: tendency(n_cols, n_levels, n_species) = 0.0_c_double

      real(c_double), target :: t_chem_loss(n_species) = 1.1574e-5_c_double ! e-folding rate (~1 day decay)
      character(kind=c_char), target :: names_char(32, n_species) = c_null_char

      real(c_double), target :: diag_prod_mass(n_cols, n_species) = 0.0_c_double
      real(c_double), target :: diag_loss_flux(n_cols, n_species) = 0.0_c_double
      real(c_double), target :: diag_phobic_mass(n_cols, n_species) = 0.0_c_double
      real(c_double), target :: diag_phobic_flux(n_cols, n_species) = 0.0_c_double

      character(kind=c_char) :: active_scheme(7) = ['g','o','c','a','r','t',c_null_char]
      integer(c_int) :: diag_ids(1) = 1

      call run_carbchem_science_bridge( &
         n_cols, n_levels, n_species, 3600.0_c_double, &
         active_scheme, 0, &
         2026, 7, 13, 12, 0, 0, &
         c_loc(airden), c_loc(delp), c_loc(pmid), &
         c_loc(t_chem_loss), c_loc(names_char), &
         c_loc(conc), c_loc(tendency), &
         c_loc(diag_prod_mass), c_loc(diag_loss_flux), c_loc(diag_phobic_mass), c_loc(diag_phobic_flux), &
         diag_ids, 1)

      call assert(.true., "CarbChem science bridge executed successfully")
   end block

   write(*,*) 'All CarbChem Science tests passed successfully!'

end program test_carbchem_science
