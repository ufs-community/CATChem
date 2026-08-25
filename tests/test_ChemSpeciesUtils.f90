!> \file test_ChemSpeciesUtils.f90
!! \brief Test program for ChemSpeciesUtils module
!!
program test_ChemSpeciesUtils
   use testing_mod, only: assert
   use ChemSpeciesUtils_Mod, only: create_species_mapping
   use StateManager_Mod, only: StateManagerType
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE

   implicit none

   type(StateManagerType) :: state_mgr
   character(len=32) :: species_names(2)
   integer :: mapping(2)
   integer :: rc

   write(*,*) 'Testing ChemSpeciesUtils module...'
   write(*,*) ''

   ! Test 1: Unassociated StateManager mapping returns CC_FAILURE
   write(*,*) 'Test 1: Unassociated StateManager mapping'
   species_names(1) = 'O3'
   species_names(2) = 'NO2'
   call create_species_mapping(state_mgr, species_names, mapping, rc)
   call assert(rc == CC_FAILURE, "Unassociated StateManager mapping should return CC_FAILURE")
   call assert(mapping(1) == -1, "Unassociated mapping(1) should be -1")
   call assert(mapping(2) == -1, "Unassociated mapping(2) should be -1")

   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   write(*,*) 'All ChemSpeciesUtils tests passed!'

end program test_ChemSpeciesUtils
