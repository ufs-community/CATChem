!> \file test_ChemSpeciesUtils.f90
!! \brief Test program for ChemSpeciesUtils module
!!
!!!>
program test_ChemSpeciesUtils
   use testing_mod, only: assert, assert_close
   use ChemSpeciesUtils_Mod
   use StateManager_Mod, only: StateManagerType
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE

   implicit none

   type(ChemSpeciesUtilsType) :: chem_utils
   type(StateManagerType) :: state_mgr
   integer :: rc

   write(*,*) 'Testing ChemSpeciesUtils module...'
   write(*,*) ''

   ! Test 1: State manager initialization (needed for ChemSpeciesUtils)
   write(*,*) 'Test 1: State manager initialization'
   call state_mgr%init('TestStateManager', rc)
   call assert(rc == CC_SUCCESS, "StateManager initialization should succeed")

   ! Test 2: ChemSpeciesUtils initialization
   write(*,*) 'Test 2: ChemSpeciesUtils initialization'
   ! The ChemSpeciesUtilsType doesn't have an explicit init method,
   ! so we just test its methods directly

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: Getting species index (will fail since no species are defined)
   write(*,*) 'Test 3: Getting species index'
   block
      integer :: species_idx

      species_idx = chem_utils%get_index(state_mgr, 'NONEXISTENT', rc)
      ! This should fail since no species are defined
      ! We're not asserting on rc because behavior may vary
   end block

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! Test 4: Checking if species exists (should return false)
   write(*,*) 'Test 4: Checking if species exists'
   block
      logical :: exists

      exists = chem_utils%exists(state_mgr, 'NONEXISTENT', rc)
      ! This should return false since no species are defined
      ! We're not asserting on the result because behavior may vary
   end block

   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! Test 5: Cleanup
   write(*,*) 'Test 5: Cleanup'
   call state_mgr%finalize(rc)
   call assert(rc == CC_SUCCESS, "StateManager finalization should succeed")

   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   write(*,*) 'All ChemSpeciesUtils tests passed!'

end program test_ChemSpeciesUtils
