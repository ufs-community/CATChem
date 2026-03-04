!> \file SO4chemProcessCreator_Mod.F90
!! \brief Factory for creating so4chem process instances
!!
!! This module provides the factory functions for creating so4chem
!! process instances following the CATChem Process Factory pattern.
!!
!! Generated on: 2026-03-03T18:15:44.325744
!! Author: Wei Li
!! Version: 1.0.0

module SO4chemProcessCreator_Mod

   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, CC_Error, CC_Warning, ErrorManagerType
   use ProcessInterface_Mod
   use ProcessSO4chemInterface_Mod

   implicit none
   private

   public :: create_so4chem_process
   public :: register_so4chem_process
   public :: get_so4chem_default_config

contains

   !> Create a new so4chem process instance
   !!
   !! This factory function creates and returns a new instance of the
   !! so4chem process. The process is not initialized - the caller
   !! must call the init() method with appropriate configuration.
   !!
   !! @param[out] process     Allocated process instance
   !! @param[out] rc          Return code
   subroutine create_so4chem_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessSO4chemInterface), allocatable :: so4chem_process
      integer :: alloc_stat

      rc = CC_SUCCESS

      ! Allocate the process instance
      allocate(so4chem_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = CC_FAILURE
         return
      end if

      ! Move to polymorphic variable
      call move_alloc(so4chem_process, process)

   end subroutine create_so4chem_process

   !> Register the so4chem process with a ProcessManager
   !!
   !! This subroutine registers the so4chem process with a ProcessManager's
   !! factory. This is the correct way to register processes for use in
   !! applications and integration tests.
   !!
   !! @param[inout] process_mgr The ProcessManager to register with
   !! @param[out] rc Return code
   subroutine register_so4chem_process(process_mgr, rc)
      use ProcessManager_Mod, only: ProcessManagerType

      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      call process_mgr%register_process( &
         name='so4chem', &
         category='chemistry', &
         description='Process for computing chemical production of sulfate from SO2 oxidation', &
         creator=create_so4chem_process, &
         rc=rc &
         )

   end subroutine register_so4chem_process

   !> Get default configuration for so4chem process
   !!
   !! This function returns a default configuration string that can be
   !! used to initialize the so4chem process with reasonable defaults.
   !!
   !! @param[out] config_data Default configuration string
   subroutine get_so4chem_default_config(config_data)
      character(len=*), intent(out) :: config_data

      ! Return default YAML configuration
      config_data = &
         '# Default so4chem process configuration' // new_line('A') // &
         'process:' // new_line('A') // &
         '  name: "so4chem"' // new_line('A') // &
         '  version: "1.0.0"' // new_line('A') // &
         '  active_scheme: ""' // new_line('A') // &
         '  is_active: true' // new_line('A') // &
         '' // new_line('A') // &
         '# Scheme configuration' // new_line('A') // &
         'schemes:' // new_line('A') // &
         '  gocart:' // new_line('A') // &
         '    description: "GOCART SO2 to SO4 production scheme"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      scale_factor: 1.0' // new_line('A') // &
         '' // new_line('A') // &
         '# Diagnostic configuration' // new_line('A') // &
         'diagnostics:' // new_line('A') // &
         '  output_frequency: 3600.0  # seconds' // new_line('A') // &
         '  output_diagnostics: true'

   end subroutine get_so4chem_default_config

end module SO4chemProcessCreator_Mod
