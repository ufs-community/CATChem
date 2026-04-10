!> \file CarbChemProcessCreator_Mod.F90
!! \brief Factory for creating carbchem process instances
!!
!! This module provides the factory functions for creating carbchem
!! process instances following the CATChem Process Factory pattern.
!!
!! Generated on: 2026-04-10T16:46:37.706128
!! Author: Wei Li
!! Version: 1.0.0

module CarbChemProcessCreator_Mod

   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, CC_Error, CC_Warning, ErrorManagerType
   use ProcessInterface_Mod
   use ProcessCarbChemInterface_Mod

   implicit none
   private

   public :: create_carbchem_process
   public :: register_carbchem_process
   public :: get_carbchem_default_config

contains

   !> Create a new carbchem process instance
   !!
   !! This factory function creates and returns a new instance of the
   !! carbchem process. The process is not initialized - the caller
   !! must call the init() method with appropriate configuration.
   !!
   !! @param[out] process     Allocated process instance
   !! @param[out] rc          Return code
   subroutine create_carbchem_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessCarbChemInterface), allocatable :: carbchem_process
      integer :: alloc_stat

      rc = CC_SUCCESS

      ! Allocate the process instance
      allocate(carbchem_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = CC_FAILURE
         return
      end if

      ! Move to polymorphic variable
      call move_alloc(carbchem_process, process)

   end subroutine create_carbchem_process

   !> Register the carbchem process with a ProcessManager
   !!
   !! This subroutine registers the carbchem process with a ProcessManager's
   !! factory. This is the correct way to register processes for use in
   !! applications and integration tests.
   !!
   !! @param[inout] process_mgr The ProcessManager to register with
   !! @param[out] rc Return code
   subroutine register_carbchem_process(process_mgr, rc)
      use ProcessManager_Mod, only: ProcessManagerType

      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      call process_mgr%register_process( &
         name='carbchem', &
         category='chemistry', &
         description='Process for computing chemical production and loss of carbon species', &
         creator=create_carbchem_process, &
         rc=rc &
         )

   end subroutine register_carbchem_process

   !> Get default configuration for carbchem process
   !!
   !! This function returns a default configuration string that can be
   !! used to initialize the carbchem process with reasonable defaults.
   !!
   !! @param[out] config_data Default configuration string
   subroutine get_carbchem_default_config(config_data)
      character(len=*), intent(out) :: config_data

      ! Return default YAML configuration
      config_data = &
         '# Default carbchem process configuration' // new_line('A') // &
         'process:' // new_line('A') // &
         '  name: "carbchem"' // new_line('A') // &
         '  version: "1.0.0"' // new_line('A') // &
         '  active_scheme: ""' // new_line('A') // &
         '  is_active: true' // new_line('A') // &
         '' // new_line('A') // &
         '# Scheme configuration' // new_line('A') // &
         'schemes:' // new_line('A') // &
         '  gocart:' // new_line('A') // &
         '    description: "GOCART carbon species chemical production and loss scheme"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      scale_factor: 1.0' // new_line('A') // &
         '' // new_line('A') // &
         '# Diagnostic configuration' // new_line('A') // &
         'diagnostics:' // new_line('A') // &
         '  output_frequency: 3600.0  # seconds' // new_line('A') // &
         '  output_diagnostics: true'

   end subroutine get_carbchem_default_config

end module CarbChemProcessCreator_Mod
