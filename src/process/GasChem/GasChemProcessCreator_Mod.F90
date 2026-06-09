!> \file GasChemProcessCreator_Mod.F90
!! \brief Factory for creating GasChem process instances
!!
!! This module provides the factory functions for creating GasChem
!! process instances following the CATChem Process Factory pattern.
!!
!! Generated on: 2026-06-05T09:03:17.458211
!! Author: Maggie Bruckner
!! Version: 1.0.0

module GasChemProcessCreator_Mod

   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, CC_Error, CC_Warning, ErrorManagerType
   use ProcessInterface_Mod
   use ProcessGasChemInterface_Mod

   implicit none
   private

   public :: create_GasChem_process
   public :: register_GasChem_process
   public :: get_GasChem_default_config

contains

   !> Create a new GasChem process instance
   !!
   !! This factory function creates and returns a new instance of the
   !! GasChem process. The process is not initialized - the caller
   !! must call the init() method with appropriate configuration.
   !!
   !! @param[out] process     Allocated process instance
   !! @param[out] rc          Return code
   subroutine create_GasChem_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessGasChemInterface), allocatable :: GasChem_process
      integer :: alloc_stat

      rc = CC_SUCCESS

      ! Allocate the process instance
      allocate(GasChem_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = CC_FAILURE
         return
      end if

      ! Move to polymorphic variable
      call move_alloc(GasChem_process, process)

   end subroutine create_GasChem_process

   !> Register the GasChem process with a ProcessManager
   !!
   !! This subroutine registers the GasChem process with a ProcessManager's
   !! factory. This is the correct way to register processes for use in
   !! applications and integration tests.
   !!
   !! @param[inout] process_mgr The ProcessManager to register with
   !! @param[out] rc Return code
   subroutine register_GasChem_process(process_mgr, rc)
      use ProcessManager_Mod, only: ProcessManagerType

      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      call process_mgr%register_process( &
         name='GasChem', &
         category='chemistry', &
         description='Process for MICM gas phase chemical solver', &
         creator=create_GasChem_process, &
         rc=rc &
      )

   end subroutine register_GasChem_process

   !> Get default configuration for GasChem process
   !!
   !! This function returns a default configuration string that can be
   !! used to initialize the GasChem process with reasonable defaults.
   !!
   !! @param[out] config_data Default configuration string
   subroutine get_GasChem_default_config(config_data)
      character(len=*), intent(out) :: config_data

      ! Return default YAML configuration
      config_data = &
         '# Default GasChem process configuration' // new_line('A') // &
         'process:' // new_line('A') // &
         '  name: "GasChem"' // new_line('A') // &
         '  version: "1.0.0"' // new_line('A') // &
         '  active_scheme: ""' // new_line('A') // &
         '  is_active: true' // new_line('A') // &
         '' // new_line('A') // &
         '# Scheme configuration' // new_line('A') // &
         'schemes:' // new_line('A') // &
         '  no_phot:' // new_line('A') // &
         '    description: "MICM gas phase chemistry solver without photolysis"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      scale_factor: 1.0' // new_line('A') // &
         '' // new_line('A') // &
         '# Diagnostic configuration' // new_line('A') // &
         'diagnostics:' // new_line('A') // &
         '  output_frequency: 3600.0  # seconds' // new_line('A') // &
         '  output_diagnostics: true'

   end subroutine get_GasChem_default_config

end module GasChemProcessCreator_Mod