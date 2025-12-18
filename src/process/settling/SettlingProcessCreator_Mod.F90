!> \file SettlingProcessCreator_Mod.F90
!! \brief Factory for creating settling process instances
!!
!! This module provides the factory functions for creating settling
!! process instances following the CATChem Process Factory pattern.
!!
!! Generated on: 2025-12-18T14:12:32.971831
!! Author: Wei Li
!! Version: 1.0.0

module SettlingProcessCreator_Mod

   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, CC_Error, CC_Warning, ErrorManagerType
   use ProcessInterface_Mod
   use ProcessSettlingInterface_Mod

   implicit none
   private

   public :: create_settling_process
   public :: register_settling_process
   public :: get_settling_default_config

contains

   !> Create a new settling process instance
   !!
   !! This factory function creates and returns a new instance of the
   !! settling process. The process is not initialized - the caller
   !! must call the init() method with appropriate configuration.
   !!
   !! @param[out] process     Allocated process instance
   !! @param[out] rc          Return code
   subroutine create_settling_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessSettlingInterface), allocatable :: settling_process
      integer :: alloc_stat

      rc = CC_SUCCESS

      ! Allocate the process instance
      allocate(settling_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = CC_FAILURE
         return
      end if

      ! Move to polymorphic variable
      call move_alloc(settling_process, process)

   end subroutine create_settling_process

   !> Register the settling process with a ProcessManager
   !!
   !! This subroutine registers the settling process with a ProcessManager's
   !! factory. This is the correct way to register processes for use in
   !! applications and integration tests.
   !!
   !! @param[inout] process_mgr The ProcessManager to register with
   !! @param[out] rc Return code
   subroutine register_settling_process(process_mgr, rc)
      use ProcessManager_Mod, only: ProcessManagerType

      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      call process_mgr%register_process( &
         name='settling', &
         category='deposition', &
         description='Process for computing gravitational settling of aerosol species', &
         creator=create_settling_process, &
         rc=rc &
         )

   end subroutine register_settling_process

   !> Get default configuration for settling process
   !!
   !! This function returns a default configuration string that can be
   !! used to initialize the settling process with reasonable defaults.
   !!
   !! @param[out] config_data Default configuration string
   subroutine get_settling_default_config(config_data)
      character(len=*), intent(out) :: config_data

      ! Return default YAML configuration
      config_data = &
         '# Default settling process configuration' // new_line('A') // &
         'process:' // new_line('A') // &
         '  name: "settling"' // new_line('A') // &
         '  version: "1.0.0"' // new_line('A') // &
         '  active_scheme: ""' // new_line('A') // &
         '  is_active: true' // new_line('A') // &
         '' // new_line('A') // &
         '# Scheme configuration' // new_line('A') // &
         'schemes:' // new_line('A') // &
         '  gocart:' // new_line('A') // &
         '    description: "GOCART gravitational settling scheme"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      scale_factor: 1.0' // new_line('A') // &
         '' // new_line('A') // &
         '# Diagnostic configuration' // new_line('A') // &
         'diagnostics:' // new_line('A') // &
         '  output_frequency: 3600.0  # seconds' // new_line('A') // &
         '  output_diagnostics: true'

   end subroutine get_settling_default_config

end module SettlingProcessCreator_Mod
