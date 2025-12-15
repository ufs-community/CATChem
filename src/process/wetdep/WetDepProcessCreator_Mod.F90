!> \file WetDepProcessCreator_Mod.F90
!! \brief Factory for creating wetdep process instances
!!
!! This module provides the factory functions for creating wetdep
!! process instances following the CATChem Process Factory pattern.
!!
!! Generated on: 2025-12-15T16:30:33.627354
!! Author: Wei Li
!! Version: 1.0.0

module WetDepProcessCreator_Mod

   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, CC_Error, CC_Warning, ErrorManagerType
   use ProcessInterface_Mod
   use ProcessWetDepInterface_Mod

   implicit none
   private

   public :: create_wetdep_process
   public :: register_wetdep_process
   public :: get_wetdep_default_config

contains

   !> Create a new wetdep process instance
   !!
   !! This factory function creates and returns a new instance of the
   !! wetdep process. The process is not initialized - the caller
   !! must call the init() method with appropriate configuration.
   !!
   !! @param[out] process     Allocated process instance
   !! @param[out] rc          Return code
   subroutine create_wetdep_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessWetDepInterface), allocatable :: wetdep_process
      integer :: alloc_stat

      rc = CC_SUCCESS

      ! Allocate the process instance
      allocate(wetdep_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = CC_FAILURE
         return
      end if

      ! Move to polymorphic variable
      call move_alloc(wetdep_process, process)

   end subroutine create_wetdep_process

   !> Register the wetdep process with a ProcessManager
   !!
   !! This subroutine registers the wetdep process with a ProcessManager's
   !! factory. This is the correct way to register processes for use in
   !! applications and integration tests.
   !!
   !! @param[inout] process_mgr The ProcessManager to register with
   !! @param[out] rc Return code
   subroutine register_wetdep_process(process_mgr, rc)
      use ProcessManager_Mod, only: ProcessManagerType

      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      call process_mgr%register_process( &
         name='wetdep', &
         category='deposition', &
         description='Process for computing wet deposition of gas and aerosol species', &
         creator=create_wetdep_process, &
         rc=rc &
         )

   end subroutine register_wetdep_process

   !> Get default configuration for wetdep process
   !!
   !! This function returns a default configuration string that can be
   !! used to initialize the wetdep process with reasonable defaults.
   !!
   !! @param[out] config_data Default configuration string
   subroutine get_wetdep_default_config(config_data)
      character(len=*), intent(out) :: config_data

      ! Return default YAML configuration
      config_data = &
         '# Default wetdep process configuration' // new_line('A') // &
         'process:' // new_line('A') // &
         '  name: "wetdep"' // new_line('A') // &
         '  version: "1.0.0"' // new_line('A') // &
         '  active_scheme: ""' // new_line('A') // &
         '  is_active: true' // new_line('A') // &
         '' // new_line('A') // &
         '# Scheme configuration' // new_line('A') // &
         'schemes:' // new_line('A') // &
         '  jacob:' // new_line('A') // &
         '    description: "Jacob et al. [2000] wet deposition scheme"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      scale_factor: 1.0' // new_line('A') // &
         '' // new_line('A') // &
         '# Diagnostic configuration' // new_line('A') // &
         'diagnostics:' // new_line('A') // &
         '  output_frequency: 3600.0  # seconds' // new_line('A') // &
         '  output_diagnostics: true'

   end subroutine get_wetdep_default_config

end module WetDepProcessCreator_Mod
