!> \file DryDepProcessCreator_Mod.F90
!! \brief Factory for creating drydep process instances
!!
!! This module provides the factory functions for creating drydep
!! process instances following the CATChem Process Factory pattern.
!!
!! Generated on: 2025-11-14T22:58:26.364759
!! Author: Wei Li
!! Version: 1.0.0

module DryDepProcessCreator_Mod

   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, CC_Error, CC_Warning, ErrorManagerType
   use ProcessInterface_Mod
   use ProcessDryDepInterface_Mod

   implicit none
   private

   public :: create_drydep_process
   public :: register_drydep_process
   public :: get_drydep_default_config

contains

   !> Create a new drydep process instance
   !!
   !! This factory function creates and returns a new instance of the
   !! drydep process. The process is not initialized - the caller
   !! must call the init() method with appropriate configuration.
   !!
   !! @param[out] process     Allocated process instance
   !! @param[out] rc          Return code
   subroutine create_drydep_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessDryDepInterface), allocatable :: drydep_process
      integer :: alloc_stat

      rc = CC_SUCCESS

      ! Allocate the process instance
      allocate(drydep_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = CC_FAILURE
         return
      end if

      ! Move to polymorphic variable
      call move_alloc(drydep_process, process)

   end subroutine create_drydep_process

   !> Register the drydep process with a ProcessManager
   !!
   !! This subroutine registers the drydep process with a ProcessManager's
   !! factory. This is the correct way to register processes for use in
   !! applications and integration tests.
   !!
   !! @param[inout] process_mgr The ProcessManager to register with
   !! @param[out] rc Return code
   subroutine register_drydep_process(process_mgr, rc)
      use ProcessManager_Mod, only: ProcessManagerType

      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      call process_mgr%register_process( &
         name='drydep', &
         category='deposition', &
         description='Process for computing dry deposition of gas and aerosol species', &
         creator=create_drydep_process, &
         rc=rc &
         )

   end subroutine register_drydep_process

   !> Get default configuration for drydep process
   !!
   !! This function returns a default configuration string that can be
   !! used to initialize the drydep process with reasonable defaults.
   !!
   !! @param[out] config_data Default configuration string
   subroutine get_drydep_default_config(config_data)
      character(len=*), intent(out) :: config_data

      ! Return default YAML configuration
      config_data = &
         '# Default drydep process configuration' // new_line('A') // &
         'process:' // new_line('A') // &
         '  name: "drydep"' // new_line('A') // &
         '  version: "1.0.0"' // new_line('A') // &
         '  active_scheme: ""' // new_line('A') // &
         '  is_active: true' // new_line('A') // &
         '' // new_line('A') // &
         '# Scheme configuration' // new_line('A') // &
         'schemes:' // new_line('A') // &
         '  wesely:' // new_line('A') // &
         '    description: "Wesely 1989 gas dry deposition scheme"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      scale_factor: 1.0' // new_line('A') // &
         '' // new_line('A') // &
         '  gocart:' // new_line('A') // &
         '    description: "GOCART-2G aerosol dry deposition scheme"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      scale_factor: 1.0' // new_line('A') // &
         '' // new_line('A') // &
         '  zhang:' // new_line('A') // &
         '    description: "Zhang et al. [2001] scheme with Emerson et al. [2020] updates"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      scale_factor: 1.0' // new_line('A') // &
         '' // new_line('A') // &
         '# Diagnostic configuration' // new_line('A') // &
         'diagnostics:' // new_line('A') // &
         '  output_frequency: 3600.0  # seconds' // new_line('A') // &
         '  output_diagnostics: true'

   end subroutine get_drydep_default_config

end module DryDepProcessCreator_Mod
