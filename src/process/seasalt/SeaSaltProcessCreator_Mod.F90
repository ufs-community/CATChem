!> \file SeaSaltProcessCreator_Mod.F90
!! \brief Factory for creating seasalt process instances
!!
!! This module provides the factory functions for creating seasalt
!! process instances following the CATChem Process Factory pattern.
!!
!! Generated on: 2025-11-14T23:01:21.603774
!! Author: Barry Baker & Wei Li
!! Version: 1.0.0

module SeaSaltProcessCreator_Mod

   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, CC_Error, CC_Warning, ErrorManagerType
   use ProcessInterface_Mod
   use ProcessSeaSaltInterface_Mod

   implicit none
   private

   public :: create_seasalt_process
   public :: register_seasalt_process
   public :: get_seasalt_default_config

contains

   !> Create a new seasalt process instance
   !!
   !! This factory function creates and returns a new instance of the
   !! seasalt process. The process is not initialized - the caller
   !! must call the init() method with appropriate configuration.
   !!
   !! @param[out] process     Allocated process instance
   !! @param[out] rc          Return code
   subroutine create_seasalt_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessSeaSaltInterface), allocatable :: seasalt_process
      integer :: alloc_stat

      rc = CC_SUCCESS

      ! Allocate the process instance
      allocate(seasalt_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = CC_FAILURE
         return
      end if

      ! Move to polymorphic variable
      call move_alloc(seasalt_process, process)

   end subroutine create_seasalt_process

   !> Register the seasalt process with a ProcessManager
   !!
   !! This subroutine registers the seasalt process with a ProcessManager's
   !! factory. This is the correct way to register processes for use in
   !! applications and integration tests.
   !!
   !! @param[inout] process_mgr The ProcessManager to register with
   !! @param[out] rc Return code
   subroutine register_seasalt_process(process_mgr, rc)
      use ProcessManager_Mod, only: ProcessManagerType

      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      call process_mgr%register_process( &
         name='seasalt', &
         category='emission', &
         description='Process for computing sea salt aerosol emissions over ocean surfaces', &
         creator=create_seasalt_process, &
         rc=rc &
         )

   end subroutine register_seasalt_process

   !> Get default configuration for seasalt process
   !!
   !! This function returns a default configuration string that can be
   !! used to initialize the seasalt process with reasonable defaults.
   !!
   !! @param[out] config_data Default configuration string
   subroutine get_seasalt_default_config(config_data)
      character(len=*), intent(out) :: config_data

      ! Return default YAML configuration
      config_data = &
         '# Default seasalt process configuration' // new_line('A') // &
         'process:' // new_line('A') // &
         '  name: "seasalt"' // new_line('A') // &
         '  version: "1.0.0"' // new_line('A') // &
         '  active_scheme: ""' // new_line('A') // &
         '  is_active: true' // new_line('A') // &
         '' // new_line('A') // &
         '# Scheme configuration' // new_line('A') // &
         'schemes:' // new_line('A') // &
         '  gong97:' // new_line('A') // &
         '    description: "Gong 1997 sea salt emission scheme"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      scale_factor: 1.0' // new_line('A') // &
         '      weibull_flag: False' // new_line('A') // &
         '' // new_line('A') // &
         '  gong03:' // new_line('A') // &
         '    description: "Gong 2003 sea salt emission scheme with improved sub- and super-micron treatment"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      scale_factor: 1.0' // new_line('A') // &
         '      weibull_flag: False' // new_line('A') // &
         '' // new_line('A') // &
         '  geos12:' // new_line('A') // &
         '    description: "GEOS-Chem 2012 sea salt emission scheme with observational constraints"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      scale_factor: 1.0' // new_line('A') // &
         '' // new_line('A') // &
         '# Diagnostic configuration' // new_line('A') // &
         'diagnostics:' // new_line('A') // &
         '  output_frequency: 3600.0  # seconds' // new_line('A') // &
         '  output_diagnostics: true'

   end subroutine get_seasalt_default_config

end module SeaSaltProcessCreator_Mod
