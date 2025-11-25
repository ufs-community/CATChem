

# File SeaSaltProcessCreator\_Mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**seasalt**](dir_3e6fd2e121e43ca7d4114b6c0b4e05b6.md) **>** [**SeaSaltProcessCreator\_Mod.F90**](_sea_salt_process_creator___mod_8_f90.md)

[Go to the documentation of this file](_sea_salt_process_creator___mod_8_f90.md)


```Fortran


module seasaltprocesscreator_mod

   use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure, cc_error, cc_warning, errormanagertype
   use processinterface_mod
   use processseasaltinterface_mod

   implicit none
   private

   public :: create_seasalt_process
   public :: register_seasalt_process
   public :: get_seasalt_default_config

contains

   subroutine create_seasalt_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessSeaSaltInterface), allocatable :: seasalt_process
      integer :: alloc_stat

      rc = cc_success

      ! Allocate the process instance
      allocate(seasalt_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = cc_failure
         return
      end if

      ! Move to polymorphic variable
      call move_alloc(seasalt_process, process)

   end subroutine create_seasalt_process

   subroutine register_seasalt_process(process_mgr, rc)
      use processmanager_mod, only: processmanagertype

      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = cc_success

      call process_mgr%register_process( &
         name='seasalt', &
         category='emission', &
         description='Process for computing sea salt aerosol emissions over ocean surfaces', &
         creator=create_seasalt_process, &
         rc=rc &
         )

   end subroutine register_seasalt_process

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

end module seasaltprocesscreator_mod
```


