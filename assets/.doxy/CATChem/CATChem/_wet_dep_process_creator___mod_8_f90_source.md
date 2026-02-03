

# File WetDepProcessCreator\_Mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**wetdep**](dir_8b9a0ce556ea4a65f6920dfb49dcd69d.md) **>** [**WetDepProcessCreator\_Mod.F90**](_wet_dep_process_creator___mod_8_f90.md)

[Go to the documentation of this file](_wet_dep_process_creator___mod_8_f90.md)


```Fortran


module wetdepprocesscreator_mod

   use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure, cc_error, cc_warning, errormanagertype
   use processinterface_mod
   use processwetdepinterface_mod

   implicit none
   private

   public :: create_wetdep_process
   public :: register_wetdep_process
   public :: get_wetdep_default_config

contains

   subroutine create_wetdep_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessWetDepInterface), allocatable :: wetdep_process
      integer :: alloc_stat

      rc = cc_success

      ! Allocate the process instance
      allocate(wetdep_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = cc_failure
         return
      end if

      ! Move to polymorphic variable
      call move_alloc(wetdep_process, process)

   end subroutine create_wetdep_process

   subroutine register_wetdep_process(process_mgr, rc)
      use processmanager_mod, only: processmanagertype

      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = cc_success

      call process_mgr%register_process( &
         name='wetdep', &
         category='deposition', &
         description='Process for computing wet deposition of gas and aerosol species', &
         creator=create_wetdep_process, &
         rc=rc &
         )

   end subroutine register_wetdep_process

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

end module wetdepprocesscreator_mod
```


