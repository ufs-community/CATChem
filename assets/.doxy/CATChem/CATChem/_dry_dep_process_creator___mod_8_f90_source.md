

# File DryDepProcessCreator\_Mod.F90

[**File List**](files.md) **>** [**drydep**](dir_57fb5aa14ddb2cd518a6d90b65ffd000.md) **>** [**DryDepProcessCreator\_Mod.F90**](_dry_dep_process_creator___mod_8_f90.md)

[Go to the documentation of this file](_dry_dep_process_creator___mod_8_f90.md)


```Fortran


module drydepprocesscreator_mod

   use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure, cc_error, cc_warning, errormanagertype
   use processinterface_mod
   use processdrydepinterface_mod

   implicit none
   private

   public :: create_drydep_process
   public :: register_drydep_process
   public :: get_drydep_default_config

contains

   subroutine create_drydep_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessDryDepInterface), allocatable :: drydep_process
      integer :: alloc_stat

      rc = cc_success

      ! Allocate the process instance
      allocate(drydep_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = cc_failure
         return
      end if

      ! Move to polymorphic variable
      call move_alloc(drydep_process, process)

   end subroutine create_drydep_process

   subroutine register_drydep_process(process_mgr, rc)
      use processmanager_mod, only: processmanagertype

      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = cc_success

      call process_mgr%register_process( &
         name='drydep', &
         category='deposition', &
         description='Process for computing dry deposition of gas and aerosol species', &
         creator=create_drydep_process, &
         rc=rc &
         )

   end subroutine register_drydep_process

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

end module drydepprocesscreator_mod
```


