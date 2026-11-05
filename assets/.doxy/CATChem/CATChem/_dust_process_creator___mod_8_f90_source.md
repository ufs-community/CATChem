

# File DustProcessCreator\_Mod.F90

[**File List**](files.md) **>** [**dust**](dir_1c14dfbaca1e3f4c2e26e74290119ebd.md) **>** [**DustProcessCreator\_Mod.F90**](_dust_process_creator___mod_8_f90.md)

[Go to the documentation of this file](_dust_process_creator___mod_8_f90.md)


```Fortran


module dustprocesscreator_mod

   use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure, cc_error, cc_warning, errormanagertype
   use processinterface_mod
   use processdustinterface_mod

   implicit none
   private

   public :: create_dust_process
   public :: register_dust_process
   public :: get_dust_default_config

contains

   subroutine create_dust_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessDustInterface), allocatable :: dust_process
      integer :: alloc_stat

      rc = cc_success

      ! Allocate the process instance
      allocate(dust_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = cc_failure
         return
      end if

      ! Move to polymorphic variable
      call move_alloc(dust_process, process)

   end subroutine create_dust_process

   subroutine register_dust_process(process_mgr, rc)
      use processmanager_mod, only: processmanagertype

      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = cc_success

      call process_mgr%register_process( &
         name='dust', &
         category='emission', &
         description='Process for computing windblown dust emissions', &
         creator=create_dust_process, &
         rc=rc &
         )

   end subroutine register_dust_process

   subroutine get_dust_default_config(config_data)
      character(len=*), intent(out) :: config_data

      ! Return default YAML configuration
      config_data = &
         '# Default dust process configuration' // new_line('A') // &
         'process:' // new_line('A') // &
         '  name: "dust"' // new_line('A') // &
         '  version: "1.0.0"' // new_line('A') // &
         '  active_scheme: ""' // new_line('A') // &
         '  is_active: true' // new_line('A') // &
         '' // new_line('A') // &
         '# Scheme configuration' // new_line('A') // &
         'schemes:' // new_line('A') // &
         '  fengsha:' // new_line('A') // &
         '    description: "Fengsha Dust emission scheme developed at NOAA ARL for use at NOAA NWS"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      scale_factor: 1.0' // new_line('A') // &
         '' // new_line('A') // &
         '  ginoux:' // new_line('A') // &
         '    description: "Ginoux dust emission scheme"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      scale_factor: 1.0' // new_line('A') // &
         '' // new_line('A') // &
         '# Diagnostic configuration' // new_line('A') // &
         'diagnostics:' // new_line('A') // &
         '  output_frequency: 3600.0  # seconds' // new_line('A') // &
         '  output_diagnostics: true'

   end subroutine get_dust_default_config

end module dustprocesscreator_mod
```


