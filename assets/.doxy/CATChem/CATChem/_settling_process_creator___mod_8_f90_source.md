

# File SettlingProcessCreator\_Mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**settling**](dir_1a0bba2ffdf6e6637fcb76856471cb75.md) **>** [**SettlingProcessCreator\_Mod.F90**](_settling_process_creator___mod_8_f90.md)

[Go to the documentation of this file](_settling_process_creator___mod_8_f90.md)


```Fortran


module settlingprocesscreator_mod

   use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure, cc_error, cc_warning, errormanagertype
   use processinterface_mod
   use processsettlinginterface_mod

   implicit none
   private

   public :: create_settling_process
   public :: register_settling_process
   public :: get_settling_default_config

contains

   subroutine create_settling_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessSettlingInterface), allocatable :: settling_process
      integer :: alloc_stat

      rc = cc_success

      ! Allocate the process instance
      allocate(settling_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = cc_failure
         return
      end if

      ! Move to polymorphic variable
      call move_alloc(settling_process, process)

   end subroutine create_settling_process

   subroutine register_settling_process(process_mgr, rc)
      use processmanager_mod, only: processmanagertype

      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = cc_success

      call process_mgr%register_process( &
         name='settling', &
         category='deposition', &
         description='Process for computing gravitational settling of aerosol species', &
         creator=create_settling_process, &
         rc=rc &
         )

   end subroutine register_settling_process

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

end module settlingprocesscreator_mod
```


