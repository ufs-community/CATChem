

# File SO4chemProcessCreator\_Mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**so4chem**](dir_fb8fc0df5ebe1b02f5e46b98d91cbc63.md) **>** [**SO4chemProcessCreator\_Mod.F90**](_s_o4chem_process_creator___mod_8_f90.md)

[Go to the documentation of this file](_s_o4chem_process_creator___mod_8_f90.md)


```Fortran


module so4chemprocesscreator_mod

   use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure, cc_error, cc_warning, errormanagertype
   use processinterface_mod
   use processso4cheminterface_mod

   implicit none
   private

   public :: create_so4chem_process
   public :: register_so4chem_process
   public :: get_so4chem_default_config

contains

   subroutine create_so4chem_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessSO4chemInterface), allocatable :: so4chem_process
      integer :: alloc_stat

      rc = cc_success

      ! Allocate the process instance
      allocate(so4chem_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = cc_failure
         return
      end if

      ! Move to polymorphic variable
      call move_alloc(so4chem_process, process)

   end subroutine create_so4chem_process

   subroutine register_so4chem_process(process_mgr, rc)
      use processmanager_mod, only: processmanagertype

      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = cc_success

      call process_mgr%register_process( &
         name='so4chem', &
         category='chemistry', &
         description='Process for computing chemical production of sulfate from SO2 oxidation', &
         creator=create_so4chem_process, &
         rc=rc &
         )

   end subroutine register_so4chem_process

   subroutine get_so4chem_default_config(config_data)
      character(len=*), intent(out) :: config_data

      ! Return default YAML configuration
      config_data = &
         '# Default so4chem process configuration' // new_line('A') // &
         'process:' // new_line('A') // &
         '  name: "so4chem"' // new_line('A') // &
         '  version: "1.0.0"' // new_line('A') // &
         '  active_scheme: ""' // new_line('A') // &
         '  is_active: true' // new_line('A') // &
         '' // new_line('A') // &
         '# Scheme configuration' // new_line('A') // &
         'schemes:' // new_line('A') // &
         '  gocart:' // new_line('A') // &
         '    description: "GOCART SO2 to SO4 production scheme"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      scale_factor: 1.0' // new_line('A') // &
         '' // new_line('A') // &
         '# Diagnostic configuration' // new_line('A') // &
         'diagnostics:' // new_line('A') // &
         '  output_frequency: 3600.0  # seconds' // new_line('A') // &
         '  output_diagnostics: true'

   end subroutine get_so4chem_default_config

end module so4chemprocesscreator_mod
```


