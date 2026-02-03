

# File drydep\_example.F90

[**File List**](files.md) **>** [**drydep**](dir_57fb5aa14ddb2cd518a6d90b65ffd000.md) **>** [**examples**](dir_3740d38ba4ad4417edd3d1220f13e03f.md) **>** [**drydep\_example.F90**](drydep__example_8_f90.md)

[Go to the documentation of this file](drydep__example_8_f90.md)


```Fortran


program drydep_example

   use precision_mod, only: fp
   use iso_fortran_env, only: output_unit, error_unit
   use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure
   use processinterface_mod
   use processdrydepinterface_mod
   use drydepcommon_mod
   use drydepprocesscreator_mod
   use statemanager_mod

   implicit none

   ! Process and state management
   class(ProcessInterface), allocatable :: process
   type(StateManagerType) :: state_manager
   integer :: rc

   ! Configuration
   character(len=256) :: config_file = "drydep_config.yaml"

   ! Simulation parameters
   integer, parameter :: n_columns = 10
   integer, parameter :: n_levels = 50
   integer, parameter :: n_time_steps = 24
   real(fp), parameter :: dt = 3600.0_fp  ! 1 hour time step

   ! Working variables
   integer :: i_time, i_col
   real(fp) :: total_time

   write(output_unit, '(A)') "=== DryDep Process Example ==="
   write(output_unit, '(A)') "Process for computing dry deposition of gas and aerosol species"
   write(output_unit, '(A)') "Author: Wei Li"
   write(output_unit, '(A)') ""

   ! Step 1: Initialize state manager
   call state_manager%init(n_levels, n_columns, 1, rc)
   if (rc /= cc_success) then
      write(error_unit, *) 'ERROR: Failed to initialize state manager'
      stop 1
   end if

   ! Step 2: Create process instance
   call create_process(process, rc)
   if (rc /= cc_success) then
      write(error_unit, *) 'ERROR: Failed to create process'
      stop 1
   end if

   ! Step 3: Initialize process
   call process%init(state_manager, rc)
   if (rc /= cc_success) then
      write(error_unit, *) 'ERROR: Process initialization failed'
      stop 1
   end if

   write(output_unit, '(A)') "Process initialized successfully"

   ! Step 4: Time loop simulation
   write(output_unit, '(A,I0,A)') "Starting ", n_time_steps, " time step simulation"

   do i_time = 1, n_time_steps
      total_time = real(i_time - 1, fp) * dt

      write(output_unit, '(A,I0,A,F8.1,A)') "Time step ", i_time, " (t=", total_time/3600.0_fp, " hours)"

      ! Run process for this time step
      call process%run(state_manager, rc)
      if (rc /= cc_success) then
         write(error_unit, *) 'ERROR: Process execution failed at time step ', i_time
         stop 1
      end if
   end do

   write(output_unit, '(A)') "Simulation completed successfully"

   ! Step 5: Finalize process
   call process%finalize(rc)
   if (rc /= cc_success) then
      write(error_unit, *) 'ERROR: Process finalization failed'
      stop 1
   end if

   write(output_unit, '(A)') "Process finalized successfully"
   write(output_unit, '(A)') "Example completed!"

end program drydep_example
write(error_unit, '(A)') "Error creating drydep process"
call error_handler%print_errors()
stop 1
end if
write(output_unit, '(A)') "5. DryDep process created"

 ! Step 6: Load configuration
call load_configuration(config_data, error_handler)
if (error_handler%has_error()) then
   write(error_unit, '(A)') "Error loading configuration"
   call error_handler%print_errors()
   stop 1
end if
write(output_unit, '(A)') "6. Configuration loaded"

 ! Step 7: Initialize process
call process%init(state_manager, config_data, error_handler)
if (error_handler%has_error()) then
   write(error_unit, '(A)') "Error initializing drydep process"
   call error_handler%print_errors()
   stop 1
end if
write(output_unit, '(A)') "7. DryDep process initialized"

 ! Step 8: Print process information
call print_process_info(process)

 ! Step 9: Set up initial conditions
call setup_initial_conditions(state_manager, error_handler)
if (error_handler%has_error()) then
   write(error_unit, '(A)') "Error setting up initial conditions"
   call error_handler%print_errors()
   stop 1
end if
write(output_unit, '(A)') "8. Initial conditions set"

 ! Step 10: Time integration loop
write(output_unit, '(A)') "9. Starting time integration..."
write(output_unit, '(A,I0,A,F0.1,A)') "   Running ", n_time_steps, &
   " time steps with dt = ", dt, " seconds"

total_time = 0.0_fp
do i_time = 1, n_time_steps

   ! Run process for one time step
   call process%run(state_manager, dt, error_handler)
   if (error_handler%has_error()) then
      write(error_unit, '(A,I0)') "Error during time step ", i_time
      call error_handler%print_errors()
      exit
   end if

   total_time = total_time + dt

   ! Print progress
   if (mod(i_time, max(n_time_steps/10, 1)) == 0) then
      write(output_unit, '(A,I0,A,I0,A,F0.1,A)') "   Step ", i_time, &
         "/", n_time_steps, " (t = ", total_time/3600.0_fp, " hours)"
   end if

end do

if (.not. error_handler%has_error()) then
   write(output_unit, '(A)') "10. Time integration completed successfully"
end if

 ! Step 11: Print diagnostic summary
diag_mgr => state_manager%get_diagnostic_manager()
call print_diagnostic_summary(diag_mgr)

 ! Step 12: Clean up
call process%finalize(error_handler)
call state_manager%finalize(error_handler)

write(output_unit, '(A)') "11. Cleanup completed"
write(output_unit, '(A)') ""
write(output_unit, '(A)') "=== Example completed successfully ==="

contains

subroutine setup_chemical_species(state_manager, error_handler)
   type(StateManagerType), intent(inout) :: state_manager
   type(ErrorHandler), intent(inout) :: error_handler

   ! Add generic species
   call state_manager%add_species('GENERIC_SPECIES', error_handler)

   ! Initialize species concentrations
   call state_manager%allocate_species_arrays(error_handler)

end subroutine setup_chemical_species

subroutine setup_meteorological_fields(state_manager, error_handler)
   type(StateManagerType), intent(inout) :: state_manager
   type(ErrorHandler), intent(inout) :: error_handler

   ! Add required meteorological fields

   ! Add optional meteorological fields

   ! Initialize meteorological data with test values
   call initialize_test_meteorology(state_manager, error_handler)

end subroutine setup_meteorological_fields

subroutine initialize_test_meteorology(state_manager, error_handler)
   type(StateManagerType), intent(inout) :: state_manager
   type(ErrorHandler), intent(inout) :: error_handler

   integer :: i_col, i_lev
   real(fp) :: height, latitude, longitude

   ! Set test meteorological values
   do i_col = 1, n_columns
      do i_lev = 1, n_levels

         height = real(i_lev - 1, fp) * 1000.0_fp  ! Height in meters
         latitude = 45.0_fp + real(i_col - 1, fp) * 1.0_fp  ! Latitude
         longitude = -120.0_fp + real(i_col - 1, fp) * 1.0_fp  ! Longitude


      end do
   end do

end subroutine initialize_test_meteorology

subroutine load_configuration(config_data, error_handler)
   character(len=*), intent(out) :: config_data
   type(ErrorHandler), intent(inout) :: error_handler

   ! For this example, use inline configuration
   ! In practice, this would be loaded from a file
   config_data = &
      'process:' // new_line('A') // &
      '  name: "drydep"' // new_line('A') // &
      '  version: "1.0.0"' // new_line('A') // &
      '  active_scheme: ""' // new_line('A') // &
      '  is_active: true' // new_line('A') // &
      'diagnostics:' // new_line('A') // &
      '  output_frequency: 3600.0' // new_line('A') // &
      '  output_diagnostics: true'

end subroutine load_configuration

subroutine setup_initial_conditions(state_manager, error_handler)
   type(StateManagerType), intent(inout) :: state_manager
   type(ErrorHandler), intent(inout) :: error_handler

   integer :: i_col, i_lev, i_spec
   real(fp) :: initial_concentration

   ! Set initial chemical concentrations

end subroutine setup_initial_conditions

subroutine print_process_info(process)
   class(ProcessInterface), intent(in) :: process

   character(len=32), allocatable :: species_list(:)
   character(len=32), allocatable :: required_fields(:)
   integer :: i

   write(output_unit, '(A)') ""
   write(output_unit, '(A)') "Process Information:"
   write(output_unit, '(A,A)') "  Name: ", trim(process%get_name())
   write(output_unit, '(A,A)') "  Version: ", trim(process%get_version())
   write(output_unit, '(A,A)') "  Description: ", trim(process%get_description())

   ! Print species list
   species_list = process%get_species_list()
   write(output_unit, '(A,I0)') "  Number of species: ", size(species_list)
   if (size(species_list) > 0) then
      write(output_unit, '(A)', advance='no') "  Species: "
      do i = 1, size(species_list)
         write(output_unit, '(A)', advance='no') trim(species_list(i))
         if (i < size(species_list)) write(output_unit, '(A)', advance='no') ", "
      end do
      write(output_unit, '(A)') ""
   end if

   ! Print required fields
   required_fields = process%get_required_met_fields()
   write(output_unit, '(A,I0)') "  Required met fields: ", size(required_fields)
   if (size(required_fields) > 0) then
      do i = 1, size(required_fields)
         write(output_unit, '(A,A)') "    - ", trim(required_fields(i))
      end do
   end if
   write(output_unit, '(A)') ""

end subroutine print_process_info

subroutine print_diagnostic_summary(diag_mgr)
   type(DiagnosticManager), intent(in) :: diag_mgr

   integer :: n_diagnostics, i
   character(len=64), allocatable :: diag_names(:)
   real(fp), allocatable :: diag_values(:)

   write(output_unit, '(A)') ""
   write(output_unit, '(A)') "Diagnostic Summary:"

   n_diagnostics = diag_mgr%get_n_diagnostics()
   write(output_unit, '(A,I0)') "  Number of diagnostics: ", n_diagnostics

   if (n_diagnostics > 0) then
      allocate(diag_names(n_diagnostics))
      allocate(diag_values(n_diagnostics))

      call diag_mgr%get_diagnostic_names(diag_names)
      call diag_mgr%get_diagnostic_values(diag_values)

      do i = 1, n_diagnostics
         write(output_unit, '(A,A,A,ES12.4)') "    ", &
            trim(diag_names(i)), ": ", diag_values(i)
      end do
   end if
   write(output_unit, '(A)') ""

end subroutine print_diagnostic_summary

end program drydep_example
```


