!> \file modern_initialization_example.F90
!! \brief Example of modern configuration-driven CATChem initialization
!! \date 2025
!!
!! This example demonstrates the modern, configuration-driven approach to
!! initializing CATChem where nspecies and other parameters are read from
!! configuration files rather than hardcoded in the API.
!!
program modern_initialization_example
   use CATChemAPI_Mod
   implicit none

   type(CATChemInstanceType) :: catchem
   type(CATChemConfigType) :: config
   type(CATChemDataType) :: input_data, output_data
   integer :: rc, i
   character(len=256) :: error_msg

   ! ================================================================
   ! MODERN CONFIGURATION-DRIVEN INITIALIZATION
   ! ================================================================

   print *, "=== CATChem Modern Initialization Example ==="
   print *, ""

   ! 1. Setup basic grid configuration
   config%nx = 144
   config%ny = 91
   config%nz = 72
   ! NOTE: nspecies is NOT set here - it will be read from species config file

   ! 2. Specify configuration files (read in this order)
   config%config_file = 'tests/CATChem_config.yml'           ! Main config (first)
   config%species_file = 'tests/Configs/Default/CATChem_species.yml'  ! Species (second)
   config%emission_file = 'external_emission_config.yaml'    ! Emissions (third)

   ! 3. Configure run phases (optional)
   config%enable_run_phases = .true.
   allocate(config%run_phase_names(4))
   config%run_phase_names = ['emissions ', 'chemistry', 'transport', 'deposition']

   ! 4. Configure performance options
   config%use_column_processing = .true.
   config%enable_diagnostics = .true.

   print *, "Configuration setup complete"
   print *, "  Grid: ", config%nx, "x", config%ny, "x", config%nz
   print *, "  Main config file: ", trim(config%config_file)
   print *, "  Species file: ", trim(config%species_file)
   print *, "  Emission file: ", trim(config%emission_file)
   print *, "  Run phases enabled: ", config%enable_run_phases
   print *, ""

   ! ================================================================
   ! INITIALIZATION USING MODERN METHOD
   ! ================================================================

   print *, "Initializing CATChem from configuration files..."

   ! Use modern initialization method that reads config files
   call catchem%init_from_config_files(config, rc)
   if (rc /= CATCHEM_SUCCESS) then
      call catchem%get_error_message(error_msg)
      print *, "ERROR: Initialization failed: ", trim(error_msg)
      print *, "Note: This example requires actual config files to run"
      print *, "      For demo purposes, using legacy initialization..."

      ! Fallback to legacy initialization for demo
      config%nspecies = 25  ! Would normally come from species file
      call catchem%init(config, rc)
      if (rc /= CATCHEM_SUCCESS) then
         print *, "ERROR: Legacy initialization also failed"
         stop 1
      end if
   end if

   print *, "✓ CATChem initialized successfully"
   print *, ""

   ! ================================================================
   ! DEMONSTRATE RUN PHASE CAPABILITIES
   ! ================================================================

   if (config%enable_run_phases) then
      print *, "Setting up run phases..."
      call catchem%setup_run_phases(config%run_phase_names, rc)
      if (rc == CATCHEM_SUCCESS) then
         print *, "✓ Run phases configured successfully"

         ! Show configured phases
         character(len=64), allocatable :: phase_names(:)
         call catchem%get_phase_names(phase_names, rc)
         if (rc == CATCHEM_SUCCESS) then
            print *, "  Configured phases:"
            do i = 1, size(phase_names)
               print *, "    ", i, ": ", trim(phase_names(i))
            end do
         end if
         print *, ""
      else
         print *, "⚠ Run phase setup failed, using single-phase mode"
      end if
   end if

   ! ================================================================
   ! DEMONSTRATE CONFIGURATION-DRIVEN PARAMETERS
   ! ================================================================

   print *, "Configuration-driven parameters:"

   ! Get species names (would come from species config file)
   character(len=64), allocatable :: species_names(:)
   call catchem%get_species_names(species_names, rc)
   if (rc == CATCHEM_SUCCESS .and. allocated(species_names)) then
      print *, "  Number of species: ", size(species_names)
      print *, "  Species read from: ", trim(config%species_file)
   else
      print *, "  Species configuration not available (demo mode)"
   end if

   print *, ""

   ! ================================================================
   ! DEMONSTRATE PROCESS MANAGER CONFIGURATION
   ! ================================================================

   print *, "Configuring process manager..."
   call catchem%configure_process_manager(max_processes=50, &
      enable_column_batching=.true., rc=rc)
   if (rc == CATCHEM_SUCCESS) then
      print *, "✓ Process manager configured"
   else
      print *, "⚠ Process manager configuration failed"
   end if
   print *, ""

   ! ================================================================
   ! ADD PROCESSES AND DEMONSTRATE EXECUTION
   ! ================================================================

   print *, "Adding processes..."

   ! Add some example processes
   call catchem%add_process('dust', rc=rc)
   if (rc == CATCHEM_SUCCESS) then
      print *, "✓ Added dust process"
   else
      print *, "⚠ Failed to add dust process"
   end if

   call catchem%add_process('drydep', rc=rc)
   if (rc == CATCHEM_SUCCESS) then
      print *, "✓ Added dry deposition process"
   else
      print *, "⚠ Failed to add dry deposition process"
   end if

   print *, ""

   ! ================================================================
   ! DEMONSTRATE DIFFERENT EXECUTION MODES
   ! ================================================================

   if (catchem%is_ready_to_run()) then
      print *, "CATChem is ready to run!"
      print *, ""
      print *, "Execution modes available:"

      if (config%enable_run_phases) then
         print *, "  1. Single phase execution:"
         print *, "     call catchem%run_phase('emissions', rc)"
         print *, ""
         print *, "  2. All phases in sequence:"
         print *, "     call catchem%run_all_phases(rc)"
         print *, ""
      end if

      print *, "  3. Traditional timestep (all processes):"
      print *, "     call catchem%run_timestep(input_data, output_data, rc)"
      print *, ""
   else
      print *, "CATChem is not ready to run (missing grid setup or processes)"
   end if

   ! ================================================================
   ! CLEANUP
   ! ================================================================

   print *, "Cleaning up..."
   call catchem%finalize(rc)
   if (rc == CATCHEM_SUCCESS) then
      print *, "✓ CATChem finalized successfully"
   else
      print *, "⚠ Finalization issues detected"
   end if

   print *, ""
   print *, "=== Modern Initialization Example Complete ==="
   print *, ""
   print *, "Key benefits of modern approach:"
   print *, "  • nspecies determined from species config file"
   print *, "  • Multi-phase execution support"
   print *, "  • Configuration-driven initialization"
   print *, "  • Flexible process manager configuration"
   print *, "  • Clear separation of concerns"

end program modern_initialization_example
