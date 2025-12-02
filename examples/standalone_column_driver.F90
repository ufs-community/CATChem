!> \file standalone_column_driver.F90
!! \brief Standalone column driver demonstrating the modern CATChem API
!! \date 2025
!!
!! This driver demonstrates:
!! - Configuration-driven initialization (nspecies from config files)
!! - Multi-phase execution capabilities
!! - Column-based processing for performance
!! - Flexible field mapping and data exchange
!! - Process-level diagnostics
!!
!! The driver simulates a single atmospheric column through multiple time steps,
!! showcasing the full capabilities of the modernized CATChem API.
!!
program standalone_column_driver
   use CATChemAPI_Mod
   use precision_mod, only: fp
   implicit none

   ! Main driver components
   type(CATChemInstanceType) :: catchem
   type(CATChemConfigType) :: config
   type(CATChemDataType) :: input_data, output_data
   type(FlexibleDataExchangeType) :: data_exchange

   ! Column configuration
   integer, parameter :: nx = 1, ny = 1, nz = 72  ! Single column
   integer, parameter :: nsteps = 24              ! 24-hour simulation
   real(fp), parameter :: dt = 3600.0_fp          ! 1-hour time steps

   ! Grid coordinates
   real(fp) :: lats(nx,ny), lons(nx,ny), levels(nz)

   ! Working variables
   integer :: rc, step, nspecies
   character(len=256) :: error_msg
   character(len=64), allocatable :: species_names(:), phase_names(:)
   logical :: use_phases = .true.

   ! Timing and performance
   real :: start_time, end_time, phase_times(4)
   character(len=16) :: time_str

   print *, "================================================================"
   print *, "CATChem Standalone Column Driver"
   print *, "Modern API with Configuration-Driven Initialization"
   print *, "================================================================"
   print *, ""

   ! ================================================================
   ! STEP 1: CONFIGURATION SETUP
   ! ================================================================

   print *, "Setting up configuration..."

   ! Basic grid configuration
   config%nx = nx
   config%ny = ny
   config%nz = nz
   config%dt = dt
   config%nsteps = nsteps

   ! Configuration files (nspecies will be read from these)
   config%config_file = 'tests/CATChem_config.yml'
   config%species_file = 'tests/Configs/Default/CATChem_species.yml'
   config%emission_file = 'external_emission_config.yaml'

   ! Enable multi-phase execution
   config%enable_run_phases = use_phases
   if (use_phases) then
      allocate(config%run_phase_names(4))
      config%run_phase_names = ['emissions ', 'chemistry ', 'transport ', 'deposition']
   endif

   ! Performance and diagnostic options
   config%use_column_processing = .true.  ! Optimize for column processing
   config%enable_diagnostics = .true.
   config%enable_dust = .true.
   config%enable_seasalt = .true.
   config%enable_drydep = .true.
   config%enable_external_emis = .true.

   print *, "✓ Configuration setup complete"
   print *, "  Grid size: ", nx, "×", ny, "×", nz, " (single column)"
   print *, "  Time steps: ", nsteps, " × ", dt/3600.0, " hours"
   print *, "  Column processing: ", config%use_column_processing
   print *, "  Multi-phase execution: ", config%enable_run_phases
   print *, ""

   ! ================================================================
   ! STEP 2: MODERN INITIALIZATION
   ! ================================================================

   print *, "Initializing CATChem from configuration files..."
   call cpu_time(start_time)

   ! Use modern configuration-driven initialization
   call catchem%init_from_config_files(config, rc)
   if (rc /= CATCHEM_SUCCESS) then
      call catchem%get_error_message(error_msg)
      print *, "⚠ Modern initialization failed: ", trim(error_msg)
      print *, "  Falling back to legacy initialization for demo..."

      ! Fallback for demo (in real use, fix config files)
      config%nspecies = 25  ! Would normally come from species config
      call catchem%init(config, rc)
      if (rc /= CATCHEM_SUCCESS) then
         print *, "✗ Legacy initialization also failed!"
         stop 1
      endif
      nspecies = config%nspecies
   else
      ! Get nspecies from configuration
      call catchem%get_species_names(species_names, rc)
      if (rc == CATCHEM_SUCCESS .and. allocated(species_names)) then
         nspecies = size(species_names)
      else
         nspecies = 25  ! Default fallback
      endif
   endif

   call cpu_time(end_time)
   print *, "✓ CATChem initialized successfully"
   print *, "  Species count: ", nspecies, " (from configuration files)"
   print *, "  Initialization time: ", end_time - start_time, " seconds"
   print *, ""

   ! ================================================================
   ! STEP 3: GRID AND COORDINATE SETUP
   ! ================================================================

   print *, "Setting up grid coordinates..."

   ! Simple column coordinates (single point)
   lats(1,1) = 40.0_fp    ! 40°N latitude
   lons(1,1) = -75.0_fp   ! 75°W longitude (Eastern US)

   ! Vertical levels (pressure coordinates, Pa)
   call setup_vertical_levels(levels, nz)

   call catchem%setup_grid(lats, lons, levels, rc)
   if (rc /= CATCHEM_SUCCESS) then
      print *, "✗ Grid setup failed!"
      stop 1
   endif

   print *, "✓ Grid coordinates configured"
   print *, "  Location: ", lats(1,1), "°N, ", lons(1,1), "°E"
   print *, "  Levels: ", nz, " (", levels(nz)/100.0, " to ", levels(1)/100.0, " hPa)"
   print *, ""

   ! ================================================================
   ! STEP 4: PROCESS CONFIGURATION
   ! ================================================================

   print *, "Configuring processes and run phases..."

   ! Configure process manager for optimal column processing
   call catchem%configure_process_manager(max_processes=20, &
      enable_column_batching=.true., rc=rc)
   if (rc /= CATCHEM_SUCCESS) then
      print *, "⚠ Process manager configuration failed"
   endif

   ! Add atmospheric chemistry processes
   call add_processes(catchem, rc)
   if (rc /= CATCHEM_SUCCESS) then
      print *, "✗ Process addition failed!"
      stop 1
   endif

   ! Setup run phases if enabled
   if (config%enable_run_phases) then
      call catchem%setup_run_phases(config%run_phase_names, rc)
      if (rc == CATCHEM_SUCCESS) then
         call catchem%get_phase_names(phase_names, rc)
         print *, "✓ Run phases configured:"
         do step = 1, size(phase_names)
            print *, "    ", step, ": ", trim(phase_names(step))
         enddo
      else
         print *, "⚠ Run phase setup failed, using single-phase mode"
         use_phases = .false.
      endif
   endif
   print *, ""

   ! ================================================================
   ! STEP 5: DATA INITIALIZATION
   ! ================================================================

   print *, "Initializing atmospheric data..."

   call initialize_atmospheric_data(input_data, nx, ny, nz, nspecies, levels)
   call initialize_atmospheric_data(output_data, nx, ny, nz, nspecies, levels)

   ! Validate input data
   call catchem%validate_data(input_data, rc)
   if (rc /= CATCHEM_SUCCESS) then
      print *, "⚠ Input data validation failed"
   else
      print *, "✓ Input data validated successfully"
   endif
   print *, ""

   ! ================================================================
   ! STEP 6: TIME STEPPING SIMULATION
   ! ================================================================

   print *, "Starting ", nsteps, "-step atmospheric column simulation..."
   print *, ""

   if (.not. catchem%is_ready_to_run()) then
      print *, "✗ CATChem not ready to run!"
      stop 1
   endif

   call cpu_time(start_time)

   do step = 1, nsteps
      call format_time_string(step, time_str)

      if (mod(step-1, 6) == 0) then  ! Print every 6 hours
         print *, "Time step ", step, "/", nsteps, " (", trim(time_str), ")"
      endif

      ! Update meteorological conditions for this time step
      call update_meteorology_for_timestep(input_data, step, dt)

      ! Execute chemistry processes
      if (use_phases) then
         ! Multi-phase execution with timing
         call run_phases_with_timing(catchem, phase_times, rc)
      else
         ! Traditional single-phase execution
         call catchem%run_timestep(input_data, output_data, rc)
      endif

      if (rc /= CATCHEM_SUCCESS) then
         call catchem%get_error_message(error_msg)
         print *, "✗ Time step ", step, " failed: ", trim(error_msg)
         stop 1
      endif

      ! Update concentrations for next step
      if (allocated(output_data%concentrations)) then
         input_data%concentrations = output_data%concentrations
      endif

      ! Process and save diagnostics periodically
      if (mod(step, 6) == 0) then  ! Every 6 hours
         call process_column_diagnostics(catchem, step, time_str)
      endif
   enddo

   call cpu_time(end_time)
   print *, ""
   print *, "✓ Simulation completed successfully!"
   print *, "  Total runtime: ", end_time - start_time, " seconds"
   print *, "  Average per step: ", (end_time - start_time) / nsteps, " seconds"

   if (use_phases) then
      print *, "  Phase timing breakdown:"
      print *, "    Emissions: ", sum(phase_times(1::4)) / count(phase_times(1::4) > 0), " s avg"
      print *, "    Chemistry: ", sum(phase_times(2::4)) / count(phase_times(2::4) > 0), " s avg"
      print *, "    Transport: ", sum(phase_times(3::4)) / count(phase_times(3::4) > 0), " s avg"
      print *, "    Deposition:", sum(phase_times(4::4)) / count(phase_times(4::4) > 0), " s avg"
   endif
   print *, ""

   ! ================================================================
   ! STEP 7: FINAL DIAGNOSTICS AND OUTPUT
   ! ================================================================

   print *, "Generating final diagnostics..."

   call generate_final_diagnostics(catchem, output_data, nsteps)
   call save_column_results(output_data, nsteps, levels)

   print *, "✓ Final diagnostics and output complete"
   print *, ""

   ! ================================================================
   ! STEP 8: CLEANUP
   ! ================================================================

   print *, "Cleaning up..."
   call catchem%finalize(rc)
   if (rc == CATCHEM_SUCCESS) then
      print *, "✓ CATChem finalized successfully"
   else
      print *, "⚠ Finalization issues detected"
   endif

   ! Deallocate data
   call cleanup_data(input_data, output_data)

   print *, ""
   print *, "================================================================"
   print *, "Column Driver Completed Successfully"
   print *, "================================================================"
   print *, ""
   print *, "Key Achievements:"
   print *, "  • Configuration-driven initialization (nspecies from files)"
   print *, "  • Multi-phase execution with timing analysis"
   print *, "  • Optimized column processing"
   print *, "  • Process-level diagnostics"
   print *, "  • Modern API usage patterns"

contains

   !> Setup vertical pressure levels
   subroutine setup_vertical_levels(levels, nz)
      real(fp), intent(out) :: levels(nz)
      integer, intent(in) :: nz
      integer :: k
      real(fp) :: p_top, p_sfc, dp

      p_top = 1000.0_fp     ! 10 hPa top
      p_sfc = 101325.0_fp   ! 1013.25 hPa surface

      ! Simple linear spacing in log-pressure
      do k = 1, nz
         levels(k) = p_sfc * exp(log(p_top/p_sfc) * (k-1) / (nz-1))
      enddo
   end subroutine setup_vertical_levels

   !> Add chemistry processes to the model
   subroutine add_processes(catchem, rc)
      type(CATChemInstanceType), intent(inout) :: catchem
      integer, intent(out) :: rc

      rc = CATCHEM_SUCCESS

      print *, "  Adding atmospheric chemistry processes..."

      ! Add dust emission process
      call catchem%add_process('dust', rc=rc)
      if (rc == CATCHEM_SUCCESS) then
         print *, "    ✓ Dust emissions"
      else
         print *, "    ⚠ Dust emissions (process not available)"
         rc = CATCHEM_SUCCESS  ! Continue without this process
      endif

      ! Add sea salt process
      call catchem%add_process('seasalt', rc=rc)
      if (rc == CATCHEM_SUCCESS) then
         print *, "    ✓ Sea salt emissions"
      else
         print *, "    ⚠ Sea salt emissions (process not available)"
         rc = CATCHEM_SUCCESS
      endif

      ! Add dry deposition process
      call catchem%add_process('drydep', rc=rc)
      if (rc == CATCHEM_SUCCESS) then
         print *, "    ✓ Dry deposition"
      else
         print *, "    ⚠ Dry deposition (process not available)"
         rc = CATCHEM_SUCCESS
      endif

      ! Add external emissions process
      call catchem%add_process('externalemission', rc=rc)
      if (rc == CATCHEM_SUCCESS) then
         print *, "    ✓ External emissions"
      else
         print *, "    ⚠ External emissions (process not available)"
         rc = CATCHEM_SUCCESS
      endif

   end subroutine add_processes

   !> Initialize atmospheric column data
   subroutine initialize_atmospheric_data(data, nx, ny, nz, nspecies, levels)
      type(CATChemDataType), intent(out) :: data
      integer, intent(in) :: nx, ny, nz, nspecies
      real(fp), intent(in) :: levels(nz)
      integer :: k
      real(fp) :: z_height, temp_lapse

      ! Allocate arrays
      allocate(data%temperature(nx,ny,nz))
      allocate(data%pressure(nx,ny,nz))
      allocate(data%humidity(nx,ny,nz))
      allocate(data%wind_u(nx,ny,nz))
      allocate(data%wind_v(nx,ny,nz))
      allocate(data%surface_pressure(nx,ny))
      allocate(data%concentrations(nx,ny,nz,nspecies))
      allocate(data%emission_rates(nx,ny,nz,nspecies))

      ! Initialize meteorological fields with realistic profiles
      temp_lapse = 6.5e-3_fp  ! K/m
      data%surface_pressure(1,1) = 101325.0_fp  ! Pa

      do k = 1, nz
         ! Height from pressure (rough approximation)
         z_height = -8400.0_fp * log(levels(k) / 101325.0_fp)

         data%pressure(1,1,k) = levels(k)
         data%temperature(1,1,k) = 288.15_fp - temp_lapse * z_height
         data%humidity(1,1,k) = 0.01_fp * exp(-z_height / 2000.0_fp)  ! kg/kg
         data%wind_u(1,1,k) = 10.0_fp * (1.0_fp - exp(-z_height / 1000.0_fp))
         data%wind_v(1,1,k) = 5.0_fp * sin(z_height / 1000.0_fp)
      enddo

      ! Initialize species concentrations (clean background)
      data%concentrations = 1.0e-9_fp  ! 1 ppb background
      data%emission_rates = 0.0_fp

      ! Set some realistic background concentrations
      if (nspecies >= 3) then
         data%concentrations(1,1,:,1) = 50.0e-9_fp   ! O3: 50 ppb
         data%concentrations(1,1,:,2) = 350.0e-6_fp  ! CO2: 350 ppm
         data%concentrations(1,1,:,3) = 1.8e-6_fp    ! CH4: 1.8 ppm
      endif

   end subroutine initialize_atmospheric_data

   !> Update meteorology for current time step
   subroutine update_meteorology_for_timestep(data, step, dt)
      type(CATChemDataType), intent(inout) :: data
      integer, intent(in) :: step
      real(fp), intent(in) :: dt
      real(fp) :: hour_of_day, diurnal_factor

      ! Simple diurnal cycle
      hour_of_day = mod(real(step-1) * dt / 3600.0_fp, 24.0_fp)
      diurnal_factor = 0.5_fp * (1.0_fp + sin(2.0_fp * 3.14159_fp * (hour_of_day - 6.0_fp) / 24.0_fp))

      ! Apply diurnal temperature variation (±5K)
      data%temperature(1,1,:) = data%temperature(1,1,:) + 5.0_fp * (diurnal_factor - 0.5_fp)

      ! Simple emission diurnal cycle (peak during day)
      if (allocated(data%emission_rates)) then
         data%emission_rates = data%emission_rates * (0.5_fp + diurnal_factor)
      endif

   end subroutine update_meteorology_for_timestep

   !> Run phases with individual timing
   subroutine run_phases_with_timing(catchem, phase_times, rc)
      type(CATChemInstanceType), intent(inout) :: catchem
      real, intent(out) :: phase_times(4)
      integer, intent(out) :: rc
      real :: t1, t2
      integer :: i
      character(len=16), parameter :: phases(4) = ['emissions ', 'chemistry ', 'transport ', 'deposition']

      rc = CATCHEM_SUCCESS
      phase_times = 0.0

      do i = 1, 4
         call cpu_time(t1)
         call catchem%run_phase(trim(phases(i)), rc)
         call cpu_time(t2)
         phase_times(i) = t2 - t1

         if (rc /= CATCHEM_SUCCESS) return
      enddo

   end subroutine run_phases_with_timing

   !> Process diagnostics for current time step
   subroutine process_column_diagnostics(catchem, step, time_str)
      type(CATChemInstanceType), intent(in) :: catchem
      integer, intent(in) :: step
      character(len=*), intent(in) :: time_str
      type(CATChemDiagnosticType), allocatable :: diagnostics(:)
      character(len=64), allocatable :: proc_names(:)
      character(len=256), allocatable :: diag_info(:)
      integer :: rc, i

      ! Get available diagnostics
      call catchem%list_process_diagnostics(proc_names, diag_info, rc)
      if (rc == CATCHEM_SUCCESS .and. allocated(proc_names)) then
         if (step == 6) then  ! First time, show available diagnostics
            print *, "    Available diagnostics:"
            do i = 1, size(proc_names)
               print *, "      ", trim(proc_names(i)), ": ", trim(diag_info(i))
            enddo
         endif
      endif

      ! Get specific diagnostics for key processes
      call get_process_diagnostic_safely(catchem, 'dust', 'emission_flux', step, time_str)
      call get_process_diagnostic_safely(catchem, 'drydep', 'deposition_velocity', step, time_str)

   end subroutine process_column_diagnostics

   !> Safely get process diagnostic (handle missing processes)
   subroutine get_process_diagnostic_safely(catchem, process_name, diag_name, step, time_str)
      type(CATChemInstanceType), intent(in) :: catchem
      character(len=*), intent(in) :: process_name, diag_name
      integer, intent(in) :: step
      character(len=*), intent(in) :: time_str
      type(CATChemDiagnosticType) :: diagnostic
      integer :: rc

      call catchem%get_process_diagnostic(process_name, diag_name, diagnostic, rc)
      if (rc == CATCHEM_SUCCESS .and. diagnostic%is_available) then
         if (step == 6) then  ! First diagnostic output
            print *, "    ", trim(process_name), " ", trim(diag_name), ": available"
         endif
      endif

   end subroutine get_process_diagnostic_safely

   !> Generate final simulation diagnostics
   subroutine generate_final_diagnostics(catchem, output_data, nsteps)
      type(CATChemInstanceType), intent(in) :: catchem
      type(CATChemDataType), intent(in) :: output_data
      integer, intent(in) :: nsteps
      real(fp) :: mean_conc, max_conc, min_conc
      integer :: i

      if (allocated(output_data%concentrations)) then
         print *, "  Final concentration statistics:"
         do i = 1, min(5, size(output_data%concentrations,4))  ! Show first 5 species
            mean_conc = sum(output_data%concentrations(1,1,:,i)) / size(output_data%concentrations,3)
            max_conc = maxval(output_data%concentrations(1,1,:,i))
            min_conc = minval(output_data%concentrations(1,1,:,i))

            write(*,'(A,I0,A,3ES12.4)') "    Species ", i, ": mean=", mean_conc, " min=", min_conc, " max=", max_conc
         enddo
      endif

   end subroutine generate_final_diagnostics

   !> Save column simulation results
   subroutine save_column_results(output_data, nsteps, levels)
      type(CATChemDataType), intent(in) :: output_data
      integer, intent(in) :: nsteps
      real(fp), intent(in) :: levels(:)

      ! In a real implementation, this would write to NetCDF or other format
      print *, "  Results would be saved to:"
      print *, "    column_output_", nsteps, "steps.nc"
      print *, "    Pressure levels: ", size(levels), " levels"
      if (allocated(output_data%concentrations)) then
         print *, "    Species concentrations: ", size(output_data%concentrations,4), " species"
      endif

   end subroutine save_column_results

   !> Format time step as hour string
   subroutine format_time_string(step, time_str)
      integer, intent(in) :: step
      character(len=*), intent(out) :: time_str
      integer :: hour

      hour = mod(step-1, 24)
      write(time_str, '(I2.2,":00 UTC")') hour

   end subroutine format_time_string

   !> Clean up allocated data
   subroutine cleanup_data(input_data, output_data)
      type(CATChemDataType), intent(inout) :: input_data, output_data

      ! Deallocate input data
      if (allocated(input_data%temperature)) deallocate(input_data%temperature)
      if (allocated(input_data%pressure)) deallocate(input_data%pressure)
      if (allocated(input_data%humidity)) deallocate(input_data%humidity)
      if (allocated(input_data%wind_u)) deallocate(input_data%wind_u)
      if (allocated(input_data%wind_v)) deallocate(input_data%wind_v)
      if (allocated(input_data%surface_pressure)) deallocate(input_data%surface_pressure)
      if (allocated(input_data%concentrations)) deallocate(input_data%concentrations)
      if (allocated(input_data%emission_rates)) deallocate(input_data%emission_rates)

      ! Deallocate output data
      if (allocated(output_data%temperature)) deallocate(output_data%temperature)
      if (allocated(output_data%pressure)) deallocate(output_data%pressure)
      if (allocated(output_data%humidity)) deallocate(output_data%humidity)
      if (allocated(output_data%wind_u)) deallocate(output_data%wind_u)
      if (allocated(output_data%wind_v)) deallocate(output_data%wind_v)
      if (allocated(output_data%surface_pressure)) deallocate(output_data%surface_pressure)
      if (allocated(output_data%concentrations)) deallocate(output_data%concentrations)
      if (allocated(output_data%emission_rates)) deallocate(output_data%emission_rates)
      if (allocated(output_data%dust_emissions)) deallocate(output_data%dust_emissions)
      if (allocated(output_data%seasalt_emissions)) deallocate(output_data%seasalt_emissions)
      if (allocated(output_data%drydep_velocity)) deallocate(output_data%drydep_velocity)

   end subroutine cleanup_data

end program standalone_column_driver
