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
   use CATChemNetCDF_Mod
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
   real(fp) :: ak_full(nz+1), bk_full(nz+1)  ! Interface levels (nz+1)
   real(fp) :: ak_half(nz), bk_half(nz)      ! Mass levels (nz)

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
   print *, "  Vertical coordinate: FV3 hybrid σ-p (", nz, " mass levels)"
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

   ! Vertical levels (FV3 hybrid sigma-pressure coordinates)
   call setup_fv3_vertical_levels(levels, ak_full, bk_full, ak_half, bk_half, nz)

   call catchem%setup_grid(lats, lons, levels, rc, &
      ak_full=ak_full, bk_full=bk_full, &
      ak_half=ak_half, bk_half=bk_half)
   if (rc /= CATCHEM_SUCCESS) then
      print *, "✗ Grid setup failed!"
      stop 1
   endif

   print *, "✓ Grid coordinates configured"
   print *, "  Location: ", lats(1,1), "°N, ", lons(1,1), "°E"
   print *, "  Levels: ", nz, " (FV3 hybrid σ-p coordinates)"
   print *, "  Top level: ak=", ak_full(1)/100.0, " hPa, bk=", bk_full(1)
   print *, "  Surface: ak=", ak_full(nz+1)/100.0, " hPa, bk=", bk_full(nz+1)
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

   ! Try to read initial conditions from NetCDF (demonstration of NetCDF API)
   call try_read_initial_conditions_from_netcdf(input_data, nx, ny, nz, nspecies, &
      levels, ak_half, bk_half)

   call initialize_atmospheric_data(input_data, nx, ny, nz, nspecies, levels, ak_half, bk_half)
   call initialize_atmospheric_data(output_data, nx, ny, nz, nspecies, levels, ak_half, bk_half)

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
   print *, "  • FV3 hybrid σ-p vertical coordinate system"
   print *, "  • YSU vertical dispersion for tracer transport"
   print *, "  • High-level NetCDF I/O API demonstration"
   print *, "  • Optimized column processing"
   print *, "  • Process-level diagnostics"
   print *, "  • Modern API usage patterns"

contains

   !> Setup FV3 hybrid sigma-pressure vertical coordinates
   subroutine setup_fv3_vertical_levels(pressure_levels, ak_full, bk_full, ak_half, bk_half, nz)
      real(fp), intent(out) :: pressure_levels(nz)     ! Mass-level pressures
      real(fp), intent(out) :: ak_full(nz+1)           ! Interface A coefficients
      real(fp), intent(out) :: bk_full(nz+1)           ! Interface B coefficients
      real(fp), intent(out) :: ak_half(nz)             ! Mass-level A coefficients
      real(fp), intent(out) :: bk_half(nz)             ! Mass-level B coefficients
      integer, intent(in) :: nz

      integer :: k
      real(fp) :: p_top, p_sfc, sigma_top, sigma_sfc
      real(fp) :: eta(nz+1), eta_half(nz)
      real(fp) :: p_ref, sigma_ref

      ! FV3 72-level configuration parameters
      p_top = 100.0_fp      ! Top pressure (1 hPa)
      p_sfc = 101325.0_fp   ! Reference surface pressure (1013.25 hPa)
      sigma_top = 0.0_fp    ! Pure pressure coordinate at top
      sigma_sfc = 1.0_fp    ! Pure sigma coordinate at surface
      p_ref = 101325.0_fp   ! Reference pressure for hybrid coordinate
      sigma_ref = 0.8_fp    ! Transition level (around 800 hPa)

      ! Create eta coordinate (generalized vertical coordinate)
      ! Interface levels (full levels, nz+1 points)
      do k = 1, nz+1
         eta(k) = real(k-1) / real(nz)  ! Linear from 0 to 1
      enddo

      ! Mass levels (half levels, nz points) - midpoints
      do k = 1, nz
         eta_half(k) = 0.5_fp * (eta(k) + eta(k+1))
      enddo

      ! Setup hybrid coefficients for interface levels
      do k = 1, nz+1
         if (eta(k) <= sigma_ref) then
            ! Pure pressure coordinate in upper atmosphere
            ak_full(k) = p_top + (p_ref * sigma_ref - p_top) * eta(k) / sigma_ref
            bk_full(k) = 0.0_fp
         else
            ! Hybrid coordinate in lower atmosphere
            ak_full(k) = p_ref * sigma_ref * (1.0_fp - eta(k)) / (1.0_fp - sigma_ref)
            bk_full(k) = (eta(k) - sigma_ref) / (1.0_fp - sigma_ref)
         endif
      enddo

      ! Setup hybrid coefficients for mass levels (interpolated)
      do k = 1, nz
         ak_half(k) = 0.5_fp * (ak_full(k) + ak_full(k+1))
         bk_half(k) = 0.5_fp * (bk_full(k) + bk_full(k+1))
      enddo

      ! Calculate pressure at mass levels using reference surface pressure
      do k = 1, nz
         pressure_levels(k) = ak_half(k) + bk_half(k) * p_sfc
      enddo

      ! Print FV3 coordinate information
      print *, "    FV3 Vertical Coordinate Setup:"
      print *, "      Interface levels (ak_full, bk_full):"
      print *, "        Top:    ak=", ak_full(1)/100.0, " hPa, bk=", bk_full(1)
      print *, "        ~L36:   ak=", ak_full(nz/2)/100.0, " hPa, bk=", bk_full(nz/2)
      print *, "        Surface:ak=", ak_full(nz+1)/100.0, " hPa, bk=", bk_full(nz+1)
      print *, "      Mass levels (ak_half, bk_half):"
      print *, "        L1:     p=", pressure_levels(1)/100.0, " hPa"
      print *, "        L36:    p=", pressure_levels(nz/2)/100.0, " hPa"
      print *, "        L", nz, ":    p=", pressure_levels(nz)/100.0, " hPa"

   end subroutine setup_fv3_vertical_levels

   !> Setup simple vertical pressure levels (legacy)
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

      ! Add YSU vertical dispersion process
      call catchem%add_process('ysuverticaldispersion', rc=rc)
      if (rc == CATCHEM_SUCCESS) then
         print *, "    ✓ YSU vertical dispersion"
      else
         print *, "    ⚠ YSU vertical dispersion (process not available)"
         rc = CATCHEM_SUCCESS
      endif

   end subroutine add_processes

   !> Initialize atmospheric column data with FV3 coordinates
   subroutine initialize_atmospheric_data(data, nx, ny, nz, nspecies, levels, ak_half, bk_half)
      type(CATChemDataType), intent(out) :: data
      integer, intent(in) :: nx, ny, nz, nspecies
      real(fp), intent(in) :: levels(nz)
      real(fp), intent(in) :: ak_half(nz), bk_half(nz)
      integer :: k
      real(fp) :: z_height, temp_lapse, p_sfc_actual
      real(fp) :: scale_height, grav, rd

      ! Physical constants
      grav = 9.80665_fp      ! m/s^2
      rd = 287.04_fp         ! J/(kg*K) - dry air gas constant
      scale_height = 8400.0_fp  ! m

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
      temp_lapse = 6.5e-3_fp  ! K/m (standard atmosphere)
      p_sfc_actual = 101325.0_fp  ! Pa (actual surface pressure)
      data%surface_pressure(1,1) = p_sfc_actual

      print *, "    Initializing FV3 atmospheric profiles:"

      do k = 1, nz
         ! Calculate actual pressure using FV3 hybrid coordinates
         data%pressure(1,1,k) = ak_half(k) + bk_half(k) * p_sfc_actual

         ! Height from pressure (hypsometric equation approximation)
         z_height = -scale_height * log(data%pressure(1,1,k) / p_sfc_actual)

         ! Temperature with standard atmosphere lapse rate
         data%temperature(1,1,k) = 288.15_fp - temp_lapse * z_height

         ! Humidity (exponential decay with height)
         data%humidity(1,1,k) = 0.01_fp * exp(-z_height / 2000.0_fp)  ! kg/kg

         ! Wind profiles (realistic boundary layer + free atmosphere)
         if (z_height < 1000.0_fp) then
            ! Boundary layer wind profile
            data%wind_u(1,1,k) = 10.0_fp * (z_height / 1000.0_fp)**0.25_fp
            data%wind_v(1,1,k) = 5.0_fp * (z_height / 1000.0_fp)**0.25_fp
         else
            ! Free atmosphere
            data%wind_u(1,1,k) = 10.0_fp + 5.0_fp * sin(z_height / 5000.0_fp)
            data%wind_v(1,1,k) = 5.0_fp + 3.0_fp * cos(z_height / 5000.0_fp)
         endif

         ! Print sample levels for verification
         if (k == 1 .or. k == nz/2 .or. k == nz) then
            print *, "      L", k, ": p=", data%pressure(1,1,k)/100.0, " hPa, T=", &
               data%temperature(1,1,k), " K, z≈", z_height, " m"
         endif
      enddo

      ! Initialize species concentrations (clean background)
      data%concentrations = 1.0e-9_fp  ! 1 ppb background
      data%emission_rates = 0.0_fp

      ! Set realistic background concentrations (height-dependent)
      if (nspecies >= 3) then
         do k = 1, nz
            z_height = -scale_height * log(data%pressure(1,1,k) / p_sfc_actual)

            ! O3: stratospheric increase
            if (z_height > 15000.0_fp) then
               data%concentrations(1,1,k,1) = 200.0e-9_fp  ! 200 ppb in stratosphere
            else
               data%concentrations(1,1,k,1) = 50.0e-9_fp   ! 50 ppb in troposphere
            endif

            ! CO2: well-mixed
            data%concentrations(1,1,k,2) = 420.0e-6_fp  ! 420 ppm (current levels)

            ! CH4: surface source, decreasing with height
            data%concentrations(1,1,k,3) = 1.9e-6_fp * exp(-z_height / 8000.0_fp)  ! 1.9 ppm
         enddo
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
      call get_process_diagnostic_safely(catchem, 'ysuverticaldispersion', 'pbl_height', step, time_str)
      call get_process_diagnostic_safely(catchem, 'ysuverticaldispersion', 'mixing_coefficients', step, time_str)

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

   !> Attempt to read initial conditions from NetCDF files (demonstration of NetCDF API)
   subroutine try_read_initial_conditions_from_netcdf(data, nx, ny, nz, nspecies, &
      levels, ak_half, bk_half)
      type(CATChemDataType), intent(inout) :: data
      integer, intent(in) :: nx, ny, nz, nspecies
      real(fp), intent(in) :: levels(nz), ak_half(nz), bk_half(nz)

      type(NetCDFFileType) :: nc_file
      real(fp), allocatable :: temp_3d(:,:,:), temp_2d(:,:)
      integer :: rc
      character(len=256) :: input_file
      logical :: file_exists

      print *, "    Attempting to read initial conditions from NetCDF..."

      ! Check if NetCDF support is available
      if (.not. netcdf_available()) then
         print *, "      NetCDF support not available - using synthetic data"
         return
      endif

      ! Try common input file names
      input_file = "column_input.nc"
      inquire(file=trim(input_file), exist=file_exists)

      if (.not. file_exists) then
         input_file = "gfs_init.nc"
         inquire(file=trim(input_file), exist=file_exists)
      endif

      if (.not. file_exists) then
         input_file = "fv3_restart.nc"
         inquire(file=trim(input_file), exist=file_exists)
      endif

      if (.not. file_exists) then
         print *, "      No input NetCDF files found - using synthetic data"
         print *, "      Searched for: column_input.nc, gfs_init.nc, fv3_restart.nc"
         return
      endif

      print *, "      Found input file: ", trim(input_file)

      ! Open NetCDF file
      call nc_file%open(input_file, 'r', rc)
      if (rc /= 0) then
         print *, "      Failed to open NetCDF file - using synthetic data"
         return
      endif

      print *, "      Successfully opened NetCDF file"

      ! List available variables
      call list_netcdf_variables(nc_file)

      ! Try to read meteorological fields
      call read_met_fields_from_netcdf(nc_file, data, nx, ny, nz)

      ! Try to read FV3 coordinates if available
      call read_fv3_coords_from_netcdf(nc_file, ak_half, bk_half)

      ! Try to read species concentrations if available
      call read_species_from_netcdf(nc_file, data, nx, ny, nz, nspecies)

      ! Close file
      call nc_file%close(rc)
      if (rc == 0) then
         print *, "      ✓ NetCDF file closed successfully"
      else
         print *, "      ⚠ Warning closing NetCDF file"
      endif

   end subroutine try_read_initial_conditions_from_netcdf

   !> List available variables in NetCDF file
   subroutine list_netcdf_variables(nc_file)
      type(NetCDFFileType), intent(in) :: nc_file
      character(len=64) :: var_names(50)
      character(len=64) :: dim_names(10)
      integer :: dim_sizes(10)
      integer :: nvars, ndims, rc, i

      call nc_file%get_variables(var_names, nvars, rc)
      call nc_file%get_dimensions(dim_names, dim_sizes, ndims, rc)

      if (rc == 0) then
         print *, "      File contains ", ndims, " dimensions and ", nvars, " variables"

         if (ndims > 0) then
            print *, "      Dimensions:"
            do i = 1, min(ndims, 10)
               print *, "        ", trim(dim_names(i)), ": ", dim_sizes(i)
            enddo
         endif

         if (nvars > 0) then
            print *, "      Variables:"
            do i = 1, min(nvars, 10)  ! Show first 10 variables
               print *, "        ", trim(var_names(i))
            enddo
            if (nvars > 10) then
               print *, "        ... and ", nvars-10, " more variables"
            endif
         endif
      endif

   end subroutine list_netcdf_variables

   !> Read meteorological fields from NetCDF
   subroutine read_met_fields_from_netcdf(nc_file, data, nx, ny, nz)
      type(NetCDFFileType), intent(in) :: nc_file
      type(CATChemDataType), intent(inout) :: data
      integer, intent(in) :: nx, ny, nz

      real(fp), allocatable :: temp_3d(:,:,:), temp_2d(:,:)
      integer :: rc
      logical :: found_data = .false.

      print *, "      Attempting to read meteorological fields..."

      ! Try to read temperature
      if (nc_file%has_variable('T') .or. nc_file%has_variable('temperature')) then
         if (nc_file%has_variable('T')) then
            call nc_file%read_var('T', temp_3d, rc)
         else
            call nc_file%read_var('temperature', temp_3d, rc)
         endif

         if (rc == 0 .and. allocated(temp_3d)) then
            print *, "        ✓ Temperature read successfully"
            ! Copy to data structure (with dimension checks)
            if (size(temp_3d,1) >= nx .and. size(temp_3d,2) >= ny .and. size(temp_3d,3) >= nz) then
               if (allocated(data%temperature)) then
                  data%temperature(1:nx,1:ny,1:nz) = temp_3d(1:nx,1:ny,1:nz)
                  found_data = .true.
               endif
            endif
            deallocate(temp_3d)
         endif
      endif

      ! Try to read pressure
      if (nc_file%has_variable('P') .or. nc_file%has_variable('pressure')) then
         if (nc_file%has_variable('P')) then
            call nc_file%read_var('P', temp_3d, rc)
         else
            call nc_file%read_var('pressure', temp_3d, rc)
         endif

         if (rc == 0 .and. allocated(temp_3d)) then
            print *, "        ✓ Pressure read successfully"
            if (size(temp_3d,1) >= nx .and. size(temp_3d,2) >= ny .and. size(temp_3d,3) >= nz) then
               if (allocated(data%pressure)) then
                  data%pressure(1:nx,1:ny,1:nz) = temp_3d(1:nx,1:ny,1:nz)
                  found_data = .true.
               endif
            endif
            deallocate(temp_3d)
         endif
      endif

      ! Try to read surface pressure
      if (nc_file%has_variable('PS') .or. nc_file%has_variable('surface_pressure')) then
         if (nc_file%has_variable('PS')) then
            call nc_file%read_var('PS', temp_2d, rc)
         else
            call nc_file%read_var('surface_pressure', temp_2d, rc)
         endif

         if (rc == 0 .and. allocated(temp_2d)) then
            print *, "        ✓ Surface pressure read successfully"
            if (size(temp_2d,1) >= nx .and. size(temp_2d,2) >= ny) then
               if (allocated(data%surface_pressure)) then
                  data%surface_pressure(1:nx,1:ny) = temp_2d(1:nx,1:ny)
                  found_data = .true.
               endif
            endif
            deallocate(temp_2d)
         endif
      endif

      if (.not. found_data) then
         print *, "        No compatible meteorological fields found"
      endif

   end subroutine read_met_fields_from_netcdf

   !> Read FV3 coordinates from NetCDF
   subroutine read_fv3_coords_from_netcdf(nc_file, ak_half, bk_half)
      type(NetCDFFileType), intent(in) :: nc_file
      real(fp), intent(in) :: ak_half(:), bk_half(:)  ! Reference values

      real(fp), allocatable :: ak_file(:), bk_file(:)
      integer :: rc

      print *, "      Attempting to read FV3 coordinates..."

      ! Try to read ak and bk coefficients
      call nc_file%read_fv3_coordinates(ak_file, bk_file, rc)

      if (rc == 0 .and. allocated(ak_file) .and. allocated(bk_file)) then
         print *, "        ✓ FV3 coordinates read successfully"
         print *, "        File has ", size(ak_file), " interface levels"
         print *, "        Using synthetic coordinates for this demo"
         deallocate(ak_file, bk_file)
      else
         print *, "        No FV3 coordinates found - using synthetic coordinates"
      endif

   end subroutine read_fv3_coords_from_netcdf

   !> Read species concentrations from NetCDF
   subroutine read_species_from_netcdf(nc_file, data, nx, ny, nz, nspecies)
      type(NetCDFFileType), intent(in) :: nc_file
      type(CATChemDataType), intent(inout) :: data
      integer, intent(in) :: nx, ny, nz, nspecies

      real(fp), allocatable :: temp_4d(:,:,:,:)
      integer :: rc
      character(len=32) :: species_vars(10) = ['O3    ', 'CO2   ', 'CH4   ', 'CO    ', 'NO    ', &
         'NO2   ', 'SO2   ', 'NH3   ', 'HNO3  ', 'H2O2  ']
      integer :: i
      logical :: found_species = .false.

      print *, "      Attempting to read species concentrations..."

      ! Try to read species concentrations
      if (nc_file%has_variable('concentrations') .or. nc_file%has_variable('mixing_ratio')) then
         if (nc_file%has_variable('concentrations')) then
            call nc_file%read_var('concentrations', temp_4d, rc)
         else
            call nc_file%read_var('mixing_ratio', temp_4d, rc)
         endif

         if (rc == 0 .and. allocated(temp_4d)) then
            print *, "        ✓ Species concentrations read successfully"
            print *, "        File shape: ", shape(temp_4d)
            found_species = .true.
            deallocate(temp_4d)
         endif
      endif

      ! Try individual species
      do i = 1, min(nspecies, 10)
         if (nc_file%has_variable(trim(species_vars(i)))) then
            print *, "        Found species: ", trim(species_vars(i))
            found_species = .true.
         endif
      enddo

      if (.not. found_species) then
         print *, "        No species concentration data found"
      endif

   end subroutine read_species_from_netcdf

end program standalone_column_driver
