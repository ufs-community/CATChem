!> \file simple_integration_example.F90
!! \brief Simple example of CATChem integration
!!
!! This example demonstrates the basic usage of the CATChem high-level API
!! for a simple atmospheric chemistry simulation.
!!
program simple_integration_example
   use CATChemAPI_Mod
   use catchem_bridge_precision
   implicit none

   ! CATChem components
   type(CATChemInstanceType) :: catchem
   type(CATChemConfigType) :: config
   type(CATChemDataType) :: input_data, output_data

   ! Grid and data
   integer, parameter :: nx = 144, ny = 91, nz = 72, nspecies = 10
   real(fp), allocatable :: lats(:,:), lons(:,:), levels(:)
   real(fp), allocatable :: temperature(:,:,:), pressure(:,:,:)
   real(fp), allocatable :: concentrations(:,:,:,:)

   ! Control variables
   integer :: rc, step, nsteps
   character(len=256) :: error_msg

   print *, "=== CATChem Simple Integration Example ==="

   ! ====================================================================
   ! 1. CONFIGURATION
   ! ====================================================================

   print *, "Setting up configuration..."

   ! Basic grid configuration
   config%nx = nx
   config%ny = ny
   config%nz = nz
   config%nspecies = nspecies

   ! Time stepping
   config%dt = 3600.0_fp    ! 1 hour time step
   nsteps = 24              ! 24 hour simulation

   ! Enable desired processes
   config%enable_dust = .true.
   config%enable_seasalt = .true.
   config%enable_drydep = .true.
   config%enable_external_emis = .false.

   ! Performance and diagnostics
   config%use_column_processing = .true.
   config%enable_diagnostics = .true.

   ! Configuration files (optional)
   config%config_file = 'examples/simple_config.yml'
   config%species_file = 'examples/species_list.yml'

   ! ====================================================================
   ! 2. INITIALIZATION
   ! ====================================================================

   print *, "Initializing CATChem..."

   call catchem%init(config, rc)
   if (rc /= CATCHEM_SUCCESS) then
      call catchem%get_error_message(error_msg)
      print *, "ERROR: Failed to initialize CATChem: ", trim(error_msg)
      stop 1
   end if

   ! Setup grid coordinates
   call setup_example_grid(lats, lons, levels)
   call catchem%setup_grid(lats, lons, levels, rc)
   if (rc /= CATCHEM_SUCCESS) then
      print *, "ERROR: Failed to setup grid"
      stop 1
   end if

   ! ====================================================================
   ! 3. ADD PROCESSES
   ! ====================================================================

   print *, "Adding processes..."

   ! Add dust emissions process
   call catchem%add_process('dust', rc=rc)
   if (rc /= CATCHEM_SUCCESS) then
      print *, "WARNING: Failed to add dust process"
   else
      print *, "  - Dust emissions enabled"
   end if

   ! Add sea salt emissions process
   call catchem%add_process('seasalt', rc=rc)
   if (rc /= CATCHEM_SUCCESS) then
      print *, "WARNING: Failed to add seasalt process"
   else
      print *, "  - Sea salt emissions enabled"
   end if

   ! Add dry deposition process
   call catchem%add_process('drydep', rc=rc)
   if (rc /= CATCHEM_SUCCESS) then
      print *, "WARNING: Failed to add drydep process"
   else
      print *, "  - Dry deposition enabled"
   end if

   ! Check if ready to run
   if (.not. catchem%is_ready_to_run()) then
      print *, "ERROR: CATChem is not ready to run"
      stop 1
   end if

   print *, "CATChem successfully initialized and ready to run"

   ! ====================================================================
   ! 4. PREPARE INPUT DATA
   ! ====================================================================

   print *, "Preparing input data..."

   ! Allocate input data arrays
   allocate(input_data%temperature(nx, ny, nz))
   allocate(input_data%pressure(nx, ny, nz))
   allocate(input_data%humidity(nx, ny, nz))
   allocate(input_data%wind_u(nx, ny, nz))
   allocate(input_data%wind_v(nx, ny, nz))
   allocate(input_data%surface_pressure(nx, ny))
   allocate(input_data%concentrations(nx, ny, nz, nspecies))

   ! Initialize with example data
   call setup_example_meteorology(input_data)
   call setup_example_chemistry(input_data)

   ! Validate input data
   call catchem%validate_data(input_data, rc)
   if (rc /= CATCHEM_SUCCESS) then
      print *, "ERROR: Input data validation failed"
      stop 1
   end if

   ! ====================================================================
   ! 5. TIME STEPPING LOOP
   ! ====================================================================

   print *, "Starting time integration..."
   print *, "Number of time steps: ", nsteps
   print *, "Time step: ", config%dt, " seconds"

   do step = 1, nsteps

      print *, "  Step ", step, " of ", nsteps

      ! Update meteorology for this time step (if needed)
      call update_meteorology_for_step(input_data, step)

      ! Run CATChem for this time step
      call catchem%run_timestep(input_data, output_data, rc)
      if (rc /= CATCHEM_SUCCESS) then
         call catchem%get_error_message(error_msg)
         print *, "ERROR: Time step failed: ", trim(error_msg)
         stop 1
      end if

      ! Update concentrations for next step
      if (allocated(output_data%concentrations)) then
         input_data%concentrations = output_data%concentrations
      end if

      ! Process diagnostics (every 6 hours)
      if (mod(step, 6) == 0) then
         call process_diagnostics(catchem, step)
      end if

      ! Optional: Save output data
      if (mod(step, 12) == 0) then
         call save_output_data(output_data, step)
      end if

   end do

   ! ====================================================================
   ! 6. FINAL DIAGNOSTICS AND CLEANUP
   ! ====================================================================

   print *, "Processing final diagnostics..."
   call process_final_diagnostics(catchem)

   print *, "Finalizing CATChem..."
   call catchem%finalize(rc)
   if (rc /= CATCHEM_SUCCESS) then
      print *, "WARNING: Finalization had issues"
   end if

   print *, "=== Simulation completed successfully ==="

contains

   !> Setup example grid coordinates
   subroutine setup_example_grid(lats, lons, levels)
      real(fp), allocatable, intent(out) :: lats(:,:), lons(:,:), levels(:)
      integer :: i, j, k

      allocate(lats(nx, ny), lons(nx, ny), levels(nz))

      ! Simple regular lat-lon grid
      do j = 1, ny
         do i = 1, nx
            lons(i, j) = -180.0_fp + (i-1) * 360.0_fp / real(nx-1, fp)
            lats(i, j) = -90.0_fp + (j-1) * 180.0_fp / real(ny-1, fp)
         end do
      end do

      ! Simple pressure levels (Pa)
      do k = 1, nz
         levels(k) = 1000.0_fp * (1.0_fp - real(k-1, fp) / real(nz-1, fp))
      end do

   end subroutine setup_example_grid

   !> Setup example meteorological data
   subroutine setup_example_meteorology(data)
      type(CATChemDataType), intent(inout) :: data
      integer :: i, j, k
      real(fp) :: lat, lon, pressure_hpa, lapse_rate

      lapse_rate = 6.5e-3_fp  ! K/m

      do k = 1, nz
         pressure_hpa = levels(k) / 100.0_fp  ! Convert Pa to hPa
         do j = 1, ny
            do i = 1, nx
               lat = lats(i, j)
               lon = lons(i, j)

               ! Temperature with simple lapse rate
               data%temperature(i, j, k) = 288.15_fp - lapse_rate * (1000.0_fp - pressure_hpa) * 10.0_fp

               ! Pressure
               data%pressure(i, j, k) = levels(k)

               ! Humidity (simple pattern)
               data%humidity(i, j, k) = 0.01_fp * exp(-pressure_hpa/500.0_fp)

               ! Winds (simple pattern)
               data%wind_u(i, j, k) = 10.0_fp * sin(lat * 3.14159_fp/180.0_fp)
               data%wind_v(i, j, k) = 5.0_fp * cos(lon * 3.14159_fp/180.0_fp)

            end do
         end do
      end do

      ! Surface pressure
      data%surface_pressure = levels(nz)

   end subroutine setup_example_meteorology

   !> Setup example chemical initial conditions
   subroutine setup_example_chemistry(data)
      type(CATChemDataType), intent(inout) :: data
      integer :: i, j, k, s

      ! Initialize all species to small background values
      data%concentrations = 1.0e-12_fp

      ! Set some realistic background concentrations (mol/mol)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               ! Example: O3 profile
               if (nspecies >= 1) then
                  data%concentrations(i, j, k, 1) = 30.0e-9_fp * exp(-levels(k)/20000.0_fp)
               end if

               ! Example: CO background
               if (nspecies >= 2) then
                  data%concentrations(i, j, k, 2) = 100.0e-9_fp
               end if

               ! Add other species as needed...
            end do
         end do
      end do

   end subroutine setup_example_chemistry

   !> Update meteorology for each time step
   subroutine update_meteorology_for_step(data, step)
      type(CATChemDataType), intent(inout) :: data
      integer, intent(in) :: step
      real(fp) :: time_factor

      ! Simple diurnal cycle
      time_factor = sin(real(step, fp) * 3.14159_fp / 12.0_fp)

      ! Adjust temperature with diurnal cycle (±5K)
      data%temperature = data%temperature + 5.0_fp * time_factor

   end subroutine update_meteorology_for_step

   !> Process and display diagnostics
   subroutine process_diagnostics(catchem, step)
      type(CATChemInstanceType), intent(inout) :: catchem
      integer, intent(in) :: step

      character(len=64), allocatable :: diagnostic_names(:)
      type(CATChemDiagnosticType) :: diagnostic
      integer :: rc, i

      print *, "    Diagnostics for step ", step, ":"

      ! Get dust diagnostics
      call catchem%get_available_diagnostics('dust', diagnostic_names, rc)
      if (rc == CATCHEM_SUCCESS .and. size(diagnostic_names) > 0) then
         call catchem%get_process_diagnostic('dust', diagnostic_names(1), diagnostic, rc)
         if (diagnostic%is_available) then
            print *, "      Dust diagnostic available: ", trim(diagnostic%field_name)
         end if
      end if

      ! Get dry deposition diagnostics
      call catchem%get_available_diagnostics('drydep', diagnostic_names, rc)
      if (rc == CATCHEM_SUCCESS .and. size(diagnostic_names) > 0) then
         print *, "      Dry deposition diagnostics: ", size(diagnostic_names), " fields"
      end if

   end subroutine process_diagnostics

   !> Process final diagnostics
   subroutine process_final_diagnostics(catchem)
      type(CATChemInstanceType), intent(inout) :: catchem

      character(len=32), allocatable :: process_names(:)
      character(len=256), allocatable :: diagnostic_info(:)
      integer :: rc, i

      ! List all available diagnostics
      call catchem%list_process_diagnostics(process_names, diagnostic_info, rc)
      if (rc == CATCHEM_SUCCESS) then
         print *, "Final diagnostic summary:"
         print *, "  Processes with diagnostics: ", size(process_names)
         print *, "  Total diagnostic fields: ", size(diagnostic_info)

         do i = 1, min(5, size(diagnostic_info))  ! Show first 5
            print *, "    ", trim(diagnostic_info(i))
         end do

         if (size(diagnostic_info) > 5) then
            print *, "    ... and ", size(diagnostic_info) - 5, " more"
         end if
      end if

   end subroutine process_final_diagnostics

   !> Save output data (placeholder)
   subroutine save_output_data(data, step)
      type(CATChemDataType), intent(in) :: data
      integer, intent(in) :: step

      print *, "      Saving output for step ", step, " (placeholder)"
      ! In a real application, you would save to NetCDF, HDF5, etc.

   end subroutine save_output_data

end program simple_integration_example
