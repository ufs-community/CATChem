!> \file CATChemAPI_Example.F90
!! \brief Example usage of the high-level CATChem API
!! \ingroup catchem_api
!!
!! This file demonstrates how to use the high-level CATChem API for
!! easy integration with different modeling architectures.
!!
program CATChemAPI_Example
   use CATChemAPI_Mod
   use precision_mod

   implicit none

   ! =======================================================================
   ! EXAMPLE 1: Basic Usage - Simple dust and sea salt simulation
   ! =======================================================================
   call example_basic_usage()

   ! =======================================================================
   ! EXAMPLE 2: Advanced Usage - Custom processes and diagnostics
   ! =======================================================================
   call example_advanced_usage()

   ! =======================================================================
   ! EXAMPLE 3: Integration Example - WRF-Chem style integration
   ! =======================================================================
   call example_host_model_integration()

contains

   !> Example 1: Basic usage with sensible defaults
   subroutine example_basic_usage()
      type(CATChemInstanceType) :: catchem
      type(CATChemConfigType) :: config
      type(CATChemDataType) :: input_data, output_data
      integer :: rc
      character(len=256) :: error_msg

      write(*,*) '=== EXAMPLE 1: Basic CATChem Usage ==='

      ! ===================================================================
      ! Step 1: Setup configuration with sensible defaults
      ! ===================================================================
      config%nx = 10           ! 10x10 horizontal grid
      config%ny = 10
      config%nz = 20           ! 20 vertical levels
      config%nspecies = 50     ! Number of chemical species
      config%dt = 3600.0_fp    ! 1-hour time step

      ! Enable basic processes
      config%enable_dust = .true.
      config%enable_seasalt = .true.
      config%enable_drydep = .true.

      ! Use performance optimizations
      config%use_column_processing = .true.
      config%enable_diagnostics = .true.

      ! ===================================================================
      ! Step 2: Initialize CATChem - just one call!
      ! ===================================================================
      call catchem%init(config, rc)
      if (rc /= CATCHEM_SUCCESS) then
         call catchem%get_error_message(error_msg)
         write(*,*) 'ERROR: Failed to initialize CATChem: ', trim(error_msg)
         return
      endif

      write(*,*) 'CATChem initialized successfully!'

      ! ===================================================================
      ! Step 3: Setup grid coordinates (optional)
      ! ===================================================================
      block
         real(fp) :: lats(10,10), lons(10,10), levels(20)
         integer :: i, j, k

         ! Simple grid setup
         do j = 1, 10
            do i = 1, 10
               lats(i,j) = 30.0_fp + real(j-1, fp) * 0.1_fp  ! 30-31N
               lons(i,j) = -90.0_fp + real(i-1, fp) * 0.1_fp ! 90-89W
            end do
         end do

         do k = 1, 20
            levels(k) = 1000.0_fp - real(k-1, fp) * 45.0_fp  ! Pressure levels
         end do

         call catchem%setup_grid(lats, lons, levels, rc)
         if (rc /= CATCHEM_SUCCESS) then
            write(*,*) 'WARNING: Grid setup failed, using defaults'
         endif
      end block

      ! ===================================================================
      ! Step 4: Add processes - simple calls
      ! ===================================================================
      call catchem%add_process('dust', rc=rc)
      if (rc /= CATCHEM_SUCCESS) then
         write(*,*) 'ERROR: Failed to add dust process'
         return
      endif

      call catchem%add_process('seasalt', rc=rc)
      if (rc /= CATCHEM_SUCCESS) then
         write(*,*) 'ERROR: Failed to add seasalt process'
         return
      endif

      call catchem%add_process('drydep', rc=rc)
      if (rc /= CATCHEM_SUCCESS) then
         write(*,*) 'ERROR: Failed to add drydep process'
         return
      endif

      write(*,*) 'All processes added successfully!'

      ! ===================================================================
      ! Step 5: Setup input data
      ! ===================================================================
      call setup_example_input_data(input_data, config)

      ! ===================================================================
      ! Step 6: Run simulation - just one call!
      ! ===================================================================
      call catchem%run_timestep(input_data, output_data, rc)
      if (rc /= CATCHEM_SUCCESS) then
         call catchem%get_error_message(error_msg)
         write(*,*) 'ERROR: Simulation failed: ', trim(error_msg)
         return
      endif

      write(*,*) 'Simulation completed successfully!'

      ! ===================================================================
      ! Step 7: Process results
      ! ===================================================================
      if (allocated(output_data%dust_emissions)) then
         write(*,*) 'Dust emissions calculated for ', &
            size(output_data%dust_emissions, 3), ' size bins'
         write(*,*) 'Max dust emission: ', maxval(output_data%dust_emissions), ' kg/m2/s'
      endif

      if (allocated(output_data%concentrations)) then
         write(*,*) 'Species concentrations updated'
         write(*,*) 'Total species concentration: ', sum(output_data%concentrations)
      endif

      ! ===================================================================
      ! Step 8: Clean up
      ! ===================================================================
      call catchem%finalize(rc)
      write(*,*) 'CATChem finalized successfully!'

   end subroutine example_basic_usage

   !> Example 2: Advanced usage with custom configuration
   subroutine example_advanced_usage()
      type(CATChemInstanceType) :: catchem
      type(CATChemConfigType) :: config
      type(CATChemDataType) :: input_data(24), output_data(24)  ! 24-hour simulation
      integer :: rc, step
      character(len=256) :: error_msg

      write(*,*) '=== EXAMPLE 2: Advanced CATChem Usage ==='

      ! ===================================================================
      ! Advanced configuration
      ! ===================================================================
      config%nx = 50
      config%ny = 50
      config%nz = 30
      config%nspecies = 100
      config%dt = 3600.0_fp
      config%nsteps = 24

      ! Enable external emissions
      config%enable_external_emis = .true.
      config%enable_dust = .true.
      config%enable_seasalt = .true.
      config%enable_drydep = .true.

      ! High-performance options
      config%use_column_processing = .true.
      config%enable_diagnostics = .true.

      ! Configuration files
      config%config_file = 'configs/advanced_config.yml'
      config%species_file = 'configs/species_list.yml'

      ! ===================================================================
      ! Initialize with advanced configuration
      ! ===================================================================
      call catchem%init(config, rc)
      if (rc /= CATCHEM_SUCCESS) then
         call catchem%get_error_message(error_msg)
         write(*,*) 'ERROR: ', trim(error_msg)
         return
      endif

      ! ===================================================================
      ! Add processes with custom configurations
      ! ===================================================================
      call catchem%add_process('external_emission', 'configs/external_emis.yml', rc)
      call catchem%add_process('dust', 'configs/dust_fengsha.yml', rc)
      call catchem%add_process('seasalt', 'configs/seasalt_gong.yml', rc)
      call catchem%add_process('drydep', rc=rc)

      ! ===================================================================
      ! Setup time-varying input data
      ! ===================================================================
      do step = 1, 24
         call setup_hourly_input_data(input_data(step), config, step)
      end do

      ! ===================================================================
      ! Run multi-step simulation
      ! ===================================================================
      call catchem%run_multiple_steps(24, input_data, output_data, rc)
      if (rc /= CATCHEM_SUCCESS) then
         call catchem%get_error_message(error_msg)
         write(*,*) 'ERROR: Multi-step simulation failed: ', trim(error_msg)
         return
      endif

      write(*,*) '24-hour simulation completed successfully!'

      ! ===================================================================
      ! Process time series results
      ! ===================================================================
      call analyze_time_series_output(output_data)

      call catchem%finalize(rc)

   end subroutine example_advanced_usage

   !> Example 3: Integration with a host model (e.g., WRF-Chem style)
   subroutine example_host_model_integration()
      type(CATChemInstanceType) :: catchem

      write(*,*) '=== EXAMPLE 3: Host Model Integration ==='

      ! This example shows how a host model might integrate CATChem
      call host_model_main()

   end subroutine example_host_model_integration

   !> Simulated host model main routine
   subroutine host_model_main()
      type(CATChemInstanceType) :: catchem
      type(CATChemConfigType) :: catchem_config
      type(CATChemDataType) :: chem_data

      ! Host model's own data structures
      real(fp), allocatable :: host_temperature(:,:,:)
      real(fp), allocatable :: host_concentrations(:,:,:,:)
      real(fp), allocatable :: host_emissions(:,:,:,:)

      integer :: rc, time_step, max_steps
      logical :: chemistry_enabled

      ! ===================================================================
      ! Host model initialization
      ! ===================================================================
      chemistry_enabled = .true.  ! Could be from namelist
      max_steps = 100

      if (chemistry_enabled) then
         ! Setup CATChem configuration based on host model grid
         call setup_catchem_from_host_model(catchem_config)

         ! Initialize CATChem
         call catchem%init(catchem_config, rc)
         if (rc /= CATCHEM_SUCCESS) then
            write(*,*) 'ERROR: CATChem initialization failed'
            chemistry_enabled = .false.
         else
            ! Add required processes
            call catchem%add_process('dust', rc=rc)
            call catchem%add_process('seasalt', rc=rc)
            call catchem%add_process('external_emission', rc=rc)
            call catchem%add_process('drydep', rc=rc)
         endif
      endif

      ! ===================================================================
      ! Host model time loop
      ! ===================================================================
      do time_step = 1, max_steps

         ! Host model physics (dynamics, microphysics, etc.)
         call host_model_physics(time_step)

         ! Run chemistry if enabled
         if (chemistry_enabled) then
            ! Transfer data from host model to CATChem
            call transfer_host_to_catchem(host_temperature, host_concentrations, &
               host_emissions, chem_data)

            ! Run CATChem for this time step
            call catchem%run_timestep(chem_data, chem_data, rc)
            if (rc /= CATCHEM_SUCCESS) then
               write(*,*) 'WARNING: Chemistry failed at step ', time_step
               cycle
            endif

            ! Transfer results back to host model
            call transfer_catchem_to_host(chem_data, host_concentrations)
         endif

         ! Host model output and diagnostics
         if (mod(time_step, 10) == 0) then
            call host_model_output(time_step, host_concentrations)
         endif

      end do

      ! ===================================================================
      ! Cleanup
      ! ===================================================================
      if (chemistry_enabled) then
         call catchem%finalize(rc)
      endif

      write(*,*) 'Host model integration example completed!'

   end subroutine host_model_main

   ! =======================================================================
   ! Helper routines (simplified implementations)
   ! =======================================================================

   subroutine setup_example_input_data(data, config)
      type(CATChemDataType), intent(out) :: data
      type(CATChemConfigType), intent(in) :: config

      integer :: i, j, k, n

      ! Allocate arrays
      allocate(data%temperature(config%nx, config%ny, config%nz))
      allocate(data%pressure(config%nx, config%ny, config%nz))
      allocate(data%humidity(config%nx, config%ny, config%nz))
      allocate(data%wind_u(config%nx, config%ny, config%nz))
      allocate(data%wind_v(config%nx, config%ny, config%nz))
      allocate(data%concentrations(config%nx, config%ny, config%nz, config%nspecies))

      ! Fill with example data
      do k = 1, config%nz
         do j = 1, config%ny
            do i = 1, config%nx
               data%temperature(i,j,k) = 288.0_fp + real(k-1, fp) * (-6.5e-3_fp * 100.0_fp)
               data%pressure(i,j,k) = 1013.25e2_fp * (data%temperature(i,j,k)/288.0_fp)**5.26_fp
               data%humidity(i,j,k) = 0.01_fp
               data%wind_u(i,j,k) = 5.0_fp + real(i-1, fp) * 0.1_fp
               data%wind_v(i,j,k) = 2.0_fp + real(j-1, fp) * 0.1_fp

               do n = 1, config%nspecies
                  data%concentrations(i,j,k,n) = 1.0e-9_fp * real(n, fp)
               end do
            end do
         end do
      end do

   end subroutine setup_example_input_data

   subroutine setup_hourly_input_data(data, config, hour)
      type(CATChemDataType), intent(out) :: data
      type(CATChemConfigType), intent(in) :: config
      integer, intent(in) :: hour

      ! Setup base data
      call setup_example_input_data(data, config)

      ! Add hourly variations (simplified)
      data%temperature = data%temperature + 5.0_fp * sin(real(hour, fp) * 3.14159_fp / 12.0_fp)
      data%wind_u = data%wind_u * (1.0_fp + 0.2_fp * sin(real(hour, fp) * 3.14159_fp / 12.0_fp))

   end subroutine setup_hourly_input_data

   subroutine analyze_time_series_output(output_data)
      type(CATChemDataType), intent(in) :: output_data(:)
      integer :: step

      write(*,*) 'Time series analysis:'
      do step = 1, size(output_data)
         if (allocated(output_data(step)%dust_emissions)) then
            write(*,'(A,I3,A,E12.4)') 'Hour ', step, ': Max dust emission = ', &
               maxval(output_data(step)%dust_emissions)
         endif
      end do

   end subroutine analyze_time_series_output

   subroutine setup_catchem_from_host_model(config)
      type(CATChemConfigType), intent(out) :: config

      ! Host model would set these based on its own configuration
      config%nx = 100
      config%ny = 100
      config%nz = 40
      config%nspecies = 80
      config%dt = 300.0_fp  ! 5-minute time step

      config%enable_dust = .true.
      config%enable_seasalt = .true.
      config%enable_external_emis = .true.
      config%enable_drydep = .true.

      config%use_column_processing = .true.
      config%enable_diagnostics = .true.

   end subroutine setup_catchem_from_host_model

   subroutine host_model_physics(time_step)
      integer, intent(in) :: time_step
      ! Placeholder for host model physics
   end subroutine host_model_physics

   subroutine transfer_host_to_catchem(host_temp, host_conc, host_emis, chem_data)
      real(fp), intent(in) :: host_temp(:,:,:)
      real(fp), intent(in) :: host_conc(:,:,:,:)
      real(fp), intent(in) :: host_emis(:,:,:,:)
      type(CATChemDataType), intent(out) :: chem_data

      ! Transfer data from host model arrays to CATChem data structure
      ! This would include unit conversions, array copying, etc.
      ! Simplified implementation here

   end subroutine transfer_host_to_catchem

   subroutine transfer_catchem_to_host(chem_data, host_conc)
      type(CATChemDataType), intent(in) :: chem_data
      real(fp), intent(out) :: host_conc(:,:,:,:)

      ! Transfer results back to host model
      if (allocated(chem_data%concentrations)) then
         host_conc = chem_data%concentrations
      endif

   end subroutine transfer_catchem_to_host

   subroutine host_model_output(time_step, concentrations)
      integer, intent(in) :: time_step
      real(fp), intent(in) :: concentrations(:,:,:,:)

      write(*,'(A,I6,A,E12.4)') 'Step ', time_step, ': Total concentration = ', &
         sum(concentrations)

   end subroutine host_model_output

end program CATChemAPI_Example
