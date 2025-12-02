!> \file host_model_integration_template.F90
!! \brief Template for integrating CATChem into a host atmospheric model
!!
!! This template provides a practical framework for integrating CATChem
!! into an existing atmospheric model. Adapt the interfaces and data
!! structures to match your specific host model.
!!
module HostCATChemInterface_Mod
   use CATChemAPI_Mod
   use precision_mod
   implicit none
   private

   public :: HostCATChemInterfaceType
   public :: HOST_CATCHEM_SUCCESS, HOST_CATCHEM_FAILURE

   ! Return codes
   integer, parameter :: HOST_CATCHEM_SUCCESS = 0
   integer, parameter :: HOST_CATCHEM_FAILURE = -1

   !> Host model state (adapt to your model's data structures)
   type :: HostModelStateType
      ! Grid dimensions
      integer :: nx, ny, nz, nspecies

      ! Meteorological fields
      real(fp), allocatable :: temperature(:,:,:)       ! [K]
      real(fp), allocatable :: pressure(:,:,:)          ! [Pa]
      real(fp), allocatable :: specific_humidity(:,:,:) ! [kg/kg]
      real(fp), allocatable :: u_wind(:,:,:)            ! [m/s]
      real(fp), allocatable :: v_wind(:,:,:)            ! [m/s]
      real(fp), allocatable :: surface_pressure(:,:)    ! [Pa]
      real(fp), allocatable :: surface_roughness(:,:)   ! [m]

      ! Chemical species (adapt names to your chemical mechanism)
      real(fp), allocatable :: species_concentrations(:,:,:,:) ! [mol/mol]
      character(len=32), allocatable :: species_names(:)

      ! Emission fields (if used)
      real(fp), allocatable :: surface_emissions(:,:,:) ! [mol/m²/s]
      real(fp), allocatable :: volume_emissions(:,:,:,:) ! [mol/m³/s]

      ! Grid coordinates
      real(fp), allocatable :: latitudes(:,:)           ! [degrees]
      real(fp), allocatable :: longitudes(:,:)          ! [degrees]
      real(fp), allocatable :: vertical_levels(:)       ! [Pa]

   end type HostModelStateType

   !> Configuration for host-CATChem integration
   type :: HostCATChemConfigType
      ! Process selection
      logical :: enable_dust_emissions = .false.
      logical :: enable_seasalt_emissions = .false.
      logical :: enable_dry_deposition = .false.
      logical :: enable_external_emissions = .false.
      logical :: enable_chemistry = .false.

      ! Performance options
      logical :: use_column_processing = .true.
      logical :: enable_diagnostics = .true.
      integer :: diagnostic_frequency = 1  ! Every N time steps

      ! Configuration files
      character(len=256) :: catchem_config_file = ''
      character(len=256) :: field_mapping_file = ''
      character(len=256) :: species_mapping_file = ''

      ! Integration options
      logical :: update_host_chemistry = .true.
      logical :: output_process_diagnostics = .false.
      real(fp) :: chemistry_timestep = 3600.0_fp  ! [seconds]

   end type HostCATChemConfigType

   !> Main interface type for host model integration
   type :: HostCATChemInterfaceType
      private

      ! CATChem components
      type(CATChemInstanceType) :: catchem
      type(FlexibleDataExchangeType) :: data_exchange

      ! Configuration
      type(HostCATChemConfigType) :: config
      logical :: is_initialized = .false.

      ! State tracking
      integer :: current_step = 0
      real(fp) :: current_time = 0.0_fp

      ! Diagnostics
      type(CATChemDiagnosticType), allocatable :: diagnostics(:)

   contains
      ! Public interface methods
      procedure :: initialize => interface_initialize
      procedure :: run_chemistry => interface_run_chemistry
      procedure :: update_emissions => interface_update_emissions
      procedure :: get_diagnostics => interface_get_diagnostics
      procedure :: finalize => interface_finalize

      ! Configuration methods
      procedure :: set_processes => interface_set_processes
      procedure :: load_field_mappings => interface_load_field_mappings
      procedure :: validate_configuration => interface_validate_configuration

      ! Utility methods
      procedure :: is_ready => interface_is_ready
      procedure :: get_species_mapping => interface_get_species_mapping
      procedure :: get_last_error => interface_get_last_error

      ! Private helper methods
      procedure, private :: map_host_to_catchem => interface_map_host_to_catchem
      procedure, private :: map_catchem_to_host => interface_map_catchem_to_host
      procedure, private :: setup_species_mapping => interface_setup_species_mapping

   end type HostCATChemInterfaceType

contains

   !> Initialize the CATChem interface
   subroutine interface_initialize(this, host_state, host_config, rc)
      class(HostCATChemInterfaceType), intent(inout) :: this
      type(HostModelStateType), intent(in) :: host_state
      type(HostCATChemConfigType), intent(in) :: host_config
      integer, intent(out) :: rc

      type(CATChemConfigType) :: catchem_config
      character(len=256) :: error_msg

      rc = HOST_CATCHEM_SUCCESS
      this%config = host_config

      print *, "Initializing CATChem interface..."

      ! ================================================================
      ! 1. SETUP CATCHEM CONFIGURATION
      ! ================================================================

      ! Map host configuration to CATChem configuration
      catchem_config%nx = host_state%nx
      catchem_config%ny = host_state%ny
      catchem_config%nz = host_state%nz
      catchem_config%nspecies = host_state%nspecies
      catchem_config%dt = host_config%chemistry_timestep

      ! Process configuration
      catchem_config%enable_dust = host_config%enable_dust_emissions
      catchem_config%enable_seasalt = host_config%enable_seasalt_emissions
      catchem_config%enable_drydep = host_config%enable_dry_deposition
      catchem_config%enable_external_emis = host_config%enable_external_emissions

      ! Performance options
      catchem_config%use_column_processing = host_config%use_column_processing
      catchem_config%enable_diagnostics = host_config%enable_diagnostics

      ! Configuration files
      catchem_config%config_file = host_config%catchem_config_file

      ! ================================================================
      ! 2. INITIALIZE CATCHEM
      ! ================================================================

      call this%catchem%init(catchem_config, rc)
      if (rc /= CATCHEM_SUCCESS) then
         call this%catchem%get_error_message(error_msg)
         print *, "ERROR: CATChem initialization failed: ", trim(error_msg)
         rc = HOST_CATCHEM_FAILURE
         return
      end if

      ! Setup grid
      call this%catchem%setup_grid(host_state%latitudes, &
         host_state%longitudes, &
         host_state%vertical_levels, rc)
      if (rc /= CATCHEM_SUCCESS) then
         print *, "ERROR: CATChem grid setup failed"
         rc = HOST_CATCHEM_FAILURE
         return
      end if

      ! ================================================================
      ! 3. SETUP FIELD MAPPINGS
      ! ================================================================

      call this%load_field_mappings(rc)
      if (rc /= HOST_CATCHEM_SUCCESS) then
         print *, "WARNING: Field mapping setup had issues"
      end if

      call this%setup_species_mapping(host_state, rc)
      if (rc /= HOST_CATCHEM_SUCCESS) then
         print *, "WARNING: Species mapping setup had issues"
      end if

      ! ================================================================
      ! 4. ADD PROCESSES
      ! ================================================================

      call this%set_processes(rc)
      if (rc /= HOST_CATCHEM_SUCCESS) then
         print *, "ERROR: Process setup failed"
         rc = HOST_CATCHEM_FAILURE
         return
      end if

      ! ================================================================
      ! 5. VALIDATE CONFIGURATION
      ! ================================================================

      call this%validate_configuration(rc)
      if (rc /= HOST_CATCHEM_SUCCESS) then
         print *, "ERROR: Configuration validation failed"
         rc = HOST_CATCHEM_FAILURE
         return
      end if

      this%is_initialized = .true.
      this%current_step = 0
      this%current_time = 0.0_fp

      print *, "CATChem interface initialized successfully"

   end subroutine interface_initialize

   !> Run chemistry for one time step
   subroutine interface_run_chemistry(this, host_state, dt, rc)
      class(HostCATChemInterfaceType), intent(inout) :: this
      type(HostModelStateType), intent(inout) :: host_state
      real(fp), intent(in) :: dt
      integer, intent(out) :: rc

      type(CATChemDataType) :: input_data, output_data
      character(len=256) :: error_msg

      rc = HOST_CATCHEM_SUCCESS

      if (.not. this%is_initialized) then
         print *, "ERROR: CATChem interface not initialized"
         rc = HOST_CATCHEM_FAILURE
         return
      end if

      this%current_step = this%current_step + 1
      this%current_time = this%current_time + dt

      ! ================================================================
      ! 1. MAP HOST DATA TO CATCHEM FORMAT
      ! ================================================================

      call this%map_host_to_catchem(host_state, input_data, rc)
      if (rc /= HOST_CATCHEM_SUCCESS) then
         print *, "ERROR: Failed to map host data to CATChem format"
         return
      end if

      ! ================================================================
      ! 2. RUN CATCHEM
      ! ================================================================

      call this%catchem%run_timestep(input_data, output_data, rc)
      if (rc /= CATCHEM_SUCCESS) then
         call this%catchem%get_error_message(error_msg)
         print *, "ERROR: CATChem time step failed: ", trim(error_msg)
         rc = HOST_CATCHEM_FAILURE
         return
      end if

      ! ================================================================
      ! 3. MAP RESULTS BACK TO HOST FORMAT
      ! ================================================================

      if (this%config%update_host_chemistry) then
         call this%map_catchem_to_host(output_data, host_state, rc)
         if (rc /= HOST_CATCHEM_SUCCESS) then
            print *, "WARNING: Failed to map CATChem results back to host"
         end if
      end if

      ! ================================================================
      ! 4. COLLECT DIAGNOSTICS
      ! ================================================================

      if (this%config%enable_diagnostics .and. &
         mod(this%current_step, this%config%diagnostic_frequency) == 0) then
         call this%get_diagnostics(rc)
         if (rc /= HOST_CATCHEM_SUCCESS) then
            print *, "WARNING: Diagnostic collection failed"
         end if
      end if

   end subroutine interface_run_chemistry

   !> Update emission fields
   subroutine interface_update_emissions(this, host_state, rc)
      class(HostCATChemInterfaceType), intent(inout) :: this
      type(HostModelStateType), intent(in) :: host_state
      integer, intent(out) :: rc

      type(CATChemDataType) :: emission_data

      rc = HOST_CATCHEM_SUCCESS

      if (.not. this%config%enable_external_emissions) return

      ! Map host emissions to CATChem format
      if (allocated(host_state%surface_emissions)) then
         allocate(emission_data%emission_rates(host_state%nx, host_state%ny, 1, host_state%nspecies))
         emission_data%emission_rates(:, :, 1, :) = host_state%surface_emissions
      end if

      ! Set emissions in CATChem
      call this%catchem%set_emissions(emission_data, rc)
      if (rc /= CATCHEM_SUCCESS) then
         print *, "WARNING: Failed to set emissions in CATChem"
         rc = HOST_CATCHEM_FAILURE
      end if

   end subroutine interface_update_emissions

   !> Get process diagnostics
   subroutine interface_get_diagnostics(this, rc)
      class(HostCATChemInterfaceType), intent(inout) :: this
      integer, intent(out) :: rc

      character(len=32), allocatable :: process_names(:)
      character(len=256), allocatable :: diagnostic_info(:)
      integer :: i

      rc = HOST_CATCHEM_SUCCESS

      if (.not. this%config%enable_diagnostics) return

      ! Get list of all available diagnostics
      call this%catchem%list_process_diagnostics(process_names, diagnostic_info, rc)
      if (rc /= CATCHEM_SUCCESS) then
         print *, "WARNING: Failed to list process diagnostics"
         return
      end if

      ! Store diagnostics for output
      if (this%config%output_process_diagnostics) then
         print *, "Step ", this%current_step, " diagnostics:"
         do i = 1, min(5, size(diagnostic_info))
            print *, "  ", trim(diagnostic_info(i))
         end do
      end if

   end subroutine interface_get_diagnostics

   !> Set up processes based on configuration
   subroutine interface_set_processes(this, rc)
      class(HostCATChemInterfaceType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = HOST_CATCHEM_SUCCESS

      ! Add dust emissions
      if (this%config%enable_dust_emissions) then
         call this%catchem%add_process('dust', rc=rc)
         if (rc == CATCHEM_SUCCESS) then
            print *, "  Dust emissions process added"
         else
            print *, "  WARNING: Failed to add dust process"
         end if
      end if

      ! Add sea salt emissions
      if (this%config%enable_seasalt_emissions) then
         call this%catchem%add_process('seasalt', rc=rc)
         if (rc == CATCHEM_SUCCESS) then
            print *, "  Sea salt emissions process added"
         else
            print *, "  WARNING: Failed to add seasalt process"
         end if
      end if

      ! Add dry deposition
      if (this%config%enable_dry_deposition) then
         call this%catchem%add_process('drydep', rc=rc)
         if (rc == CATCHEM_SUCCESS) then
            print *, "  Dry deposition process added"
         else
            print *, "  WARNING: Failed to add drydep process"
         end if
      end if

      ! Add external emissions
      if (this%config%enable_external_emissions) then
         call this%catchem%add_process('external_emissions', rc=rc)
         if (rc == CATCHEM_SUCCESS) then
            print *, "  External emissions process added"
         else
            print *, "  WARNING: Failed to add external emissions process"
         end if
      end if

      ! Check if ready
      if (.not. this%catchem%is_ready_to_run()) then
         print *, "ERROR: CATChem is not ready to run after adding processes"
         rc = HOST_CATCHEM_FAILURE
      end if

   end subroutine interface_set_processes

   !> Load field mappings from configuration file
   subroutine interface_load_field_mappings(this, rc)
      class(HostCATChemInterfaceType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = HOST_CATCHEM_SUCCESS

      ! Initialize data exchange
      call this%data_exchange%setup_flexible_exchange(this%config%nx, this%config%ny, this%config%nz, rc)
      if (rc /= CATCHEM_SUCCESS) then
         print *, "ERROR: Failed to setup flexible data exchange"
         rc = HOST_CATCHEM_FAILURE
         return
      end if

      ! Load from file if specified
      if (len_trim(this%config%field_mapping_file) > 0) then
         ! call this%catchem%setup_field_mapping(this%config%field_mapping_file, rc)
         print *, "  Field mappings loaded from: ", trim(this%config%field_mapping_file)
      else
         ! Setup default mappings
         call this%setup_default_field_mappings(rc)
      end if

   end subroutine interface_load_field_mappings

   !> Setup default field mappings
   subroutine interface_setup_default_field_mappings(this, rc)
      class(HostCATChemInterfaceType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = HOST_CATCHEM_SUCCESS

      ! Meteorological mappings
      call this%catchem%add_field_mapping('met', 'temperature', 'temperature', 'real64', rc=rc)
      call this%catchem%add_field_mapping('met', 'pressure', 'pressure', 'real64', rc=rc)
      call this%catchem%add_field_mapping('met', 'humidity', 'specific_humidity', 'real64', rc=rc)
      call this%catchem%add_field_mapping('met', 'wind_u', 'u_wind', 'real64', rc=rc)
      call this%catchem%add_field_mapping('met', 'wind_v', 'v_wind', 'real64', rc=rc)

      print *, "  Default field mappings configured"

   end subroutine interface_setup_default_field_mappings

   !> Setup species mapping between host and CATChem
   subroutine interface_setup_species_mapping(this, host_state, rc)
      class(HostCATChemInterfaceType), intent(inout) :: this
      type(HostModelStateType), intent(in) :: host_state
      integer, intent(out) :: rc

      character(len=64), allocatable :: catchem_species(:)
      integer :: i

      rc = HOST_CATCHEM_SUCCESS

      ! Get CATChem species names
      allocate(catchem_species(host_state%nspecies))
      call this%catchem%get_species_names(catchem_species, rc)
      if (rc /= CATCHEM_SUCCESS) then
         print *, "WARNING: Failed to get CATChem species names"
         return
      end if

      ! Map species (simplified - in practice, use a mapping file)
      do i = 1, min(size(host_state%species_names), size(catchem_species))
         call this%catchem%add_field_mapping('chem', trim(catchem_species(i)), &
            trim(host_state%species_names(i)), 'real64', rc=rc)
      end do

      print *, "  Species mappings configured for ", size(catchem_species), " species"

   end subroutine interface_setup_species_mapping

   !> Map host model data to CATChem format
   subroutine interface_map_host_to_catchem(this, host_state, catchem_data, rc)
      class(HostCATChemInterfaceType), intent(in) :: this
      type(HostModelStateType), intent(in) :: host_state
      type(CATChemDataType), intent(out) :: catchem_data
      integer, intent(out) :: rc

      rc = HOST_CATCHEM_SUCCESS

      ! Allocate CATChem data arrays
      allocate(catchem_data%temperature(host_state%nx, host_state%ny, host_state%nz))
      allocate(catchem_data%pressure(host_state%nx, host_state%ny, host_state%nz))
      allocate(catchem_data%humidity(host_state%nx, host_state%ny, host_state%nz))
      allocate(catchem_data%wind_u(host_state%nx, host_state%ny, host_state%nz))
      allocate(catchem_data%wind_v(host_state%nx, host_state%ny, host_state%nz))
      allocate(catchem_data%surface_pressure(host_state%nx, host_state%ny))
      allocate(catchem_data%concentrations(host_state%nx, host_state%ny, host_state%nz, host_state%nspecies))

      ! Copy data
      catchem_data%temperature = host_state%temperature
      catchem_data%pressure = host_state%pressure
      catchem_data%humidity = host_state%specific_humidity
      catchem_data%wind_u = host_state%u_wind
      catchem_data%wind_v = host_state%v_wind
      catchem_data%surface_pressure = host_state%surface_pressure
      catchem_data%concentrations = host_state%species_concentrations

      ! Add surface roughness if available
      if (allocated(host_state%surface_roughness)) then
         allocate(catchem_data%roughness_length(host_state%nx, host_state%ny))
         catchem_data%roughness_length = host_state%surface_roughness
      end if

   end subroutine interface_map_host_to_catchem

   !> Map CATChem results back to host format
   subroutine interface_map_catchem_to_host(this, catchem_data, host_state, rc)
      class(HostCATChemInterfaceType), intent(in) :: this
      type(CATChemDataType), intent(in) :: catchem_data
      type(HostModelStateType), intent(inout) :: host_state
      integer, intent(out) :: rc

      rc = HOST_CATCHEM_SUCCESS

      ! Update chemical concentrations
      if (allocated(catchem_data%concentrations)) then
         host_state%species_concentrations = catchem_data%concentrations
      end if

      ! Note: Meteorological fields typically don't change in CATChem
      ! but could be updated if needed

   end subroutine interface_map_catchem_to_host

   !> Validate the interface configuration
   subroutine interface_validate_configuration(this, rc)
      class(HostCATChemInterfaceType), intent(in) :: this
      integer, intent(out) :: rc

      rc = HOST_CATCHEM_SUCCESS

      ! Validate field mappings
      call this%catchem%validate_field_mappings(rc)
      if (rc /= CATCHEM_SUCCESS) then
         print *, "ERROR: Field mapping validation failed"
         rc = HOST_CATCHEM_FAILURE
         return
      end if

      ! Check if CATChem is ready
      if (.not. this%catchem%is_ready_to_run()) then
         print *, "ERROR: CATChem is not ready to run"
         rc = HOST_CATCHEM_FAILURE
         return
      end if

      print *, "  Configuration validation passed"

   end subroutine interface_validate_configuration

   !> Check if interface is ready
   function interface_is_ready(this) result(ready)
      class(HostCATChemInterfaceType), intent(in) :: this
      logical :: ready

      ready = this%is_initialized .and. this%catchem%is_ready_to_run()

   end function interface_is_ready

   !> Get species mapping information
   subroutine interface_get_species_mapping(this, host_species, catchem_species, rc)
      class(HostCATChemInterfaceType), intent(in) :: this
      character(len=*), intent(out) :: host_species(:)
      character(len=*), intent(out) :: catchem_species(:)
      integer, intent(out) :: rc

      rc = HOST_CATCHEM_SUCCESS

      ! Get CATChem species names
      call this%catchem%get_species_names(catchem_species, rc)

      ! Map to host species names (implementation specific)
      ! This would depend on your species mapping configuration

   end subroutine interface_get_species_mapping

   !> Get last error message
   subroutine interface_get_last_error(this, error_message)
      class(HostCATChemInterfaceType), intent(in) :: this
      character(len=*), intent(out) :: error_message

      call this%catchem%get_error_message(error_message)

   end subroutine interface_get_last_error

   !> Finalize the interface
   subroutine interface_finalize(this, rc)
      class(HostCATChemInterfaceType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = HOST_CATCHEM_SUCCESS

      ! Clean up CATChem
      call this%catchem%finalize(rc)
      if (rc /= CATCHEM_SUCCESS) then
         print *, "WARNING: CATChem finalization had issues"
         rc = HOST_CATCHEM_FAILURE
      end if

      this%is_initialized = .false.
      print *, "CATChem interface finalized"

   end subroutine interface_finalize

end module HostCATChemInterface_Mod

!> Example program demonstrating the host model integration
program host_model_integration_example
   use HostCATChemInterface_Mod
   use precision_mod
   implicit none

   type(HostCATChemInterfaceType) :: catchem_interface
   type(HostModelStateType) :: model_state
   type(HostCATChemConfigType) :: config
   integer :: rc, step
   real(fp) :: dt
   character(len=256) :: error_msg

   print *, "=== Host Model CATChem Integration Example ==="

   ! Setup host model state (simplified)
   call setup_host_model_state(model_state)

   ! Configure CATChem integration
   config%enable_dust_emissions = .true.
   config%enable_dry_deposition = .true.
   config%enable_diagnostics = .true.
   config%chemistry_timestep = 3600.0_fp
   config%diagnostic_frequency = 6

   ! Initialize CATChem interface
   call catchem_interface%initialize(model_state, config, rc)
   if (rc /= HOST_CATCHEM_SUCCESS) then
      call catchem_interface%get_last_error(error_msg)
      print *, "ERROR: ", trim(error_msg)
      stop 1
   end if

   ! Time stepping loop
   dt = 3600.0_fp  ! 1 hour
   do step = 1, 24  ! 24 hour simulation

      print *, "Step ", step, " of 24"

      ! Update host model meteorology (your model's routine)
      call update_host_meteorology(model_state, step)

      ! Run chemistry
      call catchem_interface%run_chemistry(model_state, dt, rc)
      if (rc /= HOST_CATCHEM_SUCCESS) then
         print *, "WARNING: Chemistry step failed"
      end if

      ! Continue with host model physics, dynamics, etc.

   end do

   ! Clean up
   call catchem_interface%finalize(rc)
   print *, "Integration completed"

contains

   subroutine setup_host_model_state(state)
      type(HostModelStateType), intent(out) :: state

      ! Set dimensions
      state%nx = 144
      state%ny = 91
      state%nz = 72
      state%nspecies = 10

      ! Allocate arrays
      allocate(state%temperature(state%nx, state%ny, state%nz))
      allocate(state%pressure(state%nx, state%ny, state%nz))
      allocate(state%specific_humidity(state%nx, state%ny, state%nz))
      allocate(state%u_wind(state%nx, state%ny, state%nz))
      allocate(state%v_wind(state%nx, state%ny, state%nz))
      allocate(state%surface_pressure(state%nx, state%ny))
      allocate(state%species_concentrations(state%nx, state%ny, state%nz, state%nspecies))
      allocate(state%latitudes(state%nx, state%ny))
      allocate(state%longitudes(state%nx, state%ny))
      allocate(state%vertical_levels(state%nz))
      allocate(state%species_names(state%nspecies))

      ! Initialize with example values
      state%temperature = 288.15_fp
      state%pressure = 50000.0_fp
      state%specific_humidity = 0.01_fp
      state%u_wind = 0.0_fp
      state%v_wind = 0.0_fp
      state%surface_pressure = 101325.0_fp
      state%species_concentrations = 1.0e-12_fp

      ! Species names
      state%species_names = ['O3  ', 'NO2 ', 'CO  ', 'SO2 ', 'NH3 ', &
         'PM25', 'PM10', 'BC  ', 'OC  ', 'DUST']

   end subroutine setup_host_model_state

   subroutine update_host_meteorology(state, step)
      type(HostModelStateType), intent(inout) :: state
      integer, intent(in) :: step

      ! Simple diurnal temperature cycle
      real(fp) :: time_factor
      time_factor = sin(real(step, fp) * 3.14159_fp / 12.0_fp)
      state%temperature = state%temperature + 5.0_fp * time_factor

   end subroutine update_host_meteorology

end program host_model_integration_example
