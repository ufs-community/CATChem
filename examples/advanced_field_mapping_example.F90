!> \file advanced_field_mapping_example.F90
!! \brief Advanced example demonstrating flexible field mapping
!!
!! This example shows how to use the flexible field mapping system
!! for bidirectional data exchange between a host model and CATChem.
!!
program advanced_field_mapping_example
   use CATChemAPI_Mod
   use catchem_bridge_precision
   implicit none

   ! CATChem components
   type(CATChemInstanceType) :: catchem
   type(CATChemConfigType) :: config
   type(FlexibleDataExchangeType) :: data_exchange

   ! Host model data (simulated)
   type :: HostModelDataType
      real(fp), allocatable :: temp_3d(:,:,:)           ! Temperature [K]
      real(fp), allocatable :: pres_3d(:,:,:)           ! Pressure [hPa]
      real(fp), allocatable :: humidity_3d(:,:,:)       ! Relative humidity [%]
      real(fp), allocatable :: ozone_vmr(:,:,:)         ! Ozone volume mixing ratio [ppbv]
      real(fp), allocatable :: no2_vmr(:,:,:)           ! NO2 volume mixing ratio [ppbv]
      real(fp), allocatable :: dust_emissions_2d(:,:)   ! Dust emissions [μg/m²/s]
      real(fp), allocatable :: surface_pres(:,:)        ! Surface pressure [hPa]
   end type HostModelDataType

   type(HostModelDataType) :: host_data

   ! Grid and control
   integer, parameter :: nx = 72, ny = 46, nz = 25, nspecies = 5
   real(fp), allocatable :: lats(:,:), lons(:,:), levels(:)
   integer :: rc, step
   character(len=256) :: error_msg

   print *, "=== CATChem Advanced Field Mapping Example ==="

   ! ====================================================================
   ! 1. INITIALIZE CATCHEM
   ! ====================================================================

   print *, "Initializing CATChem..."

   config%nx = nx
   config%ny = ny
   config%nz = nz
   config%nspecies = nspecies
   config%dt = 3600.0_fp
   config%enable_dust = .true.
   config%enable_drydep = .true.
   config%enable_diagnostics = .true.

   call catchem%init(config, rc)
   if (rc /= CATCHEM_SUCCESS) then
      call catchem%get_error_message(error_msg)
      print *, "ERROR: ", trim(error_msg)
      stop 1
   end if

   ! Setup grid
   call setup_grid(lats, lons, levels)
   call catchem%setup_grid(lats, lons, levels, rc)

   ! Add processes
   call catchem%add_process('dust', rc=rc)
   call catchem%add_process('drydep', rc=rc)

   ! ====================================================================
   ! 2. SETUP FIELD MAPPINGS
   ! ====================================================================

   print *, "Setting up field mappings..."

   ! Initialize flexible data exchange
   call data_exchange%setup_flexible_exchange(nx, ny, nz, rc)
   if (rc /= CATCHEM_SUCCESS) then
      print *, "ERROR: Failed to setup flexible data exchange"
      stop 1
   end if

   ! Meteorological field mappings
   call catchem%add_field_mapping('met', 'temperature', 'temp_3d', 'real64', rc=rc)
   if (rc /= CATCHEM_SUCCESS) print *, "WARNING: Temperature mapping failed"

   ! Pressure mapping with unit conversion (hPa to Pa)
   call catchem%add_field_mapping('met', 'pressure', 'pres_3d', 'real64', 100.0_fp, rc)
   if (rc /= CATCHEM_SUCCESS) print *, "WARNING: Pressure mapping failed"

   ! Humidity mapping with unit conversion (% to kg/kg)
   call catchem%add_field_mapping('met', 'humidity', 'humidity_3d', 'real64', 0.01_fp, rc)
   if (rc /= CATCHEM_SUCCESS) print *, "WARNING: Humidity mapping failed"

   ! Chemical species mappings with unit conversion (ppbv to mol/mol)
   call catchem%add_field_mapping('chem', 'O3', 'ozone_vmr', 'real64', 1.0e-9_fp, rc)
   call catchem%add_field_mapping('chem', 'NO2', 'no2_vmr', 'real64', 1.0e-9_fp, rc)

   ! Diagnostic mappings
   call catchem%add_field_mapping('diag', 'dust_emission_flux', 'dust_emissions_2d', 'real64', rc=rc)

   ! Validate all mappings
   call catchem%validate_field_mappings(rc)
   if (rc /= CATCHEM_SUCCESS) then
      print *, "WARNING: Field mapping validation had issues"
   else
      print *, "  All field mappings validated successfully"
   end if

   ! ====================================================================
   ! 3. PREPARE HOST MODEL DATA
   ! ====================================================================

   print *, "Preparing host model data..."

   call allocate_host_data(host_data)
   call initialize_host_data(host_data)

   ! Map host data to flexible exchange format
   call map_host_to_exchange(host_data, data_exchange)

   ! ====================================================================
   ! 4. BIDIRECTIONAL DATA EXCHANGE DEMO
   ! ====================================================================

   print *, "Demonstrating bidirectional data exchange..."

   do step = 1, 5

      print *, "  Step ", step

      ! 4a. Fill CATChem states from host data
      print *, "    Filling CATChem states from host data..."
      call catchem%fill_met_state_from_host(data_exchange, rc)
      if (rc /= CATCHEM_SUCCESS) then
         print *, "    WARNING: Met state filling had issues"
      end if

      call catchem%fill_chem_state_from_host(data_exchange, rc)
      if (rc /= CATCHEM_SUCCESS) then
         print *, "    WARNING: Chem state filling had issues"
      end if

      ! 4b. Run CATChem processes (simplified)
      print *, "    Running CATChem processes..."
      ! In a real application, you would use the full run_timestep method
      ! Here we just simulate the process execution

      ! 4c. Extract results back to host format
      print *, "    Extracting results to host format..."
      call catchem%extract_chem_state_to_host(data_exchange, rc)
      if (rc /= CATCHEM_SUCCESS) then
         print *, "    WARNING: Chem state extraction had issues"
      end if

      call catchem%extract_diagnostics_to_host(data_exchange, rc)
      if (rc /= CATCHEM_SUCCESS) then
         print *, "    WARNING: Diagnostic extraction had issues"
      end if

      ! 4d. Map results back to host data structures
      call map_exchange_to_host(data_exchange, host_data)

      ! 4e. Process results
      call process_host_results(host_data, step)

   end do

   ! ====================================================================
   ! 5. DEMONSTRATE FIELD MAPPING FLEXIBILITY
   ! ====================================================================

   print *, "Demonstrating field mapping flexibility..."

   ! Show mapping details
   call show_mapping_details(data_exchange)

   ! Demonstrate unit conversions
   call demonstrate_unit_conversions()

   ! Show diagnostic field access
   call show_diagnostic_access(catchem)

   ! ====================================================================
   ! 6. CLEANUP
   ! ====================================================================

   print *, "Cleaning up..."
   call catchem%finalize(rc)

   print *, "=== Advanced field mapping example completed ==="

contains

   !> Setup example grid
   subroutine setup_grid(lats, lons, levels)
      real(fp), allocatable, intent(out) :: lats(:,:), lons(:,:), levels(:)
      integer :: i, j, k

      allocate(lats(nx, ny), lons(nx, ny), levels(nz))

      do j = 1, ny
         do i = 1, nx
            lons(i, j) = real(i-1, fp) * 360.0_fp / real(nx, fp)
            lats(i, j) = -90.0_fp + real(j-1, fp) * 180.0_fp / real(ny-1, fp)
         end do
      end do

      do k = 1, nz
         levels(k) = 1000.0_fp * (1.0_fp - real(k-1, fp) / real(nz-1, fp))
      end do

   end subroutine setup_grid

   !> Allocate host model data arrays
   subroutine allocate_host_data(data)
      type(HostModelDataType), intent(inout) :: data

      allocate(data%temp_3d(nx, ny, nz))
      allocate(data%pres_3d(nx, ny, nz))
      allocate(data%humidity_3d(nx, ny, nz))
      allocate(data%ozone_vmr(nx, ny, nz))
      allocate(data%no2_vmr(nx, ny, nz))
      allocate(data%dust_emissions_2d(nx, ny))
      allocate(data%surface_pres(nx, ny))

   end subroutine allocate_host_data

   !> Initialize host model data with example values
   subroutine initialize_host_data(data)
      type(HostModelDataType), intent(inout) :: data
      integer :: i, j, k

      ! Temperature in Kelvin
      data%temp_3d = 288.15_fp

      ! Pressure in hPa (will be converted to Pa)
      do k = 1, nz
         data%pres_3d(:, :, k) = 1013.25_fp * (1.0_fp - real(k-1, fp) / real(nz, fp))
      end do

      ! Relative humidity in %
      data%humidity_3d = 60.0_fp

      ! Ozone in ppbv (will be converted to mol/mol)
      data%ozone_vmr = 50.0_fp

      ! NO2 in ppbv
      data%no2_vmr = 10.0_fp

      ! Initialize emissions to zero
      data%dust_emissions_2d = 0.0_fp

      ! Surface pressure in hPa
      data%surface_pres = 1013.25_fp

   end subroutine initialize_host_data

   !> Map host data to flexible exchange format
   subroutine map_host_to_exchange(host_data, exchange)
      type(HostModelDataType), intent(in) :: host_data
      type(FlexibleDataExchangeType), intent(inout) :: exchange

      ! This is a simplified mapping - in practice, the flexible exchange
      ! system would handle this automatically based on the field mappings

      print *, "    Mapping host data to exchange format..."
      print *, "      Temperature range: ", minval(host_data%temp_3d), " to ", maxval(host_data%temp_3d), " K"
      print *, "      Pressure range: ", minval(host_data%pres_3d), " to ", maxval(host_data%pres_3d), " hPa"
      print *, "      Ozone range: ", minval(host_data%ozone_vmr), " to ", maxval(host_data%ozone_vmr), " ppbv"

   end subroutine map_host_to_exchange

   !> Map exchange format back to host data
   subroutine map_exchange_to_host(exchange, host_data)
      type(FlexibleDataExchangeType), intent(in) :: exchange
      type(HostModelDataType), intent(inout) :: host_data

      print *, "    Mapping exchange format back to host data..."
      ! This would extract the processed data back to host format

   end subroutine map_exchange_to_host

   !> Process host model results
   subroutine process_host_results(data, step)
      type(HostModelDataType), intent(in) :: data
      integer, intent(in) :: step

      print *, "    Processing results for step ", step
      print *, "      Updated ozone range: ", minval(data%ozone_vmr), " to ", maxval(data%ozone_vmr), " ppbv"
      print *, "      Dust emissions range: ", minval(data%dust_emissions_2d), " to ", maxval(data%dust_emissions_2d), " μg/m²/s"

   end subroutine process_host_results

   !> Show detailed mapping information
   subroutine show_mapping_details(exchange)
      type(FlexibleDataExchangeType), intent(in) :: exchange

      print *, "  Field mapping details:"
      print *, "    Meteorological mappings: ", exchange%field_mapping_registry%num_met_mappings
      print *, "    Chemical mappings: ", exchange%field_mapping_registry%num_chem_mappings
      print *, "    Emission mappings: ", exchange%field_mapping_registry%num_emis_mappings
      print *, "    Diagnostic mappings: ", exchange%field_mapping_registry%num_diag_mappings

   end subroutine show_mapping_details

   !> Demonstrate unit conversion capabilities
   subroutine demonstrate_unit_conversions()
      real(fp) :: pressure_hpa, pressure_pa
      real(fp) :: humidity_percent, humidity_kgkg
      real(fp) :: ozone_ppbv, ozone_molmol

      print *, "  Unit conversion examples:"

      ! Pressure conversion
      pressure_hpa = 1013.25_fp
      pressure_pa = pressure_hpa * 100.0_fp
      print *, "    Pressure: ", pressure_hpa, " hPa = ", pressure_pa, " Pa"

      ! Humidity conversion
      humidity_percent = 60.0_fp
      humidity_kgkg = humidity_percent * 0.01_fp  ! Simplified conversion
      print *, "    Humidity: ", humidity_percent, " % ≈ ", humidity_kgkg, " kg/kg (simplified)"

      ! Ozone conversion
      ozone_ppbv = 50.0_fp
      ozone_molmol = ozone_ppbv * 1.0e-9_fp
      print *, "    Ozone: ", ozone_ppbv, " ppbv = ", ozone_molmol, " mol/mol"

   end subroutine demonstrate_unit_conversions

   !> Show diagnostic field access
   subroutine show_diagnostic_access(catchem)
      type(CATChemInstanceType), intent(inout) :: catchem

      character(len=64), allocatable :: diagnostic_names(:)
      type(CATChemDiagnosticType) :: diagnostic
      integer :: rc, i

      print *, "  Diagnostic field access:"

      ! Get dust diagnostics
      call catchem%get_available_diagnostics('dust', diagnostic_names, rc)
      if (rc == CATCHEM_SUCCESS) then
         print *, "    Dust process diagnostics: ", size(diagnostic_names)
         do i = 1, min(3, size(diagnostic_names))
            print *, "      ", trim(diagnostic_names(i))
         end do
      end if

      ! Get dry deposition diagnostics
      call catchem%get_available_diagnostics('drydep', diagnostic_names, rc)
      if (rc == CATCHEM_SUCCESS) then
         print *, "    Dry deposition diagnostics: ", size(diagnostic_names)
         do i = 1, min(3, size(diagnostic_names))
            print *, "      ", trim(diagnostic_names(i))
         end do
      end if

   end subroutine show_diagnostic_access

end program advanced_field_mapping_example
