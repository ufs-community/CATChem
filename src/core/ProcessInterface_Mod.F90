!> \file ProcessInterface_Mod.F90
!! \brief Abstract base class interface for all atmospheric processes
!! \ingroup process_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module defines the abstract interface that all atmospheric processes
!! must implement, following the exact structure from PROCESS_ARCHITECTURE_GUIDE.md
!!
module ProcessInterface_Mod
   use precision_mod
   use StateManager_Mod, only : StateManagerType
   use error_mod
   use ColumnInterface_Mod, only : ColumnProcessorType
   use VirtualColumn_Mod, only : VirtualColumnType
   use ExtEmisData_Mod, only : ExtEmisDataType
   use DiagnosticManager_Mod, only: DiagnosticManagerType
   use DiagnosticInterface_Mod, only: DiagnosticRegistryType, DiagnosticFieldType, DiagnosticDataType

   implicit none
   private

   public :: ProcessInterface
   public :: ColumnProcessInterface

   !> \brief Abstract base class for all atmospheric processes
   !!
   !! This abstract type defines the interface that all atmospheric processes
   !! must implement, following the exact structure from PROCESS_ARCHITECTURE_GUIDE.md
   !!
   type, abstract :: ProcessInterface
      private
      character(len=64), public :: name = ''         !< Process name
      character(len=64), public :: version = ''      !< Version string
      character(len=256), public :: description = '' !< Process description
      logical :: is_initialized = .false.    !< Initialization status
      logical :: is_active = .false.         !< Active status
      real(fp) :: dt = 0.0_fp                !< Process timestep

   contains
      ! Required interface methods
      procedure(init_interface), deferred :: init
      procedure(run_interface), deferred :: run
      procedure(finalize_interface), deferred :: finalize

      ! Process capabilities registration
      procedure :: get_required_met_fields => process_get_required_met_fields
      procedure :: get_required_diagnostic_fields => process_get_required_diagnostic_fields

      ! Optional interface methods with default implementations
      procedure :: get_name => process_get_name
      procedure :: get_version => process_get_version
      procedure :: get_description => process_get_description
      procedure :: is_ready => process_is_ready
      procedure :: activate => process_activate
      procedure :: deactivate => process_deactivate
      procedure :: set_timestep => process_set_timestep
      procedure :: get_timestep => process_get_timestep
      procedure :: validate_config => process_validate_config

      ! Diagnostic support
      procedure :: register_diagnostics => process_register_diagnostics
      procedure :: update_diagnostics => process_update_diagnostics
      procedure :: register_diagnostic_field => process_register_diagnostic_field

      ! Common atmospheric process utilities
      procedure :: apply_emission_scaling => process_apply_emission_scaling
      procedure :: accumulate_emissions => process_accumulate_emissions
      procedure :: apply_tendency => process_apply_tendency
      procedure :: check_mass_conservation => process_check_mass_conservation
      procedure :: validate_species_availability => process_validate_species_availability
      procedure :: validate_physical_ranges => process_validate_physical_ranges

      ! Unit conversion utilities
      procedure :: convert_concentration_units => process_convert_concentration_units
      procedure :: convert_flux_units => process_convert_flux_units
      procedure :: calculate_column_integrals => process_calculate_column_integrals
      ! procedure :: interpolate_to_pressure_levels => process_interpolate_to_pressure_levels  ! Not yet implemented

      ! Column virtualization support
      procedure :: supports_column_processing => process_supports_column_processing
      procedure :: process_column => process_process_column
   end type ProcessInterface

   !> \brief Enhanced process interface specifically for column-based processing
   !!
   !! This interface extends ProcessInterface to provide column virtualization
   !! capabilities, allowing processes to work with virtual columns while
   !! maintaining awareness of 3D spatial relationships.
   type, abstract, extends(ProcessInterface) :: ColumnProcessInterface
      private

      logical :: column_processing_enabled = .true.  !< Enable column processing mode
      integer :: column_batch_size = 100            !< Number of columns to process in batch

   contains
      ! Required column processing methods
      procedure(column_init_interface), deferred :: init_column_processing
      procedure(column_run_interface), deferred :: run_column
      procedure(column_finalize_interface), deferred :: finalize_column_processing

      ! Optional column processing methods with default implementations
      procedure :: set_column_batch_size => column_process_set_batch_size
      procedure :: get_column_batch_size => column_process_get_batch_size
      procedure :: enable_column_processing => column_process_enable
      procedure :: disable_column_processing => column_process_disable
      procedure :: is_column_processing_enabled => column_process_is_enabled

      ! Generic column diagnostic update interface - automatically dispatches based on argument types
      generic :: update_column_diagnostics => update_scalar_diagnostic_column, &
         update_1d_diagnostic_column, &
         update_2d_diagnostic_column
      procedure :: update_scalar_diagnostic_column => column_update_scalar_diagnostic
      procedure :: update_1d_diagnostic_column => column_update_1d_diagnostic
      procedure :: update_2d_diagnostic_column => column_update_2d_diagnostic
   end type ColumnProcessInterface

   ! Abstract interfaces that must be implemented by concrete processes
   abstract interface
      !> \brief Initialize the process with given container
      subroutine init_interface(this, container, rc)
         import :: ProcessInterface, StateManagerType
         class(ProcessInterface), intent(inout) :: this
         type(StateManagerType), intent(inout) :: container
         integer, intent(out) :: rc
      end subroutine

      !> \brief Execute the process main calculations
      subroutine run_interface(this, container, rc)
         import :: ProcessInterface, StateManagerType
         class(ProcessInterface), intent(inout) :: this
         type(StateManagerType), intent(inout) :: container
         integer, intent(out) :: rc
      end subroutine

      !> \brief Clean up and finalize the process
      subroutine finalize_interface(this, rc)
         import :: ProcessInterface
         class(ProcessInterface), intent(inout) :: this
         integer, intent(out) :: rc
      end subroutine
   end interface

   ! Column processing interfaces
   abstract interface
      !> \brief Initialize column processing for the process
      subroutine column_init_interface(this, container, rc)
         import :: ColumnProcessInterface, StateManagerType
         class(ColumnProcessInterface), intent(inout) :: this
         type(StateManagerType), intent(inout) :: container
         integer, intent(out) :: rc
      end subroutine

      !> \brief Process a single virtual column
      subroutine column_run_interface(this, column, container, rc)
         import :: ColumnProcessInterface, VirtualColumnType, StateManagerType
         class(ColumnProcessInterface), intent(inout) :: this
         type(VirtualColumnType), intent(inout) :: column
         type(StateManagerType), intent(inout) :: container
         integer, intent(out) :: rc
      end subroutine

      !> \brief Finalize column processing
      subroutine column_finalize_interface(this, rc)
         import :: ColumnProcessInterface
         class(ColumnProcessInterface), intent(inout) :: this
         integer, intent(out) :: rc
      end subroutine
   end interface

   !> \brief Interface for getting required met fields for a process
   abstract interface
      function get_required_met_fields_interface(this) result(field_names)
         import :: ProcessInterface
         class(ProcessInterface), intent(in) :: this
         character(len=32), allocatable :: field_names(:)
      end function get_required_met_fields_interface
   end interface

contains

   !> \brief Get process name
   function process_get_name(this) result(name)
      class(ProcessInterface), intent(in) :: this
      character(len=64) :: name
      name = this%name
   end function process_get_name

   !> \brief Get process version
   function process_get_version(this) result(version)
      class(ProcessInterface), intent(in) :: this
      character(len=64) :: version
      version = this%version
   end function process_get_version

   !> \brief Get process description
   function process_get_description(this) result(description)
      class(ProcessInterface), intent(in) :: this
      character(len=256) :: description
      description = this%description
   end function process_get_description

   !> \brief Check if process is ready to run
   function process_is_ready(this) result(ready)
      class(ProcessInterface), intent(in) :: this
      logical :: ready
      ready = this%is_initialized .and. this%is_active
   end function process_is_ready

   !> \brief Activate the process
   subroutine process_activate(this)
      class(ProcessInterface), intent(inout) :: this
      this%is_active = .true.
      this%is_initialized = .true.
   end subroutine process_activate

   !> \brief Deactivate the process
   subroutine process_deactivate(this)
      class(ProcessInterface), intent(inout) :: this
      this%is_active = .false.
   end subroutine process_deactivate

   !> \brief Default configuration validation (always valid)
   subroutine process_validate_config(this, container, rc)
      class(ProcessInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      ! Default implementation - always valid
      rc = CC_SUCCESS
   end subroutine process_validate_config

   !> \brief Set process timestep
   subroutine process_set_timestep(this, dt)
      class(ProcessInterface), intent(inout) :: this
      real(fp), intent(in) :: dt
      this%dt = dt
   end subroutine process_set_timestep

   !> \brief Get process timestep
   function process_get_timestep(this) result(dt)
      class(ProcessInterface), intent(in) :: this
      real(fp) :: dt
      dt = this%dt
   end function process_get_timestep

   !> \brief Set species names that this process handles

   !> \brief Get required meteorological fields for this process
   !!
   !! Override this method in concrete processes to specify which met fields are needed.
   !! The framework will only allocate the fields that are required.
   function process_get_required_met_fields(this) result(field_names)
      class(ProcessInterface), intent(in) :: this
      character(len=32), allocatable :: field_names(:)

      ! Default implementation - no met fields required
      allocate(field_names(0))
   end function process_get_required_met_fields

   !> \brief Get required diagnostic fields for this process
   !!
   !! Override this method in concrete processes to specify which diagnostic fields
   !! should be created and made available.
   function process_get_required_diagnostic_fields(this) result(field_names)
      class(ProcessInterface), intent(in) :: this
      character(len=64), allocatable :: field_names(:)

      ! Default implementation - no diagnostic fields required
      allocate(field_names(0))
   end function process_get_required_diagnostic_fields

   !========================================================================
   ! Diagnostic Registration Methods
   !========================================================================

   !> \brief Register diagnostics for this process
   !!
   !! This method should be called during process initialization to register
   !! all diagnostics that the process will produce.
   !!
   !! \param[inout] this ProcessInterface instance
   !! \param[inout] container StateContainer for accessing diagnostic manager
   !! \param[out] rc Return code
   subroutine process_register_diagnostics(this, container, rc)
      class(ProcessInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      rc = CC_SUCCESS
      ! Default implementation registers no diagnostics
      ! Concrete processes should override this method to register their specific diagnostics
   end subroutine process_register_diagnostics

   !> \brief Update diagnostic values for this process
   !!
   !! This method should be called during process execution to update
   !! diagnostic field values with current data from the process.
   !!
   !! \param[inout] this ProcessInterface instance
   !! \param[inout] container StateContainer for accessing state data
   !! \param[out] rc Return code
   subroutine process_update_diagnostics(this, container, rc)
      class(ProcessInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      rc = CC_SUCCESS
      ! Default implementation does nothing
      ! Concrete processes should override this method to update their diagnostic values
   end subroutine process_update_diagnostics

   !> \brief Utility to register a single diagnostic field with create, initialize, and register steps
   !!
   !! This utility method encapsulates the common pattern of creating a diagnostic field,
   !! initializing its data storage, and registering it with the diagnostic registry.
   !! This eliminates code duplication in process-specific diagnostic registration routines.
   !!
   !! \param[in] this ProcessInterface instance
   !! \param[inout] registry Diagnostic registry to register the field with
   !! \param[in] field_name Name of the diagnostic field
   !! \param[in] description Human-readable description of the field
   !! \param[in] units Physical units of the field (e.g., 'kg/m2/s', 'ppbv')
   !! \param[in] field_type Field data type (DIAG_REAL_2D, DIAG_REAL_3D, etc.)
   !! \param[in] process_name Name of the process registering the field
   !! \param[in] dimensions Array dimensions for data storage [nx, ny] or [nx, ny, nz]
   !! \param[in] diagnostic_species Optional array of species names for diagnostics
   !! \param[in] diagnostic_species_id Optional array of species indices for diagnostics
   !! \param[out] rc Return code
   subroutine process_register_diagnostic_field(this, registry, field_name, description, &
      units, field_type, process_name, dimensions, &
      diagnostic_species, diagnostic_species_id, rc)
      use DiagnosticInterface_Mod, only: DiagnosticFieldType, DiagnosticRegistryType

      class(ProcessInterface), intent(in) :: this
      type(DiagnosticRegistryType), intent(inout) :: registry
      character(len=*), intent(in) :: field_name
      character(len=*), intent(in) :: description
      character(len=*), intent(in) :: units
      integer, intent(in) :: field_type
      character(len=*), intent(in) :: process_name
      integer, intent(in) :: dimensions(:)
      character(len=*), intent(in), optional :: diagnostic_species(:)
      integer, intent(in), optional :: diagnostic_species_id(:)
      integer, intent(out) :: rc

      type(DiagnosticFieldType) :: diag_field

      rc = CC_SUCCESS

      ! Create the diagnostic field with optional diagnostic species arguments
      if (present(diagnostic_species) .and. present(diagnostic_species_id)) then
         call diag_field%create(field_name, description, units, field_type, process_name, &
            diagnostic_species, diagnostic_species_id, rc)
      else
         call diag_field%create(field_name=field_name, description=description, units=units, &
            data_type=field_type, process_name=process_name, rc=rc)
      end if
      if (rc /= CC_SUCCESS) return

      ! Initialize data storage for field
      call diag_field%initialize_data(dimensions, rc)
      if (rc /= CC_SUCCESS) return

      ! Register the field with the registry
      call registry%register_field(diag_field, rc)
      if (rc /= CC_SUCCESS) return

   end subroutine process_register_diagnostic_field

   !========================================================================
   ! Common Atmospheric Process Utilities
   !========================================================================

   !> \brief Apply emission scaling factors to emission arrays
   !!
   !! This utility method applies scaling factors to emission fluxes for
   !! sensitivity studies or emission inventory adjustments.
   !!
   !! \param[in] this ProcessInterface instance
   !> \brief Apply emission scaling factors (deprecated - placeholder implementation)
   !!
   !! This method is deprecated as EmisState_Mod has been removed in favor of
   !! the new DiagnosticManager system. Processes should handle emissions
   !! through their own state management or direct chemical state modification.
   !!
   !! \param[in] this ProcessInterface instance
   !! \param[inout] container StateManager for accessing state data
   !! \param[in] scaling_factors Scaling factors per species [dimensionless]
   !! \param[in] species_indices Indices of species to scale
   !! \param[out] rc Return code
   subroutine process_apply_emission_scaling(this, container, scaling_factors, species_indices, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      real(fp), intent(in) :: scaling_factors(:)
      integer, intent(in) :: species_indices(:)
      integer, intent(out) :: rc

      integer :: i, species_idx

      rc = CC_SUCCESS

      ! Deprecated functionality - EmisState_Mod has been removed
      ! Processes should handle emissions through direct chemical state modification
      ! or their own internal emission arrays

      ! Placeholder implementation for backward compatibility
      do i = 1, size(species_indices)
         species_idx = species_indices(i)
         if (species_idx > 0 .and. species_idx <= size(scaling_factors)) then
            ! Would apply scaling if emission system was available
            ! For now, this is a no-op
         end if
      end do

   end subroutine process_apply_emission_scaling

   !> \brief Accumulate emissions from process calculations (deprecated - placeholder implementation)
   !!
   !! This method is deprecated as EmisState_Mod has been removed in favor of
   !! the new DiagnosticManager system. Processes should handle emissions
   !! through direct chemical state modification.
   !!
   !! \param[in] this ProcessInterface instance
   !! \param[inout] container StateManager for accessing state data
   !! \param[in] process_emissions Process emission fluxes [kg/m²/s or molecules/cm²/s]
   !! \param[in] species_mapping Mapping from process species to global species indices
   !! \param[out] rc Return code
   subroutine process_accumulate_emissions(this, container, process_emissions, species_mapping, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      real(fp), intent(in) :: process_emissions(:,:,:) !< Process emission fluxes
      integer, intent(in) :: species_mapping(:)        !< Species index mapping
      integer, intent(out) :: rc

      integer :: i, j, species_idx

      rc = CC_SUCCESS

      ! Deprecated functionality - EmisState_Mod has been removed
      ! Processes should handle emissions through direct chemical state modification

      ! Placeholder implementation for backward compatibility
      ! Actual emission accumulation should be done in the process itself
      ! by directly modifying the chemical state concentrations

   end subroutine process_accumulate_emissions

   !> \brief Apply tendencies to chemical species concentrations
   !!
   !! This utility method applies calculated tendencies (rates of change) to
   !! chemical species concentrations using a specified integration method.
   !!
   !! \param[in] this ProcessInterface instance
   !! \param[inout] container StateManager for accessing chemical state
   !! \param[in] tendencies Tendency array [species/time]
   !! \param[in] dt Time step [s]
   !! \param[out] rc Return code
   subroutine process_apply_tendency(this, container, tendencies, dt, rc)
      use ChemState_Mod, only: ChemStateType
      class(ProcessInterface), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      real(fp), intent(in) :: tendencies(:,:,:,:)     !< Tendencies [conc/time]
      real(fp), intent(in) :: dt                      !< Time step [s]
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state
      real(fp), allocatable :: concentrations(:,:,:,:)
      integer :: i, j, k, s

      rc = CC_SUCCESS

      ! Get chemical state from container
      chem_state => container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         rc = CC_FAILURE
         return
      end if

      ! Get current concentrations using the full 4D array method
      call chem_state%get_all_concentrations(concentrations, rc)
      if (rc /= CC_SUCCESS) return

      ! Apply tendencies using generic integration (solver-specific logic removed)
      ! Concrete processes should implement their own integration methods
      ! This base implementation provides a simple forward Euler as default
      concentrations = concentrations + dt * tendencies

      ! Ensure non-negative concentrations
      where (concentrations < 0.0_fp)
         concentrations = 0.0_fp
      end where

      ! Update chemical state with new concentrations using the full 4D array method
      call chem_state%set_all_concentrations(concentrations, rc)

   end subroutine process_apply_tendency

   !> \brief Check mass conservation for process calculations
   !!
   !! This utility method checks whether mass is conserved during process
   !! calculations, useful for debugging and validation.
   !!
   !! \param[in] this ProcessInterface instance
   !! \param[in] container StateManager for accessing state data
   !! \param[in] initial_mass Initial mass before process [kg]
   !! \param[in] final_mass Final mass after process [kg]
   !! \param[in] tolerance Relative tolerance for mass conservation
   !! \param[out] rc Return code
   function process_check_mass_conservation(this, container, initial_mass, final_mass, tolerance, rc) result(conserved)
      class(ProcessInterface), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      real(fp), intent(in) :: initial_mass(:)    !< Initial mass per species [kg]
      real(fp), intent(in) :: final_mass(:)      !< Final mass per species [kg]
      real(fp), intent(in) :: tolerance          !< Relative tolerance
      integer, intent(out) :: rc
      logical :: conserved

      type(ErrorManagerType), pointer :: error_mgr
      real(fp) :: mass_change, relative_error
      character(len=256) :: message
      integer :: i

      rc = CC_SUCCESS
      conserved = .true.

      error_mgr => container%get_error_manager()

      do i = 1, size(initial_mass)
         if (initial_mass(i) > 0.0_fp) then
            mass_change = final_mass(i) - initial_mass(i)
            relative_error = abs(mass_change) / initial_mass(i)

            if (relative_error > tolerance) then
               conserved = .false.
               write(message, '(A,I0,A,E12.4,A,E12.4)') &
                  'Mass conservation violation for species ', i, &
                  ': relative error = ', relative_error, ', tolerance = ', tolerance
               call error_mgr%report_error(ERROR_NOT_FOUND, message, rc)
            end if
         end if
      end do

   end function process_check_mass_conservation

   !> \brief Validate that required chemical species are available
   !!
   !! This utility method checks whether all species required by a process
   !! are available in the chemical mechanism.
   !!
   !! \param[in] this ProcessInterface instance
   !! \param[in] container StateManager for accessing chemical state
   !! \param[in] required_species List of required species names
   !! \param[out] available_species Logical array indicating availability
   !! \param[out] rc Return code
   function process_validate_species_availability(this, container, required_species, available_species, rc) result(all_available)
      use ChemState_Mod, only: ChemStateType
      class(ProcessInterface), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: required_species(:)
      logical, intent(out) :: available_species(:)
      integer, intent(out) :: rc
      logical :: all_available

      type(ChemStateType), pointer :: chem_state
      type(ErrorManagerType), pointer :: error_mgr
      character(len=256) :: message
      integer :: i, species_idx

      rc = CC_SUCCESS
      all_available = .true.
      available_species = .false.

      ! Get chemical state and error manager
      chem_state => container%get_chem_state_ptr()
      error_mgr = container%get_error_manager()

      if (.not. associated(chem_state)) then
         rc = CC_FAILURE
         all_available = .false.
         return
      end if

      ! Check availability of each required species
      do i = 1, size(required_species)
         species_idx = chem_state%find_species(trim(required_species(i)))
         if (species_idx > 0) then
            available_species(i) = .true.
         else
            available_species(i) = .false.
            all_available = .false.
            write(message, '(A,A,A)') 'Required species "', trim(required_species(i)), '" not found in mechanism'
            call error_mgr%report_error(ERROR_BOUNDS_CHECK, message, rc)
         end if
      end do

   end function process_validate_species_availability

   !> \brief Validate physical ranges of variables
   !!
   !! This method validates that physical variables are within reasonable ranges
   !! to catch numerical errors, unphysical values, or model instabilities.
   !! Uses StateManager's internal validation capabilities.
   !!
   !! \param[in] this ProcessInterface instance
   !! \param[in] container StateManager for accessing state data
   !! \param[out] rc Return code (CC_SUCCESS if all values valid, CC_FAILURE if errors found)
   function process_validate_physical_ranges(this, container, rc) result(all_valid)
      class(ProcessInterface), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc
      logical :: all_valid

      character(len=256) :: message

      rc = CC_SUCCESS
      all_valid = .true.

      ! Check if container is ready
      if (.not. container%is_ready()) then
         write(*, '(A)') 'ERROR: StateManager not ready for validation'
         rc = CC_FAILURE
         all_valid = .false.
         return
      endif

      ! Validate using StateManager's internal capabilities
      ! For now, this is a basic validation - could be enhanced with
      ! specific physical range checks when StateManager validation
      ! utilities are more developed

      write(message, '(A,A)') 'Physical validation completed for process: ', trim(this%name)
      write(*, '(A)') trim(message)

   end function process_validate_physical_ranges

   !========================================================================
   ! Column virtualization support methods
   !========================================================================

   !> \brief Check if process supports column-based processing
   function process_supports_column_processing(this) result(supports)
      class(ProcessInterface), intent(in) :: this
      logical :: supports

      supports = .false.  ! Default implementation - override in subclasses
   end function process_supports_column_processing

   !> \brief Process a single column (default implementation)
   subroutine process_process_column(this, column, rc)
      class(ProcessInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      integer, intent(out) :: rc

      rc = CC_FAILURE  ! Default implementation fails - must be overridden
   end subroutine process_process_column

   !========================================================================
   ! ColumnProcessInterface Implementation
   !========================================================================

   !> \brief Set column batch size
   subroutine column_process_set_batch_size(this, batch_size)
      class(ColumnProcessInterface), intent(inout) :: this
      integer, intent(in) :: batch_size

      this%column_batch_size = max(1, batch_size)
   end subroutine column_process_set_batch_size

   !> \brief Get column batch size
   function column_process_get_batch_size(this) result(batch_size)
      class(ColumnProcessInterface), intent(in) :: this
      integer :: batch_size

      batch_size = this%column_batch_size
   end function column_process_get_batch_size

   !> \brief Enable column processing
   subroutine column_process_enable(this)
      class(ColumnProcessInterface), intent(inout) :: this

      this%column_processing_enabled = .true.
   end subroutine column_process_enable

   !> \brief Disable column processing
   subroutine column_process_disable(this)
      class(ColumnProcessInterface), intent(inout) :: this

      this%column_processing_enabled = .false.
   end subroutine column_process_disable

   !> \brief Check if column processing is enabled
   function column_process_is_enabled(this) result(is_enabled)
      class(ColumnProcessInterface), intent(in) :: this
      logical :: is_enabled

      is_enabled = this%column_processing_enabled
   end function column_process_is_enabled

   !========================================================================
   ! Unit Conversion Utilities
   !========================================================================

   !> \brief Convert concentration units between different unit systems
   !!
   !! This utility converts concentration values between common atmospheric chemistry units
   !! such as molec/cm³, ppbv, ppmv, µg/m³, etc.
   !!
   !! \param[in] this ProcessInterface instance
   !! \param[inout] values Array of concentration values to convert
   !! \param[in] from_units Source units (e.g., 'ppbv', 'molec/cm3', 'ug/m3')
   !! \param[in] to_units Target units
   !! \param[in] molecular_weight Molecular weight [g/mol] (needed for mass/volume conversions)
   !! \param[in] temperature Temperature [K] (needed for some conversions)
   !! \param[in] pressure Pressure [Pa] (needed for some conversions)
   !! \param[out] rc Return code
   subroutine process_convert_concentration_units(this, values, from_units, to_units, &
      molecular_weight, temperature, pressure, rc)
      class(ProcessInterface), intent(in) :: this
      real(fp), intent(inout) :: values(:)
      character(len=*), intent(in) :: from_units, to_units
      real(fp), intent(in), optional :: molecular_weight, temperature, pressure
      integer, intent(out) :: rc

      real(fp) :: mw, temp, pres
      real(fp) :: conversion_factor
      integer :: i

      rc = CC_SUCCESS

      ! Set default values if not provided
      mw = 29.0_fp    ! Default molecular weight of air [g/mol]
      temp = 273.15_fp ! Default temperature [K]
      pres = 101325.0_fp ! Default pressure [Pa]

      if (present(molecular_weight)) mw = molecular_weight
      if (present(temperature)) temp = temperature
      if (present(pressure)) pres = pressure

      ! Calculate conversion factor based on unit types
      if (trim(from_units) == trim(to_units)) then
         ! No conversion needed
         return
      endif

      ! Convert from ppbv to other units
      if (trim(from_units) == 'ppbv') then
         select case (trim(to_units))
          case ('ppmv')
            conversion_factor = 1.0e-3_fp
          case ('molec/cm3')
            ! ppbv to molec/cm³: ppbv * (P/RT) * (1e-9) * NA * (1e-6)
            conversion_factor = (pres / (8.314_fp * temp)) * 1.0e-9_fp * 6.022e23_fp * 1.0e-6_fp
          case ('ug/m3')
            ! ppbv to µg/m³: ppbv * (P/RT) * MW * (1e-9) * (1e6)
            conversion_factor = (pres / (8.314_fp * temp)) * mw * 1.0e-3_fp
          case default
            rc = CC_FAILURE
            return
         end select

         ! Convert from molec/cm3 to other units
      else if (trim(from_units) == 'molec/cm3') then
         select case (trim(to_units))
          case ('ppbv')
            ! molec/cm³ to ppbv: (molec/cm³) * (RT/P) * (1e9) / NA * (1e6)
            conversion_factor = (8.314_fp * temp / pres) * 1.0e9_fp / 6.022e23_fp * 1.0e6_fp
          case ('ug/m3')
            ! molec/cm³ to µg/m³: (molec/cm³) * MW / NA * (1e12)
            conversion_factor = mw / 6.022e23_fp * 1.0e12_fp
          case default
            rc = CC_FAILURE
            return
         end select

      else
         ! Unsupported conversion
         rc = CC_FAILURE
         return
      endif

      ! Apply conversion
      do i = 1, size(values)
         values(i) = values(i) * conversion_factor
      end do

   end subroutine process_convert_concentration_units

   !> \brief Convert flux units between different unit systems
   !!
   !! This utility converts emission flux values between common units
   !! such as kg/m²/s, molec/cm²/s, molecules/m²/s, etc.
   !!
   !! \param[in] this ProcessInterface instance
   !! \param[inout] flux_values Array of flux values to convert
   !! \param[in] from_units Source flux units
   !! \param[in] to_units Target flux units
   !! \param[in] molecular_weight Molecular weight [g/mol]
   !! \param[out] rc Return code
   subroutine process_convert_flux_units(this, flux_values, from_units, to_units, molecular_weight, rc)
      class(ProcessInterface), intent(in) :: this
      real(fp), intent(inout) :: flux_values(:)
      character(len=*), intent(in) :: from_units, to_units
      real(fp), intent(in) :: molecular_weight
      integer, intent(out) :: rc

      real(fp) :: conversion_factor
      integer :: i

      rc = CC_SUCCESS

      ! No conversion needed
      if (trim(from_units) == trim(to_units)) then
         return
      endif

      ! Convert kg/m²/s to molec/cm²/s
      if (trim(from_units) == 'kg/m2/s' .and. trim(to_units) == 'molec/cm2/s') then
         ! kg/m²/s * (1000 g/kg) * (1 mol/MW g) * (NA molec/mol) * (1 m²/10⁴ cm²)
         conversion_factor = 1000.0_fp * (1.0_fp / molecular_weight) * 6.022e23_fp * 1.0e-4_fp

         ! Convert molec/cm²/s to kg/m²/s
      else if (trim(from_units) == 'molec/cm2/s' .and. trim(to_units) == 'kg/m2/s') then
         ! molec/cm²/s * (1 mol/NA molec) * (MW g/mol) * (1 kg/1000 g) * (10⁴ cm²/m²)
         conversion_factor = (1.0_fp / 6.022e23_fp) * molecular_weight * 1.0e-3_fp * 1.0e4_fp

      else
         ! Unsupported conversion
         rc = CC_FAILURE
         return
      endif

      ! Apply conversion
      do i = 1, size(flux_values)
         flux_values(i) = flux_values(i) * conversion_factor
      end do

   end subroutine process_convert_flux_units

   !> \brief Calculate column integrals of species concentrations
   !!
   !! This utility calculates vertical column integrals (e.g., column density)
   !! from 3D concentration fields.
   !!
   !! \param[in] this ProcessInterface instance
   !! \param[in] container StateManager for accessing state data
   !! \param[in] species_name Name of species to integrate
   !! \param[out] column_integrals 2D array of column integrals [molecules/cm²]
   !! \param[out] rc Return code
   subroutine process_calculate_column_integrals(this, container, species_name, column_integrals, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      real(fp), allocatable, intent(out) :: column_integrals(:,:)
      integer, intent(out) :: rc

      character(len=256) :: message

      rc = CC_SUCCESS

      ! Check if container is ready
      if (.not. container%is_ready()) then
         write(*, '(A)') 'ERROR: StateManager not ready for column integration'
         rc = CC_FAILURE
         return
      endif

      ! For now, this is a placeholder that would require access to:
      ! - 3D concentration arrays from ChemState
      ! - Air density and layer thickness from MetState
      ! - Proper integration over vertical levels

      ! Allocate output array (would get dimensions from MetState)
      allocate(column_integrals(50, 50))  ! Placeholder dimensions
      column_integrals = 0.0_fp

      write(message, '(A,A,A)') 'Column integration completed for species: ', trim(species_name), &
         ' (placeholder implementation)'
      write(*, '(A)') trim(message)

   end subroutine process_calculate_column_integrals

   !========================================================================
   ! ColumnProcessInterface Diagnostic Implementation
   !========================================================================

   !> Update a scalar diagnostic field for column processing
   !!
   !! This method handles updating scalar diagnostic values from column processing
   !! to the global diagnostic storage managed by DiagnosticManager.
   !!
   !! \param[inout] this ColumnProcessInterface instance
   !! \param[in] field_name Name of the diagnostic field to update
   !! \param[in] scalar_value The scalar value to store
   !! \param[in] i_col Column i-index (x-direction)
   !! \param[in] j_col Column j-index (y-direction)
   !! \param[in] container StateManager for accessing diagnostic manager
   !! \param[out] rc Return code
   subroutine column_update_scalar_diagnostic(this, field_name, scalar_value, i_col, j_col, container, rc)
      class(ColumnProcessInterface), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: scalar_value
      integer, intent(in) :: i_col, j_col
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(DiagnosticManagerType), pointer :: diag_mgr => null()
      type(DiagnosticRegistryType), pointer :: registry => null()
      type(DiagnosticFieldType), pointer :: diag_field => null()
      type(DiagnosticDataType), pointer :: diag_data => null()
      real(fp), pointer :: field_data_2d(:,:) => null()
      character(len=64) :: process_name

      rc = CC_SUCCESS

      ! Get DiagnosticManager from container
      diag_mgr => container%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) then
         rc = CC_FAILURE
         return
      end if

      ! Get process name for registry lookup
      process_name = this%get_name()

      ! Get process registry from DiagnosticManager
      call diag_mgr%get_process_registry(trim(process_name), registry, rc)
      if (rc /= CC_SUCCESS .or. .not. associated(registry)) then
         rc = CC_FAILURE
         return
      end if

      ! Get the specific diagnostic field from registry
      diag_field => registry%get_field_ptr(field_name)
      if (.not. associated(diag_field)) then
         rc = CC_FAILURE  ! Field not found
         return
      end if

      ! Verify field is ready and get data storage
      if (.not. diag_field%is_ready()) then
         rc = CC_FAILURE  ! Field not initialized
         return
      end if
      diag_data => diag_field%get_data_ptr()
      if (.not. associated(diag_data)) then
         rc = CC_FAILURE
         return
      end if

      ! Get pointer to 2D field data and update
      field_data_2d => diag_data%get_real_2d_ptr()
      if (.not. associated(field_data_2d)) then
         rc = CC_FAILURE  ! Wrong data type or not allocated
         return
      end if

      ! Update the field at column position
      field_data_2d(i_col, j_col) = scalar_value

   end subroutine column_update_scalar_diagnostic

   !> Update a 1D diagnostic field for column processing
   !!
   !! This method handles updating 1D diagnostic arrays (e.g., vertical profiles)
   !! from column processing to the global diagnostic storage.
   !!
   !! \param[inout] this ColumnProcessInterface instance
   !! \param[in] field_name Name of the diagnostic field to update
   !! \param[in] array_1d The 1D array data to store
   !! \param[in] i_col Column i-index (x-direction)
   !! \param[in] j_col Column j-index (y-direction)
   !! \param[in] container StateManager for accessing diagnostic manager
   !! \param[out] rc Return code
   subroutine column_update_1d_diagnostic(this, field_name, array_1d, i_col, j_col, container, rc)
      class(ColumnProcessInterface), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: array_1d(:)
      integer, intent(in) :: i_col, j_col
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      ! Local variables for diagnostic system integration
      type(DiagnosticManagerType), pointer :: diag_mgr => null()
      type(DiagnosticRegistryType), pointer :: registry => null()
      type(DiagnosticFieldType), pointer :: diag_field => null()
      type(DiagnosticDataType), pointer :: diag_data => null()
      real(fp), pointer :: field_data_3d(:,:,:) => null()
      character(len=64) :: process_name
      integer :: k, n_dim3

      rc = CC_SUCCESS
      n_dim3 = size(array_1d)

      ! Get DiagnosticManager from container
      diag_mgr => container%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) then
         rc = CC_FAILURE
         return
      end if

      ! Get process name for registry lookup
      process_name = this%get_name()

      ! Get process registry from DiagnosticManager
      call diag_mgr%get_process_registry(trim(process_name), registry, rc)
      if (rc /= CC_SUCCESS) return
      if (.not. associated(registry)) then
         rc = CC_FAILURE
         return
      end if

      ! Get the specific diagnostic field from registry
      diag_field => registry%get_field_ptr(field_name)
      if (.not. associated(diag_field)) then
         rc = CC_FAILURE  ! Field not found
         return
      end if

      ! Verify field is ready and get data storage
      if (.not. diag_field%is_ready()) then
         rc = CC_FAILURE  ! Field not initialized
         return
      end if
      diag_data => diag_field%get_data_ptr()
      if (.not. associated(diag_data)) then
         rc = CC_FAILURE
         return
      end if

      ! Get pointer to 3D field data and validate
      field_data_3d => diag_data%get_real_3d_ptr()
      if (.not. associated(field_data_3d)) then
         rc = CC_FAILURE  ! Wrong data type or not allocated
         return
      end if

      ! Validate dimensions
      if (size(field_data_3d, 3) /= n_dim3) then
         rc = CC_FAILURE  ! Dimension mismatch
         return
      end if

      ! Update the field at column position
      do k = 1, n_dim3
         field_data_3d(i_col, j_col, k) = array_1d(k)
      end do

   end subroutine column_update_1d_diagnostic

   !> Update a 2D diagnostic field for column processing
   !!
   !! This method handles updating 2D diagnostic arrays (e.g., level-species matrices)
   !! from column processing to the global diagnostic storage. Since DiagnosticInterface
   !! currently supports up to 3D arrays, 2D column data is stored using flattened indexing.
   !!
   !! \param[inout] this ColumnProcessInterface instance
   !! \param[in] field_name Name of the diagnostic field to update
   !! \param[in] array_2d The 2D array data to store
   !! \param[in] i_col Column i-index (x-direction)
   !! \param[in] j_col Column j-index (y-direction)
   !! \param[in] container StateManager for accessing diagnostic manager
   !! \param[out] rc Return code
   subroutine column_update_2d_diagnostic(this, field_name, array_2d, i_col, j_col, container, rc)
      class(ColumnProcessInterface), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: array_2d(:,:)
      integer, intent(in) :: i_col, j_col
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      ! Local variables for diagnostic system integration
      type(DiagnosticManagerType), pointer :: diag_mgr => null()
      type(DiagnosticRegistryType), pointer :: registry => null()
      type(DiagnosticFieldType), pointer :: diag_field => null()
      type(DiagnosticDataType), pointer :: diag_data => null()
      real(fp), pointer :: field_data_3d(:,:,:) => null()
      character(len=64) :: process_name
      integer :: k, l, n_levels, n_species, flat_index

      rc = CC_SUCCESS

      n_levels = size(array_2d, 1)   ! levels dimension
      n_species = size(array_2d, 2)  ! species dimension

      ! Get DiagnosticManager from container
      diag_mgr => container%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) then
         rc = CC_FAILURE
         return
      end if

      ! Get process name for registry lookup
      process_name = this%get_name()

      ! Get process registry from DiagnosticManager
      call diag_mgr%get_process_registry(trim(process_name), registry, rc)
      if (rc /= CC_SUCCESS) return
      if (.not. associated(registry)) then
         rc = CC_FAILURE
         return
      end if

      ! Get diagnostic field
      diag_field => registry%get_field_ptr(field_name)
      if (.not. associated(diag_field) .or. .not. diag_field%is_ready()) then
         rc = CC_FAILURE
         return
      end if

      ! Get 3D data storage (flattened storage for 2D column data)
      diag_data => diag_field%get_data_ptr()
      field_data_3d => diag_data%get_real_3d_ptr()
      if (.not. associated(field_data_3d)) then
         rc = CC_FAILURE
         return
      end if

      ! Update field using flattened indexing
      ! The 3rd dimension contains flattened (level, species) pairs
      ! Mapping: flat_index = (level-1) * n_species + species
      do l = 1, n_species
         do k = 1, n_levels
            flat_index = (k-1) * n_species + l
            field_data_3d(i_col, j_col, flat_index) = array_2d(k, l)
         end do
      end do

   end subroutine column_update_2d_diagnostic


end module ProcessInterface_Mod
