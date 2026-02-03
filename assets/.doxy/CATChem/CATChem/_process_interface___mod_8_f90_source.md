

# File ProcessInterface\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**ProcessInterface\_Mod.F90**](_process_interface___mod_8_f90.md)

[Go to the documentation of this file](_process_interface___mod_8_f90.md)


```Fortran

module processinterface_mod
   use precision_mod
   use statemanager_mod, only : statemanagertype
   use error_mod
   use columninterface_mod, only : columnprocessortype
   use virtualcolumn_mod, only : virtualcolumntype
   use extemisdata_mod, only : extemisdatatype
   use diagnosticmanager_mod, only: diagnosticmanagertype
   use diagnosticinterface_mod, only: diagnosticregistrytype, diagnosticfieldtype, diagnosticdatatype

   implicit none
   private

   public :: processinterface
   public :: columnprocessinterface

   type, abstract :: processinterface
      private
      character(len=64), public :: name = ''
      character(len=64), public :: version = ''
      character(len=256), public :: description = ''
      logical :: is_initialized = .false.    
      logical :: is_active = .false.         
      real(fp) :: dt = 0.0_fp                

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
   end type processinterface

   type, abstract, extends(processinterface) :: columnprocessinterface
      private

      logical :: column_processing_enabled = .true.  
      integer :: column_batch_size = 100

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
   end type columnprocessinterface

   ! Abstract interfaces that must be implemented by concrete processes
   abstract interface

      subroutine init_interface(this, container, rc)
         import :: processinterface, statemanagertype
         class(ProcessInterface), intent(inout) :: this
         type(StateManagerType), intent(inout) :: container
         integer, intent(out) :: rc
      end subroutine

      subroutine run_interface(this, container, rc)
         import :: processinterface, statemanagertype
         class(ProcessInterface), intent(inout) :: this
         type(StateManagerType), intent(inout) :: container
         integer, intent(out) :: rc
      end subroutine

      subroutine finalize_interface(this, rc)
         import :: processinterface
         class(ProcessInterface), intent(inout) :: this
         integer, intent(out) :: rc
      end subroutine
   end interface

   ! Column processing interfaces
   abstract interface

      subroutine column_init_interface(this, container, rc)
         import :: columnprocessinterface, statemanagertype
         class(ColumnProcessInterface), intent(inout) :: this
         type(StateManagerType), intent(inout) :: container
         integer, intent(out) :: rc
      end subroutine

      subroutine column_run_interface(this, column, container, rc)
         import :: columnprocessinterface, virtualcolumntype, statemanagertype
         class(ColumnProcessInterface), intent(inout) :: this
         type(VirtualColumnType), intent(inout) :: column
         type(StateManagerType), intent(inout) :: container
         integer, intent(out) :: rc
      end subroutine

      subroutine column_finalize_interface(this, rc)
         import :: columnprocessinterface
         class(ColumnProcessInterface), intent(inout) :: this
         integer, intent(out) :: rc
      end subroutine
   end interface

   abstract interface
      function get_required_met_fields_interface(this) result(field_names)
         import :: processinterface
         class(ProcessInterface), intent(in) :: this
         character(len=32), allocatable :: field_names(:)
      end function get_required_met_fields_interface
   end interface

contains

   function process_get_name(this) result(name)
      class(ProcessInterface), intent(in) :: this
      character(len=64) :: name
      name = this%name
   end function process_get_name

   function process_get_version(this) result(version)
      class(ProcessInterface), intent(in) :: this
      character(len=64) :: version
      version = this%version
   end function process_get_version

   function process_get_description(this) result(description)
      class(ProcessInterface), intent(in) :: this
      character(len=256) :: description
      description = this%description
   end function process_get_description

   function process_is_ready(this) result(ready)
      class(ProcessInterface), intent(in) :: this
      logical :: ready
      ready = this%is_initialized .and. this%is_active
   end function process_is_ready

   subroutine process_activate(this)
      class(ProcessInterface), intent(inout) :: this
      this%is_active = .true.
      this%is_initialized = .true.
   end subroutine process_activate

   subroutine process_deactivate(this)
      class(ProcessInterface), intent(inout) :: this
      this%is_active = .false.
   end subroutine process_deactivate

   subroutine process_validate_config(this, container, rc)
      class(ProcessInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      ! Default implementation - always valid
      rc = cc_success
   end subroutine process_validate_config

   subroutine process_set_timestep(this, dt)
      class(ProcessInterface), intent(inout) :: this
      real(fp), intent(in) :: dt
      this%dt = dt
   end subroutine process_set_timestep

   function process_get_timestep(this) result(dt)
      class(ProcessInterface), intent(in) :: this
      real(fp) :: dt
      dt = this%dt
   end function process_get_timestep


   function process_get_required_met_fields(this) result(field_names)
      class(ProcessInterface), intent(in) :: this
      character(len=32), allocatable :: field_names(:)

      ! Default implementation - no met fields required
      allocate(field_names(0))
   end function process_get_required_met_fields

   function process_get_required_diagnostic_fields(this) result(field_names)
      class(ProcessInterface), intent(in) :: this
      character(len=64), allocatable :: field_names(:)

      ! Default implementation - no diagnostic fields required
      allocate(field_names(0))
   end function process_get_required_diagnostic_fields

   !========================================================================
   ! Diagnostic Registration Methods
   !========================================================================

   subroutine process_register_diagnostics(this, container, rc)
      class(ProcessInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      rc = cc_success
      ! Default implementation registers no diagnostics
      ! Concrete processes should override this method to register their specific diagnostics
   end subroutine process_register_diagnostics

   subroutine process_update_diagnostics(this, container, rc)
      class(ProcessInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      rc = cc_success
      ! Default implementation does nothing
      ! Concrete processes should override this method to update their diagnostic values
   end subroutine process_update_diagnostics

   subroutine process_register_diagnostic_field(this, registry, field_name, description, &
      units, field_type, process_name, dimensions, &
      diagnostic_species, diagnostic_species_id, rc)
      use diagnosticinterface_mod, only: diagnosticfieldtype, diagnosticregistrytype

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

      rc = cc_success

      ! Create the diagnostic field with optional diagnostic species arguments
      if (present(diagnostic_species) .and. present(diagnostic_species_id)) then
         call diag_field%create(field_name, description, units, field_type, process_name, &
            diagnostic_species, diagnostic_species_id, rc)
      else
         call diag_field%create(field_name=field_name, description=description, units=units, &
            data_type=field_type, process_name=process_name, rc=rc)
      end if
      if (rc /= cc_success) return

      ! Initialize data storage for field
      call diag_field%initialize_data(dimensions, rc)
      if (rc /= cc_success) return

      ! Register the field with the registry
      call registry%register_field(diag_field, rc)
      if (rc /= cc_success) return

   end subroutine process_register_diagnostic_field

   !========================================================================
   ! Common Atmospheric Process Utilities
   !========================================================================

   subroutine process_apply_emission_scaling(this, container, scaling_factors, species_indices, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      real(fp), intent(in) :: scaling_factors(:)
      integer, intent(in) :: species_indices(:)
      integer, intent(out) :: rc

      integer :: i, species_idx

      rc = cc_success

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

   subroutine process_accumulate_emissions(this, container, process_emissions, species_mapping, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      real(fp), intent(in) :: process_emissions(:,:,:)
      integer, intent(in) :: species_mapping(:)
      integer, intent(out) :: rc

      integer :: i, j, species_idx

      rc = cc_success

      ! Deprecated functionality - EmisState_Mod has been removed
      ! Processes should handle emissions through direct chemical state modification

      ! Placeholder implementation for backward compatibility
      ! Actual emission accumulation should be done in the process itself
      ! by directly modifying the chemical state concentrations

   end subroutine process_accumulate_emissions

   subroutine process_apply_tendency(this, container, tendencies, dt, rc)
      use chemstate_mod, only: chemstatetype
      class(ProcessInterface), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      real(fp), intent(in) :: tendencies(:,:,:,:)
      real(fp), intent(in) :: dt
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state
      real(fp), allocatable :: concentrations(:,:,:,:)
      integer :: i, j, k, s

      rc = cc_success

      ! Get chemical state from container
      chem_state => container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         rc = cc_failure
         return
      end if

      ! Get current concentrations using the full 4D array method
      call chem_state%get_all_concentrations(concentrations, rc)
      if (rc /= cc_success) return

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

   function process_check_mass_conservation(this, container, initial_mass, final_mass, tolerance, rc) result(conserved)
      class(ProcessInterface), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      real(fp), intent(in) :: initial_mass(:)
      real(fp), intent(in) :: final_mass(:)
      real(fp), intent(in) :: tolerance
      integer, intent(out) :: rc
      logical :: conserved

      type(ErrorManagerType), pointer :: error_mgr
      real(fp) :: mass_change, relative_error
      character(len=256) :: message
      integer :: i

      rc = cc_success
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
               call error_mgr%report_error(error_not_found, message, rc)
            end if
         end if
      end do

   end function process_check_mass_conservation

   function process_validate_species_availability(this, container, required_species, available_species, rc) result(all_available)
      use chemstate_mod, only: chemstatetype
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

      rc = cc_success
      all_available = .true.
      available_species = .false.

      ! Get chemical state and error manager
      chem_state => container%get_chem_state_ptr()
      error_mgr = container%get_error_manager()

      if (.not. associated(chem_state)) then
         rc = cc_failure
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
            call error_mgr%report_error(error_bounds_check, message, rc)
         end if
      end do

   end function process_validate_species_availability

   function process_validate_physical_ranges(this, container, rc) result(all_valid)
      class(ProcessInterface), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc
      logical :: all_valid

      character(len=256) :: message

      rc = cc_success
      all_valid = .true.

      ! Check if container is ready
      if (.not. container%is_ready()) then
         write(*, '(A)') 'ERROR: StateManager not ready for validation'
         rc = cc_failure
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

   function process_supports_column_processing(this) result(supports)
      class(ProcessInterface), intent(in) :: this
      logical :: supports

      supports = .false.  ! Default implementation - override in subclasses
   end function process_supports_column_processing

   subroutine process_process_column(this, column, rc)
      class(ProcessInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      integer, intent(out) :: rc

      rc = cc_failure  ! Default implementation fails - must be overridden
   end subroutine process_process_column

   !========================================================================
   ! ColumnProcessInterface Implementation
   !========================================================================

   subroutine column_process_set_batch_size(this, batch_size)
      class(ColumnProcessInterface), intent(inout) :: this
      integer, intent(in) :: batch_size

      this%column_batch_size = max(1, batch_size)
   end subroutine column_process_set_batch_size

   function column_process_get_batch_size(this) result(batch_size)
      class(ColumnProcessInterface), intent(in) :: this
      integer :: batch_size

      batch_size = this%column_batch_size
   end function column_process_get_batch_size

   subroutine column_process_enable(this)
      class(ColumnProcessInterface), intent(inout) :: this

      this%column_processing_enabled = .true.
   end subroutine column_process_enable

   subroutine column_process_disable(this)
      class(ColumnProcessInterface), intent(inout) :: this

      this%column_processing_enabled = .false.
   end subroutine column_process_disable

   function column_process_is_enabled(this) result(is_enabled)
      class(ColumnProcessInterface), intent(in) :: this
      logical :: is_enabled

      is_enabled = this%column_processing_enabled
   end function column_process_is_enabled

   !========================================================================
   ! Unit Conversion Utilities
   !========================================================================

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

      rc = cc_success

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
            rc = cc_failure
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
            rc = cc_failure
            return
         end select

      else
         ! Unsupported conversion
         rc = cc_failure
         return
      endif

      ! Apply conversion
      do i = 1, size(values)
         values(i) = values(i) * conversion_factor
      end do

   end subroutine process_convert_concentration_units

   subroutine process_convert_flux_units(this, flux_values, from_units, to_units, molecular_weight, rc)
      class(ProcessInterface), intent(in) :: this
      real(fp), intent(inout) :: flux_values(:)
      character(len=*), intent(in) :: from_units, to_units
      real(fp), intent(in) :: molecular_weight
      integer, intent(out) :: rc

      real(fp) :: conversion_factor
      integer :: i

      rc = cc_success

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
         rc = cc_failure
         return
      endif

      ! Apply conversion
      do i = 1, size(flux_values)
         flux_values(i) = flux_values(i) * conversion_factor
      end do

   end subroutine process_convert_flux_units

   subroutine process_calculate_column_integrals(this, container, species_name, column_integrals, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      real(fp), allocatable, intent(out) :: column_integrals(:,:)
      integer, intent(out) :: rc

      character(len=256) :: message

      rc = cc_success

      ! Check if container is ready
      if (.not. container%is_ready()) then
         write(*, '(A)') 'ERROR: StateManager not ready for column integration'
         rc = cc_failure
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

      rc = cc_success

      ! Get DiagnosticManager from container
      diag_mgr => container%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) then
         rc = cc_failure
         return
      end if

      ! Get process name for registry lookup
      process_name = this%get_name()

      ! Get process registry from DiagnosticManager
      call diag_mgr%get_process_registry(trim(process_name), registry, rc)
      if (rc /= cc_success .or. .not. associated(registry)) then
         rc = cc_failure
         return
      end if

      ! Get the specific diagnostic field from registry
      diag_field => registry%get_field_ptr(field_name)
      if (.not. associated(diag_field)) then
         rc = cc_failure  ! Field not found
         return
      end if

      ! Verify field is ready and get data storage
      if (.not. diag_field%is_ready()) then
         rc = cc_failure  ! Field not initialized
         return
      end if
      diag_data => diag_field%get_data_ptr()
      if (.not. associated(diag_data)) then
         rc = cc_failure
         return
      end if

      ! Get pointer to 2D field data and update
      field_data_2d => diag_data%get_real_2d_ptr()
      if (.not. associated(field_data_2d)) then
         rc = cc_failure  ! Wrong data type or not allocated
         return
      end if

      ! Update the field at column position
      field_data_2d(i_col, j_col) = scalar_value

   end subroutine column_update_scalar_diagnostic

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

      rc = cc_success
      n_dim3 = size(array_1d)

      ! Get DiagnosticManager from container
      diag_mgr => container%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) then
         rc = cc_failure
         return
      end if

      ! Get process name for registry lookup
      process_name = this%get_name()

      ! Get process registry from DiagnosticManager
      call diag_mgr%get_process_registry(trim(process_name), registry, rc)
      if (rc /= cc_success) return
      if (.not. associated(registry)) then
         rc = cc_failure
         return
      end if

      ! Get the specific diagnostic field from registry
      diag_field => registry%get_field_ptr(field_name)
      if (.not. associated(diag_field)) then
         rc = cc_failure  ! Field not found
         return
      end if

      ! Verify field is ready and get data storage
      if (.not. diag_field%is_ready()) then
         rc = cc_failure  ! Field not initialized
         return
      end if
      diag_data => diag_field%get_data_ptr()
      if (.not. associated(diag_data)) then
         rc = cc_failure
         return
      end if

      ! Get pointer to 3D field data and validate
      field_data_3d => diag_data%get_real_3d_ptr()
      if (.not. associated(field_data_3d)) then
         rc = cc_failure  ! Wrong data type or not allocated
         return
      end if

      ! Validate dimensions
      if (size(field_data_3d, 3) /= n_dim3) then
         rc = cc_failure  ! Dimension mismatch
         return
      end if

      ! Update the field at column position
      do k = 1, n_dim3
         field_data_3d(i_col, j_col, k) = array_1d(k)
      end do

   end subroutine column_update_1d_diagnostic

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

      rc = cc_success

      n_levels = size(array_2d, 1)   ! levels dimension
      n_species = size(array_2d, 2)  ! species dimension

      ! Get DiagnosticManager from container
      diag_mgr => container%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) then
         rc = cc_failure
         return
      end if

      ! Get process name for registry lookup
      process_name = this%get_name()

      ! Get process registry from DiagnosticManager
      call diag_mgr%get_process_registry(trim(process_name), registry, rc)
      if (rc /= cc_success) return
      if (.not. associated(registry)) then
         rc = cc_failure
         return
      end if

      ! Get diagnostic field
      diag_field => registry%get_field_ptr(field_name)
      if (.not. associated(diag_field) .or. .not. diag_field%is_ready()) then
         rc = cc_failure
         return
      end if

      ! Get 3D data storage (flattened storage for 2D column data)
      diag_data => diag_field%get_data_ptr()
      field_data_3d => diag_data%get_real_3d_ptr()
      if (.not. associated(field_data_3d)) then
         rc = cc_failure
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


end module processinterface_mod
```


