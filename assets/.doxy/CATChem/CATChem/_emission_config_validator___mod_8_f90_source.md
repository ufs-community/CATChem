

# File EmissionConfigValidator\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**EmissionConfigValidator\_Mod.F90**](_emission_config_validator___mod_8_f90.md)

[Go to the documentation of this file](_emission_config_validator___mod_8_f90.md)


```Fortran

module emissionconfigvalidator_mod
   use precision_mod
   use error_mod
   use statemanager_mod, only : statemanagertype
   use chemstate_mod, only : chemstatetype
   ! use logging_mod, only : log_message, LOG_INFO

   implicit none
   private

   public :: emissionconfigvalidatortype
   public :: validationresulttype
   public :: validation_success, validation_warning, validation_error

   ! Validation result constants
   integer, parameter :: VALIDATION_SUCCESS = 0
   integer, parameter :: VALIDATION_WARNING = 1
   integer, parameter :: VALIDATION_ERROR = 2

   type :: validationresulttype
      integer :: status = validation_success
      character(len=512) :: message = ''
      character(len=32) :: species_name = ''
      character(len=32) :: emission_source = ''
   end type validationresulttype

   type :: emissionconfigvalidatortype
      private

      ! Validation settings
      logical :: strict_mode = .false.
      logical :: require_mass_conservation = .true.
      real(fp) :: mass_conservation_tolerance = 1.0e-6_fp
      logical :: check_species_existence = .true.
      logical :: validate_units = .true.
      logical :: check_scale_factors = .true.

      ! Validation statistics
      integer :: n_warnings = 0
      integer :: n_errors = 0
      integer :: n_species_validated = 0

   contains
      procedure :: validate_emission_config
      procedure :: validate_species_mapping
      procedure :: validate_emission_sources
      procedure :: validate_vertical_distributions
      procedure :: validate_mass_conservation

      ! Utility methods
      procedure :: reset_counters
      procedure :: get_validation_summary
      procedure :: set_validation_mode

      ! Private validation helpers
      procedure, private :: validate_single_species
      procedure, private :: check_species_in_mechanism
      procedure, private :: validate_scale_factors
   end type emissionconfigvalidatortype

contains

   subroutine validate_emission_config(this, config_file, container, results, rc)
      class(EmissionConfigValidatorType), intent(inout) :: this
      character(len=*), intent(in) :: config_file
      type(StateManagerType), intent(inout) :: container
      type(ValidationResultType), allocatable, intent(out) :: results(:)
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      character(len=256) :: message
      integer :: max_results = 1000
      integer :: n_results = 0

      rc = cc_success
      error_mgr => container%get_error_manager()

      ! Reset validation counters
      call this%reset_counters()

      ! Allocate results array
      allocate(results(max_results))

      write(message, '(A,A)') 'Starting validation of emission config: ', trim(config_file)
      print *, trim(message)

      ! TODO: Parse YAML configuration file
      ! For now, we'll validate the example configuration

      ! Validate species mappings
      call this%validate_species_mapping(container, results, n_results, rc)
      if (rc /= cc_success) return

      ! Validate emission sources
      call this%validate_emission_sources(results, n_results, rc)
      if (rc /= cc_success) return

      ! Validate vertical distributions
      call this%validate_vertical_distributions(results, n_results, rc)
      if (rc /= cc_success) return

      ! Resize results array to actual size
      if (n_results > 0) then
         results = results(1:n_results)
      else
         deallocate(results)
         allocate(results(0))
      end if

      ! Report validation summary
      call this%get_validation_summary(message)
      print *, trim(message)

   end subroutine validate_emission_config

   subroutine validate_species_mapping(this, container, results, n_results, rc)
      class(EmissionConfigValidatorType), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      type(ValidationResultType), intent(inout) :: results(:)
      integer, intent(inout) :: n_results
      integer, intent(out) :: rc

      ! Example validation for common species
      character(len=32), parameter :: test_species(5) = &
         ['NO  ', 'SO2 ', 'CO  ', 'ISOP', 'NH3 ']
      integer :: i

      rc = cc_success

      do i = 1, size(test_species)
         call this%validate_single_species(trim(test_species(i)), container, &
            results, n_results, rc)
         if (rc /= cc_success) return
      end do

   end subroutine validate_species_mapping

   subroutine validate_single_species(this, species_name, container, &
      results, n_results, rc)
      class(EmissionConfigValidatorType), intent(inout) :: this
      character(len=*), intent(in) :: species_name
      type(StateManagerType), intent(inout) :: container
      type(ValidationResultType), intent(inout) :: results(:)
      integer, intent(inout) :: n_results
      integer, intent(out) :: rc

      type(ValidationResultType) :: result
      logical :: species_exists
      character(len=256) :: message

      rc = cc_success
      this%n_species_validated = this%n_species_validated + 1

      ! Check if species exists in chemical mechanism
      call this%check_species_in_mechanism(species_name, container, species_exists, rc)
      if (rc /= cc_success) return

      if (.not. species_exists) then
         result%status = validation_error
         result%species_name = species_name
         write(result%message, '(A,A,A)') 'Species "', trim(species_name), &
            '" not found in chemical mechanism'

         n_results = n_results + 1
         results(n_results) = result
         this%n_errors = this%n_errors + 1

         if (this%strict_mode) then
            rc = error_invalid_config
            return
         end if
      end if

      ! TODO: Add more validations for this species
      ! - Scale factors
      ! - Units compatibility
      ! - Vertical distribution validity

   end subroutine validate_single_species

   subroutine check_species_in_mechanism(this, species_name, container, &
      species_exists, rc)
      class(EmissionConfigValidatorType), intent(in) :: this
      character(len=*), intent(in) :: species_name
      type(StateManagerType), intent(inout) :: container
      logical, intent(out) :: species_exists
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state
      integer :: species_idx

      rc = cc_success
      species_exists = .false.

      if (.not. this%check_species_existence) then
         species_exists = .true.  ! Skip check if disabled
         return
      end if

      chem_state => container%get_chem_state_ptr()

      species_idx = chem_state%find_species(species_name)
      species_exists = (species_idx > 0)

      ! Reset rc since this is just a check
      rc = cc_success

   end subroutine check_species_in_mechanism

   subroutine validate_scale_factors(this, scale_factors, n_factors, &
      species_name, results, n_results, rc)
      class(EmissionConfigValidatorType), intent(inout) :: this
      real(fp), intent(in) :: scale_factors(:)
      integer, intent(in) :: n_factors
      character(len=*), intent(in) :: species_name
      type(ValidationResultType), intent(inout) :: results(:)
      integer, intent(inout) :: n_results
      integer, intent(out) :: rc

      type(ValidationResultType) :: result
      real(fp) :: total_scale
      character(len=256) :: message

      rc = cc_success

      if (.not. this%check_scale_factors) return

      total_scale = sum(scale_factors(1:n_factors))

      if (this%require_mass_conservation) then
         if (abs(total_scale - 1.0_fp) > this%mass_conservation_tolerance) then
            result%status = validation_warning
            result%species_name = species_name
            write(result%message, '(A,A,A,F8.5,A)') 'Scale factors for "', &
               trim(species_name), '" sum to ', total_scale, &
               ' (expected 1.0 for mass conservation)'

            n_results = n_results + 1
            results(n_results) = result
            this%n_warnings = this%n_warnings + 1
         end if
      end if

      ! Check for negative scale factors
      if (any(scale_factors(1:n_factors) < 0.0_fp)) then
         result%status = validation_error
         result%species_name = species_name
         write(result%message, '(A,A,A)') 'Negative scale factors found for "', &
            trim(species_name), '"'

         n_results = n_results + 1
         results(n_results) = result
         this%n_errors = this%n_errors + 1

         if (this%strict_mode) then
            rc = error_invalid_config
            return
         end if
      end if

   end subroutine validate_scale_factors

   subroutine validate_emission_sources(this, results, n_results, rc)
      class(EmissionConfigValidatorType), intent(inout) :: this
      type(ValidationResultType), intent(inout) :: results(:)
      integer, intent(inout) :: n_results
      integer, intent(out) :: rc

      ! TODO: Validate that emission source files exist and are accessible
      rc = cc_success

   end subroutine validate_emission_sources

   subroutine validate_vertical_distributions(this, results, n_results, rc)
      class(EmissionConfigValidatorType), intent(inout) :: this
      type(ValidationResultType), intent(inout) :: results(:)
      integer, intent(inout) :: n_results
      integer, intent(out) :: rc

      ! TODO: Validate vertical distribution parameters
      rc = cc_success

   end subroutine validate_vertical_distributions

   subroutine validate_mass_conservation(this, results, n_results, rc)
      class(EmissionConfigValidatorType), intent(inout) :: this
      type(ValidationResultType), intent(inout) :: results(:)
      integer, intent(inout) :: n_results
      integer, intent(out) :: rc

      ! TODO: Check mass conservation across all mappings
      rc = cc_success

   end subroutine validate_mass_conservation

   subroutine reset_counters(this)
      class(EmissionConfigValidatorType), intent(inout) :: this

      this%n_warnings = 0
      this%n_errors = 0
      this%n_species_validated = 0

   end subroutine reset_counters

   subroutine get_validation_summary(this, summary)
      class(EmissionConfigValidatorType), intent(in) :: this
      character(len=*), intent(out) :: summary

      write(summary, '(A,I0,A,I0,A,I0,A)') &
         'Emission config validation complete: ', &
         this%n_species_validated, ' species validated, ', &
         this%n_warnings, ' warnings, ', &
         this%n_errors, ' errors'

   end subroutine get_validation_summary

   subroutine set_validation_mode(this, strict_mode, require_mass_conservation, &
      check_species_existence)
      class(EmissionConfigValidatorType), intent(inout) :: this
      logical, intent(in), optional :: strict_mode
      logical, intent(in), optional :: require_mass_conservation
      logical, intent(in), optional :: check_species_existence

      if (present(strict_mode)) this%strict_mode = strict_mode
      if (present(require_mass_conservation)) &
         this%require_mass_conservation = require_mass_conservation
      if (present(check_species_existence)) &
         this%check_species_existence = check_species_existence

   end subroutine set_validation_mode

end module emissionconfigvalidator_mod
```


