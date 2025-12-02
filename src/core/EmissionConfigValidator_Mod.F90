!> \file EmissionConfigValidator_Mod.F90
!! \brief Configuration validation for emission species mapping
!! \ingroup process_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 1.0
!!
!! This module provides comprehensive validation for emission configuration files,
!! ensuring compatibility between emission inputs and chemical mechanisms.
!!
module EmissionConfigValidator_Mod
   use Precision_Mod
   use Error_Mod
   use StateManager_Mod, only : StateManagerType
   use ChemState_Mod, only : ChemStateType
   ! use logging_mod, only : log_message, LOG_INFO

   implicit none
   private

   public :: EmissionConfigValidatorType
   public :: ValidationResultType
   public :: VALIDATION_SUCCESS, VALIDATION_WARNING, VALIDATION_ERROR

   ! Validation result constants
   integer, parameter :: VALIDATION_SUCCESS = 0
   integer, parameter :: VALIDATION_WARNING = 1
   integer, parameter :: VALIDATION_ERROR = 2

   !> Validation result type
   type :: ValidationResultType
      integer :: status = VALIDATION_SUCCESS
      character(len=512) :: message = ''
      character(len=32) :: species_name = ''
      character(len=32) :: emission_source = ''
   end type ValidationResultType

   !> Configuration validator type
   type :: EmissionConfigValidatorType
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
      !> \copydoc validate_emission_config
      procedure :: validate_emission_config
      !> \copydoc validate_species_mapping
      procedure :: validate_species_mapping
      !> \copydoc validate_emission_sources
      procedure :: validate_emission_sources
      !> \copydoc validate_vertical_distributions
      procedure :: validate_vertical_distributions
      !> \copydoc validate_mass_conservation
      procedure :: validate_mass_conservation

      ! Utility methods
      !> \copydoc reset_counters
      procedure :: reset_counters
      !> \copydoc get_validation_summary
      procedure :: get_validation_summary
      !> \copydoc set_validation_mode
      procedure :: set_validation_mode

      ! Private validation helpers
      !> \copydoc validate_single_species
      procedure, private :: validate_single_species
      !> \copydoc check_species_in_mechanism
      procedure, private :: check_species_in_mechanism
      !> \copydoc validate_scale_factors
      procedure, private :: validate_scale_factors
   end type EmissionConfigValidatorType

contains

   !> Validate complete emission configuration
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

      rc = CC_SUCCESS
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
      if (rc /= CC_SUCCESS) return

      ! Validate emission sources
      call this%validate_emission_sources(results, n_results, rc)
      if (rc /= CC_SUCCESS) return

      ! Validate vertical distributions
      call this%validate_vertical_distributions(results, n_results, rc)
      if (rc /= CC_SUCCESS) return

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

   !> Validate species mapping configuration
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

      rc = CC_SUCCESS

      do i = 1, size(test_species)
         call this%validate_single_species(trim(test_species(i)), container, &
            results, n_results, rc)
         if (rc /= CC_SUCCESS) return
      end do

   end subroutine validate_species_mapping

   !> Validate single species mapping
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

      rc = CC_SUCCESS
      this%n_species_validated = this%n_species_validated + 1

      ! Check if species exists in chemical mechanism
      call this%check_species_in_mechanism(species_name, container, species_exists, rc)
      if (rc /= CC_SUCCESS) return

      if (.not. species_exists) then
         result%status = VALIDATION_ERROR
         result%species_name = species_name
         write(result%message, '(A,A,A)') 'Species "', trim(species_name), &
            '" not found in chemical mechanism'

         n_results = n_results + 1
         results(n_results) = result
         this%n_errors = this%n_errors + 1

         if (this%strict_mode) then
            rc = ERROR_INVALID_CONFIG
            return
         end if
      end if

      ! TODO: Add more validations for this species
      ! - Scale factors
      ! - Units compatibility
      ! - Vertical distribution validity

   end subroutine validate_single_species

   !> Check if species exists in chemical mechanism
   subroutine check_species_in_mechanism(this, species_name, container, &
      species_exists, rc)
      class(EmissionConfigValidatorType), intent(in) :: this
      character(len=*), intent(in) :: species_name
      type(StateManagerType), intent(inout) :: container
      logical, intent(out) :: species_exists
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state
      integer :: species_idx

      rc = CC_SUCCESS
      species_exists = .false.

      if (.not. this%check_species_existence) then
         species_exists = .true.  ! Skip check if disabled
         return
      end if

      chem_state => container%get_chem_state_ptr()

      species_idx = chem_state%find_species(species_name)
      species_exists = (species_idx > 0)

      ! Reset rc since this is just a check
      rc = CC_SUCCESS

   end subroutine check_species_in_mechanism

   !> Validate scale factors for mass conservation
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

      rc = CC_SUCCESS

      if (.not. this%check_scale_factors) return

      total_scale = sum(scale_factors(1:n_factors))

      if (this%require_mass_conservation) then
         if (abs(total_scale - 1.0_fp) > this%mass_conservation_tolerance) then
            result%status = VALIDATION_WARNING
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
         result%status = VALIDATION_ERROR
         result%species_name = species_name
         write(result%message, '(A,A,A)') 'Negative scale factors found for "', &
            trim(species_name), '"'

         n_results = n_results + 1
         results(n_results) = result
         this%n_errors = this%n_errors + 1

         if (this%strict_mode) then
            rc = ERROR_INVALID_CONFIG
            return
         end if
      end if

   end subroutine validate_scale_factors

   !> Validate emission sources
   subroutine validate_emission_sources(this, results, n_results, rc)
      class(EmissionConfigValidatorType), intent(inout) :: this
      type(ValidationResultType), intent(inout) :: results(:)
      integer, intent(inout) :: n_results
      integer, intent(out) :: rc

      ! TODO: Validate that emission source files exist and are accessible
      rc = CC_SUCCESS

   end subroutine validate_emission_sources

   !> Validate vertical distribution configurations
   subroutine validate_vertical_distributions(this, results, n_results, rc)
      class(EmissionConfigValidatorType), intent(inout) :: this
      type(ValidationResultType), intent(inout) :: results(:)
      integer, intent(inout) :: n_results
      integer, intent(out) :: rc

      ! TODO: Validate vertical distribution parameters
      rc = CC_SUCCESS

   end subroutine validate_vertical_distributions

   !> Validate mass conservation for all species
   subroutine validate_mass_conservation(this, results, n_results, rc)
      class(EmissionConfigValidatorType), intent(inout) :: this
      type(ValidationResultType), intent(inout) :: results(:)
      integer, intent(inout) :: n_results
      integer, intent(out) :: rc

      ! TODO: Check mass conservation across all mappings
      rc = CC_SUCCESS

   end subroutine validate_mass_conservation

   !> Reset validation counters
   subroutine reset_counters(this)
      class(EmissionConfigValidatorType), intent(inout) :: this

      this%n_warnings = 0
      this%n_errors = 0
      this%n_species_validated = 0

   end subroutine reset_counters

   !> Get validation summary message
   subroutine get_validation_summary(this, summary)
      class(EmissionConfigValidatorType), intent(in) :: this
      character(len=*), intent(out) :: summary

      write(summary, '(A,I0,A,I0,A,I0,A)') &
         'Emission config validation complete: ', &
         this%n_species_validated, ' species validated, ', &
         this%n_warnings, ' warnings, ', &
         this%n_errors, ' errors'

   end subroutine get_validation_summary

   !> Set validation mode (strict or permissive)
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

end module EmissionConfigValidator_Mod
