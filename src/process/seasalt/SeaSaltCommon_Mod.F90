!> \file SeaSaltCommon_Mod.F90
!! \brief Common types and utilities for seasalt process
!!
!! This module defines the configuration types used by the
!! seasalt process and its schemes.
!!
!! Generated on: 2025-11-14T23:01:21.563867
!! Author: Barry Baker & Wei Li
!! Version: 1.0.0

module SeaSaltCommon_Mod

   use precision_mod, only: fp
   ! use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, CC_Error, CC_Warning, ErrorManagerType, &
      ERROR_INVALID_CONFIG, ERROR_INVALID_STATE, ERROR_NOT_FOUND
   use ConfigManager_Mod, only: ConfigManagerType  ! ConfigManager integration
   use StateManager_Mod, only: StateManagerType  ! Add StateManager integration

   implicit none
   private

   ! Export types
   public :: SeaSaltProcessConfig  ! New unified process config
   public :: SeaSaltConfig
   public :: SeaSaltSchemeGONG97Config
   public :: SeaSaltSchemeGONG03Config
   public :: SeaSaltSchemeGEOS12Config

   ! Export utility functions
   public :: int_to_string

   !> Main configuration type for seasalt process
   type :: SeaSaltConfig

      ! Process settings
      character(len=32) :: scheme = 'gong97'
      logical :: is_active = .true.
      logical :: diagnostics = .false.  ! Diagnostic switch

      ! Diagnostic species configuration
      integer :: n_diagnostic_species = 0
      character(len=32), allocatable :: diagnostic_species(:)  ! User-defined species for diagnostics
      integer, allocatable :: diagnostic_species_id(:)  ! Indices mapping diagnostic_species to species_names
      real(fp) :: dt_min = 1.0_fp     ! Minimum time step (seconds)
      real(fp) :: dt_max = 3600.0_fp  ! Maximum time step (seconds)

      ! Species configuration
      integer :: n_species = 0
      character(len=32), allocatable :: species_names(:)
      integer, allocatable :: species_indices(:)  ! Indices of seasalt species in ChemState



      ! Species properties
      real(fp), allocatable :: species_density(:)      ! density for each species
      real(fp), allocatable :: species_lower_radius(:)      ! lower_radius for each species
      real(fp), allocatable :: species_radius(:)      ! radius for each species
      real(fp), allocatable :: species_upper_radius(:)      ! upper_radius for each species

      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   contains
      procedure, public :: validate => validate_seasalt_config
      procedure, public :: finalize => finalize_seasalt_config
      procedure, public :: print_summary => print_seasalt_config_summary
   end type SeaSaltConfig

   !> Configuration type for gong97 scheme
   type :: SeaSaltSchemeGONG97Config

      ! Scheme metadata
      character(len=64) :: scheme_name = 'gong97'
      character(len=256) :: description = 'Gong 1997 sea salt emission scheme'
      character(len=64) :: author = 'Barry Baker'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .false.  ! Surface-only processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Emission scale factor
      logical :: weibull_flag = .false.  ! Apply Weibull distribution for particle size

      ! Required meteorological fields
      integer :: n_required_met_fields = 5
      character(len=32) :: required_met_fields(5)

   contains
      procedure, public :: validate => validate_gong97_config
      procedure, public :: finalize => finalize_gong97_config
   end type SeaSaltSchemeGONG97Config

   ! gong97 scheme uses local variables only - no persistent state type needed

   !> Configuration type for gong03 scheme
   type :: SeaSaltSchemeGONG03Config

      ! Scheme metadata
      character(len=64) :: scheme_name = 'gong03'
      character(len=256) :: description = 'Gong 2003 sea salt emission scheme with improved sub- and super-micron treatment'
      character(len=64) :: author = 'Barry Baker'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .false.  ! Surface-only processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Emission scale factor
      logical :: weibull_flag = .false.  ! Apply Weibull distribution for particle size

      ! Required meteorological fields
      integer :: n_required_met_fields = 5
      character(len=32) :: required_met_fields(5)

   contains
      procedure, public :: validate => validate_gong03_config
      procedure, public :: finalize => finalize_gong03_config
   end type SeaSaltSchemeGONG03Config

   ! gong03 scheme uses local variables only - no persistent state type needed

   !> Configuration type for geos12 scheme
   type :: SeaSaltSchemeGEOS12Config

      ! Scheme metadata
      character(len=64) :: scheme_name = 'geos12'
      character(len=256) :: description = 'GEOS-Chem 2012 sea salt emission scheme with observational constraints'
      character(len=64) :: author = 'Barry Baker'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .false.  ! Surface-only processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Emission scale factor

      ! Required meteorological fields
      integer :: n_required_met_fields = 4
      character(len=32) :: required_met_fields(4)

   contains
      procedure, public :: validate => validate_geos12_config
      procedure, public :: finalize => finalize_geos12_config
   end type SeaSaltSchemeGEOS12Config

   ! geos12 scheme uses local variables only - no persistent state type needed


   !> Unified process configuration type that bridges ConfigManager and process-specific configs
   !! This is the main configuration type that ProcessInterface should use
   type :: SeaSaltProcessConfig

      ! Process metadata
      character(len=64) :: process_name = 'seasalt'
      character(len=16) :: process_version = '1.0.0'
      logical :: is_active = .true.

      ! Process-specific configuration (delegate to SeaSaltConfig)
      type(SeaSaltConfig) :: seasalt_config

      ! Scheme configurations
      type(SeaSaltSchemeGONG97Config) :: gong97_config
      type(SeaSaltSchemeGONG03Config) :: gong03_config
      type(SeaSaltSchemeGEOS12Config) :: geos12_config

   contains
      procedure, public :: load_from_config => seasalt_process_load_config
      procedure, public :: load_species_from_chem_state => load_species_from_chem_state
      procedure, public :: validate => seasalt_process_validate
      procedure, public :: finalize => seasalt_process_finalize
      procedure, public :: get_active_scheme_config => get_active_scheme_config
      procedure, public :: load_gong97_config
      procedure, public :: load_gong03_config
      procedure, public :: load_geos12_config
      procedure, public :: map_diagnostic_species_indices
   end type SeaSaltProcessConfig

contains

   !> Validate seasalt configuration
   subroutine validate_seasalt_config(this, error_handler)
      class(SeaSaltConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      character(len=256) :: error_msg
      integer :: rc

      ! Validate time step bounds
      if (this%dt_min <= 0.0_fp) then
         call error_handler%report_error(ERROR_INVALID_CONFIG, &
            "Minimum time step must be positive", rc)
         return
      end if

      if (this%dt_max < this%dt_min) then
         call error_handler%report_error(ERROR_INVALID_CONFIG, &
            "Maximum time step must be >= minimum time step", rc)
         return
      end if

      ! Validate active scheme(s)
      ! Validate scheme
      if (trim(this%scheme) /= 'gong97' .and. &
         trim(this%scheme) /= 'gong03' .and. &
         trim(this%scheme) /= 'geos12' .and. &
         .true.) then
         write(error_msg, '(A)') "Invalid scheme: " // trim(this%scheme)
         call error_handler%report_error(ERROR_INVALID_CONFIG, error_msg, rc)
         return
      end if

   end subroutine validate_seasalt_config

   !> Print configuration summary
   subroutine print_seasalt_config_summary(this)
      class(SeaSaltConfig), intent(in) :: this

      write(*, '(A)') "=== SeaSalt Process Configuration ==="
      write(*, '(A,A)') "  Active scheme: ", trim(this%scheme)
      write(*, '(A,I0)') "  Number of species: ", this%n_species
      write(*, '(A,F0.1,A)') "  Minimum time step: ", this%dt_min, " s"
      write(*, '(A,F0.1,A)') "  Maximum time step: ", this%dt_max, " s"
      write(*, '(A,L1)') "  Output diagnostics: ", this%output_diagnostics
      write(*, '(A)') "============================================="

   end subroutine print_seasalt_config_summary

   !> Finalize seasalt configuration
   subroutine finalize_seasalt_config(this)
      class(SeaSaltConfig), intent(inout) :: this

      ! Deallocate species names array
      if (allocated(this%species_names)) then
         deallocate(this%species_names)
      end if

      ! Deallocate species indices array
      if (allocated(this%species_indices)) then
         deallocate(this%species_indices)
      end if

      ! Deallocate species properties arrays
      if (allocated(this%species_density)) then
         deallocate(this%species_density)
      end if
      if (allocated(this%species_lower_radius)) then
         deallocate(this%species_lower_radius)
      end if
      if (allocated(this%species_radius)) then
         deallocate(this%species_radius)
      end if
      if (allocated(this%species_upper_radius)) then
         deallocate(this%species_upper_radius)
      end if


      ! Deallocate diagnostic species array
      if (allocated(this%diagnostic_species)) then
         deallocate(this%diagnostic_species)
      end if

      ! Deallocate diagnostic species indices array
      if (allocated(this%diagnostic_species_id)) then
         deallocate(this%diagnostic_species_id)
      end if

   end subroutine finalize_seasalt_config

   !> Validate gong97 scheme configuration
   subroutine validate_gong97_config(this, error_handler)
      class(SeaSaltSchemeGONG97Config), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_gong97_config

   !> Finalize gong97 scheme configuration
   subroutine finalize_gong97_config(this)
      class(SeaSaltSchemeGONG97Config), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_gong97_config


   !> Validate gong03 scheme configuration
   subroutine validate_gong03_config(this, error_handler)
      class(SeaSaltSchemeGONG03Config), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_gong03_config

   !> Finalize gong03 scheme configuration
   subroutine finalize_gong03_config(this)
      class(SeaSaltSchemeGONG03Config), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_gong03_config


   !> Validate geos12 scheme configuration
   subroutine validate_geos12_config(this, error_handler)
      class(SeaSaltSchemeGEOS12Config), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_geos12_config

   !> Finalize geos12 scheme configuration
   subroutine finalize_geos12_config(this)
      class(SeaSaltSchemeGEOS12Config), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_geos12_config



   !> Convert integer to string (utility function)
   function int_to_string(int_val) result(str_val)
      integer, intent(in) :: int_val
      character(len=32) :: str_val

      write(str_val, '(I0)') int_val
      str_val = adjustl(str_val)

   end function int_to_string

   !> Load unified process configuration from ConfigManager
   !! This is the main function that ProcessInterface.parse_process_config should call
   !! Process reads its configuration directly from the master YAML via ConfigManager
   subroutine seasalt_process_load_config(this, config_manager, error_handler)
      class(SeaSaltProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      character(len=256) :: scheme_name
      integer :: ierr, rc

      ! Process reads directly from master YAML structure: processes.seasalt
      ! ConfigManager provides generic YAML access, process handles its own configuration

      ! Load process metadata
      call config_manager%get_string("processes/seasalt/name", this%process_name, rc, "seasalt")
      if (rc /= CC_SUCCESS) this%process_name = "seasalt"  ! default

      call config_manager%get_string("processes/seasalt/version", this%process_version, rc, "1.0.0")
      if (rc /= CC_SUCCESS) this%process_version = "1.0.0"  ! default

      call config_manager%get_logical("processes/seasalt/activate", this%is_active, rc, .true.)
      if (rc /= CC_SUCCESS) this%is_active = .true.  ! default

      ! Load process-specific configuration directly from master YAML
      call config_manager%get_string("processes/seasalt/scheme", this%seasalt_config%scheme, rc, "gong97")
      if (rc /= CC_SUCCESS) then
         call error_handler%report_error(ERROR_INVALID_CONFIG, &
            "Missing required 'scheme' in processes/seasalt configuration", rc)
         return
      end if

      ! Load diagnostic switch
      call config_manager%get_logical("processes/seasalt/diagnostics", this%seasalt_config%diagnostics, rc, .false.)
      if (rc /= CC_SUCCESS) this%seasalt_config%diagnostics = .false.  ! Default

      ! Load diagnostic species list
      call config_manager%get_array("processes/seasalt/diag_species", this%seasalt_config%diagnostic_species, &
         rc, default_values=["All"])
      if (rc /= CC_SUCCESS) then
         ! Default to all species if not specified
         allocate(this%seasalt_config%diagnostic_species(1))
         this%seasalt_config%diagnostic_species(1) = "All"
         this%seasalt_config%n_diagnostic_species = 1
      else
         ! Set the count based on the returned array size
         if (allocated(this%seasalt_config%diagnostic_species)) then
            this%seasalt_config%n_diagnostic_species = size(this%seasalt_config%diagnostic_species)
         else
            this%seasalt_config%n_diagnostic_species = 0
         end if
      end if

      ! Species configuration is loaded from ChemState in load_species_from_chem_state
      ! The species come from the master species YAML file (CATChem_species.yml)
      ! and are filtered by is_seasalt property


      ! Load scheme-specific configuration from master YAML
      scheme_name = trim(this%seasalt_config%scheme)
      select case (scheme_name)
       case ('gong97')
         call this%load_gong97_config(config_manager, error_handler)
       case ('gong03')
         call this%load_gong03_config(config_manager, error_handler)
       case ('geos12')
         call this%load_geos12_config(config_manager, error_handler)
       case default
         call error_handler%report_error(ERROR_INVALID_STATE, &
            "Unknown seasalt scheme: " // trim(scheme_name), rc)
         return
      end select

   end subroutine seasalt_process_load_config


   !> Load species from ChemState
   !! This function is used for dynamic species discovery (by_metadata or all_species)
   !! For 'all_species' mode: loads all species using nSpecies and SpeciesIndex
   !! For 'by_metadata' mode: loads species by type using nSpeciesSeaSalt and SeaSaltIndex
   subroutine load_species_from_chem_state(this, chem_state, error_handler)
      use ChemState_Mod, only: ChemStateType

      class(SeaSaltProcessConfig), intent(inout) :: this
      type(ChemStateType), pointer, intent(in) :: chem_state
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: i, rc
      integer :: species_idx

      if (.not. associated(chem_state)) then
         call error_handler%report_error(ERROR_INVALID_STATE, &
            "ChemState not associated in load_species_from_chem_state", rc)
         return
      end if

      ! by_metadata mode: Load species by type from ChemState using dynamic metadata flag mapping
      ! Dynamic mapping: is_seasalt -> nSpeciesSeaSalt and SeaSaltIndex
      this%seasalt_config%n_species = chem_state%nSpeciesSeaSalt

      if (this%seasalt_config%n_species <= 0) then
         call error_handler%report_error(ERROR_INVALID_STATE, &
            "No seasalt species found in ChemState", rc)
         return
      end if

      ! Check if SeaSaltIndex is allocated and has correct size
      if (.not. allocated(chem_state%SeaSaltIndex)) then
         call error_handler%report_error(ERROR_INVALID_STATE, &
            "SeaSaltIndex not allocated in ChemState", rc)
         return
      end if

      if (size(chem_state%SeaSaltIndex) < this%seasalt_config%n_species) then
         call error_handler%report_error(ERROR_INVALID_STATE, &
            "SeaSaltIndex size inconsistent with nSpeciesSeaSalt", rc)
         return
      end if

      ! Deallocate existing arrays if allocated
      if (allocated(this%seasalt_config%species_names)) then
         deallocate(this%seasalt_config%species_names)
      end if
      if (allocated(this%seasalt_config%species_indices)) then
         deallocate(this%seasalt_config%species_indices)
      end if

      ! Allocate arrays
      allocate(this%seasalt_config%species_names(this%seasalt_config%n_species))
      allocate(this%seasalt_config%species_indices(this%seasalt_config%n_species))

      ! Allocate species properties arrays
      allocate(this%seasalt_config%species_density(this%seasalt_config%n_species))
      allocate(this%seasalt_config%species_lower_radius(this%seasalt_config%n_species))
      allocate(this%seasalt_config%species_radius(this%seasalt_config%n_species))
      allocate(this%seasalt_config%species_upper_radius(this%seasalt_config%n_species))

      ! by_metadata mode: Copy indices from metadata-specific index array using dynamic mapping
      ! Dynamic mapping: is_seasalt -> SeaSaltIndex
      this%seasalt_config%species_indices(1:this%seasalt_config%n_species) = &
         chem_state%SeaSaltIndex(1:this%seasalt_config%n_species)

      ! Get species names using the indices
      do i = 1, this%seasalt_config%n_species
         if (this%seasalt_config%species_indices(i) > 0 .and. &
            this%seasalt_config%species_indices(i) <= size(chem_state%SpeciesNames)) then
            this%seasalt_config%species_names(i) = &
               trim(chem_state%SpeciesNames(this%seasalt_config%species_indices(i)))
         else
            call error_handler%report_error(ERROR_INVALID_STATE, &
               "Invalid species index in species index array", rc)
            return
         end if
      end do

      ! Load species properties from ChemState
      do i = 1, this%seasalt_config%n_species
         species_idx = this%seasalt_config%species_indices(i)
         this%seasalt_config%species_density(i) = chem_state%ChemSpecies(species_idx)%density
         this%seasalt_config%species_lower_radius(i) = chem_state%ChemSpecies(species_idx)%lower_radius
         this%seasalt_config%species_radius(i) = chem_state%ChemSpecies(species_idx)%radius
         this%seasalt_config%species_upper_radius(i) = chem_state%ChemSpecies(species_idx)%upper_radius
      end do

   end subroutine load_species_from_chem_state


   !> Load gong97 scheme configuration from master YAML
   subroutine load_gong97_config(this, config_manager, error_handler)
      class(SeaSaltProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: ierr, rc

      ! Load scheme parameters directly from processes/seasalt/gong97/ in master YAML
      call config_manager%get_real("processes/seasalt/gong97/scale_factor", &
         this%gong97_config%scale_factor, rc, 1.0_fp)
      if (rc /= CC_SUCCESS) this%gong97_config%scale_factor = 1.0_fp
      call config_manager%get_logical("processes/seasalt/gong97/weibull_flag", &
         this%gong97_config%weibull_flag, rc, .false.)
      if (rc /= CC_SUCCESS) this%gong97_config%weibull_flag = .false.


   end subroutine load_gong97_config

   !> Load gong03 scheme configuration from master YAML
   subroutine load_gong03_config(this, config_manager, error_handler)
      class(SeaSaltProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: ierr, rc

      ! Load scheme parameters directly from processes/seasalt/gong03/ in master YAML
      call config_manager%get_real("processes/seasalt/gong03/scale_factor", &
         this%gong03_config%scale_factor, rc, 1.0_fp)
      if (rc /= CC_SUCCESS) this%gong03_config%scale_factor = 1.0_fp
      call config_manager%get_logical("processes/seasalt/gong03/weibull_flag", &
         this%gong03_config%weibull_flag, rc, .false.)
      if (rc /= CC_SUCCESS) this%gong03_config%weibull_flag = .false.


   end subroutine load_gong03_config

   !> Load geos12 scheme configuration from master YAML
   subroutine load_geos12_config(this, config_manager, error_handler)
      class(SeaSaltProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: ierr, rc

      ! Load scheme parameters directly from processes/seasalt/geos12/ in master YAML
      call config_manager%get_real("processes/seasalt/geos12/scale_factor", &
         this%geos12_config%scale_factor, rc, 1.0_fp)
      if (rc /= CC_SUCCESS) this%geos12_config%scale_factor = 1.0_fp


   end subroutine load_geos12_config


   !> Validate unified process configuration
   subroutine seasalt_process_validate(this, state_manager, error_handler)
      class(SeaSaltProcessConfig), intent(inout) :: this
      type(StateManagerType), intent(in) :: state_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      ! Validate main config
      call this%seasalt_config%validate(error_handler)

      ! Validate scheme-specific config
      select case (trim(this%seasalt_config%scheme))
       case ('gong97')
         call this%gong97_config%validate(error_handler)
       case ('gong03')
         call this%gong03_config%validate(error_handler)
       case ('geos12')
         call this%geos12_config%validate(error_handler)
      end select

   end subroutine seasalt_process_validate

   !> Finalize unified process configuration
   subroutine seasalt_process_finalize(this)
      class(SeaSaltProcessConfig), intent(inout) :: this

      call this%seasalt_config%finalize()
      call this%gong97_config%finalize()
      call this%gong03_config%finalize()
      call this%geos12_config%finalize()

   end subroutine seasalt_process_finalize

   !> Get active scheme configuration (polymorphic return)
   function get_active_scheme_config(this) result(scheme_config)
      class(SeaSaltProcessConfig), intent(in) :: this
      class(*), allocatable :: scheme_config

      select case (trim(this%seasalt_config%scheme))
       case ('gong97')
         allocate(scheme_config, source=this%gong97_config)
       case ('gong03')
         allocate(scheme_config, source=this%gong03_config)
       case ('geos12')
         allocate(scheme_config, source=this%geos12_config)
       case default
         ! Return null
      end select

   end function get_active_scheme_config

   !> Map diagnostic species names to indices in the species_names array
   !! This function creates the diagnostic_species_id array that maps each diagnostic species
   !! to its corresponding index in the full species_names array
   subroutine map_diagnostic_species_indices(this, error_handler)
      class(SeaSaltProcessConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: i, j, rc
      character(len=256) :: error_msg
      logical :: found_species

      ! Only proceed if diagnostic species are defined
      if (this%seasalt_config%n_diagnostic_species == 0) return

      ! Handle "All" case - map all available species
      if (this%seasalt_config%n_diagnostic_species == 1 .and. &
         trim(this%seasalt_config%diagnostic_species(1)) == "All") then

         ! Deallocate and reallocate for all species
         if (allocated(this%seasalt_config%diagnostic_species_id)) deallocate(this%seasalt_config%diagnostic_species_id)
         allocate(this%seasalt_config%diagnostic_species_id(this%seasalt_config%n_species))
         if (allocated(this%seasalt_config%diagnostic_species)) deallocate(this%seasalt_config%diagnostic_species)
         allocate(this%seasalt_config%diagnostic_species(this%seasalt_config%n_species))
         this%seasalt_config%n_diagnostic_species = this%seasalt_config%n_species
         this%seasalt_config%diagnostic_species = this%seasalt_config%species_names

         ! Map all species indices (1:n_species)
         do i = 1, this%seasalt_config%n_species
            this%seasalt_config%diagnostic_species_id(i) = i
         end do

         return
      end if

      ! Allocate diagnostic species indices array
      if (allocated(this%seasalt_config%diagnostic_species_id)) deallocate(this%seasalt_config%diagnostic_species_id)
      allocate(this%seasalt_config%diagnostic_species_id(this%seasalt_config%n_diagnostic_species))

      ! Map each diagnostic species name to its index in species_names
      do i = 1, this%seasalt_config%n_diagnostic_species
         found_species = .false.

         do j = 1, this%seasalt_config%n_species
            if (trim(this%seasalt_config%diagnostic_species(i)) == trim(this%seasalt_config%species_names(j))) then
               this%seasalt_config%diagnostic_species_id(i) = j
               found_species = .true.
               exit
            end if
         end do

         if (.not. found_species) then
            write(error_msg, '(A,A,A)') "Diagnostic species '", &
               trim(this%seasalt_config%diagnostic_species(i)), &
               "' not found in process species list"
            call error_handler%report_error(ERROR_NOT_FOUND, error_msg, rc)
            return
         end if
      end do

   end subroutine map_diagnostic_species_indices

end module SeaSaltCommon_Mod
