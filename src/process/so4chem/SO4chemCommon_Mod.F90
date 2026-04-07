!> \file SO4chemCommon_Mod.F90
!! \brief Common types and utilities for so4chem process
!!
!! This module defines the configuration types used by the
!! so4chem process and its schemes.
!!
!! Generated on: 2026-03-03T18:15:44.287860
!! Author: Wei Li
!! Version: 1.0.0

module SO4chemCommon_Mod

   use precision_mod, only: fp
   ! use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, CC_Error, CC_Warning, ErrorManagerType, &
      ERROR_INVALID_CONFIG, ERROR_INVALID_STATE, ERROR_NOT_FOUND
   use ConfigManager_Mod, only: ConfigManagerType  ! ConfigManager integration
   use StateManager_Mod, only: StateManagerType  ! Add StateManager integration

   implicit none
   private

   ! Export types
   public :: SO4chemProcessConfig  ! New unified process config
   public :: SO4chemConfig
   public :: SO4chemSchemeGOCARTConfig

   ! Export utility functions
   public :: int_to_string

   !> Main configuration type for so4chem process
   type :: SO4chemConfig

      ! Process settings
      character(len=32) :: scheme = 'gocart'
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
      integer, allocatable :: species_indices(:)  ! Indices of so4chem species in ChemState



      ! Species properties
      real(fp), allocatable :: species_mw_g(:)      ! mw_g for each species

      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   contains
      procedure, public :: validate => validate_so4chem_config
      procedure, public :: finalize => finalize_so4chem_config
      procedure, public :: print_summary => print_so4chem_config_summary
   end type SO4chemConfig

   !> Configuration type for gocart scheme
   type :: SO4chemSchemeGOCARTConfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'gocart'
      character(len=256) :: description = 'GOCART SO2 to SO4 production scheme'
      character(len=64) :: author = 'Wei Li'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .true.  ! Full column processing

      ! Scheme parameters
      logical :: update_so2 = .true.  ! whether to update SO2 concentration based on chemical production/loss

      ! Required meteorological fields
      integer :: n_required_met_fields = 16
      character(len=32) :: required_met_fields(16)

   contains
      procedure, public :: validate => validate_gocart_config
      procedure, public :: finalize => finalize_gocart_config
   end type SO4chemSchemeGOCARTConfig

   !> Persistent state type for gocart scheme
   !! Contains variables that persist across time steps for each column
   type :: SO4chemGOCARTPersistentState
      logical :: firsttime  ! flag for first time step
      integer :: nymd_last  ! last day of H2O2 update
      integer :: nhms_last_recycle  ! last time step of H2O2 recycle
      real(fp), allocatable :: xh2o2_init(:)  ! H2O2 column initialization
   end type SO4chemGOCARTPersistentState


   !> Unified process configuration type that bridges ConfigManager and process-specific configs
   !! This is the main configuration type that ProcessInterface should use
   type :: SO4chemProcessConfig

      ! Process metadata
      character(len=64) :: process_name = 'so4chem'
      character(len=16) :: process_version = '1.0.0'
      logical :: is_active = .true.

      ! Process-specific configuration (delegate to SO4chemConfig)
      type(SO4chemConfig) :: so4chem_config

      ! Scheme configurations
      type(SO4chemSchemeGOCARTConfig) :: gocart_config

      ! Persistent state arrays for column processing
      type(SO4chemGocartPersistentState), allocatable :: gocart_persistent_state(:)  ! Per-column state for gocart scheme
      integer :: total_columns = 0  ! Total number of columns for state allocation

   contains
      procedure, public :: load_from_config => so4chem_process_load_config
      procedure, public :: load_species_from_list => load_species_from_list
      procedure, public :: validate => so4chem_process_validate
      procedure, public :: finalize => so4chem_process_finalize
      procedure, public :: initialize_persistent_state => so4chem_initialize_persistent_state
      procedure, public :: get_active_scheme_config => get_active_scheme_config
      procedure, public :: load_gocart_config
      procedure, public :: map_diagnostic_species_indices
   end type SO4chemProcessConfig

contains

   !> Validate so4chem configuration
   subroutine validate_so4chem_config(this, error_handler)
      class(SO4chemConfig), intent(inout) :: this
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
      if (trim(this%scheme) /= 'gocart' .and. &
         .true.) then
         write(error_msg, '(A)') "Invalid scheme: " // trim(this%scheme)
         call error_handler%report_error(ERROR_INVALID_CONFIG, error_msg, rc)
         return
      end if

   end subroutine validate_so4chem_config

   !> Print configuration summary
   subroutine print_so4chem_config_summary(this)
      class(SO4chemConfig), intent(in) :: this

      write(*, '(A)') "=== SO4chem Process Configuration ==="
      write(*, '(A,A)') "  Active scheme: ", trim(this%scheme)
      write(*, '(A,I0)') "  Number of species: ", this%n_species
      write(*, '(A,F0.1,A)') "  Minimum time step: ", this%dt_min, " s"
      write(*, '(A,F0.1,A)') "  Maximum time step: ", this%dt_max, " s"
      write(*, '(A,L1)') "  Output diagnostics: ", this%output_diagnostics
      write(*, '(A)') "============================================="

   end subroutine print_so4chem_config_summary

   !> Finalize so4chem configuration
   subroutine finalize_so4chem_config(this)
      class(SO4chemConfig), intent(inout) :: this

      ! Deallocate species names array
      if (allocated(this%species_names)) then
         deallocate(this%species_names)
      end if

      ! Deallocate species indices array
      if (allocated(this%species_indices)) then
         deallocate(this%species_indices)
      end if

      ! Deallocate species properties arrays
      if (allocated(this%species_mw_g)) then
         deallocate(this%species_mw_g)
      end if


      ! Deallocate diagnostic species array
      if (allocated(this%diagnostic_species)) then
         deallocate(this%diagnostic_species)
      end if

      ! Deallocate diagnostic species indices array
      if (allocated(this%diagnostic_species_id)) then
         deallocate(this%diagnostic_species_id)
      end if

   end subroutine finalize_so4chem_config

   !> Validate gocart scheme configuration
   subroutine validate_gocart_config(this, error_handler)
      class(SO4chemSchemeGOCARTConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_gocart_config

   !> Finalize gocart scheme configuration
   subroutine finalize_gocart_config(this)
      class(SO4chemSchemeGOCARTConfig), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_gocart_config



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
   subroutine so4chem_process_load_config(this, config_manager, error_handler)
      class(SO4chemProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      character(len=256) :: scheme_name
      integer :: ierr, rc

      ! Process reads directly from master YAML structure: processes.so4chem
      ! ConfigManager provides generic YAML access, process handles its own configuration

      ! Load process metadata
      call config_manager%get_string("processes/so4chem/name", this%process_name, rc, "so4chem")
      if (rc /= CC_SUCCESS) this%process_name = "so4chem"  ! default

      call config_manager%get_string("processes/so4chem/version", this%process_version, rc, "1.0.0")
      if (rc /= CC_SUCCESS) this%process_version = "1.0.0"  ! default

      call config_manager%get_logical("processes/so4chem/activate", this%is_active, rc, .true.)
      if (rc /= CC_SUCCESS) this%is_active = .true.  ! default

      ! Load process-specific configuration directly from master YAML
      call config_manager%get_string("processes/so4chem/scheme", this%so4chem_config%scheme, rc, "gocart")
      if (rc /= CC_SUCCESS) then
         call error_handler%report_error(ERROR_INVALID_CONFIG, &
            "Missing required 'scheme' in processes/so4chem configuration", rc)
         return
      end if

      ! Load diagnostic switch
      call config_manager%get_logical("processes/so4chem/diagnostics", this%so4chem_config%diagnostics, rc, .false.)
      if (rc /= CC_SUCCESS) this%so4chem_config%diagnostics = .false.  ! Default

      ! Load diagnostic species list
      call config_manager%get_array("processes/so4chem/diag_species", this%so4chem_config%diagnostic_species, &
         rc, default_values=["All"])
      if (rc /= CC_SUCCESS) then
         ! Default to all species if not specified
         allocate(this%so4chem_config%diagnostic_species(1))
         this%so4chem_config%diagnostic_species(1) = "All"
         this%so4chem_config%n_diagnostic_species = 1
      else
         ! Set the count based on the returned array size
         if (allocated(this%so4chem_config%diagnostic_species)) then
            this%so4chem_config%n_diagnostic_species = size(this%so4chem_config%diagnostic_species)
         else
            this%so4chem_config%n_diagnostic_species = 0
         end if
      end if

      ! Species configuration is loaded from ChemState in load_species_from_chem_state
      ! The species come from the master species YAML file (CATChem_species.yml)
      ! and are filtered by is_so4chem property


      ! Load scheme-specific configuration from master YAML
      scheme_name = trim(this%so4chem_config%scheme)
      select case (scheme_name)
       case ('gocart')
         call this%load_gocart_config(config_manager, error_handler)
       case default
         call error_handler%report_error(ERROR_INVALID_STATE, &
            "Unknown so4chem scheme: " // trim(scheme_name), rc)
         return
      end select

   end subroutine so4chem_process_load_config

   !> Load species from explicit list (from species_filter)
   !! This function is used when explicit species are provided via species_filter
   subroutine load_species_from_list(this, chem_state, error_handler, rc)
      use ChemState_Mod, only: ChemStateType

      class(SO4chemProcessConfig), intent(inout) :: this
      type(ChemStateType), pointer, intent(in) :: chem_state
      type(ErrorManagerType), intent(inout) :: error_handler
      integer, intent(inout) :: rc

      integer :: i, species_idx
      character(len=32) :: species_name
      logical :: found

      ! Explicit species list from configuration
      character(len=32), parameter :: EXPLICIT_SPECIES(8) = [ &
         'h2o2                            ', &
         'oh                              ', &
         'no3                             ', &
         'dms                             ', &
         'so2                             ', &
         'so4                             ', &
         'msa                             ', &
         'dms_in                          ' ]

      if (.not. associated(chem_state)) then
         call error_handler%report_error(ERROR_INVALID_STATE, &
            "ChemState not associated in load_species_from_list", rc)
         return
      end if

      ! Set number of species from explicit list
      this%so4chem_config%n_species = 8

      ! Deallocate existing arrays if allocated
      if (allocated(this%so4chem_config%species_names)) then
         deallocate(this%so4chem_config%species_names)
      end if
      if (allocated(this%so4chem_config%species_indices)) then
         deallocate(this%so4chem_config%species_indices)
      end if
      if (allocated(this%so4chem_config%species_mw_g)) then
         deallocate(this%so4chem_config%species_mw_g)
      end if

      ! Allocate arrays
      allocate(this%so4chem_config%species_names(this%so4chem_config%n_species))
      allocate(this%so4chem_config%species_indices(this%so4chem_config%n_species))

      ! Allocate species properties arrays
      allocate(this%so4chem_config%species_mw_g(this%so4chem_config%n_species))

      ! Find indices for each explicit species
      do i = 1, this%so4chem_config%n_species
         species_name = trim(EXPLICIT_SPECIES(i))
         found = .false.

         ! Search for species in ChemState
         do species_idx = 1, size(chem_state%SpeciesNames)
            if (trim(chem_state%SpeciesNames(species_idx)) == species_name) then
               this%so4chem_config%species_names(i) = species_name
               this%so4chem_config%species_indices(i) = species_idx
               found = .true.
               exit
            end if
         end do

         if (.not. found) then
            call error_handler%report_error(ERROR_INVALID_STATE, &
               "Species not found in ChemState: " // trim(species_name), rc)
            return
         end if
      end do

      ! Load species properties from ChemState using found indices
      do i = 1, this%so4chem_config%n_species
         species_idx = this%so4chem_config%species_indices(i)
         this%so4chem_config%species_mw_g(i) = chem_state%ChemSpecies(species_idx)%mw_g
      end do

   end subroutine load_species_from_list



   !> Load gocart scheme configuration from master YAML
   subroutine load_gocart_config(this, config_manager, error_handler)
      class(SO4chemProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: ierr, rc

      ! Load scheme parameters directly from processes/so4chem/gocart/ in master YAML
      call config_manager%get_logical("processes/so4chem/gocart/update_so2", &
         this%gocart_config%update_so2, rc, .true.)
      if (rc /= CC_SUCCESS) this%gocart_config%update_so2 = .true.


   end subroutine load_gocart_config


   !> Validate unified process configuration
   subroutine so4chem_process_validate(this, state_manager, error_handler)
      class(SO4chemProcessConfig), intent(inout) :: this
      type(StateManagerType), intent(in) :: state_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      ! Validate main config
      call this%so4chem_config%validate(error_handler)

      ! Validate scheme-specific config
      select case (trim(this%so4chem_config%scheme))
       case ('gocart')
         call this%gocart_config%validate(error_handler)
      end select

   end subroutine so4chem_process_validate

   !> Finalize unified process configuration
   subroutine so4chem_process_finalize(this)
      class(SO4chemProcessConfig), intent(inout) :: this
      integer :: i

      ! Deallocate persistent state arrays
      if (allocated(this%gocart_persistent_state)) then
         ! Deallocate allocatable components in each column
         do i = 1, size(this%gocart_persistent_state)
            if (allocated(this%gocart_persistent_state(i)%xh2o2_init)) then
               deallocate(this%gocart_persistent_state(i)%xh2o2_init)
            end if
         end do
         deallocate(this%gocart_persistent_state)
      end if

      call this%so4chem_config%finalize()
      call this%gocart_config%finalize()

   end subroutine so4chem_process_finalize

   !> Initialize persistent state arrays for column processing
   subroutine so4chem_initialize_persistent_state(this, grid_manager)
      use GridManager_Mod, only: GridManagerType
      class(SO4chemProcessConfig), intent(inout) :: this
      type(GridManagerType), intent(in) :: grid_manager

      integer :: i

      ! Get total number of columns from GridManager
      this%total_columns = grid_manager%get_total_columns()

      ! Allocate persistent state arrays for each scheme that needs them
      if (.not. allocated(this%gocart_persistent_state)) then
         allocate(this%gocart_persistent_state(this%total_columns))

         ! Initialize each column's state with default values
         do i = 1, this%total_columns
            this%gocart_persistent_state(i)%firsttime = .true.
            this%gocart_persistent_state(i)%nymd_last = -1
            this%gocart_persistent_state(i)%nhms_last_recycle = -1
            ! xh2o2_init is allocatable - will be allocated when needed
         end do
      end if

   end subroutine so4chem_initialize_persistent_state

   !> Get active scheme configuration (polymorphic return)
   function get_active_scheme_config(this) result(scheme_config)
      class(SO4chemProcessConfig), intent(in) :: this
      class(*), allocatable :: scheme_config

      select case (trim(this%so4chem_config%scheme))
       case ('gocart')
         allocate(scheme_config, source=this%gocart_config)
       case default
         ! Return null
      end select

   end function get_active_scheme_config

   !> Map diagnostic species names to indices in the species_names array
   !! This function creates the diagnostic_species_id array that maps each diagnostic species
   !! to its corresponding index in the full species_names array
   subroutine map_diagnostic_species_indices(this, error_handler)
      class(SO4chemProcessConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: i, j, rc
      character(len=256) :: error_msg
      logical :: found_species

      ! Only proceed if diagnostic species are defined
      if (this%so4chem_config%n_diagnostic_species == 0) return

      ! Handle "All" case - map all available species
      if (this%so4chem_config%n_diagnostic_species == 1 .and. &
         trim(this%so4chem_config%diagnostic_species(1)) == "All") then

         ! Deallocate and reallocate for all species
         if (allocated(this%so4chem_config%diagnostic_species_id)) deallocate(this%so4chem_config%diagnostic_species_id)
         allocate(this%so4chem_config%diagnostic_species_id(this%so4chem_config%n_species))
         if (allocated(this%so4chem_config%diagnostic_species)) deallocate(this%so4chem_config%diagnostic_species)
         allocate(this%so4chem_config%diagnostic_species(this%so4chem_config%n_species))
         this%so4chem_config%n_diagnostic_species = this%so4chem_config%n_species
         this%so4chem_config%diagnostic_species = this%so4chem_config%species_names

         ! Map all species indices (1:n_species)
         do i = 1, this%so4chem_config%n_species
            this%so4chem_config%diagnostic_species_id(i) = i
         end do

         return
      end if

      ! Allocate diagnostic species indices array
      if (allocated(this%so4chem_config%diagnostic_species_id)) deallocate(this%so4chem_config%diagnostic_species_id)
      allocate(this%so4chem_config%diagnostic_species_id(this%so4chem_config%n_diagnostic_species))

      ! Map each diagnostic species name to its index in species_names
      do i = 1, this%so4chem_config%n_diagnostic_species
         found_species = .false.

         do j = 1, this%so4chem_config%n_species
            if (trim(this%so4chem_config%diagnostic_species(i)) == trim(this%so4chem_config%species_names(j))) then
               this%so4chem_config%diagnostic_species_id(i) = j
               found_species = .true.
               exit
            end if
         end do

         if (.not. found_species) then
            write(error_msg, '(A,A,A)') "Diagnostic species '", &
               trim(this%so4chem_config%diagnostic_species(i)), &
               "' not found in process species list"
            call error_handler%report_error(ERROR_NOT_FOUND, error_msg, rc)
            !return !do not return and the diagnostics for this unspecified species will be zero in the output
         end if
      end do

   end subroutine map_diagnostic_species_indices

end module SO4chemCommon_Mod
