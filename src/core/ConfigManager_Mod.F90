!> \file ConfigManager_Mod.F90
!! \brief Enhanced configuration management for CATChem
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module provides a modern, flexible configuration management system
!! for CATChem, supporting modular configuration loading, validation, and
!! runtime configuration updates.
!!
!! \details
!! The ConfigManager consolidates and modernizes CATChem configuration management
!! by replacing two previous modules:
!!
!! **Replaces config_opt_mod.F90:**
!! - Legacy config support removed, now using modern YAML-based ConfigDataType
!! - All configuration parameters and options
!! - Process-specific configuration flags
!!
!! **Replaces config_mod.F90:**
!! - Read_Input_File subroutine -> ConfigManager%load_from_file() method
!! - YAML parsing logic with enhanced error handling
!! - Configuration validation and initialization
!!
!! **New features:**
!! - Hierarchical configuration loading with includes and inheritance
!! - Schema-based validation with detailed error reporting
!! - Plugin-based process configuration loading
!! - Runtime configuration updates and hot-reloading
!! - Configuration templates and presets
!! - Environment variable and command-line override support
!!
!! \section config_usage Usage Example
!! \code{.f90}
!! use ConfigManager_Mod
!! type(ConfigManagerType) :: config_mgr
!! type(StateContainer) :: container
!! integer :: rc
!!
!! call config_mgr%init(rc)
!! call config_mgr%load_from_file('CATChem_config.yml', rc)
!! call config_mgr%apply_to_container(container, rc)
!! \endcode
!!
module ConfigManager_Mod
   use iso_c_binding, only: c_associated
   use Precision_Mod, only: fp, MISSING_BOOL, MISSING
   use Error_Mod, only : CC_SUCCESS, CC_FAILURE, ERROR_INVALID_CONFIG, ERROR_INVALID_INPUT, ErrorManagerType
   use Species_Mod, only: SpeciesType
   use yaml_interface_mod, only : yaml_node_t, yaml_load_file, yaml_load_string, yaml_destroy_node, &
      yaml_get_string, yaml_get_integer, yaml_get_real, yaml_get_logical, &
      yaml_has_key, yaml_get, yaml_set, yaml_is_map, yaml_is_sequence, &
      yaml_get_size, yaml_get_string_array, yaml_get_all_keys, &
      yaml_get_real_array, safe_yaml_get_real, safe_yaml_get_logical, &
      safe_yaml_get_integer

   implicit none
   private

   ! Define precision types

   public :: ConfigManagerType
   public :: ConfigDataType      ! Modern YAML-based configuration data structure
   public :: ConfigSchemaType
   public :: ConfigPresetType
   public :: RunPhaseType        ! Run phase configuration type
   public :: ProcessConfigType   ! Process configuration type
   public :: EmissionCategoryMapping  ! Emission category mapping structure
   public :: EmisSpeciesMappingEntry  ! Individual emission species mapping
   public :: EmissionMappingConfig
   public :: discover_yaml_section_items
   public :: discover_nested_yaml_section_items
   public :: CONFIG_STRATEGY_STRICT, CONFIG_STRATEGY_PERMISSIVE, CONFIG_STRATEGY_FALLBACK

   !> \brief Configuration loading strategies
   integer, parameter :: CONFIG_STRATEGY_STRICT = 1    !< Fail on any validation error
   integer, parameter :: CONFIG_STRATEGY_PERMISSIVE = 2 !< Warning on validation errors
   integer, parameter :: CONFIG_STRATEGY_FALLBACK = 3   !< Use defaults on errors

   !> \brief Runtime and MPI configuration
   type :: RuntimeConfig
      integer :: numCPUs = 1                         !< Number of MPI processes
      integer :: thisCPU = 0                         !< Local MPI process handle
      integer :: MPIComm = -1                        !< MPI Communicator Handle
      logical :: isMPI = .false.                     !< Is this an MPI simulation?
      logical :: amIRoot = .true.                    !< Is this the root cpu?
      logical :: DryRun = .false.                    !< Is this a dry run?
      character(len=255) :: SimulationName = ''      !< Name of the simulation
      logical :: DiagEnabled = .false.               !< Was diagnostic output requested?
      logical :: VerboseRequested = .false.          !< Was verbose output requested?
      character(len=10) :: VerboseOnCores = 'root'   !< Which cores should produce verbose output
      logical :: Verbose = .false.                   !< Should verbose output be produced?

      ! Simulation dimensions
      integer :: nLevs = 127                         !< Number of vertical levels
      integer :: nx = 1                              !< Number of grid points in x-direction
      integer :: ny = 1                              !< Number of grid points in y-direction
      integer :: nSpecies = 50                       !< Total number of chemical species
      integer :: maxSpecies = 500                    !< Maximum number of chemical species
      integer :: nSpecies_drydep = 20                !< Number of species with dry deposition
      integer :: nEmissionCategories = 10            !< Number of emission categories
      integer :: nEmissionSpecies = 50               !< Number of emission species per category
   end type RuntimeConfig

   !> \brief File paths and data sources
   type :: FilePathConfig
      character(len=255) :: Emission_File = ''       !< Path to emission mapping file
      character(len=255) :: Species_File = ''        !< Path to species configuration file
      character(len=255) :: Mie_Directory = ''        !< Path to Mie optics data directory
      character(len=255) :: Input_Directory = './'   !< Input data directory
      character(len=255) :: Output_Directory = './'  !< Output data directory
   end type FilePathConfig


   !> \brief External emissions configuration
   type :: ExternalEmisConfig
      logical :: activate = .false.                   !< Enable external emissions
      character(len=256) :: config_file = ''          !< External emissions configuration file
      character(len=64) :: temporal_profile = 'constant' !< Temporal profile type
      logical :: dynamic_mapping = .true.             !< Enable dynamic species mapping
      real(fp) :: global_scale_factor = 1.0_fp        !< Global scaling factor
   end type ExternalEmisConfig

   !> \brief Emission species mapping entry
   !!
   !! Defines how an emission field maps to one or more chemical species
   !! with associated scaling factors
   type :: EmisSpeciesMappingEntry
      character(len=64) :: emission_field = ''        !< Emission field name (e.g., "SeaS1")
      character(len=256) :: long_name = ''            !< Human-readable description
      character(len=64) :: units = ''                  !< units
      integer :: n_mappings = 0                       !< Number of mappings
      character(len=64), allocatable :: map(:)        !< Target chemical species names
      real(fp), allocatable :: scale(:)               !< Scaling factors for each mapping
      integer, allocatable :: index(:)                !< Chemical species indices in ChemState
      logical :: is_active = .true.                   !< Mapping enabled/disabled
   contains
      procedure :: init => emis_species_mapping_init
      procedure :: cleanup => emis_species_mapping_cleanup
      procedure :: copy => emis_species_mapping_copy
   end type EmisSpeciesMappingEntry

   !> \brief Emission category/process mapping
   !!
   !! Contains all emission-to-species mappings for a specific category or process
   type :: EmissionCategoryMapping
      character(len=64) :: category_name = ''         !< Category/process name (e.g., "seasalt")
      integer :: n_emission_species = 0               !< Number of emission species
      type(EmisSpeciesMappingEntry), allocatable :: species_mappings(:) !< Per-species mappings
      logical :: is_active = .true.                   !< Category enabled/disabled
   contains
      procedure :: init => emis_category_mapping_init
      procedure :: cleanup => emis_category_mapping_cleanup
      procedure :: copy => emis_category_mapping_copy
   end type EmissionCategoryMapping

   !> \brief Overall emission mapping configuration
   !!
   !! Contains all emission-to-species mappings for all categories and processes
   type :: EmissionMappingConfig
      integer :: n_categories = 0                     !< Number of emission categories
      type(EmissionCategoryMapping), allocatable :: categories(:) !< All categories
      character(len=256) :: config_file = ''          !< Source configuration file
      logical :: is_loaded = .false.                  !< Whether mapping has been loaded
   contains
      procedure :: init => emis_mapping_config_init
      procedure :: cleanup => emis_mapping_config_cleanup
   end type EmissionMappingConfig

   !> \brief Modernized configuration data structure
   !!
   !! This type provides a modern YAML-based configuration system
   !! with better organization, validation, and extensibility.
   !!
   !! \details
   !! Key features:
   !! - Better organization with nested types for different categories
   !! - Built-in validation and consistency checking
   !! - Support for dynamic process configuration
   !! - Enhanced debugging and introspection capabilities
   !! - Full YAML configuration support
   type :: ConfigDataType

      ! Configuration categories
      type(RuntimeConfig) :: runtime                    !< Runtime and MPI settings
      type(FilePathConfig) :: file_paths                !< File paths and data sources
      type(ExternalEmisConfig) :: external_emissions    !< External emissions configuration
      type(EmissionMappingConfig) :: emission_mapping   !< Emission-to-species mapping configuration

      ! Metadata
      character(len=64) :: config_version = '2.0'       !< Configuration version
      character(len=256) :: source_file = ''            !< Source configuration file
      logical :: is_validated = .false.                 !< Has configuration been validated?
      logical :: run_phases_enabled = .false.           !< Are run phases configured?

      type(ProcessConfigType), allocatable :: run_phase_processes(:)  !< All process configs for run phases
      type(RunPhaseType), allocatable :: run_phases(:)               !< Array of run phase configurations

   contains
      ! Initialization and cleanup
      procedure :: init => config_data_init
      procedure :: cleanup => config_data_cleanup
      procedure :: validate => config_data_validate

      ! Utility methods
      procedure :: print_summary => config_data_print_summary
      procedure :: to_yaml_string => config_data_to_yaml_string
      procedure :: copy => config_data_copy

   end type ConfigDataType

   !> \brief Process configuration for dynamic run phases
   !! This type allows for flexible, phase-based process configuration.
   type :: ProcessConfigType
      character(len=64) :: name             !< Process name
      character(len=64) :: process_type     !< Process type (e.g., 'verticaltransport')
      character(len=64) :: scheme           !< Scheme name (e.g., 'YSU')
      logical :: enabled                    !< Whether process is enabled
      integer :: priority                   !< Execution priority (lower = earlier)
      integer :: process_index              !< Index mapping to global process numbers
      character(len=16) :: timing           !< Timing type ('explicit', 'implicit')
      integer :: subcycling                 !< Number of subcycles
      !character(len=256) :: config_details  !< Additional configuration as string (YAML/JSON)
   end type ProcessConfigType

   !> \brief Run phase configuration for phase-based execution
   !! Contains the phase name, description, frequency, and the list of processes in the phase.
   type :: RunPhaseType
      character(len=64) :: name             !< Phase name
      character(len=256) :: description     !< Phase description
      character(len=32) :: frequency        !< Run frequency (e.g., 'every timestep')
      integer :: subcycling                 !< Number of subcycles for the phase
      integer :: num_processes              !< Number of processes in this phase
      type(ProcessConfigType), allocatable :: processes(:)  !< Array of process configurations
   end type RunPhaseType

   !> \brief Configuration schema definition
   !!
   !! Defines the expected structure and validation rules for configuration files.
   type :: ConfigSchemaType
      private
      character(len=256) :: name                        !< Schema name
      character(len=512) :: description                 !< Schema description
      character(len=64), allocatable :: required_fields(:) !< Required configuration fields
      character(len=64), allocatable :: optional_fields(:) !< Optional configuration fields
      logical :: strict_validation = .true.            !< Strict validation mode
   contains
      procedure :: init => schema_init
      procedure :: add_required_field => schema_add_required_field
      procedure :: add_optional_field => schema_add_optional_field
      procedure :: validate_config => schema_validate_config
      ! Add emission-specific validation methods
      procedure :: validate_emission_config => schema_validate_emission_config
      procedure :: validate_species_mapping => schema_validate_species_mapping
      procedure :: validate_scaling_factors => schema_validate_scaling_factors
   end type ConfigSchemaType

   !> \brief Configuration preset for common use cases
   !!
   !! Provides predefined configuration templates for common modeling scenarios.
   type :: ConfigPresetType
      character(len=256) :: name                        !< Preset name
      character(len=512) :: description                 !< Preset description
      character(len=1024) :: yaml_content               !< YAML configuration content
   contains
      procedure :: load_from_string => preset_load_from_string
      procedure :: save_to_file => preset_save_to_file
   end type ConfigPresetType

   !> \brief Main configuration manager
   !!
   !! Central configuration management with validation, loading, and state application.
   type :: ConfigManagerType
      private

      ! Configuration data
      type(yaml_node_t) :: yaml_data                   !< Loaded YAML configuration
      type(ConfigDataType), public :: config_data              !< Structured configuration data
      logical :: is_loaded = .false.                    !< Configuration loaded flag

      ! Configuration metadata
      character(len=512) :: config_file = ''            !< Current config file path
      character(len=256) :: config_version = ''         !< Configuration version
      integer :: load_strategy = CONFIG_STRATEGY_STRICT !< Loading strategy

      ! Schema and validation
      type(ConfigSchemaType) :: schema                  !< Configuration schema
      logical :: schema_loaded = .false.                !< Schema loaded flag

      ! Environment and overrides
      character(len=64), allocatable :: env_overrides(:)  !< Environment variable overrides
      character(len=512), allocatable :: cli_overrides(:) !< Command line overrides

   contains
      ! Initialization and cleanup
      procedure :: init => config_manager_init
      procedure :: finalize => config_manager_finalize

      ! Configuration loading
      procedure :: load_from_file => config_manager_load_from_file
      procedure :: load_from_string => config_manager_load_from_string
      procedure :: reload => config_manager_reload
      procedure :: load_preset => config_manager_load_preset

      ! Schema management
      procedure :: load_schema => config_manager_load_schema
      procedure :: validate => config_manager_validate

      ! Configuration access
      procedure :: get_string => config_manager_get_string
      procedure :: get_integer => config_manager_get_integer
      procedure :: get_real => config_manager_get_real
      procedure :: get_logical => config_manager_get_logical
      procedure :: get_array => config_manager_get_array
      procedure :: get_real_array => config_manager_get_real_array
      procedure :: get_nspecies => config_manager_get_nspecies
      procedure :: get_max_species => config_manager_get_max_species
      procedure :: get_nemission_categories => config_manager_get_nemission_categories
      procedure :: get_nemission_species => config_manager_get_nemission_species
      procedure :: get_species_file => config_manager_get_species_file
      procedure :: get_emission_file => config_manager_get_emission_file
      procedure :: get_mie_data => config_manager_get_mie_data

      ! State integration
      !procedure :: apply_to_container => config_manager_apply_to_container
      !procedure :: extract_from_container => config_manager_extract_from_container

      ! Species and emission configuration
      procedure :: load_and_init_species => config_manager_load_and_init_species
      procedure :: load_emission_config => config_manager_load_emission_config
      procedure :: load_emission_mapping => config_manager_load_emission_mapping

      ! Emission mapping methods
      procedure :: find_category_mapping => config_manager_find_category_mapping
      procedure :: apply_emission_mapping => config_manager_apply_emission_mapping
      procedure :: get_emission_mapping_for_category => config_manager_get_emission_mapping_for_category

      ! Configuration update methods
      procedure :: update_runtime_from_configs => config_manager_update_runtime_from_configs
      procedure :: parse_config_data
      procedure :: load_run_phases => config_manager_load_run_phases

      ! Utility methods
      procedure :: print_summary => config_manager_print_summary
      procedure :: save_to_file => config_manager_save_to_file
      procedure :: set_loading_strategy => config_manager_set_loading_strategy
      procedure :: add_env_override => config_manager_add_env_override
      procedure :: add_cli_override => config_manager_add_cli_override

   end type ConfigManagerType

   ! Built-in configuration presets
   ! type(ConfigPresetType), parameter :: PRESET_BASIC = ConfigPresetType( &
   !    name = 'basic', &
   !    description = 'Basic CATChem configuration for testing', &
   !    yaml_content = 'simulation: {start_date: "2023-01-01", end_date: "2023-01-02"}' &
   ! )

contains

   !========================================================================
   ! ConfigSchema Implementation
   !========================================================================

   !> \brief Initialize configuration schema
   subroutine schema_init(this, name, description, strict)
      class(ConfigSchemaType), intent(inout) :: this
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: description
      logical, intent(in), optional :: strict

      this%name = trim(name)
      this%description = trim(description)

      if (present(strict)) then
         this%strict_validation = strict
      else
         this%strict_validation = .true.
      endif

      ! Allocate empty arrays
      if (.not. allocated(this%required_fields)) allocate(this%required_fields(0))
      if (.not. allocated(this%optional_fields)) allocate(this%optional_fields(0))

   end subroutine schema_init

   !> \brief Add required field to schema
   subroutine schema_add_required_field(this, field_name)
      class(ConfigSchemaType), intent(inout) :: this
      character(len=*), intent(in) :: field_name

      character(len=64), allocatable :: temp_array(:)
      integer :: n, i

      n = size(this%required_fields)
      allocate(temp_array(n+1))

      do i = 1, n
         temp_array(i) = this%required_fields(i)
      end do
      temp_array(n+1) = trim(field_name)

      call move_alloc(temp_array, this%required_fields)

   end subroutine schema_add_required_field

   !> \brief Add optional field to schema
   subroutine schema_add_optional_field(this, field_name)
      class(ConfigSchemaType), intent(inout) :: this
      character(len=*), intent(in) :: field_name

      character(len=64), allocatable :: temp_array(:)
      integer :: n, i

      n = size(this%optional_fields)
      allocate(temp_array(n+1))

      do i = 1, n
         temp_array(i) = this%optional_fields(i)
      end do
      temp_array(n+1) = trim(field_name)

      call move_alloc(temp_array, this%optional_fields)

   end subroutine schema_add_optional_field

   !> \brief Validate configuration against schema
   subroutine schema_validate_config(this, yaml_data, rc)
      class(ConfigSchemaType), intent(in) :: this
      type(yaml_node_t), intent(in) :: yaml_data
      integer, intent(out) :: rc

      integer :: i
      logical :: key_exists
      character(len=512) :: missing_fields

      rc = CC_SUCCESS
      missing_fields = ''

      ! Basic validation - check if we have a valid node
      if (.not. c_associated(yaml_data%ptr)) then
         rc = CC_FAILURE
         return
      endif

      ! Check required fields
      if (allocated(this%required_fields)) then
         do i = 1, size(this%required_fields)
            key_exists = yaml_has_key(yaml_data, trim(this%required_fields(i)))
            if (.not. key_exists) then
               if (len_trim(missing_fields) > 0) then
                  missing_fields = trim(missing_fields) // ', ' // trim(this%required_fields(i))
               else
                  missing_fields = trim(this%required_fields(i))
               endif
            endif
         end do
      endif

      ! Report missing required fields
      if (len_trim(missing_fields) > 0) then
         write(*, '(A,A)') 'ERROR: Missing required configuration fields: ', trim(missing_fields)
         if (this%strict_validation) then
            rc = CC_FAILURE
            return
         endif
      endif

      ! Validate optional fields exist (if present)
      if (allocated(this%optional_fields)) then
         do i = 1, size(this%optional_fields)
            key_exists = yaml_has_key(yaml_data, trim(this%optional_fields(i)))
            if (key_exists) then
               write(*, '(A,A)') 'INFO: Found optional field: ', trim(this%optional_fields(i))
            endif
         end do
      endif

      write(*, '(A)') 'INFO: Configuration validation completed'

   end subroutine schema_validate_config

   !> \brief Validate emission configuration against schema
   subroutine schema_validate_emission_config(this, yaml_data, config_file, rc)
      implicit none
      class(ConfigSchemaType), intent(in) :: this
      type(yaml_node_t), intent(in) :: yaml_data
      character(len=*), intent(in) :: config_file
      integer, intent(out) :: rc

      logical :: file_exists, key_exists
      type(yaml_node_t) :: emission_config
      integer :: n_sources, n_species, local_rc
      character(len=256) :: data_directory

      ! Suppress warning for unused argument (used for interface compatibility)
      if (.false.) then
         key_exists = yaml_has_key(yaml_data, "dummy")
      endif

      rc = CC_SUCCESS

      ! Check if emission config file exists
      inquire(file=trim(config_file), exist=file_exists)
      if (.not. file_exists) then
         write(*, '(A,A)') 'WARNING: Emission configuration file not found: ', trim(config_file)
         ! Don't fail - emissions may be optional
         return
      endif

      ! Load and validate emission configuration
      emission_config = yaml_load_file(config_file)
      if (.not. c_associated(emission_config%ptr)) then
         write(*, '(A)') 'ERROR: Failed to parse emission configuration file'
         rc = CC_FAILURE
         return
      endif

      ! Validate required emission fields
      key_exists = yaml_has_key(emission_config, "emissions")
      if (.not. key_exists) then
         write(*, '(A)') 'ERROR: Missing emissions section in configuration'
         rc = CC_FAILURE
         call yaml_destroy_node(emission_config)
         return
      endif

      ! Check for data directory
      key_exists = yaml_has_key(emission_config, "emissions/data_directory")
      if (key_exists) then
         if (yaml_get_string(emission_config, "emissions/data_directory", data_directory)) then
            inquire(file=trim(data_directory), exist=file_exists)
            if (.not. file_exists) then
               write(*, '(A,A)') 'WARNING: Emission data directory does not exist: ', trim(data_directory)
            endif
         endif
      endif

      ! Validate emission sources
      call safe_yaml_get_integer(emission_config, "emissions/n_sources", n_sources, local_rc)
      if (local_rc == 0) then
         if (n_sources <= 0) then
            write(*, '(A)') 'WARNING: No emission sources configured'
         else
            write(*, '(A,I0,A)') 'INFO: Found ', n_sources, ' emission sources'
         endif
      endif

      ! Validate species mapping
      call safe_yaml_get_integer(emission_config, "emissions/n_species", n_species, local_rc)
      if (local_rc == 0) then
         if (n_species <= 0) then
            write(*, '(A)') 'WARNING: No emission species configured'
         else
            write(*, '(A,I0,A)') 'INFO: Found ', n_species, ' emission species'
         endif
      endif

      ! Clean up
      call yaml_destroy_node(emission_config)

      write(*, '(A)') 'INFO: Emission configuration validation completed'

   end subroutine schema_validate_emission_config

   !> \brief Validate species mapping consistency
   subroutine schema_validate_species_mapping(this, species_mappings, rc)
      implicit none
      class(ConfigSchemaType), intent(in) :: this
      character(len=*), intent(in) :: species_mappings(:)
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      ! Validate species mapping entries
      do i = 1, size(species_mappings)
         if (len_trim(species_mappings(i)) == 0) then
            write(*, '(A,I0)') 'WARNING: Empty species mapping at index ', i
         endif
      end do

   end subroutine schema_validate_species_mapping

   !> \brief Validate scaling factors
   subroutine schema_validate_scaling_factors(this, scale_factors, rc)
      implicit none
      class(ConfigSchemaType), intent(in) :: this
      real(fp), intent(in) :: scale_factors(:)
      integer, intent(out) :: rc

      integer :: i
      real(fp) :: total_scale

      rc = CC_SUCCESS

      ! Check for reasonable scaling factors
      do i = 1, size(scale_factors)
         if (scale_factors(i) < 0.0_fp) then
            write(*, '(A,I0,A,F8.3)') 'WARNING: Negative scaling factor at index ', i, ': ', scale_factors(i)
         endif
         if (scale_factors(i) > 10.0_fp) then
            write(*, '(A,I0,A,F8.3)') 'WARNING: Large scaling factor at index ', i, ': ', scale_factors(i)
         endif
      end do

      total_scale = sum(scale_factors)
      if (total_scale > 5.0_fp) then
         write(*, '(A,F8.3)') 'WARNING: Very large total scaling factor: ', total_scale
      endif

   end subroutine schema_validate_scaling_factors

   !========================================================================
   ! ConfigPreset Implementation
   !========================================================================

   !> \brief Load preset from string content
   subroutine preset_load_from_string(this, yaml_string)
      implicit none
      class(ConfigPresetType), intent(inout) :: this
      character(len=*), intent(in) :: yaml_string

      this%yaml_content = trim(yaml_string)

   end subroutine preset_load_from_string

   !> \brief Save preset to file
   subroutine preset_save_to_file(this, filename, rc)
      implicit none
      class(ConfigPresetType), intent(in) :: this
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      integer :: unit_num

      rc = CC_SUCCESS

      open(newunit=unit_num, file=trim(filename), status='replace', iostat=rc)
      if (rc /= 0) return

      write(unit_num, '(A)') trim(this%yaml_content)
      close(unit_num)

   end subroutine preset_save_to_file

   !========================================================================
   ! ConfigManager Implementation
   !========================================================================

   !> \brief Initialize configuration manager
   subroutine config_manager_init(this, rc)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Initialize schema with common CATChem fields
      call this%schema%init('CATChem', 'CATChem configuration schema')
      call this%schema%add_required_field('simulation')
      call this%schema%add_optional_field('processes')
      call this%schema%add_optional_field('output')
      call this%schema%add_optional_field('external_emissions')
      call this%schema%add_optional_field('species_mapping')

      this%schema_loaded = .true.

      ! Initialize override arrays
      if (.not. allocated(this%env_overrides)) allocate(this%env_overrides(0))
      if (.not. allocated(this%cli_overrides)) allocate(this%cli_overrides(0))

   end subroutine config_manager_init

   !> \brief Finalize configuration manager
   subroutine config_manager_finalize(this, rc)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: local_rc

      rc = CC_SUCCESS

      ! Clean up configuration data (including emission mapping)
      call this%config_data%cleanup(local_rc)
      if (local_rc /= CC_SUCCESS) then
         rc = local_rc
      endif

      if (this%is_loaded) then
         call yaml_destroy_node(this%yaml_data)
      endif

      this%is_loaded = .false.
      this%schema_loaded = .false.

      if (allocated(this%env_overrides)) deallocate(this%env_overrides)
      if (allocated(this%cli_overrides)) deallocate(this%cli_overrides)

   end subroutine config_manager_finalize

   !> \brief Load configuration from file
   subroutine config_manager_load_from_file(this, filename, rc)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      logical :: file_exists

      rc = CC_SUCCESS

      ! Check if file exists
      inquire(file=trim(filename), exist=file_exists)
      if (.not. file_exists) then
         rc = CC_FAILURE
         return
      endif

      ! Clean up any existing configuration
      if (this%is_loaded) then
         call yaml_destroy_node(this%yaml_data)
      endif

      ! Load YAML configuration using yaml_interface_mod
      this%yaml_data = yaml_load_file(filename)
      if (.not. c_associated(this%yaml_data%ptr)) then
         rc = CC_FAILURE
         return
      endif

      this%config_file = trim(filename)
      this%is_loaded = .true.

      ! Validate against schema if loaded
      if (this%schema_loaded) then
         call this%validate(rc)
         if (rc /= CC_SUCCESS .and. this%load_strategy == CONFIG_STRATEGY_STRICT) then
            return
         endif
      endif

      ! Parse structured configuration data
      call this%parse_config_data(rc)

      ! Parse run phases configuration
      call this%load_run_phases(rc)
      if (rc /= CC_SUCCESS) then
         write(*,*) 'Warning: Failed to load run phases or process configuration'
      endif

   end subroutine config_manager_load_from_file

   !> \brief Load configuration from string
   subroutine config_manager_load_from_string(this, yaml_string, rc)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      character(len=*), intent(in) :: yaml_string
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Clean up any existing configuration
      if (this%is_loaded) then
         call yaml_destroy_node(this%yaml_data)
      endif

      ! Load YAML configuration from string using yaml_interface_mod
      this%yaml_data = yaml_load_string(yaml_string)
      if (.not. c_associated(this%yaml_data%ptr)) then
         rc = CC_FAILURE
         return
      endif

      this%config_file = '<string>'
      this%is_loaded = .true.

      ! Validate against schema if loaded
      if (this%schema_loaded) then
         call this%validate(rc)
         if (rc /= CC_SUCCESS .and. this%load_strategy == CONFIG_STRATEGY_STRICT) then
            return
         endif
      endif

      ! Parse structured configuration data
      call this%parse_config_data(rc)

   end subroutine config_manager_load_from_string

   !> \brief Reload configuration from current file
   subroutine config_manager_reload(this, rc)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      if (len_trim(this%config_file) == 0) then
         rc = CC_FAILURE
         return
      endif

      call this%load_from_file(this%config_file, rc)

   end subroutine config_manager_reload

   !> \brief Load a predefined configuration preset
   subroutine config_manager_load_preset(this, preset, rc)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      type(ConfigPresetType), intent(in) :: preset
      integer, intent(out) :: rc

      call this%load_from_string(preset%yaml_content, rc)

   end subroutine config_manager_load_preset

   !> \brief Load configuration schema
   subroutine config_manager_load_schema(this, schema_file, rc)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      character(len=*), intent(in) :: schema_file
      integer, intent(out) :: rc

      type(yaml_node_t) :: schema_config
      logical :: file_exists, success
      character(len=64) :: required_fields(100), optional_fields(100)
      character(len=256) :: schema_name, schema_description
      integer :: n_required, n_optional, i, local_rc
      character(len=256) :: key

      rc = CC_SUCCESS

      ! Check if schema file exists
      inquire(file=trim(schema_file), exist=file_exists)
      if (.not. file_exists) then
         write(*, '(A,A)') 'INFO: Schema file not found, using default schema: ', trim(schema_file)
         ! Use default schema initialization
         this%schema_loaded = .true.
         return
      endif

      ! Load schema configuration file
      schema_config = yaml_load_file(schema_file)
      if (.not. c_associated(schema_config%ptr)) then
         write(*, '(A)') 'WARNING: Failed to parse schema file, using defaults'
         this%schema_loaded = .true.
         return
      endif

      ! Read schema metadata
      success = yaml_get_string(schema_config, "schema/name", schema_name)
      if (success) then
         write(*, '(A,A)') 'INFO: Loading schema: ', trim(schema_name)
      endif

      success = yaml_get_string(schema_config, "schema/description", &
         schema_description)

      ! Read required fields
      call safe_yaml_get_integer(schema_config, "schema/required_fields/n_fields", n_required, local_rc)
      if (local_rc == 0 .and. n_required > 0) then
         do i = 1, min(n_required, size(required_fields))
            write(key, '(A,I0)') 'schema/required_fields/field_', i
            success = yaml_get_string(schema_config, trim(key), required_fields(i))
            if (success) then
               call this%schema%add_required_field(trim(required_fields(i)))
            endif
         end do
         write(*, '(A,I0,A)') 'INFO: Added ', min(n_required, size(required_fields)), ' required fields to schema'
      endif

      ! Read optional fields
      call safe_yaml_get_integer(schema_config, "schema/optional_fields/n_fields", n_optional, local_rc)
      if (local_rc == 0 .and. n_optional > 0) then
         do i = 1, min(n_optional, size(optional_fields))
            write(key, '(A,I0)') 'schema/optional_fields/field_', i
            success = yaml_get_string(schema_config, trim(key), optional_fields(i))
            if (success) then
               call this%schema%add_optional_field(trim(optional_fields(i)))
            endif
         end do
         write(*, '(A,I0,A)') 'INFO: Added ', min(n_optional, size(optional_fields)), ' optional fields to schema'
      endif

      ! Clean up
      call yaml_destroy_node(schema_config)

      this%schema_loaded = .true.
      write(*, '(A)') 'INFO: Schema loaded successfully'

   end subroutine config_manager_load_schema

   !> \brief Validate loaded configuration
   subroutine config_manager_validate(this, rc)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      if (.not. this%is_loaded) then
         rc = CC_FAILURE
         return
      endif

      if (.not. this%schema_loaded) then
         rc = CC_SUCCESS  ! No schema to validate against
         return
      endif

      call this%schema%validate_config(this%yaml_data, rc)

   end subroutine config_manager_validate

   !> \brief Get string value from configuration
   subroutine config_manager_get_string(this, key, value, rc, default_value)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      character(len=*), intent(in) :: key
      character(len=*), intent(out) :: value
      integer, intent(out) :: rc
      character(len=*), optional, intent(in) :: default_value

      ! Try to get value from loaded YAML data using yaml_get interface
      if (this%is_loaded) then
         call yaml_get(this%yaml_data, key, value, rc, default_value)
      else
         ! Use default value if provided
         if (present(default_value)) then
            value = default_value
            rc = CC_SUCCESS
         else
            value = ''
            rc = CC_FAILURE
         endif
      endif

   end subroutine config_manager_get_string

   !> \brief Get integer value from configuration
   subroutine config_manager_get_integer(this, key, value, rc, default_value)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      character(len=*), intent(in) :: key
      integer, intent(out) :: value
      integer, intent(out) :: rc
      integer, optional, intent(in) :: default_value

      ! Try to get value from loaded YAML data using yaml_get interface
      if (this%is_loaded) then
         call yaml_get(this%yaml_data, key, value, rc, default_value)
      else
         ! Use default value if provided
         if (present(default_value)) then
            value = default_value
            rc = CC_SUCCESS
         else
            value = 0
            rc = CC_FAILURE
         endif
      endif

   end subroutine config_manager_get_integer

   !> \brief Get real value from configuration
   subroutine config_manager_get_real(this, key, value, rc, default_value)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      character(len=*), intent(in) :: key
      real(fp), intent(out) :: value
      integer, intent(out) :: rc
      real(fp), optional, intent(in) :: default_value

      ! Try to get value from loaded YAML data using safe conversion
      if (this%is_loaded) then
         call safe_yaml_get_real(this%yaml_data, key, value, rc)
         if (rc /= 0 .and. present(default_value)) then
            value = default_value
            rc = CC_SUCCESS
         endif
      else
         ! Use default value if provided
         if (present(default_value)) then
            value = default_value
            rc = CC_SUCCESS
         else
            value = 0.0_fp
            rc = CC_FAILURE
         endif
      endif

   end subroutine config_manager_get_real

   !> \brief Get logical value from configuration
   subroutine config_manager_get_logical(this, key, value, rc, default_value)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      character(len=*), intent(in) :: key
      logical, intent(out) :: value
      integer, intent(out) :: rc
      logical, optional, intent(in) :: default_value

      ! Try to get value from loaded YAML data using safe conversion
      if (this%is_loaded) then
         call safe_yaml_get_logical(this%yaml_data, key, value, rc)
         if (rc /= 0 .and. present(default_value)) then
            value = default_value
            rc = CC_SUCCESS
         endif
      else
         ! Use default value if provided
         if (present(default_value)) then
            value = default_value
            rc = CC_SUCCESS
         else
            value = .false.
            rc = CC_FAILURE
         endif
      endif

   end subroutine config_manager_get_logical

   !> \brief Get array value from configuration
   subroutine config_manager_get_array(this, key, values, rc, default_values)
      use yaml_interface_mod, only: yaml_get_string_array
      implicit none
      class(ConfigManagerType), intent(in) :: this
      character(len=*), intent(in) :: key
      character(len=*), allocatable, intent(out) :: values(:)
      integer, intent(out) :: rc
      character(len=*), optional, intent(in) :: default_values(:)

      integer :: actual_size, max_size
      logical :: success
      character(len=256) :: temp_values(100)  ! Temporary array with maximum size

      rc = CC_FAILURE

      ! Check if configuration is loaded
      if (.not. this%is_loaded) then
         ! Use default values if provided
         if (present(default_values)) then
            allocate(values(size(default_values)))
            values = default_values
            rc = CC_SUCCESS
         else
            allocate(values(0))
         endif
         return
      endif

      ! Try to get string array from loaded YAML data
      max_size = size(temp_values)
      success = yaml_get_string_array(this%yaml_data, key, temp_values, actual_size)

      if (success .and. actual_size > 0) then
         ! Allocate the output array with the correct size
         allocate(values(actual_size))
         values(1:actual_size) = temp_values(1:actual_size)
         rc = CC_SUCCESS
      else
         ! Use default values if provided
         if (present(default_values)) then
            allocate(values(size(default_values)))
            values = default_values
            rc = CC_SUCCESS
         else
            ! Return empty array if unsuccessful and no default provided
            allocate(values(0))
            rc = CC_FAILURE
         endif
      endif

   end subroutine config_manager_get_array

   !> \brief Get real array value from configuration
   subroutine config_manager_get_real_array(this, key, values, rc)
      use yaml_interface_mod, only: yaml_get_real_array
      implicit none
      class(ConfigManagerType), intent(in) :: this
      character(len=*), intent(in) :: key
      real(fp), allocatable, intent(out) :: values(:)
      integer, intent(out) :: rc

      integer :: actual_size, max_size
      logical :: success
      real(fp) :: temp_values(100)  ! Temporary array with maximum size

      rc = CC_FAILURE

      ! Check if configuration is loaded
      if (.not. this%is_loaded) then
         allocate(values(0))
         return
      endif

      ! Try to get real array from loaded YAML data
      max_size = size(temp_values)
      success = yaml_get_real_array(this%yaml_data, key, temp_values, actual_size)

      if (success .and. actual_size > 0) then
         ! Allocate the output array with the correct size
         allocate(values(actual_size))
         values(1:actual_size) = temp_values(1:actual_size)
         rc = CC_SUCCESS
      else
         ! Return empty array if unsuccessful
         allocate(values(0))
         rc = CC_FAILURE
      endif

   end subroutine config_manager_get_real_array

   !> \brief Get number of species from configuration
   function config_manager_get_nspecies(this) result(nspecies)
      implicit none
      class(ConfigManagerType), intent(in) :: this
      integer :: nspecies

      nspecies = this%config_data%runtime%nSpecies
   end function config_manager_get_nspecies

   !> \brief Get maximum number of species from configuration
   function config_manager_get_max_species(this) result(max_species)
      implicit none
      class(ConfigManagerType), intent(in) :: this
      integer :: max_species

      max_species = this%config_data%runtime%maxSpecies
   end function config_manager_get_max_species

   !> \brief Get number of emission categories from configuration
   function config_manager_get_nemission_categories(this) result(nemission_categories)
      implicit none
      class(ConfigManagerType), intent(in) :: this
      integer :: nemission_categories

      nemission_categories = this%config_data%runtime%nEmissionCategories
   end function config_manager_get_nemission_categories

   !> \brief Get number of emission species from configuration
   function config_manager_get_nemission_species(this) result(nemission_species)
      implicit none
      class(ConfigManagerType), intent(in) :: this
      integer :: nemission_species

      nemission_species = this%config_data%runtime%nEmissionSpecies
   end function config_manager_get_nemission_species

   !> \brief Get species filename from configuration
   function config_manager_get_species_file(this) result(species_file)
      implicit none
      class(ConfigManagerType), intent(in) :: this
      character(len=255) :: species_file

      species_file = this%config_data%file_paths%Species_File
   end function config_manager_get_species_file

   !> \brief Get emission filename from configuration
   function config_manager_get_emission_file(this) result(emission_file)
      implicit none
      class(ConfigManagerType), intent(in) :: this
      character(len=255) :: emission_file

      emission_file = this%config_data%file_paths%Emission_File
   end function config_manager_get_emission_file

   !> \brief Print configuration summary
   subroutine config_manager_print_summary(this)
      implicit none
      class(ConfigManagerType), intent(in) :: this

      write(*,'(A)') '=== ConfigManager Summary ==='
      write(*,'(A,A)') 'Config file: ', trim(this%config_file)
      write(*,'(A,L1)') 'Loaded: ', this%is_loaded
      write(*,'(A,L1)') 'Schema loaded: ', this%schema_loaded
      write(*,'(A,I0)') 'Environment overrides: ', size(this%env_overrides)
      write(*,'(A,I0)') 'CLI overrides: ', size(this%cli_overrides)
      write(*,'(A)') '=============================='

   end subroutine config_manager_print_summary

   !> \brief Save current configuration to file
   subroutine config_manager_save_to_file(this, filename, rc)
      implicit none
      class(ConfigManagerType), intent(in) :: this
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      integer :: unit_num, io_stat

      rc = CC_SUCCESS

      if (.not. this%is_loaded) then
         rc = CC_FAILURE
         return
      endif

      ! Open file for writing
      open(newunit=unit_num, file=trim(filename), status='replace', &
         action='write', iostat=io_stat)
      if (io_stat /= 0) then
         rc = CC_FAILURE
         return
      endif

      ! Write basic configuration info (placeholder implementation)
      write(unit_num, '(A)') '# CATChem Configuration (Generated)'
      write(unit_num, '(A)') '# This is a placeholder - full YAML serialization not yet implemented'
      write(unit_num, '(A,I0)') 'nSpecies: ', this%config_data%runtime%nSpecies
      write(unit_num, '(A,I0)') 'nEmissionCategories: ', this%config_data%runtime%nEmissionCategories
      write(unit_num, '(A,I0)') 'nEmissionSpecies: ', this%config_data%runtime%nEmissionSpecies

      close(unit_num)

   end subroutine config_manager_save_to_file

   !> \brief Set configuration loading strategy
   subroutine config_manager_set_loading_strategy(this, strategy)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      integer, intent(in) :: strategy

      this%load_strategy = strategy

   end subroutine config_manager_set_loading_strategy

   !> \brief Add environment variable override
   subroutine config_manager_add_env_override(this, env_var, rc)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      character(len=*), intent(in) :: env_var
      integer, intent(out) :: rc

      character(len=64), allocatable :: temp_array(:)
      integer :: n, i

      rc = CC_SUCCESS

      n = size(this%env_overrides)
      allocate(temp_array(n+1))

      do i = 1, n
         temp_array(i) = this%env_overrides(i)
      end do
      temp_array(n+1) = trim(env_var)

      call move_alloc(temp_array, this%env_overrides)

   end subroutine config_manager_add_env_override

   !> \brief Add command line override
   subroutine config_manager_add_cli_override(this, cli_arg, rc)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      character(len=*), intent(in) :: cli_arg
      integer, intent(out) :: rc

      character(len=512), allocatable :: temp_array(:)
      integer :: n, i

      rc = CC_SUCCESS

      n = size(this%cli_overrides)
      allocate(temp_array(n+1))

      do i = 1, n
         temp_array(i) = this%cli_overrides(i)
      end do
      temp_array(n+1) = trim(cli_arg)

      call move_alloc(temp_array, this%cli_overrides)

   end subroutine config_manager_add_cli_override

   !========================================================================
   ! ConfigDataType Procedures
   !========================================================================

   !> \brief Initialize configuration data
   subroutine config_data_init(this, rc)
      implicit none
      class(ConfigDataType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Initialize all components to defaults (already done by default initialization)
      this%is_validated = .false.
      this%source_file = ''
      this%config_version = '2.0'

      ! Initialize emission mapping
      call this%emission_mapping%init()

      ! Initialize emission mapping
      call this%emission_mapping%init()

   end subroutine config_data_init

   !> \brief Clean up configuration data
   subroutine config_data_cleanup(this, rc)
      implicit none
      class(ConfigDataType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Reset to defaults
      call this%init(rc)

   end subroutine config_data_cleanup

   !> \brief Validate configuration data
   function config_data_validate(this, error_mgr, rc) result(is_valid)
      implicit none
      class(ConfigDataType), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_mgr
      integer, intent(out) :: rc
      logical :: is_valid

      rc = CC_SUCCESS
      is_valid = .true.

      call error_mgr%push_context('config_data_validate', 'ConfigManager_Mod.F90')

      ! Basic validation checks
      if (this%runtime%numCPUs < 1) then
         call error_mgr%report_error(ERROR_INVALID_CONFIG, &
            'Number of CPUs must be positive', rc, 'config_data_validate')
         is_valid = .false.
         call error_mgr%pop_context()
         return
      endif

      if (this%runtime%nLevs < 1) then
         call error_mgr%report_error(ERROR_INVALID_CONFIG, &
            'Number of levels must be positive', rc, 'config_data_validate')
         is_valid = .false.
         call error_mgr%pop_context()
         return
      endif

      if (this%runtime%nx < 1) then
         call error_mgr%report_error(ERROR_INVALID_CONFIG, &
            'Number of nx must be positive', rc, 'config_data_validate')
         is_valid = .false.
         call error_mgr%pop_context()
         return
      endif

      if (this%runtime%ny < 1) then
         call error_mgr%report_error(ERROR_INVALID_CONFIG, &
            'Number of ny must be positive', rc, 'config_data_validate')
         is_valid = .false.
         call error_mgr%pop_context()
         return
      endif

      if (this%runtime%maxSpecies < 1) then
         call error_mgr%report_error(ERROR_INVALID_CONFIG, &
            'Maximum species must be positive', rc, 'config_data_validate')
         is_valid = .false.
         call error_mgr%pop_context()
         return
      endif

      if (len_trim(this%file_paths%Input_Directory) == 0) then
         call error_mgr%report_error(ERROR_INVALID_CONFIG, &
            'Input directory must be specified', rc, 'config_data_validate')
         is_valid = .false.
         call error_mgr%pop_context()
         return
      endif

      ! Mark as validated if all checks pass
      this%is_validated = is_valid
      call error_mgr%pop_context()

   end function config_data_validate

   !> \brief Print configuration summary
   subroutine config_data_print_summary(this)
      implicit none
      class(ConfigDataType), intent(in) :: this

      write(*, '(A)') '=== Configuration Summary ==='
      write(*, '(A,A)') 'Version: ', trim(this%config_version)
      write(*, '(A,A)') 'Source file: ', trim(this%source_file)
      write(*, '(A,L1)') 'Validated: ', this%is_validated
      write(*, '(A,I0)') 'Number of CPUs: ', this%runtime%numCPUs
      write(*, '(A,A)') 'Simulation name: ', trim(this%runtime%SimulationName)
      write(*, '(A,L1)') 'External emissions activated: ', this%external_emissions%activate
      write(*, '(A)') '============================'

   end subroutine config_data_print_summary

   !> \brief Convert configuration to YAML string
   function config_data_to_yaml_string(this, rc) result(yaml_string)
      implicit none
      class(ConfigDataType), intent(in) :: this
      integer, intent(out) :: rc
      character(len=:), allocatable :: yaml_string

      rc = CC_SUCCESS

      ! This is a placeholder - would convert to actual YAML format
      yaml_string = 'configuration:\n  version: ' // trim(this%config_version) // '\n'

   end function config_data_to_yaml_string

   !> \brief Copy configuration data
   subroutine config_data_copy(this, source, rc)
      implicit none
      class(ConfigDataType), intent(inout) :: this
      class(ConfigDataType), intent(in) :: source
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Deep copy all components
      this%runtime = source%runtime
      this%file_paths = source%file_paths
      this%external_emissions = source%external_emissions
      this%config_version = source%config_version
      this%source_file = source%source_file
      this%is_validated = source%is_validated

   end subroutine config_data_copy

   !> \brief Parse YAML data into structured configuration
   subroutine parse_config_data(this, rc)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      !local variables
      integer :: local_rc

      rc = CC_SUCCESS

      if (.not. this%is_loaded) then
         rc = CC_FAILURE
         return
      endif

      ! Parse runtime configuration - use safe conversion for numeric values
      call safe_yaml_get_integer(this%yaml_data, 'runtime/nEmissionSpecies', this%config_data%runtime%nEmissionSpecies, local_rc)
      if (local_rc /= 0) this%config_data%runtime%nEmissionSpecies = 50  ! default value

      call safe_yaml_get_logical(this%yaml_data, 'diagnostics/output/enabled', this%config_data%runtime%DiagEnabled, local_rc)
      if (local_rc /= 0) this%config_data%runtime%DiagEnabled = .false.  ! default value

      ! Parse file paths
      call yaml_get(this%yaml_data, 'diagnostics/output/directory', this%config_data%file_paths%Output_Directory, rc, './')
      call yaml_get(this%yaml_data, 'mie/directory', this%config_data%file_paths%Mie_Directory, rc, './')
      call yaml_get(this%yaml_data, 'simulation/species_filename', this%config_data%file_paths%Species_File, rc, '')
      call yaml_get(this%yaml_data, 'simulation/emission_filename', this%config_data%file_paths%Emission_File, rc, '')

      ! Parse external emissions configuration
      if (yaml_has_key(this%yaml_data, 'external_emissions')) then
         call safe_yaml_get_real(this%yaml_data, 'external_emissions/global_scale_factor', this%config_data%external_emissions%global_scale_factor, local_rc)
         if (local_rc /= 0) this%config_data%external_emissions%global_scale_factor = 1.0_fp  ! default value
      endif

      ! Mark configuration as validated
      this%config_data%is_validated = .true.
      this%config_data%source_file = this%config_file

   end subroutine parse_config_data

   !> \brief Load complete species data with all properties from YAML file and initialize ChemState
   !!
   !! This enhanced method loads species configuration including all properties
   !! such as molecular weight, type flags, background concentrations, etc.,
   !! and initializes the provided ChemState with the loaded species data.
   !! Based on the YAML structure in test_species.yml.
   !!
   !! \param[in]    filename     Species configuration file path
   !! \param[inout] chem_state   ChemState to initialize with species data
   !! \param[in]    grid         Grid geometry for concentration array initialization
   !! \param[out]   num_species  Number of species loaded
   !! \param[out]   rc          Return code
   subroutine config_manager_load_and_init_species(this, filename, chem_state, error_mgr, grid, rc, num_species)
      use ChemState_Mod, only: ChemStateType
      use Error_Mod, only: ErrorManagerType
      use GridGeometry_Mod, only: GridGeometryType
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      character(len=*), intent(in) :: filename
      type(ChemStateType), intent(inout) :: chem_state
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      type(GridGeometryType), pointer, intent(in) :: grid
      integer, intent(out) :: rc
      integer, optional, intent(out) :: num_species

      type(yaml_node_t) :: species_config
      logical :: file_exists, success
      integer :: i, j, list_size, total_keys, species_index
      character(len=256) :: species_path
      character(len=64), allocatable :: species_keys(:)
      character(len=64) :: all_yaml_keys(200)

      rc = CC_SUCCESS

      ! Check if file exists
      inquire(file=trim(filename), exist=file_exists)
      if (.not. file_exists) then
         write(*, '(A,A)') 'ERROR: Species configuration file not found: ', trim(filename)
         rc = CC_FAILURE
         return
      endif

      ! Load species configuration file
      species_config = yaml_load_file(filename)
      if (.not. c_associated(species_config%ptr)) then
         write(*, '(A)') 'ERROR: Failed to load species configuration file'
         rc = CC_FAILURE
         return
      endif

      ! Check if this is a map/dictionary structure
      if (.not. yaml_is_map(species_config)) then
         write(*, '(A)') 'ERROR: Species configuration file must be a YAML map/dictionary'
         rc = CC_FAILURE
         call yaml_destroy_node(species_config)
         return
      endif

      ! Get all top-level keys - these are the species names
      total_keys = yaml_get_size(species_config)
      success = yaml_get_all_keys(species_config, all_yaml_keys, list_size)

      if (.not. success .or. list_size <= 0) then
         write(*, '(A)') 'ERROR: Failed to read species keys from configuration file'
         rc = CC_FAILURE
         call yaml_destroy_node(species_config)
         return
      endif

      write(*, '(A,I0,A)') 'INFO: Found ', list_size, ' species keys in configuration'

      ! Initialize ChemState with the number of species and grid geometry
      call chem_state%init(list_size, error_mgr, rc, grid)
      if (rc /= CC_SUCCESS) then
         write(*, '(A)') 'ERROR: Failed to initialize ChemState'
         call yaml_destroy_node(species_config)
         return
      endif

      ! Allocate species keys array
      allocate(species_keys(list_size))

      ! Copy species keys (null character cleaning now handled in yaml_get_all_keys)
      do i = 1, list_size
         species_keys(i) = all_yaml_keys(i)
      end do

      ! Load each species and directly populate ChemState in single loop
      do i = 1, list_size
         species_path = trim(species_keys(i))

         ! Load species properties directly into ChemState%ChemSpecies array
         call load_species_properties(species_config, trim(species_path), chem_state%ChemSpecies(i), rc)
         if (rc /= CC_SUCCESS) then
            write(*, '(A,A,A,I0)') 'ERROR: Failed to load species properties for "', &
               trim(species_path), '", species #', i
            deallocate(species_keys)
            call yaml_destroy_node(species_config)
            return
         endif

         ! Directly populate ChemState arrays - species_index is just i
         species_index = i
         chem_state%nSpecies = species_index
         ! Clean the species name before storing it
         chem_state%ChemSpecies(i)%short_name = trim(adjustl(chem_state%ChemSpecies(i)%short_name))
         chem_state%SpeciesNames(species_index) = trim(chem_state%ChemSpecies(i)%short_name)
         chem_state%SpeciesIndex(species_index) = species_index

         ! Classify species into categories (allow multiple classifications)
         if (chem_state%ChemSpecies(i)%is_gas) then
            chem_state%nSpeciesGas = chem_state%nSpeciesGas + 1
            chem_state%GasIndex(chem_state%nSpeciesGas) = species_index
         endif

         if (chem_state%ChemSpecies(i)%is_aerosol) then
            chem_state%nSpeciesAero = chem_state%nSpeciesAero + 1
            chem_state%AeroIndex(chem_state%nSpeciesAero) = species_index
         endif

         if (chem_state%ChemSpecies(i)%is_dust) then
            chem_state%nSpeciesDust = chem_state%nSpeciesDust + 1
            chem_state%DustIndex(chem_state%nSpeciesDust) = species_index
         endif

         if (chem_state%ChemSpecies(i)%is_seasalt) then
            chem_state%nSpeciesSeaSalt = chem_state%nSpeciesSeaSalt + 1
            chem_state%SeaSaltIndex(chem_state%nSpeciesSeaSalt) = species_index
         endif

         if (chem_state%ChemSpecies(i)%is_drydep) then
            chem_state%nSpeciesDryDep = chem_state%nSpeciesDryDep + 1
            chem_state%DryDepIndex(chem_state%nSpeciesDryDep) = species_index

            ! Also add to aerosol dry deposition if it's an aerosol
            if (chem_state%ChemSpecies(i)%is_aerosol .and. chem_state%ChemSpecies(i)%is_drydep) then
               chem_state%nSpeciesAeroDryDep = chem_state%nSpeciesAeroDryDep + 1
               chem_state%AeroDryDepIndex(chem_state%nSpeciesAeroDryDep) = species_index
            endif
         endif

         if (chem_state%ChemSpecies(i)%is_wetdep) then
            chem_state%nSpeciesWetDep = chem_state%nSpeciesWetDep + 1
            chem_state%WetDepIndex(chem_state%nSpeciesWetDep) = species_index
         endif

         if (chem_state%ChemSpecies(i)%is_tracer) then
            chem_state%nSpeciesTracer = chem_state%nSpeciesTracer + 1
            chem_state%TracerIndex(chem_state%nSpeciesTracer) = species_index
         endif

         !print species info as a test
         write(*, '(A,A)') 'Species name: ', chem_state%ChemSpecies(i)%short_name
         ! write(*, '(A,A)') 'Description: ', chem_state%ChemSpecies(i)%description
         ! write(*, *) 'lower radius: ', chem_state%ChemSpecies(i)%lower_radius
         ! write(*, *) 'upper radius: ', chem_state%ChemSpecies(i)%upper_radius
         ! write(*, *) 'density: ', chem_state%ChemSpecies(i)%density
         ! write(*, *) 'Molecular weight: ', chem_state%ChemSpecies(i)%mw_g
         ! write(*, *) 'is gas: ', chem_state%ChemSpecies(i)%is_gas
         ! write(*, *) 'is aerosol: ', chem_state%ChemSpecies(i)%is_aerosol
         ! write(*, *) 'is dust: ', chem_state%ChemSpecies(i)%is_dust
         ! write(*, *) 'is sea salt: ', chem_state%ChemSpecies(i)%is_seasalt
         ! write(*, *) 'is dry deposition: ', chem_state%ChemSpecies(i)%is_drydep
         ! write(*, *) 'is tracer: ', chem_state%ChemSpecies(i)%is_tracer
         write(*, *) 'wd_rainouteff: ', chem_state%ChemSpecies(i)%wd_rainouteff

      enddo

      ! Clean up
      deallocate(species_keys)
      call yaml_destroy_node(species_config)

      write(*, '(A,I0,A)') 'INFO: Successfully initialized ChemState with ', list_size, ' species'
      write(*, '(A,I0)') '  Gas species: ', chem_state%nSpeciesGas
      write(*, '(A,I0)') '  Aerosol species: ', chem_state%nSpeciesAero
      write(*, '(A,I0)') '  Dust species: ', chem_state%nSpeciesDust
      write(*, '(A,I0)') '  Sea salt species: ', chem_state%nSpeciesSeaSalt
      write(*, '(A,I0)') '  Dry deposition species: ', chem_state%nSpeciesDryDep
      write(*, '(A,I0)') '  Wet deposition species: ', chem_state%nSpeciesWetDep
      write(*, '(A,I0)') '  Aerosol dry deposition species: ', chem_state%nSpeciesAeroDryDep
      write(*, '(A,I0)') '  Tracer species: ', chem_state%nSpeciesTracer

      ! Initialize Mie data if configuration is available
      call this%get_mie_data(chem_state, rc)
      if (rc /= CC_SUCCESS) then
         write(*, '(A)') 'WARNING: Mie data initialization failed or skipped'
         ! Don't fail the overall species loading for Mie issues
         rc = CC_SUCCESS
      end if

      ! Return optional output parameters
      if (present(num_species)) then
         num_species = list_size
      endif

   end subroutine config_manager_load_and_init_species

   !> \brief Helper subroutine to load species properties from YAML
   !! \param[in] yaml_root Root YAML node
   !! \param[in] species_path Path to the species data in YAML
   !! \param[inout] species Species object to populate
   !! \param[out] rc Return code
   subroutine load_species_properties(yaml_root, species_path, species, rc)
      implicit none
      type(yaml_node_t), intent(in) :: yaml_root
      character(len=*), intent(in) :: species_path
      type(SpeciesType), intent(inout) :: species
      integer, intent(out) :: rc

      logical :: success
      character(len=256) :: field_path
      character(len=64) :: species_name
      real(fp) :: temp_real  ! Using project-wide fp precision
      real(fp), allocatable :: temp_real_array(:)  ! Using project-wide fp precision
      logical :: temp_logical
      character(len=256) :: temp_string
      integer :: yaml_rc  ! Separate return code for YAML operations
      integer :: i, j, actual_size  ! Loop variables for debugging

      rc = CC_SUCCESS

      ! Initialize species with defaults
      call species%init('UNKNOWN', 'Unknown Species', 28.0_fp, rc)
      if (rc /= CC_SUCCESS) then
         write(*, '(A,A)') 'ERROR: Failed to initialize species at path: ', trim(species_path)
         return
      endif

      ! Load species name (required) - try 'name' field first, then use path as fallback
      write(field_path, '(A,A)') trim(species_path), '/name'
      success = yaml_get_string(yaml_root, trim(field_path), species_name)
      if (success) then
         species%short_name = trim(adjustl(species_name))
         species%long_name = trim(adjustl(species_name))
      else
         ! Use the key/path as the species name (for flat format like so2:, so4:)
         species_name = species_path
         species%short_name = trim(adjustl(species_name))
         species%long_name = trim(adjustl(species_name))
         ! Suppress verbose output: using key as fallback name is expected behavior
         ! write(*, '(A,A,A,A)') 'INFO: No name field found, using key "', &
         !                      trim(species_path), '" as species name'
      endif

      ! Load long_name if explicitly provided
      write(field_path, '(A,A)') trim(species_path), '/long_name'
      call yaml_get(yaml_root, trim(field_path), temp_string, yaml_rc)
      if (yaml_rc == 0) then
         species%long_name = trim(adjustl(temp_string))
      endif

      ! Load description (optional)
      write(field_path, '(A,A)') trim(species_path), '/description'
      call yaml_get(yaml_root, trim(field_path), temp_string, yaml_rc)
      if (yaml_rc == 0) then
         species%description = trim(adjustl(temp_string))
      endif

      ! Load mie name (optional)
      write(field_path, '(A,A)') trim(species_path), '/mie_name'
      call yaml_get(yaml_root, trim(field_path), temp_string, yaml_rc)
      if (yaml_rc == 0) then
         species%mie_name = trim(adjustl(temp_string))
      endif

      ! Load molecular weight (optional, but important) - use safe conversion for numeric values
      write(field_path, '(A,A)') trim(species_path), '/mw_g'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%mw_g = temp_real
      else
         ! No molecular weight specified - keep default from init
         write(*, '(A,A,A)') 'WARNING: No molecular_weight found for species ', &
            trim(species%short_name), ', using default 28.0'
      endif

      ! Load physical properties
      write(field_path, '(A,A)') trim(species_path), '/density'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%density = temp_real
      else
         species%density = MISSING
      endif

      write(field_path, '(A,A)') trim(species_path), '/radius'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%radius = temp_real
      else
         species%radius = MISSING
      endif

      write(field_path, '(A,A)') trim(species_path), '/lower_radius'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%lower_radius = temp_real
      else
         species%lower_radius = MISSING
      endif

      write(field_path, '(A,A)') trim(species_path), '/upper_radius'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%upper_radius = temp_real
      else
         species%upper_radius = MISSING
      endif

      write(field_path, '(A,A)') trim(species_path), '/viscosity'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%viscosity = temp_real
      else
         species%viscosity = MISSING
      endif

      write(field_path, '(A,A)') trim(species_path), '/dd_f0'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%dd_f0 = temp_real
      else
         species%dd_f0 = MISSING
      endif

      write(field_path, '(A,A)') trim(species_path), '/dd_hstar'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%dd_hstar = temp_real
      else
         species%dd_hstar = MISSING
      endif

      write(field_path, '(A,A)') trim(species_path), '/dd_DvzAerSnow'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%dd_DvzAerSnow = temp_real
      else
         species%dd_DvzAerSnow = MISSING
      endif

      write(field_path, '(A,A)') trim(species_path), '/dd_DvzMinVal_snow'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%dd_DvzMinVal_snow = temp_real
      else
         species%dd_DvzMinVal_snow = MISSING
      endif

      write(field_path, '(A,A)') trim(species_path), '/dd_DvzMinVal_land'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%dd_DvzMinVal_land = temp_real
      else
         species%dd_DvzMinVal_land = MISSING
      endif

      write(field_path, '(A,A)') trim(species_path), '/henry_k0'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%henry_k0 = temp_real
      else
         species%henry_k0 = MISSING
      endif

      write(field_path, '(A,A)') trim(species_path), '/henry_cr'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%henry_cr = temp_real
      else
         species%henry_cr = MISSING
      endif

      write(field_path, '(A,A)') trim(species_path), '/henry_pKa'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%henry_pKa = temp_real
      else
         species%henry_pKa = 0.0_fp  ! Default to 0.0 if not specified
      endif

      write(field_path, '(A,A)') trim(species_path), '/wd_retfactor'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%wd_retfactor = temp_real
      else
         species%wd_retfactor = MISSING
      endif

      write(field_path, '(A,A)') trim(species_path), '/wd_LiqAndGas'
      call safe_yaml_get_logical(yaml_root, trim(field_path), temp_logical, yaml_rc)
      if (yaml_rc == 0) then
         species%wd_LiqAndGas = temp_logical
      else
         species%wd_LiqAndGas = MISSING_BOOL
      endif

      write(field_path, '(A,A)') trim(species_path), '/wd_convfacI2G'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%wd_convfacI2G = temp_real
      else
         species%wd_convfacI2G = MISSING
      endif

      write(field_path, '(A,A)') trim(species_path), '/wd_rainouteff'
      allocate(temp_real_array(10))  ! Assume max size of 10 for temporary array
      success = yaml_get_real_array(yaml_root, trim(field_path), temp_real_array, actual_size)
      if (success .and. actual_size > 0) then
         species%wd_rainouteff(1:actual_size) = temp_real_array(1:actual_size)
         deallocate(temp_real_array)
      else
         ! Return missing array
         species%wd_rainouteff(:) = MISSING
         deallocate(temp_real_array)
      endif

      ! Load type flags (with proper default handling)
      write(field_path, '(A,A)') trim(species_path), '/is_gas'
      call safe_yaml_get_logical(yaml_root, trim(field_path), temp_logical, yaml_rc)
      if (yaml_rc == 0) then
         species%is_gas = temp_logical
      else
         species%is_gas = MISSING_BOOL
      endif

      write(field_path, '(A,A)') trim(species_path), '/is_aerosol'
      call safe_yaml_get_logical(yaml_root, trim(field_path), temp_logical, yaml_rc)
      if (yaml_rc == 0) then
         species%is_aerosol = temp_logical
      else
         species%is_aerosol = MISSING_BOOL
      endif

      write(field_path, '(A,A)') trim(species_path), '/is_dust'
      call safe_yaml_get_logical(yaml_root, trim(field_path), temp_logical, yaml_rc)
      if (yaml_rc == 0) then
         species%is_dust = temp_logical
      else
         species%is_dust = MISSING_BOOL
      endif

      write(field_path, '(A,A)') trim(species_path), '/is_seasalt'
      call safe_yaml_get_logical(yaml_root, trim(field_path), temp_logical, yaml_rc)
      if (yaml_rc == 0) then
         species%is_seasalt = temp_logical
      else
         species%is_seasalt = MISSING_BOOL
      endif

      write(field_path, '(A,A)') trim(species_path), '/is_tracer'
      call safe_yaml_get_logical(yaml_root, trim(field_path), temp_logical, yaml_rc)
      if (yaml_rc == 0) then
         species%is_tracer = temp_logical
      else
         species%is_tracer = MISSING_BOOL
      endif

      write(field_path, '(A,A)') trim(species_path), '/is_drydep'
      call safe_yaml_get_logical(yaml_root, trim(field_path), temp_logical, yaml_rc)
      if (yaml_rc == 0) then
         species%is_drydep = temp_logical
      else
         species%is_drydep = MISSING_BOOL
      endif

      write(field_path, '(A,A)') trim(species_path), '/is_wetdep'
      call safe_yaml_get_logical(yaml_root, trim(field_path), temp_logical, yaml_rc)
      if (yaml_rc == 0) then
         species%is_wetdep = temp_logical
      else
         species%is_wetdep = MISSING_BOOL
      endif

      write(field_path, '(A,A)') trim(species_path), '/is_photolysis'
      call safe_yaml_get_logical(yaml_root, trim(field_path), temp_logical, yaml_rc)
      if (yaml_rc == 0) then
         species%is_photolysis = temp_logical
      else
         species%is_photolysis = MISSING_BOOL
      endif

      ! Load background concentration (optional)
      write(field_path, '(A,A)') trim(species_path), '/background_vv'
      call safe_yaml_get_real(yaml_root, trim(field_path), temp_real, yaml_rc)
      if (yaml_rc == 0) then
         species%BackgroundVV = temp_real
      else
         species%BackgroundVV = MISSING
      endif

      ! Set species as valid
      species%is_valid = .true.

      ! Print species information in a single line
      write(*, '(A,A,A,F6.1,A,L1,A,L1,A,L1,A,L1,A)') &
         'INFO: Loaded species "', trim(adjustl(species%short_name)), &
         '" (MW=', species%mw_g, ', gas=', species%is_gas, &
         ', aerosol=', species%is_aerosol, ', dust=', species%is_dust, &
         ', seasalt=', species%is_seasalt, ')'

   end subroutine load_species_properties

   !> \brief Get Mie data and initialize ChemState Mie arrays
   !!
   !! This subroutine reads the mie/files section from the configuration,
   !! loads the Mie data files, and initializes the ChemState Mie arrays.
   !!
   !! \param[in] this ConfigManager object
   !! \param[inout] chem_state ChemState object to initialize with Mie data
   !! \param[out] rc Return code
   subroutine config_manager_get_mie_data(this, chem_state, rc)
      use ChemState_Mod, only: ChemStateType
      implicit none
      class(ConfigManagerType), intent(in) :: this
      type(ChemStateType), intent(inout) :: chem_state
      integer, intent(out) :: rc

      ! Temporary variables for file information
      integer :: n_mie_files
      character(len=30) :: mie_names(100)
      character(len=255) :: mie_filenames(100)
      character(len=512) :: mie_full_paths(100)
      character(len=320) :: key_value_pairs(100)  ! key:value entries
      integer :: i, colon_pos
      character(len=255) :: mie_dir
      character(len=320) :: trimmed_entry

      rc = CC_SUCCESS
      n_mie_files = 0

      ! Get Mie directory
      mie_dir = this%config_data%file_paths%Mie_Directory
      if (len_trim(mie_dir) == 0) then
         mie_dir = './'
      end if

      ! Discover key:value pairs in mie/files section using modified function
      call discover_nested_yaml_section_items(this%config_file, 'mie/files', key_value_pairs, n_mie_files, rc, 'key_value_pairs')

      if (rc /= CC_SUCCESS .or. n_mie_files == 0) then
         write(*, '(A)') 'WARNING: No Mie files found in configuration. Skipping Mie data initialization.'
         rc = CC_SUCCESS  ! Don't fail, just skip Mie initialization
         return
      end if

      ! Parse each key:value pair inline (reusing existing colon parsing logic)
      do i = 1, n_mie_files
         trimmed_entry = trim(key_value_pairs(i))

         ! Find the colon separator (reusing existing pattern)
         colon_pos = index(trimmed_entry, ':')
         if (colon_pos == 0) then
            rc = CC_FAILURE
            return
         end if

         ! Extract Mie name (before colon) and filename (after colon)
         mie_names(i) = trim(trimmed_entry(1:colon_pos-1))
         mie_filenames(i) = trim(adjustl(trimmed_entry(colon_pos+1:)))

         ! Construct full path
         if (trim(mie_dir) == './') then
            mie_full_paths(i) = trim(mie_filenames(i))
         else
            mie_full_paths(i) = trim(mie_dir) // trim(mie_filenames(i))
         end if
      end do

      ! Initialize Mie data in ChemState directly
      call chem_state%init_mie_data(n_mie_files, mie_names(1:n_mie_files), mie_full_paths(1:n_mie_files), rc)

      if (rc == CC_SUCCESS) then
         write(*, '(A,I0,A)') 'INFO: Successfully initialized Mie data with ', n_mie_files, ' files'
      else
         write(*, '(A)') 'ERROR: Failed to initialize Mie data in ChemState'
      end if

   end subroutine config_manager_get_mie_data

   !> \brief Find category mapping for a specific category
   subroutine config_manager_find_category_mapping(this, category_name, category_mapping, rc)
      implicit none
      class(ConfigManagerType), intent(in) :: this
      character(len=*), intent(in) :: category_name
      type(EmissionCategoryMapping), intent(out) :: category_mapping
      integer, intent(out) :: rc

      integer :: i

      rc = CC_FAILURE

      if (.not. this%config_data%emission_mapping%is_loaded) then
         return
      endif

      do i = 1, this%config_data%emission_mapping%n_categories
         if (trim(this%config_data%emission_mapping%categories(i)%category_name) == trim(category_name)) then
            ! Copy the found mapping
            category_mapping = this%config_data%emission_mapping%categories(i)
            rc = CC_SUCCESS
            return
         end if
      end do

   end subroutine config_manager_find_category_mapping

   !> \brief Apply emission mapping to map emission field to chemical species
   subroutine config_manager_apply_emission_mapping(this, category_name, emission_field, &
      chemical_species, scaling_factors, rc)
      implicit none
      class(ConfigManagerType), intent(in) :: this
      character(len=*), intent(in) :: category_name
      character(len=*), intent(in) :: emission_field
      character(len=64), allocatable, intent(out) :: chemical_species(:)
      real(fp), allocatable, intent(out) :: scaling_factors(:)
      integer, intent(out) :: rc

      type(EmissionCategoryMapping) :: category_mapping
      integer :: i, j, n_mappings

      rc = CC_FAILURE

      ! Find the category
      call this%find_category_mapping(category_name, category_mapping, rc)
      if (rc /= CC_SUCCESS) return

      ! Find the emission field in this category
      do i = 1, category_mapping%n_emission_species
         if (trim(category_mapping%species_mappings(i)%emission_field) == trim(emission_field)) then
            n_mappings = category_mapping%species_mappings(i)%n_mappings

            if (n_mappings > 0) then
               ! Allocate output arrays
               allocate(chemical_species(n_mappings))
               allocate(scaling_factors(n_mappings))

               ! Copy mapping data
               do j = 1, n_mappings
                  chemical_species(j) = category_mapping%species_mappings(i)%map(j)
                  scaling_factors(j) = category_mapping%species_mappings(i)%scale(j)
               end do

               rc = CC_SUCCESS
            endif

            return
         end if
      end do

   end subroutine config_manager_apply_emission_mapping

   !> \brief Get emission mapping configuration for specific category
   function config_manager_get_emission_mapping_for_category(this, category_name) result(mapping_index)
      implicit none
      class(ConfigManagerType), intent(in) :: this
      character(len=*), intent(in) :: category_name
      integer :: mapping_index

      integer :: i

      mapping_index = 0

      if (.not. this%config_data%emission_mapping%is_loaded) then
         return
      endif

      do i = 1, this%config_data%emission_mapping%n_categories
         if (trim(this%config_data%emission_mapping%categories(i)%category_name) == trim(category_name)) then
            mapping_index = i
            return
         end if
      end do

   end function config_manager_get_emission_mapping_for_category

   !> \brief Load emission configuration from file
   subroutine config_manager_load_emission_config(this, filename, rc)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      type(yaml_node_t) :: emission_config
      logical :: file_exists, success
      character(len=256) :: emission_directory, scaling_method
      integer :: n_emission_sources, local_rc
      real(fp) :: global_scaling_factor

      rc = CC_SUCCESS

      ! Check if file exists
      inquire(file=trim(filename), exist=file_exists)
      if (.not. file_exists) then
         write(*, '(A,A)') 'WARNING: Emission configuration file not found: ', trim(filename)
         ! Don't fail - emissions may be optional
         return
      endif

      ! Load emission configuration file
      emission_config = yaml_load_file(filename)
      if (.not. c_associated(emission_config%ptr)) then
         rc = CC_FAILURE
         return
      endif

      ! Read basic emission settings
      success = yaml_get_string(emission_config, "emissions/data_directory", &
         emission_directory)
      if (success) then
         write(*, '(A,A)') 'INFO: Emission data directory: ', trim(emission_directory)
      endif

      success = yaml_get_string(emission_config, "emissions/scaling_method", &
         scaling_method)
      if (success) then
         write(*, '(A,A)') 'INFO: Emission scaling method: ', trim(scaling_method)
      endif

      call safe_yaml_get_real(emission_config, "emissions/global_scaling_factor", &
         global_scaling_factor, local_rc)
      if (local_rc == 0) then
         write(*, '(A,F8.3)') 'INFO: Global emission scaling factor: ', global_scaling_factor
      endif

      call safe_yaml_get_integer(emission_config, "emissions/n_sources", n_emission_sources, local_rc)
      if (local_rc == 0) then
         write(*, '(A,I0)') 'INFO: Number of emission sources: ', n_emission_sources
      endif

      ! TODO: Parse individual emission sources, species mapping,
      ! vertical distributions, temporal profiles, etc.
      ! This would require extending the ConfigDataType to store emission data

      ! Clean up
      call yaml_destroy_node(emission_config)

      write(*, '(A)') 'INFO: Emission configuration loaded successfully'

   end subroutine config_manager_load_emission_config

   !> \brief Load emission-to-species mapping configuration from file
   !! \param[in] chem_state Optional ChemState for immediate index resolution
   subroutine config_manager_load_emission_mapping(this, filename, rc, chem_state)
      use ChemState_Mod, only: ChemStateType
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc
      type(ChemStateType), optional, intent(in) :: chem_state

      type(yaml_node_t) :: mapping_config
      logical :: file_exists, success
      integer :: n_categories, n_species, i, j, n_maps, n_scales, k, species_idx
      integer :: n_resolved, n_unresolved
      real(fp) :: single_scale
      character(len=64), allocatable :: all_categories(:), all_species(:)
      character(len=64), allocatable :: emission_fields(:)
      integer :: n_fields

      rc = CC_SUCCESS

      ! Check if file exists
      inquire(file=trim(filename), exist=file_exists)
      if (.not. file_exists) then
         write(*, '(A,A)') 'WARNING: Emission mapping configuration file not found: ', trim(filename)
         ! Don't fail - mappings may be optional
         return
      endif

      ! Load emission mapping configuration file
      mapping_config = yaml_load_file(filename)
      if (.not. c_associated(mapping_config%ptr)) then
         rc = CC_FAILURE
         write(*, '(A)') 'ERROR: Failed to load emission mapping config file.'
         return
      endif

      ! Allocate categories array for discovery
      allocate(all_categories(100))  ! Maximum 100 categories

      ! Get all top-level keys - these are the category names (seasalt, dust, etc.)
      success = yaml_get_all_keys(mapping_config, all_categories, n_categories)
      if (.not. success) then
         write(*, '(A)') 'WARNING: yaml_get_all_keys failed - no categories found'
         call yaml_destroy_node(mapping_config)
         return
      end if

      ! Initialize emission mapping configuration
      call this%config_data%emission_mapping%init()
      this%config_data%emission_mapping%n_categories = n_categories
      allocate(this%config_data%emission_mapping%categories(n_categories))
      this%config_data%emission_mapping%config_file = filename

      ! Process each category
      do i = 1, n_categories
         call this%config_data%emission_mapping%categories(i)%init()
         this%config_data%emission_mapping%categories(i)%category_name = all_categories(i)

         ! Get all emission field names for this category using character-based discovery
         ! Since yaml_get_node is not available, we'll use the text-based approach
         allocate(all_species(100))  ! Temporary array, we'll reallocate later
         call discover_yaml_section_items(filename, trim(all_categories(i)), 'emission_fields', &
            all_species, n_species, rc)

         if (rc /= CC_SUCCESS .or. n_species == 0) then
            ! Skip this category if we can't find any emission fields
            if (allocated(all_species)) deallocate(all_species)
            write(*, '(A,A)') 'INFO: No emission fields found for category: ', trim(all_categories(i))
            cycle
         end if

         this%config_data%emission_mapping%categories(i)%n_emission_species = n_species
         allocate(this%config_data%emission_mapping%categories(i)%species_mappings(n_species))

         ! Process each emission species
         do j = 1, n_species
            call this%config_data%emission_mapping%categories(i)%species_mappings(j)%init()

            ! Get emission field name
            this%config_data%emission_mapping%categories(i)%species_mappings(j)%emission_field = all_species(j)

            ! Get long name using full path
            success = yaml_get_string(mapping_config, trim(all_categories(i))//'/'//trim(all_species(j))//'/long_name', &
               this%config_data%emission_mapping%categories(i)%species_mappings(j)%long_name)
            if (.not. success) then
               write(*, '(A,A,A,A)') 'DEBUG: Failed to get long_name for ', trim(all_categories(i)), '/', trim(all_species(j))
               this%config_data%emission_mapping%categories(i)%species_mappings(j)%long_name = 'No description'
            end if

            ! Get units using full path
            success = yaml_get_string(mapping_config, trim(all_categories(i))//'/'//trim(all_species(j))//'/units', &
               this%config_data%emission_mapping%categories(i)%species_mappings(j)%units)
            if (.not. success) then
               write(*, '(A,A,A,A)') 'WARNING: Failed to get units for ', trim(all_categories(i)), '/', trim(all_species(j))
               this%config_data%emission_mapping%categories(i)%species_mappings(j)%units = 'kg/m2/s'  ! Default units
            else
               !write(*, '(A,A,A,A,A,A)') 'DEBUG: Read units "', trim(this%config_data%emission_mapping%categories(i)%species_mappings(j)%units), &
               !      '" for ', trim(all_categories(i)), '/', trim(all_species(j))
            end if

            ! Try to read mapping arrays using YAML API with full key paths
            ! First allocate arrays with maximum possible size
            allocate(this%config_data%emission_mapping%categories(i)%species_mappings(j)%map(10))
            allocate(this%config_data%emission_mapping%categories(i)%species_mappings(j)%scale(10))
            allocate(this%config_data%emission_mapping%categories(i)%species_mappings(j)%index(10))

            success = yaml_get_string_array(mapping_config, trim(all_categories(i))//'/'//trim(all_species(j))//'/map', &
               this%config_data%emission_mapping%categories(i)%species_mappings(j)%map, n_maps)

            if (success .and. n_maps > 0) then
               ! Try to read scale array
               success = yaml_get_real_array(mapping_config, trim(all_categories(i))//'/'//trim(all_species(j))//'/scale', &
                  this%config_data%emission_mapping%categories(i)%species_mappings(j)%scale, n_scales)

               if (success .and. n_scales > 0) then
                  ! Make sure both arrays have the same size (use minimum)
                  n_maps = min(n_maps, n_scales)
                  this%config_data%emission_mapping%categories(i)%species_mappings(j)%n_mappings = n_maps
                  !write(*, '(A,I0,A,A,A,A)') 'DEBUG: Successfully read ', n_maps, ' mappings for ', &
                  !      trim(all_categories(i)), '/', trim(all_species(j))
               else
                  ! If scale reading fails, use default scale of 1.0
                  this%config_data%emission_mapping%categories(i)%species_mappings(j)%scale(1:n_maps) = 1.0_fp
                  this%config_data%emission_mapping%categories(i)%species_mappings(j)%n_mappings = n_maps
                  !write(*, '(A,I0,A,A,A,A)') 'DEBUG: Read ', n_maps, ' mappings with default scale for ', &
                  !      trim(all_categories(i)), '/', trim(all_species(j))
               end if

               ! Initialize indices to zero (will be resolved later when ChemState is available)
               this%config_data%emission_mapping%categories(i)%species_mappings(j)%index(1:n_maps) = 0
            else
               ! No mapping found, cleanup arrays
               deallocate(this%config_data%emission_mapping%categories(i)%species_mappings(j)%map)
               deallocate(this%config_data%emission_mapping%categories(i)%species_mappings(j)%scale)
               deallocate(this%config_data%emission_mapping%categories(i)%species_mappings(j)%index)
               this%config_data%emission_mapping%categories(i)%species_mappings(j)%n_mappings = 0
               write(*, '(A,A,A,A)') 'WARNING: No mapping found for ', trim(all_categories(i)), '/', trim(all_species(j))
            end if
         end do

         ! Clean up temporary species array for this category
         if (allocated(all_species)) deallocate(all_species)
      end do

      this%config_data%emission_mapping%is_loaded = .true.

      ! Resolve chemical species indices if ChemState is provided
      if (present(chem_state)) then
         n_resolved = 0
         n_unresolved = 0

         ! Loop through all categories and species mappings to resolve indices
         do i = 1, this%config_data%emission_mapping%n_categories
            do j = 1, this%config_data%emission_mapping%categories(i)%n_emission_species
               ! Resolve indices for each mapped species
               do k = 1, this%config_data%emission_mapping%categories(i)%species_mappings(j)%n_mappings
                  species_idx = chem_state%find_species(trim(this%config_data%emission_mapping%categories(i)%species_mappings(j)%map(k)))

                  if (species_idx > 0) then
                     this%config_data%emission_mapping%categories(i)%species_mappings(j)%index(k) = species_idx
                     n_resolved = n_resolved + 1
                  else
                     this%config_data%emission_mapping%categories(i)%species_mappings(j)%index(k) = 0
                     n_unresolved = n_unresolved + 1
                     write(*, '(A,A,A,A,A,A)') 'WARNING: Species "', &
                        trim(this%config_data%emission_mapping%categories(i)%species_mappings(j)%map(k)), &
                        '" from emission field "', &
                        trim(this%config_data%emission_mapping%categories(i)%species_mappings(j)%emission_field), &
                        '" in category "', trim(this%config_data%emission_mapping%categories(i)%category_name), '"'
                  endif
               end do
            end do
         end do

         write(*, '(A,I0,A,I0,A)') 'INFO: Resolved ', n_resolved, ' species indices, ', n_unresolved, ' unresolved'

         if (n_unresolved > 0) then
            write(*, '(A)') 'WARNING: Some emission mapping species were not found in ChemState'
         endif
      end if

      ! Clean up
      call yaml_destroy_node(mapping_config)
      if (allocated(all_categories)) deallocate(all_categories)
      if (allocated(all_species)) deallocate(all_species)

      write(*, '(A,I0,A)') 'INFO: Loaded emission mapping for ', n_categories, ' categories'

   end subroutine config_manager_load_emission_mapping

   !> \brief Generic YAML text parser for discovering items in a section
   !! This function reads the YAML file as plain text and extracts item names
   !! by parsing the structure line by line, looking for keys with ':' under the section
   !! \param[in] filename YAML configuration file path
   !! \param[in] section_name Section name to search for items
   !! \param[in] parse_mode Parsing mode: 'simple', 'emission_fields'
   !! \param[out] item_names Array to store discovered item names
   !! \param[out] n_items Number of items discovered
   !! \param[out] rc Return code
   subroutine discover_yaml_section_items(filename, section_name, parse_mode, item_names, n_items, rc)
      implicit none
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: section_name
      character(len=*), intent(in) :: parse_mode  ! 'simple' or 'emission_fields'
      character(len=64), intent(inout) :: item_names(:)
      integer, intent(out) :: n_items
      integer, intent(out) :: rc

      integer :: unit_num, io_stat, colon_pos, indent_level, section_indent
      character(len=256) :: line, trimmed_line, field_name
      logical :: in_section, found_section
      integer :: line_number, field_indent

      ! Emission field specific variables
      logical :: field_has_map, field_has_scale, field_added_on_exit
      character(len=64) :: current_field

      ! Variables for duplicate detection
      logical :: already_exists
      integer :: check_idx, i
      rc = CC_SUCCESS
      n_items = 0
      in_section = .false.
      found_section = .false.
      section_indent = -1
      line_number = 0

      ! Initialize emission field specific variables
      field_has_map = .false.
      field_has_scale = .false.
      current_field = ''
      field_indent = -1
      field_added_on_exit = .false.

      ! Open file for reading
      open(newunit=unit_num, file=trim(filename), status='old', action='read', iostat=io_stat)
      if (io_stat /= 0) then
         write(*, '(A,A)') 'ERROR: Cannot open configuration file: ', trim(filename)
         rc = CC_FAILURE
         return
      endif

      ! Read file line by line
      do
         read(unit_num, '(A)', iostat=io_stat) line
         if (io_stat /= 0) exit  ! End of file or error

         line_number = line_number + 1
         trimmed_line = trim(adjustl(line))

         ! Skip empty lines and comments
         if (len_trim(trimmed_line) == 0 .or. trimmed_line(1:1) == '#') cycle

         ! Calculate indentation level
         do indent_level = 1, len_trim(line)
            if (line(indent_level:indent_level) /= ' ') exit
         end do
         indent_level = indent_level - 1

         ! Look for lines with colons
         if (index(trimmed_line, ':') > 0) then
            colon_pos = index(trimmed_line, ':')
            field_name = trimmed_line(1:colon_pos-1)
            field_name = trim(adjustl(field_name))

            ! Check if we found our target section
            if (trim(field_name) == trim(section_name) .and. indent_level == 0) then
               in_section = .true.
               found_section = .true.
               section_indent = indent_level
               if (trim(parse_mode) == 'emission_fields') then
                  write(*, '(A,A,A,I0)') 'INFO: Found category "', trim(section_name), '" at line ', line_number
               endif
               cycle
            endif

            ! Process items within the section
            if (in_section .and. indent_level > section_indent) then

               if (trim(parse_mode) == 'simple') then
                  ! Simple mode: direct children are items
                  if (indent_level == section_indent + 2) then
                     if (n_items < size(item_names)) then
                        n_items = n_items + 1
                        item_names(n_items) = trim(field_name)
                        write(*, '(A,A)') 'INFO: Discovered item: ', trim(field_name)
                     endif
                  endif

               elseif (trim(parse_mode) == 'emission_fields') then
                  ! Emission fields mode: validate fields have both 'map' and 'scale'
                  if (indent_level == section_indent + 2) then
                     ! Direct child of category - potential emission field
                     ! Process previous field if it was complete
                     if (current_field /= '' .and. field_has_map .and. field_has_scale) then
                        if (n_items < size(item_names)) then
                           ! Check for duplicates before adding
                           already_exists = .false.
                           do check_idx = 1, n_items
                              if (trim(item_names(check_idx)) == trim(current_field)) then
                                 already_exists = .true.
                                 write(*, '(A,A)') 'WARNING: Duplicate detected and skipped: ', trim(current_field)
                                 exit
                              endif
                           end do

                           if (.not. already_exists) then
                              n_items = n_items + 1
                              item_names(n_items) = trim(current_field)
                              !write(*, '(A,A,A,I0)') 'INFO: Added emission field: ', trim(current_field), ' (new field trigger at line ', line_number, ')'
                              field_added_on_exit = .true.
                           endif
                        endif
                     endif

                     ! Start tracking new potential field
                     current_field = trim(field_name)
                     field_indent = indent_level
                     field_has_map = .false.
                     field_has_scale = .false.
                     field_added_on_exit = .false.

                  elseif (indent_level == field_indent + 2 .and. current_field /= '') then
                     ! Properties of current field
                     if (trim(field_name) == 'map') then
                        field_has_map = .true.
                     elseif (trim(field_name) == 'scale') then
                        field_has_scale = .true.
                     endif
                  endif
               endif

            elseif (in_section .and. indent_level <= section_indent) then
               ! We've left our section - process final field
               if (trim(parse_mode) == 'emission_fields' .and. current_field /= '' .and. &
                  field_has_map .and. field_has_scale .and. .not. field_added_on_exit) then
                  if (n_items < size(item_names)) then
                     ! Check for duplicates before adding
                     already_exists = .false.
                     do check_idx = 1, n_items
                        if (trim(item_names(check_idx)) == trim(current_field)) then
                           already_exists = .true.
                           write(*, '(A,A)') 'WARNING: Duplicate detected at exit and skipped: ', trim(current_field)
                           exit
                        endif
                     end do

                     if (.not. already_exists) then
                        n_items = n_items + 1
                        item_names(n_items) = trim(current_field)
                        !write(*, '(A,A,A,I0)') 'INFO: Added emission field: ', trim(current_field), ' (section exit trigger at line ', line_number, ')'
                        field_added_on_exit = .true.
                     endif
                  endif
               endif
               exit
            endif
         endif
      end do

      ! Handle pending emission field at end of file (only for the very last category)
      if (trim(parse_mode) == 'emission_fields' .and. current_field /= '' .and. &
         field_has_map .and. field_has_scale .and. .not. field_added_on_exit) then
         if (n_items < size(item_names)) then
            ! Check for duplicates before adding
            already_exists = .false.
            do check_idx = 1, n_items
               if (trim(item_names(check_idx)) == trim(current_field)) then
                  already_exists = .true.
                  write(*, '(A,A)') 'WARNING: Duplicate detected at EOF and skipped: ', trim(current_field)
                  exit
               endif
            end do

            if (.not. already_exists) then
               n_items = n_items + 1
               item_names(n_items) = trim(current_field)
               !write(*, '(A,A,A,I0)') 'INFO: Added emission field: ', trim(current_field), ' (end of file trigger at line ', line_number, ')'
            endif
         endif
      endif

      close(unit_num)

      ! Report results
      if (.not. found_section) then
         write(*, '(A,A)') 'INFO: Section not found: ', trim(section_name)
         rc = CC_SUCCESS  ! Not an error - section may not exist
      elseif (n_items == 0) then
         if (trim(parse_mode) == 'emission_fields') then
            write(*, '(A,A)') 'INFO: No emission fields found for category: ', trim(section_name)
         else
            write(*, '(A,A)') 'INFO: No items found in section: ', trim(section_name)
         endif
         rc = CC_SUCCESS
      else
         if (trim(parse_mode) == 'emission_fields') then
            write(*, '(A,I0,A,A)') 'INFO: Found ', n_items, ' emission fields in category: ', trim(section_name)
         else
            write(*, '(A,I0,A,A)') 'INFO: Found ', n_items, ' items in section: ', trim(section_name)
         endif
      endif

   end subroutine discover_yaml_section_items

   !> \brief Discover items in a nested YAML section (supports arbitrary depth paths)
   !!
   !! This function discovers direct child items in any nested YAML section.
   !! It supports arbitrary nesting depth and flexible indentation.
   !!
   !! Examples:
   !!   - "processes/extemis" -> finds anthro, point, fire, fengsha
   !!   - "processes/extemis/anthro" -> finds activate, scale_factor, source_file, etc.
   !!   - "mie/files" with mode 'key_value_pairs' -> finds "SS: opticsBands_SS.v3_3.RRTMG.nc"
   !!
   !! \param[in]    filename      YAML file to parse
   !! \param[in]    section_path  Nested path (e.g., "processes/extemis/anthro")
   !! \param[inout] item_names    Array to store discovered item names
   !! \param[out]   n_items       Number of items found
   !! \param[out]   rc           Return code
   !! \param[in]    search_mode   Optional: 'section_headers' (default) or 'key_value_pairs'
   subroutine discover_nested_yaml_section_items(filename, section_path, item_names, n_items, rc, search_mode)
      implicit none
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: section_path
      character(len=*), intent(inout) :: item_names(:)
      integer, intent(out) :: n_items
      integer, intent(out) :: rc
      character(len=*), optional, intent(in) :: search_mode

      integer :: unit_num, io_stat, colon_pos, indent_level
      character(len=256) :: line, trimmed_line, field_name, content_after_colon
      integer :: line_number, target_indent
      character(len=20) :: mode

      ! Path navigation variables
      character(len=64) :: path_components(10)  ! Support up to 10 levels deep
      integer :: n_path_components, current_depth
      integer :: path_indents(10)  ! Track actual indentation for each path level
      logical :: path_matched, in_target_section

      rc = CC_SUCCESS
      n_items = 0
      line_number = 0
      target_indent = -1

      ! Set search mode (default to section headers for backward compatibility)
      if (present(search_mode)) then
         mode = trim(search_mode)
      else
         mode = 'section_headers'
      end if

      ! Initialize path tracking
      current_depth = 0
      path_matched = .false.
      path_indents = -1
      in_target_section = .false.

      ! Parse the section path (e.g., "processes/extemis/anthro" -> ["processes", "extemis", "anthro"])
      call parse_yaml_path(section_path, path_components, n_path_components)

      if (n_path_components == 0) then
         write(*, '(A,A)') 'ERROR: Invalid section path: ', trim(section_path)
         rc = CC_FAILURE
         return
      end if

      ! Open file for reading
      open(newunit=unit_num, file=trim(filename), status='old', action='read', iostat=io_stat)
      if (io_stat /= 0) then
         write(*, '(A,A)') 'ERROR: Cannot open configuration file: ', trim(filename)
         rc = CC_FAILURE
         return
      endif

      ! Read file line by line
      do
         read(unit_num, '(A)', iostat=io_stat) line
         if (io_stat /= 0) exit  ! End of file or error

         line_number = line_number + 1
         trimmed_line = trim(adjustl(line))

         ! Skip empty lines and comments
         if (len_trim(trimmed_line) == 0 .or. trimmed_line(1:1) == '#') cycle

         ! Calculate indentation level
         do indent_level = 1, len_trim(line)
            if (line(indent_level:indent_level) /= ' ') exit
         end do
         indent_level = indent_level - 1

         ! Look for lines with colons
         if (index(trimmed_line, ':') > 0) then
            colon_pos = index(trimmed_line, ':')
            field_name = trimmed_line(1:colon_pos-1)
            field_name = trim(adjustl(field_name))

            ! Update path tracking with flexible indentation
            call update_flexible_path_tracking(field_name, indent_level, path_components, n_path_components, &
               current_depth, path_indents, path_matched, target_indent, in_target_section)

            ! Once we've found the target section, process child items directly
            ! Don't rely on complex path tracking for children
            if (in_target_section) then
               ! We're in the target section, process items at deeper indentation
               if (indent_level > target_indent) then
                  ! Check if this is a direct child (first level below target)
                  if (indent_level == target_indent + 2) then
                     content_after_colon = adjustl(trimmed_line(colon_pos+1:))

                     ! Handle different search modes
                     if (trim(mode) == 'key_value_pairs') then
                        ! Look for entries WITH content after colon (key:value pairs)
                        if (len_trim(content_after_colon) > 0) then
                           if (n_items < size(item_names)) then
                              n_items = n_items + 1
                              item_names(n_items) = trim(trimmed_line)  ! Store full line
                           endif
                        endif
                     else
                        ! Default: look for entries WITHOUT content after colon (section headers)
                        if (len_trim(content_after_colon) == 0) then
                           if (n_items < size(item_names)) then
                              n_items = n_items + 1
                              item_names(n_items) = trim(field_name)  ! Store just field name
                           endif
                        endif
                     endif
                  endif
               endif

               ! Check exit conditions for the target section
               if (indent_level <= target_indent .and. trim(field_name) /= trim(path_components(n_path_components))) then
                  ! We've moved to a different section at the same level or higher
                  exit
               endif
            endif
         endif
      end do

      ! Close file
      close(unit_num)

      ! Report results
      if (n_items == 0) then
         write(*, '(A,A)') 'INFO: Section not found or empty: ', trim(section_path)
         rc = CC_FAILURE
      else
         write(*, '(A,I0,A,A)') 'INFO: Found ', n_items, ' categories in section: ', trim(section_path)
      endif

   end subroutine discover_nested_yaml_section_items

   !> \brief Parse a YAML path into components (e.g., "processes/extemis" -> ["processes", "extemis"])
   subroutine parse_yaml_path(path_string, components, n_components)
      implicit none
      character(len=*), intent(in) :: path_string
      character(len=64), intent(out) :: components(:)
      integer, intent(out) :: n_components

      integer :: start_pos, end_pos, slash_pos
      character(len=256) :: remaining_path

      n_components = 0
      remaining_path = trim(path_string)

      do while (len_trim(remaining_path) > 0 .and. n_components < size(components))
         slash_pos = index(remaining_path, '/')

         if (slash_pos > 0) then
            ! Found a slash, extract component
            n_components = n_components + 1
            components(n_components) = trim(remaining_path(1:slash_pos-1))
            remaining_path = remaining_path(slash_pos+1:)
         else
            ! No more slashes, this is the last component
            n_components = n_components + 1
            components(n_components) = trim(remaining_path)
            exit
         endif
      end do

   end subroutine parse_yaml_path

   !> \brief Update path tracking state while parsing YAML with flexible indentation
   subroutine update_flexible_path_tracking(field_name, indent_level, path_components, n_path_components, &
      current_depth, path_indents, path_matched, target_indent, in_target_section)
      implicit none
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: indent_level, n_path_components
      character(len=64), intent(in) :: path_components(:)
      integer, intent(inout) :: current_depth, path_indents(:)
      logical, intent(inout) :: path_matched, in_target_section
      integer, intent(inout) :: target_indent

      ! Reset tracking if indentation level suggests we've backed out
      ! Only check if we have a valid current_depth
      if (current_depth >= 1 .and. current_depth <= size(path_indents)) then
         if (path_indents(current_depth) >= 0 .and. indent_level <= path_indents(current_depth)) then
            ! Back out to appropriate depth
            do while (current_depth >= 1 .and. current_depth <= size(path_indents) .and. &
               path_indents(current_depth) >= 0 .and. path_indents(current_depth) >= indent_level)
               current_depth = current_depth - 1
               ! Exit early if we've gone to root level
               if (current_depth < 1) exit
            end do
         endif
      endif

      ! Update path_matched based on current state
      path_matched = (current_depth == n_path_components)

      ! Only set in_target_section to true when we first find the target
      ! Don't reset it unless we explicitly exit the section
      if (path_matched .and. current_depth >= 1 .and. current_depth <= size(path_indents)) then
         target_indent = path_indents(current_depth)
      endif

      ! Try to advance to the next path component
      if (current_depth >= 0 .and. current_depth < n_path_components) then
         if (current_depth + 1 <= n_path_components .and. &
            trim(field_name) == trim(path_components(current_depth + 1))) then
            current_depth = current_depth + 1
            if (current_depth >= 1 .and. current_depth <= size(path_indents)) then
               path_indents(current_depth) = indent_level

               ! Check if we've found our complete target path
               if (current_depth == n_path_components) then
                  path_matched = .true.
                  in_target_section = .true.
                  target_indent = indent_level
                  !write(*, '(A,A,A,I0)') 'INFO: Found target section: ', trim(field_name), ' at depth ', current_depth
               endif
            endif
         endif
      endif
   end subroutine update_flexible_path_tracking

   !> \brief Update runtime configuration from loaded configs
   subroutine config_manager_update_runtime_from_configs(this, rc)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Placeholder implementation - would update runtime settings

   end subroutine config_manager_update_runtime_from_configs

   !> \brief Initialize EmisSpeciesMappingEntry
   subroutine emis_species_mapping_init(this)
      implicit none
      class(EmisSpeciesMappingEntry), intent(inout) :: this

      this%emission_field = ''
      this%long_name = ''
      this%units = ''
      this%n_mappings = 0
      if (allocated(this%map)) deallocate(this%map)
      if (allocated(this%scale)) deallocate(this%scale)
      if (allocated(this%index)) deallocate(this%index)
      this%is_active = .true.

   end subroutine emis_species_mapping_init

   !> \brief Cleanup EmisSpeciesMappingEntry
   subroutine emis_species_mapping_cleanup(this)
      implicit none
      class(EmisSpeciesMappingEntry), intent(inout) :: this

      if (allocated(this%map)) deallocate(this%map)
      if (allocated(this%scale)) deallocate(this%scale)
      if (allocated(this%index)) deallocate(this%index)
      this%n_mappings = 0

   end subroutine emis_species_mapping_cleanup

   !> \brief Copy EmisSpeciesMappingEntry
   subroutine emis_species_mapping_copy(this, other)
      implicit none
      class(EmisSpeciesMappingEntry), intent(inout) :: this
      type(EmisSpeciesMappingEntry), intent(in) :: other

      integer :: i

      call this%cleanup()

      this%emission_field = other%emission_field
      this%long_name = other%long_name
      this%units = other%units
      this%n_mappings = other%n_mappings
      this%is_active = other%is_active

      if (other%n_mappings > 0) then
         allocate(this%map(other%n_mappings))
         allocate(this%scale(other%n_mappings))
         if (allocated(other%index)) then
            allocate(this%index(other%n_mappings))
         end if

         do i = 1, other%n_mappings
            this%map(i) = other%map(i)
            this%scale(i) = other%scale(i)
            if (allocated(this%index) .and. allocated(other%index)) then
               this%index(i) = other%index(i)
            end if
         end do
      end if

   end subroutine emis_species_mapping_copy

   !> \brief Initialize EmissionCategoryMapping
   subroutine emis_category_mapping_init(this)
      implicit none
      class(EmissionCategoryMapping), intent(inout) :: this

      integer :: i

      this%category_name = ''
      if (allocated(this%species_mappings)) then
         do i = 1, size(this%species_mappings)
            call this%species_mappings(i)%cleanup()
         end do
         deallocate(this%species_mappings)
      end if
      this%n_emission_species = 0
      this%is_active = .true.

   end subroutine emis_category_mapping_init

   !> \brief Cleanup EmissionCategoryMapping
   subroutine emis_category_mapping_cleanup(this)
      implicit none
      class(EmissionCategoryMapping), intent(inout) :: this

      integer :: i

      if (allocated(this%species_mappings)) then
         do i = 1, size(this%species_mappings)
            call this%species_mappings(i)%cleanup()
         end do
         deallocate(this%species_mappings)
      end if
      this%n_emission_species = 0

   end subroutine emis_category_mapping_cleanup

   !> \brief Copy EmissionCategoryMapping
   subroutine emis_category_mapping_copy(this, other)
      implicit none
      class(EmissionCategoryMapping), intent(inout) :: this
      type(EmissionCategoryMapping), intent(in) :: other

      integer :: i

      call this%cleanup()

      this%category_name = other%category_name
      this%n_emission_species = other%n_emission_species
      this%is_active = other%is_active

      if (other%n_emission_species > 0 .and. allocated(other%species_mappings)) then
         allocate(this%species_mappings(other%n_emission_species))

         do i = 1, other%n_emission_species
            call this%species_mappings(i)%copy(other%species_mappings(i))
         end do
      end if

   end subroutine emis_category_mapping_copy

   !> \brief Initialize EmissionMappingConfig
   subroutine emis_mapping_config_init(this)
      implicit none
      class(EmissionMappingConfig), intent(inout) :: this

      integer :: i

      this%config_file = ''
      if (allocated(this%categories)) then
         do i = 1, size(this%categories)
            call this%categories(i)%cleanup()
         end do
         deallocate(this%categories)
      end if
      this%n_categories = 0
      this%is_loaded = .false.

   end subroutine emis_mapping_config_init

   !> \brief Cleanup EmissionMappingConfig
   subroutine emis_mapping_config_cleanup(this)
      implicit none
      class(EmissionMappingConfig), intent(inout) :: this

      integer :: i

      if (allocated(this%categories)) then
         do i = 1, size(this%categories)
            call this%categories(i)%cleanup()
         end do
         deallocate(this%categories)
      end if
      this%n_categories = 0
      this%is_loaded = .false.

   end subroutine emis_mapping_config_cleanup

   !> \brief Load run phases configuration from YAML
   !!
   !! Parses the run_phases section from the loaded YAML configuration and populates
   !! the ConfigDataType%run_phases and ConfigDataType%run_phase_processes arrays.
   !!
   !! Expected YAML format:
   !! run_phases:
   !!   phase_name:
   !!     process1
   !!     process2
   !!     process3
   !!
   !! \param[out] rc Return code (CC_SUCCESS/CC_FAILURE)
   subroutine config_manager_load_run_phases(this, rc)
      implicit none
      class(ConfigManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      character(len=64), allocatable :: process_array(:)
      character(len=64), allocatable :: discovered_phases(:)
      character(len=64), allocatable :: unique_processes(:)  ! Track unique process names
      integer, allocatable :: unique_process_indices(:)      ! Map unique process names to indices
      character(len=64) :: phase_name, process_name, test_value
      integer :: phase_idx, process_idx, num_phases, num_processes
      integer :: total_processes, n_discovered_phases, global_process_idx
      integer :: n_unique_processes, unique_idx
      logical :: has_run_phases, has_processes, success, process_found, is_duplicate
      character(len=256) :: process_scheme, temp_string
      logical :: temp_logical
      integer :: temp_integer, valid_phases

      rc = CC_SUCCESS

      if (.not. this%is_loaded) then
         write(*,*) 'Warning: No YAML configuration loaded for run phases'
         rc = CC_FAILURE
         return
      endif

      ! Initialize arrays
      allocate(process_array(50))
      allocate(discovered_phases(20))
      allocate(unique_processes(50))
      allocate(unique_process_indices(50))
      n_unique_processes = 0

      ! Initialize run_phases_enabled to false
      this%config_data%run_phases_enabled = .false.

      ! Check if run_phases section exists by trying to discover phases
      call discover_yaml_section_items(this%config_file, 'run_phases', 'simple', discovered_phases, n_discovered_phases, rc)
      has_run_phases = (rc == CC_SUCCESS .and. n_discovered_phases > 0)

      if (has_run_phases) then
         write(*,*) 'Loading run phases from configuration...'

         if (n_discovered_phases == 0) then
            write(*,*) 'Warning: No phases found in run_phases section'
            deallocate(process_array, discovered_phases, unique_processes, unique_process_indices)
            return
         endif

         ! Set run_phases_enabled to true since we found phases
         this%config_data%run_phases_enabled = .true.

         ! First pass: count unique processes and valid phases
         n_unique_processes = 0
         valid_phases = 0

         do phase_idx = 1, n_discovered_phases
            phase_name = trim(discovered_phases(phase_idx))

            ! Get processes for this phase - try new nested structure first
            success = yaml_get_string_array(this%yaml_data, 'run_phases/' // trim(phase_name) // '/processes', process_array, num_processes)

            if (.not. success) then
               ! Try old direct structure for backward compatibility
               success = yaml_get_string_array(this%yaml_data, 'run_phases/' // trim(phase_name), process_array, num_processes)

               if (.not. success) then
                  ! Try as a single string and parse it manually (new structure)
                  success = yaml_get_string(this%yaml_data, 'run_phases/' // trim(phase_name) // '/processes', test_value)
                  if (.not. success) then
                     ! Try old structure as single string
                     success = yaml_get_string(this%yaml_data, 'run_phases/' // trim(phase_name), test_value)
                  endif

                  if (success) then
                     call parse_space_separated_string(test_value, process_array, num_processes)
                     if (num_processes > 0) then
                        success = .true.
                     endif
                  endif
               endif
            endif

            if (success .and. num_processes > 0) then
               ! Validate that all processes exist in processes section before counting
               do process_idx = 1, num_processes
                  process_name = trim(process_array(process_idx))
                  ! Try to read a value from the process to validate it exists
                  success = yaml_get_string(this%yaml_data, 'processes/' // trim(process_name) // '/scheme', temp_string)
                  if (.not. success) then
                     ! If scheme doesn't exist, try another field to double-check
                     success = yaml_get_logical(this%yaml_data, 'processes/' // trim(process_name) // '/activate', temp_logical)
                  endif
                  if (.not. success) then
                     write(*,'(A,A,A,A,A)') 'ERROR: Run phases process "', trim(process_name), &
                        '" listed in phase "', trim(phase_name), '" is not defined in the processes section of the configuration YAML file'
                     rc = CC_FAILURE
                     deallocate(process_array, discovered_phases, unique_processes, unique_process_indices)
                     return
                  endif
               end do

               valid_phases = valid_phases + 1

               ! Check for unique processes
               do process_idx = 1, num_processes
                  process_name = trim(process_array(process_idx))

                  ! Check if already in unique list
                  is_duplicate = .false.
                  do unique_idx = 1, n_unique_processes
                     if (trim(unique_processes(unique_idx)) == trim(process_name)) then
                        is_duplicate = .true.
                        exit
                     endif
                  end do

                  ! Add to unique list if not duplicate
                  if (.not. is_duplicate) then
                     n_unique_processes = n_unique_processes + 1
                     unique_processes(n_unique_processes) = trim(process_name)
                     unique_process_indices(n_unique_processes) = n_unique_processes
                  endif
               end do
            else
               write(*,'(A,A,A)') 'Info: Skipping phase "', trim(phase_name), '" - no processes found'
            endif
         end do

         write(*,'(A,I0,A,I0,A)') 'Found ', valid_phases, ' valid phases with ', n_unique_processes, ' unique processes'

         ! Allocate arrays based on valid phase count
         if (allocated(this%config_data%run_phases)) then
            deallocate(this%config_data%run_phases)
         endif
         allocate(this%config_data%run_phases(valid_phases))

         if (allocated(this%config_data%run_phase_processes)) then
            deallocate(this%config_data%run_phase_processes)
         endif
         allocate(this%config_data%run_phase_processes(n_unique_processes))

         ! Second pass: populate the data structures (only valid phases)
         valid_phases = 0

         do phase_idx = 1, n_discovered_phases
            phase_name = trim(discovered_phases(phase_idx))

            ! Get processes for this phase - try new nested structure first
            success = yaml_get_string_array(this%yaml_data, 'run_phases/' // trim(phase_name) // '/processes', process_array, num_processes)

            if (.not. success) then
               ! Try old direct structure for backward compatibility
               success = yaml_get_string_array(this%yaml_data, 'run_phases/' // trim(phase_name), process_array, num_processes)

               if (.not. success) then
                  ! Try as a single string and parse it manually (new structure)
                  success = yaml_get_string(this%yaml_data, 'run_phases/' // trim(phase_name) // '/processes', test_value)
                  if (.not. success) then
                     ! Try old structure as single string
                     success = yaml_get_string(this%yaml_data, 'run_phases/' // trim(phase_name), test_value)
                  endif

                  if (success) then
                     call parse_space_separated_string(test_value, process_array, num_processes)
                     if (num_processes > 0) then
                        success = .true.
                     endif
                  endif
               endif
            endif

            if (.not. success .or. num_processes == 0) then
               cycle  ! Skip phases with no processes
            endif

            ! Increment the valid phase counter
            valid_phases = valid_phases + 1

            ! Initialize phase metadata with flexible YAML reading
            call populate_phase_config(this, phase_name, this%config_data%run_phases(valid_phases), rc)
            if (rc /= CC_SUCCESS) then
               deallocate(process_array, discovered_phases, unique_processes, unique_process_indices)
               return
            endif
            this%config_data%run_phases(valid_phases)%num_processes = num_processes

            ! Allocate processes array for this phase
            if (allocated(this%config_data%run_phases(valid_phases)%processes)) then
               deallocate(this%config_data%run_phases(valid_phases)%processes)
            endif
            allocate(this%config_data%run_phases(valid_phases)%processes(num_processes))

            ! Process each process in this phase
            do process_idx = 1, num_processes
               process_name = trim(process_array(process_idx))

               ! Note: Process validation already done in first pass, no need to check again

               ! Find unique index for this process
               unique_idx = 0
               do global_process_idx = 1, n_unique_processes
                  if (trim(unique_processes(global_process_idx)) == trim(process_name)) then
                     unique_idx = unique_process_indices(global_process_idx)
                     exit
                  endif
               end do

               ! Initialize process in phase with flexible configuration reading
               call populate_process_config(this, process_name, phase_name, &
                  this%config_data%run_phases(valid_phases)%processes(process_idx), &
                  unique_idx, process_idx, rc)
               if (rc /= CC_SUCCESS) then
                  deallocate(process_array, discovered_phases, unique_processes, unique_process_indices)
                  return
               endif

               write(*,'(A,A,A,I0,A,A,A,A,A,I0,A)') 'Phase "', trim(phase_name), '" Process ', process_idx, &
                  ': ', trim(process_name), ' (scheme: ', &
                  trim(this%config_data%run_phases(valid_phases)%processes(process_idx)%scheme), &
                  ', index: ', unique_idx, ')'
            end do

            write(*,'(A,A,A,I0,A)') 'Completed phase "', trim(phase_name), '" with ', num_processes, ' processes'
         end do

         ! Populate global run_phase_processes array with unique processes
         do unique_idx = 1, n_unique_processes
            process_name = trim(unique_processes(unique_idx))

            call populate_process_config(this, process_name, '', &
               this%config_data%run_phase_processes(unique_idx), &
               unique_idx, unique_idx, rc)
            if (rc /= CC_SUCCESS) then
               deallocate(process_array, discovered_phases, unique_processes, unique_process_indices)
               return
            endif

            write(*,'(A,I0,A,A,A,A,A)') 'Global Process ', unique_idx, ': ', trim(process_name), &
               ' (scheme: ', trim(this%config_data%run_phase_processes(unique_idx)%scheme), ')'
         end do

         write(*,*) 'All run phases loaded successfully!'

      else
         ! No run_phases section found - try direct 'processes' section
         ! Check if processes section exists by trying to discover processes
         call discover_yaml_section_items(this%config_file, 'processes', 'simple', process_array, num_processes, rc)
         has_processes = (rc == CC_SUCCESS .and. num_processes > 0)
         if (has_processes) then
            write(*,*) 'No run_phases found, loading direct processes section...'

            ! Set run_phases_enabled to false for direct processes mode
            this%config_data%run_phases_enabled = .false.

            if (num_processes == 0) then
               write(*,*) 'Info: No processes found in direct processes section'
               deallocate(process_array, discovered_phases, unique_processes, unique_process_indices)
               return
            endif

            write(*,'(A,I0,A)') 'Found ', num_processes, ' processes in direct processes section'

            ! Allocate only run_phase_processes array (no phases structure)
            if (allocated(this%config_data%run_phase_processes)) then
               deallocate(this%config_data%run_phase_processes)
            endif
            allocate(this%config_data%run_phase_processes(num_processes))

            ! Populate run_phase_processes array directly
            do process_idx = 1, num_processes
               process_name = trim(process_array(process_idx))

               call populate_process_config(this, process_name, '', &
                  this%config_data%run_phase_processes(process_idx), &
                  process_idx, process_idx, rc)
               if (rc /= CC_SUCCESS) then
                  deallocate(process_array, discovered_phases, unique_processes, unique_process_indices)
                  return
               endif

               write(*,'(A,I0,A,A,A,A,A)') 'Process ', process_idx, ': ', trim(process_name), &
                  ' (scheme: ', trim(this%config_data%run_phase_processes(process_idx)%scheme), ')'
            end do

            write(*,*) 'Direct processes loaded successfully!'
         else
            write(*,*) 'Info: No run_phases or processes section found in configuration'
         endif
      endif

      deallocate(process_array)
      deallocate(discovered_phases)
      deallocate(unique_processes)
      deallocate(unique_process_indices)

   end subroutine config_manager_load_run_phases

   !> \brief Helper subroutine to populate ProcessConfigType with flexible YAML reading
   !! Reads configuration values from YAML first, then assigns defaults if not found
   !! \param[in] config_mgr ConfigManager instance
   !! \param[in] process_name Process name to read configuration for
   !! \param[in] phase_name Phase name (empty for direct processes)
   !! \param[inout] process_config ProcessConfigType to populate
   !! \param[in] process_index Global process index to assign
   !! \param[in] local_priority Local priority within phase/global list
   !! \param[out] rc Return code
   subroutine populate_process_config(config_mgr, process_name, phase_name, process_config, &
      process_index, local_priority, rc)
      implicit none
      class(ConfigManagerType), intent(in) :: config_mgr
      character(len=*), intent(in) :: process_name
      character(len=*), intent(in) :: phase_name
      type(ProcessConfigType), intent(inout) :: process_config
      integer, intent(in) :: process_index
      integer, intent(in) :: local_priority
      integer, intent(out) :: rc

      character(len=256) :: process_path, temp_string, gas_scheme, aero_scheme
      logical :: temp_logical, success, gas_success, aero_success
      integer :: temp_integer, local_rc

      rc = CC_SUCCESS

      ! Build base path for this process in YAML
      write(process_path, '(A,A)') 'processes/', trim(process_name)

      ! Always set basic required fields
      process_config%name = trim(process_name)
      process_config%process_index = process_index

      ! Read process_type from YAML or default to process name
      success = yaml_get_string(config_mgr%yaml_data, trim(process_path) // '/type', temp_string)
      if (success) then
         process_config%process_type = trim(temp_string)
      else
         process_config%process_type = trim(process_name)
      endif

      ! Read enabled from YAML or default to true
      success = yaml_get_logical(config_mgr%yaml_data, trim(process_path) // '/activate', temp_logical)
      if (success) then
         process_config%enabled = temp_logical
      else
         process_config%enabled = .true.
      endif

      ! Read priority from YAML or default to local priority
      call safe_yaml_get_integer(config_mgr%yaml_data, trim(process_path) // '/priority', temp_integer, local_rc)
      if (local_rc == 0) then
         process_config%priority = temp_integer
      else
         process_config%priority = local_priority
      endif

      ! Read timing from YAML or default to 'default'
      success = yaml_get_string(config_mgr%yaml_data, trim(process_path) // '/timing', temp_string)
      if (success) then
         process_config%timing = trim(temp_string)
      else
         process_config%timing = 'explicit'
      endif

      ! Read subcycling from YAML or default to 1
      call safe_yaml_get_integer(config_mgr%yaml_data, trim(process_path) // '/subcycling', temp_integer, local_rc)
      if (local_rc == 0) then
         process_config%subcycling = temp_integer
      else
         process_config%subcycling = 1
      endif

      ! Read scheme from YAML or use process name as default
      success = yaml_get_string(config_mgr%yaml_data, trim(process_path) // '/scheme', temp_string)
      if (success) then
         process_config%scheme = trim(temp_string)
      else
         ! Try to read separate gas_scheme and aero_scheme
         gas_success = yaml_get_string(config_mgr%yaml_data, trim(process_path) // '/gas_scheme', gas_scheme)
         aero_success = yaml_get_string(config_mgr%yaml_data, trim(process_path) // '/aero_scheme', aero_scheme)

         if (gas_success .and. aero_success) then
            ! Combine gas and aero schemes
            process_config%scheme = trim(gas_scheme) // ' (gas) & ' // trim(aero_scheme) // ' (aero)'
         elseif (gas_success) then
            ! Only gas scheme found
            process_config%scheme = trim(gas_scheme) // ' (gas)'
         elseif (aero_success) then
            ! Only aero scheme found
            process_config%scheme = trim(aero_scheme) // ' (aero)'
         else
            ! No scheme configuration found
            write(*,'(A,A,A)') 'Warning: Scheme of process "', trim(process_name), '" is not defined!'
            process_config%scheme = 'default'
         endif
      endif

   end subroutine populate_process_config

   !> \brief Helper subroutine to populate run phase metadata with flexible YAML reading
   !! Reads phase configuration values from YAML first, then assigns defaults if not found
   !! \param[in] config_mgr ConfigManager instance
   !! \param[in] phase_name Phase name to read configuration for
   !! \param[inout] phase_config RunPhaseType to populate
   !! \param[out] rc Return code
   subroutine populate_phase_config(config_mgr, phase_name, phase_config, rc)
      implicit none
      class(ConfigManagerType), intent(in) :: config_mgr
      character(len=*), intent(in) :: phase_name
      type(RunPhaseType), intent(inout) :: phase_config
      integer, intent(out) :: rc

      character(len=256) :: phase_path, temp_string
      logical :: success
      integer :: temp_integer, local_rc

      rc = CC_SUCCESS

      ! Build base path for this phase in YAML
      write(phase_path, '(A,A)') 'run_phases/', trim(phase_name)

      ! Always set the phase name
      phase_config%name = trim(phase_name)

      ! Read description from YAML or default to 'Phase: <phase_name>'
      success = yaml_get_string(config_mgr%yaml_data, trim(phase_path) // '/description', temp_string)
      if (success) then
         phase_config%description = trim(temp_string)
      else
         phase_config%description = 'Phase: ' // trim(phase_name)
      endif

      ! Read frequency from YAML or default to 'every timestep'
      success = yaml_get_string(config_mgr%yaml_data, trim(phase_path) // '/frequency', temp_string)
      if (success) then
         phase_config%frequency = trim(temp_string)
      else
         phase_config%frequency = 'every timestep'
      endif

      ! Read subcycling from YAML or default to 1
      call safe_yaml_get_integer(config_mgr%yaml_data, trim(phase_path) // '/subcycling', temp_integer, local_rc)
      if (local_rc == 0) then
         phase_config%subcycling = temp_integer
      else
         phase_config%subcycling = 1
      endif

   end subroutine populate_phase_config
   !> \brief Parse space-separated string into array elements
   subroutine parse_space_separated_string(input_string, output_array, num_elements)
      implicit none
      character(len=*), intent(in) :: input_string
      character(len=64), intent(out) :: output_array(:)
      integer, intent(out) :: num_elements

      character(len=len(input_string)) :: work_string
      integer :: pos, start_pos, str_len, i
      logical :: in_word

      num_elements = 0
      work_string = trim(adjustl(input_string))
      str_len = len_trim(work_string)

      if (str_len == 0) return

      ! Simple space-separated parsing
      start_pos = 1
      in_word = .false.

      do pos = 1, str_len + 1
         if (pos > str_len .or. work_string(pos:pos) == ' ') then
            ! End of word or end of string
            if (in_word) then
               num_elements = num_elements + 1
               if (num_elements <= size(output_array)) then
                  output_array(num_elements) = work_string(start_pos:pos-1)
               endif
               in_word = .false.
            endif
         else
            ! Start of new word
            if (.not. in_word) then
               start_pos = pos
               in_word = .true.
            endif
         endif
      end do

   end subroutine parse_space_separated_string

   !> \brief Remove duplicate strings from an array
   subroutine remove_duplicates_from_array(array, n_items)
      implicit none
      character(len=64), intent(inout) :: array(:)
      integer, intent(inout) :: n_items

      integer :: i, j, new_count
      character(len=64) :: temp_array(size(array))
      logical :: is_duplicate

      if (n_items <= 1) return

      new_count = 0
      do i = 1, n_items
         is_duplicate = .false.
         do j = 1, new_count
            if (trim(array(i)) == trim(temp_array(j))) then
               is_duplicate = .true.
               exit
            endif
         end do

         if (.not. is_duplicate) then
            new_count = new_count + 1
            temp_array(new_count) = array(i)
         endif
      end do

      ! Copy back to original array
      do i = 1, new_count
         array(i) = temp_array(i)
      end do
      n_items = new_count
   end subroutine remove_duplicates_from_array

end module ConfigManager_Mod
