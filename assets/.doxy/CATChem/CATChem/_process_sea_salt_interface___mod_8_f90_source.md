

# File ProcessSeaSaltInterface\_Mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**seasalt**](dir_3e6fd2e121e43ca7d4114b6c0b4e05b6.md) **>** [**ProcessSeaSaltInterface\_Mod.F90**](_process_sea_salt_interface___mod_8_f90.md)

[Go to the documentation of this file](_process_sea_salt_interface___mod_8_f90.md)


```Fortran


module processseasaltinterface_mod

   ! Core CATChem infrastructure
   use precision_mod, only: fp
   use processinterface_mod, only: processinterface, columnprocessinterface
   use statemanager_mod, only: statemanagertype
   use gridmanager_mod, only: gridmanagertype
   use error_mod, only: cc_success, cc_failure, cc_error, cc_warning, errormanagertype
   use diagnosticmanager_mod, only: diagnosticmanagertype
   use diagnosticinterface_mod, only: diagnosticregistrytype, diagnosticfieldtype, diagnosticdatatype
   use virtualcolumn_mod, only: virtualcolumntype, virtualmettype
   use constants, only: g0, airmw  ! Gravitational acceleration for tendency physics and air molecular weight for unit conversion

   ! Core utilities (leverage existing infrastructure)
   use configmanager_mod, only: configmanagertype
   use chemstate_mod, only: chemstatetype
   use metstate_mod, only: metstatetype

   ! Common utilities - unified configuration
   use seasaltcommon_mod, only: seasaltprocessconfig

   ! Scheme modules
   use seasaltscheme_gong97_mod, only: compute_gong97
   use seasaltscheme_gong03_mod, only: compute_gong03
   use seasaltscheme_geos12_mod, only: compute_geos12

   implicit none
   private

   public :: processseasaltinterface

   type, extends(columnprocessinterface) :: processseasaltinterface
      private

      ! Unified process configuration (bridges ConfigManager to process-specific config)
      type(SeaSaltProcessConfig), public :: process_config

      ! Process utilities (leverage core infrastructure)
      type(ChemStateType), pointer :: chem_state => null()
      type(MetStateType), pointer :: met_state => null()
      ! Note: state_manager pointer removed as it was never used

      ! Process-specific diagnostic indices (base class handles storage)
      integer :: diag_seasalt_mass_emission_total_idx = -1
      integer :: diag_seasalt_number_emission_total_idx = -1

      ! Column-level diagnostic storage for interfacing with DiagManager
      ! These are allocated per-column during processing and can be used to
      ! accumulate data for DiagManager updates
      real(fp), allocatable :: column_seasalt_mass_emission_total              ! scalar - per column
      real(fp), allocatable :: column_seasalt_number_emission_total              ! scalar - per column

      ! Scheme-specific diagnostic storage (shared across all schemes that use them)
      real(fp), allocatable :: column_seasalt_mass_emission_per_bin(:)           ! scheme-specific 1D: species - per column
      real(fp), allocatable :: column_seasalt_number_emission_per_bin(:)           ! scheme-specific 1D: species - per column

   contains
      ! Required ProcessInterface implementations
      procedure :: init => process_init
      procedure :: run => process_run
      procedure :: finalize => process_finalize
      procedure :: parse_process_config => parse_seasalt_config

      ! Required ColumnProcessInterface implementations
      procedure :: init_column_processing => init_column_processing
      procedure :: run_column => run_column
      procedure :: finalize_column_processing => finalize_column_processing

      ! ProcessInterface capability registration
      procedure :: get_required_met_fields => get_required_met_fields
      procedure :: get_required_diagnostic_fields => get_required_diagnostic_fields

      ! Public testing interface for scheme manipulation
      procedure :: set_scheme => set_seasalt_scheme
      procedure :: get_scheme => get_seasalt_scheme

      ! Process-specific implementations (column virtualization)
      procedure, private :: run_active_scheme_column
      procedure, private :: run_gong97_scheme_column
      procedure, private :: run_gong03_scheme_column
      procedure, private :: run_geos12_scheme_column

      ! Diagnostic procedures (override base class method)
      procedure :: register_diagnostics => register_and_allocate_diagnostics
      procedure, private :: register_and_allocate_diagnostics
      procedure, private :: calculate_and_update_diagnostics

   end type processseasaltinterface

contains

   subroutine process_init(this, container, rc)
      class(ProcessSeaSaltInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_manager

      rc = cc_success

      ! Get error manager
      error_manager => container%get_error_manager()

      ! Initialize column processing capabilities
      call this%init_column_processing(container, rc)
      if (rc /= cc_success) return

      ! Set process-specific name and info
      this%name = 'seasalt'
      this%version = '1.0.0'
      this%description = 'Process for computing sea salt aerosol emissions over ocean surfaces'

      ! Parse process-specific configuration using unified approach
      call this%parse_process_config(container, error_manager, rc)
      if (rc /= cc_success) then
         return
      end if

      ! Get state pointers from container (needed for species loading)
      this%chem_state => container%get_chem_state_ptr()
      this%met_state => container%get_met_state_ptr()

      ! Load species from ChemState based on is_seasalt property
      call this%process_config%load_species_from_chem_state(this%chem_state, error_manager)
      ! Note: Error handling managed by error_manager internally

      ! Map diagnostic species names to indices in the species array
      call this%process_config%map_diagnostic_species_indices(error_manager)

      ! Validate the configuration with StateManager
      call this%process_config%validate(container, error_manager)
      ! Note: validate doesn't return rc, but error_manager tracks errors

      ! Register diagnostics for this process (only if diagnostics enabled)
      call this%register_diagnostics(container, rc)
      if (rc /= cc_success) then
         return
      end if

      ! Mark process as initialized and active
      call this%activate()

   end subroutine process_init

   subroutine process_run(this, container, rc)
      class(ProcessSeaSaltInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      rc = cc_success

      ! Check if process is active
      if (.not. this%process_config%is_active) then
         return
      end if

      ! For ColumnProcessInterface processes, the ProcessManager handles column iteration
      ! and calls run_column() for each virtual column. This method is mainly a placeholder
      ! for any global 3D operations that need to happen before/after column processing.

      ! Currently no global 3D operations needed for seasalt process
      ! All processing happens in run_column() method

   end subroutine process_run

   subroutine process_finalize(this, rc)
      class(ProcessSeaSaltInterface), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success

      ! Finalize column processing
      call this%finalize_column_processing(rc)
      if (rc /= cc_success) return

      ! Deallocate diagnostic class members
      if (allocated(this%column_seasalt_mass_emission_total)) deallocate(this%column_seasalt_mass_emission_total)
      if (allocated(this%column_seasalt_number_emission_total)) deallocate(this%column_seasalt_number_emission_total)
      ! Deallocate scheme-specific diagnostic fields (only deallocate unique fields once)
      if (allocated(this%column_seasalt_mass_emission_per_bin)) deallocate(this%column_seasalt_mass_emission_per_bin)
      if (allocated(this%column_seasalt_number_emission_per_bin)) deallocate(this%column_seasalt_number_emission_per_bin)

      ! Finalize unified configuration
      call this%process_config%finalize()

   end subroutine process_finalize

   subroutine parse_seasalt_config(this, state_manager, error_manager, rc)
      class(ProcessSeaSaltInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: state_manager  ! Changed to inout for function call
      type(ErrorManagerType), intent(inout) :: error_manager
      integer, intent(out) :: rc

      type(ConfigManagerType), pointer :: config_manager

      rc = cc_success

      ! Get configuration manager from state manager
      config_manager => state_manager%get_config_ptr()
      if (.not. associated(config_manager)) then
         call error_manager%report_error(1003, &
            'ConfigManager not available from StateManager', rc)
         return
      end if

      ! Use the unified configuration loader from SeaSaltCommon_Mod
      ! This handles the complexity of parsing hierarchical YAML into process-specific types
      call this%process_config%load_from_config(config_manager, error_manager)
      ! Note: Error handling managed by error_manager internally

      ! Process is now configured - the unified config contains all scheme-specific settings

   end subroutine parse_seasalt_config

   !========================================================================
   ! Column Processing Interface Implementation
   !========================================================================

   subroutine init_column_processing(this, container, rc)
      class(ProcessSeaSaltInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      rc = cc_success

      ! Enable column processing and set batch size for optimal performance
      call this%enable_column_processing()
      call this%set_column_batch_size(50)  ! Process 50 columns at a time

      ! Any process-specific column processing setup would go here
      ! For seasalt, no additional setup is needed

   end subroutine init_column_processing

   subroutine run_column(this, column, container, rc)
      class(ProcessSeaSaltInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      rc = cc_success

      ! Check if process is active
      if (.not. this%process_config%is_active) return

      ! Delegate to the active scheme for column processing
      call this%run_active_scheme_column(column, rc)

      ! Calculate and update diagnostics if enabled
      if (this%process_config%seasalt_config%diagnostics .and. rc == cc_success) then
         call this%calculate_and_update_diagnostics(column, container, rc)
      end if

   end subroutine run_column

   subroutine finalize_column_processing(this, rc)
      class(ProcessSeaSaltInterface), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success

      ! Clean up column processing
      call this%disable_column_processing()

      ! Any process-specific cleanup would go here

   end subroutine finalize_column_processing

   subroutine run_active_scheme_column(this, column, rc)
      class(ProcessSeaSaltInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      integer, intent(out) :: rc

      rc = cc_success

      ! Delegate to appropriate scheme using unified config
      select case (trim(this%process_config%seasalt_config%scheme))
       case ('gong97')
         call this%run_gong97_scheme_column(column, rc)
       case ('gong03')
         call this%run_gong03_scheme_column(column, rc)
       case ('geos12')
         call this%run_geos12_scheme_column(column, rc)
       case default
         rc = cc_failure
      end select

   end subroutine run_active_scheme_column

   subroutine run_gong97_scheme_column(this, column, rc)
      class(ProcessSeaSaltInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      integer, intent(out) :: rc

      ! Local variables for scheme calculation
      type(VirtualMetType), pointer :: met => null()  ! Pointer to meteorological data
      ! Meteorological fields
      real(fp), allocatable :: frocean(:)
      real(fp), allocatable :: frseaice(:)
      real(fp), allocatable :: sst(:)
      real(fp), allocatable :: u10m(:)
      real(fp), allocatable :: v10m(:)
      ! Species properties
      real(fp), allocatable :: species_density(:)
      real(fp), allocatable :: species_radius(:)
      real(fp), allocatable :: species_lower_radius(:)
      real(fp), allocatable :: species_upper_radius(:)
      real(fp), allocatable :: species_conc(:,:)
      real(fp), allocatable :: species_tendencies(:,:)
      integer :: n_species, n_levels, n_chem, n_emis, i, k
      integer, allocatable :: species_indices(:)
      real(fp) :: dqa  ! Concentration change for additive tendencies
      real(fp) :: converter  ! Unit conversion factor for emissions

      rc = cc_success

      ! Get dimensions from virtual column
      n_levels = 1  ! Surface-only processing

      ! Get seasalt species information from process configuration
      n_species = this%process_config%seasalt_config%n_species
      if (n_species <= 0) then
         return
      end if

      ! Get species indices directly from configuration (pre-computed)
      allocate(species_indices(n_species))
      species_indices(1:n_species) = this%process_config%seasalt_config%species_indices(1:n_species)

      ! Allocate arrays
      allocate(species_conc(1, n_species))
      allocate(species_tendencies(1, n_species))
      ! Allocate meteorological field arrays based on field type and process configuration
      allocate(frocean(1))  ! Surface field - always scalar
      allocate(frseaice(1))  ! Surface field - always scalar
      allocate(sst(1))  ! Surface field - always scalar
      allocate(u10m(1))  ! Surface field - always scalar
      allocate(v10m(1))  ! Surface field - always scalar
      allocate(species_density(n_species))
      allocate(species_radius(n_species))
      allocate(species_lower_radius(n_species))
      allocate(species_upper_radius(n_species))
      species_tendencies = 0.0_fp

      ! Get meteorological data pointer from virtual column (VirtualMet pattern)
      met => column%get_met()

      ! Now allocate categorical fields using the met pointer dimensions

      ! Extract required fields from met pointer based on field type and processing mode
      frocean(1) = met%FROCEAN  ! Surface field - scalar access
      frseaice(1) = met%FRSEAICE  ! Surface field - scalar access
      sst(1) = met%SST  ! Surface field - scalar access
      u10m(1) = met%U10M  ! Surface field - scalar access
      v10m(1) = met%V10M  ! Surface field - scalar access

      ! Get species concentrations from virtual column
      ! Surface-only processing - get surface level concentrations
      do i = 1, n_species
         species_conc(1, i) = column%get_chem_field(species_indices(i), 1)
      end do

      ! Get species properties from configuration (pre-loaded during initialization)
      ! Use species properties from process configuration
      species_density(1:n_species) = this%process_config%seasalt_config%species_density(1:n_species)
      ! Use species properties from process configuration
      species_radius(1:n_species) = this%process_config%seasalt_config%species_radius(1:n_species)
      ! Use species properties from process configuration
      species_lower_radius(1:n_species) = this%process_config%seasalt_config%species_lower_radius(1:n_species)
      ! Use species properties from process configuration
      species_upper_radius(1:n_species) = this%process_config%seasalt_config%species_upper_radius(1:n_species)

      ! Call the science scheme with optional diagnostic parameters
      ! Note: gong97 uses the following diagnostic fields (if diagnostics enabled):
      ! - seasalt_mass_emission_total (Sea salt mass emission flux total)
      ! - seasalt_number_emission_total (Sea salt number emission flux total)
      if (this%process_config%seasalt_config%diagnostics) then
         ! Call with diagnostic outputs enabled
         call compute_gong97( &
            n_levels, &
            n_species, &
            this%process_config%gong97_config, &
            frocean(1), &
            frseaice(1), &
            sst(1), &
            u10m(1), &
            v10m(1)            , &
            species_density, &
            species_radius, &
            species_lower_radius, &
            species_upper_radius, &
            species_conc, &
            species_tendencies, &
            this%column_seasalt_mass_emission_total, &
            this%column_seasalt_number_emission_total, &
            this%column_seasalt_mass_emission_per_bin, &
            this%column_seasalt_number_emission_per_bin, &
            this%process_config%seasalt_config%diagnostic_species_id         )
      else
         ! Call without diagnostic outputs (optional parameters not passed)
         call compute_gong97( &
            n_levels, &
            n_species, &
            this%process_config%gong97_config, &
            frocean(1), &
            frseaice(1), &
            sst(1), &
            u10m(1), &
            v10m(1)            , &
            species_density, &
            species_radius, &
            species_lower_radius, &
            species_upper_radius, &
            species_conc, &
            species_tendencies &
            )
      end if

      ! Apply tendencies back to virtual column based on tendency_mode
      ! Surface-only processing - apply tendencies to surface level only
      do i = 1, n_species
         ! Additive tendency: convert emission flux to concentration change
         ! Step 1: Convert to mass mixing ratio change (kg/kg)
         dqa = species_tendencies(1, i) * this%get_timestep() * g0 / met%DELP(1)

         ! Step 2: Convert to final concentration units
         ! For gas species: convert kg/kg to ppmv
         ! For aerosol species: convert kg/kg to ug/kg
         if (this%chem_state%ChemSpecies(species_indices(i))%is_gas) then
            converter = airmw / this%chem_state%ChemSpecies(species_indices(i))%mw_g * 1.0e6_fp
         else
            converter = 1.0e9_fp
         end if
         dqa = dqa * converter

         call column%set_chem_field(1, species_indices(i), &
            species_conc(1, i) + dqa)

      end do

   end subroutine run_gong97_scheme_column

   subroutine run_gong03_scheme_column(this, column, rc)
      class(ProcessSeaSaltInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      integer, intent(out) :: rc

      ! Local variables for scheme calculation
      type(VirtualMetType), pointer :: met => null()  ! Pointer to meteorological data
      ! Meteorological fields
      real(fp), allocatable :: frocean(:)
      real(fp), allocatable :: frseaice(:)
      real(fp), allocatable :: sst(:)
      real(fp), allocatable :: u10m(:)
      real(fp), allocatable :: v10m(:)
      ! Species properties
      real(fp), allocatable :: species_density(:)
      real(fp), allocatable :: species_radius(:)
      real(fp), allocatable :: species_lower_radius(:)
      real(fp), allocatable :: species_upper_radius(:)
      real(fp), allocatable :: species_conc(:,:)
      real(fp), allocatable :: species_tendencies(:,:)
      integer :: n_species, n_levels, n_chem, n_emis, i, k
      integer, allocatable :: species_indices(:)
      real(fp) :: dqa  ! Concentration change for additive tendencies
      real(fp) :: converter  ! Unit conversion factor for emissions

      rc = cc_success

      ! Get dimensions from virtual column
      n_levels = 1  ! Surface-only processing

      ! Get seasalt species information from process configuration
      n_species = this%process_config%seasalt_config%n_species
      if (n_species <= 0) then
         return
      end if

      ! Get species indices directly from configuration (pre-computed)
      allocate(species_indices(n_species))
      species_indices(1:n_species) = this%process_config%seasalt_config%species_indices(1:n_species)

      ! Allocate arrays
      allocate(species_conc(1, n_species))
      allocate(species_tendencies(1, n_species))
      ! Allocate meteorological field arrays based on field type and process configuration
      allocate(frocean(1))  ! Surface field - always scalar
      allocate(frseaice(1))  ! Surface field - always scalar
      allocate(sst(1))  ! Surface field - always scalar
      allocate(u10m(1))  ! Surface field - always scalar
      allocate(v10m(1))  ! Surface field - always scalar
      allocate(species_density(n_species))
      allocate(species_radius(n_species))
      allocate(species_lower_radius(n_species))
      allocate(species_upper_radius(n_species))
      species_tendencies = 0.0_fp

      ! Get meteorological data pointer from virtual column (VirtualMet pattern)
      met => column%get_met()

      ! Now allocate categorical fields using the met pointer dimensions

      ! Extract required fields from met pointer based on field type and processing mode
      frocean(1) = met%FROCEAN  ! Surface field - scalar access
      frseaice(1) = met%FRSEAICE  ! Surface field - scalar access
      sst(1) = met%SST  ! Surface field - scalar access
      u10m(1) = met%U10M  ! Surface field - scalar access
      v10m(1) = met%V10M  ! Surface field - scalar access

      ! Get species concentrations from virtual column
      ! Surface-only processing - get surface level concentrations
      do i = 1, n_species
         species_conc(1, i) = column%get_chem_field(species_indices(i), 1)
      end do

      ! Get species properties from configuration (pre-loaded during initialization)
      ! Use species properties from process configuration
      species_density(1:n_species) = this%process_config%seasalt_config%species_density(1:n_species)
      ! Use species properties from process configuration
      species_radius(1:n_species) = this%process_config%seasalt_config%species_radius(1:n_species)
      ! Use species properties from process configuration
      species_lower_radius(1:n_species) = this%process_config%seasalt_config%species_lower_radius(1:n_species)
      ! Use species properties from process configuration
      species_upper_radius(1:n_species) = this%process_config%seasalt_config%species_upper_radius(1:n_species)

      ! Call the science scheme with optional diagnostic parameters
      ! Note: gong03 uses the following diagnostic fields (if diagnostics enabled):
      ! - seasalt_mass_emission_total (Sea salt mass emission flux total)
      ! - seasalt_number_emission_total (Sea salt number emission flux total)
      if (this%process_config%seasalt_config%diagnostics) then
         ! Call with diagnostic outputs enabled
         call compute_gong03( &
            n_levels, &
            n_species, &
            this%process_config%gong03_config, &
            frocean(1), &
            frseaice(1), &
            sst(1), &
            u10m(1), &
            v10m(1)            , &
            species_density, &
            species_radius, &
            species_lower_radius, &
            species_upper_radius, &
            species_conc, &
            species_tendencies, &
            this%column_seasalt_mass_emission_total, &
            this%column_seasalt_number_emission_total, &
            this%column_seasalt_mass_emission_per_bin, &
            this%column_seasalt_number_emission_per_bin, &
            this%process_config%seasalt_config%diagnostic_species_id         )
      else
         ! Call without diagnostic outputs (optional parameters not passed)
         call compute_gong03( &
            n_levels, &
            n_species, &
            this%process_config%gong03_config, &
            frocean(1), &
            frseaice(1), &
            sst(1), &
            u10m(1), &
            v10m(1)            , &
            species_density, &
            species_radius, &
            species_lower_radius, &
            species_upper_radius, &
            species_conc, &
            species_tendencies &
            )
      end if

      ! Apply tendencies back to virtual column based on tendency_mode
      ! Surface-only processing - apply tendencies to surface level only
      do i = 1, n_species
         ! Additive tendency: convert emission flux to concentration change
         ! Step 1: Convert to mass mixing ratio change (kg/kg)
         dqa = species_tendencies(1, i) * this%get_timestep() * g0 / met%DELP(1)

         ! Step 2: Convert to final concentration units
         ! For gas species: convert kg/kg to ppmv
         ! For aerosol species: convert kg/kg to ug/kg
         if (this%chem_state%ChemSpecies(species_indices(i))%is_gas) then
            converter = airmw / this%chem_state%ChemSpecies(species_indices(i))%mw_g * 1.0e6_fp
         else
            converter = 1.0e9_fp
         end if
         dqa = dqa * converter

         call column%set_chem_field(1, species_indices(i), &
            species_conc(1, i) + dqa)

      end do

   end subroutine run_gong03_scheme_column

   subroutine run_geos12_scheme_column(this, column, rc)
      class(ProcessSeaSaltInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      integer, intent(out) :: rc

      ! Local variables for scheme calculation
      type(VirtualMetType), pointer :: met => null()  ! Pointer to meteorological data
      ! Meteorological fields
      real(fp), allocatable :: frocean(:)
      real(fp), allocatable :: frseaice(:)
      real(fp), allocatable :: sst(:)
      real(fp), allocatable :: ustar(:)
      ! Species properties
      real(fp), allocatable :: species_density(:)
      real(fp), allocatable :: species_radius(:)
      real(fp), allocatable :: species_lower_radius(:)
      real(fp), allocatable :: species_upper_radius(:)
      real(fp), allocatable :: species_conc(:,:)
      real(fp), allocatable :: species_tendencies(:,:)
      integer :: n_species, n_levels, n_chem, n_emis, i, k
      integer, allocatable :: species_indices(:)
      real(fp) :: dqa  ! Concentration change for additive tendencies
      real(fp) :: converter  ! Unit conversion factor for emissions

      rc = cc_success

      ! Get dimensions from virtual column
      n_levels = 1  ! Surface-only processing

      ! Get seasalt species information from process configuration
      n_species = this%process_config%seasalt_config%n_species
      if (n_species <= 0) then
         return
      end if

      ! Get species indices directly from configuration (pre-computed)
      allocate(species_indices(n_species))
      species_indices(1:n_species) = this%process_config%seasalt_config%species_indices(1:n_species)

      ! Allocate arrays
      allocate(species_conc(1, n_species))
      allocate(species_tendencies(1, n_species))
      ! Allocate meteorological field arrays based on field type and process configuration
      allocate(frocean(1))  ! Surface field - always scalar
      allocate(frseaice(1))  ! Surface field - always scalar
      allocate(sst(1))  ! Surface field - always scalar
      allocate(ustar(1))  ! Surface field - always scalar
      allocate(species_density(n_species))
      allocate(species_radius(n_species))
      allocate(species_lower_radius(n_species))
      allocate(species_upper_radius(n_species))
      species_tendencies = 0.0_fp

      ! Get meteorological data pointer from virtual column (VirtualMet pattern)
      met => column%get_met()

      ! Now allocate categorical fields using the met pointer dimensions

      ! Extract required fields from met pointer based on field type and processing mode
      frocean(1) = met%FROCEAN  ! Surface field - scalar access
      frseaice(1) = met%FRSEAICE  ! Surface field - scalar access
      sst(1) = met%SST  ! Surface field - scalar access
      ustar(1) = met%USTAR  ! Surface field - scalar access

      ! Get species concentrations from virtual column
      ! Surface-only processing - get surface level concentrations
      do i = 1, n_species
         species_conc(1, i) = column%get_chem_field(species_indices(i), 1)
      end do

      ! Get species properties from configuration (pre-loaded during initialization)
      ! Use species properties from process configuration
      species_density(1:n_species) = this%process_config%seasalt_config%species_density(1:n_species)
      ! Use species properties from process configuration
      species_radius(1:n_species) = this%process_config%seasalt_config%species_radius(1:n_species)
      ! Use species properties from process configuration
      species_lower_radius(1:n_species) = this%process_config%seasalt_config%species_lower_radius(1:n_species)
      ! Use species properties from process configuration
      species_upper_radius(1:n_species) = this%process_config%seasalt_config%species_upper_radius(1:n_species)

      ! Call the science scheme with optional diagnostic parameters
      ! Note: geos12 uses the following diagnostic fields (if diagnostics enabled):
      ! - seasalt_mass_emission_total (Sea salt mass emission flux total)
      ! - seasalt_number_emission_total (Sea salt number emission flux total)
      if (this%process_config%seasalt_config%diagnostics) then
         ! Call with diagnostic outputs enabled
         call compute_geos12( &
            n_levels, &
            n_species, &
            this%process_config%geos12_config, &
            frocean(1), &
            frseaice(1), &
            sst(1), &
            ustar(1)            , &
            species_density, &
            species_radius, &
            species_lower_radius, &
            species_upper_radius, &
            species_conc, &
            species_tendencies, &
            this%column_seasalt_mass_emission_total, &
            this%column_seasalt_number_emission_total, &
            this%column_seasalt_mass_emission_per_bin, &
            this%column_seasalt_number_emission_per_bin, &
            this%process_config%seasalt_config%diagnostic_species_id         )
      else
         ! Call without diagnostic outputs (optional parameters not passed)
         call compute_geos12( &
            n_levels, &
            n_species, &
            this%process_config%geos12_config, &
            frocean(1), &
            frseaice(1), &
            sst(1), &
            ustar(1)            , &
            species_density, &
            species_radius, &
            species_lower_radius, &
            species_upper_radius, &
            species_conc, &
            species_tendencies &
            )
      end if

      ! Apply tendencies back to virtual column based on tendency_mode
      ! Surface-only processing - apply tendencies to surface level only
      do i = 1, n_species
         ! Additive tendency: convert emission flux to concentration change
         ! Step 1: Convert to mass mixing ratio change (kg/kg)
         dqa = species_tendencies(1, i) * this%get_timestep() * g0 / met%DELP(1)

         ! Step 2: Convert to final concentration units
         ! For gas species: convert kg/kg to ppmv
         ! For aerosol species: convert kg/kg to ug/kg
         if (this%chem_state%ChemSpecies(species_indices(i))%is_gas) then
            converter = airmw / this%chem_state%ChemSpecies(species_indices(i))%mw_g * 1.0e6_fp
         else
            converter = 1.0e9_fp
         end if
         dqa = dqa * converter

         call column%set_chem_field(1, species_indices(i), &
            species_conc(1, i) + dqa)

      end do

   end subroutine run_geos12_scheme_column



   function get_required_met_fields(this) result(field_names)
      class(ProcessSeaSaltInterface), intent(in) :: this
      character(len=32), allocatable :: field_names(:)
      character(len=32), allocatable :: scheme_fields(:)
      character(len=32), allocatable :: process_fields(:)
      character(len=32), allocatable :: unique_fields(:)
      integer :: total_fields, scheme_count, process_count, i, j, unique_count
      logical :: is_duplicate

      ! Process-level required fields
      process_count = 1
      allocate(process_fields(process_count))
      process_fields(1) = 'DELP'

      ! Get scheme-specific fields based on selected scheme
      select case (trim(this%process_config%seasalt_config%scheme))
       case ('gong97')
         scheme_count = 5
         allocate(scheme_fields(scheme_count))
         scheme_fields(1) = 'FROCEAN'
         scheme_fields(2) = 'FRSEAICE'
         scheme_fields(3) = 'SST'
         scheme_fields(4) = 'U10M'
         scheme_fields(5) = 'V10M'
       case ('gong03')
         scheme_count = 5
         allocate(scheme_fields(scheme_count))
         scheme_fields(1) = 'FROCEAN'
         scheme_fields(2) = 'FRSEAICE'
         scheme_fields(3) = 'SST'
         scheme_fields(4) = 'U10M'
         scheme_fields(5) = 'V10M'
       case ('geos12')
         scheme_count = 4
         allocate(scheme_fields(scheme_count))
         scheme_fields(1) = 'FROCEAN'
         scheme_fields(2) = 'FRSEAICE'
         scheme_fields(3) = 'SST'
         scheme_fields(4) = 'USTAR'
       case default
         scheme_count = 0
         allocate(scheme_fields(0))
      end select

      ! Combine process-level and scheme-specific fields and remove duplicates
      ! First estimate maximum possible size (without duplicates)
      total_fields = process_count + scheme_count
      allocate(unique_fields(total_fields))
      unique_count = 0

      ! Add process-level fields first
      do i = 1, process_count
         unique_count = unique_count + 1
         unique_fields(unique_count) = process_fields(i)
      end do

      ! Add scheme-specific fields (check for duplicates)
      do i = 1, scheme_count
         is_duplicate = .false.
         do j = 1, unique_count
            if (trim(scheme_fields(i)) == trim(unique_fields(j))) then
               is_duplicate = .true.
               exit
            end if
         end do
         if (.not. is_duplicate) then
            unique_count = unique_count + 1
            unique_fields(unique_count) = scheme_fields(i)
         end if
      end do

      ! Allocate final result array with exact size
      allocate(field_names(unique_count))
      field_names(1:unique_count) = unique_fields(1:unique_count)

      ! Clean up temporary arrays
      if (allocated(unique_fields)) deallocate(unique_fields)
      if (allocated(process_fields)) deallocate(process_fields)
      if (allocated(scheme_fields)) deallocate(scheme_fields)

   end function get_required_met_fields

   function get_required_diagnostic_fields(this) result(field_names)
      class(ProcessSeaSaltInterface), intent(in) :: this
      character(len=64), allocatable :: field_names(:)

      allocate(field_names(2))
      field_names(1) = 'seasalt_mass_emission_total'
      field_names(2) = 'seasalt_number_emission_total'

   end function get_required_diagnostic_fields


   subroutine register_and_allocate_diagnostics(this, container, rc)
      use diagnosticinterface_mod, only: diagnosticregistrytype, diag_real_2d, diag_real_3d

      class(ProcessSeaSaltInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(DiagnosticManagerType), pointer :: diag_mgr
      type(DiagnosticRegistryType), pointer :: registry
      type(GridManagerType), pointer :: grid_mgr
      character(len=256) :: field_name  ! For constructing species-specific field names
      integer :: i  ! Loop variable for diagnostic species
      integer :: nx, ny, nz
      integer :: n_species
      integer :: dims_2d(2)
      integer :: dims_3d_species(3)

      rc = cc_success

      ! Only register diagnostics if enabled in config
      if (.not. this%process_config%seasalt_config%diagnostics) then
         return
      endif

      ! Get managers
      diag_mgr => container%get_diagnostic_manager()
      grid_mgr => container%get_grid_manager()

      ! Register this process with diagnostic manager (only once per process)
      call diag_mgr%register_process('seasalt', rc)
      if (rc /= cc_success) return

      ! Get the process registry for registering individual diagnostics
      call diag_mgr%get_process_registry('seasalt', registry, rc)
      if (rc /= cc_success) return

      ! Get grid dimensions
      call grid_mgr%get_shape(nx, ny, nz)
      dims_2d = [nx, ny]

      ! Get species count for 3D diagnostics
      n_species = this%process_config%seasalt_config%n_species
      dims_3d_species = [nx, ny, n_species]

      ! Register seasalt_mass_emission_total
      ! Register single field for non-species or level-only diagnostics
      call this%register_diagnostic_field(registry, 'seasalt_mass_emission_total', &
         'Sea salt mass emission flux total', &
         'kg/m2/s', diag_real_2d, &
         'seasalt', dims_2d, rc=rc)
      if (rc /= cc_success) return

      ! Register seasalt_number_emission_total
      ! Register single field for non-species or level-only diagnostics
      call this%register_diagnostic_field(registry, 'seasalt_number_emission_total', &
         'Sea salt number emission flux total', &
         'kg/m2/s', diag_real_2d, &
         'seasalt', dims_2d, rc=rc)
      if (rc /= cc_success) return

      ! Get selected scheme(s)
      ! Register scheme-specific diagnostics based on selected scheme
      select case (trim(this%process_config%seasalt_config%scheme))

       case ('gong97')
         ! Register gong97-specific diagnostics
         ! Register individual 2D fields for each diagnostic species (species-only diagnostics)
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            do i = 1, this%process_config%seasalt_config%n_diagnostic_species
               write(field_name, '(A,A,A)') 'seasalt_mass_emission_', &
                  trim(this%process_config%seasalt_config%diagnostic_species(i))
               call this%register_diagnostic_field(registry, trim(field_name), &
                  'Sea salt mass emission flux per bin', &
                  'kg/m2/s', diag_real_2d, &
                  'seasalt', dims_2d, rc=rc)
               if (rc /= cc_success) return
            end do
         end if
         if (rc /= cc_success) return

         ! Register individual 2D fields for each diagnostic species (species-only diagnostics)
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            do i = 1, this%process_config%seasalt_config%n_diagnostic_species
               write(field_name, '(A,A,A)') 'seasalt_number_emission_', &
                  trim(this%process_config%seasalt_config%diagnostic_species(i))
               call this%register_diagnostic_field(registry, trim(field_name), &
                  'Sea salt number emission flux per bin', &
                  'kg/m2/s', diag_real_2d, &
                  'seasalt', dims_2d, rc=rc)
               if (rc /= cc_success) return
            end do
         end if
         if (rc /= cc_success) return

       case ('gong03')
         ! Register gong03-specific diagnostics
         ! Register individual 2D fields for each diagnostic species (species-only diagnostics)
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            do i = 1, this%process_config%seasalt_config%n_diagnostic_species
               write(field_name, '(A,A,A)') 'seasalt_mass_emission_', &
                  trim(this%process_config%seasalt_config%diagnostic_species(i))
               call this%register_diagnostic_field(registry, trim(field_name), &
                  'Sea salt mass emission flux per bin', &
                  'kg/m2/s', diag_real_2d, &
                  'seasalt', dims_2d, rc=rc)
               if (rc /= cc_success) return
            end do
         end if
         if (rc /= cc_success) return

         ! Register individual 2D fields for each diagnostic species (species-only diagnostics)
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            do i = 1, this%process_config%seasalt_config%n_diagnostic_species
               write(field_name, '(A,A,A)') 'seasalt_number_emission_', &
                  trim(this%process_config%seasalt_config%diagnostic_species(i))
               call this%register_diagnostic_field(registry, trim(field_name), &
                  'Sea salt number emission flux per bin', &
                  'kg/m2/s', diag_real_2d, &
                  'seasalt', dims_2d, rc=rc)
               if (rc /= cc_success) return
            end do
         end if
         if (rc /= cc_success) return

       case ('geos12')
         ! Register geos12-specific diagnostics
         ! Register individual 2D fields for each diagnostic species (species-only diagnostics)
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            do i = 1, this%process_config%seasalt_config%n_diagnostic_species
               write(field_name, '(A,A,A)') 'seasalt_mass_emission_', &
                  trim(this%process_config%seasalt_config%diagnostic_species(i))
               call this%register_diagnostic_field(registry, trim(field_name), &
                  'Sea salt mass emission flux per bin', &
                  'kg/m2/s', diag_real_2d, &
                  'seasalt', dims_2d, rc=rc)
               if (rc /= cc_success) return
            end do
         end if
         if (rc /= cc_success) return

         ! Register individual 2D fields for each diagnostic species (species-only diagnostics)
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            do i = 1, this%process_config%seasalt_config%n_diagnostic_species
               write(field_name, '(A,A,A)') 'seasalt_number_emission_', &
                  trim(this%process_config%seasalt_config%diagnostic_species(i))
               call this%register_diagnostic_field(registry, trim(field_name), &
                  'Sea salt number emission flux per bin', &
                  'kg/m2/s', diag_real_2d, &
                  'seasalt', dims_2d, rc=rc)
               if (rc /= cc_success) return
            end do
         end if
         if (rc /= cc_success) return

       case default
         ! Unknown scheme - only register common diagnostics
         ! (already done above)

      end select

      ! Now allocate diagnostic class members after successful registration
      ! First, deallocate if already allocated (for scheme switching)
      if (allocated(this%column_seasalt_mass_emission_total)) deallocate(this%column_seasalt_mass_emission_total)
      if (allocated(this%column_seasalt_number_emission_total)) deallocate(this%column_seasalt_number_emission_total)
      if (allocated(this%column_seasalt_mass_emission_per_bin)) deallocate(this%column_seasalt_mass_emission_per_bin)
      if (allocated(this%column_seasalt_number_emission_per_bin)) deallocate(this%column_seasalt_number_emission_per_bin)

      ! Allocate and initialize scheme-specific diagnostic fields based on selected scheme
      ! For non-gas/aero differentiated process, allocate diagnostics normally

      ! Allocate common diagnostic fields (used by all schemes)
      ! Scalar diagnostic
      allocate(this%column_seasalt_mass_emission_total)
      this%column_seasalt_mass_emission_total = 0.0_fp
      ! Scalar diagnostic
      allocate(this%column_seasalt_number_emission_total)
      this%column_seasalt_number_emission_total = 0.0_fp

      ! Allocate scheme-specific diagnostics
      select case (trim(this%process_config%seasalt_config%scheme))
       case ('gong97')
         ! Scheme-specific diagnostics for gong97
         ! 1D diagnostic: diagnostic species only - allocated based on n_diagnostic_species
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            allocate(this%column_seasalt_mass_emission_per_bin(this%process_config%seasalt_config%n_diagnostic_species))
         end if
         if (allocated(this%column_seasalt_mass_emission_per_bin)) this%column_seasalt_mass_emission_per_bin = 0.0_fp
         ! 1D diagnostic: diagnostic species only - allocated based on n_diagnostic_species
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            allocate(this%column_seasalt_number_emission_per_bin(this%process_config%seasalt_config%n_diagnostic_species))
         end if
         if (allocated(this%column_seasalt_number_emission_per_bin)) this%column_seasalt_number_emission_per_bin = 0.0_fp
       case ('gong03')
         ! Scheme-specific diagnostics for gong03
         ! 1D diagnostic: diagnostic species only - allocated based on n_diagnostic_species
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            allocate(this%column_seasalt_mass_emission_per_bin(this%process_config%seasalt_config%n_diagnostic_species))
         end if
         if (allocated(this%column_seasalt_mass_emission_per_bin)) this%column_seasalt_mass_emission_per_bin = 0.0_fp
         ! 1D diagnostic: diagnostic species only - allocated based on n_diagnostic_species
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            allocate(this%column_seasalt_number_emission_per_bin(this%process_config%seasalt_config%n_diagnostic_species))
         end if
         if (allocated(this%column_seasalt_number_emission_per_bin)) this%column_seasalt_number_emission_per_bin = 0.0_fp
       case ('geos12')
         ! Scheme-specific diagnostics for geos12
         ! 1D diagnostic: diagnostic species only - allocated based on n_diagnostic_species
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            allocate(this%column_seasalt_mass_emission_per_bin(this%process_config%seasalt_config%n_diagnostic_species))
         end if
         if (allocated(this%column_seasalt_mass_emission_per_bin)) this%column_seasalt_mass_emission_per_bin = 0.0_fp
         ! 1D diagnostic: diagnostic species only - allocated based on n_diagnostic_species
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            allocate(this%column_seasalt_number_emission_per_bin(this%process_config%seasalt_config%n_diagnostic_species))
         end if
         if (allocated(this%column_seasalt_number_emission_per_bin)) this%column_seasalt_number_emission_per_bin = 0.0_fp
       case default
         ! No scheme-specific diagnostics for unknown schemes
      end select

   end subroutine register_and_allocate_diagnostics

   subroutine calculate_and_update_diagnostics(this, column, container, rc)
      class(ProcessSeaSaltInterface), intent(inout) :: this
      type(VirtualColumnType), intent(in) :: column
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: i_col, j_col  ! Column grid position
      integer :: i  ! Loop variable for diagnostic species
      character(len=256) :: field_name  ! For constructing species-specific field names
      character(len=64) :: selected_scheme

      rc = cc_success

      ! Skip if diagnostics not enabled
      if (.not. this%process_config%seasalt_config%diagnostics) return

      ! Get column grid position (x, y indices)
      call column%get_position(i_col, j_col)

      ! Update common diagnostic fields (used by all schemes)
      ! Scalar diagnostic field
      call this%update_scalar_diagnostic_column('seasalt_mass_emission_total', &
         this%column_seasalt_mass_emission_total, &
         i_col, j_col, container, rc)
      if (rc /= cc_success) return
      ! Scalar diagnostic field
      call this%update_scalar_diagnostic_column('seasalt_number_emission_total', &
         this%column_seasalt_number_emission_total, &
         i_col, j_col, container, rc)
      if (rc /= cc_success) return
      ! Update scheme-specific diagnostic fields based on active scheme
      select case (trim(this%process_config%seasalt_config%scheme))
       case ("gong97")
         ! Scheme-specific diagnostics for gong97
         ! Update individual species diagnostic fields (species-only diagnostics)
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            do i = 1, this%process_config%seasalt_config%n_diagnostic_species
               write(field_name, '(A,A,A)') 'seasalt_mass_emission_', &
                  trim(this%process_config%seasalt_config%diagnostic_species(i))
               call this%update_scalar_diagnostic_column(trim(field_name), &
                  this%column_seasalt_mass_emission_per_bin(i), &
                  i_col, j_col, container, rc)
               if (rc /= cc_success) return
            end do
         end if
         ! Update individual species diagnostic fields (species-only diagnostics)
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            do i = 1, this%process_config%seasalt_config%n_diagnostic_species
               write(field_name, '(A,A,A)') 'seasalt_number_emission_', &
                  trim(this%process_config%seasalt_config%diagnostic_species(i))
               call this%update_scalar_diagnostic_column(trim(field_name), &
                  this%column_seasalt_number_emission_per_bin(i), &
                  i_col, j_col, container, rc)
               if (rc /= cc_success) return
            end do
         end if
       case ("gong03")
         ! Scheme-specific diagnostics for gong03
         ! Update individual species diagnostic fields (species-only diagnostics)
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            do i = 1, this%process_config%seasalt_config%n_diagnostic_species
               write(field_name, '(A,A,A)') 'seasalt_mass_emission_', &
                  trim(this%process_config%seasalt_config%diagnostic_species(i))
               call this%update_scalar_diagnostic_column(trim(field_name), &
                  this%column_seasalt_mass_emission_per_bin(i), &
                  i_col, j_col, container, rc)
               if (rc /= cc_success) return
            end do
         end if
         ! Update individual species diagnostic fields (species-only diagnostics)
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            do i = 1, this%process_config%seasalt_config%n_diagnostic_species
               write(field_name, '(A,A,A)') 'seasalt_number_emission_', &
                  trim(this%process_config%seasalt_config%diagnostic_species(i))
               call this%update_scalar_diagnostic_column(trim(field_name), &
                  this%column_seasalt_number_emission_per_bin(i), &
                  i_col, j_col, container, rc)
               if (rc /= cc_success) return
            end do
         end if
       case ("geos12")
         ! Scheme-specific diagnostics for geos12
         ! Update individual species diagnostic fields (species-only diagnostics)
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            do i = 1, this%process_config%seasalt_config%n_diagnostic_species
               write(field_name, '(A,A,A)') 'seasalt_mass_emission_', &
                  trim(this%process_config%seasalt_config%diagnostic_species(i))
               call this%update_scalar_diagnostic_column(trim(field_name), &
                  this%column_seasalt_mass_emission_per_bin(i), &
                  i_col, j_col, container, rc)
               if (rc /= cc_success) return
            end do
         end if
         ! Update individual species diagnostic fields (species-only diagnostics)
         if (this%process_config%seasalt_config%n_diagnostic_species > 0) then
            do i = 1, this%process_config%seasalt_config%n_diagnostic_species
               write(field_name, '(A,A,A)') 'seasalt_number_emission_', &
                  trim(this%process_config%seasalt_config%diagnostic_species(i))
               call this%update_scalar_diagnostic_column(trim(field_name), &
                  this%column_seasalt_number_emission_per_bin(i), &
                  i_col, j_col, container, rc)
               if (rc /= cc_success) return
            end do
         end if
      end select

   end subroutine calculate_and_update_diagnostics


   subroutine set_seasalt_scheme(this, scheme_name)
      class(ProcessSeaSaltInterface), intent(inout) :: this
      character(len=*), intent(in) :: scheme_name

      this%process_config%seasalt_config%scheme = trim(scheme_name)

   end subroutine set_seasalt_scheme

   function get_seasalt_scheme(this) result(scheme_name)
      class(ProcessSeaSaltInterface), intent(in) :: this
      character(len=64) :: scheme_name

      scheme_name = trim(this%process_config%seasalt_config%scheme)

   end function get_seasalt_scheme

end module processseasaltinterface_mod
```


