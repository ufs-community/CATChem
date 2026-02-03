

# File ProcessDryDepInterface\_Mod.F90

[**File List**](files.md) **>** [**drydep**](dir_57fb5aa14ddb2cd518a6d90b65ffd000.md) **>** [**ProcessDryDepInterface\_Mod.F90**](_process_dry_dep_interface___mod_8_f90.md)

[Go to the documentation of this file](_process_dry_dep_interface___mod_8_f90.md)


```Fortran


module processdrydepinterface_mod

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
   use drydepcommon_mod, only: drydepprocessconfig

   ! Scheme modules
   use drydepscheme_wesely_mod, only: compute_wesely
   use drydepscheme_gocart_mod, only: compute_gocart
   use drydepscheme_zhang_mod, only: compute_zhang

   implicit none
   private

   public :: processdrydepinterface

   type, extends(columnprocessinterface) :: processdrydepinterface
      private

      ! Unified process configuration (bridges ConfigManager to process-specific config)
      type(DryDepProcessConfig), public :: process_config

      ! Process utilities (leverage core infrastructure)
      type(ChemStateType), pointer :: chem_state => null()
      type(MetStateType), pointer :: met_state => null()
      ! Note: state_manager pointer removed as it was never used

      ! Process-specific diagnostic indices (base class handles storage)
      integer :: diag_drydep_con_per_species_idx = -1
      integer :: diag_drydep_velocity_per_species_idx = -1

      ! Column-level diagnostic storage for interfacing with DiagManager
      ! These are allocated per-column during processing and can be used to
      ! accumulate data for DiagManager updates
      real(fp), allocatable :: column_drydep_con_per_species(:)           ! 1D: species - per column
      real(fp), allocatable :: column_drydep_velocity_per_species(:)           ! 1D: species - per column

      ! Scheme-specific diagnostic storage (shared across all schemes that use them)

   contains
      ! Required ProcessInterface implementations
      procedure :: init => process_init
      procedure :: run => process_run
      procedure :: finalize => process_finalize
      procedure :: parse_process_config => parse_drydep_config

      ! Required ColumnProcessInterface implementations
      procedure :: init_column_processing => init_column_processing
      procedure :: run_column => run_column
      procedure :: finalize_column_processing => finalize_column_processing

      ! ProcessInterface capability registration
      procedure :: get_required_met_fields => get_required_met_fields
      procedure :: get_required_diagnostic_fields => get_required_diagnostic_fields

      ! Public testing interface for scheme manipulation
      procedure :: set_scheme => set_drydep_scheme
      procedure :: get_scheme => get_drydep_scheme

      ! Process-specific implementations (column virtualization)
      procedure, private :: run_active_scheme_column
      procedure, private :: run_wesely_scheme_column
      procedure, private :: run_gocart_scheme_column
      procedure, private :: run_zhang_scheme_column

      ! Diagnostic procedures (override base class method)
      procedure :: register_diagnostics => register_and_allocate_diagnostics
      procedure, private :: register_and_allocate_diagnostics
      procedure, private :: calculate_and_update_diagnostics

   end type processdrydepinterface

contains

   subroutine process_init(this, container, rc)
      class(ProcessDryDepInterface), intent(inout) :: this
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
      this%name = 'drydep'
      this%version = '1.0.0'
      this%description = 'Process for computing dry deposition of gas and aerosol species'

      ! Parse process-specific configuration using unified approach
      call this%parse_process_config(container, error_manager, rc)
      if (rc /= cc_success) then
         return
      end if

      ! Get state pointers from container (needed for species loading)
      this%chem_state => container%get_chem_state_ptr()
      this%met_state => container%get_met_state_ptr()

      ! Load species from ChemState based on is_drydep property
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
      class(ProcessDryDepInterface), intent(inout) :: this
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

      ! Currently no global 3D operations needed for drydep process
      ! All processing happens in run_column() method

   end subroutine process_run

   subroutine process_finalize(this, rc)
      class(ProcessDryDepInterface), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success

      ! Finalize column processing
      call this%finalize_column_processing(rc)
      if (rc /= cc_success) return

      ! Deallocate diagnostic class members
      if (allocated(this%column_drydep_con_per_species)) deallocate(this%column_drydep_con_per_species)
      if (allocated(this%column_drydep_velocity_per_species)) deallocate(this%column_drydep_velocity_per_species)
      ! Deallocate scheme-specific diagnostic fields (only deallocate unique fields once)

      ! Finalize unified configuration
      call this%process_config%finalize()

   end subroutine process_finalize

   subroutine parse_drydep_config(this, state_manager, error_manager, rc)
      class(ProcessDryDepInterface), intent(inout) :: this
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

      ! Use the unified configuration loader from DryDepCommon_Mod
      ! This handles the complexity of parsing hierarchical YAML into process-specific types
      call this%process_config%load_from_config(config_manager, error_manager)
      ! Note: Error handling managed by error_manager internally

      ! Process is now configured - the unified config contains all scheme-specific settings

   end subroutine parse_drydep_config

   !========================================================================
   ! Column Processing Interface Implementation
   !========================================================================

   subroutine init_column_processing(this, container, rc)
      class(ProcessDryDepInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      rc = cc_success

      ! Enable column processing and set batch size for optimal performance
      call this%enable_column_processing()
      call this%set_column_batch_size(50)  ! Process 50 columns at a time

      ! Any process-specific column processing setup would go here
      ! For drydep, no additional setup is needed

   end subroutine init_column_processing

   subroutine run_column(this, column, container, rc)
      class(ProcessDryDepInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      rc = cc_success

      ! Check if process is active
      if (.not. this%process_config%is_active) return

      ! Delegate to the active scheme for column processing
      call this%run_active_scheme_column(column, rc)

      ! Calculate and update diagnostics if enabled
      if (this%process_config%drydep_config%diagnostics .and. rc == cc_success) then
         call this%calculate_and_update_diagnostics(column, container, rc)
      end if

   end subroutine run_column

   subroutine finalize_column_processing(this, rc)
      class(ProcessDryDepInterface), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success

      ! Clean up column processing
      call this%disable_column_processing()

      ! Any process-specific cleanup would go here

   end subroutine finalize_column_processing

   subroutine run_active_scheme_column(this, column, rc)
      class(ProcessDryDepInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      integer, intent(out) :: rc

      rc = cc_success

      ! For gas/aerosol differentiated processing, run schemes based on species type
      ! Run gas scheme for gas species
      select case (trim(this%process_config%drydep_config%gas_scheme))
       case ('wesely')
         call this%run_wesely_scheme_column(column, rc)
         if (rc /= cc_success) return
       case default
         rc = cc_failure
         return
      end select

      ! Run aerosol scheme for aerosol species
      select case (trim(this%process_config%drydep_config%aero_scheme))
       case ('gocart')
         call this%run_gocart_scheme_column(column, rc)
       case ('zhang')
         call this%run_zhang_scheme_column(column, rc)
       case default
         rc = cc_failure
      end select

   end subroutine run_active_scheme_column

   subroutine run_wesely_scheme_column(this, column, rc)
      class(ProcessDryDepInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      integer, intent(out) :: rc

      ! Local variables for scheme calculation
      type(VirtualMetType), pointer :: met => null()  ! Pointer to meteorological data
      ! Meteorological fields
      real(fp), allocatable :: bxheight(:)
      real(fp), allocatable :: cldfrc(:)
      real(fp), allocatable :: frlai(:)
      real(fp), allocatable :: frlanduse(:)
      integer, allocatable :: iland(:)
      logical, allocatable :: isice(:)
      logical, allocatable :: island(:)
      logical, allocatable :: issnow(:)
      real(fp), allocatable :: lat(:)
      real(fp), allocatable :: lon(:)
      character(len=255), allocatable :: lucname(:)
      real(fp), allocatable :: obk(:)
      real(fp), allocatable :: ps(:)
      real(fp), allocatable :: salinity(:)
      real(fp), allocatable :: suncosmid(:)
      real(fp), allocatable :: swgdn(:)
      real(fp), allocatable :: ts(:)
      real(fp), allocatable :: tskin(:)
      real(fp), allocatable :: tstep(:)
      real(fp), allocatable :: ustar(:)
      real(fp), allocatable :: z0(:)
      ! Species properties
      real(fp), allocatable :: species_mw_g(:)
      real(fp), allocatable :: species_dd_f0(:)
      real(fp), allocatable :: species_dd_hstar(:)
      real(fp), allocatable :: species_dd_DvzAerSnow(:)
      real(fp), allocatable :: species_dd_DvzMinVal_snow(:)
      real(fp), allocatable :: species_dd_DvzMinVal_land(:)
      real(fp), allocatable :: species_conc(:,:)
      real(fp), allocatable :: species_tendencies(:,:)
      integer :: n_species, n_levels, n_chem, n_emis, i, k
      integer, allocatable :: species_indices(:)
      real(fp) :: loss_fraction  ! Loss fraction for multiplicative tendencies

      rc = cc_success

      ! Get dimensions from virtual column
      n_levels = 1  ! Surface-only processing

      ! Get drydep species information from process configuration
      n_species = this%process_config%drydep_config%n_species
      if (n_species <= 0) then
         return
      end if

      ! Get species indices directly from configuration (pre-computed)
      allocate(species_indices(n_species))
      species_indices(1:n_species) = this%process_config%drydep_config%species_indices(1:n_species)

      ! Allocate arrays
      allocate(species_conc(1, n_species))
      allocate(species_tendencies(1, n_species))
      ! Allocate meteorological field arrays based on field type and process configuration
      allocate(bxheight(n_levels))  ! Atmospheric field - always n_levels
      allocate(cldfrc(1))  ! Surface field - always scalar



      allocate(isice(1))  ! Surface field - always scalar
      allocate(island(1))  ! Surface field - always scalar
      allocate(issnow(1))  ! Surface field - always scalar
      allocate(lat(1))  ! Surface field - always scalar
      allocate(lon(1))  ! Surface field - always scalar
      allocate(lucname(1))  ! Surface field - always scalar
      allocate(obk(1))  ! Surface field - always scalar
      allocate(ps(1))  ! Surface field - always scalar
      allocate(salinity(1))  ! Surface field - always scalar
      allocate(suncosmid(1))  ! Surface field - always scalar
      allocate(swgdn(1))  ! Surface field - always scalar
      allocate(ts(1))  ! Surface field - always scalar
      allocate(tskin(1))  ! Surface field - always scalar
      allocate(tstep(1))  ! Special timestep field - scalar
      allocate(ustar(1))  ! Surface field - always scalar
      allocate(z0(1))  ! Surface field - always scalar
      allocate(species_mw_g(n_species))
      allocate(species_dd_f0(n_species))
      allocate(species_dd_hstar(n_species))
      allocate(species_dd_dvzaersnow(n_species))
      allocate(species_dd_dvzminval_snow(n_species))
      allocate(species_dd_dvzminval_land(n_species))
      species_tendencies = 0.0_fp

      ! Get meteorological data pointer from virtual column (VirtualMet pattern)
      met => column%get_met()

      ! Now allocate categorical fields using the met pointer dimensions
      allocate(frlai(size(met%FRLAI)))  ! Categorical field - get size from met pointer
      allocate(frlanduse(size(met%FRLANDUSE)))  ! Categorical field - get size from met pointer
      allocate(iland(size(met%ILAND)))  ! Categorical field - get size from met pointer

      ! Extract required fields from met pointer based on field type and processing mode
      bxheight(1:n_levels) = met%BXHEIGHT(1:n_levels)  ! Atmospheric field - always n_levels
      cldfrc(1) = met%CLDFRC  ! Surface field - scalar access
      frlai(:) = met%FRLAI(:)  ! Categorical field - full dimension
      frlanduse(:) = met%FRLANDUSE(:)  ! Categorical field - full dimension
      iland(:) = met%ILAND(:)  ! Categorical field - full dimension
      isice(1) = met%IsIce  ! Surface field - scalar access
      island(1) = met%IsLand  ! Surface field - scalar access
      issnow(1) = met%IsSnow  ! Surface field - scalar access
      lat(1) = met%LAT  ! Surface field - scalar access
      lon(1) = met%LON  ! Surface field - scalar access
      lucname(1) = met%LUCNAME  ! Surface field - scalar access
      obk(1) = met%OBK  ! Surface field - scalar access
      ps(1) = met%PS  ! Surface field - scalar access
      salinity(1) = met%SALINITY  ! Surface field - scalar access
      suncosmid(1) = met%SUNCOSmid  ! Surface field - scalar access
      swgdn(1) = met%SWGDN  ! Surface field - scalar access
      ts(1) = met%TS  ! Surface field - scalar access
      tskin(1) = met%TSKIN  ! Surface field - scalar access
      tstep(1) = this%get_timestep()  ! Special timestep field - retrieved from ProcessInterface
      ustar(1) = met%USTAR  ! Surface field - scalar access
      z0(1) = met%Z0  ! Surface field - scalar access

      ! Get species concentrations from virtual column
      ! Surface-only processing - get surface level concentrations
      do i = 1, n_species
         species_conc(1, i) = column%get_chem_field(species_indices(i), 1)
      end do

      ! Get species properties from configuration (pre-loaded during initialization)
      ! Use species properties from process configuration
      species_mw_g(1:n_species) = this%process_config%drydep_config%species_mw_g(1:n_species)
      ! Use species properties from process configuration
      species_dd_f0(1:n_species) = this%process_config%drydep_config%species_dd_f0(1:n_species)
      ! Use species properties from process configuration
      species_dd_hstar(1:n_species) = this%process_config%drydep_config%species_dd_hstar(1:n_species)
      ! Use species properties from process configuration
      species_dd_dvzaersnow(1:n_species) = this%process_config%drydep_config%species_dd_DvzAerSnow(1:n_species)
      ! Use species properties from process configuration
      species_dd_dvzminval_snow(1:n_species) = this%process_config%drydep_config%species_dd_DvzMinVal_snow(1:n_species)
      ! Use species properties from process configuration
      species_dd_dvzminval_land(1:n_species) = this%process_config%drydep_config%species_dd_DvzMinVal_land(1:n_species)

      ! Call the science scheme with optional diagnostic parameters
      ! Note: wesely uses the following diagnostic fields (if diagnostics enabled):
      ! - drydep_con_per_species (Dry deposition concentration per species)
      ! - drydep_velocity_per_species (Dry deposition velocity)
      if (this%process_config%drydep_config%diagnostics) then
         ! Call with diagnostic outputs enabled
         call compute_wesely( &
            n_levels, &
            n_species, &
            this%process_config%wesely_config, &
            bxheight, &
            cldfrc(1), &
            frlai, &
            frlanduse, &
            iland, &
            isice(1), &
            island(1), &
            issnow(1), &
            lat(1), &
            lon(1), &
            lucname(1), &
            obk(1), &
            ps(1), &
            salinity(1), &
            suncosmid(1), &
            swgdn(1), &
            ts(1), &
            tskin(1), &
            tstep(1), &
            ustar(1), &
            z0(1)            , &
            species_mw_g, &
            species_dd_f0, &
            this%process_config%drydep_config%species_names, &
            species_dd_hstar, &
            species_dd_dvzaersnow, &
            species_dd_dvzminval_snow, &
            species_dd_dvzminval_land, &
            species_conc, &
            species_tendencies, &
            this%process_config%drydep_config%is_gas, &
            this%column_drydep_con_per_species, &
            this%column_drydep_velocity_per_species, &
            this%process_config%drydep_config%diagnostic_species_id         )
      else
         ! Call without diagnostic outputs (optional parameters not passed)
         call compute_wesely( &
            n_levels, &
            n_species, &
            this%process_config%wesely_config, &
            bxheight, &
            cldfrc(1), &
            frlai, &
            frlanduse, &
            iland, &
            isice(1), &
            island(1), &
            issnow(1), &
            lat(1), &
            lon(1), &
            lucname(1), &
            obk(1), &
            ps(1), &
            salinity(1), &
            suncosmid(1), &
            swgdn(1), &
            ts(1), &
            tskin(1), &
            tstep(1), &
            ustar(1), &
            z0(1)            , &
            species_mw_g, &
            species_dd_f0, &
            this%process_config%drydep_config%species_names, &
            species_dd_hstar, &
            species_dd_dvzaersnow, &
            species_dd_dvzminval_snow, &
            species_dd_dvzminval_land, &
            species_conc, &
            species_tendencies, &
            this%process_config%drydep_config%is_gas &
            )
      end if

      ! Apply tendencies back to virtual column based on tendency_mode
      ! Surface-only processing - apply tendencies to surface level only
      do i = 1, n_species
         ! Skip species that don't match scheme type (gas vs aerosol)
         if (.not. this%process_config%drydep_config%is_gas(i)) cycle
         ! Multiplicative tendency: new_conc = conc * (1 - loss_fraction)
         ! where loss_fraction = 1 - exp(-tendency * dt)
         loss_fraction = max(1.0_fp - exp(-species_tendencies(1, i) * this%get_timestep()), 0.0_fp)
         call column%set_chem_field(1, species_indices(i), &
            species_conc(1, i) * (1.0_fp - loss_fraction))
      end do

   end subroutine run_wesely_scheme_column

   subroutine run_gocart_scheme_column(this, column, rc)
      class(ProcessDryDepInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      integer, intent(out) :: rc

      ! Local variables for scheme calculation
      type(VirtualMetType), pointer :: met => null()  ! Pointer to meteorological data
      ! Meteorological fields
      real(fp), allocatable :: airden(:)
      real(fp), allocatable :: frlake(:)
      real(fp), allocatable :: gwettop(:)
      real(fp), allocatable :: hflux(:)
      integer, allocatable :: lwi(:)
      real(fp), allocatable :: pblh(:)
      real(fp), allocatable :: t(:)
      real(fp), allocatable :: tstep(:)
      real(fp), allocatable :: u10m(:)
      real(fp), allocatable :: ustar(:)
      real(fp), allocatable :: v10m(:)
      real(fp), allocatable :: z(:)
      real(fp), allocatable :: z0h(:)
      ! Species properties
      real(fp), allocatable :: species_density(:)
      real(fp), allocatable :: species_radius(:)
      logical, allocatable :: species_is_seasalt(:)
      real(fp), allocatable :: species_conc(:,:)
      real(fp), allocatable :: species_tendencies(:,:)
      integer :: n_species, n_levels, n_chem, n_emis, i, k
      integer, allocatable :: species_indices(:)
      real(fp) :: loss_fraction  ! Loss fraction for multiplicative tendencies

      rc = cc_success

      ! Get dimensions from virtual column
      n_levels = 1  ! Surface-only processing

      ! Get drydep species information from process configuration
      n_species = this%process_config%drydep_config%n_species
      if (n_species <= 0) then
         return
      end if

      ! Get species indices directly from configuration (pre-computed)
      allocate(species_indices(n_species))
      species_indices(1:n_species) = this%process_config%drydep_config%species_indices(1:n_species)

      ! Allocate arrays
      allocate(species_conc(1, n_species))
      allocate(species_tendencies(1, n_species))
      ! Allocate meteorological field arrays based on field type and process configuration
      allocate(airden(n_levels))  ! Atmospheric field - always n_levels
      allocate(frlake(1))  ! Surface field - always scalar
      allocate(gwettop(1))  ! Surface field - always scalar
      allocate(hflux(1))  ! Surface field - always scalar
      allocate(lwi(1))  ! Surface field - always scalar
      allocate(pblh(1))  ! Surface field - always scalar
      allocate(t(n_levels))  ! Atmospheric field - always n_levels
      allocate(tstep(1))  ! Special timestep field - scalar
      allocate(u10m(1))  ! Surface field - always scalar
      allocate(ustar(1))  ! Surface field - always scalar
      allocate(v10m(1))  ! Surface field - always scalar
      allocate(z(n_levels+1))  ! Edge field - always n_levels+1
      allocate(z0h(1))  ! Surface field - always scalar
      allocate(species_density(n_species))
      allocate(species_radius(n_species))
      allocate(species_is_seasalt(n_species))
      species_tendencies = 0.0_fp

      ! Get meteorological data pointer from virtual column (VirtualMet pattern)
      met => column%get_met()

      ! Now allocate categorical fields using the met pointer dimensions

      ! Extract required fields from met pointer based on field type and processing mode
      airden(1:n_levels) = met%AIRDEN(1:n_levels)  ! Atmospheric field - always n_levels
      frlake(1) = met%FRLAKE  ! Surface field - scalar access
      gwettop(1) = met%GWETTOP  ! Surface field - scalar access
      hflux(1) = met%HFLUX  ! Surface field - scalar access
      lwi(1) = met%LWI  ! Surface field - scalar access
      pblh(1) = met%PBLH  ! Surface field - scalar access
      t(1:n_levels) = met%T(1:n_levels)  ! Atmospheric field - always n_levels
      tstep(1) = this%get_timestep()  ! Special timestep field - retrieved from ProcessInterface
      u10m(1) = met%U10M  ! Surface field - scalar access
      ustar(1) = met%USTAR  ! Surface field - scalar access
      v10m(1) = met%V10M  ! Surface field - scalar access
      z(1:n_levels+1) = met%Z(1:n_levels+1)  ! Edge field - always n_levels+1
      z0h(1) = met%Z0H  ! Surface field - scalar access

      ! Get species concentrations from virtual column
      ! Surface-only processing - get surface level concentrations
      do i = 1, n_species
         species_conc(1, i) = column%get_chem_field(species_indices(i), 1)
      end do

      ! Get species properties from configuration (pre-loaded during initialization)
      ! Use species properties from process configuration
      species_density(1:n_species) = this%process_config%drydep_config%species_density(1:n_species)
      ! Use species properties from process configuration
      species_radius(1:n_species) = this%process_config%drydep_config%species_radius(1:n_species)
      ! Use species properties from process configuration
      species_is_seasalt(1:n_species) = this%process_config%drydep_config%species_is_seasalt(1:n_species)

      ! Call the science scheme with optional diagnostic parameters
      ! Note: gocart uses the following diagnostic fields (if diagnostics enabled):
      ! - drydep_con_per_species (Dry deposition concentration per species)
      ! - drydep_velocity_per_species (Dry deposition velocity)
      if (this%process_config%drydep_config%diagnostics) then
         ! Call with diagnostic outputs enabled
         call compute_gocart( &
            n_levels, &
            n_species, &
            this%process_config%gocart_config, &
            airden, &
            frlake(1), &
            gwettop(1), &
            hflux(1), &
            lwi(1), &
            pblh(1), &
            t, &
            tstep(1), &
            u10m(1), &
            ustar(1), &
            v10m(1), &
            z, &
            z0h(1)            , &
            species_density, &
            species_radius, &
            species_is_seasalt, &
            species_conc, &
            species_tendencies, &
            this%process_config%drydep_config%is_gas, &
            this%column_drydep_con_per_species, &
            this%column_drydep_velocity_per_species, &
            this%process_config%drydep_config%diagnostic_species_id         )
      else
         ! Call without diagnostic outputs (optional parameters not passed)
         call compute_gocart( &
            n_levels, &
            n_species, &
            this%process_config%gocart_config, &
            airden, &
            frlake(1), &
            gwettop(1), &
            hflux(1), &
            lwi(1), &
            pblh(1), &
            t, &
            tstep(1), &
            u10m(1), &
            ustar(1), &
            v10m(1), &
            z, &
            z0h(1)            , &
            species_density, &
            species_radius, &
            species_is_seasalt, &
            species_conc, &
            species_tendencies, &
            this%process_config%drydep_config%is_gas &
            )
      end if

      ! Apply tendencies back to virtual column based on tendency_mode
      ! Surface-only processing - apply tendencies to surface level only
      do i = 1, n_species
         ! Skip species that don't match scheme type (gas vs aerosol)
         if (this%process_config%drydep_config%is_gas(i)) cycle
         ! Multiplicative tendency: new_conc = conc * (1 - loss_fraction)
         ! where loss_fraction = 1 - exp(-tendency * dt)
         loss_fraction = max(1.0_fp - exp(-species_tendencies(1, i) * this%get_timestep()), 0.0_fp)
         call column%set_chem_field(1, species_indices(i), &
            species_conc(1, i) * (1.0_fp - loss_fraction))
      end do

   end subroutine run_gocart_scheme_column

   subroutine run_zhang_scheme_column(this, column, rc)
      class(ProcessDryDepInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      integer, intent(out) :: rc

      ! Local variables for scheme calculation
      type(VirtualMetType), pointer :: met => null()  ! Pointer to meteorological data
      ! Meteorological fields
      real(fp), allocatable :: bxheight(:)
      real(fp), allocatable :: frlanduse(:)
      integer, allocatable :: iland(:)
      logical, allocatable :: isice(:)
      logical, allocatable :: issnow(:)
      character(len=255), allocatable :: lucname(:)
      real(fp), allocatable :: obk(:)
      real(fp), allocatable :: ps(:)
      real(fp), allocatable :: rh(:)
      real(fp), allocatable :: ts(:)
      real(fp), allocatable :: tstep(:)
      real(fp), allocatable :: u10m(:)
      real(fp), allocatable :: ustar(:)
      real(fp), allocatable :: v10m(:)
      real(fp), allocatable :: z0(:)
      ! Species properties
      real(fp), allocatable :: species_mw_g(:)
      real(fp), allocatable :: species_radius(:)
      real(fp), allocatable :: species_density(:)
      real(fp), allocatable :: species_dd_hstar(:)
      real(fp), allocatable :: species_dd_DvzAerSnow(:)
      real(fp), allocatable :: species_dd_DvzMinVal_snow(:)
      real(fp), allocatable :: species_dd_DvzMinVal_land(:)
      real(fp), allocatable :: species_lower_radius(:)
      real(fp), allocatable :: species_upper_radius(:)
      logical, allocatable :: species_is_dust(:)
      logical, allocatable :: species_is_seasalt(:)
      real(fp), allocatable :: species_conc(:,:)
      real(fp), allocatable :: species_tendencies(:,:)
      integer :: n_species, n_levels, n_chem, n_emis, i, k
      integer, allocatable :: species_indices(:)
      real(fp) :: loss_fraction  ! Loss fraction for multiplicative tendencies

      rc = cc_success

      ! Get dimensions from virtual column
      n_levels = 1  ! Surface-only processing

      ! Get drydep species information from process configuration
      n_species = this%process_config%drydep_config%n_species
      if (n_species <= 0) then
         return
      end if

      ! Get species indices directly from configuration (pre-computed)
      allocate(species_indices(n_species))
      species_indices(1:n_species) = this%process_config%drydep_config%species_indices(1:n_species)

      ! Allocate arrays
      allocate(species_conc(1, n_species))
      allocate(species_tendencies(1, n_species))
      ! Allocate meteorological field arrays based on field type and process configuration
      allocate(bxheight(n_levels))  ! Atmospheric field - always n_levels


      allocate(isice(1))  ! Surface field - always scalar
      allocate(issnow(1))  ! Surface field - always scalar
      allocate(lucname(1))  ! Surface field - always scalar
      allocate(obk(1))  ! Surface field - always scalar
      allocate(ps(1))  ! Surface field - always scalar
      allocate(rh(n_levels))  ! Atmospheric field - always n_levels
      allocate(ts(1))  ! Surface field - always scalar
      allocate(tstep(1))  ! Special timestep field - scalar
      allocate(u10m(1))  ! Surface field - always scalar
      allocate(ustar(1))  ! Surface field - always scalar
      allocate(v10m(1))  ! Surface field - always scalar
      allocate(z0(1))  ! Surface field - always scalar
      allocate(species_mw_g(n_species))
      allocate(species_radius(n_species))
      allocate(species_density(n_species))
      allocate(species_dd_hstar(n_species))
      allocate(species_dd_dvzaersnow(n_species))
      allocate(species_dd_dvzminval_snow(n_species))
      allocate(species_dd_dvzminval_land(n_species))
      allocate(species_lower_radius(n_species))
      allocate(species_upper_radius(n_species))
      allocate(species_is_dust(n_species))
      allocate(species_is_seasalt(n_species))
      species_tendencies = 0.0_fp

      ! Get meteorological data pointer from virtual column (VirtualMet pattern)
      met => column%get_met()

      ! Now allocate categorical fields using the met pointer dimensions
      allocate(frlanduse(size(met%FRLANDUSE)))  ! Categorical field - get size from met pointer
      allocate(iland(size(met%ILAND)))  ! Categorical field - get size from met pointer

      ! Extract required fields from met pointer based on field type and processing mode
      bxheight(1:n_levels) = met%BXHEIGHT(1:n_levels)  ! Atmospheric field - always n_levels
      frlanduse(:) = met%FRLANDUSE(:)  ! Categorical field - full dimension
      iland(:) = met%ILAND(:)  ! Categorical field - full dimension
      isice(1) = met%IsIce  ! Surface field - scalar access
      issnow(1) = met%IsSnow  ! Surface field - scalar access
      lucname(1) = met%LUCNAME  ! Surface field - scalar access
      obk(1) = met%OBK  ! Surface field - scalar access
      ps(1) = met%PS  ! Surface field - scalar access
      rh(1:n_levels) = met%RH(1:n_levels)  ! Atmospheric field - always n_levels
      ts(1) = met%TS  ! Surface field - scalar access
      tstep(1) = this%get_timestep()  ! Special timestep field - retrieved from ProcessInterface
      u10m(1) = met%U10M  ! Surface field - scalar access
      ustar(1) = met%USTAR  ! Surface field - scalar access
      v10m(1) = met%V10M  ! Surface field - scalar access
      z0(1) = met%Z0  ! Surface field - scalar access

      ! Get species concentrations from virtual column
      ! Surface-only processing - get surface level concentrations
      do i = 1, n_species
         species_conc(1, i) = column%get_chem_field(species_indices(i), 1)
      end do

      ! Get species properties from configuration (pre-loaded during initialization)
      ! Use species properties from process configuration
      species_mw_g(1:n_species) = this%process_config%drydep_config%species_mw_g(1:n_species)
      ! Use species properties from process configuration
      species_radius(1:n_species) = this%process_config%drydep_config%species_radius(1:n_species)
      ! Use species properties from process configuration
      species_density(1:n_species) = this%process_config%drydep_config%species_density(1:n_species)
      ! Use species properties from process configuration
      species_dd_hstar(1:n_species) = this%process_config%drydep_config%species_dd_hstar(1:n_species)
      ! Use species properties from process configuration
      species_dd_dvzaersnow(1:n_species) = this%process_config%drydep_config%species_dd_DvzAerSnow(1:n_species)
      ! Use species properties from process configuration
      species_dd_dvzminval_snow(1:n_species) = this%process_config%drydep_config%species_dd_DvzMinVal_snow(1:n_species)
      ! Use species properties from process configuration
      species_dd_dvzminval_land(1:n_species) = this%process_config%drydep_config%species_dd_DvzMinVal_land(1:n_species)
      ! Use species properties from process configuration
      species_lower_radius(1:n_species) = this%process_config%drydep_config%species_lower_radius(1:n_species)
      ! Use species properties from process configuration
      species_upper_radius(1:n_species) = this%process_config%drydep_config%species_upper_radius(1:n_species)
      ! Use species properties from process configuration
      species_is_dust(1:n_species) = this%process_config%drydep_config%species_is_dust(1:n_species)
      ! Use species properties from process configuration
      species_is_seasalt(1:n_species) = this%process_config%drydep_config%species_is_seasalt(1:n_species)

      ! Call the science scheme with optional diagnostic parameters
      ! Note: zhang uses the following diagnostic fields (if diagnostics enabled):
      ! - drydep_con_per_species (Dry deposition concentration per species)
      ! - drydep_velocity_per_species (Dry deposition velocity)
      if (this%process_config%drydep_config%diagnostics) then
         ! Call with diagnostic outputs enabled
         call compute_zhang( &
            n_levels, &
            n_species, &
            this%process_config%zhang_config, &
            bxheight, &
            frlanduse, &
            iland, &
            isice(1), &
            issnow(1), &
            lucname(1), &
            obk(1), &
            ps(1), &
            rh, &
            ts(1), &
            tstep(1), &
            u10m(1), &
            ustar(1), &
            v10m(1), &
            z0(1)            , &
            species_mw_g, &
            species_radius, &
            species_density, &
            this%process_config%drydep_config%species_names, &
            species_dd_hstar, &
            species_dd_dvzaersnow, &
            species_dd_dvzminval_snow, &
            species_dd_dvzminval_land, &
            species_lower_radius, &
            species_upper_radius, &
            species_is_dust, &
            species_is_seasalt, &
            species_conc, &
            species_tendencies, &
            this%process_config%drydep_config%is_gas, &
            this%column_drydep_con_per_species, &
            this%column_drydep_velocity_per_species, &
            this%process_config%drydep_config%diagnostic_species_id         )
      else
         ! Call without diagnostic outputs (optional parameters not passed)
         call compute_zhang( &
            n_levels, &
            n_species, &
            this%process_config%zhang_config, &
            bxheight, &
            frlanduse, &
            iland, &
            isice(1), &
            issnow(1), &
            lucname(1), &
            obk(1), &
            ps(1), &
            rh, &
            ts(1), &
            tstep(1), &
            u10m(1), &
            ustar(1), &
            v10m(1), &
            z0(1)            , &
            species_mw_g, &
            species_radius, &
            species_density, &
            this%process_config%drydep_config%species_names, &
            species_dd_hstar, &
            species_dd_dvzaersnow, &
            species_dd_dvzminval_snow, &
            species_dd_dvzminval_land, &
            species_lower_radius, &
            species_upper_radius, &
            species_is_dust, &
            species_is_seasalt, &
            species_conc, &
            species_tendencies, &
            this%process_config%drydep_config%is_gas &
            )
      end if

      ! Apply tendencies back to virtual column based on tendency_mode
      ! Surface-only processing - apply tendencies to surface level only
      do i = 1, n_species
         ! Skip species that don't match scheme type (gas vs aerosol)
         if (this%process_config%drydep_config%is_gas(i)) cycle
         ! Multiplicative tendency: new_conc = conc * (1 - loss_fraction)
         ! where loss_fraction = 1 - exp(-tendency * dt)
         loss_fraction = max(1.0_fp - exp(-species_tendencies(1, i) * this%get_timestep()), 0.0_fp)
         call column%set_chem_field(1, species_indices(i), &
            species_conc(1, i) * (1.0_fp - loss_fraction))
      end do

   end subroutine run_zhang_scheme_column



   function get_required_met_fields(this) result(field_names)
      class(ProcessDryDepInterface), intent(in) :: this
      character(len=32), allocatable :: field_names(:)
      character(len=32), allocatable :: scheme_fields(:)
      character(len=32), allocatable :: process_fields(:)
      character(len=32), allocatable :: gas_scheme_fields(:), aero_scheme_fields(:)
      integer :: gas_scheme_count, aero_scheme_count
      character(len=32), allocatable :: unique_fields(:)
      integer :: total_fields, scheme_count, process_count, i, j, unique_count
      logical :: is_duplicate

      ! No process-level required fields
      process_count = 0
      allocate(process_fields(0))

      ! For gas/aero differentiated processes, get fields from both schemes
      ! Get gas scheme fields
      select case (trim(this%process_config%drydep_config%gas_scheme))
       case ('wesely')
         gas_scheme_count = 21
         allocate(gas_scheme_fields(gas_scheme_count))
         gas_scheme_fields(1) = 'USTAR'
         gas_scheme_fields(2) = 'TSTEP'
         gas_scheme_fields(3) = 'TS'
         gas_scheme_fields(4) = 'SWGDN'
         gas_scheme_fields(5) = 'SUNCOSmid'
         gas_scheme_fields(6) = 'OBK'
         gas_scheme_fields(7) = 'CLDFRC'
         gas_scheme_fields(8) = 'BXHEIGHT'
         gas_scheme_fields(9) = 'Z0'
         gas_scheme_fields(10) = 'PS'
         gas_scheme_fields(11) = 'FRLAI'
         gas_scheme_fields(12) = 'ILAND'
         gas_scheme_fields(13) = 'SALINITY'
         gas_scheme_fields(14) = 'FRLANDUSE'
         gas_scheme_fields(15) = 'TSKIN'
         gas_scheme_fields(16) = 'LON'
         gas_scheme_fields(17) = 'LAT'
         gas_scheme_fields(18) = 'LUCNAME'
         gas_scheme_fields(19) = 'IsSnow'
         gas_scheme_fields(20) = 'IsIce'
         gas_scheme_fields(21) = 'IsLand'
       case default
         gas_scheme_count = 0
         allocate(gas_scheme_fields(0))
      end select

      ! Get aerosol scheme fields
      select case (trim(this%process_config%drydep_config%aero_scheme))
       case ('gocart')
         aero_scheme_count = 13
         allocate(aero_scheme_fields(aero_scheme_count))
         aero_scheme_fields(1) = 'USTAR'
         aero_scheme_fields(2) = 'TSTEP'
         aero_scheme_fields(3) = 'T'
         aero_scheme_fields(4) = 'AIRDEN'
         aero_scheme_fields(5) = 'Z'
         aero_scheme_fields(6) = 'LWI'
         aero_scheme_fields(7) = 'PBLH'
         aero_scheme_fields(8) = 'HFLUX'
         aero_scheme_fields(9) = 'Z0H'
         aero_scheme_fields(10) = 'U10M'
         aero_scheme_fields(11) = 'V10M'
         aero_scheme_fields(12) = 'FRLAKE'
         aero_scheme_fields(13) = 'GWETTOP'
       case ('zhang')
         aero_scheme_count = 15
         allocate(aero_scheme_fields(aero_scheme_count))
         aero_scheme_fields(1) = 'USTAR'
         aero_scheme_fields(2) = 'TSTEP'
         aero_scheme_fields(3) = 'TS'
         aero_scheme_fields(4) = 'OBK'
         aero_scheme_fields(5) = 'BXHEIGHT'
         aero_scheme_fields(6) = 'Z0'
         aero_scheme_fields(7) = 'RH'
         aero_scheme_fields(8) = 'PS'
         aero_scheme_fields(9) = 'U10M'
         aero_scheme_fields(10) = 'V10M'
         aero_scheme_fields(11) = 'FRLANDUSE'
         aero_scheme_fields(12) = 'ILAND'
         aero_scheme_fields(13) = 'LUCNAME'
         aero_scheme_fields(14) = 'IsSnow'
         aero_scheme_fields(15) = 'IsIce'
       case default
         aero_scheme_count = 0
         allocate(aero_scheme_fields(0))
      end select

      ! Combine all fields (process + gas_scheme + aero_scheme) and remove duplicates
      ! First estimate maximum possible size (without duplicates)
      total_fields = process_count + gas_scheme_count + aero_scheme_count
      allocate(unique_fields(total_fields))
      unique_count = 0

      ! Add process-level fields first
      do i = 1, process_count
         unique_count = unique_count + 1
         unique_fields(unique_count) = process_fields(i)
      end do

      ! Add gas scheme fields (check for duplicates)
      do i = 1, gas_scheme_count
         is_duplicate = .false.
         do j = 1, unique_count
            if (trim(gas_scheme_fields(i)) == trim(unique_fields(j))) then
               is_duplicate = .true.
               exit
            end if
         end do
         if (.not. is_duplicate) then
            unique_count = unique_count + 1
            unique_fields(unique_count) = gas_scheme_fields(i)
         end if
      end do

      ! Add aerosol scheme fields (check for duplicates)
      do i = 1, aero_scheme_count
         is_duplicate = .false.
         do j = 1, unique_count
            if (trim(aero_scheme_fields(i)) == trim(unique_fields(j))) then
               is_duplicate = .true.
               exit
            end if
         end do
         if (.not. is_duplicate) then
            unique_count = unique_count + 1
            unique_fields(unique_count) = aero_scheme_fields(i)
         end if
      end do

      ! Allocate final result array with exact size
      allocate(field_names(unique_count))
      field_names(1:unique_count) = unique_fields(1:unique_count)

      ! Clean up temporary arrays
      if (allocated(unique_fields)) deallocate(unique_fields)
      if (allocated(process_fields)) deallocate(process_fields)
      if (allocated(gas_scheme_fields)) deallocate(gas_scheme_fields)
      if (allocated(aero_scheme_fields)) deallocate(aero_scheme_fields)

   end function get_required_met_fields

   function get_required_diagnostic_fields(this) result(field_names)
      class(ProcessDryDepInterface), intent(in) :: this
      character(len=64), allocatable :: field_names(:)

      allocate(field_names(2))
      field_names(1) = 'drydep_con_per_species'
      field_names(2) = 'drydep_velocity_per_species'

   end function get_required_diagnostic_fields


   subroutine register_and_allocate_diagnostics(this, container, rc)
      use diagnosticinterface_mod, only: diagnosticregistrytype, diag_real_2d, diag_real_3d

      class(ProcessDryDepInterface), intent(inout) :: this
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
      if (.not. this%process_config%drydep_config%diagnostics) then
         return
      endif

      ! Get managers
      diag_mgr => container%get_diagnostic_manager()
      grid_mgr => container%get_grid_manager()

      ! Register this process with diagnostic manager (only once per process)
      call diag_mgr%register_process('drydep', rc)
      if (rc /= cc_success) return

      ! Get the process registry for registering individual diagnostics
      call diag_mgr%get_process_registry('drydep', registry, rc)
      if (rc /= cc_success) return

      ! Get grid dimensions
      call grid_mgr%get_shape(nx, ny, nz)
      dims_2d = [nx, ny]

      ! Get species count for 3D diagnostics
      n_species = this%process_config%drydep_config%n_species
      dims_3d_species = [nx, ny, n_species]

      ! Register drydep_con_per_species
      ! Register individual 2D fields for each diagnostic species (species-only diagnostics)
      if (this%process_config%drydep_config%n_diagnostic_species > 0) then
         do i = 1, this%process_config%drydep_config%n_diagnostic_species
            write(field_name, '(A,A,A)') 'drydep_con_', &
               trim(this%process_config%drydep_config%diagnostic_species(i))
            call this%register_diagnostic_field(registry, trim(field_name), &
               'Dry deposition concentration per species', &
               'ug/kg or ppm', diag_real_2d, &
               'drydep', dims_2d, rc=rc)
            if (rc /= cc_success) return
         end do
      end if
      if (rc /= cc_success) return

      ! Register drydep_velocity_per_species
      ! Register individual 2D fields for each diagnostic species (species-only diagnostics)
      if (this%process_config%drydep_config%n_diagnostic_species > 0) then
         do i = 1, this%process_config%drydep_config%n_diagnostic_species
            write(field_name, '(A,A,A)') 'drydep_velocity_', &
               trim(this%process_config%drydep_config%diagnostic_species(i))
            call this%register_diagnostic_field(registry, trim(field_name), &
               'Dry deposition velocity', &
               'm/s', diag_real_2d, &
               'drydep', dims_2d, rc=rc)
            if (rc /= cc_success) return
         end do
      end if
      if (rc /= cc_success) return

      ! Get selected scheme(s)
      ! For gas/aero differentiated processes, register diagnostics from both schemes
      ! Track registered diagnostics to avoid duplicates
      ! Register gas scheme diagnostics
      select case (trim(this%process_config%drydep_config%gas_scheme))
       case ('wesely')
         ! Register wesely-specific diagnostics (gas)
       case default
         ! Unknown gas scheme
      end select

      ! Register aerosol scheme diagnostics (only if not already registered)
      select case (trim(this%process_config%drydep_config%aero_scheme))
       case ('gocart')
         ! Register gocart-specific diagnostics (aerosol)
       case ('zhang')
         ! Register zhang-specific diagnostics (aerosol)
       case default
         ! Unknown aerosol scheme
      end select

      ! Now allocate diagnostic class members after successful registration
      ! First, deallocate if already allocated (for scheme switching)
      if (allocated(this%column_drydep_con_per_species)) deallocate(this%column_drydep_con_per_species)
      if (allocated(this%column_drydep_velocity_per_species)) deallocate(this%column_drydep_velocity_per_species)

      ! Allocate and initialize scheme-specific diagnostic fields based on selected scheme
      ! For gas/aero differentiated process, allocate diagnostics from both gas and aero schemes
      ! Track allocated diagnostics to avoid duplicates

      ! Allocate common diagnostic fields (used by all schemes)
      ! 1D diagnostic: diagnostic species only - allocated based on n_diagnostic_species
      if (this%process_config%drydep_config%n_diagnostic_species > 0) then
         allocate(this%column_drydep_con_per_species(this%process_config%drydep_config%n_diagnostic_species))
      end if
      if (allocated(this%column_drydep_con_per_species)) this%column_drydep_con_per_species = 0.0_fp
      ! 1D diagnostic: diagnostic species only - allocated based on n_diagnostic_species
      if (this%process_config%drydep_config%n_diagnostic_species > 0) then
         allocate(this%column_drydep_velocity_per_species(this%process_config%drydep_config%n_diagnostic_species))
      end if
      if (allocated(this%column_drydep_velocity_per_species)) this%column_drydep_velocity_per_species = 0.0_fp

      ! Gas scheme diagnostics
      select case (trim(this%process_config%drydep_config%gas_scheme))
       case ("wesely")
         ! Gas scheme-specific diagnostics for wesely
       case default
         ! No gas scheme-specific diagnostics for unknown schemes
      end select

      ! Aerosol scheme diagnostics (only allocate if not already allocated by gas scheme)
      select case (trim(this%process_config%drydep_config%aero_scheme))
       case ("gocart")
         ! Aerosol scheme-specific diagnostics for gocart
       case ("zhang")
         ! Aerosol scheme-specific diagnostics for zhang
       case default
         ! No aerosol scheme-specific diagnostics for unknown schemes
      end select

   end subroutine register_and_allocate_diagnostics

   subroutine calculate_and_update_diagnostics(this, column, container, rc)
      class(ProcessDryDepInterface), intent(inout) :: this
      type(VirtualColumnType), intent(in) :: column
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: i_col, j_col  ! Column grid position
      integer :: i  ! Loop variable for diagnostic species
      character(len=256) :: field_name  ! For constructing species-specific field names

      rc = cc_success

      ! Skip if diagnostics not enabled
      if (.not. this%process_config%drydep_config%diagnostics) return

      ! Get column grid position (x, y indices)
      call column%get_position(i_col, j_col)

      ! Update common diagnostic fields (used by all schemes)
      ! Update individual species diagnostic fields (species-only diagnostics)
      if (this%process_config%drydep_config%n_diagnostic_species > 0) then
         do i = 1, this%process_config%drydep_config%n_diagnostic_species
            write(field_name, '(A,A,A)') 'drydep_con_', &
               trim(this%process_config%drydep_config%diagnostic_species(i))
            call this%update_scalar_diagnostic_column(trim(field_name), &
               this%column_drydep_con_per_species(i), &
               i_col, j_col, container, rc)
            if (rc /= cc_success) return
         end do
      end if
      ! Update individual species diagnostic fields (species-only diagnostics)
      if (this%process_config%drydep_config%n_diagnostic_species > 0) then
         do i = 1, this%process_config%drydep_config%n_diagnostic_species
            write(field_name, '(A,A,A)') 'drydep_velocity_', &
               trim(this%process_config%drydep_config%diagnostic_species(i))
            call this%update_scalar_diagnostic_column(trim(field_name), &
               this%column_drydep_velocity_per_species(i), &
               i_col, j_col, container, rc)
            if (rc /= cc_success) return
         end do
      end if
      ! Update scheme-specific diagnostic fields based on active scheme
      ! For gas/aero differentiated process, update diagnostics from both gas and aero schemes
      ! Track updated diagnostics to avoid duplicates
      ! Update gas scheme diagnostics
      select case (trim(this%process_config%drydep_config%gas_scheme))
       case ("wesely")
         ! Gas scheme-specific diagnostics for wesely
       case default
         ! No gas scheme diagnostics for unknown schemes
      end select

      ! Update aerosol scheme diagnostics (only if not already updated)
      select case (trim(this%process_config%drydep_config%aero_scheme))
       case ("gocart")
         ! Aerosol scheme-specific diagnostics for gocart
       case ("zhang")
         ! Aerosol scheme-specific diagnostics for zhang
       case default
         ! No aerosol scheme diagnostics for unknown schemes
      end select

   end subroutine calculate_and_update_diagnostics


   subroutine set_drydep_scheme(this, scheme_name, gas_scheme)
      class(ProcessDryDepInterface), intent(inout) :: this
      character(len=*), intent(in) :: scheme_name
      logical, intent(in), optional :: gas_scheme

      logical :: is_gas

      ! Default to setting the single scheme for backward compatibility
      is_gas = .true.
      if (present(gas_scheme)) is_gas = gas_scheme

      if (is_gas) then
         this%process_config%drydep_config%gas_scheme = trim(scheme_name)
      else
         this%process_config%drydep_config%aero_scheme = trim(scheme_name)
      end if

   end subroutine set_drydep_scheme

   function get_drydep_scheme(this, gas_scheme) result(scheme_name)
      class(ProcessDryDepInterface), intent(in) :: this
      logical, intent(in), optional :: gas_scheme
      character(len=64) :: scheme_name

      logical :: is_gas

      ! Default to getting the gas scheme for backward compatibility
      is_gas = .true.
      if (present(gas_scheme)) is_gas = gas_scheme

      if (is_gas) then
         scheme_name = trim(this%process_config%drydep_config%gas_scheme)
      else
         scheme_name = trim(this%process_config%drydep_config%aero_scheme)
      end if

   end function get_drydep_scheme

end module processdrydepinterface_mod
```


