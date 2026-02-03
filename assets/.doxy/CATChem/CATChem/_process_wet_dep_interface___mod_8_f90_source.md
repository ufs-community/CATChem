

# File ProcessWetDepInterface\_Mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**wetdep**](dir_8b9a0ce556ea4a65f6920dfb49dcd69d.md) **>** [**ProcessWetDepInterface\_Mod.F90**](_process_wet_dep_interface___mod_8_f90.md)

[Go to the documentation of this file](_process_wet_dep_interface___mod_8_f90.md)


```Fortran


module processwetdepinterface_mod

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
   use wetdepcommon_mod, only: wetdepprocessconfig

   ! Scheme modules
   use wetdepscheme_jacob_mod, only: compute_jacob

   implicit none
   private

   public :: processwetdepinterface

   type, extends(columnprocessinterface) :: processwetdepinterface
      private

      ! Unified process configuration (bridges ConfigManager to process-specific config)
      type(WetDepProcessConfig), public :: process_config

      ! Process utilities (leverage core infrastructure)
      type(ChemStateType), pointer :: chem_state => null()
      type(MetStateType), pointer :: met_state => null()
      ! Note: state_manager pointer removed as it was never used

      ! Process-specific diagnostic indices (base class handles storage)
      integer :: diag_wetdep_mass_per_species_per_level_idx = -1
      integer :: diag_wetdep_flux_per_species_per_level_idx = -1

      ! Column-level diagnostic storage for interfacing with DiagManager
      ! These are allocated per-column during processing and can be used to
      ! accumulate data for DiagManager updates
      real(fp), allocatable :: column_wetdep_mass_per_species_per_level(:,:)         ! 2D: levels x species - per column
      real(fp), allocatable :: column_wetdep_flux_per_species_per_level(:,:)         ! 2D: levels x species - per column

      ! Scheme-specific diagnostic storage (shared across all schemes that use them)

   contains
      ! Required ProcessInterface implementations
      procedure :: init => process_init
      procedure :: run => process_run
      procedure :: finalize => process_finalize
      procedure :: parse_process_config => parse_wetdep_config

      ! Required ColumnProcessInterface implementations
      procedure :: init_column_processing => init_column_processing
      procedure :: run_column => run_column
      procedure :: finalize_column_processing => finalize_column_processing

      ! ProcessInterface capability registration
      procedure :: get_required_met_fields => get_required_met_fields
      procedure :: get_required_diagnostic_fields => get_required_diagnostic_fields

      ! Public testing interface for scheme manipulation
      procedure :: set_scheme => set_wetdep_scheme
      procedure :: get_scheme => get_wetdep_scheme

      ! Process-specific implementations (column virtualization)
      procedure, private :: run_active_scheme_column
      procedure, private :: run_jacob_scheme_column

      ! Diagnostic procedures (override base class method)
      procedure :: register_diagnostics => register_and_allocate_diagnostics
      procedure, private :: register_and_allocate_diagnostics
      procedure, private :: calculate_and_update_diagnostics

   end type processwetdepinterface

contains

   subroutine process_init(this, container, rc)
      class(ProcessWetDepInterface), intent(inout) :: this
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
      this%name = 'wetdep'
      this%version = '1.0.0'
      this%description = 'Process for computing wet deposition of gas and aerosol species'

      ! Parse process-specific configuration using unified approach
      call this%parse_process_config(container, error_manager, rc)
      if (rc /= cc_success) then
         return
      end if

      ! Get state pointers from container (needed for species loading)
      this%chem_state => container%get_chem_state_ptr()
      this%met_state => container%get_met_state_ptr()

      ! Load species from ChemState based on is_wetdep property
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
      class(ProcessWetDepInterface), intent(inout) :: this
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

      ! Currently no global 3D operations needed for wetdep process
      ! All processing happens in run_column() method

   end subroutine process_run

   subroutine process_finalize(this, rc)
      class(ProcessWetDepInterface), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success

      ! Finalize column processing
      call this%finalize_column_processing(rc)
      if (rc /= cc_success) return

      ! Deallocate diagnostic class members
      if (allocated(this%column_wetdep_mass_per_species_per_level)) deallocate(this%column_wetdep_mass_per_species_per_level)
      if (allocated(this%column_wetdep_flux_per_species_per_level)) deallocate(this%column_wetdep_flux_per_species_per_level)
      ! Deallocate scheme-specific diagnostic fields (only deallocate unique fields once)

      ! Finalize unified configuration
      call this%process_config%finalize()

   end subroutine process_finalize

   subroutine parse_wetdep_config(this, state_manager, error_manager, rc)
      class(ProcessWetDepInterface), intent(inout) :: this
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

      ! Use the unified configuration loader from WetDepCommon_Mod
      ! This handles the complexity of parsing hierarchical YAML into process-specific types
      call this%process_config%load_from_config(config_manager, error_manager)
      ! Note: Error handling managed by error_manager internally

      ! Process is now configured - the unified config contains all scheme-specific settings

   end subroutine parse_wetdep_config

   !========================================================================
   ! Column Processing Interface Implementation
   !========================================================================

   subroutine init_column_processing(this, container, rc)
      class(ProcessWetDepInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      rc = cc_success

      ! Enable column processing and set batch size for optimal performance
      call this%enable_column_processing()
      call this%set_column_batch_size(50)  ! Process 50 columns at a time

      ! Any process-specific column processing setup would go here
      ! For wetdep, no additional setup is needed

   end subroutine init_column_processing

   subroutine run_column(this, column, container, rc)
      class(ProcessWetDepInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      rc = cc_success

      ! Check if process is active
      if (.not. this%process_config%is_active) return

      ! Delegate to the active scheme for column processing
      call this%run_active_scheme_column(column, rc)

      ! Calculate and update diagnostics if enabled
      if (this%process_config%wetdep_config%diagnostics .and. rc == cc_success) then
         call this%calculate_and_update_diagnostics(column, container, rc)
      end if

   end subroutine run_column

   subroutine finalize_column_processing(this, rc)
      class(ProcessWetDepInterface), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success

      ! Clean up column processing
      call this%disable_column_processing()

      ! Any process-specific cleanup would go here

   end subroutine finalize_column_processing

   subroutine run_active_scheme_column(this, column, rc)
      class(ProcessWetDepInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      integer, intent(out) :: rc

      rc = cc_success

      ! Delegate to appropriate scheme using unified config
      select case (trim(this%process_config%wetdep_config%scheme))
       case ('jacob')
         call this%run_jacob_scheme_column(column, rc)
       case default
         rc = cc_failure
      end select

   end subroutine run_active_scheme_column

   subroutine run_jacob_scheme_column(this, column, rc)
      class(ProcessWetDepInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      integer, intent(out) :: rc

      ! Local variables for scheme calculation
      type(VirtualMetType), pointer :: met => null()  ! Pointer to meteorological data
      ! Meteorological fields
      real(fp), allocatable :: airden_dry(:)
      real(fp), allocatable :: mairden(:)
      real(fp), allocatable :: pedge(:)
      real(fp), allocatable :: pfilsan(:)
      real(fp), allocatable :: pfllsan(:)
      real(fp), allocatable :: reevapls(:)
      real(fp), allocatable :: t(:)
      real(fp), allocatable :: tstep(:)
      ! Species properties
      logical, allocatable :: species_is_aerosol(:)
      real(fp), allocatable :: species_henry_cr(:)
      real(fp), allocatable :: species_henry_k0(:)
      real(fp), allocatable :: species_henry_pKa(:)
      real(fp), allocatable :: species_wd_retfactor(:)
      logical, allocatable :: species_wd_LiqAndGas(:)
      real(fp), allocatable :: species_wd_convfacI2G(:)
      real(fp), allocatable :: species_wd_rainouteff(:,:)
      real(fp), allocatable :: species_radius(:)
      real(fp), allocatable :: species_mw_g(:)
      real(fp), allocatable :: species_conc(:,:)
      real(fp), allocatable :: species_tendencies(:,:)
      integer :: n_species, n_levels, n_chem, n_emis, i, k
      integer, allocatable :: species_indices(:)

      rc = cc_success

      ! Get dimensions from virtual column
      call column%get_dimensions(n_levels, n_chem, n_emis)  ! Full column processing

      ! Get wetdep species information from process configuration
      n_species = this%process_config%wetdep_config%n_species
      if (n_species <= 0) then
         return
      end if

      ! Get species indices directly from configuration (pre-computed)
      allocate(species_indices(n_species))
      species_indices(1:n_species) = this%process_config%wetdep_config%species_indices(1:n_species)

      ! Allocate arrays
      allocate(species_conc(n_levels, n_species))
      allocate(species_tendencies(n_levels, n_species))
      ! Allocate meteorological field arrays based on field type and process configuration
      allocate(airden_dry(n_levels))  ! Atmospheric field - always n_levels
      allocate(mairden(n_levels))  ! Atmospheric field - always n_levels
      allocate(pedge(n_levels+1))  ! Edge field - always n_levels+1
      allocate(pfilsan(n_levels+1))  ! Edge field - always n_levels+1
      allocate(pfllsan(n_levels+1))  ! Edge field - always n_levels+1
      allocate(reevapls(n_levels))  ! Atmospheric field - always n_levels
      allocate(t(n_levels))  ! Atmospheric field - always n_levels
      allocate(tstep(1))  ! Special timestep field - scalar
      allocate(species_is_aerosol(n_species))
      allocate(species_henry_cr(n_species))
      allocate(species_henry_k0(n_species))
      allocate(species_henry_pka(n_species))
      allocate(species_wd_retfactor(n_species))
      allocate(species_wd_liqandgas(n_species))
      allocate(species_wd_convfaci2g(n_species))
      allocate(species_wd_rainouteff(n_species, 3))
      allocate(species_radius(n_species))
      allocate(species_mw_g(n_species))
      species_tendencies = 0.0_fp

      ! Get meteorological data pointer from virtual column (VirtualMet pattern)
      met => column%get_met()

      ! Now allocate categorical fields using the met pointer dimensions

      ! Extract required fields from met pointer based on field type and processing mode
      airden_dry(1:n_levels) = met%AIRDEN_DRY(1:n_levels)  ! Atmospheric field - always n_levels
      mairden(1:n_levels) = met%MAIRDEN(1:n_levels)  ! Atmospheric field - always n_levels
      pedge(1:n_levels+1) = met%PEDGE(1:n_levels+1)  ! Edge field - always n_levels+1
      pfilsan(1:n_levels+1) = met%PFILSAN(1:n_levels+1)  ! Edge field - always n_levels+1
      pfllsan(1:n_levels+1) = met%PFLLSAN(1:n_levels+1)  ! Edge field - always n_levels+1
      reevapls(1:n_levels) = met%REEVAPLS(1:n_levels)  ! Atmospheric field - always n_levels
      t(1:n_levels) = met%T(1:n_levels)  ! Atmospheric field - always n_levels
      tstep(1) = this%get_timestep()  ! Special timestep field - retrieved from ProcessInterface

      ! Get species concentrations from virtual column
      ! Full column processing - get concentrations for all levels
      do k = 1, n_levels
         do i = 1, n_species
            species_conc(k, i) = column%get_chem_field(species_indices(i), k)
         end do
      end do

      ! Get species properties from configuration (pre-loaded during initialization)
      ! Use species properties from process configuration
      species_is_aerosol(1:n_species) = this%process_config%wetdep_config%species_is_aerosol(1:n_species)
      ! Use species properties from process configuration
      species_henry_cr(1:n_species) = this%process_config%wetdep_config%species_henry_cr(1:n_species)
      ! Use species properties from process configuration
      species_henry_k0(1:n_species) = this%process_config%wetdep_config%species_henry_k0(1:n_species)
      ! Use species properties from process configuration
      species_henry_pka(1:n_species) = this%process_config%wetdep_config%species_henry_pKa(1:n_species)
      ! Use species properties from process configuration
      species_wd_retfactor(1:n_species) = this%process_config%wetdep_config%species_wd_retfactor(1:n_species)
      ! Use species properties from process configuration
      species_wd_liqandgas(1:n_species) = this%process_config%wetdep_config%species_wd_LiqAndGas(1:n_species)
      ! Use species properties from process configuration
      species_wd_convfaci2g(1:n_species) = this%process_config%wetdep_config%species_wd_convfacI2G(1:n_species)
      ! Use species properties from process configuration
      species_wd_rainouteff(1:n_species, :) = this%process_config%wetdep_config%species_wd_rainouteff(1:n_species, :)
      ! Use species properties from process configuration
      species_radius(1:n_species) = this%process_config%wetdep_config%species_radius(1:n_species)
      ! Use species properties from process configuration
      species_mw_g(1:n_species) = this%process_config%wetdep_config%species_mw_g(1:n_species)

      ! Call the science scheme with optional diagnostic parameters
      ! Note: jacob uses the following diagnostic fields (if diagnostics enabled):
      ! - wetdep_mass_per_species_per_level (Wet deposition mass loss per species per level)
      ! - wetdep_flux_per_species_per_level (Wet deposition flux per species per level)
      if (this%process_config%wetdep_config%diagnostics) then
         ! Call with diagnostic outputs enabled
         call compute_jacob( &
            n_levels, &
            n_species, &
            this%process_config%jacob_config, &
            airden_dry, &
            mairden, &
            pedge, &
            pfilsan, &
            pfllsan, &
            reevapls, &
            t, &
            tstep(1)            , &
            species_is_aerosol, &
            this%process_config%wetdep_config%species_names, &
            species_henry_cr, &
            species_henry_k0, &
            species_henry_pka, &
            species_wd_retfactor, &
            species_wd_liqandgas, &
            species_wd_convfaci2g, &
            species_wd_rainouteff, &
            species_radius, &
            species_mw_g, &
            species_conc, &
            species_tendencies, &
            this%column_wetdep_mass_per_species_per_level, &
            this%column_wetdep_flux_per_species_per_level, &
            this%process_config%wetdep_config%diagnostic_species_id         )
      else
         ! Call without diagnostic outputs (optional parameters not passed)
         call compute_jacob( &
            n_levels, &
            n_species, &
            this%process_config%jacob_config, &
            airden_dry, &
            mairden, &
            pedge, &
            pfilsan, &
            pfllsan, &
            reevapls, &
            t, &
            tstep(1)            , &
            species_is_aerosol, &
            this%process_config%wetdep_config%species_names, &
            species_henry_cr, &
            species_henry_k0, &
            species_henry_pka, &
            species_wd_retfactor, &
            species_wd_liqandgas, &
            species_wd_convfaci2g, &
            species_wd_rainouteff, &
            species_radius, &
            species_mw_g, &
            species_conc, &
            species_tendencies &
            )
      end if

      ! Apply tendencies back to virtual column based on tendency_mode
      ! Full column processing - apply tendencies to all levels
      do k = 1, n_levels
         do i = 1, n_species
            ! Replacement tendency: new_conc = tendency (tendency is the new value)
            call column%set_chem_field(k, species_indices(i), &
               species_tendencies(k, i))
         end do
      end do

   end subroutine run_jacob_scheme_column



   function get_required_met_fields(this) result(field_names)
      class(ProcessWetDepInterface), intent(in) :: this
      character(len=32), allocatable :: field_names(:)
      character(len=32), allocatable :: scheme_fields(:)
      character(len=32), allocatable :: process_fields(:)
      character(len=32), allocatable :: unique_fields(:)
      integer :: total_fields, scheme_count, process_count, i, j, unique_count
      logical :: is_duplicate

      ! No process-level required fields
      process_count = 0
      allocate(process_fields(0))

      ! Get scheme-specific fields based on selected scheme
      select case (trim(this%process_config%wetdep_config%scheme))
       case ('jacob')
         scheme_count = 8
         allocate(scheme_fields(scheme_count))
         scheme_fields(1) = 'T'
         scheme_fields(2) = 'TSTEP'
         scheme_fields(3) = 'AIRDEN_DRY'
         scheme_fields(4) = 'MAIRDEN'
         scheme_fields(5) = 'PFLLSAN'
         scheme_fields(6) = 'PFILSAN'
         scheme_fields(7) = 'PEDGE'
         scheme_fields(8) = 'REEVAPLS'
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
      class(ProcessWetDepInterface), intent(in) :: this
      character(len=64), allocatable :: field_names(:)

      allocate(field_names(2))
      field_names(1) = 'wetdep_mass_per_species_per_level'
      field_names(2) = 'wetdep_flux_per_species_per_level'

   end function get_required_diagnostic_fields


   subroutine register_and_allocate_diagnostics(this, container, rc)
      use diagnosticinterface_mod, only: diagnosticregistrytype, diag_real_2d, diag_real_3d

      class(ProcessWetDepInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(DiagnosticManagerType), pointer :: diag_mgr
      type(DiagnosticRegistryType), pointer :: registry
      type(GridManagerType), pointer :: grid_mgr
      character(len=256) :: field_name  ! For constructing species-specific field names
      integer :: i  ! Loop variable for diagnostic species
      integer :: nx, ny, nz
      integer :: dims_2d(2)
      integer :: dims_3d_levels(3)

      rc = cc_success

      ! Only register diagnostics if enabled in config
      if (.not. this%process_config%wetdep_config%diagnostics) then
         return
      endif

      ! Get managers
      diag_mgr => container%get_diagnostic_manager()
      grid_mgr => container%get_grid_manager()

      ! Register this process with diagnostic manager (only once per process)
      call diag_mgr%register_process('wetdep', rc)
      if (rc /= cc_success) return

      ! Get the process registry for registering individual diagnostics
      call diag_mgr%get_process_registry('wetdep', registry, rc)
      if (rc /= cc_success) return

      ! Get grid dimensions
      call grid_mgr%get_shape(nx, ny, nz)
      dims_2d = [nx, ny]

      dims_3d_levels = [nx, ny, nz]

      ! Register wetdep_mass_per_species_per_level
      ! Register individual 3D fields for each diagnostic species (level + species diagnostics)
      if (this%process_config%wetdep_config%n_diagnostic_species > 0) then
         do i = 1, this%process_config%wetdep_config%n_diagnostic_species
            write(field_name, '(A,A,A)') 'wetdep_mass_', &
               trim(this%process_config%wetdep_config%diagnostic_species(i))
            call this%register_diagnostic_field(registry, trim(field_name), &
               'Wet deposition mass loss per species per level', &
               'kg/m2', diag_real_3d, &
               'wetdep', dims_3d_levels, rc=rc)
            if (rc /= cc_success) return
         end do
      end if
      if (rc /= cc_success) return

      ! Register wetdep_flux_per_species_per_level
      ! Register individual 3D fields for each diagnostic species (level + species diagnostics)
      if (this%process_config%wetdep_config%n_diagnostic_species > 0) then
         do i = 1, this%process_config%wetdep_config%n_diagnostic_species
            write(field_name, '(A,A,A)') 'wetdep_flux_', &
               trim(this%process_config%wetdep_config%diagnostic_species(i))
            call this%register_diagnostic_field(registry, trim(field_name), &
               'Wet deposition flux per species per level', &
               'kg/m2/s', diag_real_3d, &
               'wetdep', dims_3d_levels, rc=rc)
            if (rc /= cc_success) return
         end do
      end if
      if (rc /= cc_success) return

      ! Get selected scheme(s)
      ! Register scheme-specific diagnostics based on selected scheme
      select case (trim(this%process_config%wetdep_config%scheme))

       case ('jacob')
         ! Register jacob-specific diagnostics
       case default
         ! Unknown scheme - only register common diagnostics
         ! (already done above)

      end select

      ! Now allocate diagnostic class members after successful registration
      ! First, deallocate if already allocated (for scheme switching)
      if (allocated(this%column_wetdep_mass_per_species_per_level)) deallocate(this%column_wetdep_mass_per_species_per_level)
      if (allocated(this%column_wetdep_flux_per_species_per_level)) deallocate(this%column_wetdep_flux_per_species_per_level)

      ! Allocate and initialize scheme-specific diagnostic fields based on selected scheme
      ! For non-gas/aero differentiated process, allocate diagnostics normally

      ! Allocate common diagnostic fields (used by all schemes)
      ! 2D diagnostic: levels x diagnostic_species
      if (nz > 0 .and. this%process_config%wetdep_config%n_diagnostic_species > 0) then
         allocate(this%column_wetdep_mass_per_species_per_level(nz, this%process_config%wetdep_config%n_diagnostic_species))
      end if
      if (allocated(this%column_wetdep_mass_per_species_per_level)) this%column_wetdep_mass_per_species_per_level = 0.0_fp
      ! 2D diagnostic: levels x diagnostic_species
      if (nz > 0 .and. this%process_config%wetdep_config%n_diagnostic_species > 0) then
         allocate(this%column_wetdep_flux_per_species_per_level(nz, this%process_config%wetdep_config%n_diagnostic_species))
      end if
      if (allocated(this%column_wetdep_flux_per_species_per_level)) this%column_wetdep_flux_per_species_per_level = 0.0_fp

      ! Allocate scheme-specific diagnostics
      select case (trim(this%process_config%wetdep_config%scheme))
       case ('jacob')
         ! Scheme-specific diagnostics for jacob
       case default
         ! No scheme-specific diagnostics for unknown schemes
      end select

   end subroutine register_and_allocate_diagnostics

   subroutine calculate_and_update_diagnostics(this, column, container, rc)
      class(ProcessWetDepInterface), intent(inout) :: this
      type(VirtualColumnType), intent(in) :: column
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: i_col, j_col  ! Column grid position
      integer :: i  ! Loop variable for diagnostic species
      character(len=256) :: field_name  ! For constructing species-specific field names

      rc = cc_success

      ! Skip if diagnostics not enabled
      if (.not. this%process_config%wetdep_config%diagnostics) return

      ! Get column grid position (x, y indices)
      call column%get_position(i_col, j_col)

      ! Update common diagnostic fields (used by all schemes)
      ! Update individual 3D fields for each diagnostic species (level + species diagnostics)
      if (this%process_config%wetdep_config%n_diagnostic_species > 0) then
         do i = 1, this%process_config%wetdep_config%n_diagnostic_species
            write(field_name, '(A,A,A)') 'wetdep_mass_', &
               trim(this%process_config%wetdep_config%diagnostic_species(i))
            call this%update_1d_diagnostic_column(trim(field_name), &
               this%column_wetdep_mass_per_species_per_level(:,i), &
               i_col, j_col, container, rc)
            if (rc /= cc_success) return
         end do
      end if
      ! Update individual 3D fields for each diagnostic species (level + species diagnostics)
      if (this%process_config%wetdep_config%n_diagnostic_species > 0) then
         do i = 1, this%process_config%wetdep_config%n_diagnostic_species
            write(field_name, '(A,A,A)') 'wetdep_flux_', &
               trim(this%process_config%wetdep_config%diagnostic_species(i))
            call this%update_1d_diagnostic_column(trim(field_name), &
               this%column_wetdep_flux_per_species_per_level(:,i), &
               i_col, j_col, container, rc)
            if (rc /= cc_success) return
         end do
      end if
      ! Update scheme-specific diagnostic fields based on active scheme
      select case (trim(this%process_config%wetdep_config%scheme))
       case ("jacob")
         ! Scheme-specific diagnostics for jacob
      end select

   end subroutine calculate_and_update_diagnostics


   subroutine set_wetdep_scheme(this, scheme_name)
      class(ProcessWetDepInterface), intent(inout) :: this
      character(len=*), intent(in) :: scheme_name

      this%process_config%wetdep_config%scheme = trim(scheme_name)

   end subroutine set_wetdep_scheme

   function get_wetdep_scheme(this) result(scheme_name)
      class(ProcessWetDepInterface), intent(in) :: this
      character(len=64) :: scheme_name

      scheme_name = trim(this%process_config%wetdep_config%scheme)

   end function get_wetdep_scheme

end module processwetdepinterface_mod
```


