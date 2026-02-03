

# File DryDepCommon\_Mod.F90

[**File List**](files.md) **>** [**drydep**](dir_57fb5aa14ddb2cd518a6d90b65ffd000.md) **>** [**DryDepCommon\_Mod.F90**](_dry_dep_common___mod_8_f90.md)

[Go to the documentation of this file](_dry_dep_common___mod_8_f90.md)


```Fortran


module drydepcommon_mod

   use precision_mod, only: fp
   ! use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure, cc_error, cc_warning, errormanagertype, &
      error_invalid_config, error_invalid_state, error_not_found
   use configmanager_mod, only: configmanagertype  ! ConfigManager integration
   use statemanager_mod, only: statemanagertype  ! Add StateManager integration

   implicit none
   private

   ! Export types
   public :: drydepprocessconfig  ! New unified process config
   public :: drydepconfig
   public :: drydepschemeweselyconfig
   public :: drydepschemegocartconfig
   public :: drydepschemezhangconfig

   ! Export utility functions
   public :: int_to_string

   type :: drydepconfig

      ! Process settings
      character(len=32) :: gas_scheme = 'wesely'   ! Scheme for gas species
      character(len=32) :: aero_scheme = 'gocart'  ! Scheme for aerosol species
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
      integer, allocatable :: species_indices(:)  ! Indices of drydep species in ChemState
      logical, allocatable :: is_gas(:)           ! Is gas species (true for gas, false for aerosol)



      ! Species properties
      real(fp), allocatable :: species_dd_DvzAerSnow(:)      ! dd_DvzAerSnow for each species
      real(fp), allocatable :: species_dd_DvzMinVal_land(:)      ! dd_DvzMinVal_land for each species
      real(fp), allocatable :: species_dd_DvzMinVal_snow(:)      ! dd_DvzMinVal_snow for each species
      real(fp), allocatable :: species_dd_f0(:)      ! dd_f0 for each species
      real(fp), allocatable :: species_dd_hstar(:)      ! dd_hstar for each species
      real(fp), allocatable :: species_density(:)      ! density for each species
      logical, allocatable :: species_is_dust(:)      ! is_dust for each species
      logical, allocatable :: species_is_seasalt(:)      ! is_seasalt for each species
      real(fp), allocatable :: species_lower_radius(:)      ! lower_radius for each species
      real(fp), allocatable :: species_mw_g(:)      ! mw_g for each species
      real(fp), allocatable :: species_radius(:)      ! radius for each species
      real(fp), allocatable :: species_upper_radius(:)      ! upper_radius for each species

      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   contains
      procedure, public :: validate => validate_drydep_config
      procedure, public :: finalize => finalize_drydep_config
      procedure, public :: print_summary => print_drydep_config_summary
   end type drydepconfig

   type :: drydepschemeweselyconfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'wesely'
      character(len=256) :: description = 'Wesely 1989 gas dry deposition scheme'
      character(len=64) :: author = 'Wei Li'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .false.  ! Surface-only processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! DryDep velocity scale factor
      logical :: co2_effect = .true.  ! Apply CO2 effect on stomatal conductance
      real(fp) :: co2_level = 600.0  ! Ambient CO2 level for stomatal conductance adjustment
      real(fp) :: co2_reference = 380.0  ! Reference CO2 level for stomatal conductance adjustment

      ! Required meteorological fields
      integer :: n_required_met_fields = 21
      character(len=32) :: required_met_fields(21)

   contains
      procedure, public :: validate => validate_wesely_config
      procedure, public :: finalize => finalize_wesely_config
   end type drydepschemeweselyconfig

   ! wesely scheme uses local variables only - no persistent state type needed

   type :: drydepschemegocartconfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'gocart'
      character(len=256) :: description = 'GOCART-2G aerosol dry deposition scheme'
      character(len=64) :: author = 'Wei Li & Lacey Holland'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .false.  ! Surface-only processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Dry deposition velocity scale factor
      logical :: resuspension = .false.  ! Apply resuspension for dry deposition

      ! Required meteorological fields
      integer :: n_required_met_fields = 13
      character(len=32) :: required_met_fields(13)

   contains
      procedure, public :: validate => validate_gocart_config
      procedure, public :: finalize => finalize_gocart_config
   end type drydepschemegocartconfig

   ! gocart scheme uses local variables only - no persistent state type needed

   type :: drydepschemezhangconfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'zhang'
      character(len=256) :: description = 'Zhang et al. [2001] scheme with Emerson et al. [2020] updates'
      character(len=64) :: author = 'Wei Li'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .false.  ! Surface-only processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Dry deposition velocity scale factor

      ! Required meteorological fields
      integer :: n_required_met_fields = 15
      character(len=32) :: required_met_fields(15)

   contains
      procedure, public :: validate => validate_zhang_config
      procedure, public :: finalize => finalize_zhang_config
   end type drydepschemezhangconfig

   ! zhang scheme uses local variables only - no persistent state type needed


   type :: drydepprocessconfig

      ! Process metadata
      character(len=64) :: process_name = 'drydep'
      character(len=16) :: process_version = '1.0.0'
      logical :: is_active = .true.

      ! Process-specific configuration (delegate to DryDepConfig)
      type(DryDepConfig) :: drydep_config

      ! Scheme configurations
      ! Separate gas and aerosol scheme configurations
      type(DryDepSchemeWESELYConfig) :: wesely_config
      type(DryDepSchemeGOCARTConfig) :: gocart_config
      type(DryDepSchemeZHANGConfig) :: zhang_config

   contains
      procedure, public :: load_from_config => drydep_process_load_config
      procedure, public :: load_species_from_chem_state => load_species_from_chem_state
      procedure, public :: validate => drydep_process_validate
      procedure, public :: finalize => drydep_process_finalize
      procedure, public :: get_active_scheme_config => get_active_scheme_config
      procedure, public :: load_wesely_config
      procedure, public :: load_gocart_config
      procedure, public :: load_zhang_config
      procedure, public :: map_diagnostic_species_indices
   end type drydepprocessconfig

contains

   subroutine validate_drydep_config(this, error_handler)
      class(DryDepConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      character(len=256) :: error_msg
      integer :: rc

      ! Validate time step bounds
      if (this%dt_min <= 0.0_fp) then
         call error_handler%report_error(error_invalid_config, &
            "Minimum time step must be positive", rc)
         return
      end if

      if (this%dt_max < this%dt_min) then
         call error_handler%report_error(error_invalid_config, &
            "Maximum time step must be >= minimum time step", rc)
         return
      end if

      ! Validate active scheme(s)
      ! Validate gas scheme
      if (trim(this%gas_scheme) /= 'wesely' .and. &
         .true.) then
         write(error_msg, '(A)') "Invalid gas scheme: " // trim(this%gas_scheme)
         call error_handler%report_error(error_invalid_config, error_msg, rc)
         return
      end if

      ! Validate aerosol scheme
      if (trim(this%aero_scheme) /= 'gocart' .and. &
         trim(this%aero_scheme) /= 'zhang' .and. &
         .true.) then
         write(error_msg, '(A)') "Invalid aerosol scheme: " // trim(this%aero_scheme)
         call error_handler%report_error(error_invalid_config, error_msg, rc)
         return
      end if

   end subroutine validate_drydep_config

   subroutine print_drydep_config_summary(this)
      class(DryDepConfig), intent(in) :: this

      write(*, '(A)') "=== DryDep Process Configuration ==="
      write(*, '(A,A)') "  Gas scheme: ", trim(this%gas_scheme)
      write(*, '(A,A)') "  Aerosol scheme: ", trim(this%aero_scheme)
      write(*, '(A,I0)') "  Number of species: ", this%n_species
      write(*, '(A,F0.1,A)') "  Minimum time step: ", this%dt_min, " s"
      write(*, '(A,F0.1,A)') "  Maximum time step: ", this%dt_max, " s"
      write(*, '(A,L1)') "  Output diagnostics: ", this%output_diagnostics
      write(*, '(A)') "============================================="

   end subroutine print_drydep_config_summary

   subroutine finalize_drydep_config(this)
      class(DryDepConfig), intent(inout) :: this

      ! Deallocate species names array
      if (allocated(this%species_names)) then
         deallocate(this%species_names)
      end if

      ! Deallocate species indices array
      if (allocated(this%species_indices)) then
         deallocate(this%species_indices)
      end if

      ! Deallocate species properties arrays
      if (allocated(this%species_dd_DvzAerSnow)) then
         deallocate(this%species_dd_DvzAerSnow)
      end if
      if (allocated(this%species_dd_DvzMinVal_land)) then
         deallocate(this%species_dd_DvzMinVal_land)
      end if
      if (allocated(this%species_dd_DvzMinVal_snow)) then
         deallocate(this%species_dd_DvzMinVal_snow)
      end if
      if (allocated(this%species_dd_f0)) then
         deallocate(this%species_dd_f0)
      end if
      if (allocated(this%species_dd_hstar)) then
         deallocate(this%species_dd_hstar)
      end if
      if (allocated(this%species_density)) then
         deallocate(this%species_density)
      end if
      if (allocated(this%species_is_dust)) then
         deallocate(this%species_is_dust)
      end if
      if (allocated(this%species_is_seasalt)) then
         deallocate(this%species_is_seasalt)
      end if
      if (allocated(this%species_lower_radius)) then
         deallocate(this%species_lower_radius)
      end if
      if (allocated(this%species_mw_g)) then
         deallocate(this%species_mw_g)
      end if
      if (allocated(this%species_radius)) then
         deallocate(this%species_radius)
      end if
      if (allocated(this%species_upper_radius)) then
         deallocate(this%species_upper_radius)
      end if

      ! Deallocate gas/aerosol classification array
      if (allocated(this%is_gas)) then
         deallocate(this%is_gas)
      end if

      ! Deallocate diagnostic species array
      if (allocated(this%diagnostic_species)) then
         deallocate(this%diagnostic_species)
      end if

      ! Deallocate diagnostic species indices array
      if (allocated(this%diagnostic_species_id)) then
         deallocate(this%diagnostic_species_id)
      end if

   end subroutine finalize_drydep_config

   subroutine validate_wesely_config(this, error_handler)
      class(DryDepSchemeWESELYConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_wesely_config

   subroutine finalize_wesely_config(this)
      class(DryDepSchemeWESELYConfig), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_wesely_config


   subroutine validate_gocart_config(this, error_handler)
      class(DryDepSchemeGOCARTConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_gocart_config

   subroutine finalize_gocart_config(this)
      class(DryDepSchemeGOCARTConfig), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_gocart_config


   subroutine validate_zhang_config(this, error_handler)
      class(DryDepSchemeZHANGConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_zhang_config

   subroutine finalize_zhang_config(this)
      class(DryDepSchemeZHANGConfig), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_zhang_config



   function int_to_string(int_val) result(str_val)
      integer, intent(in) :: int_val
      character(len=32) :: str_val

      write(str_val, '(I0)') int_val
      str_val = adjustl(str_val)

   end function int_to_string

   subroutine drydep_process_load_config(this, config_manager, error_handler)
      class(DryDepProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      character(len=256) :: scheme_name
      integer :: ierr, rc

      ! Process reads directly from master YAML structure: processes.drydep
      ! ConfigManager provides generic YAML access, process handles its own configuration

      ! Load process metadata
      call config_manager%get_string("processes/drydep/name", this%process_name, rc, "drydep")
      if (rc /= cc_success) this%process_name = "drydep"  ! default

      call config_manager%get_string("processes/drydep/version", this%process_version, rc, "1.0.0")
      if (rc /= cc_success) this%process_version = "1.0.0"  ! default

      call config_manager%get_logical("processes/drydep/activate", this%is_active, rc, .true.)
      if (rc /= cc_success) this%is_active = .true.  ! default

      ! Load process-specific configuration directly from master YAML
      call config_manager%get_string("processes/drydep/gas_scheme", this%drydep_config%gas_scheme, rc, "wesely")
      if (rc /= cc_success) then
         call error_handler%report_error(error_invalid_config, &
            "Missing required 'gas_scheme' in processes/drydep configuration", rc)
         return
      end if

      call config_manager%get_string("processes/drydep/aero_scheme", this%drydep_config%aero_scheme, rc, "")
      if (rc /= cc_success) then
         call error_handler%report_error(error_invalid_config, &
            "Missing required 'aero_scheme' in processes/drydep configuration", rc)
         return
      end if

      ! Load diagnostic switch
      call config_manager%get_logical("processes/drydep/diagnostics", this%drydep_config%diagnostics, rc, .false.)
      if (rc /= cc_success) this%drydep_config%diagnostics = .false.  ! Default

      ! Load diagnostic species list
      call config_manager%get_array("processes/drydep/diag_species", this%drydep_config%diagnostic_species, &
         rc, default_values=["All"])
      if (rc /= cc_success) then
         ! Default to all species if not specified
         allocate(this%drydep_config%diagnostic_species(1))
         this%drydep_config%diagnostic_species(1) = "All"
         this%drydep_config%n_diagnostic_species = 1
      else
         ! Set the count based on the returned array size
         if (allocated(this%drydep_config%diagnostic_species)) then
            this%drydep_config%n_diagnostic_species = size(this%drydep_config%diagnostic_species)
         else
            this%drydep_config%n_diagnostic_species = 0
         end if
      end if

      ! Species configuration is loaded from ChemState in load_species_from_chem_state
      ! The species come from the master species YAML file (CATChem_species.yml)
      ! and are filtered by is_drydep property


      ! Load scheme-specific configuration from master YAML
      ! Load gas scheme configuration
      scheme_name = trim(this%drydep_config%gas_scheme)
      select case (scheme_name)
       case ('wesely')
         call this%load_wesely_config(config_manager, error_handler)
       case default
         call error_handler%report_error(error_invalid_state, &
            "Unknown drydep gas scheme: " // trim(scheme_name), rc)
         return
      end select

      ! Load aerosol scheme configuration
      scheme_name = trim(this%drydep_config%aero_scheme)
      select case (scheme_name)
       case ('gocart')
         call this%load_gocart_config(config_manager, error_handler)
       case ('zhang')
         call this%load_zhang_config(config_manager, error_handler)
       case default
         call error_handler%report_error(error_invalid_state, &
            "Unknown drydep aerosol scheme: " // trim(scheme_name), rc)
         return
      end select

   end subroutine drydep_process_load_config


   subroutine load_species_from_chem_state(this, chem_state, error_handler)
      use chemstate_mod, only: chemstatetype

      class(DryDepProcessConfig), intent(inout) :: this
      type(ChemStateType), pointer, intent(in) :: chem_state
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: i, rc
      integer :: species_idx

      if (.not. associated(chem_state)) then
         call error_handler%report_error(error_invalid_state, &
            "ChemState not associated in load_species_from_chem_state", rc)
         return
      end if

      ! by_metadata mode: Load species by type from ChemState using dynamic metadata flag mapping
      ! Dynamic mapping: is_drydep -> nSpeciesDryDep and DryDepIndex
      this%drydep_config%n_species = chem_state%nSpeciesDryDep

      if (this%drydep_config%n_species <= 0) then
         call error_handler%report_error(error_invalid_state, &
            "No drydep species found in ChemState", rc)
         return
      end if

      ! Check if DryDepIndex is allocated and has correct size
      if (.not. allocated(chem_state%DryDepIndex)) then
         call error_handler%report_error(error_invalid_state, &
            "DryDepIndex not allocated in ChemState", rc)
         return
      end if

      if (size(chem_state%DryDepIndex) < this%drydep_config%n_species) then
         call error_handler%report_error(error_invalid_state, &
            "DryDepIndex size inconsistent with nSpeciesDryDep", rc)
         return
      end if

      ! Deallocate existing arrays if allocated
      if (allocated(this%drydep_config%species_names)) then
         deallocate(this%drydep_config%species_names)
      end if
      if (allocated(this%drydep_config%species_indices)) then
         deallocate(this%drydep_config%species_indices)
      end if

      ! Allocate arrays
      allocate(this%drydep_config%species_names(this%drydep_config%n_species))
      allocate(this%drydep_config%species_indices(this%drydep_config%n_species))
      allocate(this%drydep_config%is_gas(this%drydep_config%n_species))

      ! Allocate species properties arrays
      allocate(this%drydep_config%species_dd_DvzAerSnow(this%drydep_config%n_species))
      allocate(this%drydep_config%species_dd_DvzMinVal_land(this%drydep_config%n_species))
      allocate(this%drydep_config%species_dd_DvzMinVal_snow(this%drydep_config%n_species))
      allocate(this%drydep_config%species_dd_f0(this%drydep_config%n_species))
      allocate(this%drydep_config%species_dd_hstar(this%drydep_config%n_species))
      allocate(this%drydep_config%species_density(this%drydep_config%n_species))
      allocate(this%drydep_config%species_is_dust(this%drydep_config%n_species))
      allocate(this%drydep_config%species_is_seasalt(this%drydep_config%n_species))
      allocate(this%drydep_config%species_lower_radius(this%drydep_config%n_species))
      allocate(this%drydep_config%species_mw_g(this%drydep_config%n_species))
      allocate(this%drydep_config%species_radius(this%drydep_config%n_species))
      allocate(this%drydep_config%species_upper_radius(this%drydep_config%n_species))

      ! by_metadata mode: Copy indices from metadata-specific index array using dynamic mapping
      ! Dynamic mapping: is_drydep -> DryDepIndex
      this%drydep_config%species_indices(1:this%drydep_config%n_species) = &
         chem_state%DryDepIndex(1:this%drydep_config%n_species)

      ! Get species names using the indices
      do i = 1, this%drydep_config%n_species
         if (this%drydep_config%species_indices(i) > 0 .and. &
            this%drydep_config%species_indices(i) <= size(chem_state%SpeciesNames)) then
            this%drydep_config%species_names(i) = &
               trim(chem_state%SpeciesNames(this%drydep_config%species_indices(i)))
         else
            call error_handler%report_error(error_invalid_state, &
               "Invalid species index in species index array", rc)
            return
         end if
      end do

      ! Load species properties from ChemState
      do i = 1, this%drydep_config%n_species
         species_idx = this%drydep_config%species_indices(i)
         this%drydep_config%species_dd_DvzAerSnow(i) = chem_state%ChemSpecies(species_idx)%dd_DvzAerSnow
         this%drydep_config%species_dd_DvzMinVal_land(i) = chem_state%ChemSpecies(species_idx)%dd_DvzMinVal_land
         this%drydep_config%species_dd_DvzMinVal_snow(i) = chem_state%ChemSpecies(species_idx)%dd_DvzMinVal_snow
         this%drydep_config%species_dd_f0(i) = chem_state%ChemSpecies(species_idx)%dd_f0
         this%drydep_config%species_dd_hstar(i) = chem_state%ChemSpecies(species_idx)%dd_hstar
         this%drydep_config%species_density(i) = chem_state%ChemSpecies(species_idx)%density
         this%drydep_config%species_is_dust(i) = chem_state%ChemSpecies(species_idx)%is_dust
         this%drydep_config%species_is_seasalt(i) = chem_state%ChemSpecies(species_idx)%is_seasalt
         this%drydep_config%species_lower_radius(i) = chem_state%ChemSpecies(species_idx)%lower_radius
         this%drydep_config%species_mw_g(i) = chem_state%ChemSpecies(species_idx)%mw_g
         this%drydep_config%species_radius(i) = chem_state%ChemSpecies(species_idx)%radius
         this%drydep_config%species_upper_radius(i) = chem_state%ChemSpecies(species_idx)%upper_radius
         this%drydep_config%is_gas(i) = chem_state%ChemSpecies(species_idx)%is_gas
      end do

   end subroutine load_species_from_chem_state


   subroutine load_wesely_config(this, config_manager, error_handler)
      class(DryDepProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: ierr, rc

      ! Load scheme parameters directly from processes/drydep/wesely/ in master YAML
      call config_manager%get_real("processes/drydep/wesely/scale_factor", &
         this%wesely_config%scale_factor, rc, 1.0_fp)
      if (rc /= cc_success) this%wesely_config%scale_factor = 1.0_fp
      call config_manager%get_logical("processes/drydep/wesely/co2_effect", &
         this%wesely_config%co2_effect, rc, .true.)
      if (rc /= cc_success) this%wesely_config%co2_effect = .true.
      call config_manager%get_real("processes/drydep/wesely/co2_level", &
         this%wesely_config%co2_level, rc, 600.0_fp)
      if (rc /= cc_success) this%wesely_config%co2_level = 600.0_fp
      call config_manager%get_real("processes/drydep/wesely/co2_reference", &
         this%wesely_config%co2_reference, rc, 380.0_fp)
      if (rc /= cc_success) this%wesely_config%co2_reference = 380.0_fp


   end subroutine load_wesely_config

   subroutine load_gocart_config(this, config_manager, error_handler)
      class(DryDepProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: ierr, rc

      ! Load scheme parameters directly from processes/drydep/gocart/ in master YAML
      call config_manager%get_real("processes/drydep/gocart/scale_factor", &
         this%gocart_config%scale_factor, rc, 1.0_fp)
      if (rc /= cc_success) this%gocart_config%scale_factor = 1.0_fp
      call config_manager%get_logical("processes/drydep/gocart/resuspension", &
         this%gocart_config%resuspension, rc, .false.)
      if (rc /= cc_success) this%gocart_config%resuspension = .false.


   end subroutine load_gocart_config

   subroutine load_zhang_config(this, config_manager, error_handler)
      class(DryDepProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: ierr, rc

      ! Load scheme parameters directly from processes/drydep/zhang/ in master YAML
      call config_manager%get_real("processes/drydep/zhang/scale_factor", &
         this%zhang_config%scale_factor, rc, 1.0_fp)
      if (rc /= cc_success) this%zhang_config%scale_factor = 1.0_fp


   end subroutine load_zhang_config


   subroutine drydep_process_validate(this, state_manager, error_handler)
      class(DryDepProcessConfig), intent(inout) :: this
      type(StateManagerType), intent(in) :: state_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      ! Validate main config
      call this%drydep_config%validate(error_handler)

      ! Validate scheme-specific config
      ! Validate gas scheme config
      select case (trim(this%drydep_config%gas_scheme))
       case ('wesely')
         call this%wesely_config%validate(error_handler)
      end select

      ! Validate aerosol scheme config
      select case (trim(this%drydep_config%aero_scheme))
       case ('gocart')
         call this%gocart_config%validate(error_handler)
       case ('zhang')
         call this%zhang_config%validate(error_handler)
      end select

   end subroutine drydep_process_validate

   subroutine drydep_process_finalize(this)
      class(DryDepProcessConfig), intent(inout) :: this

      call this%drydep_config%finalize()
      call this%wesely_config%finalize()
      call this%gocart_config%finalize()
      call this%zhang_config%finalize()

   end subroutine drydep_process_finalize

   function get_active_scheme_config(this, gas_scheme) result(scheme_config)
      class(DryDepProcessConfig), intent(in) :: this
      logical, optional, intent(in) :: gas_scheme  ! If true, return gas scheme; if false, return aerosol scheme
      class(*), allocatable :: scheme_config

      ! For gas/aero differentiated processes, return requested scheme type
      if (present(gas_scheme) .and. gas_scheme) then
         ! Return gas scheme
         select case (trim(this%drydep_config%gas_scheme))
          case ('wesely')
            allocate(scheme_config, source=this%wesely_config)
          case default
            ! Return null
         end select
      else
         ! Return aerosol scheme (default or when gas_scheme = false)
         select case (trim(this%drydep_config%aero_scheme))
          case ('gocart')
            allocate(scheme_config, source=this%gocart_config)
          case ('zhang')
            allocate(scheme_config, source=this%zhang_config)
          case default
            ! Return null
         end select
      end if

   end function get_active_scheme_config

   subroutine map_diagnostic_species_indices(this, error_handler)
      class(DryDepProcessConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: i, j, rc
      character(len=256) :: error_msg
      logical :: found_species

      ! Only proceed if diagnostic species are defined
      if (this%drydep_config%n_diagnostic_species == 0) return

      ! Handle "All" case - map all available species
      if (this%drydep_config%n_diagnostic_species == 1 .and. &
         trim(this%drydep_config%diagnostic_species(1)) == "All") then

         ! Deallocate and reallocate for all species
         if (allocated(this%drydep_config%diagnostic_species_id)) deallocate(this%drydep_config%diagnostic_species_id)
         allocate(this%drydep_config%diagnostic_species_id(this%drydep_config%n_species))
         if (allocated(this%drydep_config%diagnostic_species)) deallocate(this%drydep_config%diagnostic_species)
         allocate(this%drydep_config%diagnostic_species(this%drydep_config%n_species))
         this%drydep_config%n_diagnostic_species = this%drydep_config%n_species
         this%drydep_config%diagnostic_species = this%drydep_config%species_names

         ! Map all species indices (1:n_species)
         do i = 1, this%drydep_config%n_species
            this%drydep_config%diagnostic_species_id(i) = i
         end do

         return
      end if

      ! Allocate diagnostic species indices array
      if (allocated(this%drydep_config%diagnostic_species_id)) deallocate(this%drydep_config%diagnostic_species_id)
      allocate(this%drydep_config%diagnostic_species_id(this%drydep_config%n_diagnostic_species))

      ! Map each diagnostic species name to its index in species_names
      do i = 1, this%drydep_config%n_diagnostic_species
         found_species = .false.

         do j = 1, this%drydep_config%n_species
            if (trim(this%drydep_config%diagnostic_species(i)) == trim(this%drydep_config%species_names(j))) then
               this%drydep_config%diagnostic_species_id(i) = j
               found_species = .true.
               exit
            end if
         end do

         if (.not. found_species) then
            write(error_msg, '(A,A,A)') "Diagnostic species '", &
               trim(this%drydep_config%diagnostic_species(i)), &
               "' not found in process species list"
            call error_handler%report_error(error_not_found, error_msg, rc)
            return
         end if
      end do

   end subroutine map_diagnostic_species_indices

end module drydepcommon_mod
```


