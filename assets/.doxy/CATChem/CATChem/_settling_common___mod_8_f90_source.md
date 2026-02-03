

# File SettlingCommon\_Mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**settling**](dir_1a0bba2ffdf6e6637fcb76856471cb75.md) **>** [**SettlingCommon\_Mod.F90**](_settling_common___mod_8_f90.md)

[Go to the documentation of this file](_settling_common___mod_8_f90.md)


```Fortran


module settlingcommon_mod

   use precision_mod, only: fp
   ! use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure, cc_error, cc_warning, errormanagertype, &
      error_invalid_config, error_invalid_state, error_not_found
   use configmanager_mod, only: configmanagertype  ! ConfigManager integration
   use statemanager_mod, only: statemanagertype  ! Add StateManager integration

   implicit none
   private

   ! Export types
   public :: settlingprocessconfig  ! New unified process config
   public :: settlingconfig
   public :: settlingschemegocartconfig

   ! Export utility functions
   public :: int_to_string

   type :: settlingconfig

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
      integer, allocatable :: species_indices(:)  ! Indices of settling species in ChemState



      ! Species properties
      real(fp), allocatable :: species_density(:)      ! density for each species
      real(fp), allocatable :: species_mie_map(:)      ! mie_map for each species
      real(fp), allocatable :: species_radius(:)      ! radius for each species

      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   contains
      procedure, public :: validate => validate_settling_config
      procedure, public :: finalize => finalize_settling_config
      procedure, public :: print_summary => print_settling_config_summary
   end type settlingconfig

   type :: settlingschemegocartconfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'gocart'
      character(len=256) :: description = 'GOCART gravitational settling scheme'
      character(len=64) :: author = 'Wei Li'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .true.  ! Full column processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! settling velocity factor
      logical :: simple_scheme = .false.  ! read in mie data for wet particles if true; otherwise calculate particles wet swelling internally
      integer :: swelling_method = 1  ! method for calculating particle swelling: 1 Fitzgerald 1975; 2 for Gerber 1985
      logical :: correction_maring = .false.  ! correct the settling velocity following Maring et al, 2003

      ! Required meteorological fields
      integer :: n_required_met_fields = 7
      character(len=32) :: required_met_fields(7)

   contains
      procedure, public :: validate => validate_gocart_config
      procedure, public :: finalize => finalize_gocart_config
   end type settlingschemegocartconfig

   ! gocart scheme uses local variables only - no persistent state type needed


   type :: settlingprocessconfig

      ! Process metadata
      character(len=64) :: process_name = 'settling'
      character(len=16) :: process_version = '1.0.0'
      logical :: is_active = .true.

      ! Process-specific configuration (delegate to SettlingConfig)
      type(SettlingConfig) :: settling_config

      ! Scheme configurations
      type(SettlingSchemeGOCARTConfig) :: gocart_config

   contains
      procedure, public :: load_from_config => settling_process_load_config
      procedure, public :: load_species_from_chem_state => load_species_from_chem_state
      procedure, public :: validate => settling_process_validate
      procedure, public :: finalize => settling_process_finalize
      procedure, public :: get_active_scheme_config => get_active_scheme_config
      procedure, public :: load_gocart_config
      procedure, public :: map_diagnostic_species_indices
   end type settlingprocessconfig

contains

   subroutine validate_settling_config(this, error_handler)
      class(SettlingConfig), intent(inout) :: this
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
      ! Validate scheme
      if (trim(this%scheme) /= 'gocart' .and. &
         .true.) then
         write(error_msg, '(A)') "Invalid scheme: " // trim(this%scheme)
         call error_handler%report_error(error_invalid_config, error_msg, rc)
         return
      end if

   end subroutine validate_settling_config

   subroutine print_settling_config_summary(this)
      class(SettlingConfig), intent(in) :: this

      write(*, '(A)') "=== Settling Process Configuration ==="
      write(*, '(A,A)') "  Active scheme: ", trim(this%scheme)
      write(*, '(A,I0)') "  Number of species: ", this%n_species
      write(*, '(A,F0.1,A)') "  Minimum time step: ", this%dt_min, " s"
      write(*, '(A,F0.1,A)') "  Maximum time step: ", this%dt_max, " s"
      write(*, '(A,L1)') "  Output diagnostics: ", this%output_diagnostics
      write(*, '(A)') "============================================="

   end subroutine print_settling_config_summary

   subroutine finalize_settling_config(this)
      class(SettlingConfig), intent(inout) :: this

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
      if (allocated(this%species_mie_map)) then
         deallocate(this%species_mie_map)
      end if
      if (allocated(this%species_radius)) then
         deallocate(this%species_radius)
      end if


      ! Deallocate diagnostic species array
      if (allocated(this%diagnostic_species)) then
         deallocate(this%diagnostic_species)
      end if

      ! Deallocate diagnostic species indices array
      if (allocated(this%diagnostic_species_id)) then
         deallocate(this%diagnostic_species_id)
      end if

   end subroutine finalize_settling_config

   subroutine validate_gocart_config(this, error_handler)
      class(SettlingSchemeGOCARTConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_gocart_config

   subroutine finalize_gocart_config(this)
      class(SettlingSchemeGOCARTConfig), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_gocart_config



   function int_to_string(int_val) result(str_val)
      integer, intent(in) :: int_val
      character(len=32) :: str_val

      write(str_val, '(I0)') int_val
      str_val = adjustl(str_val)

   end function int_to_string

   subroutine settling_process_load_config(this, config_manager, error_handler)
      class(SettlingProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      character(len=256) :: scheme_name
      integer :: ierr, rc

      ! Process reads directly from master YAML structure: processes.settling
      ! ConfigManager provides generic YAML access, process handles its own configuration

      ! Load process metadata
      call config_manager%get_string("processes/settling/name", this%process_name, rc, "settling")
      if (rc /= cc_success) this%process_name = "settling"  ! default

      call config_manager%get_string("processes/settling/version", this%process_version, rc, "1.0.0")
      if (rc /= cc_success) this%process_version = "1.0.0"  ! default

      call config_manager%get_logical("processes/settling/activate", this%is_active, rc, .true.)
      if (rc /= cc_success) this%is_active = .true.  ! default

      ! Load process-specific configuration directly from master YAML
      call config_manager%get_string("processes/settling/scheme", this%settling_config%scheme, rc, "gocart")
      if (rc /= cc_success) then
         call error_handler%report_error(error_invalid_config, &
            "Missing required 'scheme' in processes/settling configuration", rc)
         return
      end if

      ! Load diagnostic switch
      call config_manager%get_logical("processes/settling/diagnostics", this%settling_config%diagnostics, rc, .false.)
      if (rc /= cc_success) this%settling_config%diagnostics = .false.  ! Default

      ! Load diagnostic species list
      call config_manager%get_array("processes/settling/diag_species", this%settling_config%diagnostic_species, &
         rc, default_values=["All"])
      if (rc /= cc_success) then
         ! Default to all species if not specified
         allocate(this%settling_config%diagnostic_species(1))
         this%settling_config%diagnostic_species(1) = "All"
         this%settling_config%n_diagnostic_species = 1
      else
         ! Set the count based on the returned array size
         if (allocated(this%settling_config%diagnostic_species)) then
            this%settling_config%n_diagnostic_species = size(this%settling_config%diagnostic_species)
         else
            this%settling_config%n_diagnostic_species = 0
         end if
      end if

      ! Species configuration is loaded from ChemState in load_species_from_chem_state
      ! The species come from the master species YAML file (CATChem_species.yml)
      ! and are filtered by is_settling property


      ! Load scheme-specific configuration from master YAML
      scheme_name = trim(this%settling_config%scheme)
      select case (scheme_name)
       case ('gocart')
         call this%load_gocart_config(config_manager, error_handler)
       case default
         call error_handler%report_error(error_invalid_state, &
            "Unknown settling scheme: " // trim(scheme_name), rc)
         return
      end select

   end subroutine settling_process_load_config


   subroutine load_species_from_chem_state(this, chem_state, error_handler)
      use chemstate_mod, only: chemstatetype

      class(SettlingProcessConfig), intent(inout) :: this
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
      ! Dynamic mapping: is_aerosol -> nSpeciesAero and AeroIndex
      this%settling_config%n_species = chem_state%nSpeciesAero

      if (this%settling_config%n_species <= 0) then
         call error_handler%report_error(error_invalid_state, &
            "No aerosol species found in ChemState", rc)
         return
      end if

      ! Check if AeroIndex is allocated and has correct size
      if (.not. allocated(chem_state%AeroIndex)) then
         call error_handler%report_error(error_invalid_state, &
            "AeroIndex not allocated in ChemState", rc)
         return
      end if

      if (size(chem_state%AeroIndex) < this%settling_config%n_species) then
         call error_handler%report_error(error_invalid_state, &
            "AeroIndex size inconsistent with nSpeciesAero", rc)
         return
      end if

      ! Deallocate existing arrays if allocated
      if (allocated(this%settling_config%species_names)) then
         deallocate(this%settling_config%species_names)
      end if
      if (allocated(this%settling_config%species_indices)) then
         deallocate(this%settling_config%species_indices)
      end if

      ! Allocate arrays
      allocate(this%settling_config%species_names(this%settling_config%n_species))
      allocate(this%settling_config%species_indices(this%settling_config%n_species))

      ! Allocate species properties arrays
      allocate(this%settling_config%species_density(this%settling_config%n_species))
      allocate(this%settling_config%species_mie_map(this%settling_config%n_species))
      allocate(this%settling_config%species_radius(this%settling_config%n_species))

      ! by_metadata mode: Copy indices from metadata-specific index array using dynamic mapping
      ! Dynamic mapping: is_aerosol -> AeroIndex
      this%settling_config%species_indices(1:this%settling_config%n_species) = &
         chem_state%AeroIndex(1:this%settling_config%n_species)

      ! Get species names using the indices
      do i = 1, this%settling_config%n_species
         if (this%settling_config%species_indices(i) > 0 .and. &
            this%settling_config%species_indices(i) <= size(chem_state%SpeciesNames)) then
            this%settling_config%species_names(i) = &
               trim(chem_state%SpeciesNames(this%settling_config%species_indices(i)))
         else
            call error_handler%report_error(error_invalid_state, &
               "Invalid species index in species index array", rc)
            return
         end if
      end do

      ! Load species properties from ChemState
      do i = 1, this%settling_config%n_species
         species_idx = this%settling_config%species_indices(i)
         this%settling_config%species_density(i) = chem_state%ChemSpecies(species_idx)%density
         if (allocated(chem_state%SpcMieMap)) then
            this%settling_config%species_mie_map(i) = chem_state%SpcMieMap(species_idx)
         else
            this%settling_config%species_mie_map(i) = -1  ! Default or error value
         end if
         this%settling_config%species_radius(i) = chem_state%ChemSpecies(species_idx)%radius
      end do

   end subroutine load_species_from_chem_state


   subroutine load_gocart_config(this, config_manager, error_handler)
      class(SettlingProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: ierr, rc

      ! Load scheme parameters directly from processes/settling/gocart/ in master YAML
      call config_manager%get_real("processes/settling/gocart/scale_factor", &
         this%gocart_config%scale_factor, rc, 1.0_fp)
      if (rc /= cc_success) this%gocart_config%scale_factor = 1.0_fp
      call config_manager%get_logical("processes/settling/gocart/simple_scheme", &
         this%gocart_config%simple_scheme, rc, .false.)
      if (rc /= cc_success) this%gocart_config%simple_scheme = .false.
      call config_manager%get_integer("processes/settling/gocart/swelling_method", &
         this%gocart_config%swelling_method, rc, 1)
      if (rc /= cc_success) this%gocart_config%swelling_method = 1
      call config_manager%get_logical("processes/settling/gocart/correction_maring", &
         this%gocart_config%correction_maring, rc, .false.)
      if (rc /= cc_success) this%gocart_config%correction_maring = .false.


   end subroutine load_gocart_config


   subroutine settling_process_validate(this, state_manager, error_handler)
      class(SettlingProcessConfig), intent(inout) :: this
      type(StateManagerType), intent(in) :: state_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      ! Validate main config
      call this%settling_config%validate(error_handler)

      ! Validate scheme-specific config
      select case (trim(this%settling_config%scheme))
       case ('gocart')
         call this%gocart_config%validate(error_handler)
      end select

   end subroutine settling_process_validate

   subroutine settling_process_finalize(this)
      class(SettlingProcessConfig), intent(inout) :: this

      call this%settling_config%finalize()
      call this%gocart_config%finalize()

   end subroutine settling_process_finalize

   function get_active_scheme_config(this) result(scheme_config)
      class(SettlingProcessConfig), intent(in) :: this
      class(*), allocatable :: scheme_config

      select case (trim(this%settling_config%scheme))
       case ('gocart')
         allocate(scheme_config, source=this%gocart_config)
       case default
         ! Return null
      end select

   end function get_active_scheme_config

   subroutine map_diagnostic_species_indices(this, error_handler)
      class(SettlingProcessConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: i, j, rc
      character(len=256) :: error_msg
      logical :: found_species

      ! Only proceed if diagnostic species are defined
      if (this%settling_config%n_diagnostic_species == 0) return

      ! Handle "All" case - map all available species
      if (this%settling_config%n_diagnostic_species == 1 .and. &
         trim(this%settling_config%diagnostic_species(1)) == "All") then

         ! Deallocate and reallocate for all species
         if (allocated(this%settling_config%diagnostic_species_id)) deallocate(this%settling_config%diagnostic_species_id)
         allocate(this%settling_config%diagnostic_species_id(this%settling_config%n_species))
         if (allocated(this%settling_config%diagnostic_species)) deallocate(this%settling_config%diagnostic_species)
         allocate(this%settling_config%diagnostic_species(this%settling_config%n_species))
         this%settling_config%n_diagnostic_species = this%settling_config%n_species
         this%settling_config%diagnostic_species = this%settling_config%species_names

         ! Map all species indices (1:n_species)
         do i = 1, this%settling_config%n_species
            this%settling_config%diagnostic_species_id(i) = i
         end do

         return
      end if

      ! Allocate diagnostic species indices array
      if (allocated(this%settling_config%diagnostic_species_id)) deallocate(this%settling_config%diagnostic_species_id)
      allocate(this%settling_config%diagnostic_species_id(this%settling_config%n_diagnostic_species))

      ! Map each diagnostic species name to its index in species_names
      do i = 1, this%settling_config%n_diagnostic_species
         found_species = .false.

         do j = 1, this%settling_config%n_species
            if (trim(this%settling_config%diagnostic_species(i)) == trim(this%settling_config%species_names(j))) then
               this%settling_config%diagnostic_species_id(i) = j
               found_species = .true.
               exit
            end if
         end do

         if (.not. found_species) then
            write(error_msg, '(A,A,A)') "Diagnostic species '", &
               trim(this%settling_config%diagnostic_species(i)), &
               "' not found in process species list"
            call error_handler%report_error(error_not_found, error_msg, rc)
            return
         end if
      end do

   end subroutine map_diagnostic_species_indices

end module settlingcommon_mod
```


