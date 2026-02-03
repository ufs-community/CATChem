

# File WetDepCommon\_Mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**wetdep**](dir_8b9a0ce556ea4a65f6920dfb49dcd69d.md) **>** [**WetDepCommon\_Mod.F90**](_wet_dep_common___mod_8_f90.md)

[Go to the documentation of this file](_wet_dep_common___mod_8_f90.md)


```Fortran


module wetdepcommon_mod

   use precision_mod, only: fp
   ! use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure, cc_error, cc_warning, errormanagertype, &
      error_invalid_config, error_invalid_state, error_not_found
   use configmanager_mod, only: configmanagertype  ! ConfigManager integration
   use statemanager_mod, only: statemanagertype  ! Add StateManager integration

   implicit none
   private

   ! Export types
   public :: wetdepprocessconfig  ! New unified process config
   public :: wetdepconfig
   public :: wetdepschemejacobconfig

   ! Export utility functions
   public :: int_to_string

   type :: wetdepconfig

      ! Process settings
      character(len=32) :: scheme = 'jacob'
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
      integer, allocatable :: species_indices(:)  ! Indices of wetdep species in ChemState



      ! Species properties
      real(fp), allocatable :: species_henry_cr(:)      ! henry_cr for each species
      real(fp), allocatable :: species_henry_k0(:)      ! henry_k0 for each species
      real(fp), allocatable :: species_henry_pKa(:)      ! henry_pKa for each species
      logical, allocatable :: species_is_aerosol(:)      ! is_aerosol for each species
      real(fp), allocatable :: species_mw_g(:)      ! mw_g for each species
      real(fp), allocatable :: species_radius(:)      ! radius for each species
      logical, allocatable :: species_wd_LiqAndGas(:)      ! wd_LiqAndGas for each species
      real(fp), allocatable :: species_wd_convfacI2G(:)      ! wd_convfacI2G for each species
      real(fp), allocatable :: species_wd_rainouteff(:,:)      ! wd_rainouteff for each species
      real(fp), allocatable :: species_wd_retfactor(:)      ! wd_retfactor for each species

      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   contains
      procedure, public :: validate => validate_wetdep_config
      procedure, public :: finalize => finalize_wetdep_config
      procedure, public :: print_summary => print_wetdep_config_summary
   end type wetdepconfig

   type :: wetdepschemejacobconfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'jacob'
      character(len=256) :: description = 'Jacob et al. [2000] wet deposition scheme'
      character(len=64) :: author = 'Wei Li'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .true.  ! Full column processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Washout tuning factor
      real(fp) :: radius_threshold = 1.0  ! Radius threshold for aerosol wet deposition (um)

      ! Required meteorological fields
      integer :: n_required_met_fields = 8
      character(len=32) :: required_met_fields(8)

   contains
      procedure, public :: validate => validate_jacob_config
      procedure, public :: finalize => finalize_jacob_config
   end type wetdepschemejacobconfig

   ! jacob scheme uses local variables only - no persistent state type needed


   type :: wetdepprocessconfig

      ! Process metadata
      character(len=64) :: process_name = 'wetdep'
      character(len=16) :: process_version = '1.0.0'
      logical :: is_active = .true.

      ! Process-specific configuration (delegate to WetDepConfig)
      type(WetDepConfig) :: wetdep_config

      ! Scheme configurations
      type(WetDepSchemeJACOBConfig) :: jacob_config

   contains
      procedure, public :: load_from_config => wetdep_process_load_config
      procedure, public :: load_species_from_chem_state => load_species_from_chem_state
      procedure, public :: validate => wetdep_process_validate
      procedure, public :: finalize => wetdep_process_finalize
      procedure, public :: get_active_scheme_config => get_active_scheme_config
      procedure, public :: load_jacob_config
      procedure, public :: map_diagnostic_species_indices
   end type wetdepprocessconfig

contains

   subroutine validate_wetdep_config(this, error_handler)
      class(WetDepConfig), intent(inout) :: this
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
      if (trim(this%scheme) /= 'jacob' .and. &
         .true.) then
         write(error_msg, '(A)') "Invalid scheme: " // trim(this%scheme)
         call error_handler%report_error(error_invalid_config, error_msg, rc)
         return
      end if

   end subroutine validate_wetdep_config

   subroutine print_wetdep_config_summary(this)
      class(WetDepConfig), intent(in) :: this

      write(*, '(A)') "=== WetDep Process Configuration ==="
      write(*, '(A,A)') "  Active scheme: ", trim(this%scheme)
      write(*, '(A,I0)') "  Number of species: ", this%n_species
      write(*, '(A,F0.1,A)') "  Minimum time step: ", this%dt_min, " s"
      write(*, '(A,F0.1,A)') "  Maximum time step: ", this%dt_max, " s"
      write(*, '(A,L1)') "  Output diagnostics: ", this%output_diagnostics
      write(*, '(A)') "============================================="

   end subroutine print_wetdep_config_summary

   subroutine finalize_wetdep_config(this)
      class(WetDepConfig), intent(inout) :: this

      ! Deallocate species names array
      if (allocated(this%species_names)) then
         deallocate(this%species_names)
      end if

      ! Deallocate species indices array
      if (allocated(this%species_indices)) then
         deallocate(this%species_indices)
      end if

      ! Deallocate species properties arrays
      if (allocated(this%species_henry_cr)) then
         deallocate(this%species_henry_cr)
      end if
      if (allocated(this%species_henry_k0)) then
         deallocate(this%species_henry_k0)
      end if
      if (allocated(this%species_henry_pKa)) then
         deallocate(this%species_henry_pKa)
      end if
      if (allocated(this%species_is_aerosol)) then
         deallocate(this%species_is_aerosol)
      end if
      if (allocated(this%species_mw_g)) then
         deallocate(this%species_mw_g)
      end if
      if (allocated(this%species_radius)) then
         deallocate(this%species_radius)
      end if
      if (allocated(this%species_wd_LiqAndGas)) then
         deallocate(this%species_wd_LiqAndGas)
      end if
      if (allocated(this%species_wd_convfacI2G)) then
         deallocate(this%species_wd_convfacI2G)
      end if
      if (allocated(this%species_wd_rainouteff)) then
         deallocate(this%species_wd_rainouteff)
      end if
      if (allocated(this%species_wd_retfactor)) then
         deallocate(this%species_wd_retfactor)
      end if


      ! Deallocate diagnostic species array
      if (allocated(this%diagnostic_species)) then
         deallocate(this%diagnostic_species)
      end if

      ! Deallocate diagnostic species indices array
      if (allocated(this%diagnostic_species_id)) then
         deallocate(this%diagnostic_species_id)
      end if

   end subroutine finalize_wetdep_config

   subroutine validate_jacob_config(this, error_handler)
      class(WetDepSchemeJACOBConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_jacob_config

   subroutine finalize_jacob_config(this)
      class(WetDepSchemeJACOBConfig), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_jacob_config



   function int_to_string(int_val) result(str_val)
      integer, intent(in) :: int_val
      character(len=32) :: str_val

      write(str_val, '(I0)') int_val
      str_val = adjustl(str_val)

   end function int_to_string

   subroutine wetdep_process_load_config(this, config_manager, error_handler)
      class(WetDepProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      character(len=256) :: scheme_name
      integer :: ierr, rc

      ! Process reads directly from master YAML structure: processes.wetdep
      ! ConfigManager provides generic YAML access, process handles its own configuration

      ! Load process metadata
      call config_manager%get_string("processes/wetdep/name", this%process_name, rc, "wetdep")
      if (rc /= cc_success) this%process_name = "wetdep"  ! default

      call config_manager%get_string("processes/wetdep/version", this%process_version, rc, "1.0.0")
      if (rc /= cc_success) this%process_version = "1.0.0"  ! default

      call config_manager%get_logical("processes/wetdep/activate", this%is_active, rc, .true.)
      if (rc /= cc_success) this%is_active = .true.  ! default

      ! Load process-specific configuration directly from master YAML
      call config_manager%get_string("processes/wetdep/scheme", this%wetdep_config%scheme, rc, "jacob")
      if (rc /= cc_success) then
         call error_handler%report_error(error_invalid_config, &
            "Missing required 'scheme' in processes/wetdep configuration", rc)
         return
      end if

      ! Load diagnostic switch
      call config_manager%get_logical("processes/wetdep/diagnostics", this%wetdep_config%diagnostics, rc, .false.)
      if (rc /= cc_success) this%wetdep_config%diagnostics = .false.  ! Default

      ! Load diagnostic species list
      call config_manager%get_array("processes/wetdep/diag_species", this%wetdep_config%diagnostic_species, &
         rc, default_values=["All"])
      if (rc /= cc_success) then
         ! Default to all species if not specified
         allocate(this%wetdep_config%diagnostic_species(1))
         this%wetdep_config%diagnostic_species(1) = "All"
         this%wetdep_config%n_diagnostic_species = 1
      else
         ! Set the count based on the returned array size
         if (allocated(this%wetdep_config%diagnostic_species)) then
            this%wetdep_config%n_diagnostic_species = size(this%wetdep_config%diagnostic_species)
         else
            this%wetdep_config%n_diagnostic_species = 0
         end if
      end if

      ! Species configuration is loaded from ChemState in load_species_from_chem_state
      ! The species come from the master species YAML file (CATChem_species.yml)
      ! and are filtered by is_wetdep property


      ! Load scheme-specific configuration from master YAML
      scheme_name = trim(this%wetdep_config%scheme)
      select case (scheme_name)
       case ('jacob')
         call this%load_jacob_config(config_manager, error_handler)
       case default
         call error_handler%report_error(error_invalid_state, &
            "Unknown wetdep scheme: " // trim(scheme_name), rc)
         return
      end select

   end subroutine wetdep_process_load_config


   subroutine load_species_from_chem_state(this, chem_state, error_handler)
      use chemstate_mod, only: chemstatetype

      class(WetDepProcessConfig), intent(inout) :: this
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
      ! Dynamic mapping: is_wetdep -> nSpeciesWetdep and WetdepIndex
      this%wetdep_config%n_species = chem_state%nSpeciesWetdep

      if (this%wetdep_config%n_species <= 0) then
         call error_handler%report_error(error_invalid_state, &
            "No wetdep species found in ChemState", rc)
         return
      end if

      ! Check if WetdepIndex is allocated and has correct size
      if (.not. allocated(chem_state%WetdepIndex)) then
         call error_handler%report_error(error_invalid_state, &
            "WetdepIndex not allocated in ChemState", rc)
         return
      end if

      if (size(chem_state%WetdepIndex) < this%wetdep_config%n_species) then
         call error_handler%report_error(error_invalid_state, &
            "WetdepIndex size inconsistent with nSpeciesWetdep", rc)
         return
      end if

      ! Deallocate existing arrays if allocated
      if (allocated(this%wetdep_config%species_names)) then
         deallocate(this%wetdep_config%species_names)
      end if
      if (allocated(this%wetdep_config%species_indices)) then
         deallocate(this%wetdep_config%species_indices)
      end if

      ! Allocate arrays
      allocate(this%wetdep_config%species_names(this%wetdep_config%n_species))
      allocate(this%wetdep_config%species_indices(this%wetdep_config%n_species))

      ! Allocate species properties arrays
      allocate(this%wetdep_config%species_henry_cr(this%wetdep_config%n_species))
      allocate(this%wetdep_config%species_henry_k0(this%wetdep_config%n_species))
      allocate(this%wetdep_config%species_henry_pKa(this%wetdep_config%n_species))
      allocate(this%wetdep_config%species_is_aerosol(this%wetdep_config%n_species))
      allocate(this%wetdep_config%species_mw_g(this%wetdep_config%n_species))
      allocate(this%wetdep_config%species_radius(this%wetdep_config%n_species))
      allocate(this%wetdep_config%species_wd_LiqAndGas(this%wetdep_config%n_species))
      allocate(this%wetdep_config%species_wd_convfacI2G(this%wetdep_config%n_species))
      allocate(this%wetdep_config%species_wd_rainouteff(this%wetdep_config%n_species, 3))
      allocate(this%wetdep_config%species_wd_retfactor(this%wetdep_config%n_species))

      ! by_metadata mode: Copy indices from metadata-specific index array using dynamic mapping
      ! Dynamic mapping: is_wetdep -> WetdepIndex
      this%wetdep_config%species_indices(1:this%wetdep_config%n_species) = &
         chem_state%WetdepIndex(1:this%wetdep_config%n_species)

      ! Get species names using the indices
      do i = 1, this%wetdep_config%n_species
         if (this%wetdep_config%species_indices(i) > 0 .and. &
            this%wetdep_config%species_indices(i) <= size(chem_state%SpeciesNames)) then
            this%wetdep_config%species_names(i) = &
               trim(chem_state%SpeciesNames(this%wetdep_config%species_indices(i)))
         else
            call error_handler%report_error(error_invalid_state, &
               "Invalid species index in species index array", rc)
            return
         end if
      end do

      ! Load species properties from ChemState
      do i = 1, this%wetdep_config%n_species
         species_idx = this%wetdep_config%species_indices(i)
         this%wetdep_config%species_henry_cr(i) = chem_state%ChemSpecies(species_idx)%henry_cr
         this%wetdep_config%species_henry_k0(i) = chem_state%ChemSpecies(species_idx)%henry_k0
         this%wetdep_config%species_henry_pKa(i) = chem_state%ChemSpecies(species_idx)%henry_pKa
         this%wetdep_config%species_is_aerosol(i) = chem_state%ChemSpecies(species_idx)%is_aerosol
         this%wetdep_config%species_mw_g(i) = chem_state%ChemSpecies(species_idx)%mw_g
         this%wetdep_config%species_radius(i) = chem_state%ChemSpecies(species_idx)%radius
         this%wetdep_config%species_wd_LiqAndGas(i) = chem_state%ChemSpecies(species_idx)%wd_LiqAndGas
         this%wetdep_config%species_wd_convfacI2G(i) = chem_state%ChemSpecies(species_idx)%wd_convfacI2G
         this%wetdep_config%species_wd_rainouteff(i, :) = chem_state%ChemSpecies(species_idx)%wd_rainouteff(:)
         this%wetdep_config%species_wd_retfactor(i) = chem_state%ChemSpecies(species_idx)%wd_retfactor
      end do

   end subroutine load_species_from_chem_state


   subroutine load_jacob_config(this, config_manager, error_handler)
      class(WetDepProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: ierr, rc

      ! Load scheme parameters directly from processes/wetdep/jacob/ in master YAML
      call config_manager%get_real("processes/wetdep/jacob/scale_factor", &
         this%jacob_config%scale_factor, rc, 1.0_fp)
      if (rc /= cc_success) this%jacob_config%scale_factor = 1.0_fp
      call config_manager%get_real("processes/wetdep/jacob/radius_threshold", &
         this%jacob_config%radius_threshold, rc, 1.0_fp)
      if (rc /= cc_success) this%jacob_config%radius_threshold = 1.0_fp


   end subroutine load_jacob_config


   subroutine wetdep_process_validate(this, state_manager, error_handler)
      class(WetDepProcessConfig), intent(inout) :: this
      type(StateManagerType), intent(in) :: state_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      ! Validate main config
      call this%wetdep_config%validate(error_handler)

      ! Validate scheme-specific config
      select case (trim(this%wetdep_config%scheme))
       case ('jacob')
         call this%jacob_config%validate(error_handler)
      end select

   end subroutine wetdep_process_validate

   subroutine wetdep_process_finalize(this)
      class(WetDepProcessConfig), intent(inout) :: this

      call this%wetdep_config%finalize()
      call this%jacob_config%finalize()

   end subroutine wetdep_process_finalize

   function get_active_scheme_config(this) result(scheme_config)
      class(WetDepProcessConfig), intent(in) :: this
      class(*), allocatable :: scheme_config

      select case (trim(this%wetdep_config%scheme))
       case ('jacob')
         allocate(scheme_config, source=this%jacob_config)
       case default
         ! Return null
      end select

   end function get_active_scheme_config

   subroutine map_diagnostic_species_indices(this, error_handler)
      class(WetDepProcessConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: i, j, rc
      character(len=256) :: error_msg
      logical :: found_species

      ! Only proceed if diagnostic species are defined
      if (this%wetdep_config%n_diagnostic_species == 0) return

      ! Handle "All" case - map all available species
      if (this%wetdep_config%n_diagnostic_species == 1 .and. &
         trim(this%wetdep_config%diagnostic_species(1)) == "All") then

         ! Deallocate and reallocate for all species
         if (allocated(this%wetdep_config%diagnostic_species_id)) deallocate(this%wetdep_config%diagnostic_species_id)
         allocate(this%wetdep_config%diagnostic_species_id(this%wetdep_config%n_species))
         if (allocated(this%wetdep_config%diagnostic_species)) deallocate(this%wetdep_config%diagnostic_species)
         allocate(this%wetdep_config%diagnostic_species(this%wetdep_config%n_species))
         this%wetdep_config%n_diagnostic_species = this%wetdep_config%n_species
         this%wetdep_config%diagnostic_species = this%wetdep_config%species_names

         ! Map all species indices (1:n_species)
         do i = 1, this%wetdep_config%n_species
            this%wetdep_config%diagnostic_species_id(i) = i
         end do

         return
      end if

      ! Allocate diagnostic species indices array
      if (allocated(this%wetdep_config%diagnostic_species_id)) deallocate(this%wetdep_config%diagnostic_species_id)
      allocate(this%wetdep_config%diagnostic_species_id(this%wetdep_config%n_diagnostic_species))

      ! Map each diagnostic species name to its index in species_names
      do i = 1, this%wetdep_config%n_diagnostic_species
         found_species = .false.

         do j = 1, this%wetdep_config%n_species
            if (trim(this%wetdep_config%diagnostic_species(i)) == trim(this%wetdep_config%species_names(j))) then
               this%wetdep_config%diagnostic_species_id(i) = j
               found_species = .true.
               exit
            end if
         end do

         if (.not. found_species) then
            write(error_msg, '(A,A,A)') "Diagnostic species '", &
               trim(this%wetdep_config%diagnostic_species(i)), &
               "' not found in process species list"
            call error_handler%report_error(error_not_found, error_msg, rc)
            return
         end if
      end do

   end subroutine map_diagnostic_species_indices

end module wetdepcommon_mod
```


