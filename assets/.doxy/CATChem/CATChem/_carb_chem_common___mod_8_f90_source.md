

# File CarbChemCommon\_Mod.F90

[**File List**](files.md) **>** [**carbchem**](dir_5dbdd03f815becc4c35f94e0692b4e09.md) **>** [**CarbChemCommon\_Mod.F90**](_carb_chem_common___mod_8_f90.md)

[Go to the documentation of this file](_carb_chem_common___mod_8_f90.md)


```Fortran


module carbchemcommon_mod

   use precision_mod, only: fp
   ! use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure, cc_error, cc_warning, errormanagertype, &
      error_invalid_config, error_invalid_state, error_not_found
   use configmanager_mod, only: configmanagertype  ! ConfigManager integration
   use statemanager_mod, only: statemanagertype  ! Add StateManager integration

   implicit none
   private

   ! Export types
   public :: carbchemprocessconfig  ! New unified process config
   public :: carbchemconfig
   public :: carbchemschemegocartconfig

   ! Export utility functions
   public :: int_to_string

   type :: carbchemconfig

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
      integer, allocatable :: species_indices(:)  ! Indices of carbchem species in ChemState



      ! Species properties
      real(fp), allocatable :: species_t_chem_loss(:)      ! t_chem_loss for each species

      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   contains
      procedure, public :: validate => validate_carbchem_config
      procedure, public :: finalize => finalize_carbchem_config
      procedure, public :: print_summary => print_carbchem_config_summary
   end type carbchemconfig

   type :: carbchemschemegocartconfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'gocart'
      character(len=256) :: description = 'GOCART carbon species chemical production and loss scheme'
      character(len=64) :: author = 'Wei Li'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .true.  ! Full column processing

      ! Scheme parameters
      real(fp) :: time_days_hydrophobic_to_hydrophilic = 2.5  ! Rate of conversion of hydrophobic to hydrophilic [days]

      ! Required meteorological fields
      integer :: n_required_met_fields = 4
      character(len=32) :: required_met_fields(4)

   contains
      procedure, public :: validate => validate_gocart_config
      procedure, public :: finalize => finalize_gocart_config
   end type carbchemschemegocartconfig

   ! gocart scheme uses local variables only - no persistent state type needed


   type :: carbchemprocessconfig

      ! Process metadata
      character(len=64) :: process_name = 'carbchem'
      character(len=16) :: process_version = '1.0.0'
      logical :: is_active = .true.

      ! Process-specific configuration (delegate to CarbChemConfig)
      type(CarbChemConfig) :: carbchem_config

      ! Scheme configurations
      type(CarbChemSchemeGOCARTConfig) :: gocart_config


   contains
      procedure, public :: load_from_config => carbchem_process_load_config
      procedure, public :: load_species_from_list => load_species_from_list
      procedure, public :: validate => carbchem_process_validate
      procedure, public :: finalize => carbchem_process_finalize
      procedure, public :: get_active_scheme_config => get_active_scheme_config
      procedure, public :: load_gocart_config
      procedure, public :: map_diagnostic_species_indices
   end type carbchemprocessconfig

contains

   subroutine validate_carbchem_config(this, error_handler)
      class(CarbChemConfig), intent(inout) :: this
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

   end subroutine validate_carbchem_config

   subroutine print_carbchem_config_summary(this)
      class(CarbChemConfig), intent(in) :: this

      write(*, '(A)') "=== CarbChem Process Configuration ==="
      write(*, '(A,A)') "  Active scheme: ", trim(this%scheme)
      write(*, '(A,I0)') "  Number of species: ", this%n_species
      write(*, '(A,F0.1,A)') "  Minimum time step: ", this%dt_min, " s"
      write(*, '(A,F0.1,A)') "  Maximum time step: ", this%dt_max, " s"
      write(*, '(A,L1)') "  Output diagnostics: ", this%output_diagnostics
      write(*, '(A)') "============================================="

   end subroutine print_carbchem_config_summary

   subroutine finalize_carbchem_config(this)
      class(CarbChemConfig), intent(inout) :: this

      ! Deallocate species names array
      if (allocated(this%species_names)) then
         deallocate(this%species_names)
      end if

      ! Deallocate species indices array
      if (allocated(this%species_indices)) then
         deallocate(this%species_indices)
      end if

      ! Deallocate species properties arrays
      if (allocated(this%species_t_chem_loss)) then
         deallocate(this%species_t_chem_loss)
      end if


      ! Deallocate diagnostic species array
      if (allocated(this%diagnostic_species)) then
         deallocate(this%diagnostic_species)
      end if

      ! Deallocate diagnostic species indices array
      if (allocated(this%diagnostic_species_id)) then
         deallocate(this%diagnostic_species_id)
      end if

   end subroutine finalize_carbchem_config

   subroutine validate_gocart_config(this, error_handler)
      class(CarbChemSchemeGOCARTConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_gocart_config

   subroutine finalize_gocart_config(this)
      class(CarbChemSchemeGOCARTConfig), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_gocart_config



   function int_to_string(int_val) result(str_val)
      integer, intent(in) :: int_val
      character(len=32) :: str_val

      write(str_val, '(I0)') int_val
      str_val = adjustl(str_val)

   end function int_to_string

   subroutine carbchem_process_load_config(this, config_manager, error_handler)
      class(CarbChemProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      character(len=256) :: scheme_name
      integer :: ierr, rc

      ! Process reads directly from master YAML structure: processes.carbchem
      ! ConfigManager provides generic YAML access, process handles its own configuration

      ! Load process metadata
      call config_manager%get_string("processes/carbchem/name", this%process_name, rc, "carbchem")
      if (rc /= cc_success) this%process_name = "carbchem"  ! default

      call config_manager%get_string("processes/carbchem/version", this%process_version, rc, "1.0.0")
      if (rc /= cc_success) this%process_version = "1.0.0"  ! default

      call config_manager%get_logical("processes/carbchem/activate", this%is_active, rc, .true.)
      if (rc /= cc_success) this%is_active = .true.  ! default

      ! Load process-specific configuration directly from master YAML
      call config_manager%get_string("processes/carbchem/scheme", this%carbchem_config%scheme, rc, "gocart")
      if (rc /= cc_success) then
         call error_handler%report_error(error_invalid_config, &
            "Missing required 'scheme' in processes/carbchem configuration", rc)
         return
      end if

      ! Load diagnostic switch
      call config_manager%get_logical("processes/carbchem/diagnostics", this%carbchem_config%diagnostics, rc, .false.)
      if (rc /= cc_success) this%carbchem_config%diagnostics = .false.  ! Default

      ! Load diagnostic species list
      call config_manager%get_array("processes/carbchem/diag_species", this%carbchem_config%diagnostic_species, &
         rc, default_values=["All"])
      if (rc /= cc_success) then
         ! Default to all species if not specified
         allocate(this%carbchem_config%diagnostic_species(1))
         this%carbchem_config%diagnostic_species(1) = "All"
         this%carbchem_config%n_diagnostic_species = 1
      else
         ! Set the count based on the returned array size
         if (allocated(this%carbchem_config%diagnostic_species)) then
            this%carbchem_config%n_diagnostic_species = size(this%carbchem_config%diagnostic_species)
         else
            this%carbchem_config%n_diagnostic_species = 0
         end if
      end if

      ! Species configuration is loaded from ChemState in load_species_from_chem_state
      ! The species come from the master species YAML file (CATChem_species.yml)
      ! and are filtered by is_carbchem property


      ! Load scheme-specific configuration from master YAML
      scheme_name = trim(this%carbchem_config%scheme)
      select case (scheme_name)
       case ('gocart')
         call this%load_gocart_config(config_manager, error_handler)
       case default
         call error_handler%report_error(error_invalid_state, &
            "Unknown carbchem scheme: " // trim(scheme_name), rc)
         return
      end select

   end subroutine carbchem_process_load_config

   subroutine load_species_from_list(this, chem_state, error_handler, rc)
      use chemstate_mod, only: chemstatetype

      class(CarbChemProcessConfig), intent(inout) :: this
      type(ChemStateType), pointer, intent(in) :: chem_state
      type(ErrorManagerType), intent(inout) :: error_handler
      integer, intent(inout) :: rc

      integer :: i, species_idx
      character(len=32) :: species_name
      logical :: found

      ! Explicit species list from configuration
      character(len=32), parameter :: EXPLICIT_SPECIES(4) = [ &
         'oc1                             ', &
         'oc2                             ', &
         'bc1                             ', &
         'bc2                             ' ]

      if (.not. associated(chem_state)) then
         call error_handler%report_error(error_invalid_state, &
            "ChemState not associated in load_species_from_list", rc)
         return
      end if

      ! Set number of species from explicit list
      this%carbchem_config%n_species = 4

      ! Deallocate existing arrays if allocated
      if (allocated(this%carbchem_config%species_names)) then
         deallocate(this%carbchem_config%species_names)
      end if
      if (allocated(this%carbchem_config%species_indices)) then
         deallocate(this%carbchem_config%species_indices)
      end if
      if (allocated(this%carbchem_config%species_t_chem_loss)) then
         deallocate(this%carbchem_config%species_t_chem_loss)
      end if

      ! Allocate arrays
      allocate(this%carbchem_config%species_names(this%carbchem_config%n_species))
      allocate(this%carbchem_config%species_indices(this%carbchem_config%n_species))

      ! Allocate species properties arrays
      allocate(this%carbchem_config%species_t_chem_loss(this%carbchem_config%n_species))

      ! Find indices for each explicit species
      do i = 1, this%carbchem_config%n_species
         species_name = trim(explicit_species(i))
         found = .false.

         ! Search for species in ChemState
         do species_idx = 1, size(chem_state%SpeciesNames)
            if (trim(chem_state%SpeciesNames(species_idx)) == species_name) then
               this%carbchem_config%species_names(i) = species_name
               this%carbchem_config%species_indices(i) = species_idx
               found = .true.
               exit
            end if
         end do

         if (.not. found) then
            call error_handler%report_error(error_invalid_state, &
               "Species not found in ChemState: " // trim(species_name), rc)
            return
         end if
      end do

      ! Load species properties from ChemState using found indices
      do i = 1, this%carbchem_config%n_species
         species_idx = this%carbchem_config%species_indices(i)
         this%carbchem_config%species_t_chem_loss(i) = chem_state%ChemSpecies(species_idx)%t_chem_loss
      end do

   end subroutine load_species_from_list



   subroutine load_gocart_config(this, config_manager, error_handler)
      class(CarbChemProcessConfig), intent(inout) :: this
      type(ConfigManagerType), intent(inout) :: config_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: ierr, rc

      ! Load scheme parameters directly from processes/carbchem/gocart/ in master YAML
      call config_manager%get_real("processes/carbchem/gocart/time_days_hydrophobic_to_hydrophilic", &
         this%gocart_config%time_days_hydrophobic_to_hydrophilic, rc, 2.5_fp)
      if (rc /= cc_success) this%gocart_config%time_days_hydrophobic_to_hydrophilic = 2.5_fp


   end subroutine load_gocart_config


   subroutine carbchem_process_validate(this, state_manager, error_handler)
      class(CarbChemProcessConfig), intent(inout) :: this
      type(StateManagerType), intent(in) :: state_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      ! Validate main config
      call this%carbchem_config%validate(error_handler)

      ! Validate scheme-specific config
      select case (trim(this%carbchem_config%scheme))
       case ('gocart')
         call this%gocart_config%validate(error_handler)
      end select

   end subroutine carbchem_process_validate

   subroutine carbchem_process_finalize(this)
      class(CarbChemProcessConfig), intent(inout) :: this


      call this%carbchem_config%finalize()
      call this%gocart_config%finalize()

   end subroutine carbchem_process_finalize


   function get_active_scheme_config(this) result(scheme_config)
      class(CarbChemProcessConfig), intent(in) :: this
      class(*), allocatable :: scheme_config

      select case (trim(this%carbchem_config%scheme))
       case ('gocart')
         allocate(scheme_config, source=this%gocart_config)
       case default
         ! Return null
      end select

   end function get_active_scheme_config

   subroutine map_diagnostic_species_indices(this, error_handler)
      class(CarbChemProcessConfig), intent(inout) :: this
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: i, j, rc
      character(len=256) :: error_msg
      logical :: found_species

      ! Only proceed if diagnostic species are defined
      if (this%carbchem_config%n_diagnostic_species == 0) return

      ! Handle "All" case - map all available species
      if (this%carbchem_config%n_diagnostic_species == 1 .and. &
         trim(this%carbchem_config%diagnostic_species(1)) == "All") then

         ! Deallocate and reallocate for all species
         if (allocated(this%carbchem_config%diagnostic_species_id)) deallocate(this%carbchem_config%diagnostic_species_id)
         allocate(this%carbchem_config%diagnostic_species_id(this%carbchem_config%n_species))
         if (allocated(this%carbchem_config%diagnostic_species)) deallocate(this%carbchem_config%diagnostic_species)
         allocate(this%carbchem_config%diagnostic_species(this%carbchem_config%n_species))
         this%carbchem_config%n_diagnostic_species = this%carbchem_config%n_species
         this%carbchem_config%diagnostic_species = this%carbchem_config%species_names

         ! Map all species indices (1:n_species)
         do i = 1, this%carbchem_config%n_species
            this%carbchem_config%diagnostic_species_id(i) = i
         end do

         return
      end if

      ! Allocate diagnostic species indices array
      if (allocated(this%carbchem_config%diagnostic_species_id)) deallocate(this%carbchem_config%diagnostic_species_id)
      allocate(this%carbchem_config%diagnostic_species_id(this%carbchem_config%n_diagnostic_species))

      ! Map each diagnostic species name to its index in species_names
      do i = 1, this%carbchem_config%n_diagnostic_species
         found_species = .false.

         do j = 1, this%carbchem_config%n_species
            if (trim(this%carbchem_config%diagnostic_species(i)) == trim(this%carbchem_config%species_names(j))) then
               this%carbchem_config%diagnostic_species_id(i) = j
               found_species = .true.
               exit
            end if
         end do

         if (.not. found_species) then
            write(error_msg, '(A,A,A)') "Diagnostic species '", &
               trim(this%carbchem_config%diagnostic_species(i)), &
               "' not found in process species list"
            call error_handler%report_error(error_not_found, error_msg, rc)
            !return !do not return and the diagnostics for this unspecified species will be zero in the output
         end if
      end do

   end subroutine map_diagnostic_species_indices

end module carbchemcommon_mod
```


