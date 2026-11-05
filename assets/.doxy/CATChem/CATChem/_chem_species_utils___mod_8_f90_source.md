

# File ChemSpeciesUtils\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**ChemSpeciesUtils\_Mod.F90**](_chem_species_utils___mod_8_f90.md)

[Go to the documentation of this file](_chem_species_utils___mod_8_f90.md)


```Fortran

module chemspeciesutils_mod
   use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure, cc_warning, errormanagertype
   use statemanager_mod, only: statemanagertype
   use chemstate_mod, only: chemstatetype
   use species_mod, only: speciestype

   implicit none
   private

   public :: chemspeciesutilstype
   public :: get_species_index
   public :: get_species_indices
   public :: check_species_exists
   public :: get_species_properties
   public :: create_species_mapping

   ! Additional utilities
   public :: get_gas_species_list
   public :: get_aerosol_species_list
   public :: get_dust_species_list
   public :: get_seasalt_species_list
   public :: get_tracer_species_list
   public :: filter_species_by_type
   public :: get_species_concentration_units
   public :: parse_species_list

   type :: chemspeciesutilstype
   contains
      procedure :: get_index => utils_get_species_index
      procedure :: get_indices => utils_get_species_indices
      procedure :: exists => utils_check_species_exists
      procedure :: get_properties => utils_get_species_properties
      procedure :: create_mapping => utils_create_species_mapping
      procedure :: get_molecular_weight => utils_get_molecular_weight
      procedure :: get_species_list => utils_get_species_list
      procedure :: validate_species_list => utils_validate_species_list

      ! Additional type-bound methods
      procedure :: get_density => utils_get_density
      procedure :: get_radius => utils_get_radius
      procedure :: is_gas => utils_is_gas
      procedure :: is_aerosol => utils_is_aerosol
      procedure :: is_dust => utils_is_dust
      procedure :: is_seasalt => utils_is_seasalt
      procedure :: is_tracer => utils_is_tracer
      procedure :: undergoes_drydep => utils_undergoes_drydep
      procedure :: undergoes_photolysis => utils_undergoes_photolysis
      procedure :: get_background_concentration => utils_get_background_concentration
      procedure :: filter_by_type => utils_filter_species_by_type
   end type chemspeciesutilstype

contains

   !========================================================================
   ! Standalone utility functions
   !========================================================================

   function get_species_index(container, species_name, rc) result(species_idx)
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      integer :: species_idx

      type(ChemStateType), pointer :: chem_state

      rc = cc_success
      species_idx = 0

      chem_state => container%get_chem_state_ptr()
      if (associated(chem_state)) then
         species_idx = chem_state%find_species(trim(species_name))
      else
         rc = cc_failure
      end if

   end function get_species_index

   subroutine get_species_indices(container, species_names, species_indices, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_names(:)
      integer, intent(out) :: species_indices(:)
      integer, intent(out) :: rc

      integer :: i

      rc = cc_success

      !initialize to zero
      species_indices = 0

      do i = 1, size(species_names)
         species_indices(i) = get_species_index(container, species_names(i), rc)
         !if (rc /= CC_SUCCESS) exit !Do not exit for now
      end do

   end subroutine get_species_indices

   function check_species_exists(container, species_name, rc) result(exists)
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: exists

      integer :: species_idx

      species_idx = get_species_index(container, species_name, rc)
      exists = (species_idx > 0)

   end function check_species_exists

   subroutine get_species_properties(container, species_name, species, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      type(SpeciesType), intent(out) :: species
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state
      integer :: species_idx

      rc = cc_success

      chem_state => container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         rc = cc_failure
         return
      end if

      species_idx = chem_state%find_species(trim(species_name))
      if (species_idx > 0) then
         ! For now, just initialize basic species properties
         ! TODO: Implement proper species property retrieval from ChemState
         call species%init(species_name, species_name, 28.0_fp, rc)
      else
         rc = cc_failure
      end if

   end subroutine get_species_properties

   subroutine create_species_mapping(container, process_species, species_mapping, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: process_species(:)
      integer, intent(out) :: species_mapping(:)
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      character(len=256) :: message
      integer :: i

      rc = cc_success
      error_mgr => container%get_error_manager()

      call get_species_indices(container, process_species, species_mapping, rc)
      if (rc /= cc_success) return

      ! Check for unmapped species and warn
      do i = 1, size(species_mapping)
         if (species_mapping(i) == 0) then
            write(message, '(A,A,A)') 'Process species "', trim(process_species(i)), &
               '" not found in chemical mechanism'
            call cc_warning(message, rc, 'create_species_mapping')
         end if
      end do

   end subroutine create_species_mapping

   !========================================================================
   ! ChemSpeciesUtilsType methods
   !========================================================================

   function utils_get_species_index(this, container, species_name, rc) result(species_idx)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      integer :: species_idx

      species_idx = get_species_index(container, species_name, rc)

   end function utils_get_species_index

   subroutine utils_get_species_indices(this, container, species_names, species_indices, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_names(:)
      integer, intent(out) :: species_indices(:)
      integer, intent(out) :: rc

      call get_species_indices(container, species_names, species_indices, rc)

   end subroutine utils_get_species_indices

   function utils_check_species_exists(this, container, species_name, rc) result(exists)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: exists

      exists = check_species_exists(container, species_name, rc)

   end function utils_check_species_exists

   subroutine utils_get_species_properties(this, container, species_name, species, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      type(SpeciesType), intent(out) :: species
      integer, intent(out) :: rc

      call get_species_properties(container, species_name, species, rc)

   end subroutine utils_get_species_properties

   subroutine utils_create_species_mapping(this, container, process_species, species_mapping, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: process_species(:)
      integer, intent(out) :: species_mapping(:)
      integer, intent(out) :: rc

      call create_species_mapping(container, process_species, species_mapping, rc)

   end subroutine utils_create_species_mapping

   subroutine utils_get_molecular_weight(this, container, species_name, molecular_weight, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      real(fp), intent(out) :: molecular_weight
      integer, intent(out) :: rc

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == cc_success) then
         molecular_weight = species%mw_g
      else
         molecular_weight = 0.0_fp
      end if

   end subroutine utils_get_molecular_weight

   subroutine utils_get_species_list(this, container, species_list, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=32), allocatable, intent(out) :: species_list(:)
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state
      integer :: n_species, i
      character(len=64), allocatable :: species_names(:)

      rc = cc_success

      chem_state => container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         rc = cc_failure
         return
      end if

      ! Get all species names using the correct method
      species_names = chem_state%get_species()
      n_species = size(species_names)

      allocate(species_list(n_species), stat=rc)
      if (rc /= 0) then
         rc = cc_failure
         return
      end if

      do i = 1, n_species
         species_list(i) = trim(species_names(i))
      end do

   end subroutine utils_get_species_list

   subroutine utils_validate_species_list(this, container, species_names, all_valid, missing_species, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_names(:)
      logical, intent(out) :: all_valid
      character(len=32), allocatable, intent(out) :: missing_species(:)
      integer, intent(out) :: rc

      logical, allocatable :: available(:)
      integer :: i, n_missing

      rc = cc_success
      all_valid = .true.

      allocate(available(size(species_names)), stat=rc)
      if (rc /= 0) then
         rc = cc_failure
         return
      end if

      ! Check availability of each species
      do i = 1, size(species_names)
         available(i) = this%exists(container, species_names(i), rc)
         if (rc /= cc_success) exit
         if (.not. available(i)) all_valid = .false.
      end do

      if (rc /= cc_success) return

      ! Collect missing species
      n_missing = count(.not. available)
      if (n_missing > 0) then
         allocate(missing_species(n_missing), stat=rc)
         if (rc /= 0) then
            rc = cc_failure
            return
         end if

         n_missing = 0
         do i = 1, size(species_names)
            if (.not. available(i)) then
               n_missing = n_missing + 1
               missing_species(n_missing) = trim(species_names(i))
            end if
         end do
      else
         allocate(missing_species(0))
      end if

   end subroutine utils_validate_species_list

   !========================================================================
   ! Additional utility functions for specific species types
   !========================================================================

   subroutine get_gas_species_list(container, gas_species_list, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=32), allocatable, intent(out) :: gas_species_list(:)
      integer, intent(out) :: rc

      call filter_species_by_type(container, 'gas', gas_species_list, rc)

   end subroutine get_gas_species_list

   subroutine get_aerosol_species_list(container, aerosol_species_list, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=32), allocatable, intent(out) :: aerosol_species_list(:)
      integer, intent(out) :: rc

      call filter_species_by_type(container, 'aerosol', aerosol_species_list, rc)

   end subroutine get_aerosol_species_list

   subroutine get_dust_species_list(container, dust_species_list, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=32), allocatable, intent(out) :: dust_species_list(:)
      integer, intent(out) :: rc

      call filter_species_by_type(container, 'dust', dust_species_list, rc)

   end subroutine get_dust_species_list

   subroutine get_seasalt_species_list(container, seasalt_species_list, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=32), allocatable, intent(out) :: seasalt_species_list(:)
      integer, intent(out) :: rc

      call filter_species_by_type(container, 'seasalt', seasalt_species_list, rc)

   end subroutine get_seasalt_species_list

   subroutine get_tracer_species_list(container, tracer_species_list, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=32), allocatable, intent(out) :: tracer_species_list(:)
      integer, intent(out) :: rc

      call filter_species_by_type(container, 'tracer', tracer_species_list, rc)

   end subroutine get_tracer_species_list

   subroutine filter_species_by_type(container, species_type, filtered_species, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_type
      character(len=32), allocatable, intent(out) :: filtered_species(:)
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state
      integer :: n_species, i, count_matching
      character(len=32), allocatable :: temp_list(:)

      rc = cc_success

      ! Suppress unused variable warning for placeholder implementation
      if (len_trim(species_type) == 0) continue

      chem_state => container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         rc = cc_failure
         return
      end if

      n_species = chem_state%get_num_species()
      allocate(temp_list(n_species), stat=rc)
      if (rc /= 0) then
         rc = cc_failure
         return
      end if

      count_matching = 0
      do i = 1, n_species
         ! TODO: Replace with proper species property access when available
         ! For now, just add all species to avoid compilation errors
         count_matching = count_matching + 1
         write(temp_list(count_matching), '(A,I0)') "Species_", i
      end do

      if (rc /= cc_success) then
         deallocate(temp_list)
         return
      end if

      ! Allocate final list with correct size
      allocate(filtered_species(count_matching), stat=rc)
      if (rc /= 0) then
         rc = cc_failure
         deallocate(temp_list)
         return
      end if

      filtered_species(1:count_matching) = temp_list(1:count_matching)
      deallocate(temp_list)

   end subroutine filter_species_by_type

   subroutine get_species_concentration_units(container, species_name, units, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      character(len=10), intent(out) :: units
      integer, intent(out) :: rc

      type(SpeciesType) :: species

      call get_species_properties(container, species_name, species, rc)
      if (rc == cc_success) then
         if (species%is_gas) then
            units = 'v/v'
         else if (species%is_aerosol .or. species%is_dust .or. species%is_seasalt) then
            units = 'kg/kg'
         else
            units = 'v/v'  ! Default for tracers
         end if
      else
         units = 'unknown'
      end if

   end subroutine get_species_concentration_units

   subroutine parse_species_list(species_string, species_array, n_species)
      character(len=*), intent(in) :: species_string
      character(len=32), intent(out) :: species_array(:)
      integer, intent(out) :: n_species

      character(len=len(species_string)) :: work_string
      integer :: comma_pos, start_pos
      character(len=32) :: temp_name

      work_string = trim(adjustl(species_string))
      n_species = 0
      start_pos = 1

      do while (n_species < size(species_array) .and. start_pos <= len_trim(work_string))
         comma_pos = index(work_string(start_pos:), ',')

         if (comma_pos == 0) then
            ! Last species in the list
            temp_name = trim(adjustl(work_string(start_pos:)))
            if (len_trim(temp_name) > 0) then
               n_species = n_species + 1
               species_array(n_species) = temp_name
            end if
            exit
         else
            ! Extract species name before comma
            temp_name = trim(adjustl(work_string(start_pos:start_pos+comma_pos-2)))
            if (len_trim(temp_name) > 0) then
               n_species = n_species + 1
               species_array(n_species) = temp_name
            end if
            start_pos = start_pos + comma_pos
         end if
      end do

   end subroutine parse_species_list

   !========================================================================
   ! Additional type-bound methods for ChemSpeciesUtilsType
   !========================================================================

   subroutine utils_get_density(this, container, species_name, density, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      real(fp), intent(out) :: density
      integer, intent(out) :: rc

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == cc_success) then
         density = species%density
      else
         density = 0.0_fp
      end if

   end subroutine utils_get_density

   subroutine utils_get_radius(this, container, species_name, radius, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      real(fp), intent(out) :: radius
      integer, intent(out) :: rc

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == cc_success) then
         radius = species%radius
      else
         radius = 0.0_fp
      end if

   end subroutine utils_get_radius

   function utils_is_gas(this, container, species_name, rc) result(is_gas)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: is_gas

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == cc_success) then
         is_gas = species%is_gas
      else
         is_gas = .false.
      end if

   end function utils_is_gas

   function utils_is_aerosol(this, container, species_name, rc) result(is_aerosol)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: is_aerosol

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == cc_success) then
         is_aerosol = species%is_aerosol
      else
         is_aerosol = .false.
      end if

   end function utils_is_aerosol

   function utils_is_dust(this, container, species_name, rc) result(is_dust)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: is_dust

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == cc_success) then
         is_dust = species%is_dust
      else
         is_dust = .false.
      end if

   end function utils_is_dust

   function utils_is_seasalt(this, container, species_name, rc) result(is_seasalt)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: is_seasalt

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == cc_success) then
         is_seasalt = species%is_seasalt
      else
         is_seasalt = .false.
      end if

   end function utils_is_seasalt

   function utils_is_tracer(this, container, species_name, rc) result(is_tracer)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: is_tracer

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == cc_success) then
         is_tracer = species%is_tracer
      else
         is_tracer = .false.
      end if

   end function utils_is_tracer

   function utils_undergoes_drydep(this, container, species_name, rc) result(undergoes_drydep)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: undergoes_drydep

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == cc_success) then
         undergoes_drydep = species%is_drydep
      else
         undergoes_drydep = .false.
      end if

   end function utils_undergoes_drydep

   function utils_undergoes_photolysis(this, container, species_name, rc) result(undergoes_photolysis)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: undergoes_photolysis

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == cc_success) then
         undergoes_photolysis = species%is_photolysis
      else
         undergoes_photolysis = .false.
      end if

   end function utils_undergoes_photolysis

   subroutine utils_get_background_concentration(this, container, species_name, background_conc, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      real(fp), intent(out) :: background_conc
      integer, intent(out) :: rc

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == cc_success) then
         background_conc = species%BackgroundVV
      else
         background_conc = 0.0_fp
      end if

   end subroutine utils_get_background_concentration

   subroutine utils_filter_species_by_type(this, container, species_type, filtered_species, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_type
      character(len=32), allocatable, intent(out) :: filtered_species(:)
      integer, intent(out) :: rc

      call filter_species_by_type(container, species_type, filtered_species, rc)

   end subroutine utils_filter_species_by_type
end module chemspeciesutils_mod
```


