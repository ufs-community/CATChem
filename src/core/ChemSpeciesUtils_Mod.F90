!> \file ChemSpeciesUtils_Mod.F90
!! \brief Utility functions for chemical species access and manipulation
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 1.0
!!
!! This module provides convenient utility functions for processes to access
!! chemical species information, indices, properties, and perform common
!! species-related operations.
!!
module ChemSpeciesUtils_Mod
   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, CC_Warning, ErrorManagerType
   use StateManager_Mod, only: StateManagerType
   use ChemState_Mod, only: ChemStateType
   use Species_Mod, only: SpeciesType

   implicit none
   private

   public :: ChemSpeciesUtilsType
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

   !> \brief Chemical species utilities helper type
   !!
   !! This type provides a collection of utility methods for working with
   !! chemical species in atmospheric processes.
   type :: ChemSpeciesUtilsType
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
   end type ChemSpeciesUtilsType

contains

   !========================================================================
   ! Standalone utility functions
   !========================================================================

   !> \brief Get the index of a single chemical species
   !!
   !! \param[in] container StateManager containing chemical state
   !! \param[in] species_name Name of the species to find
   !! \param[out] rc Return code
   !! \return Species index (> 0 if found, 0 if not found)
   function get_species_index(container, species_name, rc) result(species_idx)
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      integer :: species_idx

      type(ChemStateType), pointer :: chem_state

      rc = CC_SUCCESS
      species_idx = 0

      chem_state => container%get_chem_state_ptr()
      if (associated(chem_state)) then
         species_idx = chem_state%find_species(trim(species_name))
      else
         rc = CC_FAILURE
      end if

   end function get_species_index

   !> \brief Get indices for multiple chemical species
   !!
   !! \param[inout] container StateManager containing chemical state
   !! \param[in] species_names Array of species names to find
   !! \param[out] species_indices Array of species indices (0 if not found)
   !! \param[out] rc Return code
   subroutine get_species_indices(container, species_names, species_indices, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_names(:)
      integer, intent(out) :: species_indices(:)
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      !initialize to zero
      species_indices = 0

      do i = 1, size(species_names)
         species_indices(i) = get_species_index(container, species_names(i), rc)
         !if (rc /= CC_SUCCESS) exit !Do not exit for now
      end do

   end subroutine get_species_indices

   !> \brief Check if a chemical species exists in the mechanism
   !!
   !! \param[inout] container StateManager containing chemical state
   !! \param[in] species_name Name of the species to check
   !! \param[out] rc Return code
   !! \return True if species exists, false otherwise
   function check_species_exists(container, species_name, rc) result(exists)
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: exists

      integer :: species_idx

      species_idx = get_species_index(container, species_name, rc)
      exists = (species_idx > 0)

   end function check_species_exists

   !> \brief Get properties of a chemical species
   !!
   !! \param[in] container StateManager containing chemical state
   !! \param[in] species_name Name of the species
   !! \param[out] species Properties of the species
   !! \param[out] rc Return code
   subroutine get_species_properties(container, species_name, species, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      type(SpeciesType), intent(out) :: species
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state
      integer :: species_idx

      rc = CC_SUCCESS

      chem_state => container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         rc = CC_FAILURE
         return
      end if

      species_idx = chem_state%find_species(trim(species_name))
      if (species_idx > 0) then
         ! For now, just initialize basic species properties
         ! TODO: Implement proper species property retrieval from ChemState
         call species%init(species_name, species_name, 28.0_fp, rc)
      else
         rc = CC_FAILURE
      end if

   end subroutine get_species_properties

   !> \brief Create mapping from process species to mechanism species
   !!
   !! \param[in] container StateManager containing chemical state
   !! \param[in] process_species Array of process species names
   !! \param[out] species_mapping Array of mechanism species indices
   !! \param[out] rc Return code
   subroutine create_species_mapping(container, process_species, species_mapping, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: process_species(:)
      integer, intent(out) :: species_mapping(:)
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      character(len=256) :: message
      integer :: i

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()

      call get_species_indices(container, process_species, species_mapping, rc)
      if (rc /= CC_SUCCESS) return

      ! Check for unmapped species and warn
      do i = 1, size(species_mapping)
         if (species_mapping(i) == 0) then
            write(message, '(A,A,A)') 'Process species "', trim(process_species(i)), &
               '" not found in chemical mechanism'
            call CC_Warning(message, rc, 'create_species_mapping')
         end if
      end do

   end subroutine create_species_mapping

   !========================================================================
   ! ChemSpeciesUtilsType methods
   !========================================================================

   !> \brief Get the index of a single chemical species (type-bound method)
   function utils_get_species_index(this, container, species_name, rc) result(species_idx)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      integer :: species_idx

      species_idx = get_species_index(container, species_name, rc)

   end function utils_get_species_index

   !> \brief Get indices for multiple chemical species (type-bound method)
   subroutine utils_get_species_indices(this, container, species_names, species_indices, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_names(:)
      integer, intent(out) :: species_indices(:)
      integer, intent(out) :: rc

      call get_species_indices(container, species_names, species_indices, rc)

   end subroutine utils_get_species_indices

   !> \brief Check if a chemical species exists (type-bound method)
   function utils_check_species_exists(this, container, species_name, rc) result(exists)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: exists

      exists = check_species_exists(container, species_name, rc)

   end function utils_check_species_exists

   !> \brief Get properties of a chemical species (type-bound method)
   subroutine utils_get_species_properties(this, container, species_name, species, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      type(SpeciesType), intent(out) :: species
      integer, intent(out) :: rc

      call get_species_properties(container, species_name, species, rc)

   end subroutine utils_get_species_properties

   !> \brief Create mapping from process species to mechanism species (type-bound method)
   subroutine utils_create_species_mapping(this, container, process_species, species_mapping, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: process_species(:)
      integer, intent(out) :: species_mapping(:)
      integer, intent(out) :: rc

      call create_species_mapping(container, process_species, species_mapping, rc)

   end subroutine utils_create_species_mapping

   !> \brief Get molecular weight of a species
   !!
   !! \param[in] this ChemSpeciesUtilsType instance
   !! \param[in] container StateManager containing chemical state
   !! \param[in] species_name Name of the species
   !! \param[out] molecular_weight Molecular weight [kg/mol]
   !! \param[out] rc Return code
   subroutine utils_get_molecular_weight(this, container, species_name, molecular_weight, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      real(fp), intent(out) :: molecular_weight
      integer, intent(out) :: rc

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == CC_SUCCESS) then
         molecular_weight = species%mw_g
      else
         molecular_weight = 0.0_fp
      end if

   end subroutine utils_get_molecular_weight

   !> \brief Get list of all available species in the mechanism
   !!
   !! \param[in] this ChemSpeciesUtilsType instance
   !! \param[in] container StateManager containing chemical state
   !! \param[out] species_list List of all species names
   !! \param[out] rc Return code
   subroutine utils_get_species_list(this, container, species_list, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=32), allocatable, intent(out) :: species_list(:)
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state
      integer :: n_species, i
      character(len=64), allocatable :: species_names(:)

      rc = CC_SUCCESS

      chem_state => container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         rc = CC_FAILURE
         return
      end if

      ! Get all species names using the correct method
      species_names = chem_state%get_species()
      n_species = size(species_names)

      allocate(species_list(n_species), stat=rc)
      if (rc /= 0) then
         rc = CC_FAILURE
         return
      end if

      do i = 1, n_species
         species_list(i) = trim(species_names(i))
      end do

   end subroutine utils_get_species_list

   !> \brief Validate that a list of species are all available
   !!
   !! \param[in] this ChemSpeciesUtilsType instance
   !! \param[in] container StateManager containing chemical state
   !! \param[in] species_names List of species names to validate
   !! \param[out] all_valid True if all species are available
   !! \param[out] missing_species List of missing species (if any)
   !! \param[out] rc Return code
   subroutine utils_validate_species_list(this, container, species_names, all_valid, missing_species, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_names(:)
      logical, intent(out) :: all_valid
      character(len=32), allocatable, intent(out) :: missing_species(:)
      integer, intent(out) :: rc

      logical, allocatable :: available(:)
      integer :: i, n_missing

      rc = CC_SUCCESS
      all_valid = .true.

      allocate(available(size(species_names)), stat=rc)
      if (rc /= 0) then
         rc = CC_FAILURE
         return
      end if

      ! Check availability of each species
      do i = 1, size(species_names)
         available(i) = this%exists(container, species_names(i), rc)
         if (rc /= CC_SUCCESS) exit
         if (.not. available(i)) all_valid = .false.
      end do

      if (rc /= CC_SUCCESS) return

      ! Collect missing species
      n_missing = count(.not. available)
      if (n_missing > 0) then
         allocate(missing_species(n_missing), stat=rc)
         if (rc /= 0) then
            rc = CC_FAILURE
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

   !> \brief Get list of gas species
   !!
   !! \param[in] container StateManager containing chemical state
   !! \param[out] gas_species_list List of gas species names
   !! \param[out] rc Return code
   subroutine get_gas_species_list(container, gas_species_list, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=32), allocatable, intent(out) :: gas_species_list(:)
      integer, intent(out) :: rc

      call filter_species_by_type(container, 'gas', gas_species_list, rc)

   end subroutine get_gas_species_list

   !> \brief Get list of aerosol species
   !!
   !! \param[in] container StateManager containing chemical state
   !! \param[out] aerosol_species_list List of aerosol species names
   !! \param[out] rc Return code
   subroutine get_aerosol_species_list(container, aerosol_species_list, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=32), allocatable, intent(out) :: aerosol_species_list(:)
      integer, intent(out) :: rc

      call filter_species_by_type(container, 'aerosol', aerosol_species_list, rc)

   end subroutine get_aerosol_species_list

   !> \brief Get list of dust species
   !!
   !! \param[in] container StateManager containing chemical state
   !! \param[out] dust_species_list List of dust species names
   !! \param[out] rc Return code
   subroutine get_dust_species_list(container, dust_species_list, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=32), allocatable, intent(out) :: dust_species_list(:)
      integer, intent(out) :: rc

      call filter_species_by_type(container, 'dust', dust_species_list, rc)

   end subroutine get_dust_species_list

   !> \brief Get list of sea salt species
   !!
   !! \param[in] container StateManager containing chemical state
   !! \param[out] seasalt_species_list List of sea salt species names
   !! \param[out] rc Return code
   subroutine get_seasalt_species_list(container, seasalt_species_list, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=32), allocatable, intent(out) :: seasalt_species_list(:)
      integer, intent(out) :: rc

      call filter_species_by_type(container, 'seasalt', seasalt_species_list, rc)

   end subroutine get_seasalt_species_list

   !> \brief Get list of tracer species
   !!
   !! \param[in] container StateManager containing chemical state
   !! \param[out] tracer_species_list List of tracer species names
   !! \param[out] rc Return code
   subroutine get_tracer_species_list(container, tracer_species_list, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=32), allocatable, intent(out) :: tracer_species_list(:)
      integer, intent(out) :: rc

      call filter_species_by_type(container, 'tracer', tracer_species_list, rc)

   end subroutine get_tracer_species_list

   !> \brief Filter species by type
   !!
   !! \param[in] container StateManager containing chemical state
   !! \param[in] species_type Type of species ('gas', 'aerosol', 'dust', 'seasalt', 'tracer')
   !! \param[out] filtered_species List of species names matching the type
   !! \param[out] rc Return code
   subroutine filter_species_by_type(container, species_type, filtered_species, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_type
      character(len=32), allocatable, intent(out) :: filtered_species(:)
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state
      integer :: n_species, i, count_matching
      character(len=32), allocatable :: temp_list(:)

      rc = CC_SUCCESS

      ! Suppress unused variable warning for placeholder implementation
      if (len_trim(species_type) == 0) continue

      chem_state => container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         rc = CC_FAILURE
         return
      end if

      n_species = chem_state%get_num_species()
      allocate(temp_list(n_species), stat=rc)
      if (rc /= 0) then
         rc = CC_FAILURE
         return
      end if

      count_matching = 0
      do i = 1, n_species
         ! TODO: Replace with proper species property access when available
         ! For now, just add all species to avoid compilation errors
         count_matching = count_matching + 1
         write(temp_list(count_matching), '(A,I0)') "Species_", i
      end do

      if (rc /= CC_SUCCESS) then
         deallocate(temp_list)
         return
      end if

      ! Allocate final list with correct size
      allocate(filtered_species(count_matching), stat=rc)
      if (rc /= 0) then
         rc = CC_FAILURE
         deallocate(temp_list)
         return
      end if

      filtered_species(1:count_matching) = temp_list(1:count_matching)
      deallocate(temp_list)

   end subroutine filter_species_by_type

   !> \brief Get concentration units for a species
   !!
   !! \param[in] container StateManager containing chemical state
   !! \param[in] species_name Name of the species
   !! \param[out] units Units string ('v/v', 'kg/kg', etc.)
   !! \param[out] rc Return code
   subroutine get_species_concentration_units(container, species_name, units, rc)
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      character(len=10), intent(out) :: units
      integer, intent(out) :: rc

      type(SpeciesType) :: species

      call get_species_properties(container, species_name, species, rc)
      if (rc == CC_SUCCESS) then
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

   !> \brief Parse a comma-separated species list
   !!
   !! \param[in] species_string Comma-separated species names
   !! \param[out] species_array Array of parsed species names
   !! \param[out] n_species Number of species parsed
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

   !> \brief Get density of a species (type-bound method)
   subroutine utils_get_density(this, container, species_name, density, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      real(fp), intent(out) :: density
      integer, intent(out) :: rc

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == CC_SUCCESS) then
         density = species%density
      else
         density = 0.0_fp
      end if

   end subroutine utils_get_density

   !> \brief Get radius of a species (type-bound method)
   subroutine utils_get_radius(this, container, species_name, radius, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      real(fp), intent(out) :: radius
      integer, intent(out) :: rc

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == CC_SUCCESS) then
         radius = species%radius
      else
         radius = 0.0_fp
      end if

   end subroutine utils_get_radius

   !> \brief Check if species is a gas (type-bound method)
   function utils_is_gas(this, container, species_name, rc) result(is_gas)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: is_gas

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == CC_SUCCESS) then
         is_gas = species%is_gas
      else
         is_gas = .false.
      end if

   end function utils_is_gas

   !> \brief Check if species is an aerosol (type-bound method)
   function utils_is_aerosol(this, container, species_name, rc) result(is_aerosol)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: is_aerosol

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == CC_SUCCESS) then
         is_aerosol = species%is_aerosol
      else
         is_aerosol = .false.
      end if

   end function utils_is_aerosol

   !> \brief Check if species is dust (type-bound method)
   function utils_is_dust(this, container, species_name, rc) result(is_dust)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: is_dust

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == CC_SUCCESS) then
         is_dust = species%is_dust
      else
         is_dust = .false.
      end if

   end function utils_is_dust

   !> \brief Check if species is sea salt (type-bound method)
   function utils_is_seasalt(this, container, species_name, rc) result(is_seasalt)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: is_seasalt

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == CC_SUCCESS) then
         is_seasalt = species%is_seasalt
      else
         is_seasalt = .false.
      end if

   end function utils_is_seasalt

   !> \brief Check if species is a tracer (type-bound method)
   function utils_is_tracer(this, container, species_name, rc) result(is_tracer)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: is_tracer

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == CC_SUCCESS) then
         is_tracer = species%is_tracer
      else
         is_tracer = .false.
      end if

   end function utils_is_tracer

   !> \brief Check if species undergoes dry deposition (type-bound method)
   function utils_undergoes_drydep(this, container, species_name, rc) result(undergoes_drydep)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: undergoes_drydep

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == CC_SUCCESS) then
         undergoes_drydep = species%is_drydep
      else
         undergoes_drydep = .false.
      end if

   end function utils_undergoes_drydep

   !> \brief Check if species undergoes photolysis (type-bound method)
   function utils_undergoes_photolysis(this, container, species_name, rc) result(undergoes_photolysis)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      logical :: undergoes_photolysis

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == CC_SUCCESS) then
         undergoes_photolysis = species%is_photolysis
      else
         undergoes_photolysis = .false.
      end if

   end function utils_undergoes_photolysis

   !> \brief Get background concentration of a species (type-bound method)
   subroutine utils_get_background_concentration(this, container, species_name, background_conc, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_name
      real(fp), intent(out) :: background_conc
      integer, intent(out) :: rc

      type(SpeciesType) :: species

      call this%get_properties(container, species_name, species, rc)
      if (rc == CC_SUCCESS) then
         background_conc = species%BackgroundVV
      else
         background_conc = 0.0_fp
      end if

   end subroutine utils_get_background_concentration

   !> \brief Filter species by type (type-bound method)
   subroutine utils_filter_species_by_type(this, container, species_type, filtered_species, rc)
      class(ChemSpeciesUtilsType), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      character(len=*), intent(in) :: species_type
      character(len=32), allocatable, intent(out) :: filtered_species(:)
      integer, intent(out) :: rc

      call filter_species_by_type(container, species_type, filtered_species, rc)

   end subroutine utils_filter_species_by_type
end module ChemSpeciesUtils_Mod
