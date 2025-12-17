!> \file species_mod.F90
!! \brief Modern species definition and management for CATChem
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module contains the enhanced SpeciesType and related routines
!! for managing chemical species with integrated error handling and
!! comprehensive validation.
!!
!! \details
!! **Enhanced features in v2.0:**
!! - Integrated error handling and validation
!! - Enhanced species property management
!! - Better memory management and cleanup
!! - Comprehensive species database operations
!! - Thread-safe species operations
!!

module species_mod
   use precision_mod
   use error_mod
   use, intrinsic :: ieee_arithmetic
   implicit none
   private

   public :: SpeciesType
   public :: SpeciesManagerType
   public :: validate_species
   public :: find_species_by_name
   public :: create_species_database

   !> Derived type for chemical species
   !!
   !! This type contains all properties and data for a chemical species
   !! in the CATChem atmospheric chemistry model, including physical
   !! properties, classification flags, and concentration data.
   !!
   !! @param long_name Long descriptive name for species (NetCDF attribute)
   !! @param short_name Short identifier name for species
   !! @param description Detailed description of the species
   !! @param is_gas Logical flag: true if species is gaseous
   !! @param is_aerosol Logical flag: true if species is an aerosol
   !! @param is_tracer Logical flag: true if species is a passive tracer
   !! @param is_advected Logical flag: true if species undergoes advection
   !! @param is_drydep Logical flag: true if species undergoes dry deposition
   !! @param is_photolysis Logical flag: true if species undergoes photolysis
   !! @param is_gocart_aero Logical flag: true if species is a GOCART aerosol
   !! @param is_dust Logical flag: true if species is dust
   !! @param is_seasalt Logical flag: true if species is sea salt
   !! @param mw_g Gaseous molecular weight [g/mol]
   !! @param density Particle density [kg/m³]
   !! @param radius Mean molecular diameter [m]
   !! @param lower_radius Lower radius bound [m]
   !! @param upper_radius Upper radius bound [m]
   !! @param viscosity Kinematic viscosity [m²/s]
   !! @param BackgroundVV Background concentration [v/v]
   !! @param species_index Index in species array
   !! @param drydep_index Index in dry deposition array
   !! @param photolysis_index Index in photolysis array
   !! @param gocart_aero_index Index in GOCART aerosol array
   !! @param conc Species concentration [v/v] or [kg/kg]
   type :: SpeciesType

      ! Names
      character(len=30) :: long_name   !< Long name for species used for NetCDF attribute "long_name"
      character(len=30) :: short_name  !< Short name for species
      character(len=50) :: description !< Description of species

      ! Logical switches
      logical :: is_gas               !< If true, species is a gas and not an aerosol
      logical :: is_aerosol           !< If true, species is aerosol and not a gas
      logical :: is_tracer            !< If true, species is a tracer and not an aerosol or gas that undergoes chemistry or photolysis
      logical :: is_advected          !< If true, species is advected
      logical :: is_drydep            !< If true, species undergoes dry deposition
      logical :: is_wetdep            !< if true, species undergoes wet deposition
      logical :: is_photolysis        !< If true, species undergoes photolysis
      logical :: is_gocart_aero       !< If true, species is a GOCART aerosol species
      logical :: is_dust              !< If true, species is dust
      logical :: is_seasalt           !< If true, species is sea salt

      ! Numerical properties
      real(kind=fp) :: mw_g                 !< Gaseous molecular weight [g/mol]
      real(kind=fp) :: density              !< Particle density [kg/m³]
      real(kind=fp) :: radius               !< Mean molecular diameter [m]
      real(kind=fp) :: lower_radius         !< Lower radius [m]
      real(kind=fp) :: upper_radius         !< Upper radius [m]
      real(kind=fp) :: viscosity            !< Kinematic viscosity [m²/s]

      ! used for dry deposition
      real(kind=fp) :: dd_f0                !< reactivity factor for oxidation of biological substances
      real(kind=fp) :: dd_hstar             !< Henry’s law constant
      real(kind=fp) :: dd_DvzAerSnow        !< fix dry deposition velocity (cm/s) over ice and snow for certain aerosol species
      real(kind=fp) :: dd_DvzMinVal_snow    !< minimum dry deposition velocity (cm/s) over snow and ice
      real(kind=fp) :: dd_DvzMinVal_land    !< minimum dry deposition velocity (cm/s) over land

      ! used for wet deposition
      !real(kind=fp) :: radius_wet           !< mean molecular diameter in meters for wet conditions (use the same radius for both dry and wet deposition for now)
      real(kind=fp) :: henry_k0             !< Henry’s law solubility constant ( M / atm)
      real(kind=fp) :: henry_cr             !< Henry’s law volatility constant (K)
      real(kind=fp) :: henry_pKa            !< Henry’s Law pH correction factor (seems zeros for all species now)
      real(kind=fp) :: wd_retfactor         !< retention efficiency of species in the liquid cloud condensate as it is converted to precipitation
      logical       :: wd_LiqAndGas         !< whether the ice-to-gas ratio can be computed for this species by co-condensation
      real(kind=fp) :: wd_convfacI2G        !< conversion factor for computing the ice-to-gas ratio by co-condensation when wd_LiqAndGas = .true.
      real(kind=fp) :: wd_rainouteff(3)     !< temperature-dependent (T < 237k;  237 <= T < 258k;  T >= 258k) scale factor for the fraction of rainout.

      !used for settling
      character(len=30) :: mie_name         !< Mie data name associated with this species for settling velocity calculation

      ! Default background concentration
      real(kind=fp) :: BackgroundVV        !< Background concentration [v/v]

      ! Indices
      integer :: species_index        !< Species index in species array
      integer :: drydep_index         !< Dry deposition index in drydep array
      integer :: photolysis_index     !< Photolysis index in photolysis array
      integer :: gocart_aero_index    !< GOCART aerosol index in gocart_aero array

      ! Concentration
      real(kind=fp), POINTER :: conc(:,:,:)             !< Species concentration [v/v] or [kg/kg]

      ! Validation and status
      logical :: is_valid = .false.                     !< Validation status

   contains
      ! Enhanced methods with error handling
      procedure :: init => species_init
      procedure :: validate => species_validate
      procedure :: cleanup => species_cleanup
      procedure :: set_concentration => species_set_concentration
      procedure :: get_concentration => species_get_concentration
      procedure :: copy => species_copy
      procedure :: print_info => species_print_info

   end type SpeciesType

   !> \brief Species management system
   !!
   !! This type provides comprehensive species database management with
   !! enhanced error handling and validation capabilities.
   type :: SpeciesManagerType
      private

      type(SpeciesType), allocatable :: species_db(:)  !< Species database
      integer :: num_species = 0                       !< Number of species
      type(ErrorManagerType) :: error_mgr              !< Integrated error manager
      logical :: is_initialized = .false.              !< Initialization status

   contains
      procedure :: init => species_manager_init
      procedure :: add_species => species_manager_add_species
      procedure :: find_species => species_manager_find_species
      procedure :: validate_database => species_manager_validate_database
      procedure :: load_from_file => species_manager_load_from_file
      procedure :: cleanup => species_manager_cleanup
      procedure :: print_database => species_manager_print_database

   end type SpeciesManagerType

   !
   ! !DEFINED PARAMETERS:
   !
   !=========================================================================
   ! Missing species concentration value if not in restart file and special
   ! background value not defined
   !=========================================================================
   REAL(fp), PARAMETER, PUBLIC :: MISSING_VV  = 1.0e-20_fp !< Missing species concentration value

contains

   !> \brief Initialize a species with enhanced validation
   !!
   !! This subroutine initializes a SpeciesType object with comprehensive
   !! error checking and validation.
   !!
   !! \param[inout] this The species object to initialize
   !! \param[in] species_name Short name identifier for the species
   !! \param[in] long_name Descriptive long name
   !! \param[in] molecular_weight Molecular weight [g/mol]
   !! \param[out] rc Return code
   subroutine species_init(this, species_name, long_name, molecular_weight, rc)
      implicit none
      class(SpeciesType), intent(inout) :: this
      character(len=*), intent(in) :: species_name
      character(len=*), intent(in) :: long_name
      real(fp), intent(in) :: molecular_weight
      integer, intent(out) :: rc

      ! Validate inputs
      if (len_trim(species_name) == 0) then
         rc = ERROR_INVALID_INPUT
         return
      endif

      if (molecular_weight <= 0.0_fp) then
         rc = ERROR_INVALID_INPUT
         return
      endif

      ! Initialize species properties
      this%short_name = trim(species_name)
      this%long_name = trim(long_name)
      this%description = ''  ! Initialize description
      this%mw_g = molecular_weight

      ! Set defaults
      this%is_gas = .true.
      this%is_aerosol = .false.
      this%is_tracer = .false.
      this%is_advected = .true.
      this%is_drydep = .false.
      this%is_wetdep = .false.  ! Initialize wet deposition flag
      this%is_photolysis = .false.
      this%is_gocart_aero = .false.
      this%is_dust = .false.
      this%is_seasalt = .false.

      this%density = 1000.0_fp  ! Default density
      this%radius = 1.0e-9_fp   ! Default radius
      this%lower_radius = 0.0_fp
      this%upper_radius = 0.0_fp
      this%viscosity = 1.0e-5_fp

      ! Initialize dry deposition properties
      this%dd_f0 = 0.0_fp
      this%dd_hstar = 0.0_fp
      this%dd_DvzAerSnow = 0.0_fp
      this%dd_DvzMinVal_snow = 0.0_fp
      this%dd_DvzMinVal_land = 0.0_fp

      ! Initialize wet deposition properties
      this%henry_k0 = 0.0_fp
      this%henry_cr = 0.0_fp
      this%henry_pKa = 0.0_fp
      this%wd_retfactor = 0.0_fp
      this%wd_LiqAndGas = .false.
      this%wd_convfacI2G = 0.0_fp
      this%wd_rainouteff(:) = 0.0_fp

      this%BackgroundVV = MISSING_VV
      this%mie_name = ''  ! Initialize Mie name to empty

      this%species_index = -1
      this%drydep_index = -1
      this%photolysis_index = -1
      this%gocart_aero_index = -1

      this%is_valid = .true.
      rc = CC_SUCCESS

   end subroutine species_init

   !> \brief Validate species properties
   !!
   !! This function validates all species properties for consistency and
   !! physical validity.
   !!
   !! \param[in] this The species object to validate
   !! \param[out] rc Return code
   function species_validate(this, rc) result(is_valid)
      implicit none
      class(SpeciesType), intent(in) :: this
      integer, intent(out) :: rc
      logical :: is_valid

      is_valid = .true.
      rc = CC_SUCCESS

      ! Check basic properties
      if (len_trim(this%short_name) == 0) then
         is_valid = .false.
         rc = ERROR_INVALID_INPUT
         return
      endif

      if (this%mw_g <= 0.0_fp) then
         is_valid = .false.
         rc = ERROR_INVALID_INPUT
         return
      endif

      ! Check logical consistency
      if (this%is_gas .and. this%is_aerosol) then
         is_valid = .false.
         rc = ERROR_STATE_INCONSISTENCY
         return
      endif

      if (this%is_dust .and. .not. this%is_aerosol) then
         is_valid = .false.
         rc = ERROR_STATE_INCONSISTENCY
         return
      endif

      if (this%is_seasalt .and. .not. this%is_aerosol) then
         is_valid = .false.
         rc = ERROR_STATE_INCONSISTENCY
         return
      endif

      ! Check physical properties for aerosols
      if (this%is_aerosol) then
         if (this%density <= 0.0_fp) then
            is_valid = .false.
            rc = ERROR_INVALID_INPUT
            return
         endif

         if (this%radius <= 0.0_fp) then
            is_valid = .false.
            rc = ERROR_INVALID_INPUT
            return
         endif
      endif

      ! Check concentration array: all values must be positive and finite
      if (associated(this%conc)) then
         if (any(this%conc < 0.0_fp)) then
            is_valid = .false.
            rc = ERROR_INVALID_INPUT
            return
         endif
         if (any(.not. ieee_is_finite(this%conc))) then
            is_valid = .false.
            rc = ERROR_INVALID_INPUT
            return
         endif
      endif

   end function species_validate

   !> \brief Set species concentration with bounds checking
   !!
   !! This subroutine sets the species concentration with validation.
   !!
   !! \param[inout] this The species object
   !! \param[in] concentration New concentration value
   !! \param[in] grid_index Grid index for concentration
   !! \param[out] rc Return code
   subroutine species_set_concentration(this, concentration, grid_index, rc)
      implicit none
      class(SpeciesType), intent(inout) :: this
      real(fp), intent(in) :: concentration
      integer, intent(in), dimension(3) :: grid_index
      integer, intent(out) :: rc

      ! Check if concentration array is allocated
      if (.not. associated(this%conc)) then
         rc = ERROR_STATE_INCONSISTENCY
         return
      endif

      ! Check bounds
      if (any(grid_index < 1) .or. any(grid_index > shape(this%conc))) then
         rc = ERROR_BOUNDS_CHECK
         return
      endif

      ! Check for valid concentration
      if (concentration < 0.0_fp) then
         rc = ERROR_INVALID_INPUT
         return
      endif

      this%conc(grid_index(1), grid_index(2), grid_index(3)) = concentration
      rc = CC_SUCCESS

   end subroutine species_set_concentration

   !> \brief Cleanup species resources
   !!
   !! This subroutine deallocates all allocated arrays and resets the species.
   !!
   !! \param[inout] this The species object to cleanup
   !! \param[out] rc Return code
   subroutine species_cleanup(this, rc)
      implicit none
      class(SpeciesType), intent(inout) :: this
      integer, intent(out) :: rc

      if (associated(this%conc)) deallocate(this%conc)

      this%description = ''  ! Clear description
      this%mie_name = ''  ! Clear Mie name
      this%is_valid = .false.
      rc = CC_SUCCESS

   end subroutine species_cleanup

   !> \brief Initialize species manager
   !!
   !! This subroutine initializes the species management system.
   !!
   !! \param[inout] this The species manager to initialize
   !! \param[in] max_species Maximum number of species to support
   !! \param[out] rc Return code
   subroutine species_manager_init(this, max_species, rc)
      implicit none
      class(SpeciesManagerType), intent(inout) :: this
      integer, intent(in) :: max_species
      integer, intent(out) :: rc

      if (max_species <= 0) then
         rc = ERROR_INVALID_INPUT
         return
      endif

      allocate(this%species_db(max_species), stat=rc)
      if (rc /= 0) then
         rc = ERROR_MEMORY_ALLOCATION
         return
      endif

      call this%error_mgr%init(verbose=.true.)
      this%num_species = 0
      this%is_initialized = .true.
      rc = CC_SUCCESS

   end subroutine species_manager_init

   !> \brief Find species by name
   !!
   !! This function finds a species in the database by name.
   !!
   !! \param[in] this The species manager
   !! \param[in] species_name Name to search for
   !! \param[out] species_index Index of found species (-1 if not found)
   !! \param[out] rc Return code
   subroutine species_manager_find_species(this, species_name, species_index, rc)
      implicit none
      class(SpeciesManagerType), intent(in) :: this
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: species_index
      integer, intent(out) :: rc

      integer :: i

      species_index = -1
      rc = CC_SUCCESS

      if (.not. this%is_initialized) then
         rc = ERROR_STATE_INCONSISTENCY
         return
      endif

      do i = 1, this%num_species
         if (trim(this%species_db(i)%short_name) == trim(species_name)) then
            species_index = i
            return
         endif
      enddo

      ! Species not found
      rc = ERROR_INVALID_INPUT

   end subroutine species_manager_find_species

   !> \brief Standalone species validation function
   !!
   !! This function provides species validation outside of the class methods.
   !!
   !! \param[in] species Species to validate
   !! \param[out] rc Return code
   function validate_species(species, rc) result(is_valid)
      implicit none
      type(SpeciesType), intent(in) :: species
      integer, intent(out) :: rc
      logical :: is_valid

      is_valid = species%validate(rc)

   end function validate_species

   !> \brief Find species by name (standalone function)
   !!
   !! This function provides species lookup functionality.
   !!
   !! \param[in] species_db Array of species
   !! \param[in] num_species Number of species in database
   !! \param[in] species_name Name to search for
   !! \param[out] rc Return code
   function find_species_by_name(species_db, num_species, species_name, rc) result(species_index)
      implicit none
      type(SpeciesType), intent(in) :: species_db(:)
      integer, intent(in) :: num_species
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      integer :: species_index

      integer :: i

      species_index = -1
      rc = CC_SUCCESS

      do i = 1, num_species
         if (trim(species_db(i)%short_name) == trim(species_name)) then
            species_index = i
            return
         endif
      enddo

      rc = ERROR_INVALID_INPUT

   end function find_species_by_name

   !> \brief Create basic species database
   !!
   !! This subroutine creates a basic species database with common species.
   !!
   !! \param[out] species_mgr Initialized species manager
   !! \param[out] rc Return code
   subroutine create_species_database(species_mgr, rc)
      implicit none
      type(SpeciesManagerType), intent(out) :: species_mgr
      integer, intent(out) :: rc

      call species_mgr%init(100, rc)  ! Support up to 100 species
      if (rc /= CC_SUCCESS) return

      ! Add some common species (this would be expanded)
      ! This is a placeholder implementation
      rc = CC_SUCCESS

   end subroutine create_species_database

   !> \brief Get species concentration at grid point
   !!
   !! This function returns the concentration of the species at a specified
   !! grid point with bounds checking.
   !!
   !! \param[in] this The species object
   !! \param[in] grid_index Grid point index
   !! \param[out] rc Return code
   !! \return Species concentration [v/v] or [kg/kg]
   function species_get_concentration(this, grid_index, rc) result(concentration)
      implicit none
      class(SpeciesType), intent(in) :: this
      integer, intent(in), dimension(3) :: grid_index
      integer, intent(out) :: rc
      real(fp) :: concentration

      rc = CC_SUCCESS
      concentration = 0.0_fp

      ! Check if concentration array is allocated
      if (.not. associated(this%conc)) then
         rc = ERROR_MEMORY_ALLOCATION
         return
      endif

      ! Check bounds
      if (any(grid_index < 1) .or. any(grid_index > shape(this%conc))) then
         rc = ERROR_BOUNDS_CHECK
         return
      endif

      concentration = this%conc(grid_index(1), grid_index(2), grid_index(3))

   end function species_get_concentration

   !> \brief Copy species properties from another species
   !!
   !! This subroutine creates a deep copy of another species object.
   !!
   !! \param[inout] this The destination species object
   !! \param[in] source The source species object to copy from
   !! \param[out] rc Return code
   subroutine species_copy(this, source, rc)
      implicit none
      class(SpeciesType), intent(inout) :: this
      class(SpeciesType), intent(in) :: source
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Copy all scalar properties
      this%short_name = source%short_name
      this%long_name = source%long_name
      this%description = source%description

      ! Copy logical switches
      this%is_gas = source%is_gas
      this%is_aerosol = source%is_aerosol
      this%is_tracer = source%is_tracer
      this%is_advected = source%is_advected
      this%is_drydep = source%is_drydep
      this%is_wetdep = source%is_wetdep
      this%is_photolysis = source%is_photolysis
      this%is_gocart_aero = source%is_gocart_aero
      this%is_dust = source%is_dust
      this%is_seasalt = source%is_seasalt

      ! Copy numerical properties
      this%mw_g = source%mw_g
      this%density = source%density
      this%radius = source%radius
      this%lower_radius = source%lower_radius
      this%upper_radius = source%upper_radius
      this%viscosity = source%viscosity

      ! Copy dry deposition properties
      this%dd_f0 = source%dd_f0
      this%dd_hstar = source%dd_hstar
      this%dd_DvzAerSnow = source%dd_DvzAerSnow
      this%dd_DvzMinVal_snow = source%dd_DvzMinVal_snow
      this%dd_DvzMinVal_land = source%dd_DvzMinVal_land

      ! Copy wet deposition properties
      this%henry_k0 = source%henry_k0
      this%henry_cr = source%henry_cr
      this%henry_pKa = source%henry_pKa
      this%wd_retfactor = source%wd_retfactor
      this%wd_LiqAndGas = source%wd_LiqAndGas
      this%wd_convfacI2G = source%wd_convfacI2G
      this%wd_rainouteff = source%wd_rainouteff

      this%BackgroundVV = source%BackgroundVV
      this%mie_name = source%mie_name

      ! Copy indices
      this%species_index = source%species_index
      this%drydep_index = source%drydep_index
      this%photolysis_index = source%photolysis_index
      this%gocart_aero_index = source%gocart_aero_index

      ! Copy concentration array if allocated
      if (associated(source%conc)) then
         if (associated(this%conc)) deallocate(this%conc)
         allocate(this%conc(size(source%conc,1), size(source%conc,2), size(source%conc,3)), stat=rc)
         if (rc /= 0) then
            rc = ERROR_MEMORY_ALLOCATION
            return
         endif
         this%conc = source%conc
      endif

      this%is_valid = source%is_valid
      rc = CC_SUCCESS

   end subroutine species_copy

   !> \brief Print species information
   !!
   !! This subroutine prints detailed information about the species
   !! to standard output for debugging and diagnostics.
   !!
   !! \param[in] this The species object
   subroutine species_print_info(this)
      implicit none
      class(SpeciesType), intent(in) :: this

      write(*, '(A)') '=== Species Information ==='
      write(*, '(A,A)') 'Short name: ', trim(this%short_name)
      write(*, '(A,A)') 'Long name:  ', trim(this%long_name)
      write(*, '(A,A)') 'Description: ', trim(this%description)
      write(*, '(A,F12.6)') 'Molecular weight [g/mol]: ', this%mw_g
      write(*, '(A,F12.3)') 'Density [kg/m³]: ', this%density
      write(*, '(A,E12.5)') 'Radius [m]: ', this%radius
      write(*, '(A,L1)') 'Is gas: ', this%is_gas
      write(*, '(A,L1)') 'Is aerosol: ', this%is_aerosol
      write(*, '(A,L1)') 'Is tracer: ', this%is_tracer
      write(*, '(A,L1)') 'Is advected: ', this%is_advected
      write(*, '(A,L1)') 'Undergoes dry deposition: ', this%is_drydep
      write(*, '(A,L1)') 'Undergoes wet deposition: ', this%is_wetdep
      write(*, '(A,L1)') 'Undergoes photolysis: ', this%is_photolysis
      write(*, '(A,L1)') 'Is GOCART aerosol: ', this%is_gocart_aero
      write(*, '(A,L1)') 'Is dust: ', this%is_dust
      write(*, '(A,L1)') 'Is seasalt: ', this%is_seasalt
      write(*, '(A,A)') 'Mie data name: ', trim(this%mie_name)
      write(*, '(A,E12.5)') 'Background concentration [v/v]: ', this%BackgroundVV
      write(*, '(A,I0)') 'Species index: ', this%species_index
      if (associated(this%conc)) then
         write(*, '(A,I0)') 'Concentration grid size: ', size(this%conc)
      else
         write(*, '(A)') 'Concentration: Not allocated'
      endif
      write(*, '(A,L1)') 'Valid: ', this%is_valid
      write(*, '(A)') '=========================='

   end subroutine species_print_info

   !> \brief Add a species to the manager database
   !!
   !! This subroutine adds a new species to the species manager database
   !! with validation and error checking.
   !!
   !! \param[inout] this The species manager
   !! \param[in] species The species to add
   !! \param[out] rc Return code
   subroutine species_manager_add_species(this, species, rc)
      implicit none
      class(SpeciesManagerType), intent(inout) :: this
      type(SpeciesType), intent(in) :: species
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      if (.not. this%is_initialized) then
         rc = ERROR_PROCESS_INITIALIZATION
         return
      endif

      if (.not. allocated(this%species_db)) then
         rc = ERROR_MEMORY_ALLOCATION
         return
      endif

      if (this%num_species >= size(this%species_db)) then
         rc = ERROR_BOUNDS_CHECK
         return
      endif

      ! Add species and increment counter
      this%num_species = this%num_species + 1
      call this%species_db(this%num_species)%copy(species, rc)
      if (rc /= CC_SUCCESS) return

      ! Set the species index
      this%species_db(this%num_species)%species_index = this%num_species

   end subroutine species_manager_add_species

   !> \brief Validate the entire species database
   !!
   !! This subroutine validates all species in the database for consistency
   !! and physical reasonableness.
   !!
   !! \param[inout] this The species manager
   !! \param[out] rc Return code
   subroutine species_manager_validate_database(this, rc)
      implicit none
      class(SpeciesManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i, local_rc
      character(len=100) :: error_msg
      logical :: is_valid

      rc = CC_SUCCESS

      if (.not. this%is_initialized) then
         rc = ERROR_PROCESS_INITIALIZATION
         return
      endif

      ! Validate each species
      do i = 1, this%num_species
         is_valid = this%species_db(i)%validate(local_rc)
         if (local_rc /= CC_SUCCESS .or. .not. is_valid) then
            write(error_msg, '(A,I0)') 'Species validation failed for species index ', i
            call this%error_mgr%report_error(local_rc, error_msg, &
               rc, 'species_manager_validate_database')
            if (rc /= CC_SUCCESS) return
         endif
      enddo

   end subroutine species_manager_validate_database

   !> \brief Load species database from configuration file
   !!
   !! This subroutine loads species definitions from a configuration file.
   !! Currently a placeholder implementation.
   !!
   !! \param[inout] this The species manager
   !! \param[in] filename Configuration file name
   !! \param[out] rc Return code
   subroutine species_manager_load_from_file(this, filename, rc)
      implicit none
      class(SpeciesManagerType), intent(inout) :: this
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      if (.not. this%is_initialized) then
         rc = ERROR_PROCESS_INITIALIZATION
         return
      endif

      ! Placeholder implementation
      ! In a full implementation, this would:
      ! 1. Open and parse the configuration file
      ! 2. Create SpeciesType objects from file data
      ! 3. Add them to the database using add_species
      ! 4. Validate the loaded database

      write(*, '(A,A)') 'INFO: Loading species from file: ', trim(filename)
      write(*, '(A)') 'WARNING: species_manager_load_from_file is a placeholder implementation'

   end subroutine species_manager_load_from_file

   !> \brief Clean up species manager resources
   !!
   !! This subroutine deallocates all resources used by the species manager.
   !!
   !! \param[inout] this The species manager
   !! \param[out] rc Return code
   subroutine species_manager_cleanup(this, rc)
      implicit none
      class(SpeciesManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i, local_rc

      rc = CC_SUCCESS

      ! Clean up each species
      if (allocated(this%species_db)) then
         do i = 1, this%num_species
            call this%species_db(i)%cleanup(local_rc)
            ! Continue even if individual cleanup fails
         enddo
         deallocate(this%species_db)
      endif

      this%num_species = 0
      this%is_initialized = .false.

   end subroutine species_manager_cleanup

   !> \brief Print species database information
   !!
   !! This subroutine prints summary information about all species
   !! in the database.
   !!
   !! \param[in] this The species manager
   subroutine species_manager_print_database(this)
      implicit none
      class(SpeciesManagerType), intent(in) :: this

      integer :: i

      write(*, '(A)') '=== Species Database Summary ==='
      write(*, '(A,L1)') 'Initialized: ', this%is_initialized
      write(*, '(A,I0)') 'Number of species: ', this%num_species

      if (allocated(this%species_db)) then
         write(*, '(A,I0)') 'Database capacity: ', size(this%species_db)

         write(*, '(A)') 'Species list:'
         do i = 1, this%num_species
            write(*, '(I4,A,A,A,A)') i, ': ', trim(this%species_db(i)%short_name), &
               ' (', trim(this%species_db(i)%long_name), ')'
         enddo
      else
         write(*, '(A)') 'Database: Not allocated'
      endif

      write(*, '(A)') '==============================='

   end subroutine species_manager_print_database

   ! function get_name(this) result(species_name)
   !    character(len=30) :: species_name

   !    species_name = this%short_name
   ! end function get_name

   ! function get_atomic_number(this) result(atomic_num)
   !    class(Species), intent(in) :: this
   !    integer :: atomic_num

   !    atomic_num = this%atomic_number
   ! end function get_atomic_number

end module species_mod
