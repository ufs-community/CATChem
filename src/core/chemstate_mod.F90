!> \file chemstate_mod.F90
!! \brief Contains the `ChemStateType` data type and related subroutines and functions.
!!
!!
!! The `ChemState_Mod` module contains the chemstate_mod::chemstatetype data type
!! and related subroutines and functions for managing the state of the chemical model.
!!
!! \ingroup core_modules
!!!>
module ChemState_Mod
   !
   ! USES:
   !
   USE Error_Mod
   USE Precision_Mod
   USE species_mod, only: SpeciesType
   USE GridGeometry_Mod, only: GridGeometryType
   USE GOCART2G_MieMod, only: GOCART2G_Mie

   IMPLICIT NONE
   PRIVATE
   !
   ! !PUBLIC MEMBER FUNCTIONS:
   !PUBLIC :: ChemStateType          ! Main data type
   PUBLIC :: Find_Number_of_Species
   ! Legacy routines - commented out in modernization
   ! PUBLIC :: Find_Index_of_Species
   ! PUBLIC :: FindSpecByName
   ! PUBLIC :: GetSpecConc
   ! PUBLIC :: GetSpecConcByName
   ! PUBLIC :: GetSpecConcByIndex
   ! All legacy allocation functions removed - use modern type-bound procedures:
   ! - ChemState%init() instead of legacy allocation
   ! - ChemState%cleanup() for cleanup
   ! - ChemState%validate() for validation
   !
   ! !Private DATA MEMBERS:
   !
   !=========================================================================
   !

   !> \brief Data type for managing the state of the chemical model.
   !!
   !! \details chemStateType contains the following data members:
   !! \param State: A character string containing the name of this state.
   !! \param nSpecies: The total number of species.
   !! \param nSpeciesGas: The number of gas species.
   !! \param nSpeciesAero: The number of aerosol species.
   !! \param nSpeciesDust: The number of dust species.
   !! \param nSpeicesSeaSalt: The number of sea salt species.
   !! \param SpeciesIndex: An array containing the total species index.
   !! \param AeroIndex: An array containing the aerosol species index.
   !! \param GasIndex: An array containing the gas species index.
   !! \param DustIndex: An array containing the dust species index.
   !! \param SeaSaltIndex: An array containing the sea salt species index.
   !! \param chemSpecies: A 2-D array containing the concentration of each species.
   !!
   !! \ingroup core_modules
   !!!>
   type, public :: ChemStateType
      !---------------------------------------------------------------------
      ! Name of variables containing chemistry information
      !---------------------------------------------------------------------
      CHARACTER(LEN=4)  :: State     = 'Chem'    ! Name of this state

      !---------------------------------------------------------------------
      ! Integers
      !---------------------------------------------------------------------
      INTEGER              :: nSpecies          ! Total Number of Species
      INTEGER              :: nSpeciesGas       ! Number of Gas Species
      INTEGER              :: nSpeciesAero      ! Number of Aerosol Species
      INTEGER              :: nSpeciesAeroDryDep ! Number of Aerosol Species for Dry Dep
      INTEGER              :: nSpeciesDryDep    ! Number of DryDep Species
      INTEGER              :: nSpeciesWetDep    ! Number of WetDep Species
      INTEGER              :: nSpeciesTracer    ! Number of Tracer Species
      INTEGER              :: nSpeciesDust      ! Number of Dust Species
      INTEGER              :: nSpeciesSeaSalt   ! Number of SeaSalt Species
      INTEGER, ALLOCATABLE :: SpeciesIndex(:)   ! Total Species Index
      INTEGER, ALLOCATABLE :: TracerIndex(:)    ! Tracer Species Index
      INTEGER, ALLOCATABLE :: AeroIndex(:)      ! Aerosol Species Index
      INTEGER, ALLOCATABLE :: AeroDryDepIndex(:) ! Aerosol DryDep Species Index
      INTEGER, ALLOCATABLE :: GasIndex(:)       ! Gas Species Index
      INTEGER, ALLOCATABLE :: DustIndex(:)      ! Dust Species Index
      INTEGER, ALLOCATABLE :: SeaSaltIndex(:)   ! SeaSalt Species Index
      INTEGER, ALLOCATABLE :: DryDepIndex(:)   ! DryDep Species Index
      INTEGER, ALLOCATABLE :: WetDepIndex(:)   ! WetDep Species Index
      CHARACTER(len=50), ALLOCATABLE :: SpeciesNames(:)  ! Species Names
      type(GOCART2G_Mie), ALLOCATABLE :: MieData(:) ! Mie data for aerosols
      CHARACTER(len=50), ALLOCATABLE :: MieNames(:) ! Mie species names
      INTEGER, ALLOCATABLE :: SpcMieMap(:)   ! Mapping from species name to Mie data

      !---------------------------------------------------------------------
      ! Reals
      !---------------------------------------------------------------------
      type(SpeciesType), allocatable :: ChemSpecies(:)
      type(GridGeometryType), pointer :: Grid => null()  ! Pointer to grid geometry

   contains
      ! Type-bound procedures for modern initialization and cleanup
      procedure :: init => chemstate_init
      procedure :: cleanup => chemstate_cleanup
      procedure :: validate => chemstate_validate
      procedure :: reset => chemstate_reset
      procedure :: is_allocated => chemstate_is_allocated
      procedure :: get_memory_usage => chemstate_get_memory_usage
      procedure :: print_summary => chemstate_print_summary
      procedure :: find_species => chemstate_find_species
      procedure :: get_concentration => chemstate_get_concentration
      procedure :: set_concentration => chemstate_set_concentration
      procedure :: get_column_ptr => chemstate_get_column_ptr
      procedure :: get_num_species => chemstate_get_num_species
      procedure :: get_species => chemstate_get_species

      ! Concentration accessors - multiple interfaces for flexibility
      procedure :: get_concentrations => chemstate_get_concentrations          ! Get column at i,j
      procedure :: set_concentrations => chemstate_set_concentrations          ! Set column at i,j
      procedure :: get_all_concentrations => chemstate_get_all_concentrations  ! Get full 4D array
      procedure :: set_all_concentrations => chemstate_set_all_concentrations  ! Set full 4D array

      procedure :: has_species => chemstate_has_species
      procedure :: get_dimensions => chemstate_get_dimensions
      procedure :: init_mie_data => chemstate_init_mie_data
   end type ChemStateType

CONTAINS


   ! Legacy Chem_Allocate procedure (deprecated)
   ! This has been replaced by the modern chemstate_init procedure
   ! which takes explicit parameters instead of GridState
   !
   ! subroutine Chem_Allocate(GridState, ChemState, RC)
   !    USE GridState_Mod,  ONLY : GridStateType
   !    ...
   ! end subroutine Chem_Allocate

   !> \brief Find the number of species
   !!
   !! This subroutine finds the number of species
   !!
   !! \param ChemState The ChemState object
   !! \param RC The return code
   !!
   !! \ingroup core_modules
   !!!>
   subroutine Find_Number_of_Species(ChemState, RC)
      ! USES
      USE Species_Mod,  ONLY :SpeciesType

      IMPLICIT NONE

      ! INOUT Params
      type(ChemStateType), INTENT(inout) :: ChemState     ! chem State object
      ! OUTPUT Params
      INTEGER,             INTENT(OUT)   :: RC            ! Success or failure

      ! Error handling
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: thisLoc

      ! Local variables
      INTEGER :: i

      ! Initialize
      RC = CC_SUCCESS
      ErrMsg = ''
      thisLoc = ' -> at Find_Number_of_Species (in core/chemstate_mod.F90)'

      ! Initialize to zero before counting species
      ChemState%nSpeciesAero = 0
      ChemState%nSpeciesAeroDryDep = 0
      ChemState%nSpeciesDryDep = 0
      ChemState%nSpeciesWetDep = 0
      ChemState%nSpeciesDust = 0
      ChemState%nSpeciesGas = 0
      ChemState%nSpeciesSeaSalt = 0
      ChemState%nSpeciesTracer = 0

      ! Count number of species
      do i = 1, ChemState%nSpecies
         if (ChemState%ChemSpecies(i)%is_gas .eqv. .true.) then
            ChemState%nSpeciesGas = ChemState%nSpeciesGas + 1
         endif
         if (ChemState%ChemSpecies(i)%is_aerosol .eqv. .true.) then
            ChemState%nSpeciesAero = ChemState%nSpeciesAero + 1
         endif
         if (ChemState%ChemSpecies(i)%is_dust .eqv. .true.) then
            ChemState%nSpeciesDust = ChemState%nSpeciesDust + 1
         endif
         if (ChemState%ChemSpecies(i)%is_seasalt .eqv. .true.) then
            ChemState%nSpeciesSeaSalt = ChemState%nSpeciesSeaSalt + 1
         endif
         if (ChemState%ChemSpecies(i)%is_tracer .eqv. .true.) then
            ChemState%nSpeciesTracer = ChemState%nSpeciesTracer + 1
         endif
         if (ChemState%ChemSpecies(i)%is_drydep .eqv. .true.) then
            ChemState%nSpeciesDryDep = ChemState%nSpeciesDryDep + 1
         endif
         if (ChemState%ChemSpecies(i)%is_drydep .eqv. .true. .and. &
            ChemState%ChemSpecies(i)%is_aerosol .eqv. .true.) then
            ChemState%nSpeciesAeroDryDep = ChemState%nSpeciesAeroDryDep + 1
         endif
         if (ChemState%ChemSpecies(i)%is_wetdep .eqv. .true.) then
            ChemState%nSpeciesWetDep = ChemState%nSpeciesWetDep + 1
         endif
      enddo

   end subroutine Find_Number_of_Species

   !> \brief Find the indices of species (LEGACY - COMMENTED OUT)
   !!
   !! \param ChemState The ChemState object
   !! \param RC The return code
   !!
   !! \ingroup core_modules
   !!!>
   ! subroutine Find_Index_of_Species(ChemState, RC)
   !    ! USES
   !    USE Species_Mod,  ONLY : SpeciesType
   !
   !    IMPLICIT NONE
   !
   !    ! INOUT Params
   !    type(ChemStateType),  INTENT(INOUT) :: ChemState     ! chem State object
   !    ! OUTPUT Params
   !    INTEGER,             INTENT(OUT)   :: RC            ! Success or failure
   !
   !    ! Error handling
   !    CHARACTER(LEN=255) :: ErrMsg
   !    CHARACTER(LEN=255) :: thisLoc
   !
   !    ! Local variables
   !    integer :: n ! looping variable
   !    integer :: aero_index      ! Current Aerosol Index
   !    integer :: gas_index       ! Current Gas Index
   !    integer :: dust_index      ! Current Dust Index
   !    integer :: seasalt_index   ! Current Seas Salt Index
   !    integer :: tracer_index    ! Current Tracer Index
   !    integer :: drydep_index    ! Current DryDep Index
   !
   !
   !    ! Initialize
   !    RC = CC_SUCCESS
   !    ErrMsg = ''
   !    thisLoc = ' -> at Find_indices_of_Species (in core/chemstate_mod.F90)'
   !
   !
   !    ! Initialize to zero before counting species
   !    aero_index = 1
   !    gas_index = 1
   !    dust_index = 1
   !    seasalt_index = 1
   !    tracer_index = 1
   !    drydep_index = 1
   !
   !    ! Allocate index arrays
   !    ALLOCATE(Chemstate%AeroIndex(ChemState%nSpeciesAero), STAT=RC)
   !    IF ( RC /= CC_SUCCESS ) THEN
   !       errMsg = 'Error allocating Chemstate%AeroIndex'
   !       call CC_Error(errMsg, RC, thisLoc)
   !       RETURN
   !    ENDIF
   !
   !    ALLOCATE(Chemstate%TracerIndex(ChemState%nSpeciesTracer), STAT=RC)
   !    IF ( RC /= CC_SUCCESS ) THEN
   !       errMsg = 'Error allocating Chemstate%TracerIndex'
   !       call CC_Error(errMsg, RC, thisLoc)
   !       RETURN
   !    ENDIF
   !
   !    ALLOCATE(Chemstate%GasIndex(ChemState%nSpeciesGas), STAT=RC)
   !    IF ( RC /= CC_SUCCESS ) THEN
   !       errMsg = 'Error allocating Chemstate%GasIndex'
   !       call CC_Error(errMsg, RC, thisLoc)
   !       RETURN
   !    ENDIF
   !
   !    ALLOCATE(Chemstate%DustIndex(ChemState%nSpeciesDust), STAT=RC)
   !    IF ( RC /= CC_SUCCESS ) THEN
   !       errMsg = 'Error allocating Chemstate%DustIndex'
   !       call CC_Error(errMsg, RC, thisLoc)
   !       RETURN
   !    ENDIF
   !
   !    ALLOCATE(Chemstate%SeaSaltIndex(ChemState%nSpeciesSeaSalt), STAT=RC)
   !    IF ( RC /= CC_SUCCESS ) THEN
   !       errMsg = 'Error allocating Chemstate%SeaSaltIndex'
   !       call CC_Error(errMsg, RC, thisLoc)
   !       RETURN
   !    ENDIF
   !
   !    ALLOCATE(Chemstate%DryDepIndex(ChemState%nSpeciesAeroDryDep), STAT=RC)
   !    IF ( RC /= CC_SUCCESS ) THEN
   !       errMsg = 'Error allocating Chemstate%DryDepIndex'
   !       call CC_Error(errMsg, RC, thisLoc)
   !
   !    ENDIF
   !
   !    ! Find indices for species groups
   !    do n = 1, ChemState%nSpecies
   !       if (ChemState%ChemSpecies(n)%is_aerosol .eqv. .true.) then
   !          Chemstate%AeroIndex(aero_index) = n
   !          aero_index = aero_index + 1
   !       endif
   !       if (ChemState%ChemSpecies(n)%is_gas .eqv. .true.) then
   !          Chemstate%GasIndex(gas_index) = n
   !          gas_index = gas_index + 1
   !       endif
   !       if (ChemState%ChemSpecies(n)%is_dust .eqv. .true.) then
   !          Chemstate%DustIndex(dust_index) = n
   !          dust_index = dust_index + 1
   !       endif
   !       if (ChemState%ChemSpecies(n)%is_seasalt .eqv. .true.) then
   !          Chemstate%SeaSaltIndex(seasalt_index) = n
   !          seasalt_index = seasalt_index + 1
   !       endif
   !       if (ChemState%ChemSpecies(n)%is_tracer .eqv. .true.) then
   !          Chemstate%TracerIndex(tracer_index) = n
   !          tracer_index = tracer_index + 1
   !       endif
   !       if (ChemState%ChemSpecies(n)%is_drydep .eqv. .true.) then
   !          Chemstate%DryDepIndex(drydep_index) = n
   !          drydep_index = drydep_index + 1
   !       endif
   !    enddo
   !
   ! end subroutine Find_index_of_Species

   !> \brief Find the species by name (LEGACY - COMMENTED OUT)
   !!
   !! \param ChemState The ChemState object
   !! \param name The name of the species
   !! \param index The index of the species
   !! \param RC The return code
   !!
   !! \ingroup core_modules
   !!!>
   ! subroutine FindSpecByName(ChemState, name, index, RC)
   !
   !    type(ChemStateType),  INTENT(INOUT) :: ChemState     ! chem State object
   !    character(len=50),    INTENT(in)    :: name
   !    integer,              INTENT(out)   :: index
   !    integer,              INTENT(out)   :: RC
   !
   !    ! Error handling
   !    CHARACTER(LEN=255) :: ErrMsg
   !    CHARACTER(LEN=255) :: thisLoc
   !
   !    ! local variables
   !    integer :: n
   !
   !    ! Initialize
   !    RC = CC_SUCCESS
   !    ErrMsg = ''
   !    thisLoc = ' -> at FindSpecByName (in core/chemstate_mod.F90)'
   !
   !    index = 0
   !    do n = 1, ChemState%nSpecies
   !       if (TRIM(name) == TRIM(ChemState%SpeciesNames(n))) then
   !          index = n
   !          exit
   !       endif
   !    enddo
   !    if (index == 0) then
   !       RC = CC_FAILURE
   !       ErrMsg = 'Species not found: ' // TRIM(name)
   !       call CC_Warning(ErrMsg, RC, thisLoc)
   !    endif
   !
   ! end subroutine FindSpecByName
   !
   ! subroutine FindSpecByName(ChemState, name, index, RC)
   !
   !    type(ChemStateType),  INTENT(INOUT) :: ChemState     ! chem State object
   !    character(len=50),    INTENT(in)    :: name
   !    integer,              INTENT(out)   :: index
   !    integer,              INTENT(out)   :: RC
   !
   !    ! Error handling
   !    CHARACTER(LEN=255) :: ErrMsg
   !    CHARACTER(LEN=255) :: thisLoc
   !    integer :: n
   !
   !    ! Initialize
   !    RC = CC_SUCCESS
   !    ErrMsg = ''
   !    thisLoc = ' -> at FindSpecByName (in core/chemstate_mod.F90)'
   !
   !    index = 0
   !    do n = 1, ChemState%nSpecies
   !       if (TRIM(name) == TRIM(ChemState%SpeciesNames(n))) then
   !          index = n
   !          exit
   !       endif
   !    enddo
   !    if (index == 0) then
   !       RC = CC_FAILURE
   !       ErrMsg = 'Species not found: ' // TRIM(name)
   !       call CC_Warning(ErrMsg, RC, thisLoc)
   !    endif
   !
   ! end subroutine FindSpecByName


   !> \brief Get the concentration of a species (LEGACY - COMMENTED OUT)
   !!
   !! get the concentration of a species given either the index or the name of the species
   !!
   !! \param ChemState The ChemState object
   !! \param concentration The concentration of the species
   !! \param RC The return code
   !! \param index The index of the species - Optional
   !! \param name The name of the species - Optional
   !!
   !! \ingroup core_modules
   !!!>
   ! subroutine GetSpecConc(ChemState, concentration, RC, index, name)
   !
   !    type(ChemStateType),  INTENT(INOUT) :: ChemState     ! chem State object
   !    real(kind=fp), dimension(:), INTENT(out)   :: concentration
   !    integer,              INTENT(out)   :: RC
   !    integer, optional,    INTENT(inout)    :: index
   !    character(len=50), optional, INTENT(inout)    :: name
   !
   !    ! Error handling
   !    CHARACTER(LEN=255) :: ErrMsg
   !    CHARACTER(LEN=255) :: thisLoc
   !
   !    ! Initialize
   !    RC = CC_SUCCESS
   !    ErrMsg = ''
   !    thisLoc = ' -> at GetSpecConc (in core/chemstate_mod.F90)'
   !
   !    if (present(index)) then
   !       call GetSpecConcByIndex(ChemState, concentration, index, RC)
   !    elseif (present(name)) then
   !       call GetSpecConcByName(ChemState, concentration, name, RC)
   !    else
   !       RC = CC_FAILURE
   !    endif
   !
   !    if (RC /= CC_SUCCESS) then
   !       errMsg = 'Error in GetSpecConc'
   !       call CC_Error(errMsg, RC, thisLoc)
   !       RETURN
   !    endif
   !
   ! end subroutine GetSpecConc

   !> \brief Get the concentration of a species by index (LEGACY - COMMENTED OUT)
   !!
   !! \param ChemState The ChemState object
   !! \param concentration The concentration of the species
   !! \param RC The return code
   !! \param index The index of the species
   !!
   !! \ingroup core_modules
   !!!>
   ! subroutine GetSpecConcByIndex(ChemState, concentration, index, RC)
   !
   !    type(ChemStateType),  INTENT(INOUT) :: ChemState     ! chem State object
   !    real(kind=fp), dimension(:), INTENT(out)   :: concentration
   !    integer,              INTENT(in)    :: index
   !    integer,              INTENT(out)   :: RC
   !
   !    ! Error handling
   !    CHARACTER(LEN=255) :: ErrMsg
   !    CHARACTER(LEN=255) :: thisLoc
   !
   !    ! Initialize
   !    RC = CC_SUCCESS
   !    ErrMsg = ''
   !    thisLoc = ' -> at GetSpecConcByIndex (in core/chemstate_mod.F90)'
   !
   !    if (index < 1 .or. index > ChemState%nSpecies) then
   !       RC = CC_FAILURE
   !       errMsg = 'index out of bounds'
   !       call CC_Error(errMsg, RC, thisLoc)
   !       RETURN
   !    endif
   !
   !    concentration = ChemState%ChemSpecies(index)%conc
   !
   ! end subroutine GetSpecConcByIndex

   !> \brief Get the concentration of a species by name (LEGACY - COMMENTED OUT)
   !!
   !! \param ChemState The ChemState object
   !! \param concentration The concentration of the species
   !! \param RC The return code
   !! \param name The name of the species
   !!
   !! \ingroup core_modules
   !!!>
   ! subroutine GetSpecConcByName(ChemState, concentration, name, RC)
   !
   !    type(ChemStateType),  INTENT(INOUT) :: ChemState     ! chem State object
   !    real(kind=fp), dimension(:), INTENT(out)   :: concentration
   !    character(len=50),    INTENT(in)    :: name
   !    integer,              INTENT(out)   :: RC
   !
   !    ! Locals
   !    integer :: index
   !
   !    ! Error handling
   !    CHARACTER(LEN=255) :: ErrMsg
   !    CHARACTER(LEN=255) :: thisLoc
   !
   !    ! Initialize
   !    RC = CC_SUCCESS
   !    ErrMsg = ''
   !    thisLoc = ' -> at GetSpecConcByName (in core/chemstate_mod.F90)'
   !
   !    call FindSpecByName(ChemState, name, index, RC)
   !
   !    if (RC /= CC_SUCCESS) then
   !       errMsg = 'Error in GetSpecConcByName'
   !       call CC_Error(errMsg, RC, thisLoc)
   !       RETURN
   !    endif
   !
   !    concentration = ChemState%ChemSpecies(index)%conc
   !
   ! end subroutine GetSpecConcByName

   !========================================================================
   ! Modern ChemState Type-Bound Procedures
   !========================================================================

   !> \brief Modern initialization procedure for ChemStateType
   !!
   !! This procedure initializes a ChemStateType object with proper error handling
   !! and validation. It replaces the old allocation patterns with a cleaner approach.
   !!
   !! \param[inout] this The ChemStateType object to initialize
   !! \param[in] max_species Maximum number of species to allocate for
   !! \param[in] error_mgr Error manager for context and error reporting
   !! \param[out] rc Return code
   subroutine chemstate_init(this, max_species, error_mgr, rc, grid)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, ERROR_MEMORY_ALLOCATION
      use GridGeometry_Mod, only: GridGeometryType
      implicit none
      class(ChemStateType), intent(inout) :: this
      integer, intent(in) :: max_species
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc
      type(GridGeometryType), pointer, optional, intent(in) :: grid
      character(len=256) :: thisLoc
      integer :: allocStat, nx, ny, nz, s

      thisLoc = 'chemstate_init (in core/chemstate_mod.F90)'
      call error_mgr%push_context('chemstate_init', 'initializing chemistry state')

      rc = CC_SUCCESS

      ! Initialize basic state information
      this%State = 'Chem'
      this%nSpecies = 0
      this%nSpeciesGas = 0
      this%nSpeciesAero = 0
      this%nSpeciesAeroDryDep = 0
      this%nSpeciesDryDep = 0
      this%nSpeciesWetDep = 0
      this%nSpeciesTracer = 0
      this%nSpeciesDust = 0
      this%nSpeciesSeaSalt = 0

      ! Store grid pointer if provided
      if (present(grid)) then
         this%Grid => grid
      else
         this%Grid => null()
      endif

      ! Allocate species arrays
      if (max_species > 0) then
         allocate(this%SpeciesIndex(max_species), stat=allocStat)
         if (allocStat /= 0) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
               'Failed to allocate SpeciesIndex', rc, &
               thisLoc, 'Check available memory')
            call error_mgr%pop_context()
            return
         endif

         allocate(this%TracerIndex(max_species), stat=allocStat)
         if (allocStat /= 0) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
               'Failed to allocate TracerIndex', rc, &
               thisLoc, 'Check available memory')
            call error_mgr%pop_context()
            return
         endif

         allocate(this%AeroIndex(max_species), stat=allocStat)
         if (allocStat /= 0) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
               'Failed to allocate AeroIndex', rc, &
               thisLoc, 'Check available memory')
            call error_mgr%pop_context()
            return
         endif

         allocate(this%GasIndex(max_species), stat=allocStat)
         if (allocStat /= 0) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
               'Failed to allocate GasIndex', rc, &
               thisLoc, 'Check available memory')
            call error_mgr%pop_context()
            return
         endif

         allocate(this%DustIndex(max_species), stat=allocStat)
         if (allocStat /= 0) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
               'Failed to allocate DustIndex', rc, &
               thisLoc, 'Check available memory')
            call error_mgr%pop_context()
            return
         endif

         allocate(this%SeaSaltIndex(max_species), stat=allocStat)
         if (allocStat /= 0) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
               'Failed to allocate SeaSaltIndex', rc, &
               thisLoc, 'Check available memory')
            call error_mgr%pop_context()
            return
         endif

         allocate(this%DryDepIndex(max_species), stat=allocStat)
         if (allocStat /= 0) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
               'Failed to allocate DryDepIndex', rc, &
               thisLoc, 'Check available memory')
            call error_mgr%pop_context()
            return
         endif

         allocate(this%WetDepIndex(max_species), stat=allocStat)
         if (allocStat /= 0) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
               'Failed to allocate WetDepIndex', rc, &
               thisLoc, 'Check available memory')
            call error_mgr%pop_context()
            return
         endif

         allocate(this%AeroDryDepIndex(max_species), stat=allocStat)
         if (allocStat /= 0) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
               'Failed to allocate AeroDryDepIndex', rc, &
               thisLoc, 'Check available memory')
            call error_mgr%pop_context()
            return
         endif

         allocate(this%SpeciesNames(max_species), stat=allocStat)
         if (allocStat /= 0) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
               'Failed to allocate SpeciesNames', rc, &
               thisLoc, 'Check available memory')
            call error_mgr%pop_context()
            return
         endif

         allocate(this%ChemSpecies(max_species), stat=allocStat)
         if (allocStat /= 0) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
               'Failed to allocate ChemSpecies', rc, &
               thisLoc, 'Check available memory')
            call error_mgr%pop_context()
            return
         endif

         ! Initialize all SpeciesType objects to prevent garbage pointer values
         do s = 1, max_species
            nullify(this%ChemSpecies(s)%conc)
            this%ChemSpecies(s)%is_valid = .false.
         end do

         ! Allocate ChemSpecies(:)%conc(nx,ny,nz) if grid is present
         if (associated(this%Grid)) then
            nx = this%Grid%nx
            ny = this%Grid%ny
            nz = this%Grid%nz

            do s = 1, max_species
               ! Always nullify and reallocate to ensure proper dimensions
               ! Skip trying to deallocate potentially corrupted pointers
               if (associated(this%ChemSpecies(s)%conc)) then
                  nullify(this%ChemSpecies(s)%conc)
               endif

               allocate(this%ChemSpecies(s)%conc(nx,ny,nz), stat=allocStat)
               if (allocStat /= 0) then
                  call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
                     'Failed to allocate ChemSpecies(s)%conc', rc, thisLoc, 'Check available memory')
                  call error_mgr%pop_context()
                  return
               endif
               this%ChemSpecies(s)%conc = 0.0_fp
            end do
         else
            write(*,'(A)') 'DEBUG: Grid not associated, cannot allocate conc arrays'
         endif
      endif
      call error_mgr%pop_context()
   end subroutine chemstate_init

   !> \brief Clean up and deallocate ChemStateType
   subroutine chemstate_cleanup(this, rc)
      implicit none
      class(ChemStateType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Deallocate all allocatable arrays
      if (allocated(this%SpeciesIndex)) deallocate(this%SpeciesIndex)
      if (allocated(this%TracerIndex)) deallocate(this%TracerIndex)
      if (allocated(this%AeroIndex)) deallocate(this%AeroIndex)
      if (allocated(this%GasIndex)) deallocate(this%GasIndex)
      if (allocated(this%DustIndex)) deallocate(this%DustIndex)
      if (allocated(this%SeaSaltIndex)) deallocate(this%SeaSaltIndex)
      if (allocated(this%DryDepIndex)) deallocate(this%DryDepIndex)
      if (allocated(this%WetDepIndex)) deallocate(this%WetDepIndex)
      if (allocated(this%AeroDryDepIndex)) deallocate(this%AeroDryDepIndex)
      if (allocated(this%SpeciesNames)) deallocate(this%SpeciesNames)
      if (allocated(this%ChemSpecies)) deallocate(this%ChemSpecies)
      if (allocated(this%MieData)) deallocate(this%MieData)
      if (allocated(this%MieNames)) deallocate(this%MieNames)
      if (allocated(this%SpcMieMap)) deallocate(this%SpcMieMap)

      ! Clean up grid geometry pointer (nullify only, don't deallocate as we don't own it)
      if (associated(this%Grid)) then
         this%Grid => null()
      endif

      ! Reset counters
      this%nSpecies = 0
      this%nSpeciesGas = 0
      this%nSpeciesAero = 0
      this%nSpeciesAeroDryDep = 0
      this%nSpeciesDryDep = 0
      this%nSpeciesWetDep = 0
      this%nSpeciesTracer = 0
      this%nSpeciesDust = 0
      this%nSpeciesSeaSalt = 0
      this%State = ''

   end subroutine chemstate_cleanup

   !> \brief Validate ChemStateType for consistency
   subroutine chemstate_validate(this, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, ERROR_INVALID_INPUT

      implicit none
      class(ChemStateType), intent(in) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      character(len=256) :: thisLoc

      thisLoc = 'chemstate_validate (in core/chemstate_mod.F90)'
      call error_mgr%push_context('chemstate_validate', 'validating chemistry state')

      rc = CC_SUCCESS

      ! Check that species counts are consistent
      if (this%nSpecies < 0) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
            'Number of species cannot be negative', rc, &
            thisLoc, 'Check species initialization')
         call error_mgr%pop_context()
         return
      endif

      if (this%nSpeciesGas < 0 .or. this%nSpeciesAero < 0) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
            'Species type counts cannot be negative', rc, &
            thisLoc, 'Check species type initialization')
         call error_mgr%pop_context()
         return
      endif

      ! Check consistency of total vs individual counts
      if (this%nSpecies > 0) then
         if (this%nSpeciesGas + this%nSpeciesAero > this%nSpecies) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'Sum of species types exceeds total species', rc, &
               thisLoc, 'Check species counting logic')
            call error_mgr%pop_context()
            return
         endif
      endif

      ! Check array allocation consistency
      if (this%nSpecies > 0 .and. .not. this%is_allocated()) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
            'Species arrays not allocated but nSpecies > 0', rc, &
            thisLoc, 'Call init() to allocate arrays')
         call error_mgr%pop_context()
         return
      endif

      call error_mgr%pop_context()
   end subroutine chemstate_validate

   !> \brief Reset ChemStateType to initial values
   subroutine chemstate_reset(this, rc)
      implicit none
      class(ChemStateType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Reset species counts
      this%nSpecies = 0
      this%nSpeciesGas = 0
      this%nSpeciesAero = 0
      this%nSpeciesAeroDryDep = 0
      this%nSpeciesDryDep = 0
      this%nSpeciesWetDep = 0
      this%nSpeciesTracer = 0
      this%nSpeciesDust = 0
      this%nSpeciesSeaSalt = 0

      ! Reset arrays if allocated
      if (allocated(this%SpeciesIndex)) this%SpeciesIndex = 0
      if (allocated(this%TracerIndex)) this%TracerIndex = 0
      if (allocated(this%AeroIndex)) this%AeroIndex = 0
      if (allocated(this%GasIndex)) this%GasIndex = 0
      if (allocated(this%DustIndex)) this%DustIndex = 0
      if (allocated(this%SeaSaltIndex)) this%SeaSaltIndex = 0
      if (allocated(this%DryDepIndex)) this%DryDepIndex = 0
      if (allocated(this%WetDepIndex)) this%WetDepIndex = 0
      if (allocated(this%AeroDryDepIndex)) this%AeroDryDepIndex = 0
      if (allocated(this%SpeciesNames)) this%SpeciesNames = ''

   end subroutine chemstate_reset

   !> \brief Check if required arrays are allocated
   function chemstate_is_allocated(this) result(is_alloc)
      implicit none
      class(ChemStateType), intent(in) :: this
      logical :: is_alloc

      is_alloc = allocated(this%SpeciesIndex) .and. allocated(this%SpeciesNames) .and. &
         allocated(this%ChemSpecies)
   end function chemstate_is_allocated

   !> \brief Get approximate memory usage in bytes
   function chemstate_get_memory_usage(this) result(memory_bytes)
      implicit none
      class(ChemStateType), intent(in) :: this
      integer(kind=8) :: memory_bytes

      memory_bytes = 0

      ! Estimate based on allocated arrays
      if (allocated(this%SpeciesIndex)) then
         memory_bytes = memory_bytes + size(this%SpeciesIndex) * 4  ! integers
      endif
      if (allocated(this%TracerIndex)) then
         memory_bytes = memory_bytes + size(this%TracerIndex) * 4
      endif
      if (allocated(this%AeroIndex)) then
         memory_bytes = memory_bytes + size(this%AeroIndex) * 4
      endif
      if (allocated(this%GasIndex)) then
         memory_bytes = memory_bytes + size(this%GasIndex) * 4
      endif
      if (allocated(this%DustIndex)) then
         memory_bytes = memory_bytes + size(this%DustIndex) * 4
      endif
      if (allocated(this%SeaSaltIndex)) then
         memory_bytes = memory_bytes + size(this%SeaSaltIndex) * 4
      endif
      if (allocated(this%DryDepIndex)) then
         memory_bytes = memory_bytes + size(this%DryDepIndex) * 4
      endif
      if (allocated(this%WetDepIndex)) then
         memory_bytes = memory_bytes + size(this%WetDepIndex) * 4
      endif
      if (allocated(this%AeroDryDepIndex)) then
         memory_bytes = memory_bytes + size(this%AeroDryDepIndex) * 4
      endif
      if (allocated(this%SpeciesNames)) then
         memory_bytes = memory_bytes + size(this%SpeciesNames) * 50  ! character arrays
      endif
      if (allocated(this%ChemSpecies)) then
         memory_bytes = memory_bytes + size(this%ChemSpecies) * 1000  ! estimate for SpeciesType
      endif

   end function chemstate_get_memory_usage

   !> \brief Print summary of ChemStateType
   subroutine chemstate_print_summary(this)
      implicit none
      class(ChemStateType), intent(in) :: this

      write(*,'(A)') '=== ChemState Summary ==='
      write(*,'(A,A)') 'State: ', trim(this%State)
      write(*,'(A,I0)') 'Total species: ', this%nSpecies
      write(*,'(A,I0)') 'Gas species: ', this%nSpeciesGas
      write(*,'(A,I0)') 'Aerosol species: ', this%nSpeciesAero
      write(*,'(A,I0)') 'Dust species: ', this%nSpeciesDust
      write(*,'(A,I0)') 'Sea salt species: ', this%nSpeciesSeaSalt
      write(*,'(A,I0)') 'Tracer species: ', this%nSpeciesTracer
      write(*,'(A,I0)') 'DryDep species: ', this%nSpeciesDryDep
      write(*,'(A,I0)') 'WetDep species: ', this%nSpeciesWetDep
      write(*,'(A,L1)') 'Arrays allocated: ', this%is_allocated()
      write(*,'(A,I0,A)') 'Memory usage: ', this%get_memory_usage(), ' bytes'
      write(*,'(A)') '========================'
   end subroutine chemstate_print_summary

   !> \brief Find species index by name
   function chemstate_find_species(this, species_name) result(species_index)
      implicit none
      class(ChemStateType), intent(in) :: this
      character(len=*), intent(in) :: species_name
      integer :: species_index

      integer :: i

      species_index = 0  ! Not found

      if (allocated(this%SpeciesNames)) then
         do i = 1, this%nSpecies
            if (trim(this%SpeciesNames(i)) == trim(species_name)) then
               species_index = i
               exit
            endif
         enddo
      endif

   end function chemstate_find_species

   !> \brief Get species concentration by index
   function chemstate_get_concentration(this, species_index) result(concentration)
      use species_mod, only: SpeciesType

      implicit none
      class(ChemStateType), intent(in) :: this
      integer, intent(in) :: species_index
      type(SpeciesType) :: concentration

      if (species_index > 0 .and. species_index <= this%nSpecies) then
         if (allocated(this%ChemSpecies)) then
            concentration = this%ChemSpecies(species_index)
         endif
      endif

   end function chemstate_get_concentration

   !> \brief Set species concentration by index
   subroutine chemstate_set_concentration(this, species_index, concentration, rc)
      use species_mod, only: SpeciesType

      implicit none
      class(ChemStateType), intent(inout) :: this
      integer, intent(in) :: species_index
      type(SpeciesType), intent(in) :: concentration
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      if (species_index > 0 .and. species_index <= this%nSpecies) then
         if (allocated(this%ChemSpecies)) then
            this%ChemSpecies(species_index) = concentration
         else
            rc = CC_FAILURE
         endif
      else
         rc = CC_FAILURE
      endif

   end subroutine chemstate_set_concentration

   !> \brief Get a pointer to a vertical column for a given field name and (i,j) indices (type-safe)
   subroutine chemstate_get_column_ptr(this, field_name, i, j, col_ptr, rc)
      use GridGeometry_Mod, only: GridGeometryType
      class(ChemStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: i, j
      real(fp), pointer :: col_ptr(:)
      integer, intent(out) :: rc
      integer :: species_idx, nz

      rc = CC_FAILURE
      nullify(col_ptr)

      species_idx = -1
      if (allocated(this%SpeciesNames)) then
         do species_idx = 1, size(this%SpeciesNames)
            if (trim(this%SpeciesNames(species_idx)) == trim(field_name)) exit
         end do
         if (species_idx < 1 .or. species_idx > size(this%SpeciesNames)) species_idx = -1
      endif

      if (species_idx > 0 .and. allocated(this%ChemSpecies)) then
         if (associated(this%ChemSpecies(species_idx)%conc)) then
            if (associated(this%Grid)) then
               nz = this%Grid%nz
               if (i >= 1 .and. j >= 1 .and. i <= this%Grid%nx .and. j <= this%Grid%ny) then
                  col_ptr => this%ChemSpecies(species_idx)%conc(i,j,:)
                  if (associated(col_ptr)) then
                     rc = CC_SUCCESS
                     return
                  endif
               endif
            endif
         endif
      endif
      nullify(col_ptr)
      rc = CC_FAILURE
   end subroutine chemstate_get_column_ptr

   !> \brief Get the number of species in the chemical state
   function chemstate_get_num_species(this) result(num_species)
      class(ChemStateType), intent(in) :: this
      integer :: num_species

      if (allocated(this%SpeciesNames)) then
         num_species = size(this%SpeciesNames)
      else
         num_species = 0
      end if
   end function chemstate_get_num_species

   !> \brief Get species names array
   function chemstate_get_species(this) result(species_names)
      class(ChemStateType), intent(in) :: this
      character(len=64), allocatable :: species_names(:)

      if (allocated(this%SpeciesNames)) then
         allocate(species_names(size(this%SpeciesNames)))
         species_names = this%SpeciesNames
      else
         allocate(species_names(0))
      end if
   end function chemstate_get_species

   !> \brief Get concentrations for all species at a grid point
   function chemstate_get_concentrations(this, i, j) result(concentrations)
      class(ChemStateType), intent(in) :: this
      integer, intent(in) :: i, j
      real(fp), allocatable :: concentrations(:,:)
      integer :: s

      if (allocated(this%ChemSpecies) .and. associated(this%Grid)) then
         allocate(concentrations(size(this%ChemSpecies), this%Grid%nz))
         do s = 1, size(this%ChemSpecies)
            if (associated(this%ChemSpecies(s)%conc)) then
               concentrations(s, :) = this%ChemSpecies(s)%conc(i, j, :)
            else
               concentrations(s, :) = 0.0_fp
            end if
         end do
      else
         allocate(concentrations(0, 0))
      end if
   end function chemstate_get_concentrations

   !> \brief Set concentrations for all species at a grid point
   subroutine chemstate_set_concentrations(this, i, j, concentrations, rc)
      class(ChemStateType), intent(inout) :: this
      integer, intent(in) :: i, j
      real(fp), intent(in) :: concentrations(:,:)
      integer, intent(out) :: rc
      integer :: s

      rc = CC_FAILURE
      if (.not. allocated(this%ChemSpecies) .or. .not. associated(this%Grid)) return
      if (size(concentrations, 1) /= size(this%ChemSpecies)) return
      if (size(concentrations, 2) /= this%Grid%nz) return

      do s = 1, size(this%ChemSpecies)
         if (associated(this%ChemSpecies(s)%conc)) then
            this%ChemSpecies(s)%conc(i, j, :) = concentrations(s, :)
         end if
      end do
      rc = CC_SUCCESS
   end subroutine chemstate_set_concentrations

   !> \brief Check if a species exists in the chemical state
   function chemstate_has_species(this, species_name) result(has_it)
      class(ChemStateType), intent(in) :: this
      character(len=*), intent(in) :: species_name
      logical :: has_it

      has_it = (this%find_species(species_name) > 0)
   end function chemstate_has_species

   !> \brief Get grid dimensions
   function chemstate_get_dimensions(this) result(dims)
      class(ChemStateType), intent(in) :: this
      integer :: dims(3)

      if (associated(this%Grid)) then
         dims = [this%Grid%nx, this%Grid%ny, this%Grid%nz]
      else
         dims = [0, 0, 0]
      end if
   end function chemstate_get_dimensions

   !> \brief Get all concentrations as a full 4D array
   !!
   !! Returns the complete concentration array for all species and grid points.
   !! This is useful for processes that need to work with the entire domain.
   !!
   !! \param[in] this ChemStateType instance
   !! \param[out] concentrations Full 4D concentration array (nx, ny, nz, n_species)
   !! \param[out] rc Return code
   subroutine chemstate_get_all_concentrations(this, concentrations, rc)
      class(ChemStateType), intent(in) :: this
      real(fp), allocatable, intent(out) :: concentrations(:,:,:,:)
      integer, intent(out) :: rc
      integer :: s

      rc = CC_FAILURE
      if (.not. allocated(this%ChemSpecies) .or. .not. associated(this%Grid)) return

      ! Allocate the full 4D array
      allocate(concentrations(this%Grid%nx, this%Grid%ny, this%Grid%nz, size(this%ChemSpecies)), stat=rc)
      if (rc /= 0) then
         rc = CC_FAILURE
         return
      end if

      ! Copy data from individual species arrays
      do s = 1, size(this%ChemSpecies)
         if (associated(this%ChemSpecies(s)%conc)) then
            concentrations(:, :, :, s) = this%ChemSpecies(s)%conc(:, :, :)
         else
            concentrations(:, :, :, s) = 0.0_fp
         end if
      end do

      rc = CC_SUCCESS
   end subroutine chemstate_get_all_concentrations

   !> \brief Set all concentrations from a full 4D array
   !!
   !! Updates the complete concentration arrays for all species and grid points.
   !! This is useful for processes that compute tendencies for the entire domain.
   !!
   !! \param[inout] this ChemStateType instance
   !! \param[in] concentrations Full 4D concentration array (nx, ny, nz, n_species)
   !! \param[out] rc Return code
   subroutine chemstate_set_all_concentrations(this, concentrations, rc)
      class(ChemStateType), intent(inout) :: this
      real(fp), intent(in) :: concentrations(:,:,:,:)
      integer, intent(out) :: rc
      integer :: s

      rc = CC_FAILURE
      if (.not. allocated(this%ChemSpecies) .or. .not. associated(this%Grid)) return

      ! Validate array dimensions
      if (size(concentrations, 1) /= this%Grid%nx .or. &
         size(concentrations, 2) /= this%Grid%ny .or. &
         size(concentrations, 3) /= this%Grid%nz .or. &
         size(concentrations, 4) /= size(this%ChemSpecies)) then
         return
      end if

      ! Copy data to individual species arrays
      do s = 1, size(this%ChemSpecies)
         if (associated(this%ChemSpecies(s)%conc)) then
            this%ChemSpecies(s)%conc(:, :, :) = concentrations(:, :, :, s)
         end if
      end do

      rc = CC_SUCCESS
   end subroutine chemstate_set_all_concentrations

   !> \brief Initialize Mie data for aerosol optical properties
   !!
   !! This subroutine allocates and initializes Mie scattering data based on
   !! the configuration file information and species Mie name mappings.
   !!
   !! \param[inout] this ChemStateType object
   !! \param[in] n_mie_files Number of Mie files
   !! \param[in] mie_names Array of Mie type names (e.g., 'SS', 'DU', 'BC')
   !! \param[in] mie_full_paths Array of full file paths to Mie data files
   !! \param[out] rc Return code
   subroutine chemstate_init_mie_data(this, n_mie_files, mie_names, mie_full_paths, rc)
      implicit none
      class(ChemStateType), intent(inout) :: this
      integer, intent(in) :: n_mie_files
      character(len=30), intent(in) :: mie_names(:)
      character(len=512), intent(in) :: mie_full_paths(:)
      integer, intent(out) :: rc

      integer :: i, j, local_rc
      !integer :: channels(4) = [470, 550, 670, 870]  ! Example channels: 470, 550, 670, 870 nm
      character(len=255) :: err_msg
      character(len=255) :: this_loc

      rc = CC_SUCCESS
      this_loc = ' -> at chemstate_init_mie_data (in core/chemstate_mod.F90)'

      ! Allocate MieData and MieNames arrays
      if (allocated(this%MieData)) deallocate(this%MieData)
      if (allocated(this%MieNames)) deallocate(this%MieNames)

      allocate(this%MieData(n_mie_files), stat=rc)
      if (rc /= CC_SUCCESS) then
         err_msg = 'Error allocating MieData array'
         call CC_Error(err_msg, rc, this_loc)
         return
      end if

      allocate(this%MieNames(n_mie_files), stat=rc)
      if (rc /= CC_SUCCESS) then
         err_msg = 'Error allocating MieNames array'
         call CC_Error(err_msg, rc, this_loc)
         return
      end if

      ! Copy Mie names and load Mie data files
      do i = 1, n_mie_files
         this%MieNames(i) = mie_names(i)

         ! Initialize Mie data from file [470 550 670 870] nm for diagnostics
         !this%MieData(i) = GOCART2G_Mie(trim(mie_full_paths(i)), channels*1.e-9, nmom=0, rc=local_rc) !This is for diagMie
         this%MieData(i) = GOCART2G_Mie(trim(mie_full_paths(i)), rc=local_rc)
         if (local_rc /= 0) then
            err_msg = 'Error initializing Mie data for ' // trim(mie_names(i)) // &
               ' from file: ' // trim(mie_full_paths(i))
            rc = local_rc
            call CC_Error(err_msg, rc, this_loc)
            return
         end if
      end do

      ! Allocate and compute species-to-Mie mapping
      if (allocated(this%SpcMieMap)) deallocate(this%SpcMieMap)
      allocate(this%SpcMieMap(this%nSpecies), stat=rc)
      if (rc /= CC_SUCCESS) then
         err_msg = 'Error allocating SpcMieMap array'
         call CC_Error(err_msg, rc, this_loc)
         return
      end if

      ! Initialize mapping to zero (no Mie data)
      this%SpcMieMap(:) = 0

      ! Map species to Mie data based on species mie_name field
      do i = 1, this%nSpecies
         if (len_trim(this%ChemSpecies(i)%mie_name) > 0) then
            ! Find matching Mie data
            do j = 1, n_mie_files
               if (trim(this%ChemSpecies(i)%mie_name) == trim(this%MieNames(j))) then
                  this%SpcMieMap(i) = j
                  exit
               end if
            end do

            ! Warn if no matching Mie data found
            if (this%SpcMieMap(i) == 0) then
               err_msg = 'Warning: No Mie data found for species ' // &
                  trim(this%ChemSpecies(i)%short_name) // ' with mie_name: ' // &
                  trim(this%ChemSpecies(i)%mie_name)
               call CC_Warning(err_msg, rc, this_loc)
            end if
         end if
      end do

      rc = CC_SUCCESS
   end subroutine chemstate_init_mie_data

end module ChemState_Mod
