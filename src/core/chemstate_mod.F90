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
   use species_mod, only: SpeciesType

   IMPLICIT NONE
   PRIVATE
   !
   ! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Chem_Allocate
   PUBLIC :: Find_Number_of_Species
   PUBLIC :: Find_Index_of_Species
   PUBLIC :: Find_SeaSalt_Bin
   PUBLIC :: FindSpecByName
   PUBLIC :: GetSpecConc
   PUBLIC :: GetSpecConcByName
   PUBLIC :: GetSpecConcByIndex
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
      INTEGER              :: nSpecies          !< Total Number of Species
      INTEGER              :: nSpeciesGas       !< Number of Gas Species
      INTEGER              :: nSpeciesAero      !< Number of Aerosol Species
      INTEGER              :: nSpeciesAeroDryDep !< Number of Aerosol Species for Dry Dep
      INTEGER              :: nSpeciesDryDep     !< Number of all Species for Dry Dep
      INTEGER              :: nSpeciesWetDep     !< Number of all Species for Wet Dep
      INTEGER              :: nSpeciesTracer    !< Number of Tracer Species
      INTEGER              :: nSpeciesDust      !< Number of Dust Species
      INTEGER              :: nSpeciesSeaSalt   !< Number of SeaSalt Species
      INTEGER              :: nSpeciesBin       !< Number of SeaSalt Species Bin
      INTEGER, ALLOCATABLE :: SpeciesIndex(:)   !< Total Species Index
      INTEGER, ALLOCATABLE :: TracerIndex(:)    !< Tracer Species Index
      INTEGER, ALLOCATABLE :: AeroIndex(:)      !< Aerosol Species Index
      INTEGER, ALLOCATABLE :: GasIndex(:)       !< Gas Species Index
      INTEGER, ALLOCATABLE :: DustIndex(:)      !< Dust Species Index
      INTEGER, ALLOCATABLE :: SeaSaltIndex(:)   !< SeaSalt Species Index
      real(fp),ALLOCATABLE :: SeaSaltBinLower(:)  !< SeaSalt Species Bin Lower edge
      real(fp),ALLOCATABLE :: SeaSaltBinUpper(:) !< SeaSalt Species Bin upper edge
      INTEGER, ALLOCATABLE :: AeroDryDepIndex(:) !< Aerosol DryDep Species Index for Dry Dep
      INTEGER, ALLOCATABLE :: DryDepIndex(:)   !< All DryDep Species Index
      INTEGER, ALLOCATABLE :: WetDepIndex(:)   !< All WetDep Species Index
      CHARACTER(len=50), ALLOCATABLE :: SpeciesNames(:)  !< Species Names

      !---------------------------------------------------------------------
      ! Reals
      !---------------------------------------------------------------------
      type(SpeciesType), allocatable :: ChemSpecies(:)

   end type ChemStateType

CONTAINS


   !> \brief Allocate the chem state
   !!
   !! \details Allocate the chem state.
   !!
   !! \ingroup core_modules
   !!
   !! \param GridState Grid State
   !! \param ChemState Chem State
   !! \param RC Return code
   !!
   !!!>
   subroutine Chem_Allocate(GridState, ChemState, RC)

      ! USES
      USE GridState_Mod,  ONLY : GridStateType
      USE Species_Mod,    Only : SpeciesType

      IMPLICIT NONE

      ! INOUT Params
      type(GridStateType), INTENT(in)    :: GridState ! Grid State object
      !type(MetStateType), INTENT(in)    :: MetState   ! Met State object
      type(ChemStateType), INTENT(inout) :: ChemState ! chem State object
      ! type(SpeciesType),   POINTER       :: Species   ! Species object
      ! OUTPUT Params
      INTEGER,             INTENT(OUT)   :: RC            ! Success or failure

      ! Local
      integer :: i    ! Looping Variable

      ! Error handling
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: thisLoc

      ! Initialize
      RC = CC_SUCCESS
      ErrMsg = ''
      thisLoc = ' -> at chem_Allocate (in core/chemstate_mod.F90)'

      ! Allocate
      ALLOCATE( ChemState%ChemSpecies( ChemState%nSpecies ), STAT=RC )
      IF ( RC /= CC_SUCCESS ) THEN
         ErrMsg = 'Could not allocate ChemState%chemSpecies'
         CALL CC_Error( ErrMsg, RC, thisLoc )
      ENDIF
      CALL CC_CheckVar( 'ChemState%chemSpecies', 0, RC )
      IF ( RC /= CC_SUCCESS ) RETURN

      do i=0, ChemState%nSpecies
         ALLOCATE(ChemState%ChemSpecies(i)%conc(GridState%number_of_levels), STAT=RC)
         IF ( RC /= CC_SUCCESS ) THEN
            ErrMsg = 'Could not Allocate ChemState%ChemSpecies(i)%conc'
            CALL CC_Error( ErrMsg, RC, thisLoc )
         ENDIF
         ChemState%ChemSpecies(i)%conc = TINY_
      end do

   end subroutine Chem_Allocate

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
         if ( (ChemState%ChemSpecies(i)%is_drydep .eqv. .true.) .and. &
            (ChemState%ChemSpecies(i)%is_aerosol .eqv. .true.)   ) then
            ChemState%nSpeciesAeroDryDep = ChemState%nSpeciesAeroDryDep + 1
         endif
         if (ChemState%ChemSpecies(i)%is_drydep .eqv. .true.) then
            ChemState%nSpeciesDryDep = ChemState%nSpeciesDryDep + 1
         endif
         if (ChemState%ChemSpecies(i)%is_wetdep .eqv. .true.) then
            ChemState%nSpeciesWetDep = ChemState%nSpeciesWetDep + 1
         endif
      enddo

   end subroutine Find_Number_of_Species

   !> \brief Find the indices of species
   !!
   !! \param ChemState The ChemState object
   !! \param RC The return code
   !!
   !! \ingroup core_modules
   !!!>
   subroutine Find_Index_of_Species(ChemState, RC)
      ! USES
      USE Species_Mod,  ONLY : SpeciesType

      IMPLICIT NONE

      ! INOUT Params
      type(ChemStateType),  INTENT(INOUT) :: ChemState     ! chem State object
      ! OUTPUT Params
      INTEGER,             INTENT(OUT)   :: RC            ! Success or failure

      ! Error handling
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: thisLoc

      ! Local variables
      integer :: n ! looping variable
      integer :: aero_index      ! Current Aerosol Index
      integer :: gas_index       ! Current Gas Index
      integer :: dust_index      ! Current Dust Index
      integer :: seasalt_index   ! Current Seas Salt Index
      integer :: tracer_index    ! Current Tracer Index
      integer :: aero_drydep_index ! Current Aerosol DryDep Index
      integer :: drydep_index    ! Current DryDep Index
      integer :: wetdep_index    ! Current WetDep Index


      ! Initialize
      RC = CC_SUCCESS
      ErrMsg = ''
      thisLoc = ' -> at Find_indices_of_Species (in core/chemstate_mod.F90)'


      ! Initialize to zero before counting species
      aero_index = 1
      gas_index = 1
      dust_index = 1
      seasalt_index = 1
      tracer_index = 1
      aero_drydep_index = 1
      drydep_index = 1
      wetdep_index = 1

      ! Allocate index arrays
      ALLOCATE(Chemstate%AeroIndex(ChemState%nSpeciesAero), STAT=RC)
      IF ( RC /= CC_SUCCESS ) THEN
         errMsg = 'Error allocating Chemstate%AeroIndex'
         call CC_Error(errMsg, RC, thisLoc)
         RETURN
      ENDIF

      ALLOCATE(Chemstate%TracerIndex(ChemState%nSpeciesTracer), STAT=RC)
      IF ( RC /= CC_SUCCESS ) THEN
         errMsg = 'Error allocating Chemstate%TracerIndex'
         call CC_Error(errMsg, RC, thisLoc)
         RETURN
      ENDIF

      ALLOCATE(Chemstate%GasIndex(ChemState%nSpeciesGas), STAT=RC)
      IF ( RC /= CC_SUCCESS ) THEN
         errMsg = 'Error allocating Chemstate%GasIndex'
         call CC_Error(errMsg, RC, thisLoc)
         RETURN
      ENDIF

      ALLOCATE(Chemstate%DustIndex(ChemState%nSpeciesDust), STAT=RC)
      IF ( RC /= CC_SUCCESS ) THEN
         errMsg = 'Error allocating Chemstate%DustIndex'
         call CC_Error(errMsg, RC, thisLoc)
         RETURN
      ENDIF

      ALLOCATE(Chemstate%SeaSaltIndex(ChemState%nSpeciesSeaSalt), STAT=RC)
      IF ( RC /= CC_SUCCESS ) THEN
         errMsg = 'Error allocating Chemstate%SeaSaltIndex'
         call CC_Error(errMsg, RC, thisLoc)
         RETURN
      ENDIF

      ALLOCATE(Chemstate%AeroDryDepIndex(ChemState%nSpeciesAeroDryDep), STAT=RC)
      IF ( RC /= CC_SUCCESS ) THEN
         errMsg = 'Error allocating Chemstate%AeroDryDepIndex'
         call CC_Error(errMsg, RC, thisLoc)
         RETURN
      ENDIF

      ALLOCATE(Chemstate%DryDepIndex(ChemState%nSpeciesDryDep), STAT=RC)
      IF ( RC /= CC_SUCCESS ) THEN
         errMsg = 'Error allocating Chemstate%DryDepIndex'
         call CC_Error(errMsg, RC, thisLoc)
         RETURN
      ENDIF

      ALLOCATE(Chemstate%WetDepIndex(ChemState%nSpeciesWetDep), STAT=RC)
      IF ( RC /= CC_SUCCESS ) THEN
         errMsg = 'Error allocating Chemstate%WetDepIndex'
         call CC_Error(errMsg, RC, thisLoc)
         RETURN
      ENDIF

      ! Find indices for species groups
      do n = 1, ChemState%nSpecies
         if (ChemState%ChemSpecies(n)%is_aerosol .eqv. .true.) then
            Chemstate%AeroIndex(aero_index) = n
            aero_index = aero_index + 1
         endif
         if (ChemState%ChemSpecies(n)%is_gas .eqv. .true.) then
            Chemstate%GasIndex(gas_index) = n
            gas_index = gas_index + 1
         endif
         if (ChemState%ChemSpecies(n)%is_dust .eqv. .true.) then
            Chemstate%DustIndex(dust_index) = n
            dust_index = dust_index + 1
         endif
         if (ChemState%ChemSpecies(n)%is_seasalt .eqv. .true.) then
            Chemstate%SeaSaltIndex(seasalt_index) = n
            seasalt_index = seasalt_index + 1
         endif
         if (ChemState%ChemSpecies(n)%is_tracer .eqv. .true.) then
            Chemstate%TracerIndex(tracer_index) = n
            tracer_index = tracer_index + 1
         endif
         if ( (ChemState%ChemSpecies(n)%is_drydep  .eqv. .true.) .and. &
            (ChemState%ChemSpecies(n)%is_aerosol .eqv. .true.)    ) then
            Chemstate%AeroDryDepIndex(aero_drydep_index) = n
            aero_drydep_index = aero_drydep_index + 1
         endif
         if (ChemState%ChemSpecies(n)%is_drydep .eqv. .true.) then
            Chemstate%DryDepIndex(drydep_index) = n
            drydep_index = drydep_index + 1
         endif
         if (ChemState%ChemSpecies(n)%is_wetdep .eqv. .true.) then
            Chemstate%WetDepIndex(wetdep_index) = n
            wetdep_index = wetdep_index + 1
         endif
      enddo

   end subroutine Find_index_of_Species

   !> \brief Find the bins of sea salt species
   !!
   !! \param ChemState The ChemState object
   !! \param RC The return code
   !!
   !! \ingroup core_modules
   !!!>
   subroutine Find_SeaSalt_Bin(ChemState, RC)
      ! USES
      !USE Species_Mod,  ONLY : SpeciesType

      IMPLICIT NONE

      ! INOUT Params
      type(ChemStateType),  INTENT(INOUT) :: ChemState     ! chem State object
      ! OUTPUT Params
      INTEGER,             INTENT(OUT)   :: RC            ! Success or failure

      ! Error handling
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: thisLoc

      ! Local variables
      integer   :: n                  ! looping variable
      integer   :: n_bin              ! number of sea salt bins
      !real(fp)  :: radius0           ! initial radius of sea salt bin
      real(fp)  :: lower_radius(10)   ! lower radius of sea salt bin holder
      real(fp)  :: upper_radius(10)   ! upper radius of sea salt bin holder
      logical   :: mask(10) = .FALSE. ! flag to for sorting bins by radius

      ! Initialize
      RC = CC_SUCCESS
      ErrMsg = ''
      thisLoc = ' -> at Find_SeaSalt_Bin (in core/chemstate_mod.F90)'


      ! Initialize to zero before counting species
      n_bin = 0
      !radius0 = 0.0_fp

      ! Find possible sea salt bins
      do n = 1, ChemState%nSpeciesSeaSalt
         if (n == 1) then
            n_bin = 1
            lower_radius(n_bin) = ChemState%ChemSpecies(Chemstate%SeaSaltIndex(n))%lower_radius
            upper_radius(n_bin) = ChemState%ChemSpecies(Chemstate%SeaSaltIndex(n))%upper_radius
         else
            if (  ALL( ABS(lower_radius(1:n_bin) - ChemState%ChemSpecies(Chemstate%SeaSaltIndex(n))%lower_radius) > 0.0_fp )) then
               n_bin = n_bin + 1
               lower_radius(n_bin) = ChemState%ChemSpecies(Chemstate%SeaSaltIndex(n))%lower_radius
               upper_radius(n_bin) = ChemState%ChemSpecies(Chemstate%SeaSaltIndex(n))%upper_radius
            endif
         endif
      enddo

      ! Allocate index arrays
      ALLOCATE(Chemstate%SeaSaltBinLower(n_bin), STAT=RC)
      IF ( RC /= CC_SUCCESS ) THEN
         errMsg = 'Error allocating Chemstate%SeaSaltBinLower'
         call CC_Error(errMsg, RC, thisLoc)
         RETURN
      ENDIF

      ALLOCATE(Chemstate%SeaSaltBinUpper(n_bin), STAT=RC)
      IF ( RC /= CC_SUCCESS ) THEN
         errMsg = 'Error allocating Chemstate%SeaSaltBinUpper'
         call CC_Error(errMsg, RC, thisLoc)
         RETURN
      ENDIF

      !sort bins by radius from low to high for lower_radius
      mask(1:n_bin) = .TRUE.
      do n = 1, n_bin
         Chemstate%SeaSaltBinLower(n) =  MINVAL(lower_radius,mask)
         mask(MINLOC(lower_radius,mask)) = .FALSE.
      enddo

      !sort bins by radius from low to high for upper_radius
      mask(1:n_bin) = .TRUE.
      do n = 1, n_bin
         Chemstate%SeaSaltBinUpper(n) =  MINVAL(upper_radius,mask)
         mask(MINLOC(upper_radius,mask)) = .FALSE.
      enddo

      !check if the bins are continuous
      do n = 1, n_bin-1
         if ( .not. rae(Chemstate%SeaSaltBinUpper(n), Chemstate%SeaSaltBinLower(n+1)) ) then
            errMsg = 'Sea Salt Bins are not continuous'
            call CC_Error(errMsg, RC, thisLoc)
            RETURN
         endif
      enddo

   end subroutine Find_SeaSalt_Bin

   !> \brief Find the species by name
   !!
   !! \param ChemState The ChemState object
   !! \param name The name of the species
   !! \param index The index of the species
   !! \param RC The return code
   !!
   !! \ingroup core_modules
   !!!>
   subroutine FindSpecByName(ChemState, name, index, RC)

      type(ChemStateType),  INTENT(INOUT) :: ChemState     ! chem State object
      character(len=*),    INTENT(in)    :: name
      integer,              INTENT(out)   :: index
      integer,              INTENT(out)   :: RC

      ! Error handling
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: thisLoc

      ! local variables
      integer :: n

      ! Initialize
      RC = CC_SUCCESS
      ErrMsg = ''
      thisLoc = ' -> at FindSpecByName (in core/chemstate_mod.F90)'

      index = 0
      do n = 1, ChemState%nSpecies
         if (TRIM(name) == TRIM(ChemState%SpeciesNames(n))) then
            index = n
            exit
         endif
      enddo

   end subroutine FindSpecByName


   !> \brief Get the concentration of a species
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
   subroutine GetSpecConc(ChemState, concentration, RC, index, name)

      type(ChemStateType),  INTENT(INOUT) :: ChemState     ! chem State object
      real(kind=fp), dimension(:), INTENT(out)   :: concentration
      integer,              INTENT(out)   :: RC
      integer, optional,    INTENT(inout)    :: index
      character(len=50), optional, INTENT(inout)    :: name

      ! Error handling
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: thisLoc

      ! Initialize
      RC = CC_SUCCESS
      ErrMsg = ''
      thisLoc = ' -> at GetSpecConc (in core/chemstate_mod.F90)'

      if (present(index)) then
         call GetSpecConcByIndex(ChemState, concentration, index, RC)
      elseif (present(name)) then
         call GetSpecConcByName(ChemState, concentration, name, RC)
      else
         RC = CC_FAILURE
      endif

      if (RC /= CC_SUCCESS) then
         errMsg = 'Error in GetSpecConc'
         call CC_Error(errMsg, RC, thisLoc)
         RETURN
      endif

   end subroutine GetSpecConc

   !> \brief Get the concentration of a species by index
   !!
   !! \param ChemState The ChemState object
   !! \param concentration The concentration of the species
   !! \param RC The return code
   !! \param index The index of the species
   !!
   !! \ingroup core_modules
   !!!>
   subroutine GetSpecConcByIndex(ChemState, concentration, index, RC)

      type(ChemStateType),  INTENT(INOUT) :: ChemState     ! chem State object
      real(kind=fp), dimension(:), INTENT(out)   :: concentration
      integer,              INTENT(in)    :: index
      integer,              INTENT(out)   :: RC

      ! Error handling
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: thisLoc

      ! Initialize
      RC = CC_SUCCESS
      ErrMsg = ''
      thisLoc = ' -> at GetSpecConcByIndex (in core/chemstate_mod.F90)'

      if (index < 1 .or. index > ChemState%nSpecies) then
         RC = CC_FAILURE
         errMsg = 'index out of bounds'
         call CC_Error(errMsg, RC, thisLoc)
         RETURN
      endif

      concentration = ChemState%ChemSpecies(index)%conc

   end subroutine GetSpecConcByIndex

   !> \brief Get the concentration of a species by name
   !!
   !! \param ChemState The ChemState object
   !! \param concentration The concentration of the species
   !! \param RC The return code
   !! \param name The name of the species
   !!
   !! \ingroup core_modules
   !!!>
   subroutine GetSpecConcByName(ChemState, concentration, name, RC)

      type(ChemStateType),  INTENT(INOUT) :: ChemState     ! chem State object
      real(kind=fp), dimension(:), INTENT(out)   :: concentration
      character(len=50),    INTENT(in)    :: name
      integer,              INTENT(out)   :: RC

      ! Locals
      integer :: index

      ! Error handling
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: thisLoc

      ! Initialize
      RC = CC_SUCCESS
      ErrMsg = ''
      thisLoc = ' -> at GetSpecConcByName (in core/chemstate_mod.F90)'

      call FindSpecByName(ChemState, name, index, RC)

      if (RC /= CC_SUCCESS) then
         errMsg = 'Error in GetSpecConcByName'
         call CC_Error(errMsg, RC, thisLoc)
         RETURN
      endif

      concentration = ChemState%ChemSpecies(index)%conc

   end subroutine GetSpecConcByName

end module ChemState_Mod
