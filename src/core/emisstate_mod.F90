!>
!! \file emisstate_mod.F90
!! \brief Emission Species State
!!
!! \details The Emission Species States contain information about the emitted species and
!!          the emitted fluxes, scale factors, and speciation to concentration mappings.
!!
!! \ingroup core_modules
!!!>
module EmisState_Mod

   ! Uses

   USE Error_Mod
   USE Precision_Mod

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: Emis_Allocate
   PUBLIC :: Emis_Find_Chem_Map_Index
   PUBLIC :: EmisState_CleanUp
   PUBLIC :: Apply_Emis_to_Chem

   !> \brief Emission Species State
   !!
   !! \details The Emission Species States contain information about the emitted species and
   !!          the emitted fluxes, scale factors, and speciation to concentration mappings.
   !!
   !! \ingroup core_modules
   !!
   !! \param name Name of the species
   !! \param long_name Long name of the species
   !! \param units Units of the species
   !! \param nEmisMap Number of Emission mappings per emitted species
   !! \param EmisMapIndex Emission mapping to concentration index
   !! \param Scale Scale factor
   !! \param Flux Emission flux
   !!
   !!!>
   TYPE :: EmisSpeciesType
      ! Character
      character(len=10)               :: name             !< Name of the species
      character(len=50)               :: long_name        !< Long name of the species
      character(len=10)               :: units            !< Units of the species
      character(len=100), allocatable :: EmisMapName(:)   !< Emission mapping to concentration species

      ! Integers
      integer              :: nEmisMap         !< Number of Emission mappings per emitted species
      integer, ALLOCATABLE :: EmisMapIndex(:)  !< Emission mapping to concentration index
      integer              :: plumerise        !< Plumerise option (default = 0 for no plumerise) sofiev = 1 briggs = 2 simple linear weighted thickness = 3
      integer              :: EmisLayer        !< Emission layer
      integer              :: nPlmSrc         !< Plumerise number of sources

      ! Real
      real(fp), ALLOCATABLE :: Scale(:)        !< Scale factor
      real(fp), ALLOCATABLE :: Flux(:)         !< Emission flux
      real(fp)              :: EmisHeight      !< Emission Height [m] - Simple emission height or -1 for PBLH
      real(fp), ALLOCATABLE :: PlmSrcFlx(:)    !< Plumerise source emission flux [kg/m2/s]
      real(fp), ALLOCATABLE :: PlmRiseHgt(:)   !< Height of Plume rise [m]
      real(fp), ALLOCATABLE :: FRP(:)          !< Fire Radiative Power (W/m^2)
      real(fp), ALLOCATABLE :: STKDM(:)        !< Briggs stack diameter [m] (array of all point sources in grid cell)
      real(fp), ALLOCATABLE :: STKHT(:)        !< Briggs stack thickness [m] (array of all point sources in grid cell)
      real(fp), ALLOCATABLE :: STKTK(:)        !< Briggs stack thickness [m] (array of all point sources in grid cell)
      real(fp), ALLOCATABLE :: STKVE(:)        !< Briggs stack velocity [m/s] (array of all point sources in grid cell)

   END TYPE EmisSpeciesType

   !> \brief Emission Category State
   !!
   !! \details The Emission Category State contains all the Species types for a given emission
   !!          category. It includes the name of the category, the number of emitted species, and
   !!          the emitted species containers.
   !!
   !! \param name Name of the category
   !! \param nSpecies Number of emitted species
   !! \param Species Emitted Species container
   !!
   !! \ingroup core_modules
   !!!>
   TYPE :: EmisCategoryType

      ! Character
      character(len=20)    :: name                         !< Name of the species

      ! Integers
      integer              :: nSpecies                 !< Number of emitted species per process
      integer              :: nPlumerise               !< Number of Species undergoing a Plumerise

      ! Types
      type(EmisSpeciesType), ALLOCATABLE :: Species(:) !< Emitted species container

   END TYPE EmisCategoryType

   !>
   !! \brief Emission State
   !!
   !! \details The Emission State contains all the emission category containers for the simulation and
   !!          the number of emission categories.
   !!
   !! \ingroup core_modules
   !!
   !! \param State Name of the state
   !! \param nCats Number of emission categories
   !! \param Cats Emission category containers
   !!
   !!!>
   TYPE, PUBLIC :: EmisStateType

      ! Character
      CHARACTER(LEN=4) :: State = 'Emis' !< Name of this state
      CHARACTER(LEN=10), ALLOCATABLE :: TotEmisNames(:)         !< Name of the state

      ! Integers
      integer :: nCats         !< Number of emission categories
      integer :: nEmisTotal              !< Total number of emitted species
      integer :: nEmisTotalPlumerise     !< Total number of plume rise categories

      ! Types
      type(EmisSpeciesType), ALLOCATABLE :: TotSpecies(:) !< Emitted species container
      type(EmisCategoryType), ALLOCATABLE :: Cats(:) !< Emission categories container

   END TYPE EmisStateType

CONTAINS

   !> \brief Allocate the emission state
   !!
   !! \details Allocate the emission state.
   !!
   !! \ingroup core_modules
   !!
   !! \param GridState Grid State
   !! \param EmisState Emission State
   !! \param RC Return code
   !!
   !!!>
   subroutine Emis_Allocate(GridState, EmisState, RC)

      ! Uses
      USE GridState_Mod, ONLY : GridStateType

      IMPLICIT NONE

      !Input
      TYPE(GridStateType), INTENT(IN) :: GridState

      ! Input/Output
      TYPE(EmisStateType), INTENT(INOUT) :: EmisState

      ! Output
      INTEGER, INTENT(OUT) :: RC

      ! Locals
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: ThisLoc

      integer :: c ! Loop counter for emission Cats
      integer :: s ! Loop counter for emitted species
      ! integer :: nPlumes ! temporary variable for number of plumes

      ! Initialize return code
      RC = CC_SUCCESS

      ! Initialize local variables
      ErrMsg = ''
      ThisLoc = ' -> at Emis_Allocate (in core/emisstate_mod.F90)'

      if (EmisState%nCats > 0) then
         print*, 'Number of Categories = ', EmisState%nCats
         do c = 1, EmisState%nCats
            print*, 'Number of Species for CAT=', EmisState%Cats(c)%name, ' = ', EmisState%Cats(c)%nSpecies
            do s = 1, EmisState%Cats(c)%nSpecies
               print*, 'Allocating ', EmisState%Cats(c)%Species(s)%name

               ALLOCATE(EmisState%Cats(c)%Species(s)%Flux(GridState%number_of_levels), STAT=RC)
               if (RC /= CC_SUCCESS) then
                  ErrMsg = '  Error allocating "EmisState%Cats%Species%Flux"!'
                  call CC_Error(ErrMsg, RC, ThisLoc)
                  return
               endif
               print*, 'Flux allocated for ', EmisState%Cats(c)%Species(s)%name

               ALLOCATE(EmisState%Cats(c)%Species(s)%EmisMapIndex(EmisState%Cats(c)%Species(s)%nEmisMap), STAT=RC)
               if (RC /= CC_SUCCESS) then
                  ErrMsg = '  Error allocating "EmisState%Cats%Species%EmisMapIndex"!'
                  call CC_Error(ErrMsg, RC, ThisLoc)
                  return
               endif

               print*, '  EmisMapIndex allocated for ', EmisState%Cats(c)%Species(s)%name

               ! if (EmisState%Cats(c)%Species(s)%nPlmSrc > 0) then ! if there are plume sources
               !    nPlumes = EmisState%Cats(c)%Species(s)%nPlmSrc ! temporary variable for number of plumes
               !    ALLOCATE(EmisState%Cats(c)%Species(s)%PlmSrcFlx(nPlumes), STAT=RC)
               !    if (RC /= CC_SUCCESS) then
               !       ErrMsg = 'Error allocating "EmisState%Cats%Species%PlmSrcFlx"!'
               !       call CC_Error(ErrMsg, RC, ThisLoc)
               !       return
               !    endif
               !    print*, '  PlmSrcFlx allocated for ', EmisState%Cats(c)%Species(s)%name

               !    ALLOCATE(EmisState%Cats(c)%Species(s)%FRP(nPlumes), STAT=RC)
               !    if (RC /= CC_SUCCESS) then
               !       ErrMsg = 'Error allocating "EmisState%Cats%Species%FRP"!'
               !       call CC_Error(ErrMsg, RC, ThisLoc)
               !       return
               !    endif
               !    print*, '  FRP allocated for ', EmisState%Cats(c)%Species(s)%name

               !    ALLOCATE(EmisState%Cats(c)%Species(s)%PlmRiseHgt(nPlumes), STAT=RC)
               !    if (RC /= CC_SUCCESS) then
               !       ErrMsg = 'Error allocating "EmisState%Cats%Species%PlmRiseHgt"!'
               !       call CC_Error(ErrMsg, RC, ThisLoc)
               !       return
               !    endif
               !    print*, '  PlmRiseHgt allocated for ', EmisState%Cats(c)%Species(s)%name

               ! endif
            end do
         end do

         ! call TotEmisSpecies_Allocate(GridState, EmisState, RC)
         ! IF (RC /= CC_success) THEN
         !    ErrMsg = 'Error in "TotEmisSpecies_Allocate"'
         !    call CC_Error(ErrMsg, RC, ThisLoc)
         !    return
         ! ENDIF

      end if

   end subroutine Emis_Allocate

   ! !> \brief Allocate the total emitted species
   ! !!
   ! !! \details Allocate the total emitted species
   ! !!
   ! !! \ingroup core_modules
   ! !!
   ! !! \param GridState Grid State
   ! !! \param EmisState Emission State
   ! !! \param RC Return code
   ! !!
   ! !!!>
   ! subroutine TotEmisSpecies_Allocate(GridState, EmisState, RC)

   !    USE GridState_Mod, ONLY : GridStateType

   !    IMPLICIT NONE

   !    TYPE(GridStateType), INTENT(IN) :: GridState

   !    TYPE(EmisStateType), INTENT(INOUT) :: EmisState

   !    INTEGER, INTENT(OUT) :: RC

   !    integer :: s ! Loop counter for emitted species
   !    integer :: c ! Loop counter for emission Cats
   !    integer :: t ! loop counter for total emitted species

   !    ! integer :: nT ! total number of emitted species counter
   !    ! integer :: nS ! number of the same emitted species in different categories

   !    ! logical :: IsIn

   !    character(len=255) :: ErrMsg
   !    character(len=255) :: ThisLoc

   !    character(len=10) ::currName

   !    ! Initialize return code
   !    RC = CC_SUCCESS

   !    ! Initialize local variables
   !    ErrMsg = ''
   !    ThisLoc = ' -> at TotEmisSpecies_Allocate (in core/emisstate_mod.F90)'

   !    ALLOCATE(EmisState%TotEmisNames(0), STAT=RC)
   !    if (EmisState%nCats > 0) then
   !       EmisState%nEmisTotal = 0
   !       do c = 1, EmisState%nCats
   !          do s = 1, EmisState%Cats(c)%nSpecies
   !             currName = EmisState%Cats(c)%Species(s)%name
   !             if ( ANY( EmisState%TotEmisNames == TRIM(currName) ) .eqv. .False. ) THEN
   !                EmisState%nEmisTotal = EmisState%nEmisTotal + 1
   !                EmisState%TotEmisNames = [EmisState%TotEmisNames, currName]
   !             endif
   !          end do
   !       end do

   !       ALLOCATE(EmisState%TotSpecies(EmisState%nEmisTotal), STAT=RC)
   !       IF (RC /= CC_SUCCESS) THEN
   !          ErrMsg = 'Error allocating "EmisState%TotSpecies"!'
   !          call CC_Error(ErrMsg, RC, ThisLoc)
   !          return
   !       ENDIF

   !       do t = 1, EmisState%nEmisTotal
   !          ! Fill Total Species fluxes and properties
   !          EmisState%TotSpecies(t)%name = EmisState%TotEmisNames(s)
   !          EmisState%TotSpecies(t)%Flux(:) = 0. ! set all fluxes to zero
   !          currName = EmisState%TotEmisNames(t)
   !          do c = 1, EmisState%nCats
   !             do s = 1, EmisState%Cats(c)%nSpecies
   !                if (currName == EmisState%Cats(c)%Species(s)%name) then
   !                   EmisState%TotSpecies(t)%units = EmisState%Cats(c)%Species(s)%units
   !                   EmisState%TotSpecies(t)%long_name = EmisState%Cats(c)%Species(s)%long_name
   !                endif
   !             end do
   !          end do
   !       end do
   !    endif

   ! end subroutine TotEmisSpecies_Allocate

   !> \brief Find the index of the mapped species
   !!
   !! \details Loops through the emitted species and finds the index to all of the mapped emitted species to the chemical species
   !!
   !! \ingroup core_modules
   !!
   !! \param EmisState Emission State
   !! \param ChemState Chemical State
   !! \param RC Return code
   !!
   !!!>
   subroutine Emis_Find_Chem_Map_Index(EmisState, ChemState, RC)

      USE ChemState_Mod, ONLY : ChemStateType, FindSpecByName

      TYPE(EmisStateType), INTENT(INOUT) :: EmisState
      TYPE(ChemStateType), INTENT(INOUT) :: ChemState
      INTEGER, INTENT(OUT) :: RC

      ! Local
      integer :: c ! Loop counter for emission Cats
      integer :: s ! Loop counter for emitted species
      integer :: n ! Loop counter for mapped species

      integer :: index ! index of the chemical species species

      ! Error handling
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: ThisLoc

      ! Initialize
      RC = CC_SUCCESS
      ErrMsg = ''
      ThisLoc = ' -> at Emis_find_chem_map_indexs (in core/emisstate_mod.F90)'

      do c = 1, EmisState%nCats
         do s = 1, EmisState%Cats(c)%nSpecies
            do n = 1, EmisState%Cats(c)%Species(s)%nEmisMap
               call FindSpecByName(ChemState, EmisState%Cats(c)%Species(s)%EmisMapName(n), index, RC)
               if (RC /= CC_SUCCESS) then
                  ErrMsg = 'Error in find_species_by_name'
                  call CC_Error(ErrMsg, RC, ThisLoc)
                  return
               endif
               EmisState%Cats(c)%Species(s)%EmisMapIndex(n) = index
            enddo
         enddo
      enddo

   end subroutine Emis_find_chem_map_index

   !> \brief Apply emissions to the chemical state
   !!
   !! \details TODO: Apply emissions to the chemical state
   !! \ingroup core_modules
   !!
   !! \param EmisState Emission State
   !! \param ChemState Chemical State
   !! \param RC Return code
   !!
   !!!>
   subroutine Apply_Emis_to_Chem(EmisState, MetState, ChemState, RC, surface_only_flag)
      USE ChemState_Mod, ONLY : ChemStateType
      USE MetState_Mod, ONLY : MetStateType
      use constants, only : g0

      TYPE(EmisStateType), INTENT(IN)    :: EmisState
      TYPE(MetStateType),  INTENT(IN)    :: MetState
      TYPE(ChemStateType), INTENT(INOUT) :: ChemState
      LOGICAL, INTENT(IN), OPTIONAL :: surface_only_flag
      INTEGER, INTENT(OUT) :: RC

      ! Local
      integer :: c !< Loop counter for emission Cats
      integer :: s !< Loop counter for emitted species
      integer :: n !< Loop counter for mapped species
      integer :: k !< Loop counter for height levels
      real(kind=fp) :: scale  !< scaling factor
      integer :: index        !< index of the mapped chemical species
      real(kind=fp) :: cdt    !< time step
      real(kind=fp) :: dqa    !< change in concentration
      logical :: surface_only !< flag for surface only

      ! Error handling
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: ThisLoc

      ! Initialize
      RC = CC_SUCCESS
      ErrMsg = ''
      ThisLoc = ' -> at Apply_Emis_to_Chem (in core/emisstate_mod.F90)'

      cdt = MetState%TSTEP

      if (present(surface_only_flag)) then
         surface_only = surface_only_flag
      else
         surface_only = .false.
      endif

      cats: do c = 1, EmisState%nCats
         species: do s = 1, EmisState%Cats(c)%nSpecies
            mapping: do n = 1, EmisState%Cats(c)%Species(s)%nEmisMap
               ! get the emis mapping and scale
               scale = EmisState%Cats(c)%Species(s)%Scale(n)
               index = EmisState%Cats(c)%Species(s)%EmisMapIndex(n)

               associate( &
               ! current flux of the emitted species
                  emis => EmisState%Cats(c)%Species(s)%Flux, &
               ! current concentration of the species at `index`
                  conc => ChemState%ChemSpecies(index)%conc)

                  levs: do k = 1, MetState%NLEVS
                     dqa = emis(k) * scale * cdt * g0 / MetState%DELP(k)
                     conc(k) = conc(k) + dqa
                  end do levs

               end associate

            enddo mapping
         enddo species
      enddo cats

   end subroutine Apply_Emis_to_Chem

   !> \brief Clean up the emission state
   !!
   !! \ingroup core_modules
   !!
   !! \param EmisState Emission State
   !! \param RC Return code
   !!
   !!!>
   subroutine EmisState_CleanUp(EmisState, RC)

      TYPE(EmisStateType), INTENT(INOUT) :: EmisState
      INTEGER, INTENT(OUT) :: RC

      character(len=255) :: ErrMsg
      character(len=255) :: ThisLoc

      integer :: c ! Loop counter for emission Cats
      integer :: s ! Loop counter for emitted species

      ! Initialize return code
      RC = CC_SUCCESS

      ! Initialize local variables
      ErrMsg = ''
      ThisLoc = ' -> at EmisState_CleanUp (in core/emisstate_mod.F90)'

      ! Deallocate total variables
      if (allocated(EmisState%TotEmisNames)) deallocate(EmisState%TotEmisNames)
      do c = 1, EmisState%nEmisTotal
         if (allocated(EmisState%TotSpecies(c)%Flux)) deallocate(EmisState%TotSpecies(c)%Flux)
         if (allocated(EmisState%TotSpecies(c)%Flux)) deallocate(EmisState%TotSpecies(c)%Flux)
      end do

      ! Deallocate emission variables in each category
      cats: do c = 1, EmisState%nCats
         species: do s = 1, EmisState%Cats(c)%nSpecies
            if (allocated(EmisState%Cats(c)%Species(s)%Flux)) &
               deallocate(EmisState%Cats(c)%Species(s)%Flux)
            if (allocated(EmisState%Cats(c)%Species(s)%EmisMapIndex)) &
               deallocate(EmisState%Cats(c)%Species(s)%EmisMapIndex)
            if (allocated(EmisState%Cats(c)%Species(s)%Scale)) &
               deallocate(EmisState%Cats(c)%Species(s)%Scale)
            if (allocated(EmisState%Cats(c)%Species(s)%EmisMapName)) &
               deallocate(EmisState%Cats(c)%Species(s)%EmisMapName)
            if (allocated(EmisState%Cats(c)%Species(s)%PlmRiseHgt)) &
               deallocate(EmisState%Cats(c)%Species(s)%PlmRiseHgt)
            if (allocated(EmisState%Cats(c)%Species(s)%PlmSrcFlx)) &
               deallocate(EmisState%Cats(c)%Species(s)%PlmSrcFlx)
            if (allocated(EmisState%Cats(c)%Species(s)%FRP)) &
               deallocate(EmisState%Cats(c)%Species(s)%FRP)
            if (allocated(EmisState%Cats(c)%Species(s)%STKDM)) &
               deallocate(EmisState%Cats(c)%Species(s)%STKDM)
            if (allocated(EmisState%Cats(c)%Species(s)%STKHT)) &
               deallocate(EmisState%Cats(c)%Species(s)%STKHT)
            if (allocated(EmisState%Cats(c)%Species(s)%STKTK)) &
               deallocate(EmisState%Cats(c)%Species(s)%STKTK)
            if (allocated(EmisState%Cats(c)%Species(s)%STKVE)) &
               deallocate(EmisState%Cats(c)%Species(s)%STKVE)
         end do species
      end do cats
      if (allocated(EmisState%Cats)) deallocate(EmisState%Cats)
      if (allocated(EmisState%Cats)) deallocate(EmisState%Cats)

   end subroutine EmisState_CleanUp

END MODULE EmisState_Mod
