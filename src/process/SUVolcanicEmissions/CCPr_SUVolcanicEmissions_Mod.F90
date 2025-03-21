!> \brief CCPR suvolcanicemissions state types
!!
!! \defgroup catchem_suvolcanicemissions_process
!!
!! \author Lacey Holland and Wei Li
!! \date 10/2024
!!!>
MODULE CCPR_SUVolcanicEmissions_mod
   USE Precision_mod
   USE Error_Mod
   USE DiagState_Mod, Only : DiagStateType
   USE MetState_Mod,  Only : MetStateType
   USE ChemState_Mod, Only : ChemStateType
   USE Config_Opt_Mod, Only : ConfigType
   USE EmisState_Mod, Only: EmisStateType

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CCPR_SUVolcanicEmissions_Init
   PUBLIC :: CCPR_SUVolcanicEmissions_Run
   PUBLIC :: CCPR_SUVolcanicEmissions_Finalize
   PUBLIC :: SUVolcanicStateType


   !> \brief SUVolcanicStateType
   !!
   !! SUVolcanicStateType is the process-specific derived type.
   !!
   !! \param Activate Activate Process (True/False)
   !! \param Scheme Scheme Option
   !! \param nSUVolcanicSpecies # of SUVolcanic species
   !! \param SUVolcanicSpeciesIndex Index of SUVolcanic species
   !! \param SUVolcanicSpeciesName Name of SUVolcanic species
   !! \param SpcIDs CATChem species IDs
   !! \param CatIndex Index of emission category in EmisState
   !! \param TotalEmission Total emission of all species at each level [kg/m^2/s]
   !! \param EmissionPerSpecies Emission per species at each level [kg/m^2/s]
   !! \param FileDir Input file directory for reading in emissions
   !!
   !! \ingroup core_modules
   !!!>
   TYPE :: SUVolcanicStateType

      ! Generic Variables for Every Process
      LOGICAL                         :: Activate              ! Activate Process (True/False)
      INTEGER                         :: SchemeOpt             ! Scheme Option (if there is only one SchemeOpt always = 1)
      integer                         :: nSUVolcanicSpecies           !< Number of SUVolcanic species
      integer, allocatable            :: SUVolcanicSpeciesIndex(:)    !< Index of SUVolcanic species
      character(len=31), allocatable  :: SUVolcanicSpeciesName(:)     !< name of SUVolcanic species
      integer, allocatable            :: SpcIDs(:)               !< CATChem species IDs
      integer                         :: CatIndex                !< Index of emission category in EmisState

      ! Process Specific Parameters
      real(fp), allocatable           :: TotalEmission(:)          !< Total emission of all species at each level [kg/m^2/s]
      real(fp), allocatable           :: EmissionPerSpecies(:,:)   !< Emission per species at each level          [kg/m^2/s]
      character(len=1055)             :: FileDir                  !< Input file directory for reading in emissions


   END TYPE SUVolcanicStateType


CONTAINS

   !>
   !! \brief Initialize the CATChem Volcanic module
   !!
   !! \param Config       CATCHem configuration options
   !! \param SUVolcanicState   CATCHem PROCESS state
   !! \param EmisState         CATCHem emission state
   !! \param RC               Error return code
   !!
   !! \ingroup catchem_suvolcanicemissions_process
   !!
   !!!>
   SUBROUTINE CCPR_SUVolcanicEmissions_Init( Config, SUVolcanicState, EmisState, RC )
      ! USE


      IMPLICIT NONE
      ! INPUT PARAMETERS
      !-----------------
      TYPE(ConfigType)    :: Config    ! Module options
      TYPE(EmisStateType) :: EmisState ! Chemical state

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      TYPE(SUVolcanicStateType)    :: SUVolcanicState ! Volcanic state
      INTEGER,         INTENT(INOUT) :: RC       ! Success or failure

      ! Error handling
      !---------------
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc

      ! LOCAL VARIABLES
      !----------------
      INTEGER  :: c

      !=================================================================
      ! CCPR_DryDep_Init begins here!
      !=================================================================
      ErrMsg = ''
      ThisLoc = ' -> at CCPR_SUVolcanicEmissions_INIT (in process/SUVolcanicEmissions/CCPr_SUVolcanicEmissions_mod.F90)'

      ! First check if process is activated in config | if not don't allocate arrays or pointers
      if (Config%suvolcanic_activate) then

         ! Activate Process
         !------------------
         SUVolcanicState%Activate = .true.

         ! Set scheme option
         !------------------
         ! For now, the only option is SchemeOpt = 1
         SUVolcanicState%SchemeOpt = Config%suvolcanic_scheme

         !Find SUVOLCANIC caterory index in EmisState for future use
         !--------------------------------------------
         do c = 1, EmisState%nCats
            if (EmisState%Cats(c)%name == 'SUVOLCANIC') then
               SUVolcanicState%CatIndex = c
               exit
            endif
         end do

         ! Set number of species from EmisState
         !----------------------
         SUVolcanicState%nSUVolcanicSpecies = EmisState%Cats(SUVolcanicState%CatIndex)%nSpecies

         !------------------------------------
         ! Allocate emission species index
         ALLOCATE( SUVolcanicState%SUVolcanicSpeciesIndex(SUVolcanicState%nSUVolcanicSpecies), STAT=RC )
         CALL CC_CheckVar('SUVolcanicState%SUVolcanicSpeciesIndex', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
         SUVolcanicState%SUVolcanicSpeciesIndex = -1

         ! Allocate emission speceis names
         ALLOCATE( SUVolcanicState%SUVolcanicSpeciesName(SUVolcanicState%nSUVolcanicSpecies), STAT=RC )
         CALL CC_CheckVar('SUVolcanicState%SUVolcanicSpeciesName', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
         SUVolcanicState%SUVolcanicSpeciesName = ''

         ! Allocate CatChem species index
         ALLOCATE( SUVolcanicState%SpcIDs(SUVolcanicState%nSUVolcanicSpecies), STAT=RC )
         CALL CC_CheckVar('SUVolcanicState%SpcIDs', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
         SUVolcanicState%SpcIDs = -1

         ! Allocate emission flux
         ALLOCATE( SUVolcanicState%EmissionPerSpecies(SUVolcanicState%nSUVolcanicSpecies, &
            SIZE(EmisState%Cats(SUVolcanicState%CatIndex)%Species(1)%Flux)), STAT=RC )
         CALL CC_CheckVar('SUVolcanicState%EmissionPerSpecies', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
         SUVolcanicState%EmissionPerSpecies = ZERO

         ! Allocate total emissions
         ALLOCATE( SUVolcanicState%TotalEmission(SIZE(EmisState%Cats(SUVolcanicState%CatIndex)%Species(1)%Flux)), STAT=RC )
         CALL CC_CheckVar('SUVolcanicState%TotalEmission', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
         SUVolcanicState%TotalEmission = ZERO

         ! Set the file directory
         SUVolcanicState%FileDir = TRIM(Config%suvolcanic_filedir)

      else
         SUVolcanicState%Activate = .false.
      end if

   end subroutine CCPR_SUVolcanicEmissions_Init

   !>
   !! \brief Run the SUVolcanicEmissions
   !!
   !! \param [IN] MetState - The MetState object
   !! \param [INOUT] SUVolcanicState - The SUVolcanicState object
   !! \param [INOUT] EmisState - The EmisState object
   !! \param [OUT] RC Return code
   !!
   !! \ingroup catchem_suvolcanicemissions_process
   !!!>
   SUBROUTINE CCPr_SUVolcanicEmissions_Run( MetState, SUVolcanicState, EmisState, RC )

      ! USE
      USE constants, only : g0
      USE ReadEmissions, only:  ReadASCIIPointEmissions, VolcanicEmissionData
      use CCPr_Scheme_GOCART_SUVolcanicEmissions_Mod, only : CCPr_Scheme_GOCART_SUVolcanicEmissions

      IMPLICIT NONE
      ! INPUT PARAMETERS
      TYPE(MetStateType),  INTENT(IN) :: MetState       !< MetState Instance

      ! INPUT/OUTPUT PARAMETERS
      !TYPE(DiagStateType), INTENT(INOUT)      :: DiagState       !< DiagState Instance
      TYPE(SUVolcanicStateType), INTENT(INOUT)  :: SUVolcanicState  !< SUVolcanicState Instance
      !TYPE(ChemStateType), INTENT(INOUT)     :: ChemState       !< ChemState Instance
      TYPE(EmisStateType), INTENT(INOUT)     :: EmisState       !< ChemState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                                 ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc
      CHARACTER(LEN=1055) :: fname
      type(VolcanicEmissionData), allocatable :: VolcanicEmis(:)
      CHARACTER(len=7), parameter :: label='volcano'
      INTEGER :: i                          ! loop index
      INTEGER :: hms                        ! Model time [secs] TODO: format is right?
      INTEGER :: ymd                        ! Model date [YYYYMMDD] TODO: format is right?
      CHARACTER(LEN=18) :: ymd_str          ! Model date in string

      INTEGER, dimension(:), allocatable   :: vStart      ! Emissions Start time [sec]
      INTEGER, dimension(:), allocatable   :: vEnd        ! Emissions end time [sec]
      INTEGER                              :: nVolc       ! number of volcanic sources
      INTEGER, dimension(:), allocatable   :: iPoint, jPoint ! sub-domain - we only run this at the place/time of eruption??
      INTEGER, parameter                   :: nSO2 =1     ! index of SO2 relative to other sulfate tracers

      REAL, dimension(:), allocatable      :: vSO2   ! volcanic emissions  [kg]
      REAL, dimension(:,:,:),pointer  :: SO2       ! SO2 [kg kg-1]
      !REAL, dimension(:,:,:),pointer  :: SU_emis   ! SU emissions, kg/m2/s
      REAL, dimension(:), allocatable        :: vCloud    ! top elevation of emissions [m]
      REAL, dimension(:), allocatable        :: vElev     ! bottom elevation of emissions [m]
      REAL, dimension(:), allocatable        :: vLat     ! latitude specified in file [degree]
      REAL, dimension(:), allocatable        :: VLon     ! longitude specified in file [degree]
      REAL, dimension(:,:), allocatable      :: area     ! area of current grid cell [m^2]

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_SUVolcanicEmissions_Run &
      & (in process/SUVolcanicEmissions/ccpr_SUVolcanicEmissions_mod.F90)'

      ! Run  SUVolcanic
      !-------------------------
      if (SUVolcanicState%Activate) then
         ! Run the GOCART SUVolcanic Scheme
         !-------------------------
         if (SUVolcanicState%SchemeOpt == 1) then
            ! Run the SU Volcanic GOCART Scheme
            !-------------------------
            if (SUVolcanicState%nSUVolcanicSpecies  > 0) then

               !TODO: read file name from config file based on ymd?
               ymd = MetState%YMD; hms = MetState%HMS
               write(ymd_str, "(i0)") ymd ! converting integer to string
               fname = TRIM(SUVolcanicState%FileDir) // '.' // TRIM(ymd_str) // '.rc'
               call ReadASCIIPointEmissions (fname, label, VolcanicEmis, RC)
               nVolc = VolcanicEmis(1)%nPts
               allocate(vSO2(nVolc), vCloud(nVolc), vElev(nVolc), vLat(nVolc), VLon(nVolc))
               vSO2 = VolcanicEmis(:)%VEmis
               vCloud = VolcanicEmis(:)%Vtop
               vElev = VolcanicEmis(:)%Vbase
               vLon = VolcanicEmis(:)%Vlon
               vLat = VolcanicEmis(:)%Vlat

               !TODO: iPoint and jPoint needs to be determined; currently gives all one. In real run,
               !we could pre-select sources for the current grid cell so giving all ones can work fine.
               allocate(iPoint(nVolc), jPoint(nVolc))
               iPoint = 1; jPoint = 1
               !TODO: Not sure how to get Vstart and VEnd format and values.
               allocate(vStart(nVolc), vEnd(nVolc))
               vStart = ymd  + 000000
               vEnd =   ymd  + 235959
               hms = ymd + hms
               !set area and SO2
               allocate(area(1,1), SO2(1,1, MetState%NLEVS))
               area(1,1) = MetState%AREA_M2
               SO2 = ZERO

               ! loop through all species. Right now, GOCART only has SO2
               do i = 1, SUVolcanicState%nSUVolcanicSpecies

                  !Need to look up which is the index for SO2 concentrations
                  !TODO: is level index reversed in the GOCART???
                  call CCPr_Scheme_GOCART_SUVolcanicEmissions( MetState%NLEVS,   &
                     MetState%TSTEP, &
                     VStart, &
                     VEnd, &
                     nVolc, &
                     iPoint, &
                     jPoint, &
                     hms, &
                     g0, &
                     MetState%BXHEIGHT, &
                     MetState%DELP, &
                     area, &
                     vSO2, &    !volcanic contribution to so2 emissions
                     nSO2, &    !tracer number for so2 within sulfur trace
                     SO2, &     !total so2 concentration intent(inout)
                  !SU_emis, & !total emission rate for each sulfur species, !SU_emis(:,:,nSO2), nSO2=2, nDMS=1, nSO4=3, nMSA=4
                     vCloud, &
                     vElev, &
                     vLat, &
                     VLon, &
                     RC)
                  !  nso2 is used to define SU_emis:  SU_emis(:,:,nSO2). We only use SO2 for now and assign nSO2=1

                  !put it back to SUVolcanicState
                  SUVolcanicState%SUVolcanicSpeciesIndex(i) = i
                  SUVolcanicState%SUVolcanicSpeciesName(i) = EmisState%Cats(SUVolcanicState%CatIndex)%Species(i)%name
                  !TODO: convert unit from kg kg-1 to kg m-2 s-1;
                  !TODO: The test only has 8 levels and we give EmissionPerSpecies 28 levels from GridState. In real run,
                  ! "1:8" should be changed to ":" for all the vertical levels.
                  SUVolcanicState%EmissionPerSpecies(i,1:8) = SO2(1, 1, :) * MetState%DELP / g0 / MetState%TSTEP
                  SUVolcanicState%TotalEmission(:) = SUVolcanicState%TotalEmission(:) + SUVolcanicState%EmissionPerSpecies(i,:)

               end do ! do i = 1, SUVolcanicState%nSUVolcanicSpecies

            endif  ! if (SUVolcanicState%nSUVolcanicSpecies  > 0)

         endif  ! if (SUVolcanicState%SchemeOpt == 1)

         write(*,*) 'TODO: Need to figure out how to add back to the chemical species state '

      endif   !  if (SUVolcanicState%Activate)

   end subroutine CCPr_SUVolcanicEmissions_Run

   !>
   !! \brief Finalize the DryDep
   !!
   !! \param [INOUT] SUVolcanicState
   !! \param [OUT] RC Return code
   !!!>
   SUBROUTINE CCPr_SUVolcanicEmissions_Finalize( SUVolcanicState, RC )

      ! USE
      !----

      IMPLICIT NONE

      ! INPUT/OUTPUT PARAMETERS
      TYPE(SUVolcanicStateType), INTENT(INOUT) :: SUVolcanicState  ! SUVolcanicState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                                  ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_SUVolcanicEmissions_Finalize &
      &(in process/SUVolcanicEmissions/ccpr_SUVolcanicEmissions_mod.F90)'

      ! Deallocate any arrays here
      IF (ALLOCATED(SUVolcanicState%SpcIDs)) THEN
         DEALLOCATE( SUVolcanicState%SpcIDs, STAT=RC )
         CALL CC_CheckVar('SUVolcanicState%SpcIDs', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      END IF

      IF (ALLOCATED(SUVolcanicState%SUVolcanicSpeciesIndex)) THEN
         DEALLOCATE( SUVolcanicState%SUVolcanicSpeciesIndex, STAT=RC )
         CALL CC_CheckVar('SUVolcanicState%SUVolcanicSpeciesIndex', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      END IF
      
      IF (ALLOCATED(SUVolcanicState%SUVolcanicSpeciesName)) THEN
         DEALLOCATE( SUVolcanicState%SUVolcanicSpeciesName, STAT=RC )
         CALL CC_CheckVar('SUVolcanicState%SUVolcanicSpeciesName', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      END IF
      
      IF (ALLOCATED(SUVolcanicState%EmissionPerSpecies)) THEN
         DEALLOCATE( SUVolcanicState%EmissionPerSpecies, STAT=RC )
         CALL CC_CheckVar('SUVolcanicState%EmissionPerSpecies', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      END IF

      IF (ALLOCATED(SUVolcanicState%TotalEmission)) THEN
         DEALLOCATE( SUVolcanicState%TotalEmission, STAT=RC )
         CALL CC_CheckVar('SUVolcanicState%TotalEmission', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      END IF
   
   end subroutine CCPr_SUVolcanicEmissions_Finalize


END MODULE CCPR_SUVolcanicEmissions_Mod
