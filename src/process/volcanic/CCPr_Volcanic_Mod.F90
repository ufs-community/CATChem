!> \brief CCPR Volcanice state types
!!
!! \defgroup catchem_Volcanic_process
!!
!! \author Lacey Holland and Wei Li
!! \date 10/2024
!!!>
MODULE CCPR_Volcanic_mod
   USE Precision_mod
   USE Error_Mod
   USE DiagState_Mod, Only : DiagStateType
   USE MetState_Mod,  Only : MetStateType
   USE ChemState_Mod, Only : ChemStateType
   USE Config_Opt_Mod, Only : ConfigType
   USE EmisState_Mod, Only: EmisStateType

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CCPR_Volcanic_Init
   PUBLIC :: CCPR_Volcanic_Run
   PUBLIC :: CCPR_Volcanic_Finalize
   PUBLIC :: VolcanicStateType


   !> \brief VolcanicStateType
   !!
   !! VolcanicStateType is the process-specific derived type.
   !!
   !! \param Activate Activate Process (True/False)
   !! \param Scheme Scheme Option
   !! \param nVolcanicSpecies # of Volcanic species
   !! \param VolcanicSpeciesIndex Index of Volcanic species
   !! \param VolcanicSpeciesName Name of Volcanic species
   !! \param SpcIDs CATChem species IDs
   !! \param CatIndex Index of emission category in EmisState
   !! \param TotalEmission Total emission of all species at each level [kg/m^2/s]
   !! \param EmissionPerSpecies Emission per species at each level [kg/m^2/s]
   !! \param FileDir Input file directory for reading in emissions
   !!
   !! \ingroup catchem_Volcanic_process
   !!!>
   TYPE :: VolcanicStateType

      ! Generic Variables for Every Process
      LOGICAL                         :: Activate              ! Activate Process (True/False)
      INTEGER                         :: SchemeOpt             ! Scheme Option (if there is only one SchemeOpt always = 1)
      integer                         :: nVolcanicSpecies           !< Number of Volcanic species
      integer, allocatable            :: VolcanicSpeciesIndex(:)    !< Index of Volcanic species
      character(len=31), allocatable  :: VolcanicSpeciesName(:)     !< name of Volcanic species
      integer, allocatable            :: SpcIDs(:)               !< CATChem species IDs
      integer                         :: CatIndex                !< Index of emission category in EmisState

      ! Process Specific Parameters
      real(fp), allocatable           :: TotalEmission(:)          !< Total emission of all species at each level [kg/m^2/s]
      real(fp), allocatable           :: EmissionPerSpecies(:,:)   !< Emission per species at each level          [kg/m^2/s]
      character(len=1055)             :: FileDir                  !< Input file directory for reading in emissions


   END TYPE VolcanicStateType


CONTAINS

   !>
   !! \brief Initialize the CATChem Volcanic module
   !!
   !! \param Config       CATCHem configuration options
   !! \param VolcanicState   CATCHem PROCESS state
   !! \param EmisState         CATCHem emission state
   !! \param RC               Error return code
   !!
   !! \ingroup catchem_Volcanicemissions_process
   !!
   !!!>
   SUBROUTINE CCPR_Volcanic_Init( Config, VolcanicState, EmisState, RC )
      ! USE


      IMPLICIT NONE
      ! INPUT PARAMETERS
      !-----------------
      TYPE(ConfigType)    :: Config    ! Module options
      TYPE(EmisStateType) :: EmisState ! Chemical state

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      TYPE(VolcanicStateType)    :: VolcanicState ! Volcanic state
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
      ThisLoc = ' -> at CCPR_Volcanic_INIT (in process/Volcanic/CCPr_Volcanic_mod.F90)'

      ! First check if process is activated in config | if not don't allocate arrays or pointers
      if (Config%volcanic_activate) then

         ! Activate Process
         !------------------
         VolcanicState%Activate = .true.

         ! Set scheme option
         !------------------
         ! For now, the only option is SchemeOpt = 1
         VolcanicState%SchemeOpt = Config%volcanic_scheme

         !Find VOLCANIC caterory index in EmisState for future use
         !--------------------------------------------
         do c = 1, EmisState%nCats
            if (EmisState%Cats(c)%name == 'VOLCANIC') then
               VolcanicState%CatIndex = c
               exit
            endif
         end do

         ! Set number of species from EmisState
         !----------------------
         VolcanicState%nVolcanicSpecies = EmisState%Cats(VolcanicState%CatIndex)%nSpecies

         !------------------------------------
         ! Allocate emission species index
         ALLOCATE( VolcanicState%VolcanicSpeciesIndex(VolcanicState%nVolcanicSpecies), STAT=RC )
         CALL CC_CheckVar('VolcanicState%VolcanicSpeciesIndex', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
         VolcanicState%VolcanicSpeciesIndex = -1

         ! Allocate emission speceis names
         ALLOCATE( VolcanicState%VolcanicSpeciesName(VolcanicState%nVolcanicSpecies), STAT=RC )
         CALL CC_CheckVar('VolcanicState%VolcanicSpeciesName', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
         VolcanicState%VolcanicSpeciesName = ''

         ! Allocate CatChem species index
         ALLOCATE( VolcanicState%SpcIDs(VolcanicState%nVolcanicSpecies), STAT=RC )
         CALL CC_CheckVar('VolcanicState%SpcIDs', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
         VolcanicState%SpcIDs = -1

         ! Allocate emission flux
         ALLOCATE( VolcanicState%EmissionPerSpecies(VolcanicState%nVolcanicSpecies, &
            SIZE(EmisState%Cats(VolcanicState%CatIndex)%Species(1)%Flux)), STAT=RC )
         CALL CC_CheckVar('VolcanicState%EmissionPerSpecies', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
         VolcanicState%EmissionPerSpecies = ZERO

         ! Allocate total emissions
         ALLOCATE( VolcanicState%TotalEmission(SIZE(EmisState%Cats(VolcanicState%CatIndex)%Species(1)%Flux)), STAT=RC )
         CALL CC_CheckVar('VolcanicState%TotalEmission', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
         VolcanicState%TotalEmission = ZERO

         ! Set the file directory
         VolcanicState%FileDir = TRIM(Config%Volcanic_filedir)

      else
         VolcanicState%Activate = .false.
      end if

   end subroutine CCPR_Volcanic_Init

   !>
   !! \brief Run the VolcanicEmissions
   !!
   !! \param [IN] MetState - The MetState object
   !! \param [INOUT] VolcanicState - The VolcanicState object
   !! \param [INOUT] EmisState - The EmisState object
   !! \param [INOUT] RC Return code
   !!
   !! \ingroup catchem_Volcanicemissions_process
   !!!>
   SUBROUTINE CCPr_Volcanic_Run( MetState, VolcanicState, EmisState, RC )

      ! USE
      USE constants, only : g0
      use CCPr_Scheme_Volcanic_GOCART_Mod, only : CCPr_Scheme_Volcanic_GOCART, VolcanicEmisData, ReadASCIIPointEmis

      IMPLICIT NONE
      ! INPUT PARAMETERS
      TYPE(MetStateType),  INTENT(IN) :: MetState       !< MetState Instance

      ! INPUT/OUTPUT PARAMETERS
      !TYPE(DiagStateType), INTENT(INOUT)      :: DiagState       !< DiagState Instance
      TYPE(VolcanicStateType), INTENT(INOUT)  :: VolcanicState  !< VolcanicState Instance
      !TYPE(ChemStateType), INTENT(INOUT)     :: ChemState       !< ChemState Instance
      TYPE(EmisStateType), INTENT(INOUT)     :: EmisState       !< ChemState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(INOUT) :: RC                                 ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc
      CHARACTER(LEN=1055) :: fname
      type(VolcanicEmisData), allocatable :: VolcanicEmis(:)
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
      thisLoc = ' -> at CCPr_VolcanicEmissions_Run &
      & (in process/VolcanicEmissions/ccpr_VolcanicEmissions_mod.F90)'

      ! Run  Volcanic
      !-------------------------
      if (VolcanicState%Activate) then
         ! Run the GOCART Volcanic Scheme
         !-------------------------
         if (VolcanicState%SchemeOpt == 1) then
            ! Run the SU Volcanic GOCART Scheme
            !-------------------------
            if (VolcanicState%nVolcanicSpecies  > 0) then

               !TODO: read file name from config file based on ymd?
               ymd = MetState%YMD; hms = MetState%HMS
               write(ymd_str, "(i0)") ymd ! converting integer to string
               fname = TRIM(VolcanicState%FileDir) // '.' // TRIM(ymd_str) // '.rc'
               call ReadASCIIPointEmis (fname, label, VolcanicEmis, RC)
               nVolc = VolcanicEmis(1)%nPts
               allocate(vSO2(nVolc), vCloud(nVolc), vElev(nVolc), vLat(nVolc), VLon(nVolc))
               vSO2 = VolcanicEmis(:)%VEmis
               vCloud = VolcanicEmis(:)%Vtop
               vElev = VolcanicEmis(:)%Vbase
               vLon = VolcanicEmis(:)%Vlon
               vLat = VolcanicEmis(:)%Vlat
               if (allocated(VolcanicEmis)) deallocate(VolcanicEmis)

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
               do i = 1, VolcanicState%nVolcanicSpecies

                  !Need to look up which is the index for SO2 concentrations
                  !TODO: is level index reversed in the GOCART???
                  call CCPr_Scheme_Volcanic_GOCART( MetState%NLEVS,   &
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
                     vCloud, &
                     vElev, &
                     vLat, &
                     VLon, &
                     RC)
                  !  nso2 is used to define SU_emis:  SU_emis(:,:,nSO2). We only use SO2 for now and assign nSO2=1

                  !put it back to VolcanicState
                  VolcanicState%VolcanicSpeciesIndex(i) = i
                  VolcanicState%VolcanicSpeciesName(i) = EmisState%Cats(VolcanicState%CatIndex)%Species(i)%name
                  !TODO: convert unit from kg kg-1 to kg m-2 s-1;
                  !TODO: The test only has 8 levels and we give EmissionPerSpecies 28 levels from GridState. In real run,
                  ! "1:8" should be changed to ":" for all the vertical levels.
                  VolcanicState%EmissionPerSpecies(i,1:8) = SO2(1, 1, :) * MetState%DELP / g0 / MetState%TSTEP
                  VolcanicState%TotalEmission(:) = VolcanicState%TotalEmission(:) + VolcanicState%EmissionPerSpecies(i,:)

               end do ! do i = 1, VolcanicState%nVolcanicSpecies

            endif  ! if (VolcanicState%nVolcanicSpecies  > 0)

         endif  ! if (VolcanicState%SchemeOpt == 1)

         write(*,*) 'TODO: Need to figure out how to add back to the chemical species state '

      endif   !  if (VolcanicState%Activate)

   end subroutine CCPr_Volcanic_Run

   !>
   !! \brief Finalize the DryDep
   !!
   !! \param [INOUT] VolcanicState
   !! \param [INOUT] RC Return code
   !!!>
   SUBROUTINE CCPr_Volcanic_Finalize( VolcanicState, RC )

      ! USE
      !----

      IMPLICIT NONE

      ! INPUT/OUTPUT PARAMETERS
      TYPE(VolcanicStateType), INTENT(INOUT) :: VolcanicState  ! VolcanicState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(INOUT) :: RC                                  ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_VolcanicEmissions_Finalize &
      &(in process/VolcanicEmissions/ccpr_VolcanicEmissions_mod.F90)'

      ! Deallocate any arrays here
      IF (ALLOCATED(VolcanicState%SpcIDs)) THEN
         DEALLOCATE( VolcanicState%SpcIDs, STAT=RC )
         CALL CC_CheckVar('VolcanicState%SpcIDs', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      END IF

      IF (ALLOCATED(VolcanicState%VolcanicSpeciesIndex)) THEN
         DEALLOCATE( VolcanicState%VolcanicSpeciesIndex, STAT=RC )
         CALL CC_CheckVar('VolcanicState%VolcanicSpeciesIndex', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      END IF

      IF (ALLOCATED(VolcanicState%VolcanicSpeciesName)) THEN
         DEALLOCATE( VolcanicState%VolcanicSpeciesName, STAT=RC )
         CALL CC_CheckVar('VolcanicState%VolcanicSpeciesName', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      END IF

      IF (ALLOCATED(VolcanicState%EmissionPerSpecies)) THEN
         DEALLOCATE( VolcanicState%EmissionPerSpecies, STAT=RC )
         CALL CC_CheckVar('VolcanicState%EmissionPerSpecies', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      END IF

      IF (ALLOCATED(VolcanicState%TotalEmission)) THEN
         DEALLOCATE( VolcanicState%TotalEmission, STAT=RC )
         CALL CC_CheckVar('VolcanicState%TotalEmission', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      END IF

   end subroutine CCPr_Volcanic_Finalize


END MODULE CCPR_Volcanic_Mod
