!> \brief Driver for CATChem DMS process
!!
!!\defgroup catchem_dms_process
!! The CATChem DMS Process group holds all the CATCHem DMS processes.
!!
!! \author Lacey Holland and Wei Li
!! \date 01/2025
!!!>
MODULE CCPR_DMS_mod
   USE Precision_mod
   USE Error_Mod
   USE DiagState_Mod, Only : DiagStateType
   USE MetState_Mod,  Only : MetStateType
   USE ChemState_Mod, Only : ChemStateType
   USE EmisState_Mod,  Only : EmisStateType
   USE GridState_Mod,  Only : GridStateType
   USE Config_Opt_Mod, Only : ConfigType

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CCPR_DMS_Init
   PUBLIC :: CCPR_DMS_Run
   PUBLIC :: CCPR_DMS_Finalize
   PUBLIC :: DMSStateType

   !> \brief DMSStateType
   !!
   !! \details Contains all the information needed to run DMS Process
   !!
   !! This type contains the following variables:
   !! \param Activate Activate Process (True/False)
   !! \param SchemeOpt Scheme Option
   !! \param nDMSSpecies Number of DMS species
   !! \param DMSSpeciesIndex Index of DMS species
   !! \param DMSSpeciesName Name of DMS species
   !! \param SpcIDs CATChem species IDs
   !! \param CatIndex Index of emission category in EmisState
   !! \param TotalEmission Total emission [kg/m^2/s]
   !! \param EmissionPerSpecies Emission per species [kg/m^2/s]
   !!
   !! \ingroup catchem_dms_process
   !!!>

   TYPE :: DMSStateType
      LOGICAL                         :: Activate              ! Activate Process (True/False)
      INTEGER                         :: SchemeOpt             ! SchemeOption
      integer                         :: nDMSSpecies           !< Number of DMS species
      integer, allocatable            :: DMSSpeciesIndex(:)    !< Index of DMS species
      character(len=31), allocatable  :: DMSSpeciesName(:)     !< name of DMS species
      integer, allocatable            :: SpcIDs(:)             !< CATChem species IDs
      integer                         :: CatIndex              !< Index of emission category in EmisState

      ! Process Specific Parameters
      real(fp)                        :: TotalEmission         !< Total emission [kg/m^2/s]
      real(fp), allocatable           :: EmissionPerSpecies(:) !< Emission per species [kg/m^2/s]

   END TYPE DMSStateType

CONTAINS

   !>
   !! \brief Initialize the CATChem DMS process
   !!
   !! \param Config     CATCHem configuration options
   !! \param DMSState   CATCHem DMS state
   !! \param ChemState  CATCHem chemical state
   !! \param RC         Error return code
   !!
   !!!>
   SUBROUTINE CCPR_DMS_Init( Config, DMSState, EmisState, RC )
      ! USE

      IMPLICIT NONE
      ! INPUT PARAMETERS
      !-----------------
      TYPE(ConfigType)       :: Config    ! Module options
      TYPE(EmisStateType)    :: EmisState ! Chemical state

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      TYPE(DMSStateType)     :: DMSState  ! DMS state
      INTEGER, INTENT(INOUT) :: RC        ! Success or failure

      ! Error handling
      !---------------
      CHARACTER(LEN=255)     :: ErrMsg
      CHARACTER(LEN=255)     :: ThisLoc

      ! LOCAL VARIABLES
      !----------------
      INTEGER  :: c

      !=================================================================
      ! CCPR_DMS_Init begins here!
      !=================================================================
      ErrMsg = ''
      ThisLoc = ' -> at CCPR_DMS_INIT (in process/DMS/ccpr_dms_mod.F90)'

      ! First check if process is activated in config
      if (Config%DMS_activate) then

         ! Activate Process
         !------------------
         DMSState%Activate = .true.

         ! Set scheme option
         !------------------
         DMSState%SchemeOpt = config%DMS_Scheme

         !Find DMS caterory index in EmisState for future use
         !--------------------------------------------
         do c = 1, EmisState%nCats
            if (EmisState%Cats(c)%name == 'DMSO') then
               DMSState%CatIndex = c
               exit
            endif
         end do

         ! Set number of species from EmisState
         !----------------------
         DMSState%nDMSSpecies = EmisState%Cats(DMSState%CatIndex)%nSpecies

         !------------------------------------
         ! Allocate emission species index
         ALLOCATE( DMSState%DMSSpeciesIndex(DMSState%nDMSSpecies), STAT=RC )
         CALL CC_CheckVar('DMSState%DMSSpeciesIndex', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
         DMSState%DMSSpeciesIndex = -1

         ! Allocate emission speceis names
         ALLOCATE( DMSState%DMSSpeciesName(DMSState%nDMSSpecies), STAT=RC )
         CALL CC_CheckVar('DMSState%DMSSpeciesName', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
         DMSState%DMSSpeciesName = ''

         ! Allocate CatChem species index
         ALLOCATE( DMSState%SpcIDs(DMSState%nDMSSpecies), STAT=RC )
         CALL CC_CheckVar('DMSState%SpcIDs', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
         DMSState%SpcIDs = -1

         ! Allocate emission flux
         ALLOCATE( DMSState%EmissionPerSpecies(DMSState%nDMSSpecies), STAT=RC )
         CALL CC_CheckVar('DMSState%EmissionPerSpecies', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
         DMSState%EmissionPerSpecies = ZERO

         !initialize total emissions
         DMSState%TotalEmission = ZERO

      else

         DMSState%Activate = .false.

      endif

   end subroutine CCPR_DMS_Init

   !>
   !! \brief Run the DMS emission scheme
   !!
   !! \param [IN] MetState - The MetState object
   !! \param [INOUT] DMSState - The DMSState object
   !! \param [INOUT] EmisState - The EmisState object
   !! \param [OUT] RC Return code
   !!!>
   SUBROUTINE CCPr_DMS_Run( MetState, DMSState, EmisState, RC )

      ! USE
      USE constants, only : g0
      USE CCPr_Scheme_GOCART_DMS_Mod, only : CCPR_Scheme_GOCART_DMS

      IMPLICIT NONE

      ! INPUT PARAMETERS
      TYPE(MetStateType),  INTENT(IN)      :: MetState     ! MetState Instance

      ! INPUT/OUTPUT PARAMETERS
      TYPE(DMSStateType), INTENT(INOUT)    :: DMSState     ! DMSState Instance
      TYPE(EmisStateType),  INTENT(INOUT)  :: EmisState    ! EmisState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT)                 :: RC           ! Return Code

      ! LOCAL VARIABLES
      INTEGER :: s
      CHARACTER(LEN=255) :: ErrMsg, thisLoc
      REAL, dimension(:,:,:),pointer  :: SU_emis   ! DMS emissions in kg/m2/s

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_DMS_Run (in process/DMS/ccpr_DMS_mod.F90)'

      ! If DMS is activated
      !-------------------------
      if (DMSState%Activate) then
         ! Run the DMS Scheme
         !-------------------------
         do s = 1, DMSState%nDMSSpecies !only one species for now

            if (DMSState%SchemeOpt == 1) then

               !Dimensions (lon,lat,nSpecies) to be consistent with GOCART;
               !we all have one since it is a column model with only DMS emission
               allocate(SU_emis(1,1,1)); SU_emis = ZERO

               call CCPr_Scheme_GOCART_DMS(MetState%NLEVS, &
                  MetState%TSTEP, &
                  g0, &
                  MetState%T, &
                  MetState%U10M, &
                  MetState%V10M, &
                  MetState%LWI, &
                  MetState%DELP, &
                  MetState%DMSO_CONC, &
                  SU_emis, &
                  RC)

               if (RC /= CC_SUCCESS) then
                  errMsg = 'Error in CCPr_Scheme_GOCART_DMS'
                  CALL CC_Error( errMsg, RC, thisLoc )
               endif

               !put it back to DMSState
               DMSState%DMSSpeciesIndex(s)  = s
               DMSState%DMSSpeciesName(s)   = EmisState%Cats(DMSState%CatIndex)%Species(s)%name
               DMSState%EmissionPerSpecies(s) = SU_emis(1,1,1)
               DMSState%TotalEmission = DMSState%TotalEmission + DMSState%EmissionPerSpecies(s)
               if (associated(SU_emis)) nullify(SU_emis)

            else
               errMsg =  'ERROR: Unknown DMS scheme option'
               RC = CC_FAILURE
               CALL CC_Error( errMsg, RC, thisLoc )
               return

            endif !DMS scheme option
         end do !for each species

      endif

   end subroutine CCPr_DMS_Run


   !>
   !! \brief Finalize the DMS emission scheme
   !!
   !! \param [INOUT] DMSState
   !! \param [OUT] RC Return code
   !!!>
   SUBROUTINE CCPr_DMS_Finalize( DMSState, RC )

      ! USE
      !----

      IMPLICIT NONE

      ! INPUT/OUTPUT PARAMETERS
      TYPE(DMSStateType), INTENT(INOUT) :: DMSState  ! DMSState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT)              :: RC        ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_DMS_Finalize (in process/DMS/ccpr_DMS_mod.F90)'

      !Deallocate DMSState
      IF (ALLOCATED(DMSState%DMSSpeciesIndex)) THEN
         DEALLOCATE(DMSState%DMSSpeciesIndex, STAT=RC)
         CALL CC_CheckVar('DMSState%DMSSpeciesIndex', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      ENDIF
      IF (ALLOCATED(DMSState%DMSSpeciesName)) THEN
         DEALLOCATE(DMSState%DMSSpeciesName, STAT=RC)
         CALL CC_CheckVar('DMSState%DMSSpeciesName', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      ENDIF
      IF (ALLOCATED(DMSState%SpcIDs)) THEN
         DEALLOCATE(DMSState%SpcIDs, STAT=RC)
         CALL CC_CheckVar('DMSState%SpcIDs', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      ENDIF
      IF (ALLOCATED(DMSState%EmissionPerSpecies)) THEN
         DEALLOCATE(DMSState%EmissionPerSpecies, STAT=RC)
         CALL CC_CheckVar('DMSState%EmissionPerSpecies', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      ENDIF

   end subroutine CCPr_DMS_Finalize


END MODULE CCPR_DMS_Mod
