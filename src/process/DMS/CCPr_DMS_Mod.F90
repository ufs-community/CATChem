!> \brief CCPR DMS state types
!!
!!
!!
!! \author Lacey Holland
!! \date 01/2025
!!!>
MODULE CCPR_DMS_mod
   USE Precision_mod
   USE Error_Mod
   USE DiagState_Mod, Only : DiagStateType
   USE MetState_Mod,  Only : MetStateType
   USE ChemState_Mod, Only : ChemStateType
   USE Config_Opt_Mod, Only : ConfigType

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CCPR_DMS_Init
   PUBLIC :: CCPR_DMS_Run
   PUBLIC :: CCPR_DMS_Finalize
   PUBLIC :: DMSStateType

   !> \brief DMSStateType
   !!
   !! DMSStateType is the process-specific derived type. It should hold all module
   !! variables and arrays that are required to compute the emissions.
   !! For instance, if the process relies on an input field read through the
   !! CATChem configuration file (e.g. MY_INPUT_FIELD), the data array pointer
   !! to that field should be listed within the instance and NOT outside of it.
   !! This ensures that the same process can be invoked in various instances,
   !! all of them potentially pointing to different data fields.
   !!
   !! \param Activate Activate Process (True/False)
   !! \param SchemeOpt Scheme Option
   !!!>

   TYPE :: DMSStateType
      LOGICAL                         :: Activate              ! Activate Process (True/False)
      INTEGER                         :: SchemeOpt              ! SchemeOption (True/False)

      ! Process Specific Parameters

      ! Namelist parameters for specific DMS goes here as well
      !=================================================================
      ! Module specific variables/arrays/data pointers come below
      !=================================================================

   END TYPE DMSStateType

CONTAINS

   !>
   !! \brief Initialize the CATChem DMS module
   !!
   !! \param Config       CATCHem configuration options
   !! \param DMSState   CATCHem PROCESS state
   !! \param ChemState         CATCHem chemical state
   !! \param RC               Error return code
   !!
   !!!>


   SUBROUTINE CCPR_DMS_Init( Config, DMSState, ChemState, RC )
      ! USE


      IMPLICIT NONE
      ! INPUT PARAMETERS
      !-----------------
      TYPE(ConfigType), POINTER       :: Config    ! Module options
      TYPE(ChemStateType), POINTER    :: ChemState ! Chemical state

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      TYPE(DMSStateType), POINTER  :: DMSState ! DMS state
      INTEGER,         INTENT(INOUT)    :: RC       ! Success or failure

      ! Error handling
      !---------------
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc

      ! LOCAL VARIABLES
      !----------------


      ! Put any local variables here

      !=================================================================
      ! CCPR_DMS_Init begins here!
      !=================================================================
      ErrMsg = ''
      ThisLoc = ' -> at CCPR_DMS_INIT (in process/DMSemissions/ccpr_dms_mod.F90)'

      ! First check if process is activated in config | if not don't allocate arrays or pointers
      if (Config%DMS_activate) then

         ! Activate Process
         !------------------
         DMSState%Activate = .true.


         ! Set scheme option
         !------------------
         DMSState%SchemeOpt = config%DMS_Scheme

      else

         DMSState%Activate = .false.

      endif

   end subroutine CCPR_DMS_Init

   !>
   !! \brief Run the DMS emission scheme
   !!
   !! \param [IN] MetState - The MetState object
   !! \param [INOUT] DiagState - The DiagState object
   !! \param [INOUT] DMSState - The DMSState object
   !! \param [INOUT] ChemState - The ChemState object
   !! \param [OUT] RC Return code
   !!!>
   SUBROUTINE CCPr_DMS_Run( MetState, DiagState, DMSState, ChemState, RC )

      ! USE
      USE constants, only : g0
      USE CCPr_Scheme_GOCART_DMS_Mod, only : CCPR_Scheme_GOCART_DMS

      IMPLICIT NONE
      ! INPUT PARAMETERS
      TYPE(MetStateType),  INTENT(IN) :: MetState       ! MetState Instance

      ! INPUT/OUTPUT PARAMETERS
      TYPE(DiagStateType), INTENT(INOUT)   :: DiagState       ! DiagState Instance
      TYPE(DMSStateType), INTENT(INOUT)    :: DMSState     ! DMSState Instance
      TYPE(ChemStateType),  INTENT(INOUT)  :: ChemState       ! ChemState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                                 ! Return Code


      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc
      INTEGER            :: NDMS
      REAL, dimension(:,:), pointer   :: dmso_conc   ! concentration of DMS
      REAL, dimension(:,:,:),pointer  :: DMS       ! DMS [kg kg-1]
      REAL, dimension(:,:,:),pointer  :: SU_emis   ! SU emissions, kg/m2/s

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_DMS_Run (in process/DMSemissions/ccpr_DMS_mod.F90)'

      ! Run the DMS Scheme
      !-------------------------
      if (DMSState%Activate) then
         ! Run the DMS Scheme
         !-------------------------
         if (DMSState%SchemeOpt == 1) then
            ! Run the DMS Scheme
            !-------------------------

            call CCPr_Scheme_GOCART_DMS(MetState%NLEVS, &
               MetState%TSTEP, &
               g0, &
               MetState%T, &
               MetState%U10M, &
               MetState%V10M, &
               MetState%LWI, &
               MetState%DELP, &
               dmso_conc, &
               dms, &
               SU_emis, &
               ndms, &
               RC)

         endif

      endif


   end subroutine CCPr_DMS_Run


   !>
   !! \brief Finalize the DMSemissions
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
      INTEGER, INTENT(OUT) :: RC                                  ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_DMS_Finalize (in process/DMSemissions/ccpr_DMS_mod.F90)'


   end subroutine CCPr_DMS_Finalize


END MODULE CCPR_DMS_Mod
