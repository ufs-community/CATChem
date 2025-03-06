!> \brief Driver for the CATCHem Process: BVOC
!!
!!
!! \defgroup catchem_bvoc_process
!!
!! \author Wei Li
!! \date 07/2024
!!!>
MODULE CCPR_BVOC_mod
   USE Precision_mod, only : fp
   USE Error_Mod,   Only : CC_Error, CC_SUCCESS, CC_FAILURE, CC_CheckVar
   USE constants, only : PI_180
   USE DiagState_Mod, Only : DiagStateType
   USE MetState_Mod,  Only : MetStateType
   USE ChemState_Mod, Only : ChemStateType
   USE Config_Opt_Mod,    Only : ConfigType
   USE CCPr_BVOC_Common_Mod, Only : BvocStateType
   USE EmisState_Mod, Only : EmisStateType

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CCPR_BVOC_Init
   PUBLIC :: CCPR_BVOC_Run
   PUBLIC :: CCPR_BVOC_Final

CONTAINS

   !>
   !! \brief Initialize the CATChem BVOC module
   !!
   !! \param Config_Opt       CATCHem configuration options
   !! \param BvocState        CATCHem Bvoc state
   !! \param EmisState        CATCHem Emission state
   !! \param ChmState         CATCHem chemical state
   !! \param RC               Error return code
   !!
   !!!>
   !SUBROUTINE CCPR_BVOC_Init( Config, ChemState, EmisState, BvocState, RC )
   SUBROUTINE CCPR_BVOC_Init( Config, EmisState, BvocState, RC )
      ! USE

      IMPLICIT NONE
      ! INPUT PARAMETERS
      !-----------------
      TYPE(ConfigType),  intent(in)    :: Config        ! Module options
      !TYPE(ChemStateType),  intent(in)    :: ChemState  ! Chemical state
      TYPE(EmisStateType),  intent(in)    :: EmisState  ! Emission state

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      TYPE(BvocStateType), intent(inout) :: BvocState   ! Bvoc state
      INTEGER,              intent(inout) :: RC         ! Success or failure

      ! Error handling
      !---------------
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc

      ! LOCAL VARIABLES
      !----------------
      INTEGER               :: c

      ! Put any local variables here

      !=================================================================
      ! CCPR_BVOC_Init begins here!
      !=================================================================
      RC = CC_SUCCESS
      ThisLoc = ' -> at CCPR_BVOC_INIT (in process/bvoc/ccpr_bvoc_mod.F90)'

      ! First check if process is activated in config | if not don't allocate arrays or pointers
      if (Config%Bvoc_activate) then

         ! Activate Process
         !------------------
         BvocState%Activate = .true.

         ! Set Scheme Options
         !-------------------
         if (Config%bvoc_scheme < 0) then ! not listed in config
            BvocState%SchemeOpt = 1
         else
            BvocState%SchemeOpt = Config%bvoc_scheme
         endif

         ! CO2 inhibition option
         !TODO: what if it is not given in the configuration file properly
         !------------------
         BvocState%CO2Inhib = Config%megan_co2_inhib

         ! Set CO2 concentration (ppm)
         ! In case it is given as a negative value in config
         !----------------------------
         if (Config%megan_co2_conc_ppm < 0) then
            BvocState%CO2conc = 390.0_fp
         else
            BvocState%CO2conc = Config%megan_co2_conc_ppm
         endif

         ! Check GLOBCO2 if CO2 inhibition is turned on (LISOPCO2 = .TRUE.)
         ! GLOBCO2 should be between 150-1250 ppmv. Isoprene response to
         ! CO2 outside this range has no empirical basis.
         if ( BvocState%CO2Inhib ) then
            if ( BvocState%CO2conc <  150.0_fp .or. &
               BvocState%CO2conc > 1250.0_fp     ) then
               RC = CC_FAILURE
               ErrMsg = 'Global CO2 outside valid range of 150-1250 ppmv!'
               call CC_Error( errMsg, RC, thisLoc )
               return
            endif
         endif

         !--------------------------------------------
         !Find bvoc caterory index in EmisState for future use
         do c = 1, EmisState%nCats
            if (EmisState%Cats(c)%name == 'BVOC') then
               BvocState%CatIndex = c
               exit
            endif
         end do

         ! Set number of species from EmisState
         !----------------------
         BvocState%nBvocSpecies = EmisState%Cats(BvocState%CatIndex)%nSpecies

         !------------------------------------
         ! Allocate emission species index
         ALLOCATE( BvocState%BvocSpeciesIndex(BvocState%nBvocSpecies), STAT=RC )
         CALL CC_CheckVar('BvocState%BvocSpeciesIndex', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN

         ! Allocate emission speceis names
         ALLOCATE( BvocState%BvocSpeciesName(BvocState%nBvocSpecies), STAT=RC )
         CALL CC_CheckVar('BvocState%BvocSpeciesName', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN

         ! Allocate emission flux
         ALLOCATE( BvocState%EmissionPerSpecies(BvocState%nBvocSpecies), STAT=RC )
         CALL CC_CheckVar('BvocState%EmissionPerSpecies', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN

         ! Allocate normalized factor
         ! There should be a different normalization factor for each compound, but
         ! we calculate only 1 normalization factor for all compounds
         ALLOCATE( BvocState%EmisNormFactor(1) , STAT=RC)
         CALL CC_CheckVar('BvocState%EmisNormFactor', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN

         !TODO: emission factor from 7 speceis are read from files and put in MetSate by now
         !      Others are calculated using 'PFT_16', which is also added in MetState
         !      Some met values (last 15 day T average) may need another function and be saved to restart file.

         !TODO: emission species name and ID should read from a namelist. Give them values for now
         !      They are not really used since EmisState controls the species needed now. Keep it as comments for now.
         !BvocState%BvocSpeciesIndex = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21/)
         !BvocState%BvocSpeciesName(/'ISOP','APIN','BPIN','LIMO','SABI','MYRC','CARE',
         !'OCIM','OMON','ALD2','MOH', 'EOH', 'MBOX','FAXX',
         !'AAXX','ACET','PRPE','C2H4','FARN','BCAR','OSQT' /)

      else

         BvocState%Activate = .false.

      endif

   end subroutine CCPR_BVOC_Init

   !>
   !! \brief Run the BVOC process
   !!
   !! \param [IN] MetState The MetState object
   !! \param [INOUT] DiagState The DiagState object
   !! \param [INOUT] BvocState The BvocState object
   !! \param [INOUT] ChemState The ChemState object
   !! \param [INOUT] EmisState The EmisState object
   !! \param [OUT] RC Return code
   !!!>
   !SUBROUTINE CCPr_BVOC_Run( MetState, EmisState, DiagState, BvocState, ChemState, RC )
   SUBROUTINE CCPr_BVOC_Run( MetState, EmisState, DiagState, BvocState, RC )

      ! USE
      USE CCPr_Scheme_MeganV21_Mod, ONLY: CCPr_Scheme_MeganV21  ! Megan scheme
      USE CCPr_BVOC_Common_Mod, Only : CALC_NORM_FAC

      IMPLICIT NONE
      ! INPUT PARAMETERS
      TYPE(MetStatetype),  INTENT(INOUT) :: MetState       ! MetState Instance

      ! INPUT/OUTPUT PARAMETERS
      TYPE(EmisStateType), INTENT(INOUT) :: EmisState   ! Emission Instance
      TYPE(DiagStatetype), INTENT(INOUT) :: DiagState   ! DiagState Instance
      TYPE(BvocStateType), INTENT(INOUT) :: BvocState   ! Bvoc State Instance
      !TYPE(ChemStatetype), INTENT(INOUT) :: ChemState   ! ChemState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                         ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc, MSG
      REAL(fp),PARAMETER :: D2RAD = PI_180
      integer            :: s !, c, cat_bvoc_idx
      ! e-folding time to be applied to long-term past conditions (in days)
      REAL(fp), PARAMETER  :: TAU_DAYS  = 5.0_fp
      ! e-folding time to be applied to short-term past conditions (in hours)
      REAL(fp), PARAMETER  :: TAU_HOURS = 12.0_fp
      REAL(fp)           :: TS_EMIS  !emission time step
      REAL(fp)           :: DNEWFRAC
      REAL(fp)           :: DOLDFRAC
      REAL(fp)           :: HNEWFRAC
      REAL(fp)           :: HOLDFRAC
      logical, save      :: FIRST = .TRUE.

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_BVOC_Run (in process/bvoc/ccpr_bvoc_mod.F90)'

      ! Run the BVOC Scheme if activated
      !----------------------------------
      if (BvocState%Activate) then

         if (FIRST) then
            ! Calculate normalization factor
            ! Not really used now based on Sam's method in which 0.21 is used
            CALL CALC_NORM_FAC( D2RAD, BvocState%EmisNormFactor(1), RC )
            if (RC /= CC_SUCCESS) then
               MSG = 'call on CALC_NORM_FAC failed!'
               call CC_Error( MSG, RC , thisLoc)
            endif
            FIRST = .FALSE.
         endif

         ! Set to 1 day since we know that LAI is updated every day.
         MetState%D_BTW_M=  1.0_fp

         ! Run the BVOC Scheme (only MEGANv2.1 for now)
         ! Put the scheme function in a loop based on EmisSate%cat%nSpecies and calculate
         ! the flux we need and do not need to map the flux back to EmisState
         !-------------------------
         do s = 1, EmisState%Cats(BvocState%CatIndex)%nSpecies

            if (BvocState%SchemeOpt == 1) then ! MEGANv2.1
               call CCPr_Scheme_MeganV21(                                &
                  EmisState%Cats(BvocState%CatIndex)%Species(s)%name,    &
                  EmisState%Cats(BvocState%CatIndex)%Species(s)%Flux(1), &
                  MetState%LAI,                 &
                  MetState%PFT_16,              &
                  MetState%PMISOLAI,            &
                  MetState%Q_DIR_2,             &
                  MetState%Q_DIFF_2,            &
                  MetState%PARDR_LASTXDAYS,     &
                  MetState%PARDF_LASTXDAYS,     &
                  MetState%TS,                  &
                  MetState%T_LASTXDAYS,         &
                  MetState%T_LAST24H,           &
                  MetState%GWETROOT,            &
                  BvocState%CO2Inhib,           &
                  BvocState%CO2conc,            &
                  MetState%SUNCOS,              &
                  MetState%LAT,                 &
                  MetState%DOY,                 &
                  MetState%LocalHour,           &
                  MetState%D_BTW_M,             &
                  RC)
               if (RC /= CC_SUCCESS) then
                  errMsg = 'Error in CCPr_Scheme_MeganV21'
                  CALL CC_Error( errMsg, RC, thisLoc )
               endif
            else
               errMsg =  'ERROR: Unknown BVOC scheme option'
               RC = CC_FAILURE
               CALL CC_Error( errMsg, RC, thisLoc )
               return

            endif  !end if scheme option

            !put it back to BvocState (may not be necessary)
            BvocState%BvocSpeciesIndex(s)  = s
            BvocState%BvocSpeciesName(s)   = EmisState%Cats(BvocState%CatIndex)%Species(s)%name
            BvocState%EmissionPerSpecies(s) = EmisState%Cats(BvocState%CatIndex)%Species(s)%Flux(1)
            BvocState%TotalEmission = BvocState%TotalEmission + BvocState%EmissionPerSpecies(s)

         end do ! for each species requested in EmisState

         !-----------------------------------------------------------------
         ! Update historical temperature / radiation values
         !-----------------------------------------------------------------
         ! Calculate weights for running means of historic variables
         ! DNEWFRAC and DOLDFRAC are the weights given to the current
         ! and existing value, respectively, when updating running means
         ! over the last X days. HNEWFRAC and HOLDFRAC are the same but
         ! for the 24H means. (following MEGAN in HEMCO)

         !TODO: use the same time step with MET for now; emission may have its own in the future; make sure the unit is seconds
         TS_EMIS  = MetState%TSTEP
         DNEWFRAC = TS_EMIS / ( TAU_DAYS * 24.0_fp * 3600.0_fp )
         DOLDFRAC = 1.0_fp - DNEWFRAC
         HNEWFRAC = TS_EMIS / ( TAU_HOURS * 3600.0_fp )
         HOLDFRAC = 1.0_fp - HNEWFRAC

         ! Updated temperature of last 24 hours
         MetState%T_LAST24H = ( HOLDFRAC * MetState%T_LAST24H ) + ( HNEWFRAC * MetState%TS )
         ! Updated temperature of last NUM_DAYS
         MetState%T_LASTXDAYS = ( DOLDFRAC * MetState%T_LASTXDAYS ) + ( DNEWFRAC * MetState%TS )
         ! Updated direct radiation of last NUM_DAYS
         MetState%PARDR_LASTXDAYS = ( DOLDFRAC * MetState%PARDR_LASTXDAYS ) + ( DNEWFRAC * MetState%Q_DIR_2 )
         ! Updated diffuse radiation of last NUM_DAYS
         MetState%PARDF_LASTXDAYS = ( DOLDFRAC * MetState%PARDF_LASTXDAYS ) + ( DNEWFRAC * MetState%Q_DIFF_2 )
         ! Updated LAI of last 24 hours (named LAI_PREVDAY in HEMCO and given to PMISOLAI; here we use PMISOLAI directly)
         MetState%PMISOLAI = ( HOLDFRAC * MetState%PMISOLAI ) + ( HNEWFRAC * MetState%LAI )

         ! Fill Diagnostic Variables
         !TODO: This is better to be done in Restart files if we have them in the future
         !TODO: There also should be some initial values for these historical variables if it is a 'cold' start. The initial values used in HEMCO is:
         !TODO: T_LAST24H =288.15_fp; T_LASTXDAYS =288.15; PARDR_LASTXDAYS=30.0; PARDF_LASTXDAYS=48.0; PMISOLAI =current LAI
         !--------------------------
         DiagState%T_LAST24H        = MetState%T_LAST24H
         DiagState%T_LASTXDAYS      = MetState%T_LASTXDAYS
         DiagState%PARDR_LASTXDAYS  = MetState%PARDR_LASTXDAYS
         DiagState%PARDF_LASTXDAYS  = MetState%PARDF_LASTXDAYS
         DiagState%PMISOLAI         = MetState%PMISOLAI

      endif !if BOVC is activated or not

   end subroutine CCPr_BVOC_Run

   !>
   !! \brief Finalize BVOC
   !!
   !! \param [INOUT] BvocState
   !! \param [OUT] RC Return code
   !!!>
   SUBROUTINE CCPr_BVOC_Final( BVOCState, RC )

      ! USE
      !----

      IMPLICIT NONE

      ! INPUT/OUTPUT PARAMETERS
      TYPE(BvocStateType), INTENT(INOUT) :: BvocState  ! BvocState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                                  ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_BVOC_Final (in process/bvoc/ccpr_BVOC_mod.F90)'

      ! Deallocate any arrays here
      IF ( ASSOCIATED( BvocState%SpcIDs ) ) THEN
         DEALLOCATE( BvocState%SpcIDs, STAT=RC )
         CALL CC_CheckVar('BvocState%SpcIDs', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      ENDIF

      IF ( ASSOCIATED( BvocState%BvocSpeciesIndex ) ) THEN
         DEALLOCATE( BvocState%BvocSpeciesIndex, STAT=RC )
         CALL CC_CheckVar('BvocState%BvocSpeciesIndex', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      ENDIF

      IF ( ASSOCIATED( BvocState%BvocSpeciesName ) ) THEN
         DEALLOCATE( BvocState%BvocSpeciesName, STAT=RC )
         CALL CC_CheckVar('BvocState%BvocSpeciesName', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      ENDIF

      IF ( ASSOCIATED( BvocState%EmissionPerSpecies ) ) THEN
         DEALLOCATE( BvocState%EmissionPerSpecies, STAT=RC )
         CALL CC_CheckVar('BvocState%EmissionPerSpecies', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      ENDIF

      IF ( ASSOCIATED( BvocState%EmisNormFactor ) ) THEN
         DEALLOCATE( BvocState%EmisNormFactor, STAT=RC )
         CALL CC_CheckVar('BvocState%EmisNormFactor', 0, RC)
         IF (RC /= CC_SUCCESS) RETURN
      ENDIF

   end subroutine CCPr_BVOC_Final

END MODULE CCPR_BVOC_Mod
