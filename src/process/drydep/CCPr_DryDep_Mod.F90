!> \brief CCPR drydep state types
!!
!! \defgroup catchem_drydep_process
!!
!! \author Lacey Holland
!! \date 07/2024
!!!>
MODULE CCPR_DryDep_mod
   USE Precision_mod
   USE Error_Mod
   USE DiagState_Mod, Only : DiagStateType
   USE MetState_Mod,  Only : MetStateType
   USE ChemState_Mod, Only : ChemStateType
   USE Config_Opt_Mod, Only : ConfigType
   USE CCPr_drydep_Common_Mod   !to initialize INIT_WEIGHTSS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CCPR_DryDep_Init
   PUBLIC :: CCPR_DryDep_Run
   PUBLIC :: CCPR_DryDep_Finalize
   PUBLIC :: DryDepStateType


   !> \brief DryDepStateType
   !!
   !! DryDepStateType is the process-specific derived type.
   !!
   !! \param Activate Activate Process (True/False)
   !! \param Scheme Scheme Option
   !! \param DryDepSpeciesIndex Effected Chemical Species from DryDep
   !! \param nSpc # of species
   !! \param SpcIDs CATChem species IDs
   !! \param ScaleFactor Scale Factor
   !! \param Resuspension Activate resuspension  (True/False)
   !!
   !! \ingroup core_modules
   !!!>
   TYPE :: DryDepStateType

      ! Process Specific Parameters

      ! Namelist parameters for specific DryDep goes here as well
      !=================================================================
      ! Module specific variables/arrays/data pointers come below
      !=================================================================

      LOGICAL                         :: Activate              ! Activate Process (True/False)
      LOGICAL                         :: Resuspension          ! Activate resuspension  (True/False)
      INTEGER                         :: SchemeOpt             ! Scheme Option (if there is only one SchemeOpt always = 1)
      real                            :: particleradius        ! Particle radius (m)
      real                            :: particledensity       ! Particle density (kg/m^3)
      LOGICAL                         :: co2_effect            ! CO2 effect on drydep
      real(fp)                        :: co2_level             ! CO2 level (ppm)
      real(fp)                        :: co2_reference         ! Reference CO2 level (ppm)
      real, allocatable               :: drydep_frequency(:)   ! could have one per chem species, revisit later
      real, allocatable               :: drydep_vel(:)         ! could have one per chem species, revisit later


   END TYPE DryDepStateType


CONTAINS

   !>
   !! \brief Initialize the CATChem DryDep module
   !!
   !! \param Config       CATCHem configuration options
   !! \param DryDepState   CATCHem PROCESS state
   !! \param ChemState         CATCHem chemical state
   !! \param RC               Error return code
   !!
   !! \ingroup catchem_drydep_process
   !!
   !!!>
   SUBROUTINE CCPR_DryDep_Init( Config, DryDepState, ChemState, RC )
      ! USE

      IMPLICIT NONE
      ! INPUT PARAMETERS
      !-----------------
      TYPE(ConfigType)    :: Config    ! Module options
      TYPE(ChemStateType) :: ChemState ! Chemical state

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      TYPE(DryDepStateType)          :: DryDepState ! DryDep state
      INTEGER,         INTENT(INOUT) :: RC       ! Success or failure

      ! Error handling
      !---------------
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc

      ! LOCAL VARIABLES
      !----------------


      ! Put any local variables here

      !=================================================================
      ! CCPR_DryDep_Init begins here!
      !=================================================================
      ErrMsg = ''
      ThisLoc = ' -> at CCPR_DryDep_INIT (in process/drydep/ccpr_drydep_mod.F90)'

      ! First check if process is activated in config | if not don't allocate arrays or pointers
      if (Config%drydep_activate) then

         ! Activate Process
         !------------------
         DryDepState%Activate = .true.

         ! Set scheme option
         !------------------
         ! For now, the only option is SchemeOpt = 1,2; default is 1
         DryDepState%SchemeOpt = Config%drydep_scheme

         if (DryDepState%SchemeOpt == 1) then !only aerosol drydep

            allocate(DryDepState%drydep_frequency(ChemState%nSpeciesAeroDryDep), STAT=RC)
            IF ( RC /= CC_SUCCESS ) THEN
               ErrMsg = 'Could not Allocate DryDepState%drydep_frequency(ChemState%nSpeciesAeroDryDep)'
               CALL CC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
            DryDepState%drydep_frequency(1:ChemState%nSpeciesAeroDryDep)=ZERO

            allocate(DryDepState%drydep_vel(ChemState%nSpeciesAeroDryDep), STAT=RC)
            IF ( RC /= CC_SUCCESS ) THEN
               ErrMsg = 'Could not Allocate DryDepState%drydep_vel(ChemState%nSpeciesAeroDryDep)'
               CALL CC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
            DryDepState%drydep_vel(1:ChemState%nSpeciesAeroDryDep)=ZERO

         else if (DryDepState%SchemeOpt == 2) then  !aerosol and gas drydep wesely

            allocate(DryDepState%drydep_frequency(ChemState%nSpeciesDryDep), STAT=RC)
            IF ( RC /= CC_SUCCESS ) THEN
               ErrMsg = 'Could not Allocate DryDepState%drydep_frequency(ChemState%nSpeciesDryDep)'
               CALL CC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
            DryDepState%drydep_frequency(1:ChemState%nSpeciesDryDep)=ZERO

            allocate(DryDepState%drydep_vel(ChemState%nSpeciesDryDep), STAT=RC)
            IF ( RC /= CC_SUCCESS ) THEN
               ErrMsg = 'Could not Allocate DryDepState%drydep_vel(ChemState%nSpeciesDryDep)'
               CALL CC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
            DryDepState%drydep_vel(1:ChemState%nSpeciesDryDep)=ZERO

            !calculate the volume distribution of sea salt aerosols (only need to do this once)
            !TODO: The bin is hard coded in ccpr_drydep_common_mod.F90
            CALL INIT_WEIGHTSS(SALA_REDGE_um(1), SALC_REDGE_um(2), RC)
            IF ( RC /= CC_SUCCESS ) THEN
               ErrMsg = 'Could not Allocate arrays in INIT_WEIGHTSS'
               CALL CC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
         end if  ! if (DryDepState%SchemeOpt == 1)

         ! Set other scheme-related  options
         !-----------------------------------
         DryDepState%Resuspension = Config%drydep_resuspension
         DryDepState%co2_effect = Config%drydep_co2_effect
         DryDepState%co2_level = Config%drydep_co2_level
         DryDepState%co2_reference = Config%drydep_co2_reference
      else
         DryDepState%Activate = .false.
      end if

   end subroutine CCPR_DryDep_Init

   !>
   !! \brief Run the DryDep
   !!
   !! \param [IN] MetState - The MetState object
   !! \param [INOUT] DiagState - The DiagState object
   !! \param [INOUT] DryDepState - The DryDepState object
   !! \param [INOUT] ChemState - The ChemState object
   !! \param [OUT] RC Return code
   !!
   !! \ingroup catchem_drydep_process
   !!!>
   SUBROUTINE CCPr_DryDep_Run( MetState, DiagState, DryDepState, ChemState, RC )

      ! USE
      USE constants, only : Cp, g0, VON_KARMAN
      use CCPr_Scheme_GOCART_DryDep_Mod, only : CCPr_Scheme_GOCART_DryDep
      use CCPr_Scheme_Wesely_Mod, only : CCPr_Scheme_Wesely

      IMPLICIT NONE
      ! INPUT PARAMETERS
      TYPE(MetStateType),  INTENT(IN) :: MetState       !< MetState Instance

      ! INPUT/OUTPUT PARAMETERS
      TYPE(DiagStateType), INTENT(INOUT)      :: DiagState       !< DiagState Instance
      TYPE(DryDepStateType), INTENT(INOUT)    :: DryDepState     !< DryDepState Instance
      TYPE(ChemStateType),  INTENT(INOUT)     :: ChemState       !< ChemState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                                 ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc
      INTEGER :: km
      INTEGER :: i !< counter
      real :: radius
      real :: rhop
      real(fp) :: W10, F0   !calculated 10m wind speed from U10M and V10M
      real(fp) :: THIK   !codespell:ignore
      real(fp) :: VD, DDFreq
      real :: drydepf(1,1)
      REAL(fp) :: dqa                                    ! Change in Species due to drydep
      REAL(fp) :: SpecConc                               ! Temporary Species concentration

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_DryDep_Run (in process/drydep/ccpr_DryDep_mod.F90)'

      km = MetState%NLEVS

      ! Run the DryDep Scheme
      !-------------------------
      if (DryDepState%Activate) then
         ! Run the DryDep Scheme
         !-------------------------
         if (DryDepState%SchemeOpt == 1) then
            ! Run the DryDep Scheme - Only Applicable to AEROSOL species
            !-------------------------
            if (ChemState%nSpeciesAeroDryDep > 0) then

               ! loop through aerosol species
               do i = 1, ChemState%nSpeciesAeroDryDep

                  radius = ChemState%chemSpecies(ChemState%AeroDryDepIndex(i))%radius
                  rhop = ChemState%chemSpecies(ChemState%AeroDryDepIndex(i))%density

                  call CCPr_Scheme_GOCART_DryDep( MetState%NLEVS,   &
                     MetState%T,       &
                     MetState%AIRDEN,  &
                     MetState%ZMID,    &
                     MetState%LWI,     &
                     MetState%USTAR,   &
                     MetSTate%PBLH,    &
                     MetState%HFLUX,   &
                     VON_KARMAN,       &
                     Cp,               &
                     g0,               &
                     MetState%Z0H,     &
                     drydepf,          &
                     DryDepState%Resuspension, &
                     radius,           &
                     rhop,             &
                     MetState%U10M,    &
                     MetSTate%V10M,    &
                     MetState%FRLAKE,  &
                     MetState%GWETTOP, &
                     RC)


                  if (RC /= 0) then
                     errMsg = 'Error in GOCART DryDeposition'
                     CALL CC_Error( errMsg, RC, thisLoc )
                  endif  !if (RC /= CC_SUCCESS)

                  ! Fill Diagnostic Variables
                  !--------------------------
                  DryDepState%drydep_frequency(i) = drydepf(1,1)
                  DryDepState%drydep_vel(i) = MetState%ZMID(1) * drydepf(1,1)
                  DiagState%drydep_frequency(i)= drydepf(1,1)
                  DiagState%drydep_vel(i) = MetState%ZMID(1) * drydepf(1,1)

                  ! apply drydep velocities/freq to chem species
                  dqa = 0.
                  SpecConc = ChemState%chemSpecies(ChemState%AeroDryDepIndex(i))%conc(1)
                  dqa = MAX(0.0_fp, SpecConc * (1.-exp(-1*drydepf(1,1) * MetState%TSTEP)))
                  ChemState%chemSpecies(ChemState%AeroDryDepIndex(i))%conc(1) = SpecConc - dqa

               end do ! do i = 1, ChemState%nSpeciesAeroDryDep

            endif  ! if (ChemState%nSpeciesAeroDryDep > 0)

         else if (DryDepState%SchemeOpt == 2) then
            ! Run the DryDep Scheme - Wesely scheme
            !-------------------------
            if (ChemState%nSpeciesDryDep > 0) then

               W10 = sqrt(MetState%U10M**2 + MetState%V10M**2)
               ! loop through aerosol species
               !TODO: nSpeciesAeroDryDep is actually all the drydep species, not just aerosols
               do i = 1, ChemState%nSpeciesDryDep

                  radius = ChemState%chemSpecies(ChemState%DryDepIndex(i))%radius
                  !if (radius > 0.1_fp) radius = radius * 1.e-6_fp !TODO: This is temporary solution if A_RADI is in um in the input
                  rhop = ChemState%chemSpecies(ChemState%DryDepIndex(i))%density
                  !These two can be changed in the function so as not to modify the original values in the States
                  THIK = MetState%BXHEIGHT(1) !codespell:ignore
                  F0 = ChemState%chemSpecies(ChemState%DryDepIndex(i))%dd_f0

                  call CCPr_Scheme_Wesely( &
                     MetState%SWGDN, &
                     MetState%TS,       &
                     MetState%SUNCOSmid,  &
                     F0, &
                     ChemState%chemSpecies(ChemState%DryDepIndex(i))%dd_hstar, &
                     ChemState%chemSpecies(ChemState%DryDepIndex(i))%mw_g/1000.0_fp, &
                     radius,           &
                     rhop,             &
                     MetState%USTAR,   &
                     MetState%OBK,     & !TODO: Need to add Obukhov length to met state
                     MetState%CLDFRC,   &
                     MetState%PBLH,    &
                     THIK,  & !codespell:ignore
                     MetState%Z0,     &
                     MetState%RH(1)/100.0_fp,     & !TODO: input is percent & RH is a array
                     MetState%PS * 100.0_fp,     & !TODO: input is hPa
                     W10,     &
                     ChemState%chemSpecies(ChemState%DryDepIndex(i))%short_name,     &
                     MetState%FRLAI,     & !TODO: whether LAI is separated to each land type?
                     MetState%ILAND,     & !TODO: Need to add land use type to met state
                     MetState%FRLANDUSE,     &
                     MetState%SALINITY,     & !TODO: Need to add salinity to met state
                     MetState%TSKIN,     &
                     MetState%IODIDE,   & !TODO: Need to read from ChemState in the future
                     MetState%LON,     & !TODO: Need to add longitude to met state
                     MetState%LAT,     &
                     DryDepState%co2_effect,     &
                     DryDepState%co2_level,     &
                     DryDepState%co2_reference,     &
                     MetState%LNLPBL,     & !TODO: Need to add PBL option somewhere
                     ChemState%chemSpecies(ChemState%DryDepIndex(i))%is_gas,     &
                     ChemState%chemSpecies(ChemState%DryDepIndex(i))%is_dust,     &
                     ChemState%chemSpecies(ChemState%DryDepIndex(i))%is_seasalt,     &
                     MetState%IsSnow, MetState%IsIce, MetState%IsLand, &
                     ChemState%chemSpecies(ChemState%DryDepIndex(i))%dd_DvzAerSnow,     &
                     ChemState%chemSpecies(ChemState%DryDepIndex(i))%dd_DvzMinVal_snow,     &
                     ChemState%chemSpecies(ChemState%DryDepIndex(i))%dd_DvzMinVal_land,     &
                     VD, DDFreq, RC )

                  if (RC /= CC_SUCCESS ) then
                     errMsg = 'Error in Wesely DryDeposition'
                     CALL CC_Error( errMsg, RC, thisLoc )
                     RETURN
                  endif

                  ! Fill Diagnostic Variables
                  !--------------------------
                  DryDepState%drydep_frequency(i) = DDFreq
                  DryDepState%drydep_vel(i) = VD
                  DiagState%drydep_frequency(i)= DDFreq
                  DiagState%drydep_vel(i) = VD

                  ! apply drydep velocities/freq to chem species (TODO: need to see if the mapping is right)
                  dqa = 0.
                  SpecConc = ChemState%chemSpecies(ChemState%DryDepIndex(i))%conc(1)
                  dqa = MAX(0.0_fp, SpecConc * (1.-exp(-1*drydepf(1,1) * MetState%TSTEP)))
                  ChemState%chemSpecies(ChemState%DryDepIndex(i))%conc(1) = SpecConc - dqa

               end do ! do i = 1, ChemState%nSpeciesDryDep

            endif  ! if (ChemState%nSpeciesDryDep > 0)

         endif  ! if (DryDepState%SchemeOpt == 1 or 2)

         ! TO DO:  apply dry dep velocities/freq to chem species
         write(*,*) 'TODO: Need to figure out how to add back to the chemical species state '

      endif   !  if (DryDepState%Activate)


   end subroutine CCPr_DryDep_Run

   !>
   !! \brief Finalize the DryDep
   !!
   !! \param [INOUT] DryDepState
   !! \param [OUT] RC Return code
   !!!>
   SUBROUTINE CCPr_DryDep_Finalize( DryDepState, RC )

      ! USE
      !----

      IMPLICIT NONE

      ! INPUT/OUTPUT PARAMETERS
      TYPE(DryDepStateType), INTENT(INOUT) :: DryDepState  ! DryDepState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                                  ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_DryDep_Finalize (in process/drydep/ccpr_DryDep_mod.F90)'

      DEALLOCATE( DryDepState%drydep_frequency, STAT=RC )
      IF ( RC /= CC_SUCCESS ) THEN
         ErrMsg = 'Could not Deallocate DryDepState%drydep_frequency'
         CALL CC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      DEALLOCATE( DryDepState%drydep_vel, STAT=RC )
      IF ( RC /= CC_SUCCESS ) THEN
         ErrMsg = 'Could not Deallocate DryDepState%drydep_vel'
         CALL CC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      IF ( ALLOCATED( DMID     ) ) DEALLOCATE( DMID,     STAT=RC )
      IF ( RC /= CC_SUCCESS ) THEN
         ErrMsg = 'Could not Deallocate DMID'
         CALL CC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      IF ( ALLOCATED( SALT_V    ) ) DEALLOCATE( SALT_V,    STAT=RC )
      IF ( RC /= CC_SUCCESS ) THEN
         ErrMsg = 'Could not Deallocate SALT_V'
         CALL CC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF


   end subroutine CCPr_DryDep_Finalize


END MODULE CCPR_DryDep_Mod
