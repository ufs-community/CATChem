!> \file CCPr_WetDep_Mod.F90
!! \brief Driver for the CCPR large scale wet deposition process
!!
!! \defgroup catchem_wetdep_process
!! The CATChem WetDep_process group holds all the modules for wet deposition.
!!
!! \authors Wei Li
!! \date 04/2025
!!!>
MODULE CCPR_WetDep_mod
   USE Precision_mod
   USE Error_Mod
   USE DiagState_Mod, Only : DiagStateType
   USE MetState_Mod,  Only : MetStateType
   USE ChemState_Mod, Only : ChemStateType, FindSpecByName
   USE Config_Opt_Mod, Only : ConfigType
   USE CCPr_wetdep_Common_Mod

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CCPR_WetDep_Init
   PUBLIC :: CCPR_WetDep_Run
   PUBLIC :: CCPR_WetDep_Finalize
   PUBLIC :: WetDepStateType


CONTAINS

   !>
   !! \brief Initialize the CATChem WetDep module
   !!
   !! \param Config        CATCHem configuration options
   !! \param WetDepState   CATCHem PROCESS state
   !! \param ChemState     CATCHem chemical state
   !! \param RC            Error return code
   !!
   !! \ingroup catchem_wetdep_process
   !!
   !!!>
   SUBROUTINE CCPR_WetDep_Init( Config, WetDepState, MetState, ChemState, RC )
      ! USE

      IMPLICIT NONE
      ! INPUT PARAMETERS
      !-----------------
      TYPE(ConfigType)    :: Config    ! Module options
      TYPE(ChemStateType) :: ChemState ! Chemical state
      TYPE(MetStateType)  :: MetState  ! Meteorological state

      ! INPUT/OUTPUT PARAMETERS
      !------------------------
      TYPE(WetDepStateType)          :: WetDepState ! WetDep state
      INTEGER,         INTENT(INOUT) :: RC       ! Success or failure

      ! Error handling
      !---------------
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc

      ! LOCAL VARIABLES
      !----------------

      !=================================================================
      ! CCPR_WetDep_Init begins here!
      !=================================================================
      ErrMsg = ''
      ThisLoc = ' -> at CCPR_WetDep_INIT (in process/wetdep/ccpr_wetdep_mod.F90)'

      ! First check if process is activated in config | if not don't allocate arrays or pointers
      if (Config%wetdep_activate) then

         ! Activate Process
         !------------------
         WetDepState%Activate = .true.

         ! Set scheme option
         !------------------
         ! For now, the only option is 1
         WetDepState%SchemeOpt = Config%wetdep_scheme

         allocate(WetDepState%WetDep_Flux(ChemState%nSpeciesWetDep, MetState%NLEVS), STAT=RC)
         IF ( RC /= CC_SUCCESS ) THEN
            ErrMsg = 'Could not Allocate WetDepState%WetDep_Flux(ChemState%nSpeciesWetDep)'
            CALL CC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
         WetDepState%WetDep_Flux(1:ChemState%nSpeciesWetDep, 1:MetState%NLEVS)=ZERO

         ! Set other scheme-related  options if having any
         !------------------------------------------------

      else
         WetDepState%Activate = .false.
      end if

   end subroutine CCPR_WetDep_Init

   !>
   !! \brief Run the WetDep
   !!
   !! \param [IN] MetState - The MetState object
   !! \param [INOUT] WetDepState - The WetDepState object
   !! \param [INOUT] ChemState - The ChemState object
   !! \param [OUT] RC Return code
   !!
   !! \ingroup catchem_wetdep_process
   !!!>
   SUBROUTINE CCPr_WetDep_Run( MetState, WetDepState, ChemState, RC )

      ! USE
      USE constants, only : g0
      use CCPr_Scheme_Jacob_WetDep_Mod, only : CCPr_Scheme_Jacob_WetDep


      IMPLICIT NONE
      ! INPUT PARAMETERS
      TYPE(MetStateType),  INTENT(IN) :: MetState       !< MetState Instance

      ! INPUT/OUTPUT PARAMETERS
      !TYPE(DiagStateType), INTENT(INOUT)      :: DiagState       !< DiagState Instance
      TYPE(WetDepStateType), INTENT(INOUT)    :: WetDepState     !< WetDepState Instance
      TYPE(ChemStateType),  INTENT(INOUT)     :: ChemState       !< ChemState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                                 ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc
      INTEGER  :: i !< counter
      INTEGER  :: H2O2_id = 1, SO4_id = 1 !update these ids only for SO2
      real(fp), dimension(:), allocatable :: fluxout

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_WetDep_Run (in process/wetdep/ccpr_WetDep_mod.F90)'

      ! Run the WetDep Scheme
      !-------------------------
      if (WetDepState%Activate) then
         ! Run the WetDep Scheme
         !-------------------------
         if (ChemState%nSpeciesWetDep > 0) then

            ! loop through all wetdep species
            do i = 1, ChemState%nSpeciesWetDep

               if (WetDepState%SchemeOpt == 1) then
                  ! Run the WetDep Scheme - Jacob scheme
                  !-------------------------
                  !initialize the fluxout array
                  allocate(fluxout(1:MetState%NLEVS)); fluxout = ZERO
                  !get H2O2 and SO4 id
                  if (ChemState%chemSpecies(ChemState%WetDepIndex(i))%short_name == 'SO2') then
                     call FindSpecByName(ChemState, 'H2O2', H2O2_id, RC)
                     call FindSpecByName(ChemState, 'SO4', SO4_id, RC)

                     if (RC /= CC_SUCCESS) then
                        errMsg = 'Error in finding H2O2 or SO4 id'
                        CALL CC_Error( errMsg, RC, thisLoc )
                        RETURN
                     endif
                  end if

                  !run the wetdep Jacob scheme
                  call CCPr_Scheme_Jacob_WetDep( &
                     MetState%NLEVS, &
                     MetState%TSTEP,       &
                     ChemState%chemSpecies(ChemState%WetDepIndex(i))%short_name,   &
                     ChemState%chemSpecies(ChemState%WetDepIndex(i))%is_aerosol,   &
                     ChemState%chemSpecies(ChemState%WetDepIndex(i))%wd_LiqAndGas,   &
                     ChemState%chemSpecies(ChemState%WetDepIndex(i))%henry_k0,     &
                     ChemState%chemSpecies(ChemState%WetDepIndex(i))%henry_cr,     &
                     ChemState%chemSpecies(ChemState%WetDepIndex(i))%henry_pKa,    &
                     ChemState%chemSpecies(ChemState%WetDepIndex(i))%wd_retfactor,    &
                     ChemState%chemSpecies(ChemState%WetDepIndex(i))%wd_convfacI2G,   &
                     g0,           &
                     ChemState%chemSpecies(ChemState%WetDepIndex(i))%radius_wet * 1.0e+6_fp,   & !Note radius is in m in the species yaml file
                     ChemState%chemSpecies(ChemState%WetDepIndex(i))%wd_rainouteff,  &
                     1.0_fp,    &  ! washout tuning factor; use 1.0 for now
                     1.0_fp,    &  ! radius_thr; use 1.0 um to be consistent with GC
                     MetState%PEDGE_DRY, &  !TODO: is this the right variable for ple?
                     MetState%T, &
                     MetState%MAIRDEN, &
                     MetState%PFLLSAN, &
                     MetState%PFILSAN, &
                     MetState%REEVAPLS, &
                     MetState%AIRDEN, &
                     ChemState%ChemSpecies(ChemState%WetDepIndex(i))%conc, & !make sure unit is kg/kg
                     ChemState%ChemSpecies(H2O2_id)%conc, & !make sure unit is kg/kg
                     ChemState%ChemSpecies(SO4_id)%conc, & !make sure unit is kg/kg
                     fluxout, &
                     RC )

                  if (RC /= CC_SUCCESS ) then
                     errMsg = 'Error in Jacob WetDeposition'
                     CALL CC_Error( errMsg, RC, thisLoc )
                     RETURN
                  endif

                  ! Fill WetDepState Variables
                  !--------------------------
                  WetDepState%WetDep_Flux(i,:) = fluxout
                  !TODO: is this necessary to save to diag?
                  !DiagState%wetdep_flux(i, :) = fluxout
                  deallocate(fluxout)

               end if  !if (WetDepState%SchemeOpt == 1)

            end do ! do i = 1, ChemState%nSpeciesWetDep

         endif  ! if (ChemState%nSpeciesWetDep > 0)

      endif   !  if (WetDepState%Activate)


   end subroutine CCPr_WetDep_Run

   !>
   !! \brief Finalize the WetDep
   !!
   !! \param [INOUT] WetDepState
   !! \param [OUT] RC Return code
   !!!>
   SUBROUTINE CCPr_WetDep_Finalize( WetDepState, RC )

      ! USE
      !----

      IMPLICIT NONE

      ! INPUT/OUTPUT PARAMETERS
      TYPE(WetDepStateType), INTENT(INOUT) :: WetDepState  ! WetDepState Instance

      ! OUTPUT PARAMETERS
      INTEGER, INTENT(OUT) :: RC                                  ! Return Code

      ! LOCAL VARIABLES
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at CCPr_WetDep_Finalize (in process/wetdep/ccpr_WetDep_mod.F90)'


      IF ( ALLOCATED( WetDepState%WetDep_Flux ) ) DEALLOCATE( WetDepState%WetDep_Flux, STAT=RC )
      IF ( RC /= CC_SUCCESS ) THEN
         ErrMsg = 'Could not Deallocate WetDepState%WetDep_Flux'
         CALL CC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

   end subroutine CCPr_WetDep_Finalize


END MODULE CCPR_WetDep_Mod
