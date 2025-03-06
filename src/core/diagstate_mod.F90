!> \file diagstate_mod.F90
!! \brief Contains the DiagStateType and Diag_Allocate subroutine
!!
!! \ingroup core_modules
!!
!! \details This module contains subroutines and functions related to the DiagState instance of CATChem.
!! It includes subroutines for initializing of the DiagState.
!!!>
module DiagState_Mod
   ! Uses
   USE Precision_Mod
   USE Error_Mod
   USE ChemState_Mod,  only : ChemStateType

   IMPLICIT NONE
   private

   ! PUBLIC :: Zero_DiagState
   PUBLIC :: Diag_Allocate

   !> \brief Data type for storing diagnostic state variables
   !!
   !! \ingroup core_modules
   !!
   !! \param dust_total_flux      Total flux of dust particles [kg m-2 s-1]
   !! \param sea_salt_total_flux  Total flux of sea salt particles [kg m-2 s-1]
   !! \param AOD550               Total AOD at 550nm [1]
   !! \param AOD380               Total AOD at 380nm [1]
   !! \param TOMSAI               TOMS Aerosol Index [1]
   !! \param briggs_plumerise_height Effective plume rise height from Briggs algorithm [m]
   !! \param sofiev_plumerise_height Effective plume rise height from Sofiev algorithm [m]
   !!
   !!!>
   type, public :: DiagStateType

      ! Surface or single-level variables
      REAL(fp) :: dust_total_flux      !< Total flux of dust particles [kg m-2 s-1]
      REAL(fp) :: sea_salt_total_flux  !< Total flux of sea salt particles [kg m-2 s-1]

      ! Aerosol properties
      ! TODO: Add support for multiple aerosol types / wavelengths and more aerosol optical properties
      real(fp), allocatable :: AOD550(:)  !< Total AOD at 550nm [1]
      real(fp), allocatable :: AOD380(:)  !< Total AOD at 380nm [1]
      real(fp), allocatable :: TOMSAI(:)  !< TOMS Aerosol Index [1]


      real(fp) :: briggs_plumerise_height !< Effective plume rise height from Briggs algorithm [m]
      real(fp) :: sofiev_plumerise_height !< Effective plume rise height from Sofiev algorithm [m]


      real(fp), allocatable :: drydep_frequency(:)
      real(fp), allocatable :: drydep_vel(:)

      ! MEGAN historical variables
      ! TODO: this is better to be done in restart files if possible in the future
      real(fp), allocatable :: T_LAST24H         !< temperature of last 24 hours
      real(fp), allocatable :: T_LASTXDAYS       !< temperature of last NUM_DAYS
      real(fp), allocatable :: PARDR_LASTXDAYS   !< direct radiation of last NUM_DAYS
      real(fp), allocatable :: PARDF_LASTXDAYS   !< diffuse radiation of last NUM_DAYS
      real(fp), allocatable :: PMISOLAI          !< LAI of last 24 hours


      ! Species Specific Variables


   end type DiagStateType

CONTAINS

   !> \brief Allocate memory for the diagnostic state variables
   !!
   !! This subroutine allocates memory for the diagnostic state variables.
   !!
   !! \param Config The configuration options
   !! \param GridState The grid state containing information about the grid
   !! \param DiagState The diagnostic state to be allocated
   !! \param RC The return code
   !! \ingroup core_modules
   !!!>
   subroutine Diag_Allocate(Config, DiagState, ChemState, RC)
      ! USES
      ! USE GridState_Mod, ONLY : GridStateType
      USE Config_Opt_Mod, ONLY : ConfigType
      use ChemState_Mod, only : ChemStateType

      ! Arguments
      type(ConfigType),    INTENT(IN)    :: Config
      ! type(GridStateType), INTENT(IN)    :: GridState ! Grid State object
      type(DiagStateType), INTENT(INOUT) :: DiagState ! Diag State object
      type(ChemStateType), INTENT(INOUT) :: ChemState ! Chem State object
      ! OUTPUT Params
      INTEGER,             INTENT(OUT)   :: RC        ! Success or failure

      ! Error handling
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: thisLoc

      ! Initialize
      RC = 0
      ErrMsg = ''
      thisLoc = ' -> at Diag_Allocate (in core/diagstate_mod.F90)'

      ! Nullify all fields for safety's sake before allocating them
      ! This can prevent compilation errors caused by uninitialized values
      ! DiagState%x => NULL()

      ! If dust process is activated then allocate dust related diagnostics
      if (Config%dust_activate) then
         DiagState%dust_total_flux = ZERO
      endif

      ! If sea salt process is activated then allocate sea salt related diagnostics
      if (Config%seasalt_activate) then
         DiagState%sea_salt_total_flux = ZERO
      endif

      ! If dry deposition process is activated then allocate dry dep related diagnostics
      !write (*,*) "ChemState%nSpeciesAeroDryDep=", ChemState%nSpeciesAeroDryDep
      if (Config%drydep_activate) then
         Allocate(DiagState%drydep_frequency(ChemState%nSpeciesAeroDryDep), STAT=RC)
         IF ( RC /= CC_SUCCESS ) THEN
            ErrMsg = 'Could not Allocate DiagState%drydep_frequency(ChemState%nSpeciesAeroDryDep)'
            CALL CC_Error( ErrMsg, RC, thisLoc )
         ENDIF
         DiagState%drydep_frequency(ChemState%nSpeciesAeroDryDep)= ZERO

         Allocate(DiagState%drydep_vel(ChemState%nSpeciesAeroDryDep), STAT=RC)
         IF ( RC /= CC_SUCCESS ) THEN
            ErrMsg = 'Could not Allocate DiagState%drydep_vel(ChemState%nSpeciesAeroDryDep)'
            CALL CC_Error( ErrMsg, RC, thisLoc )
         ENDIF
         DiagState%drydep_vel(ChemState%nSpeciesAeroDryDep)= ZERO

      endif

      ! If bvoc process is activated then allocate bvoc related diagnostics
      if (Config%bvoc_activate) then
         DiagState%T_LAST24H        = ZERO
         DiagState%T_LASTXDAYS      = ZERO
         DiagState%PARDR_LASTXDAYS  = ZERO
         DiagState%PARDF_LASTXDAYS  = ZERO
         DiagState%PMISOLAI         = ZERO
      endif


   end subroutine Diag_Allocate

end module DiagState_Mod
