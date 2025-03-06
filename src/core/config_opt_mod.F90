!> \file config_opt_mod.F90
!! \brief This file contains the module for catchem configuration options
!!
!! This file contains the module for catchem configuration options
!!
!! \ingroup core_modules
!!!>
MODULE Config_Opt_Mod
   !
   ! !USES:
   !
   USE PRECISION_MOD    ! For CATChem Precision (fp)

   IMPLICIT NONE
   PRIVATE
   !
   ! !PUBLIC MEMBER FUNCTIONS:
   !
   PUBLIC :: Set_Config
   PUBLIC :: Cleanup_Config
   !
   ! !PUBLIC DATA MEMBERS:
   !
   !=========================================================================
   ! Derived type for Input Options
   !=========================================================================
   !> \brief Derived type for Input Options
   !!
   !! ConfigType contains the following variables:
   !! - `numCPUs` : Number of MPI procs
   !! - `thisCPU` : Local MPI process handle
   !! - `MPIComm` : MPI Communicator Handle
   !! - `isMPI` : Is this an MPI sim?
   !! - `amIRoot` : Are we on the root CPU?
   !! - `DryRun` : Is this a dry run?
   !! - `SimulationName` : Name of the simulation
   !! - `SpcDatabaseFile` : Name of the species database file
   !! - `VerboseRequested` : Is the user requesting verbose mode
   !! - `VerboseOnCores` : Which cores should be verbose
   !! - `Verbose` : Is verbose mode on?
   !! - `dust_activate` : Activate dust process
   !! - `dust_scheme_opt` : Scheme option for dust process
   !! - `seasalt_activate` : Activate seasalt process
   !! - `seasalt_scheme_opt` : Scheme option for seasalt process
   !! - `drydep_activate` : Activate drydep process
   !! - `drydep_scheme` : Scheme option for drydep process
   !! - `drydep_resuspension` : Activate resuspension
   !! - `bvoc_activate` : Activate BVOC process
   !! - `bvoc_scheme` : Scheme option for BVOC process
   !! - `megan_co2_inhib` : use CO2 inhibition for isoprene?
   !! - `megan_co2_conc_ppm` : If so, provide CO2 concentrations
   !!
   !! \ingroup core_modules
   !!!>
   TYPE, PUBLIC :: ConfigType

      !----------------------------------------
      ! General Runtime & Distributed Comp Info
      !----------------------------------------
      INTEGER                     :: numCPUs    ! Number of MPI procs
      INTEGER                     :: thisCPU    ! Local MPI process handle
      INTEGER                     :: MPIComm    ! MPI Communicator Handle
      LOGICAL                     :: isMPI      ! Is this an MPI sim?
      LOGICAL                     :: amIRoot    ! Is this the root cpu?

      !----------------------------------------
      ! Dry run info (print out file names)
      !----------------------------------------
      LOGICAL                     :: DryRun     ! Is this a dry run?

      !----------------------------------------
      ! SIMULATION MENU fields
      !----------------------------------------
      CHARACTER(LEN=255)          :: SimulationName
      CHARACTER(LEN=255)          :: Emission_File
      CHARACTER(LEN=255)          :: Species_File
      LOGICAL                     :: VerboseRequested
      CHARACTER(LEN=10)           :: VerboseOnCores
      LOGICAL                     :: Verbose

      !-----------------------------------------
      ! PROCESSING MENU fields
      !-----------------------------------------

      ! Dust Process
      LOGICAL                     :: dust_activate
      INTEGER                     :: dust_scheme
      INTEGER                     :: dust_drag_opt  ! Fengsha Option for drag Parameterization (1 MB95; 2 Input Value)
      INTEGER                     :: dust_moist_opt ! Fengsha Option for moisture Parameterization (1 Fecan; 2 Shao)
      INTEGER                     :: dust_horizflux_opt ! Horizontal Flux Calculation Option
      real(fp)                    :: dust_alpha
      real(fp)                    :: dust_beta

      ! SeaSalt Process
      LOGICAL                     :: seasalt_activate
      LOGICAL                     :: seasalt_weibull
      LOGICAL                     :: seasalt_hoppel
      INTEGER                     :: seasalt_scheme
      real(fp)                    :: seasalt_scalefactor

      ! BVOC Process
      LOGICAL                     :: bvoc_activate
      INTEGER                     :: bvoc_scheme
      LOGICAL                     :: megan_co2_inhib
      real(fp)                    :: megan_co2_conc_ppm

      ! Plumerise Process
      LOGICAL                     :: plumerise_activate

      ! DryDeposition Process
      LOGICAL                     :: drydep_activate
      INTEGER                     :: drydep_scheme
      LOGICAL                     :: drydep_resuspension  !< Turn on resuspension

   END TYPE ConfigType

CONTAINS

   !> \brief Initialize the Config options
   !!
   !! This subroutine initializes the Config options
   !!
   !! \param am_I_Root  Are we on the root CPU?
   !! \param Config     The Config object
   !! \param RC         The return code
   !!
   !! \ingroup core_modules
   !!!>
   SUBROUTINE Set_Config( am_I_Root, Config, RC )
      !
      ! !USES:
      !
      USE Error_Mod
      !
      ! !INPUT PARAMETERS:
      !
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      !
      ! !INPUT/OUTPUT PARAMETERS:
      !
      TYPE(ConfigType), INTENT(INOUT) :: Config   ! Input Options object
      !
      ! !OUTPUT PARAMETERS:
      !
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
      !
      ! !LOCAL VARIABLES:
      !

      !----------------------------------------
      ! Initialize
      ! Set pointers to NULL for safety's sake
      !----------------------------------------
      RC                               =  CC_SUCCESS

      !----------------------------------------
      ! General Runtime & Distributed Comp Info
      !----------------------------------------
      Config%amIRoot                = am_I_Root
      Config%isMPI                  = .FALSE.
      Config%numCPUs                = 1
      Config%thisCPU                = -1
      Config%MPIComm                = -1

      !----------------------------------------
      ! Dry run info (print out file names)
      !----------------------------------------
      Config%DryRun                 = .FALSE.

      !-----------------------------------------
      ! PROCESSING MENU fields
      !-----------------------------------------
      ! Dust Process
      Config%dust_activate = .FALSE.
      Config%dust_scheme = 1
      Config%dust_drag_opt = 1
      Config%dust_moist_opt = 1
      Config%dust_horizflux_opt = 1

      ! SeaSalt Process
      Config%seasalt_activate = .FALSE.
      Config%seasalt_scheme = 1

      ! Dry Dep Process
      Config%drydep_activate = .FALSE.
      Config%drydep_scheme = 1
      Config%drydep_resuspension = .FALSE.

      !BVOC Process
      Config%bvoc_activate = .FALSE.
      Config%bvoc_scheme = 1
      Config%megan_co2_inhib = .TRUE.
      Config%megan_co2_conc_ppm = 390.0_fp

   END SUBROUTINE Set_Config
   !> \brief Cleanup the Config options
   !!
   !! This subroutine cleans up the Config options
   !!
   !! \param Config     The Config object
   !! \param RC         The return code
   !!
   !! \ingroup core_modules
   !!!>
   SUBROUTINE Cleanup_Config( RC )
      !
      ! !USES:
      !
      USE Error_Mod
      !
      ! !INPUT/OUTPUT PARAMETERS:
      !
      ! TYPE(ConfigType), INTENT(INOUT) :: Config   ! Input Options object
      !
      ! !OUTPUT PARAMETERS:
      !
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure

      ! Assume success
      RC = CC_SUCCESS

      !======================================================================
      ! Deallocate fields of the Input Options object
      !======================================================================

      ! Nothing to do yet

   END SUBROUTINE Cleanup_Config

END MODULE Config_Opt_Mod
