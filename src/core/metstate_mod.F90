! \file metstate_mod.F90
!! \brief Module for meteorology state variables
!!
!! This module contains subroutines and functions related to the MetStateType instance of CATChem.
!! It includes subroutines for initializing of the MetStateType.
!!
!! \ingroup core_modules
!!!>
MODULE MetState_Mod
   !
   ! USES:
   !
   USE Error_Mod
   USE Precision_Mod
   USE GridGeometry_Mod
   USE Met_Utilities_Mod
   USE TimeState_Mod, only: TimeStateType



   IMPLICIT NONE
   PRIVATE

   !PUBLIC :: MetStateType           ! Main data type

   !=========================================================================
   ! Derived type for Meteorology State
   !=========================================================================

   ! \brief Derived type for Meteorology State
   !!
   !! Contains all meteorological state variables for CATChem including
   !! land, radiation, flux, cloud, and state-related fields. Use type-bound
   !! procedures for initialization, cleanup, validation, and memory usage.
   !!
   !! \ingroup core_modules
   !!!>
   TYPE, PUBLIC :: MetStateType
      CHARACTER(LEN=3)             :: State     = 'MET'    !< Name of this state
      INTEGER                      :: NLEVS     = 127      !< Number of vertical levels (default)
      TYPE(GridGeometryType) :: geometry
      INTEGER                      :: NSURFTYPE = 20       !< Number of surface types (default)
      ! Grid flags (2D: nx, ny)
      LOGICAL, ALLOCATABLE         :: IsLand(:,:)       !< Is this a land grid box?
      LOGICAL, ALLOCATABLE         :: IsWater(:,:)      !< Is this a water grid box?
      LOGICAL, ALLOCATABLE         :: IsIce(:,:)        !< Is this an ice grid box?
      LOGICAL, ALLOCATABLE         :: IsSnow(:,:)       !< Is this a snow grid box?
      ! Vertical flags and arrays (3D: nx, ny, nz)
      LOGICAL,  ALLOCATABLE        :: InStratMeso(:,:,:)    !< Are we in the stratosphere or mesosphere?
      LOGICAL,  ALLOCATABLE        :: InStratosphere(:,:,:) !< Are we in the stratosphere?
      LOGICAL,  ALLOCATABLE        :: InTroposphere(:,:,:)  !< Are we in the troposphere?
      LOGICAL,  ALLOCATABLE        :: InPbl(:,:,:)          !< Are we in the PBL?
      LOGICAL,  ALLOCATABLE        :: IsLocalNoon(:,:)      !< Is it local noon (between 11 and 13 local solar time)?
      ! Surface properties (2D: nx, ny)
      REAL(fp), ALLOCATABLE        :: AREA_M2(:,:)      !< Grid box surface area [m2]
      INTEGER,  ALLOCATABLE        :: LWI(:,:)          !< Land water ice mask (0-sea, 1-land, 2-ice)
      INTEGER,  ALLOCATABLE        :: DLUSE(:,:)        !< Dominant land-use type
      REAL(fp), ALLOCATABLE        :: FRVEG(:,:)        !< Fraction of veg [1]
      REAL(fp), ALLOCATABLE        :: FRLAKE(:,:)       !< Fraction of lake [1]
      REAL(fp), ALLOCATABLE        :: FRLAND(:,:)       !< Fraction of land [1]
      REAL(fp), ALLOCATABLE        :: FRLANDIC(:,:)     !< Fraction of land ice [1]
      REAL(fp), ALLOCATABLE        :: FROCEAN(:,:)      !< Fraction of ocean [1]
      REAL(fp), ALLOCATABLE        :: FRSEAICE(:,:)     !< Sfc sea ice fraction
      REAL(fp), ALLOCATABLE        :: FRSNO(:,:)        !< Sfc snow fraction
      REAL(fp), ALLOCATABLE        :: LAI(:,:)          !< Leaf area index [m2/m2] (online) Dominant
      REAL(fp), ALLOCATABLE        :: GVF(:,:)          !< Green Vegetative Fraction
      ! Dust Only Variables
      REAL(fp), ALLOCATABLE        :: RDRAG(:,:)        !< Drag Partition [1]
      REAL(fp), ALLOCATABLE        :: USTAR_THRESHOLD(:,:) !< Threshold friction velocity [m/s]
      REAL(fp), ALLOCATABLE        :: SSM(:,:)          !< Sediment Supply Map [1]
      ! Surface and ice properties (2D: nx, ny)
      REAL(fp), ALLOCATABLE        :: SEAICE00(:,:)     !< Sea ice coverage 00-10%
      REAL(fp), ALLOCATABLE        :: SEAICE10(:,:)     !< Sea ice coverage 10-20%
      REAL(fp), ALLOCATABLE        :: SEAICE20(:,:)     !< Sea ice coverage 20-30%
      REAL(fp), ALLOCATABLE        :: SEAICE30(:,:)     !< Sea ice coverage 30-40%
      REAL(fp), ALLOCATABLE        :: SEAICE40(:,:)     !< Sea ice coverage 40-50%
      REAL(fp), ALLOCATABLE        :: SEAICE50(:,:)     !< Sea ice coverage 50-60%
      REAL(fp), ALLOCATABLE        :: SEAICE60(:,:)     !< Sea ice coverage 60-70%
      REAL(fp), ALLOCATABLE        :: SEAICE70(:,:)     !< Sea ice coverage 70-80%
      REAL(fp), ALLOCATABLE        :: SEAICE80(:,:)     !< Sea ice coverage 80-90%
      REAL(fp), ALLOCATABLE        :: SEAICE90(:,:)     !< Sea ice coverage 90-100%
      REAL(fp), ALLOCATABLE        :: SNODP(:,:)        !< Snow depth [m]
      REAL(fp), ALLOCATABLE        :: SNOMAS(:,:)       !< Snow mass [kg/m2]

      ! Soil and land use arrays (2D for counts, 3D for fractions)
      INTEGER,  ALLOCATABLE        :: DSOILTYPE(:,:)    !< Dominant soil type
      REAL(fp), ALLOCATABLE        :: CLAYFRAC(:,:)     !< Fraction of clay [1]
      REAL(fp), ALLOCATABLE        :: SANDFRAC(:,:)     !< Fraction of sand [1]
      INTEGER,  ALLOCATABLE        :: nLNDTYPE(:,:)     !< # of landtypes in box (I,J)
      REAL(fp), ALLOCATABLE        :: GWETTOP(:,:)      !< Top soil moisture [1]
      REAL(fp), ALLOCATABLE        :: GWETROOT(:,:)     !< Root Zone soil moisture [1]
      REAL(fp), ALLOCATABLE        :: WILT(:,:)         !< Wilt point [1]
      INTEGER                      :: nSOIL             !< # number of soil layers
      INTEGER                      :: nSOILTYPE         !< # number of soil types
      REAL(fp), ALLOCATABLE        :: SOILM(:,:,:)      !< Volumetric Soil moisture [m3/m3] (nx,ny,nsoil)
      REAL(fp), ALLOCATABLE        :: SOILT(:,:,:)      !< Temperature of soil layer [K] (nx,ny,nsoil)
      REAL(fp), ALLOCATABLE        :: FRLANDUSE(:,:,:)  !< Fractional Land Use (nx,ny,nlanduse)
      REAL(fp), ALLOCATABLE        :: FRSOIL(:,:,:)     !< Fractional Soil (nx,ny,nsoil)
      REAL(fp), ALLOCATABLE        :: FRLAI(:,:,:)      !< LAI in each Fractional Land use type [m2/m2] (nx,ny,nlanduse)
      INTEGER, ALLOCATABLE         :: ILAND(:,:,:)      !< Land type ID in current grid box (nx,ny,nlanduse)
      ! Location arrays (1D for single point, 2D for grid)
      real(fp), ALLOCATABLE        :: LAT(:,:)         !< Latitude
      real(fp), ALLOCATABLE        :: LON(:,:)         !< Longitude
      character(len=20)            :: LUCNAME          !< name of land use category
      ! Surface meteorological properties (2D: nx, ny)
      REAL(fp), ALLOCATABLE        :: ALBD_VIS(:,:)     !< Visible surface albedo [1]
      REAL(fp), ALLOCATABLE        :: ALBD_NIR(:,:)     !< Near-IR surface albedo [1]
      REAL(fp), ALLOCATABLE        :: ALBD_UV(:,:)      !< UV surface albedo [1]
      REAL(fp), ALLOCATABLE        :: PARDR(:,:)        !< Direct photsynthetically active radiation [W/m2]
      REAL(fp), ALLOCATABLE        :: PARDF(:,:)        !< Diffuse photsynthetically active radiation [W/m2]
      REAL(fp), ALLOCATABLE        :: SUNCOS(:,:)       !< COS(solar zenith angle) at current time
      REAL(fp), ALLOCATABLE        :: SUNCOSmid(:,:)    !< COS(solar zenith angle) at midpoint of chem timestep
      REAL(fp), ALLOCATABLE        :: SUNCOSsum(:,:)    !< Sum of COS(SZA) for HEMCO OH diurnal variability
      REAL(fp), ALLOCATABLE        :: SZAFACT(:,:)      !< Diurnal scale factor for HEMCO OH diurnal variability (computed) [1]
      REAL(fp), ALLOCATABLE        :: SWGDN(:,:)        !< Incident radiation @ ground [W/m2]
      REAL(fp), ALLOCATABLE        :: EFLUX(:,:)        !< Latent heat flux [W/m2]
      REAL(fp), ALLOCATABLE        :: HFLUX(:,:)        !< Sensible heat flux [W/m2]
      REAL(fp), ALLOCATABLE        :: U10M(:,:)         !< E/W wind speed @ 10m ht [m/s]
      REAL(fp), ALLOCATABLE        :: USTAR(:,:)        !< Friction velocity [m/s]
      REAL(fp), ALLOCATABLE        :: V10M(:,:)         !< N/S wind speed @ 10m ht [m/s]
      REAL(fp), ALLOCATABLE        :: Z0(:,:)           !< Surface roughness height [m]
      REAL(fp), ALLOCATABLE        :: Z0H(:,:)          !< Surface roughness height, for heat (thermal roughness) [m]
      REAL(fp), ALLOCATABLE        :: FRZ0(:,:,:)       !< Aerodynamic Roughness Length per FRLANDUSE (nx,ny,nlanduse)
      REAL(fp), ALLOCATABLE        :: PBLH(:,:)         !< PBL height [m]
      REAL(fp), ALLOCATABLE        :: SALINITY(:,:)     !< Salinity of the ocean [part per thousand]
      REAL(fp), ALLOCATABLE        :: CMM(:,:)          !< Aerodynamic conductance [m/s]
      REAL(fp), ALLOCATABLE        :: ORO(:,:)          !< surface height above sea level [m]
      REAL(fp), ALLOCATABLE        :: RCA(:,:)          !< Aerodynamic resistance in canopy [s/m]
      REAL(fp), ALLOCATABLE        :: WCA(:,:)          ! canopy water amount [kg/m2]
      ! 3D volumetric fields (3D: nx, ny, nz)
      REAL(fp), ALLOCATABLE        :: F_OF_PBL(:,:,:)       !< Fraction of box within PBL [1]
      REAL(fp), ALLOCATABLE        :: F_UNDER_PBLTOP(:,:,:) !< Fraction of box under PBL top
      real(fp), ALLOCATABLE        :: OBK(:,:)          !< Monin-Obhukov length [m]
      ! Cloud and precipitation properties (2D for surface, 3D for volumetric)
      REAL(fp), ALLOCATABLE        :: CLDFRC(:,:)       !< Column cloud fraction [1]
      REAL(fp), ALLOCATABLE        :: CONV_DEPTH(:,:)   !< Convective cloud depth [m]
      REAL(fp), ALLOCATABLE        :: FLASH_DENS(:,:)   !< Lightning flash density [#/km2/s]
      REAL(fp), ALLOCATABLE        :: CNV_FRC(:,:)      !< Convective fraction [1]
      REAL(fp), ALLOCATABLE        :: CLDF(:,:,:)       !< 3-D cloud fraction [1]
      REAL(fp), ALLOCATABLE        :: CMFMC(:,:,:)      !< Cloud mass flux [kg/m2/s]
      REAL(fp), ALLOCATABLE        :: DQRCU(:,:,:)      !< Conv precip production rate [kg/kg/s] (assume per dry air)
      REAL(fp), ALLOCATABLE        :: DQRLSAN(:,:,:)    !< LS precip prod rate [kg/kg/s] (assume per dry air)
      REAL(fp), ALLOCATABLE        :: DTRAIN(:,:,:)     !< Detrainment flux [kg/m2/s]
      REAL(fp), ALLOCATABLE        :: PRECANV(:,:)      !< Anvil previp @ ground [kg/m2/s] -> [mm/day]
      REAL(fp), ALLOCATABLE        :: PRECCON(:,:)      !< Conv  precip @ ground [kg/m2/s] -> [mm/day]
      REAL(fp), ALLOCATABLE        :: PRECLSC(:,:)      !< Large-scale precip @ ground kg/m2/s] -> [mm/day]
      real(fp), ALLOCATABLE        :: REEVAPLS(:,:,:)   !< Evap of precip LS+anvil [kg/kg/s] (assume per dry air)
      ! 3D cloud and precipitation arrays
      REAL(fp), ALLOCATABLE        :: QI(:,:,:)         !< Mass fraction of cloud ice water [kg/kg dry air]
      REAL(fp), ALLOCATABLE        :: QL(:,:,:)         !< Mass fraction of cloud liquid water [kg/kg dry air]
      REAL(fp), ALLOCATABLE        :: PFICU(:,:,:)      !< Dwn flux ice prec:conv [kg/m2/s]
      REAL(fp), ALLOCATABLE        :: PFILSAN(:,:,:)    !< Dwn flux ice prec:LS+anv [kg/m2/s] (nx,ny,nz+1)
      REAL(fp), ALLOCATABLE        :: PFLCU(:,:,:)      !< Dwn flux liq prec:conv [kg/m2/s]
      REAL(fp), ALLOCATABLE        :: PFLLSAN(:,:,:)    !< Dwn flux liq prec:LS+anv [kg/m2/s] (nx,ny,nz+1)
      REAL(fp), ALLOCATABLE        :: TAUCLI(:,:,:)     !< Opt depth of ice clouds [1]
      REAL(fp), ALLOCATABLE        :: TAUCLW(:,:,:)     !< Opt depth of H2O clouds [1]
      ! Surface scalars (now 2D: nx, ny)
      REAL(fp), ALLOCATABLE        :: PHIS(:,:)         !< Surface geopotential height [m2/s2]
      REAL(fp), ALLOCATABLE        :: PS_WET(:,:)       !< Wet surface pressure at start of timestep [Pa]
      REAL(fp), ALLOCATABLE        :: PS_DRY(:,:)       !< Dry surface pressure at start of timestep [Pa]
      REAL(fp), ALLOCATABLE        :: QV2M(:,:)         !< Specific Humidity at 2m [kg/kg]
      REAL(fp), ALLOCATABLE        :: T2M(:,:)          !< Temperature 2m [K]
      REAL(fp), ALLOCATABLE        :: TS(:,:)           !< Surface temperature [K]
      REAL(fp), ALLOCATABLE        :: TSKIN(:,:)        !< Surface skin temperature [K]
      REAL(fp), ALLOCATABLE        :: SST(:,:)          !< Sea surface temperature [K]
      REAL(fp), ALLOCATABLE        :: SLP(:,:)          !< Sea level pressure [Pa]
      REAL(fp), ALLOCATABLE        :: PS(:,:)           !< Surface Pressure [Pa]
      REAL(fp), ALLOCATABLE        :: TO3(:,:)          !< Total overhead O3 column [DU]
      REAL(fp), ALLOCATABLE        :: TROPP(:,:)        !< Tropopause pressure [Pa]
      INTEGER,  ALLOCATABLE        :: TropLev(:,:)      !< Tropopause level [1]
      REAL(fp), ALLOCATABLE        :: TropHt(:,:)       !< Tropopause height [km]
      ! 3D atmospheric variables (3D: nx, ny, nz)
      REAL(fp), ALLOCATABLE        :: Z(:,:,:)          !< Geopotential Height @ level edges [m] (nx,ny,nz+1)
      REAL(fp), ALLOCATABLE        :: ZMID(:,:,:)       !< Mid Layer Geopotential Height [m]
      REAL(fp), ALLOCATABLE        :: BXHEIGHT(:,:,:)   !< Grid box height [m] (dry air)
      REAL(fp), ALLOCATABLE        :: QV(:,:,:)         !< Specific Humidity [kg/kg]
      REAL(fp), ALLOCATABLE        :: T(:,:,:)          !< Temperature [K]
      REAL(fp), ALLOCATABLE        :: THETA(:,:,:)      !< Potential temperature [K]
      REAL(fp), ALLOCATABLE        :: TV(:,:,:)         !< Virtual temperature [K]
      REAL(fp), ALLOCATABLE        :: V(:,:,:)          !< N/S component of wind [m s-1]
      REAL(fp), ALLOCATABLE        :: U(:,:,:)          !< E/W component of wind [m s-1]
      REAL(fp), ALLOCATABLE        :: OMEGA(:,:,:)      !< Updraft velocity [Pa/s]
      REAL(fp), ALLOCATABLE        :: RH(:,:,:)         !< Relative humidity [fraction, not %]
      REAL(fp), ALLOCATABLE        :: SPHU(:,:,:)       !< Specific humidity [g H2O/kg tot air]
      REAL(fp), ALLOCATABLE        :: AIRDEN(:,:,:)     !< Wet air density [kg/m3]
      REAL(fp), ALLOCATABLE        :: AIRDEN_DRY(:,:,:) !< Dry air density [kg/m3]
      REAL(fp), ALLOCATABLE        :: AIRNUMDEN(:,:,:)  !< Dry air density [molec/cm3]
      REAL(fp), ALLOCATABLE        :: MAIRDEN(:,:,:)    !< Moist air density (same as AIRDEN to cover possible use cases) [kg/m3]
      REAL(fp), ALLOCATABLE        :: AVGW(:,:,:)       !< Water vapor volume mixing ratio [vol H2O/vol dry air]
      REAL(fp), ALLOCATABLE        :: DELP(:,:,:)       !< Delta-P (wet) across box [Pa]
      REAL(fp), ALLOCATABLE        :: DELP_DRY(:,:,:)   !< Delta-P (dry) across box [Pa]
      REAL(fp), ALLOCATABLE        :: DAIRMASS(:,:,:)   !< Dry air mass [kg] in grid box
      REAL(fp), ALLOCATABLE        :: AIRVOL(:,:,:)     !< Grid box volume [m3] (dry air)
      REAL(fp), ALLOCATABLE        :: PEDGE_DRY(:,:,:)  !< Dry air partial pressure @ level edges [Pa] (nx,ny,nz+1)
      REAL(fp), ALLOCATABLE        :: PEDGE(:,:,:)      !< Air partial pressure @ level edges [Pa] (nx,ny,nz+1)
      REAL(fp), ALLOCATABLE        :: PMID(:,:,:)       !< Average wet air pressure [Pa] defined as arithmetic average of edge pressures
      REAL(fp), ALLOCATABLE        :: PMID_DRY(:,:,:)   !< Dry air partial pressure [Pa] defined as arithmetic avg of edge pressures
   contains
      procedure :: init => metstate_init
      procedure :: cleanup => metstate_cleanup
      procedure :: validate => metstate_validate
      procedure :: reset => metstate_reset
      procedure :: is_allocated => metstate_is_allocated
      procedure :: get_memory_usage => metstate_get_memory_usage
      procedure :: print_summary => metstate_print_summary
      procedure :: get_dimensions => metstate_get_dimensions
      procedure :: get_field_ptr => metstate_get_field_ptr
      procedure :: get_field_ptr_int => metstate_get_field_ptr_int
      procedure :: get_field_ptr_logical => metstate_get_field_ptr_logical
      procedure, public :: get_column_ptr_func => metstate_get_column_ptr_func
      procedure, public :: get_column_ptr_func_int => metstate_get_column_ptr_func_int
      procedure, public :: get_column_ptr_func_logical => metstate_get_column_ptr_func_logical
      procedure, public :: get_column_ptr => metstate_get_column_ptr_subroutine
      procedure, public :: get_2Dto0D_value => metstate_get_2Dto0D_value
      procedure, public :: get_2Dto0D_value_int => metstate_get_2Dto0D_value_int
      procedure, public :: get_2Dto0D_value_logical => metstate_get_2Dto0D_value_logical
      procedure, public :: get_scalar_value => metstate_get_scalar_value
      procedure, public :: get_scalar_value_int => metstate_get_scalar_value_int
      procedure, public :: get_scalar_value_logical => metstate_get_scalar_value_logical
      ! Generic interface for setting fields with proper dimensions
      generic, public :: set_field => metstate_set_field_scalar_real, &
         metstate_set_field_scalar_int, &
         metstate_set_field_scalar_logical, &
         metstate_set_field_2d_real, &
         metstate_set_field_2d_int, &
         metstate_set_field_2d_logical, &
         metstate_set_field_3d_real, &
         metstate_set_field_3d_int, &
         metstate_set_field_3d_logical
      procedure, public :: metstate_set_field_scalar_real
      procedure, public :: metstate_set_field_scalar_int
      procedure, public :: metstate_set_field_scalar_logical
      procedure, public :: metstate_set_field_2d_real
      procedure, public :: metstate_set_field_2d_int
      procedure, public :: metstate_set_field_2d_logical
      procedure, public :: metstate_set_field_3d_real
      procedure, public :: metstate_set_field_3d_int
      procedure, public :: metstate_set_field_3d_logical
      procedure, public :: set_multiple_fields => metstate_set_multiple_fields
      procedure, public :: derive_field => metstate_derive_field
      procedure :: allocate_field => metstate_allocate_field
      procedure :: deallocate_field => metstate_deallocate_field
      procedure, private :: allocate_arrays => allocate_metstate_arrays
   end type MetStateType

CONTAINS

   !> \brief Initialize a MetStateType object
   !!
   !! Initializes the meteorological state object, sets default values, and allocates required arrays.
   !!
   !! \param[inout] this      MetStateType object to initialize
   !! \param[in]    nx        Number of grid points in x direction
   !! \param[in]    ny        Number of grid points in y direction
   !! \param[in]    nlevs     Number of vertical levels
   !! \param[in]    nsoil     Number of soil layers
   !! \param[in]    nsoiltype Number of soil types
   !! \param[in]    nsurftype Number of surface types
   !! \param[inout] error_mgr Error manager for context and error reporting
   !! \param[out]   rc        Return code (CC_SUCCESS or error code)
   subroutine metstate_init(this, nx, ny, nlevs, nsoil, nsoiltype, nsurftype, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, ERROR_MEMORY_ALLOCATION

      implicit none
      class(MetStateType), intent(inout) :: this
      integer, intent(in) :: nx, ny, nlevs
      integer, intent(in), optional :: nsoil, nsoiltype, nsurftype
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc
      character(len=256) :: thisLoc

      thisLoc = 'metstate_init (in core/metstate_mod.F90)'
      call error_mgr%push_context('metstate_init', 'initializing meteorological state')

      rc = CC_SUCCESS

      ! Initialize default values for integer parameters
      this%NLEVS = nlevs

      call this%geometry%set(nx, ny, nlevs) ! Add a set() method to GridGeometryType

      this%State = 'MET'

      ! Set soil and surface parameters if provided
      if (present(nsurftype)) then
         this%NSURFTYPE = nsurftype
      else
         this%NSURFTYPE = 0  ! Will prevent allocation of surface arrays
      end if

      ! Set soil parameters if provided
      if (present(nsoil)) then
         this%nSOIL = nsoil
      else
         this%nSOIL = 0  ! Will prevent allocation of soil arrays
      end if

      if (present(nsoiltype)) then
         this%nSOILTYPE = nsoiltype
      else
         this%nSOILTYPE = 0  ! Will prevent allocation of soil type arrays
      end if

      ! Call helper procedure to allocate arrays
      call this%allocate_arrays('ALL', error_mgr, rc)

      call error_mgr%pop_context()
   end subroutine metstate_init

   !> \brief Allocate all arrays for MetStateType (optionally, only a specific field)
   !!
   !! Helper procedure to allocate and initialize arrays in the meteorological state.
   !!
   !! \param[inout] this      MetStateType object
   !! \param[in]    field_name Name of the field to allocate (or 'ALL' for all fields)
   !! \param[inout] error_mgr Error manager for context and error reporting
   !! \param[out]   rc        Return code
   subroutine allocate_metstate_arrays(this, field_name, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, ERROR_MEMORY_ALLOCATION

      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc
      character(len=256) :: thisLoc
      integer :: nx, ny, nz, nsoil, nsoiltype, nSURFTYPE

      thisLoc = 'allocate_metstate_arrays (in core/metstate_mod.F90)'
      rc = CC_SUCCESS

      call this%geometry%get_dimensions(nx, ny, nz)

      ! Use the properly initialized values (no more defaults needed)
      nsoil = this%nSOIL
      nsoiltype = this%nSOILTYPE
      nSURFTYPE = this%NSURFTYPE

      ! Auto-allocate all REAL(fp), ALLOCATABLE fields (generated, now field-by-field)
      select case (trim(field_name))
#include "metstate_allocate_fields.inc"
      end select

      ! Initialize to safe defaults
      this%T = 288.15_fp
      this%U = 0.0_fp
      this%V = 0.0_fp
      this%QV = 0.001_fp
      this%RH = 0.50_fp
      this%AIRDEN = 1.2_fp
      this%BXHEIGHT = 100.0_fp
      this%PS = 101300.25_fp
      this%SLP = 101300.25_fp
      this%T2M = 288.15_fp
      this%TS = 288.15_fp

   end subroutine allocate_metstate_arrays

   !> \brief Deallocate and clean up all arrays in MetStateType (field-by-field)
   !!
   !! Deallocates the specified allocatable array and resets scalar values in the meteorological state.
   !!
   !! \param[inout] this MetStateType object
   !! \param[in]    field_name Name of the field to deallocate (or 'ALL' for all fields)
   !! \param[out]   rc   Return code (CC_SUCCESS)
   subroutine metstate_cleanup(this, field_name, rc)
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Deallocate only the requested field (or all if 'ALL')
      select case (trim(field_name))
#include "metstate_deallocate_fields.inc"
      end select

      this%State = ''
      this%NLEVS = 72  ! Reset to default
      this%NSURFTYPE = 1  ! Reset to default

   end subroutine metstate_cleanup

   !> \brief Validate the MetStateType object for consistency and physical reasonableness
   !!
   !! Checks that the number of levels is positive, temperatures and pressures are within physical ranges,
   !! and that required arrays are allocated.
   !!
   !! \param[in]    this      MetStateType object
   !! \param[inout] error_mgr Error manager for context and error reporting
   !! \param[out]   rc        Return code (CC_SUCCESS or error code)
   subroutine metstate_validate(this, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, ERROR_INVALID_INPUT
      use utilities_mod, only: is_valid_temperature, is_valid_pressure

      implicit none
      class(MetStateType), intent(in) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      character(len=256) :: thisLoc

      thisLoc = 'metstate_validate (in core/metstate_mod.F90)'
      call error_mgr%push_context('metstate_validate', 'validating meteorological state')

      rc = CC_SUCCESS

      ! Check basic state
      if (this%NLEVS <= 0) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
            'Number of levels must be positive', rc, &
            thisLoc, 'Set NLEVS to a positive integer')
         call error_mgr%pop_context()
         return
      endif

      ! Validate temperatures (use maxval/minval for array validation)
      if (allocated(this%T2M)) then
         if (maxval(this%T2M) > 400.0_fp .or. minval(this%T2M) < 100.0_fp) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               '2m temperature out of physical range', rc, &
               thisLoc, 'Check temperature units and values')
            call error_mgr%pop_context()
            return
         endif
      endif

      if (allocated(this%TS)) then
         if (maxval(this%TS) > 400.0_fp .or. minval(this%TS) < 100.0_fp) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'Surface temperature out of physical range', rc, &
               thisLoc, 'Check temperature units and values')
            call error_mgr%pop_context()
            return
         endif
      endif

      ! Validate pressures (use maxval/minval for array validation)
      if (allocated(this%PS)) then
         if (maxval(this%PS) > 120000.0_fp .or. minval(this%PS) < 1000.0_fp) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'Surface pressure out of physical range', rc, &
               thisLoc, 'Check pressure units and values')
            call error_mgr%pop_context()
            return
         endif
      endif

      if (allocated(this%SLP)) then
         if (maxval(this%SLP) > 120000.0_fp .or. minval(this%SLP) < 50000.0_fp) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'Sea level pressure out of physical range', rc, &
               thisLoc, 'Check pressure units and values')
            call error_mgr%pop_context()
            return
         endif
      endif

      ! Check array allocation
      if (.not. this%is_allocated()) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
            'Required arrays not allocated', rc, &
            thisLoc, 'Call init() before using MetState')
         call error_mgr%pop_context()
         return
      endif

      call error_mgr%pop_context()
   end subroutine metstate_validate

   !> \brief Reset MetStateType to initial values
   !!
   !! Resets time, surface fields, and arrays to standard atmosphere values.
   !!
   !! \param[inout] this MetStateType object
   !! \param[out]   rc   Return code (CC_SUCCESS)
   subroutine metstate_reset(this, rc)
      implicit none
      class(MetStateType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Reset to standard atmosphere values
      if (allocated(this%T2M)) this%T2M = 288.15_fp
      if (allocated(this%TS)) this%TS = 288.15_fp
      if (allocated(this%TSKIN)) this%TSKIN = 288.15_fp
      if (allocated(this%PS)) this%PS = 101300.25_fp
      if (allocated(this%SLP)) this%SLP = 101300.25_fp
      if (allocated(this%SST)) this%SST = 288.15_fp

      ! Reset arrays if allocated
      if (allocated(this%T)) this%T = 288.15_fp
      if (allocated(this%U)) this%U = 0.0_fp
      if (allocated(this%V)) this%V = 0.0_fp
      if (allocated(this%QV)) this%QV = 0.01_fp
      if (allocated(this%RH)) this%RH = 0.5_fp

   end subroutine metstate_reset

   !> \brief Check if required arrays are allocated in MetStateType
   !!
   !! Returns .true. if all required arrays are allocated, .false. otherwise.
   !!
   !! \param[in] this MetStateType object
   !! \return    Logical flag indicating allocation status
   function metstate_is_allocated(this) result(is_alloc)
      implicit none
      class(MetStateType), intent(in) :: this
      logical :: is_alloc

      is_alloc = allocated(this%T) .and. allocated(this%U) .and. allocated(this%V) .and. &
         allocated(this%QV) .and. allocated(this%PMID) .and. allocated(this%DELP)
   end function metstate_is_allocated

   !> \brief Get approximate memory usage of MetStateType in bytes
   !!
   !! Estimates the memory usage of all arrays in the meteorological state object.
   !!
   !! \param[in] this MetStateType object
   !! \return    Estimated memory usage in bytes (integer(kind=8))
   function metstate_get_memory_usage(this) result(memory_bytes)
      implicit none
      class(MetStateType), intent(in) :: this
      integer(kind=8) :: memory_bytes

      integer :: nlevs

      memory_bytes = 0
      nlevs = this%NLEVS

      if (nlevs > 0) then
         ! Estimate based on number of allocated arrays and precision
         ! Each real(fp) array: nlevs * 8 bytes (assuming fp = real64)
         ! Each logical array: nlevs * 1 byte
         memory_bytes = nlevs * 8 * 26  ! 26 real arrays
         memory_bytes = memory_bytes + nlevs * 1 * 4  ! 4 logical arrays
         memory_bytes = memory_bytes + (nlevs+1) * 8 * 2  ! 2 edge arrays
      endif
   end function metstate_get_memory_usage

   !> \brief Print a summary of the MetStateType object to standard output
   !!
   !! Prints key fields, allocation status, and memory usage for diagnostics.
   !!
   !! \param[in] this MetStateType object
   subroutine metstate_print_summary(this)
      implicit none
      class(MetStateType), intent(in) :: this

      write(*,'(A)') '=== MetState Summary ==='
      write(*,'(A,A)') 'State: ', trim(this%State)
      write(*,'(A,I0)') 'Number of levels: ', this%NLEVS
      if (allocated(this%TS)) then
         write(*,'(A,F8.2,A)') 'Surface temperature: ', this%TS(1,1), ' K'
      endif
      if (allocated(this%PS)) then
         write(*,'(A,F8.2,A)') 'Surface pressure: ', this%PS(1,1), ' Pa'
      endif
      if (allocated(this%SLP)) then
         write(*,'(A,F8.2,A)') 'Sea level pressure: ', this%SLP(1,1), ' Pa'
      endif
      write(*,'(A,L1)') 'Arrays allocated: ', this%is_allocated()
      write(*,'(A,I0,A)') 'Memory usage: ', this%get_memory_usage(), ' bytes'
      write(*,'(A)') '======================='
   end subroutine metstate_print_summary

   !========================================================================
   !! Get grid dimensions for column interface support
   !========================================================================
   subroutine metstate_get_dimensions(this, nx, ny, nlev)
      class(MetStateType), intent(in) :: this
      integer, intent(out) :: nx, ny, nlev

      ! Get actual dimensions from geometry
      call this%geometry%get_dimensions(nx, ny, nlev)

   end subroutine metstate_get_dimensions

!    !========================================================================
!    !! Get pointer to a meteorological field array by name
!    !!
!    !! Returns a pointer to the requested 1D field array (e.g., temperature, pressure) for the column model.
!    !! Returns null() if the field is not found or not allocated.
!    !!
!    !! \param[in]  this       MetStateType object
!    !! \param[in]  field_name Name of the field (e.g., 'T', 'temperature')
!    !! \return     Pointer to the requested field array, or null()
!    function metstate_get_field_ptr(this, field_name) result(field_ptr)
!       implicit none
!       class(MetStateType), intent(in), target :: this
!       character(len=*), intent(in) :: field_name
!       real(fp), pointer :: field_ptr(:)

!       ! Return pointer to 1D column data
!       field_ptr => null()

!       select case (trim(field_name))
! #include "metstate_accessor.inc"
!       end select

!    end function metstate_get_field_ptr

   !> \brief Allocate a specific field in MetStateType by name
   !!
   !! Calls the generated select-case macro to allocate only the requested field.
   !!
   !! \param[inout] this      MetStateType object
   !! \param[in]    field_name Name of the field to allocate
   !! \param[out]   rc         Return code (CC_SUCCESS or error code)
   subroutine metstate_allocate_field(this, field_name, rc)
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(out) :: rc
      integer :: nx, ny, nz, nsoil, nsoiltype, nSURFTYPE
      rc = CC_SUCCESS
      call this%geometry%get_dimensions(nx, ny, nz)
      ! Use the properly initialized values (no more defaults needed)
      nsoil = this%nSOIL
      nsoiltype = this%nSOILTYPE
      nSURFTYPE = this%NSURFTYPE
      ! Only allocate the requested field
      select case (trim(field_name))
#include "metstate_allocate_fields.inc"
      end select
   end subroutine metstate_allocate_field

   !> \brief Deallocate a specific field in MetStateType by name
   !!
   !! Calls the generated select-case macro to deallocate only the requested field.
   !!
   !! \param[inout] this      MetStateType object
   !! \param[in]    field_name Name of the field to deallocate
   !! \param[out]   rc         Return code (CC_SUCCESS or error code)
   subroutine metstate_deallocate_field(this, field_name, rc)
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(out) :: rc
      rc = CC_SUCCESS
      ! Only deallocate the requested field
      select case (trim(field_name))
#include "metstate_deallocate_fields.inc"
      end select
   end subroutine metstate_deallocate_field

   !> \brief Get a pointer to a vertical column for a given field name and (i,j) indices (type-safe)
   !!
   !! Calls the type-safe function version to obtain a real(fp) pointer, then assigns it to a polymorphic pointer for generic access.
   !! If the field is not found or not allocated, col_ptr is set to null and rc is set to CC_FAILURE.
   !!
   !! \param[inout] this      MetStateType object
   !! \param[in]    field_name Name of the field (e.g., 'T', 'temperature')
   !! \param[in]    i          Grid column index (1-based)
   !! \param[in]    j          Grid row index (1-based)
   !! \param[out]   col_ptr    Pointer to the vertical column data (polymorphic pointer)
   !! \param[out]   rc         Return code (CC_SUCCESS if found, CC_FAILURE otherwise)
   subroutine metstate_get_3Dto1D_ptr(this, field_name, i, j, col_ptr, rc)
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: i, j
      real(fp), pointer :: col_ptr(:)
      integer, intent(out) :: rc
      col_ptr => this%get_column_ptr_func(field_name, i, j)
      if (associated(col_ptr)) then
         rc = 0
      else
         rc = 1
      endif
   end subroutine metstate_get_3Dto1D_ptr

   !========================================================================
   !! Get pointer to a vertical column for a given field at (i,j)
   !! Returns a pointer to the vertical profile for a given field at grid location (i,j).
   !! For column models, returns the full 1D array.
   !! \param[in]  this       MetStateType object
   !! \param[in]  field_name Name of the field (e.g., 'T', 'temperature')
   !! \param[in]  i          Grid column index (optional, default 1)
   !! \param[in]  j          Grid row index (optional, default 1)
   !! \return     Pointer to vertical profile (1D)
   function metstate_get_column_ptr_func(this, field_name, i, j) result(column_ptr)
      implicit none
      class(MetStateType), intent(in), target :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in), optional :: i, j
      real(fp), pointer :: column_ptr(:)
      integer :: nx, ny, nlev, col_i, col_j
      column_ptr => null()
      call this%get_dimensions(nx, ny, nlev)
      col_i = 1; col_j = 1
      if (present(i)) col_i = max(1, min(i, nx))
      if (present(j)) col_j = max(1, min(j, ny))
      select case (trim(field_name))
#include "metstate_column_accessor.inc"
      end select
   end function metstate_get_column_ptr_func

   !> Get a scalar value from a 2D field at (i,j)
   function metstate_get_2Dto0D_value(this, field_name, i, j) result(scalar_val)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: i, j
      real(fp) :: scalar_val
      integer :: col_i, col_j
      col_i = i
      col_j = j
      select case (trim(field_name))
#include "metstate_2d_scalar_accessor.inc"
      end select
   end function metstate_get_2Dto0D_value

   !> Get a scalar value from a scalar field
   function metstate_get_scalar_value(this, field_name) result(scalar_val)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      real(fp) :: scalar_val
      select case (trim(field_name))
#include "metstate_scalar_accessor.inc"
      end select
   end function metstate_get_scalar_value

   !> INTEGER versions of accessor functions
   function metstate_get_column_ptr_func_int(this, field_name, i, j) result(column_ptr)
      implicit none
      class(MetStateType), intent(in), target :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in), optional :: i, j
      integer, pointer :: column_ptr(:)
      integer :: nx, ny, nlev, col_i, col_j
      column_ptr => null()
      call this%get_dimensions(nx, ny, nlev)
      col_i = 1; col_j = 1
      if (present(i)) col_i = max(1, min(i, nx))
      if (present(j)) col_j = max(1, min(j, ny))
      select case (trim(field_name))
#include "metstate_column_accessor_int.inc"
      end select
   end function metstate_get_column_ptr_func_int

   function metstate_get_2Dto0D_value_int(this, field_name, i, j) result(scalar_val)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: i, j
      integer :: scalar_val
      integer :: col_i, col_j
      col_i = i
      col_j = j
      select case (trim(field_name))
#include "metstate_2d_scalar_accessor_int.inc"
      end select
   end function metstate_get_2Dto0D_value_int

   function metstate_get_scalar_value_int(this, field_name) result(scalar_val)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      integer :: scalar_val
      select case (trim(field_name))
#include "metstate_scalar_accessor_int.inc"
      end select
   end function metstate_get_scalar_value_int

   !> LOGICAL versions of accessor functions
   function metstate_get_column_ptr_func_logical(this, field_name, i, j) result(column_ptr)
      implicit none
      class(MetStateType), intent(in), target :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in), optional :: i, j
      logical, pointer :: column_ptr(:)
      integer :: nx, ny, nlev, col_i, col_j
      column_ptr => null()
      call this%get_dimensions(nx, ny, nlev)
      col_i = 1; col_j = 1
      if (present(i)) col_i = max(1, min(i, nx))
      if (present(j)) col_j = max(1, min(j, ny))
      select case (trim(field_name))
#include "metstate_column_accessor_logical.inc"
      end select
   end function metstate_get_column_ptr_func_logical

   function metstate_get_2Dto0D_value_logical(this, field_name, i, j) result(scalar_val)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: i, j
      logical :: scalar_val
      integer :: col_i, col_j
      col_i = i
      col_j = j
      select case (trim(field_name))
#include "metstate_2d_scalar_accessor_logical.inc"
      end select
   end function metstate_get_2Dto0D_value_logical

   function metstate_get_scalar_value_logical(this, field_name) result(scalar_val)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      logical :: scalar_val
      select case (trim(field_name))
#include "metstate_scalar_accessor_logical.inc"
      end select
   end function metstate_get_scalar_value_logical

   !> High-level interface: get any field (column, 2D, or scalar)
   subroutine metstate_get_field_ptr(this, field_name, i, j, col_ptr, scalar_val, rc)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in), optional :: i, j
      real(fp), pointer, optional :: col_ptr(:)
      real(fp), optional :: scalar_val
      integer, intent(out) :: rc
      ! Try 3D column first
      if (present(col_ptr) .and. present(i) .and. present(j)) then
         col_ptr => this%get_column_ptr_func(field_name, i, j)
         if (associated(col_ptr)) then
            rc = 0
            return
         end if
      end if
      ! Try 2D scalar
      if (present(scalar_val) .and. present(i) .and. present(j)) then
         scalar_val = this%get_2Dto0D_value(field_name, i, j)
         rc = 0
         return
      end if
      ! Try scalar field
      if (present(scalar_val)) then
         scalar_val = this%get_scalar_value(field_name)
         rc = 0
         return
      end if
      rc = 1 ! Not found
   end subroutine metstate_get_field_ptr

   !> Integer version of get_field_ptr
   subroutine metstate_get_field_ptr_int(this, field_name, i, j, col_ptr, scalar_val, rc)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in), optional :: i, j
      integer, pointer, optional :: col_ptr(:)
      integer, optional :: scalar_val
      integer, intent(out) :: rc
      ! Try 3D column first
      if (present(col_ptr) .and. present(i) .and. present(j)) then
         col_ptr => this%get_column_ptr_func_int(field_name, i, j)
         if (associated(col_ptr)) then
            rc = 0
            return
         end if
      end if
      ! Try 2D scalar
      if (present(scalar_val) .and. present(i) .and. present(j)) then
         scalar_val = this%get_2Dto0D_value_int(field_name, i, j)
         rc = 0
         return
      end if
      ! Try scalar field
      if (present(scalar_val)) then
         scalar_val = this%get_scalar_value_int(field_name)
         rc = 0
         return
      end if
      rc = 1 ! Not found
   end subroutine metstate_get_field_ptr_int

   !> Logical version of get_field_ptr
   subroutine metstate_get_field_ptr_logical(this, field_name, i, j, col_ptr, scalar_val, rc)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in), optional :: i, j
      logical, pointer, optional :: col_ptr(:)
      logical, optional :: scalar_val
      integer, intent(out) :: rc
      ! Try 3D column first
      if (present(col_ptr) .and. present(i) .and. present(j)) then
         col_ptr => this%get_column_ptr_func_logical(field_name, i, j)
         if (associated(col_ptr)) then
            rc = 0
            return
         end if
      end if
      ! Try 2D scalar
      if (present(scalar_val) .and. present(i) .and. present(j)) then
         scalar_val = this%get_2Dto0D_value_logical(field_name, i, j)
         rc = 0
         return
      end if
      ! Try scalar field
      if (present(scalar_val)) then
         scalar_val = this%get_scalar_value_logical(field_name)
         rc = 0
         return
      end if
      rc = 1 ! Not found
   end subroutine metstate_get_field_ptr_logical

   !> \brief Get a pointer to a vertical column for a given field name and (i,j) indices (subroutine version)
   !!
   !! This subroutine version provides the interface expected by StateManager_Mod.
   !! It handles different variable types (2D fields, 3D fields, scalar values) and returns
   !! a 1D pointer to the appropriate data.
   !!
   !! \param[inout] this      MetStateType object
   !! \param[in]    field_name Name of the field (e.g., 'T', 'temperature', 'PS')
   !! \param[in]    i          Grid column index (1-based)
   !! \param[in]    j          Grid row index (1-based)
   !! \param[out]   col_ptr    Pointer to the vertical column data (1D array)
   !! \param[out]   rc         Return code (CC_SUCCESS if found, CC_FAILURE otherwise)
   subroutine metstate_get_column_ptr_subroutine(this, field_name, i, j, col_ptr, rc)
      use error_mod, only: CC_SUCCESS, CC_FAILURE

      implicit none
      class(MetStateType), intent(inout), target :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: i, j
      real(fp), pointer :: col_ptr(:)
      integer, intent(out) :: rc

      ! Local variables for handling different field types
      real(fp), pointer :: temp_col_ptr(:)
      real(fp) :: scalar_val
      integer :: nx, ny, nlev

      rc = CC_FAILURE
      nullify(col_ptr)

      call this%get_dimensions(nx, ny, nlev)

      ! First try to get as a 3D field (vertical column)
      temp_col_ptr => this%get_column_ptr_func(field_name, i, j)
      if (associated(temp_col_ptr)) then
         col_ptr => temp_col_ptr
         rc = CC_SUCCESS
         return
      endif

      ! If not found as 3D field, try as 2D field and create a single-element array
      ! For 2D fields, we return a pointer to a single-element array containing the scalar value
      select case (trim(field_name))
       case ('PS', 'SLP', 'TS', 'T2M', 'TSKIN', 'SST', 'PHIS', 'PS_WET', 'PS_DRY', &
          'QV2M', 'AREA_M2', 'ALBD_VIS', 'ALBD_NIR', 'ALBD_UV', 'PARDR', 'PARDF', &
          'SUNCOS', 'SUNCOSmid', 'SWGDN', 'EFLUX', 'HFLUX', 'U10M', 'V10M', &
          'USTAR', 'Z0', 'Z0H', 'PBLH', 'OBK', 'CLDFRC', 'CONV_DEPTH', &
          'FLASH_DENS', 'CNV_FRC', 'PRECANV', 'PRECCON', 'PRECLSC', &
          'LAI', 'GVF', 'RDRAG', 'CLAYFRAC', 'SANDFRAC', 'FRVEG', 'FRLAKE', &
          'FRLAND', 'FRLANDIC', 'FROCEAN', 'FRSEAICE', 'FRSNO', 'SNODP', &
          'SNOMAS', 'SSM', 'USTAR_THRESHOLD', 'GWETTOP', 'GWETROOT', 'WILT', &
          'TO3', 'TROPP', 'TropHt', 'LAT', 'LON')

         scalar_val = this%get_2Dto0D_value(field_name, i, j)

         ! For 2D fields, we need to return a pointer to a single element
         ! This is a limitation of the current interface - we can't easily return a scalar as a 1D pointer
         ! For now, we'll return null and indicate failure
         rc = CC_FAILURE
         return

       case default
         ! Try as scalar field - similar limitation
         scalar_val = this%get_scalar_value(field_name)
         rc = CC_FAILURE
         return
      end select

   end subroutine metstate_get_column_ptr_subroutine


   !---------------------------------------------------------------------------
   !                 Dimensional MetState Set Field Subroutines
   !---------------------------------------------------------------------------

   !> @brief Set a scalar REAL field
   subroutine metstate_set_field_scalar_real(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: field_data
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Handle scalar REAL fields directly
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_scalar_real.inc"
       case default
         ! If not a scalar field, try broadcasting to 2D REAL arrays
         select case (trim(adjustl(field_name)))
#include "metstate_set_field_2d_real.inc"
          case default
            ! Try broadcasting to 3D REAL arrays
            select case (trim(adjustl(field_name)))
#include "metstate_set_field_3d_real.inc"
             case default
               call error_mgr%report_error(ERROR_NOT_FOUND, &
                  'Unknown REAL field name: ' // trim(field_name), rc)
               rc = CC_FAILURE
            end select
         end select
      end select
   end subroutine metstate_set_field_scalar_real

   !> @brief Set a scalar INTEGER field
   subroutine metstate_set_field_scalar_int(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: field_data
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Handle scalar INTEGER fields directly
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_scalar_int.inc"
       case default
         ! If not a scalar field, try broadcasting to 2D INTEGER arrays
         select case (trim(adjustl(field_name)))
#include "metstate_set_field_2d_int.inc"
          case default
            ! Try broadcasting to 3D INTEGER arrays
            select case (trim(adjustl(field_name)))
#include "metstate_set_field_3d_int.inc"
             case default
               call error_mgr%report_error(ERROR_NOT_FOUND, &
                  'Unknown INTEGER field name: ' // trim(field_name), rc)
               rc = CC_FAILURE
            end select
         end select
      end select
   end subroutine metstate_set_field_scalar_int

   !> @brief Set a scalar LOGICAL field
   subroutine metstate_set_field_scalar_logical(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      logical, intent(in) :: field_data
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Handle scalar LOGICAL fields directly
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_scalar_logical.inc"
       case default
         ! If not a scalar field, try broadcasting to 2D LOGICAL arrays
         select case (trim(adjustl(field_name)))
#include "metstate_set_field_2d_logical.inc"
          case default
            ! Try broadcasting to 3D LOGICAL arrays
            select case (trim(adjustl(field_name)))
#include "metstate_set_field_3d_logical.inc"
             case default
               call error_mgr%report_error(ERROR_NOT_FOUND, &
                  'Unknown LOGICAL field name: ' // trim(field_name), rc)
               rc = CC_FAILURE
            end select
         end select
      end select
   end subroutine metstate_set_field_scalar_logical

   !> @brief Set a 2D REAL field
   subroutine metstate_set_field_2d_real(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: field_data(:,:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Generated include file for 2D REAL field assignment
      rc = CC_SUCCESS
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_2d_real.inc"
       case default
         call error_mgr%report_error(ERROR_NOT_FOUND, &
            "Unknown field name: " // trim(field_name), rc)
         rc = CC_FAILURE
      end select
   end subroutine metstate_set_field_2d_real

   !> @brief Set a 2D INTEGER field
   subroutine metstate_set_field_2d_int(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: field_data(:,:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Generated include file for 2D INTEGER field assignment
      rc = CC_SUCCESS
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_2d_int.inc"
       case default
         call error_mgr%report_error(ERROR_NOT_FOUND, &
            "Unknown field name: " // trim(field_name), rc)
         rc = CC_FAILURE
      end select
   end subroutine metstate_set_field_2d_int

   !> @brief Set a 2D LOGICAL field
   subroutine metstate_set_field_2d_logical(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      logical, intent(in) :: field_data(:,:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Generated include file for 2D LOGICAL field assignment
      rc = CC_SUCCESS
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_2d_logical.inc"
       case default
         call error_mgr%report_error(ERROR_NOT_FOUND, &
            "Unknown field name: " // trim(field_name), rc)
         rc = CC_FAILURE
      end select
   end subroutine metstate_set_field_2d_logical

   !> @brief Set a 3D REAL field
   subroutine metstate_set_field_3d_real(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: field_data(:,:,:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Generated include file for 3D REAL field assignment
      rc = CC_SUCCESS
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_3d_real.inc"
       case default
         call error_mgr%report_error(ERROR_NOT_FOUND, &
            "Unknown field name: " // trim(field_name), rc)
         rc = CC_FAILURE
      end select
   end subroutine metstate_set_field_3d_real

   !> @brief Set a 3D INTEGER field
   subroutine metstate_set_field_3d_int(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: field_data(:,:,:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Generated include file for 3D INTEGER field assignment
      rc = CC_SUCCESS
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_3d_int.inc"
       case default
         call error_mgr%report_error(ERROR_NOT_FOUND, &
            "Unknown field name: " // trim(field_name), rc)
         rc = CC_FAILURE
      end select
   end subroutine metstate_set_field_3d_int

   !> @brief Set a 3D LOGICAL field
   subroutine metstate_set_field_3d_logical(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      logical, intent(in) :: field_data(:,:,:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Generated include file for 3D LOGICAL field assignment
      rc = CC_SUCCESS
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_3d_logical.inc"
       case default
         call error_mgr%report_error(ERROR_NOT_FOUND, &
            "Unknown field name: " // trim(field_name), rc)
         rc = CC_FAILURE
      end select
   end subroutine metstate_set_field_3d_logical

! Include the auto-generated multiple fields interface
#include "metstate_multiple_fields_interface.inc"

   !> \brief Derive meteorological fields from existing data
   !!
   !! Calculates derived fields using existing meteorological variables.
   !! Supports common derived quantities like air density, virtual temperature, etc.
   !!
   !! \param[inout] this        MetStateType object
   !! \param[in]    field_name  Name of the field to derive
   !! \param[inout] error_mgr   Error manager for context and error reporting
   !! \param[out]   rc          Return code (CC_SUCCESS or error code)
   subroutine metstate_derive_field(this, field_name, error_mgr, time_state, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE, ERROR_INVALID_INPUT, ERROR_NOT_FOUND
      use constants, only: g0, Rd, Rdg0, AIRMW, H2OMW

      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      type(TimeStateType), pointer,intent(inout) :: time_state
      integer, intent(out) :: rc

      character(len=256) :: thisLoc
      integer :: nx, ny, nz, i, j, k, nlanduse
      real(fp) :: airden
      real(fp) :: avgw ! Water vapor volume mixing ratio [v/v dry air]
      real(fp) :: xh2o ! Water vapor mole fraction [mol (H2O) / mol (moist air)]

      thisLoc = 'metstate_derive_field (in core/metstate_mod.F90)'
      call error_mgr%push_context('metstate_derive_field', 'deriving field: ' // trim(field_name))

      rc = CC_SUCCESS
      call this%get_dimensions(nx, ny, nz)

      select case (trim(adjustl(field_name)))

       case ('MAIRDEN', 'mairden', 'AIRDEN', 'airden')
         ! Calculate dry air density from pressure and temperature
         ! ρ = P / (R_specific * T) where R_specific = R / MW
         if (.not. allocated(this%PMID) .or. .not. allocated(this%T)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'PMID and T fields required for MAIRDEN/AIRDEN calculation', rc, &
               thisLoc, 'Ensure pressure and temperature are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate MAIRDEN if not already allocated
         if (.not. allocated(this%MAIRDEN) .or. .not. allocated(this%AIRDEN)) then
            call error_mgr%report_error(rc, 'MAIRDEN/AIRDEN fields need to be allocated first!', rc, thisLoc)
            call error_mgr%pop_context()
            return
         endif

         ! Calculate dry air density: ρ = P / (R_dry * T)
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  this%MAIRDEN(i, j, k) = this%PMID(i, j, k) / rd / this%T(i, j, k)
                  this%AIRDEN(i, j, k) = this%PMID(i, j, k) / rd / this%T(i, j, k)
               enddo
            enddo
         enddo

       case ('AIRDEN_DRY', 'airden_dry', 'PMID_DRY', 'pmid_dry', 'PEDGE_DRY', 'pedge_dry', 'DELP_DRY', 'delp_dry')
         ! Calculate dry air density from pressure and temperature
         ! ρ = P / (R_specific * T) where R_specific = R / MW
         if (.not. allocated(this%PMID) .or. .not. allocated(this%T)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'PMID and T fields required for AIRDEN_DRY calculation', rc, &
               thisLoc, 'Ensure pressure and temperature are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate AIRDEN_DRY if not already allocated
         if (.not. allocated(this%AIRDEN_DRY) .or. .not. allocated(this%PMID_DRY) .or. &
            .not. allocated(this%PEDGE_DRY) .or. .not. allocated(this%DELP_DRY)) then
            call error_mgr%report_error(rc, 'AIRDEN_DRY/PMID_DRY/PEDGE_DRY/DELP_DRY fields need to be allocated first!', rc, thisLoc)
            call error_mgr%pop_context()
            return
         endif

         ! Calculate dry air density: ρ = P / (R_dry * T)
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  avgw = AIRMW * this%QV(i,j,k) / ( H2OMW * (1.0e+0_fp - this%QV(i,j,k)) )
                  xh2o = avgw / (1.0e+0_fp + avgw)
                  this%PMID_DRY(i, j, k) = this%PMID(i, j, k) * ( 1.e+0_fp - xh2o )
                  this%AIRDEN_DRY(i, j, k) = this%PMID_DRY(i, j, k) / rd / this%T(i, j, k)
                  this%PEDGE_DRY(i, j, k) = this%PEDGE(i, j, k) * ( 1.e+0_fp - xh2o )
                  if (k == nz) then
                     this%PEDGE_DRY(i, j, k+1) = this%PEDGE(i, j, k+1) * ( 1.e+0_fp - xh2o )
                  end if
                  this%DELP_DRY(i, j, k) = this%PEDGE_DRY(i, j, k) - this%PEDGE_DRY(i, j, k+1)
               enddo
            enddo
         enddo

       case ('RH', 'rh')
         ! Calculate virtual temperature from temperature and humidity
         if (.not. allocated(this%T) .or. .not. allocated(this%QV) .or. .not. allocated(this%PMID)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'T, PMID and QV fields required for RH calculation', rc, &
               thisLoc, 'Ensure temperature, pressure and humidity are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate RH if not already allocated
         if (.not. allocated(this%RH)) then
            call error_mgr%report_error(rc, 'RH field needs to be allocated first!', rc, thisLoc)
            call error_mgr%pop_context()
            return
         endif

         ! Calculate relative humidity from met_utility module
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  this%RH(i, j, k) = relative_humidity(this%T(i, j, k), this%QV(i, j, k), this%PMID(i, j, k))
               enddo
            enddo
         enddo

       case ('TV', 'tv')
         ! Calculate virtual temperature from temperature and humidity
         if (.not. allocated(this%T) .or. .not. allocated(this%QV)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'T and QV fields required for TV calculation', rc, &
               thisLoc, 'Ensure temperature and humidity are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate TV if not already allocated
         if (.not. allocated(this%TV)) then
            call error_mgr%report_error(rc, 'TV field needs to be allocated first!', rc, thisLoc)
            call error_mgr%pop_context()
            return
         endif

         ! Calculate virtual temperature: Tv = T * (1 + 0.608 * qv)
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  this%TV(i, j, k) = this%T(i, j, k) * (1.0_fp + 0.608_fp * this%QV(i, j, k))
               enddo
            enddo
         enddo

       case ('OBK', 'obk')
         ! Calculate OBK from sensible heat flux and air density
         if (.not. allocated(this%HFLUX) .or. .not. allocated(this%AIRDEN) .or. .not. allocated(this%TS) .or. &
            .not. allocated(this%USTAR)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'TS, USTAR, AIRDEN and HFLUX fields required for OBK calculation', rc, &
               thisLoc, 'Ensure temperature, ustar, air density, and sensible heat flux are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate OBK if not already allocated
         if (.not. allocated(this%OBK)) then
            call error_mgr%report_error(rc, 'OBK field needs to be allocated first!', rc, thisLoc)
            call error_mgr%pop_context()
            return
         endif

         ! Calculate OBK from met_utility module
         do j = 1, ny
            do i = 1, nx
               airden = this%PMID(i, j, 1) / rd / this%T(i, j, 1)
               !!!! Note we cannot use this%AIRDEN here because it may not be calculated yet
               this%OBK(i, j) = monin_obukhov_length(this%USTAR(i, j), this%TS(i, j), this%HFLUX(i, j), airden)
            enddo
         enddo

       case ('SUNCOS', 'suncos')
         ! Calculate SUNCOS
         if (.not. allocated(this%LAT) .or. .not. allocated(this%LON)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'LAT and LON fields required for SUNCOS calculation', rc, &
               thisLoc, 'Ensure latitude and longitude are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate OBK if not already allocated
         if (.not. allocated(this%SUNCOS)) then
            call error_mgr%report_error(rc, 'SUNCOS field needs to be allocated first!', rc, thisLoc)
            call error_mgr%pop_context()
            return
         endif

         ! Calculate OBK from met_utility module
         do j = 1, ny
            do i = 1, nx
               !make sure lat[-90 - 90] and lon[-180 - 180] are in degrees
               this%SUNCOS(i, j) = time_state%get_cos_sza(this%LAT(i, j), this%LON(i, j))
            enddo
         enddo

       case ('SUNCOSmid', 'suncosmid')
         ! Calculate SUNCOSmid
         if (.not. allocated(this%LAT) .or. .not. allocated(this%LON)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'LAT and LON fields required for SUNCOSmid calculation', rc, &
               thisLoc, 'Ensure latitude and longitude are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate OBK if not already allocated
         if (.not. allocated(this%SUNCOSmid)) then
            call error_mgr%report_error(rc, 'SUNCOSmid field needs to be allocated first!', rc, thisLoc)
            call error_mgr%pop_context()
            return
         endif

         ! Calculate OBK from met_utility module
         do j = 1, ny
            do i = 1, nx
               !make sure lat[-90 - 90] and lon[-180 - 180] are in degrees
               this%SUNCOSmid(i, j) = time_state%get_cos_sza(this%LAT(i, j), this%LON(i, j), .true.)
            enddo
         enddo

       case ('DELP', 'delp')
         ! Calculate box height from geopotential heights
         if (.not. allocated(this%PEDGE)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'PEDGE field required for DELP calculation', rc, &
               thisLoc, 'Ensure pressure edges are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate BXHEIGHT if not already allocated
         if (.not. allocated(this%DELP)) then
            call error_mgr%report_error(rc, 'BXHEIGHT field needs to be allocated first!', rc, thisLoc)
            call error_mgr%pop_context()
            return
         endif

         ! Calculate box height as difference between edge heights
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  ! lower edge - upper edge
                  this%DELP(i, j, k) = this%PEDGE(i, j, k) - this%PEDGE(i, j, k+1)
               enddo
            enddo
         enddo

       case ('BXHEIGHT', 'bxheight')
         ! Calculate box height from geopotential heights
         if (.not. allocated(this%PEDGE)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'PEDGE field required for BXHEIGHT calculation', rc, &
               thisLoc, 'Ensure pressure edges are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate BXHEIGHT if not already allocated
         if (.not. allocated(this%BXHEIGHT)) then
            call error_mgr%report_error(rc, 'BXHEIGHT field needs to be allocated first!', rc, thisLoc)
            call error_mgr%pop_context()
            return
         endif

         ! Calculate box height as difference between edge heights
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  ! Refer to https://github.com/geoschem/geos-chem/GeosCore/calc_met_mod.F90
                  this%BXHEIGHT(i, j, k) = Rdg0 * virtual_temperature(this%T(i, j, k), this%QV(i, j, k)) * &
                     LOG(this%PEDGE(i, j, k) / this%PEDGE(i, j, k+1))
               enddo
            enddo
         enddo

       case ('SST', 'sst')
         this%SST(:,:) = this%TS(:,:)  !just copy TS to SST

       case ('TSKIN', 'tskin')
         this%TSKIN(:,:) = this%TS(:,:)  !just copy TS to TSKIN

       case ('Z0H', 'z0h')
         this%Z0H(:,:) = this%Z0(:,:)  !just copy Z0 to Z0H

       case ('CLDFRC', 'cldfrc')
         this%CLDFRC(:,:) = this%CLDF(:,:, 1)  !just copy surface CLDF to CLDFRC

       case ('IsLand', 'island', 'ISLAND')
         do j = 1, ny
            do i = 1, nx
               this%IsLand(i, j) = ( abs(this%LWI(i, j) - 1.0_fp) < 0.5_fp ) ! Land if LWI = 1.0
            enddo
         enddo

       case ('IsIce', 'isice', 'ISICE')
         do j = 1, ny
            do i = 1, nx
               this%IsIce(i, j) = ( abs(this%LWI(i, j) - 2.0_fp) < 0.5_fp ) ! Ice if LWI = 2.0
            enddo
         enddo

       case ('IsWater', 'iswater', 'ISWATER')
         do j = 1, ny
            do i = 1, nx
               this%IsWater(i, j) = ( abs(this%LWI(i, j) - 0.0_fp) < 0.5_fp ) ! sea if LWI = 0.0
            enddo
         enddo

       case ('IsSnow', 'issnow', 'ISSNOW')
         do j = 1, ny
            do i = 1, nx
               !geos-chem has a different method: https://github.com/geoschem/geos-chem/GeosCore/calc_met_mod.F90#L324
               this%IsSnow(i, j) = ( this%FRSNO(i, j) >= 0.5_fp ) ! Snow fraction is read in
            enddo
         enddo

       case ('LUCNAME', 'lucname')
         this%LUCNAME = 'NOAH'
       case ('nLNDTYPE', 'nlndtype', 'NLNDTYPE')
         nlanduse = 20  !set to 20 for now; later we can read from a config file or pass in from outside
         this%nLNDTYPE(:,:) = nlanduse  !manually set to 20 for now; not sure if NUOPC can get it
       case ('FRLANDUSE', 'frlanduse')
         !Note that FRLANDUSE is not allocated yet in met_sate%init phase because we don't know nlanduse yet
         nlanduse = 20  !set to 20 for now; later we can read from a config file or pass in from outside
         if (.not. allocated(this%FRLANDUSE)) allocate(this%FRLANDUSE(nx, ny, nlanduse))
         this%FRLANDUSE(:,:,:) = 0.0_fp
         do j = 1, ny
            do i = 1, nx
               do k = 1, nlanduse
                  if (this%DLUSE(i, j) == k) this%FRLANDUSE(i, j, k) = 1.0_fp
                  !We receive DLUSE = 0 over water but it should be 17th type
                  if (this%DLUSE(i, j) == 0 .and. k == 17) this%FRLANDUSE(i, j, k) = 1.0_fp
               enddo
            enddo
         enddo
       case ('ILAND', 'iland')
         !Note that ILAND is not allocated yet in met_sate%init phase because we don't know nlanduse yet
         nlanduse = 20  !set to 20 for now; later we can read from a config file or pass in from outside
         if (.not. allocated(this%ILAND)) allocate(this%ILAND(nx, ny, nlanduse))
         this%ILAND(:,:,:) = 0
         do j = 1, ny
            do i = 1, nx
               do k = 1, nlanduse
                  this%ILAND(i, j, k) = k
               enddo
            enddo
         enddo
       case ('FRLAI', 'frlai')
         !Note that FRLAI is not allocated yet in met_sate%init phase because we don't know nlanduse yet
         nlanduse = 20  !set to 20 for now; later we can read from a config file or pass in from outside
         if (.not. allocated(this%FRLAI)) allocate(this%FRLAI(nx, ny, nlanduse))
         this%FRLAI(:,:,:) = 0.0_fp
         do j = 1, ny
            do i = 1, nx
               do k = 1, nlanduse
                  if (this%DLUSE(i, j) == k) this%FRLAI(i, j, k) = this%LAI(i, j) !TODO: should times fraclanduse but here is 1.0
               enddo
               this%FRLAI(i, j, 15:17) = 0.0 !manually give index 15(snow and ice), 16(barren), 17(water) zeros
            enddo
         enddo
       case ('SALINITY', 'salinity')
         this%SALINITY(:,:) = 0.0_fp  !set to zero for now, which will turn off O3 dry deposition over ocean with iodine.

       case ('REEVAPLS', 'reevapls')
         this%REEVAPLS(:,:,:) = 0.0_fp  !set to zero for now because I did not find data from GFS. This will overestimate the washout of aerosols.

       case default
         call error_mgr%report_error(ERROR_NOT_FOUND, &
            'Unknown derived field: ' // trim(field_name), rc, &
            thisLoc, 'Supported fields: AIRDEN,  TV,  BXHEIGHT')
         rc = CC_FAILURE
      end select

      call error_mgr%pop_context()
   end subroutine metstate_derive_field

END MODULE MetState_Mod
