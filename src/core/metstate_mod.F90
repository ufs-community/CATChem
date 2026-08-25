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
   USE iso_c_binding, only: c_ptr, c_null_ptr, c_associated



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
      type(c_ptr) :: cpp_ptr = c_null_ptr
      CHARACTER(LEN=3)             :: State     = 'MET'    !< Name of this state
      INTEGER                      :: NLEVS     = 127      !< Number of vertical levels (default)
      TYPE(GridGeometryType) :: geometry
      INTEGER                      :: NSURFTYPE = 20       !< Number of surface types (default)
      ! Grid flags (2D: nx, ny)
      LOGICAL, POINTER         :: IsLand(:,:) => null()       !< Is this a land grid box?
      LOGICAL, POINTER         :: IsWater(:,:) => null()      !< Is this a water grid box?
      LOGICAL, POINTER         :: IsIce(:,:) => null()        !< Is this an ice grid box?
      LOGICAL, POINTER         :: IsSnow(:,:) => null()       !< Is this a snow grid box?
      ! Vertical flags and arrays (3D: nx, ny, nz)
      LOGICAL, POINTER        :: InStratMeso(:,:,:) => null()    !< Are we in the stratosphere or mesosphere?
      LOGICAL, POINTER        :: InStratosphere(:,:,:) => null() !< Are we in the stratosphere?
      LOGICAL, POINTER        :: InTroposphere(:,:,:) => null()  !< Are we in the troposphere?
      LOGICAL, POINTER        :: InPbl(:,:,:) => null()          !< Are we in the PBL?
      LOGICAL, POINTER        :: IsLocalNoon(:,:) => null()      !< Is it local noon (between 11 and 13 local solar time)?
      ! Surface properties (2D: nx, ny)
      REAL(fp), POINTER        :: AREA_M2(:,:) => null()      !< Grid box surface area [m2]
      INTEGER, POINTER        :: LWI(:,:) => null()          !< Land water ice mask (0-sea, 1-land, 2-ice)
      INTEGER, POINTER        :: DLUSE(:,:) => null()        !< Dominant land-use type
      REAL(fp), POINTER        :: FRVEG(:,:) => null()        !< Fraction of veg [1]
      REAL(fp), POINTER        :: FRLAKE(:,:) => null()       !< Fraction of lake [1]
      REAL(fp), POINTER        :: FRLAND(:,:) => null()       !< Fraction of land [1]
      REAL(fp), POINTER        :: FRLANDIC(:,:) => null()     !< Fraction of land ice [1]
      REAL(fp), POINTER        :: FROCEAN(:,:) => null()      !< Fraction of ocean [1]
      REAL(fp), POINTER        :: FRSEAICE(:,:) => null()     !< Sfc sea ice fraction
      REAL(fp), POINTER        :: FRSNO(:,:) => null()        !< Sfc snow fraction
      REAL(fp), POINTER        :: LAI(:,:) => null()          !< Leaf area index [m2/m2] (online) Dominant
      REAL(fp), POINTER        :: GVF(:,:) => null()          !< Green Vegetative Fraction
      ! Dust Only Variables
      REAL(fp), POINTER        :: RDRAG(:,:) => null()        !< Drag Partition [1]
      REAL(fp), POINTER        :: USTAR_THRESHOLD(:,:) => null() !< Threshold friction velocity [m/s]
      REAL(fp), POINTER        :: SSM(:,:) => null()          !< Sediment Supply Map [1]
      ! Surface and ice properties (2D: nx, ny)
      REAL(fp), POINTER        :: SEAICE00(:,:) => null()     !< Sea ice coverage 00-10%
      REAL(fp), POINTER        :: SEAICE10(:,:) => null()     !< Sea ice coverage 10-20%
      REAL(fp), POINTER        :: SEAICE20(:,:) => null()     !< Sea ice coverage 20-30%
      REAL(fp), POINTER        :: SEAICE30(:,:) => null()     !< Sea ice coverage 30-40%
      REAL(fp), POINTER        :: SEAICE40(:,:) => null()     !< Sea ice coverage 40-50%
      REAL(fp), POINTER        :: SEAICE50(:,:) => null()     !< Sea ice coverage 50-60%
      REAL(fp), POINTER        :: SEAICE60(:,:) => null()     !< Sea ice coverage 60-70%
      REAL(fp), POINTER        :: SEAICE70(:,:) => null()     !< Sea ice coverage 70-80%
      REAL(fp), POINTER        :: SEAICE80(:,:) => null()     !< Sea ice coverage 80-90%
      REAL(fp), POINTER        :: SEAICE90(:,:) => null()     !< Sea ice coverage 90-100%
      REAL(fp), POINTER        :: SNODP(:,:) => null()        !< Snow depth [m]
      REAL(fp), POINTER        :: SNOMAS(:,:) => null()       !< Snow mass [kg/m2]

      ! Soil and land use arrays (2D for counts, 3D for fractions)
      INTEGER, POINTER        :: DSOILTYPE(:,:) => null()    !< Dominant soil type
      REAL(fp), POINTER        :: CLAYFRAC(:,:) => null()     !< Fraction of clay [1]
      REAL(fp), POINTER        :: SANDFRAC(:,:) => null()     !< Fraction of sand [1]
      INTEGER, POINTER        :: nLNDTYPE(:,:) => null()     !< # of landtypes in box (I,J)
      REAL(fp), POINTER        :: GWETTOP(:,:) => null()      !< Top soil moisture [1]
      REAL(fp), POINTER        :: GWETROOT(:,:) => null()     !< Root Zone soil moisture [1]
      REAL(fp), POINTER        :: WILT(:,:) => null()         !< Wilt point [1]
      INTEGER                      :: nSOIL             !< # number of soil layers
      INTEGER                      :: nSOILTYPE         !< # number of soil types
      REAL(fp), POINTER        :: SOILM(:,:,:) => null()      !< Volumetric Soil moisture [m3/m3] (nx,ny,nsoil)
      REAL(fp), POINTER        :: SOILT(:,:,:) => null()      !< Temperature of soil layer [K] (nx,ny,nsoil)
      REAL(fp), POINTER        :: FRLANDUSE(:,:,:) => null()  !< Fractional Land Use (nx,ny,nlanduse)
      REAL(fp), POINTER        :: FRSOIL(:,:,:) => null()     !< Fractional Soil (nx,ny,nsoil)
      REAL(fp), POINTER        :: FRLAI(:,:,:) => null()      !< LAI in each Fractional Land use type [m2/m2] (nx,ny,nlanduse)
      INTEGER, POINTER         :: ILAND(:,:,:) => null()      !< Land type ID in current grid box (nx,ny,nlanduse)
      ! Location arrays (1D for single point, 2D for grid)
      real(fp), POINTER        :: LAT(:,:) => null()         !< Latitude
      real(fp), POINTER        :: LON(:,:) => null()         !< Longitude
      character(len=20)            :: LUCNAME          !< name of land use category
      ! Surface meteorological properties (2D: nx, ny)
      REAL(fp), POINTER        :: ALBD_VIS(:,:) => null()     !< Visible surface albedo [1]
      REAL(fp), POINTER        :: ALBD_NIR(:,:) => null()     !< Near-IR surface albedo [1]
      REAL(fp), POINTER        :: ALBD_UV(:,:) => null()      !< UV surface albedo [1]
      REAL(fp), POINTER        :: PARDR(:,:) => null()        !< Direct photsynthetically active radiation [W/m2]
      REAL(fp), POINTER        :: PARDF(:,:) => null()        !< Diffuse photsynthetically active radiation [W/m2]
      REAL(fp), POINTER        :: SUNCOS(:,:) => null()       !< COS(solar zenith angle) at current time
      REAL(fp), POINTER        :: SUNCOSmid(:,:) => null()    !< COS(solar zenith angle) at midpoint of chem timestep
      REAL(fp), POINTER        :: SUNCOSsum(:,:) => null()    !< Sum of COS(SZA) for HEMCO OH diurnal variability
      REAL(fp), POINTER        :: SZAFACT(:,:) => null()      !< Diurnal scale factor for HEMCO OH diurnal variability (computed) [1]
      REAL(fp), POINTER        :: SWGDN(:,:) => null()        !< Incident radiation @ ground [W/m2]
      REAL(fp), POINTER        :: EFLUX(:,:) => null()        !< Latent heat flux [W/m2]
      REAL(fp), POINTER        :: HFLUX(:,:) => null()        !< Sensible heat flux [W/m2]
      REAL(fp), POINTER        :: HFLUX_UP(:,:) => null()     !< Sensible upward heat flux [W/m2]
      REAL(fp), POINTER        :: U10M(:,:) => null()         !< E/W wind speed @ 10m ht [m/s]
      REAL(fp), POINTER        :: USTAR(:,:) => null()        !< Friction velocity [m/s]
      REAL(fp), POINTER        :: V10M(:,:) => null()         !< N/S wind speed @ 10m ht [m/s]
      REAL(fp), POINTER        :: Z0(:,:) => null()           !< Surface roughness height [m]
      REAL(fp), POINTER        :: Z0H(:,:) => null()          !< Surface roughness height, for heat (thermal roughness) [m]
      REAL(fp), POINTER        :: FRZ0(:,:,:) => null()       !< Aerodynamic Roughness Length per FRLANDUSE (nx,ny,nlanduse)
      REAL(fp), POINTER        :: PBLH(:,:) => null()         !< PBL height [m]
      REAL(fp), POINTER        :: SALINITY(:,:) => null()     !< Salinity of the ocean [part per thousand]
      REAL(fp), POINTER        :: CMM(:,:) => null()          !< Aerodynamic conductance [m/s]
      REAL(fp), POINTER        :: ORO(:,:) => null()          !< surface height above sea level [m]
      REAL(fp), POINTER        :: RCA(:,:) => null()          !< Aerodynamic resistance in canopy [s/m]
      REAL(fp), ALLOCATABLE        :: WCA(:,:)          ! canopy water amount [kg/m2]
      ! 3D volumetric fields (3D: nx, ny, nz)
      REAL(fp), POINTER        :: F_OF_PBL(:,:,:) => null()       !< Fraction of box within PBL [1]
      REAL(fp), POINTER        :: F_UNDER_PBLTOP(:,:,:) => null() !< Fraction of box under PBL top
      real(fp), POINTER        :: OBK(:,:) => null()          !< Monin-Obhukov length [m]
      ! Cloud and precipitation properties (2D for surface, 3D for volumetric)
      REAL(fp), POINTER        :: CLDFRC(:,:) => null()       !< Column cloud fraction [1]
      REAL(fp), POINTER        :: CONV_DEPTH(:,:) => null()   !< Convective cloud depth [m]
      REAL(fp), POINTER        :: FLASH_DENS(:,:) => null()   !< Lightning flash density [#/km2/s]
      REAL(fp), POINTER        :: CNV_FRC(:,:) => null()      !< Convective fraction [1]
      REAL(fp), POINTER        :: CLDF(:,:,:) => null()       !< 3-D cloud fraction [1]
      REAL(fp), POINTER        :: CMFMC(:,:,:) => null()      !< Cloud mass flux [kg/m2/s]
      REAL(fp), POINTER        :: DQRCU(:,:,:) => null()      !< Conv precip production rate [kg/kg/s] (assume per dry air)
      REAL(fp), POINTER        :: DQRLSAN(:,:,:) => null()    !< LS precip prod rate [kg/kg/s] (assume per dry air)
      REAL(fp), POINTER        :: DTRAIN(:,:,:) => null()     !< Detrainment flux [kg/m2/s]
      REAL(fp), POINTER        :: PRECANV(:,:) => null()      !< Anvil previp @ ground [kg/m2/s] -> [mm/day]
      REAL(fp), POINTER        :: PRECCON(:,:) => null()      !< Conv  precip @ ground [kg/m2/s] -> [mm/day]
      REAL(fp), POINTER        :: PRECLSC(:,:) => null()      !< Large-scale precip @ ground kg/m2/s] -> [mm/day]
      real(fp), POINTER        :: REEVAPLS(:,:,:) => null()   !< Evap of precip LS+anvil [kg/kg/s] (assume per dry air)
      ! 3D cloud and precipitation arrays
      REAL(fp), POINTER        :: QI(:,:,:) => null()         !< Mass fraction of cloud ice water [kg/kg dry air]
      REAL(fp), POINTER        :: QL(:,:,:) => null()         !< Mass fraction of cloud liquid water [kg/kg dry air]
      REAL(fp), POINTER        :: PFICU(:,:,:) => null()      !< Dwn flux ice prec:conv [kg/m2/s]
      REAL(fp), POINTER        :: PFILSAN(:,:,:) => null()    !< Dwn flux ice prec:LS+anv [kg/m2/s] (nx,ny,nz+1)
      REAL(fp), POINTER        :: PFLCU(:,:,:) => null()      !< Dwn flux liq prec:conv [kg/m2/s]
      REAL(fp), POINTER        :: PFLLSAN(:,:,:) => null()    !< Dwn flux liq prec:LS+anv [kg/m2/s] (nx,ny,nz+1)
      REAL(fp), POINTER        :: TAUCLI(:,:,:) => null()     !< Opt depth of ice clouds [1]
      REAL(fp), POINTER        :: TAUCLW(:,:,:) => null()     !< Opt depth of H2O clouds [1]
      ! Surface scalars (now 2D: nx, ny)
      REAL(fp), POINTER        :: PHIS(:,:) => null()         !< Surface geopotential height [m2/s2]
      REAL(fp), POINTER        :: PS_WET(:,:) => null()       !< Wet surface pressure at start of timestep [Pa]
      REAL(fp), POINTER        :: PS_DRY(:,:) => null()       !< Dry surface pressure at start of timestep [Pa]
      REAL(fp), POINTER        :: QV2M(:,:) => null()         !< Specific Humidity at 2m [kg/kg]
      REAL(fp), POINTER        :: T2M(:,:) => null()          !< Temperature 2m [K]
      REAL(fp), POINTER        :: TS(:,:) => null()           !< Surface temperature [K]
      REAL(fp), POINTER        :: TSKIN(:,:) => null()        !< Surface skin temperature [K]
      REAL(fp), POINTER        :: SST(:,:) => null()          !< Sea surface temperature [K]
      REAL(fp), POINTER        :: SLP(:,:) => null()          !< Sea level pressure [Pa]
      REAL(fp), POINTER        :: PS(:,:) => null()           !< Surface Pressure [Pa]
      REAL(fp), POINTER        :: TO3(:,:) => null()          !< Total overhead O3 column [DU]
      REAL(fp), POINTER        :: TROPP(:,:) => null()        !< Tropopause pressure [Pa]
      INTEGER, POINTER        :: TropLev(:,:) => null()      !< Tropopause level [1]
      REAL(fp), POINTER        :: TropHt(:,:) => null()       !< Tropopause height [km]
      ! 3D atmospheric variables (3D: nx, ny, nz)
      REAL(fp), POINTER        :: Z(:,:,:) => null()          !< Geopotential Height @ level edges [m] (nx,ny,nz+1)
      REAL(fp), POINTER        :: ZMID(:,:,:) => null()       !< Mid Layer Geopotential Height [m]
      REAL(fp), POINTER        :: BXHEIGHT(:,:,:) => null()   !< Grid box height [m] (dry air)
      REAL(fp), POINTER        :: QV(:,:,:) => null()         !< Specific Humidity [kg/kg]
      REAL(fp), POINTER        :: T(:,:,:) => null()          !< Temperature [K]
      REAL(fp), POINTER        :: THETA(:,:,:) => null()      !< Potential temperature [K]
      REAL(fp), POINTER        :: TV(:,:,:) => null()         !< Virtual temperature [K]
      REAL(fp), POINTER        :: V(:,:,:) => null()          !< N/S component of wind [m s-1]
      REAL(fp), POINTER        :: U(:,:,:) => null()          !< E/W component of wind [m s-1]
      REAL(fp), POINTER        :: OMEGA(:,:,:) => null()      !< Updraft velocity [Pa/s]
      REAL(fp), POINTER        :: RH(:,:,:) => null()         !< Relative humidity [fraction, not %]
      REAL(fp), POINTER        :: SPHU(:,:,:) => null()       !< Specific humidity [g H2O/kg tot air]
      REAL(fp), POINTER        :: AIRDEN(:,:,:) => null()     !< Wet air density [kg/m3]
      REAL(fp), POINTER        :: AIRDEN_DRY(:,:,:) => null() !< Dry air density [kg/m3]
      REAL(fp), POINTER        :: AIRNUMDEN(:,:,:) => null()  !< Dry air density [molec/cm3]
      REAL(fp), POINTER        :: MAIRDEN(:,:,:) => null()    !< Moist air density (same as AIRDEN to cover possible use cases) [kg/m3]
      REAL(fp), POINTER        :: AVGW(:,:,:) => null()       !< Water vapor volume mixing ratio [vol H2O/vol dry air]
      REAL(fp), POINTER        :: DELP(:,:,:) => null()       !< Delta-P (wet) across box [Pa]
      REAL(fp), POINTER        :: DELP_DRY(:,:,:) => null()   !< Delta-P (dry) across box [Pa]
      REAL(fp), POINTER        :: DAIRMASS(:,:,:) => null()   !< Dry air mass [kg] in grid box
      REAL(fp), POINTER        :: AIRVOL(:,:,:) => null()     !< Grid box volume [m3] (dry air)
      REAL(fp), POINTER        :: PEDGE_DRY(:,:,:) => null()  !< Dry air partial pressure @ level edges [Pa] (nx,ny,nz+1)
      REAL(fp), POINTER        :: PEDGE(:,:,:) => null()      !< Air partial pressure @ level edges [Pa] (nx,ny,nz+1)
      REAL(fp), POINTER        :: PMID(:,:,:) => null()       !< Average wet air pressure [Pa] defined as arithmetic average of edge pressures
      REAL(fp), POINTER        :: PMID_DRY(:,:,:) => null()   !< Dry air partial pressure [Pa] defined as arithmetic avg of edge pressures
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
      use error_mod, only: ErrorManagerType, CC_SUCCESS

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
      use error_mod, only: ErrorManagerType, CC_SUCCESS

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

      nsoil = this%nSOIL
      nsoiltype = this%nSOILTYPE
      nSURFTYPE = this%NSURFTYPE

      select case (to_upper(trim(adjustl(field_name))))
       case ('ALL', 'all')
         ! Allocate core meteorological pointer arrays only if not already associated (prevents leaks/orphaning)
         if (.not. associated(this%T)) allocate(this%T(nx, ny, nz))
         if (.not. associated(this%QV)) allocate(this%QV(nx, ny, nz))
         if (.not. associated(this%RH)) allocate(this%RH(nx, ny, nz))
         if (.not. associated(this%PMID)) allocate(this%PMID(nx, ny, nz))
         if (.not. associated(this%PEDGE)) allocate(this%PEDGE(nx, ny, nz+1))
         if (.not. associated(this%AIRDEN)) allocate(this%AIRDEN(nx, ny, nz))
         if (.not. associated(this%AIRDEN_DRY)) allocate(this%AIRDEN_DRY(nx, ny, nz))
         if (.not. associated(this%BXHEIGHT)) allocate(this%BXHEIGHT(nx, ny, nz))
         if (.not. associated(this%DELP)) allocate(this%DELP(nx, ny, nz))
         if (.not. associated(this%DELP_DRY)) allocate(this%DELP_DRY(nx, ny, nz))
         if (.not. associated(this%PS)) allocate(this%PS(nx, ny))
         if (.not. associated(this%TS)) allocate(this%TS(nx, ny))
         if (.not. associated(this%PBLH)) allocate(this%PBLH(nx, ny))
         if (.not. associated(this%USTAR)) allocate(this%USTAR(nx, ny))
         if (.not. associated(this%HFLUX)) allocate(this%HFLUX(nx, ny))
         if (.not. associated(this%OBK)) allocate(this%OBK(nx, ny))
         if (.not. associated(this%LAT)) allocate(this%LAT(nx, ny))
         if (.not. associated(this%LON)) allocate(this%LON(nx, ny))
         if (.not. associated(this%IsLand)) allocate(this%IsLand(nx, ny))
         if (.not. associated(this%IsWater)) allocate(this%IsWater(nx, ny))
         if (.not. associated(this%IsIce)) allocate(this%IsIce(nx, ny))
         if (.not. associated(this%IsSnow)) allocate(this%IsSnow(nx, ny))
         if (.not. associated(this%LWI)) allocate(this%LWI(nx, ny))
         if (.not. associated(this%DLUSE)) allocate(this%DLUSE(nx, ny))
         if (.not. associated(this%FRVEG)) allocate(this%FRVEG(nx, ny))
         if (.not. associated(this%AREA_M2)) allocate(this%AREA_M2(nx, ny))
         if (.not. associated(this%U)) allocate(this%U(nx, ny, nz))
         if (.not. associated(this%V)) allocate(this%V(nx, ny, nz))

       case ('T', 't')
         allocate(this%T(nx, ny, nz))
       case ('QV', 'qv')
         allocate(this%QV(nx, ny, nz))
       case ('RH', 'rh')
         allocate(this%RH(nx, ny, nz))
       case ('PMID', 'pmid')
         allocate(this%PMID(nx, ny, nz))
       case ('PEDGE', 'pedge')
         allocate(this%PEDGE(nx, ny, nz+1))
       case ('AIRDEN', 'airden')
         allocate(this%AIRDEN(nx, ny, nz))
       case ('AIRDEN_DRY', 'airden_dry')
         allocate(this%AIRDEN_DRY(nx, ny, nz))
       case ('BXHEIGHT', 'bxheight')
         allocate(this%BXHEIGHT(nx, ny, nz))
       case ('DELP', 'delp')
         allocate(this%DELP(nx, ny, nz))
       case ('DELP_DRY', 'delp_dry')
         allocate(this%DELP_DRY(nx, ny, nz))
       case ('PS', 'ps')
         allocate(this%PS(nx, ny))
       case ('TS', 'ts')
         allocate(this%TS(nx, ny))
       case ('PBLH', 'pblh')
         allocate(this%PBLH(nx, ny))
       case ('USTAR', 'ustar')
         allocate(this%USTAR(nx, ny))
       case ('HFLUX', 'hflux')
         allocate(this%HFLUX(nx, ny))
       case ('OBK', 'obk')
         allocate(this%OBK(nx, ny))
       case ('LAT', 'lat')
         allocate(this%LAT(nx, ny))
       case ('LON', 'lon')
         allocate(this%LON(nx, ny))
       case ('U', 'u')
         allocate(this%U(nx, ny, nz))
       case ('V', 'v')
         allocate(this%V(nx, ny, nz))
       case ('IsLand', 'island')
         allocate(this%IsLand(nx, ny))
       case ('IsWater', 'iswater')
         allocate(this%IsWater(nx, ny))
       case ('LWI', 'lwi')
         allocate(this%LWI(nx, ny))
       case ('FRVEG', 'frveg')
         allocate(this%FRVEG(nx, ny))
       case ('AREA_M2', 'area_m2')
         allocate(this%AREA_M2(nx, ny))
       case ('InPbl', 'inpbl')
         allocate(this%InPbl(nx, ny, nz))
       case ('TropLev', 'troplev')
         allocate(this%TropLev(nx, ny))
      end select

      ! Initialize to safe defaults if associated
      if (associated(this%T)) this%T = 288.15_fp
      if (associated(this%U)) this%U = 0.0_fp
      if (associated(this%V)) this%V = 0.0_fp
      if (associated(this%QV)) this%QV = 0.001_fp
      if (associated(this%RH)) this%RH = 0.50_fp
      if (associated(this%AIRDEN)) this%AIRDEN = 1.2_fp
      if (associated(this%BXHEIGHT)) this%BXHEIGHT = 100.0_fp
      if (associated(this%PS)) this%PS = 101300.25_fp
      if (associated(this%SLP)) this%SLP = 101300.25_fp
      if (associated(this%T2M)) this%T2M = 288.15_fp
      if (associated(this%TS)) this%TS = 288.15_fp

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

      select case (to_upper(trim(adjustl(field_name))))
       case ('ALL', 'all')
         if (associated(this%T)) deallocate(this%T)
         if (associated(this%QV)) deallocate(this%QV)
         if (associated(this%RH)) deallocate(this%RH)
         if (associated(this%PMID)) deallocate(this%PMID)
         if (associated(this%PEDGE)) deallocate(this%PEDGE)
         if (associated(this%AIRDEN)) deallocate(this%AIRDEN)
         if (associated(this%AIRDEN_DRY)) deallocate(this%AIRDEN_DRY)
         if (associated(this%BXHEIGHT)) deallocate(this%BXHEIGHT)
         if (associated(this%DELP)) deallocate(this%DELP)
         if (associated(this%DELP_DRY)) deallocate(this%DELP_DRY)
         if (associated(this%PS)) deallocate(this%PS)
         if (associated(this%TS)) deallocate(this%TS)
         if (associated(this%PBLH)) deallocate(this%PBLH)
         if (associated(this%USTAR)) deallocate(this%USTAR)
         if (associated(this%HFLUX)) deallocate(this%HFLUX)
         if (associated(this%OBK)) deallocate(this%OBK)
         if (associated(this%LAT)) deallocate(this%LAT)
         if (associated(this%LON)) deallocate(this%LON)
         if (associated(this%U)) deallocate(this%U)
         if (associated(this%V)) deallocate(this%V)
         if (associated(this%IsLand)) deallocate(this%IsLand)
         if (associated(this%IsWater)) deallocate(this%IsWater)
         if (associated(this%LWI)) deallocate(this%LWI)
         if (associated(this%FRVEG)) deallocate(this%FRVEG)
         if (associated(this%AREA_M2)) deallocate(this%AREA_M2)
         if (associated(this%InPbl)) deallocate(this%InPbl)
         if (associated(this%TropLev)) deallocate(this%TropLev)

         ! Nullify pointers
         nullify(this%T, this%QV, this%RH, this%PMID, this%PEDGE, this%AIRDEN, this%AIRDEN_DRY)
         nullify(this%BXHEIGHT, this%DELP, this%DELP_DRY, this%PS, this%TS, this%PBLH)
         nullify(this%USTAR, this%HFLUX, this%OBK, this%LAT, this%LON, this%U, this%V)
         nullify(this%IsLand, this%IsWater, this%LWI, this%FRVEG, this%AREA_M2, this%InPbl, this%TropLev)

       case ('T', 't')
         if (associated(this%T)) deallocate(this%T)
         nullify(this%T)
       case ('QV', 'qv')
         if (associated(this%QV)) deallocate(this%QV)
         nullify(this%QV)
       case ('RH', 'rh')
         if (associated(this%RH)) deallocate(this%RH)
         nullify(this%RH)
       case ('PMID', 'pmid')
         if (associated(this%PMID)) deallocate(this%PMID)
         nullify(this%PMID)
       case ('PEDGE', 'pedge')
         if (associated(this%PEDGE)) deallocate(this%PEDGE)
         nullify(this%PEDGE)
       case ('AIRDEN', 'airden')
         if (associated(this%AIRDEN)) deallocate(this%AIRDEN)
         nullify(this%AIRDEN)
       case ('AIRDEN_DRY', 'airden_dry')
         if (associated(this%AIRDEN_DRY)) deallocate(this%AIRDEN_DRY)
         nullify(this%AIRDEN_DRY)
       case ('BXHEIGHT', 'bxheight')
         if (associated(this%BXHEIGHT)) deallocate(this%BXHEIGHT)
         nullify(this%BXHEIGHT)
       case ('DELP', 'delp')
         if (associated(this%DELP)) deallocate(this%DELP)
         nullify(this%DELP)
       case ('DELP_DRY', 'delp_dry')
         if (associated(this%DELP_DRY)) deallocate(this%DELP_DRY)
         nullify(this%DELP_DRY)
       case ('PS', 'ps')
         if (associated(this%PS)) deallocate(this%PS)
         nullify(this%PS)
       case ('TS', 'ts')
         if (associated(this%TS)) deallocate(this%TS)
         nullify(this%TS)
       case ('PBLH', 'pblh')
         if (associated(this%PBLH)) deallocate(this%PBLH)
         nullify(this%PBLH)
       case ('USTAR', 'ustar')
         if (associated(this%USTAR)) deallocate(this%USTAR)
         nullify(this%USTAR)
       case ('HFLUX', 'hflux')
         if (associated(this%HFLUX)) deallocate(this%HFLUX)
         nullify(this%HFLUX)
       case ('OBK', 'obk')
         if (associated(this%OBK)) deallocate(this%OBK)
         nullify(this%OBK)
       case ('LAT', 'lat')
         if (associated(this%LAT)) deallocate(this%LAT)
         nullify(this%LAT)
       case ('LON', 'lon')
         if (associated(this%LON)) deallocate(this%LON)
         nullify(this%LON)
       case ('U', 'u')
         if (associated(this%U)) deallocate(this%U)
         nullify(this%U)
       case ('V', 'v')
         if (associated(this%V)) deallocate(this%V)
         nullify(this%V)
       case ('IsLand', 'island')
         if (associated(this%IsLand)) deallocate(this%IsLand)
         nullify(this%IsLand)
       case ('IsWater', 'iswater')
         if (associated(this%IsWater)) deallocate(this%IsWater)
         nullify(this%IsWater)
       case ('LWI', 'lwi')
         if (associated(this%LWI)) deallocate(this%LWI)
         nullify(this%LWI)
       case ('FRVEG', 'frveg')
         if (associated(this%FRVEG)) deallocate(this%FRVEG)
         nullify(this%FRVEG)
       case ('AREA_M2', 'area_m2')
         if (associated(this%AREA_M2)) deallocate(this%AREA_M2)
         nullify(this%AREA_M2)
       case ('InPbl', 'inpbl')
         if (associated(this%InPbl)) deallocate(this%InPbl)
         nullify(this%InPbl)
       case ('TropLev', 'troplev')
         if (associated(this%TropLev)) deallocate(this%TropLev)
         nullify(this%TropLev)
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
      if (associated(this%T2M)) then
         if (maxval(this%T2M) > 400.0_fp .or. minval(this%T2M) < 100.0_fp) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               '2m temperature out of physical range', rc, &
               thisLoc, 'Check temperature units and values')
            call error_mgr%pop_context()
            return
         endif
      endif

      if (associated(this%TS)) then
         if (maxval(this%TS) > 400.0_fp .or. minval(this%TS) < 100.0_fp) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'Surface temperature out of physical range', rc, &
               thisLoc, 'Check temperature units and values')
            call error_mgr%pop_context()
            return
         endif
      endif

      ! Validate pressures (use maxval/minval for array validation)
      if (associated(this%PS)) then
         if (maxval(this%PS) > 120000.0_fp .or. minval(this%PS) < 1000.0_fp) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'Surface pressure out of physical range', rc, &
               thisLoc, 'Check pressure units and values')
            call error_mgr%pop_context()
            return
         endif
      endif

      if (associated(this%SLP)) then
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
      if (associated(this%T2M)) this%T2M = 288.15_fp
      if (associated(this%TS)) this%TS = 288.15_fp
      if (associated(this%TSKIN)) this%TSKIN = 288.15_fp
      if (associated(this%PS)) this%PS = 101300.25_fp
      if (associated(this%SLP)) this%SLP = 101300.25_fp
      if (associated(this%SST)) this%SST = 288.15_fp

      ! Reset arrays if allocated
      if (associated(this%T)) this%T = 288.15_fp
      if (associated(this%U)) this%U = 0.0_fp
      if (associated(this%V)) this%V = 0.0_fp
      if (associated(this%QV)) this%QV = 0.01_fp
      if (associated(this%RH)) this%RH = 0.5_fp

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

      is_alloc = associated(this%T) .and. associated(this%U) .and. associated(this%V) .and. &
         associated(this%QV) .and. associated(this%PMID) .and. associated(this%DELP)
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
      if (associated(this%TS)) then
         write(*,'(A,F8.2,A)') 'Surface temperature: ', this%TS(1,1), ' K'
      endif
      if (associated(this%PS)) then
         write(*,'(A,F8.2,A)') 'Surface pressure: ', this%PS(1,1), ' Pa'
      endif
      if (associated(this%SLP)) then
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
      nsoil = this%nSOIL
      nsoiltype = this%nSOILTYPE
      nSURFTYPE = this%NSURFTYPE
      select case (to_upper(trim(field_name)))
       case ('T')
         if (.not. associated(this%T)) allocate(this%T(nx, ny, nz))
       case ('U')
         if (.not. associated(this%U)) allocate(this%U(nx, ny, nz))
       case ('V')
         if (.not. associated(this%V)) allocate(this%V(nx, ny, nz))
       case ('QV')
         if (.not. associated(this%QV)) allocate(this%QV(nx, ny, nz))
       case ('RH')
         if (.not. associated(this%RH)) allocate(this%RH(nx, ny, nz))
       case ('PMID')
         if (.not. associated(this%PMID)) allocate(this%PMID(nx, ny, nz))
       case ('PEDGE')
         if (.not. associated(this%PEDGE)) allocate(this%PEDGE(nx, ny, nz+1))
       case ('AIRDEN')
         if (.not. associated(this%AIRDEN)) allocate(this%AIRDEN(nx, ny, nz))
       case ('AIRDEN_DRY')
         if (.not. associated(this%AIRDEN_DRY)) allocate(this%AIRDEN_DRY(nx, ny, nz))
       case ('BXHEIGHT')
         if (.not. associated(this%BXHEIGHT)) allocate(this%BXHEIGHT(nx, ny, nz))
       case ('DELP')
         if (.not. associated(this%DELP)) allocate(this%DELP(nx, ny, nz))
       case ('DELP_DRY')
         if (.not. associated(this%DELP_DRY)) allocate(this%DELP_DRY(nx, ny, nz))
       case ('PS')
         if (.not. associated(this%PS)) allocate(this%PS(nx, ny))
       case ('TS')
         if (.not. associated(this%TS)) allocate(this%TS(nx, ny))
       case ('PBLH')
         if (.not. associated(this%PBLH)) allocate(this%PBLH(nx, ny))
       case ('USTAR')
         if (.not. associated(this%USTAR)) allocate(this%USTAR(nx, ny))
       case ('HFLUX')
         if (.not. associated(this%HFLUX)) allocate(this%HFLUX(nx, ny))
       case ('OBK')
         if (.not. associated(this%OBK)) allocate(this%OBK(nx, ny))
       case ('LAT')
         if (.not. associated(this%LAT)) allocate(this%LAT(nx, ny))
       case ('LON')
         if (.not. associated(this%LON)) allocate(this%LON(nx, ny))
       case ('FRVEG')
         if (.not. associated(this%FRVEG)) allocate(this%FRVEG(nx, ny))
       case ('AREA_M2')
         if (.not. associated(this%AREA_M2)) allocate(this%AREA_M2(nx, ny))
       case ('LWI')
         if (.not. associated(this%LWI)) allocate(this%LWI(nx, ny))
       case ('DLUSE')
         if (.not. associated(this%DLUSE)) allocate(this%DLUSE(nx, ny))
       case ('DSOILTYPE')
         if (.not. associated(this%DSOILTYPE)) allocate(this%DSOILTYPE(nx, ny))
       case ('NLNDTYPE')
         if (.not. associated(this%nLNDTYPE)) allocate(this%nLNDTYPE(nx, ny))
       case ('TROPLEV')
         if (.not. associated(this%TropLev)) allocate(this%TropLev(nx, ny))
       case ('ISLAND')
         if (.not. associated(this%IsLand)) allocate(this%IsLand(nx, ny))
       case ('ISWATER')
         if (.not. associated(this%IsWater)) allocate(this%IsWater(nx, ny))
       case ('ISICE')
         if (.not. associated(this%IsIce)) allocate(this%IsIce(nx, ny))
       case ('ISSNOW')
         if (.not. associated(this%IsSnow)) allocate(this%IsSnow(nx, ny))
       case ('ISLOCALNOON')
         if (.not. associated(this%IsLocalNoon)) allocate(this%IsLocalNoon(nx, ny))
       case ('INPBL')
         if (.not. associated(this%InPbl)) allocate(this%InPbl(nx, ny, nz))
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
      select case (to_upper(trim(field_name)))
       case ('T')
         if (associated(this%T)) deallocate(this%T)
       case ('U')
         if (associated(this%U)) deallocate(this%U)
       case ('V')
         if (associated(this%V)) deallocate(this%V)
       case ('QV')
         if (associated(this%QV)) deallocate(this%QV)
       case ('RH')
         if (associated(this%RH)) deallocate(this%RH)
       case ('PMID')
         if (associated(this%PMID)) deallocate(this%PMID)
       case ('PEDGE')
         if (associated(this%PEDGE)) deallocate(this%PEDGE)
       case ('AIRDEN')
         if (associated(this%AIRDEN)) deallocate(this%AIRDEN)
       case ('AIRDEN_DRY')
         if (associated(this%AIRDEN_DRY)) deallocate(this%AIRDEN_DRY)
       case ('BXHEIGHT')
         if (associated(this%BXHEIGHT)) deallocate(this%BXHEIGHT)
       case ('DELP')
         if (associated(this%DELP)) deallocate(this%DELP)
       case ('DELP_DRY')
         if (associated(this%DELP_DRY)) deallocate(this%DELP_DRY)
       case ('PS')
         if (associated(this%PS)) deallocate(this%PS)
       case ('TS')
         if (associated(this%TS)) deallocate(this%TS)
       case ('PBLH')
         if (associated(this%PBLH)) deallocate(this%PBLH)
       case ('USTAR')
         if (associated(this%USTAR)) deallocate(this%USTAR)
       case ('HFLUX')
         if (associated(this%HFLUX)) deallocate(this%HFLUX)
       case ('OBK')
         if (associated(this%OBK)) deallocate(this%OBK)
       case ('LAT')
         if (associated(this%LAT)) deallocate(this%LAT)
       case ('LON')
         if (associated(this%LON)) deallocate(this%LON)
       case ('FRVEG')
         if (associated(this%FRVEG)) deallocate(this%FRVEG)
       case ('AREA_M2')
         if (associated(this%AREA_M2)) deallocate(this%AREA_M2)
       case ('LWI')
         if (associated(this%LWI)) deallocate(this%LWI)
       case ('DLUSE')
         if (associated(this%DLUSE)) deallocate(this%DLUSE)
       case ('DSOILTYPE')
         if (associated(this%DSOILTYPE)) deallocate(this%DSOILTYPE)
       case ('NLNDTYPE')
         if (associated(this%nLNDTYPE)) deallocate(this%nLNDTYPE)
       case ('TROPLEV')
         if (associated(this%TropLev)) deallocate(this%TropLev)
       case ('ISLAND')
         if (associated(this%IsLand)) deallocate(this%IsLand)
       case ('ISWATER')
         if (associated(this%IsWater)) deallocate(this%IsWater)
       case ('ISICE')
         if (associated(this%IsIce)) deallocate(this%IsIce)
       case ('ISSNOW')
         if (associated(this%IsSnow)) deallocate(this%IsSnow)
       case ('ISLOCALNOON')
         if (associated(this%IsLocalNoon)) deallocate(this%IsLocalNoon)
       case ('INPBL')
         if (associated(this%InPbl)) deallocate(this%InPbl)
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
      use Interop_Mod, only: get_cpp_field
      implicit none
      class(MetStateType), intent(in), target :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in), optional :: i, j
      real(fp), pointer :: column_ptr(:)
      integer :: nx, ny, nlev, col_i, col_j, rc
      real(fp), pointer :: full_3d(:,:,:)
      column_ptr => null()
      full_3d => null()
      call this%get_dimensions(nx, ny, nlev)
      col_i = 1; col_j = 1
      if (present(i)) col_i = max(1, min(i, nx))
      if (present(j)) col_j = max(1, min(j, ny))

      if (c_associated(this%cpp_ptr)) then
         call get_cpp_field(this%cpp_ptr, field_name, full_3d, [nx, ny, nlev], rc)
         if (rc == 0 .and. associated(full_3d)) then
            column_ptr => full_3d(col_i, col_j, :)
            return
         end if
      end if

      select case (to_upper(trim(field_name)))
       case ('T')
         if (associated(this%T)) full_3d => this%T
       case ('QV')
         if (associated(this%QV)) full_3d => this%QV
       case ('RH')
         if (associated(this%RH)) full_3d => this%RH
       case ('PMID')
         if (associated(this%PMID)) full_3d => this%PMID
       case ('PEDGE')
         if (associated(this%PEDGE)) full_3d => this%PEDGE
       case ('AIRDEN')
         if (associated(this%AIRDEN)) full_3d => this%AIRDEN
       case ('AIRDEN_DRY')
         if (associated(this%AIRDEN_DRY)) full_3d => this%AIRDEN_DRY
       case ('BXHEIGHT')
         if (associated(this%BXHEIGHT)) full_3d => this%BXHEIGHT
       case ('DELP')
         if (associated(this%DELP)) full_3d => this%DELP
       case ('DELP_DRY')
         if (associated(this%DELP_DRY)) full_3d => this%DELP_DRY
       case ('U')
         if (associated(this%U)) full_3d => this%U
       case ('V')
         if (associated(this%V)) full_3d => this%V
       case ('PFILSAN')
         if (associated(this%PFILSAN)) full_3d => this%PFILSAN
       case ('PFLLSAN')
         if (associated(this%PFLLSAN)) full_3d => this%PFLLSAN
      end select

      if (associated(full_3d)) then
         column_ptr => full_3d(col_i, col_j, :)
      endif
   end function metstate_get_column_ptr_func

   !> Get a scalar value from a 2D field at (i,j)
   function metstate_get_2Dto0D_value(this, field_name, i, j) result(scalar_val)
      use Interop_Mod, only: get_cpp_field
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: i, j
      real(fp) :: scalar_val
      integer :: col_i, col_j, rc, nx, ny, nz
      real(fp), pointer :: full_2d(:,:)
      col_i = i
      col_j = j
      scalar_val = 0.0_fp

      if (c_associated(this%cpp_ptr)) then
         call this%geometry%get_dimensions(nx, ny, nz)
         call get_cpp_field(this%cpp_ptr, field_name, full_2d, [nx, ny], rc)
         if (rc == 0 .and. associated(full_2d)) then
            scalar_val = full_2d(col_i, col_j)
            return
         end if
      end if

      select case (to_upper(trim(field_name)))
       case ('PS')
         if (associated(this%PS)) scalar_val = this%PS(col_i, col_j)
       case ('TS')
         if (associated(this%TS)) scalar_val = this%TS(col_i, col_j)
       case ('PBLH')
         if (associated(this%PBLH)) scalar_val = this%PBLH(col_i, col_j)
       case ('USTAR')
         if (associated(this%USTAR)) scalar_val = this%USTAR(col_i, col_j)
       case ('HFLUX')
         if (associated(this%HFLUX)) scalar_val = this%HFLUX(col_i, col_j)
       case ('OBK')
         if (associated(this%OBK)) scalar_val = this%OBK(col_i, col_j)
       case ('LAT')
         if (associated(this%LAT)) scalar_val = this%LAT(col_i, col_j)
       case ('LON')
         if (associated(this%LON)) scalar_val = this%LON(col_i, col_j)
       case ('FRVEG')
         if (associated(this%FRVEG)) scalar_val = this%FRVEG(col_i, col_j)
       case ('AREA_M2')
         if (associated(this%AREA_M2)) scalar_val = this%AREA_M2(col_i, col_j)
       case ('FROCEAN')
         if (associated(this%FROCEAN)) scalar_val = this%FROCEAN(col_i, col_j)
       case ('FRSEAICE')
         if (associated(this%FRSEAICE)) scalar_val = this%FRSEAICE(col_i, col_j)
       case ('SST')
         if (associated(this%SST)) scalar_val = this%SST(col_i, col_j)
      end select
   end function metstate_get_2Dto0D_value

   !> Get a scalar value from a scalar field
   function metstate_get_scalar_value(this, field_name) result(scalar_val)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      real(fp) :: scalar_val
      scalar_val = 0.0_fp
   end function metstate_get_scalar_value

   !> INTEGER versions of accessor functions
   function metstate_get_column_ptr_func_int(this, field_name, i, j) result(column_ptr)
      implicit none
      class(MetStateType), intent(in), target :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in), optional :: i, j
      integer, pointer :: column_ptr(:)
      integer :: nx, ny, nlev, col_i, col_j
      integer, pointer :: full_3d(:,:,:)
      column_ptr => null()
      full_3d => null()
      call this%get_dimensions(nx, ny, nlev)
      col_i = 1; col_j = 1
      if (present(i)) col_i = max(1, min(i, nx))
      if (present(j)) col_j = max(1, min(j, ny))

      select case (to_upper(trim(field_name)))
       case ('ILAND')
         if (associated(this%ILAND)) full_3d => this%ILAND
      end select

      if (associated(full_3d)) then
         column_ptr => full_3d(col_i, col_j, :)
      endif
   end function metstate_get_column_ptr_func_int

   function metstate_get_2Dto0D_value_int(this, field_name, i, j) result(scalar_val)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: i, j
      integer :: scalar_val
      integer :: col_i, col_j
      col_i = i
      col_j = j
      scalar_val = 0
      select case (to_upper(trim(field_name)))
       case ('LWI')
         if (associated(this%LWI)) scalar_val = this%LWI(col_i, col_j)
       case ('DLUSE')
         if (associated(this%DLUSE)) scalar_val = this%DLUSE(col_i, col_j)
       case ('DSOILTYPE')
         if (associated(this%DSOILTYPE)) scalar_val = this%DSOILTYPE(col_i, col_j)
       case ('NLNDTYPE')
         if (associated(this%nLNDTYPE)) scalar_val = this%nLNDTYPE(col_i, col_j)
       case ('TROPLEV')
         if (associated(this%TropLev)) scalar_val = this%TropLev(col_i, col_j)
      end select
   end function metstate_get_2Dto0D_value_int

   function metstate_get_scalar_value_int(this, field_name) result(scalar_val)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      integer :: scalar_val
      scalar_val = 0
      select case (to_upper(trim(field_name)))
       case ('NLEVS')
         scalar_val = this%NLEVS
       case ('NSOIL')
         scalar_val = this%nSOIL
       case ('NSOILTYPE')
         scalar_val = this%nSOILTYPE
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
      logical, pointer :: full_3d(:,:,:)
      column_ptr => null()
      full_3d => null()
      call this%get_dimensions(nx, ny, nlev)
      col_i = 1; col_j = 1
      if (present(i)) col_i = max(1, min(i, nx))
      if (present(j)) col_j = max(1, min(j, ny))

      select case (to_upper(trim(field_name)))
       case ('INSTRATMESO')
         if (associated(this%InStratMeso)) full_3d => this%InStratMeso
       case ('INSTRATOSPHERE')
         if (associated(this%InStratosphere)) full_3d => this%InStratosphere
       case ('INTROPOSPHERE')
         if (associated(this%InTroposphere)) full_3d => this%InTroposphere
       case ('INPBL')
         if (associated(this%InPbl)) full_3d => this%InPbl
      end select

      if (associated(full_3d)) then
         column_ptr => full_3d(col_i, col_j, :)
      endif
   end function metstate_get_column_ptr_func_logical

   function metstate_get_2Dto0D_value_logical(this, field_name, i, j) result(scalar_val)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: i, j
      logical :: scalar_val
      integer :: col_i, col_j
      col_i = i
      col_j = j
      scalar_val = .false.
      select case (to_upper(trim(field_name)))
       case ('ISLAND')
         if (associated(this%IsLand)) scalar_val = this%IsLand(col_i, col_j)
       case ('ISWATER')
         if (associated(this%IsWater)) scalar_val = this%IsWater(col_i, col_j)
       case ('ISICE')
         if (associated(this%IsIce)) scalar_val = this%IsIce(col_i, col_j)
       case ('ISSNOW')
         if (associated(this%IsSnow)) scalar_val = this%IsSnow(col_i, col_j)
       case ('ISLOCALNOON')
         if (associated(this%IsLocalNoon)) scalar_val = this%IsLocalNoon(col_i, col_j)
      end select
   end function metstate_get_2Dto0D_value_logical

   function metstate_get_scalar_value_logical(this, field_name) result(scalar_val)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      logical :: scalar_val
      scalar_val = .false.
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

      real(fp), pointer :: temp_col_ptr(:)

      rc = CC_FAILURE
      nullify(col_ptr)

      ! Call the targeted getter function to fetch the slice
      temp_col_ptr => this%get_column_ptr_func(field_name, i, j)
      if (associated(temp_col_ptr)) then
         col_ptr => temp_col_ptr
         rc = CC_SUCCESS
      end if
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
      select case (to_upper(trim(adjustl(field_name))))
       case ('T')
         if (associated(this%T)) this%T = field_data
       case ('QV')
         if (associated(this%QV)) this%QV = field_data
       case ('RH')
         if (associated(this%RH)) this%RH = field_data
       case ('PMID')
         if (associated(this%PMID)) this%PMID = field_data
       case ('PEDGE')
         if (associated(this%PEDGE)) this%PEDGE = field_data
       case ('AIRDEN')
         if (associated(this%AIRDEN)) this%AIRDEN = field_data
       case ('AIRDEN_DRY')
         if (associated(this%AIRDEN_DRY)) this%AIRDEN_DRY = field_data
       case ('BXHEIGHT')
         if (associated(this%BXHEIGHT)) this%BXHEIGHT = field_data
       case ('DELP')
         if (associated(this%DELP)) this%DELP = field_data
       case ('DELP_DRY')
         if (associated(this%DELP_DRY)) this%DELP_DRY = field_data
       case ('PS')
         if (associated(this%PS)) this%PS = field_data
       case ('TS')
         if (associated(this%TS)) this%TS = field_data
       case ('PBLH')
         if (associated(this%PBLH)) this%PBLH = field_data
       case ('USTAR')
         if (associated(this%USTAR)) this%USTAR = field_data
       case ('HFLUX')
         if (associated(this%HFLUX)) this%HFLUX = field_data
       case ('OBK')
         if (associated(this%OBK)) this%OBK = field_data
       case ('LAT')
         if (associated(this%LAT)) this%LAT = field_data
       case ('LON')
         if (associated(this%LON)) this%LON = field_data
       case ('FRVEG')
         if (associated(this%FRVEG)) this%FRVEG = field_data
       case ('AREA_M2')
         if (associated(this%AREA_M2)) this%AREA_M2 = field_data
       case default
         rc = CC_FAILURE
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
      select case (to_upper(trim(adjustl(field_name))))
       case ('NLEVS')
         this%NLEVS = field_data
       case ('NSOIL')
         this%nSOIL = field_data
       case ('NSOILTYPE')
         this%nSOILTYPE = field_data
       case ('LWI')
         if (associated(this%LWI)) this%LWI = field_data
       case ('DLUSE')
         if (associated(this%DLUSE)) this%DLUSE = field_data
       case ('DSOILTYPE')
         if (associated(this%DSOILTYPE)) this%DSOILTYPE = field_data
       case ('NLNDTYPE')
         if (associated(this%nLNDTYPE)) this%nLNDTYPE = field_data
       case ('TROPLEV')
         if (associated(this%TropLev)) this%TropLev = field_data
       case default
         rc = CC_FAILURE
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
      select case (to_upper(trim(adjustl(field_name))))
       case ('ISLAND')
         if (associated(this%IsLand)) this%IsLand = field_data
       case ('ISWATER')
         if (associated(this%IsWater)) this%IsWater = field_data
       case ('ISICE')
         if (associated(this%IsIce)) this%IsIce = field_data
       case ('ISSNOW')
         if (associated(this%IsSnow)) this%IsSnow = field_data
       case ('ISLOCALNOON')
         if (associated(this%IsLocalNoon)) this%IsLocalNoon = field_data
       case ('INPBL')
         if (associated(this%InPbl)) this%InPbl = field_data
       case default
         rc = CC_FAILURE
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

      rc = CC_SUCCESS
      select case (to_upper(trim(adjustl(field_name))))
       case ('PS')
         if (associated(this%PS)) this%PS = field_data
       case ('TS')
         if (associated(this%TS)) this%TS = field_data
       case ('PBLH')
         if (associated(this%PBLH)) this%PBLH = field_data
       case ('USTAR')
         if (associated(this%USTAR)) this%USTAR = field_data
       case ('HFLUX')
         if (associated(this%HFLUX)) this%HFLUX = field_data
       case ('OBK')
         if (associated(this%OBK)) this%OBK = field_data
       case ('LAT')
         if (associated(this%LAT)) this%LAT = field_data
       case ('LON')
         if (associated(this%LON)) this%LON = field_data
       case ('FRVEG')
         if (associated(this%FRVEG)) this%FRVEG = field_data
       case ('AREA_M2')
         if (associated(this%AREA_M2)) this%AREA_M2 = field_data
       case default
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

      rc = CC_SUCCESS
      select case (to_upper(trim(adjustl(field_name))))
       case ('LWI')
         if (associated(this%LWI)) this%LWI = field_data
       case ('DLUSE')
         if (associated(this%DLUSE)) this%DLUSE = field_data
       case ('DSOILTYPE')
         if (associated(this%DSOILTYPE)) this%DSOILTYPE = field_data
       case ('NLNDTYPE')
         if (associated(this%nLNDTYPE)) this%nLNDTYPE = field_data
       case ('TROPLEV')
         if (associated(this%TropLev)) this%TropLev = field_data
       case default
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

      rc = CC_SUCCESS
      select case (to_upper(trim(adjustl(field_name))))
       case ('ISLAND')
         if (associated(this%IsLand)) this%IsLand = field_data
       case ('ISWATER')
         if (associated(this%IsWater)) this%IsWater = field_data
       case ('ISICE')
         if (associated(this%IsIce)) this%IsIce = field_data
       case ('ISSNOW')
         if (associated(this%IsSnow)) this%IsSnow = field_data
       case ('ISLOCALNOON')
         if (associated(this%IsLocalNoon)) this%IsLocalNoon = field_data
       case default
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

      rc = CC_SUCCESS
      select case (to_upper(trim(adjustl(field_name))))
       case ('T')
         if (associated(this%T)) this%T = field_data
       case ('QV')
         if (associated(this%QV)) this%QV = field_data
       case ('RH')
         if (associated(this%RH)) this%RH = field_data
       case ('PMID')
         if (associated(this%PMID)) this%PMID = field_data
       case ('PEDGE')
         if (associated(this%PEDGE)) this%PEDGE = field_data
       case ('AIRDEN')
         if (associated(this%AIRDEN)) this%AIRDEN = field_data
       case ('AIRDEN_DRY')
         if (associated(this%AIRDEN_DRY)) this%AIRDEN_DRY = field_data
       case ('BXHEIGHT')
         if (associated(this%BXHEIGHT)) this%BXHEIGHT = field_data
       case ('DELP')
         if (associated(this%DELP)) this%DELP = field_data
       case ('DELP_DRY')
         if (associated(this%DELP_DRY)) this%DELP_DRY = field_data
       case ('U')
         if (associated(this%U)) this%U = field_data
       case ('V')
         if (associated(this%V)) this%V = field_data
       case default
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

      rc = CC_SUCCESS
      select case (to_upper(trim(adjustl(field_name))))
       case ('ILAND')
         if (associated(this%ILAND)) this%ILAND = field_data
       case default
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

      rc = CC_SUCCESS
      select case (to_upper(trim(adjustl(field_name))))
       case ('INPBL')
         if (associated(this%InPbl)) this%InPbl = field_data
       case ('INSTRATMESO')
         if (associated(this%InStratMeso)) this%InStratMeso = field_data
       case ('INSTRATOSPHERE')
         if (associated(this%InStratosphere)) this%InStratosphere = field_data
       case ('INTROPOSPHERE')
         if (associated(this%InTroposphere)) this%InTroposphere = field_data
       case default
         rc = CC_FAILURE
      end select
   end subroutine metstate_set_field_3d_logical

   subroutine metstate_set_multiple_fields(this, field_names, error_mgr, rc, &
      AREA_M2_data, FRVEG_data, LAT_data, LON_data, IsLand_data, IsWater_data, NLEVS_data, LWI_data, &
      T_data, U_data, V_data, InPbl_data, nSOIL_data, nSOILTYPE_data, QV_data, RH_data, TropLev_data)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_names(:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      real(fp), optional, intent(in) :: AREA_M2_data(:,:)
      real(fp), optional, intent(in) :: FRVEG_data(:,:)
      real(fp), optional, intent(in) :: LAT_data(:,:)
      real(fp), optional, intent(in) :: LON_data(:,:)
      logical, optional, intent(in) :: IsLand_data(:,:)
      logical, optional, intent(in) :: IsWater_data(:,:)
      integer, optional, intent(in) :: NLEVS_data
      integer, optional, intent(in) :: LWI_data(:,:)
      real(fp), optional, intent(in) :: T_data(:,:,:)
      real(fp), optional, intent(in) :: U_data(:,:,:)
      real(fp), optional, intent(in) :: V_data(:,:,:)
      logical, optional, intent(in) :: InPbl_data(:,:,:)
      integer, optional, intent(in) :: nSOIL_data
      integer, optional, intent(in) :: nSOILTYPE_data
      real(fp), optional, intent(in) :: QV_data(:,:,:)
      real(fp), optional, intent(in) :: RH_data(:,:,:)
      integer, optional, intent(in) :: TropLev_data(:,:)

      integer :: idx
      character(len=32) :: name

      rc = CC_SUCCESS

      do idx = 1, size(field_names)
         name = trim(adjustl(field_names(idx)))
         select case (name)
          case ('AREA_M2', 'area_m2')
            if (present(AREA_M2_data)) then
               call this%set_field(name, AREA_M2_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('FRVEG', 'frveg')
            if (present(FRVEG_data)) then
               call this%set_field(name, FRVEG_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('LAT', 'lat')
            if (present(LAT_data)) then
               call this%set_field(name, LAT_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('LON', 'lon')
            if (present(LON_data)) then
               call this%set_field(name, LON_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('IsLand', 'island')
            if (present(IsLand_data)) then
               call this%set_field(name, IsLand_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('IsWater', 'iswater')
            if (present(IsWater_data)) then
               call this%set_field(name, IsWater_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('NLEVS', 'nlevs')
            if (present(NLEVS_data)) then
               call this%set_field(name, NLEVS_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('LWI', 'lwi')
            if (present(LWI_data)) then
               call this%set_field(name, LWI_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('T', 't')
            if (present(T_data)) then
               call this%set_field(name, T_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('U', 'u')
            if (present(U_data)) then
               call this%set_field(name, U_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('V', 'v')
            if (present(V_data)) then
               call this%set_field(name, V_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('InPbl', 'inpbl')
            if (present(InPbl_data)) then
               call this%set_field(name, InPbl_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('nSOIL', 'nsoil')
            if (present(nSOIL_data)) then
               call this%set_field(name, nSOIL_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('nSOILTYPE', 'nsoiltype')
            if (present(nSOILTYPE_data)) then
               call this%set_field(name, nSOILTYPE_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('QV', 'qv')
            if (present(QV_data)) then
               call this%set_field(name, QV_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('RH', 'rh')
            if (present(RH_data)) then
               call this%set_field(name, RH_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case ('TropLev', 'troplev')
            if (present(TropLev_data)) then
               call this%set_field(name, TropLev_data, error_mgr, rc)
            else
               rc = CC_FAILURE
            end if
          case default
            rc = CC_FAILURE
         end select
         if (rc /= CC_SUCCESS) return
      end do
   end subroutine metstate_set_multiple_fields

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
      real(fp) :: airden, rh, air_mass
      real(fp) :: avgw ! Water vapor volume mixing ratio [v/v dry air]
      real(fp) :: xh2o ! Water vapor mole fraction [mol (H2O) / mol (moist air)]
      !some variables used for reevaporation calculations
      real(fp) :: flux_liq, flux_ice, flux_tot, reevap_liq, reevap_ice,C_evap,RH_term, frac_liq
      real(fp), parameter :: C_evap_liq = 2.0e-5_fp  ! liquid evap coefficient
      real(fp), parameter :: C_evap_ice = 0.5e-5_fp  ! ice sublimation coefficient
      real(fp), parameter :: b0         = 0.9_fp     ! Sundqvist RH threshold
      real(fp), parameter :: T_liq      = 273.15_fp  ! K - liquid threshold
      real(fp), parameter :: T_ice      = 258.15_fp  ! K - ice threshold

      thisLoc = 'metstate_derive_field (in core/metstate_mod.F90)'
      call error_mgr%push_context('metstate_derive_field', 'deriving field: ' // trim(field_name))

      rc = CC_SUCCESS
      call this%get_dimensions(nx, ny, nz)

      select case (to_upper(trim(adjustl(field_name))))

       case ('MAIRDEN', 'mairden', 'AIRDEN', 'airden')
         ! Calculate dry air density from pressure and temperature
         ! ρ = P / (R_specific * T) where R_specific = R / MW
         if (.not. associated(this%PMID) .or. .not. associated(this%T)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'PMID and T fields required for MAIRDEN/AIRDEN calculation', rc, &
               thisLoc, 'Ensure pressure and temperature are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate MAIRDEN if not already allocated
         if (.not. associated(this%MAIRDEN) .or. .not. associated(this%AIRDEN)) then
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
         if (.not. associated(this%PMID) .or. .not. associated(this%T)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'PMID and T fields required for AIRDEN_DRY calculation', rc, &
               thisLoc, 'Ensure pressure and temperature are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate AIRDEN_DRY if not already allocated
         if (.not. associated(this%AIRDEN_DRY) .or. .not. associated(this%PMID_DRY) .or. &
            .not. associated(this%PEDGE_DRY) .or. .not. associated(this%DELP_DRY)) then
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
                  !if (k == nz) then
                  this%PEDGE_DRY(i, j, k+1) = this%PEDGE(i, j, k+1) * ( 1.e+0_fp - xh2o )
                  !end if
                  this%DELP_DRY(i, j, k) = this%PEDGE_DRY(i, j, k) - this%PEDGE_DRY(i, j, k+1)
               enddo
            enddo
         enddo

       case ('RH', 'rh')
         ! Calculate virtual temperature from temperature and humidity
         if (.not. associated(this%T) .or. .not. associated(this%QV) .or. .not. associated(this%PMID)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'T, PMID and QV fields required for RH calculation', rc, &
               thisLoc, 'Ensure temperature, pressure and humidity are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate RH if not already allocated
         if (.not. associated(this%RH)) then
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
         if (.not. associated(this%T) .or. .not. associated(this%QV)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'T and QV fields required for TV calculation', rc, &
               thisLoc, 'Ensure temperature and humidity are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate TV if not already allocated
         if (.not. associated(this%TV)) then
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
         if (.not. associated(this%HFLUX) .or. .not. associated(this%AIRDEN) .or. .not. associated(this%TS) .or. &
            .not. associated(this%USTAR)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'TS, USTAR, AIRDEN and HFLUX fields required for OBK calculation', rc, &
               thisLoc, 'Ensure temperature, ustar, air density, and sensible heat flux are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate OBK if not already allocated
         if (.not. associated(this%OBK)) then
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
         if (.not. associated(this%LAT) .or. .not. associated(this%LON)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'LAT and LON fields required for SUNCOS calculation', rc, &
               thisLoc, 'Ensure latitude and longitude are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate OBK if not already allocated
         if (.not. associated(this%SUNCOS)) then
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
         if (.not. associated(this%LAT) .or. .not. associated(this%LON)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'LAT and LON fields required for SUNCOSmid calculation', rc, &
               thisLoc, 'Ensure latitude and longitude are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate OBK if not already allocated
         if (.not. associated(this%SUNCOSmid)) then
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
         if (.not. associated(this%PEDGE)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'PEDGE field required for DELP calculation', rc, &
               thisLoc, 'Ensure pressure edges are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate BXHEIGHT if not already allocated
         if (.not. associated(this%DELP)) then
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
         if (.not. associated(this%PEDGE)) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
               'PEDGE field required for BXHEIGHT calculation', rc, &
               thisLoc, 'Ensure pressure edges are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate BXHEIGHT if not already allocated
         if (.not. associated(this%BXHEIGHT)) then
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

       case ('FRLANDIC', 'frlandic')
         this%FRLANDIC(:,:) = 0.0_fp !set to zero if IsIce is false
         do j = 1, ny
            do i = 1, nx
               if (abs(this%LWI(i, j) - 2.0_fp) < 0.5_fp) this%FRLANDIC(i, j) = 1.0_fp
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
         if (.not. associated(this%FRLANDUSE)) allocate(this%FRLANDUSE(nx, ny, nlanduse))
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
         if (.not. associated(this%ILAND)) allocate(this%ILAND(nx, ny, nlanduse))
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
         if (.not. associated(this%FRLAI)) allocate(this%FRLAI(nx, ny, nlanduse))
         this%FRLAI(:,:,:) = 0.0_fp
         do j = 1, ny
            do i = 1, nx
               do k = 1, nlanduse
                  if (this%DLUSE(i, j) == k) this%FRLAI(i, j, k) = this%LAI(i, j) !TODO: should times fraclanduse but here is 1.0
               enddo
               this%FRLAI(i, j, 15:17) = 0.0 !manually give index 15(snow and ice), 16(barren), 17(water) zeros
            enddo
         enddo
       case ('CLAYFRAC', 'clayfrac', 'SANDFRAC', 'sandfrac', 'SSM', 'ssm', 'RDRAG', 'rdrag', 'USTAR_THRESHOLD', 'ustar_threshold')
         !place holder. These are read in from emission read module for now. Here is to make sure required_met is all set.
         write(*,'(A)') 'Warning: Some Fengsha related met fields are read in from emission module, which will be disabled in the future!'
       case ('SALINITY', 'salinity')
         this%SALINITY(:,:) = 0.0_fp  !set to zero for now, which will turn off O3 dry deposition over ocean with iodine.

       case ('REEVAPLS', 'reevapls')
         this%REEVAPLS(:,:,:) = 0.0_fp  !I did not find data from GFS. Try to calculate it here.
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  ! 1. GET PRECIPITATION FLUXES
                  flux_liq = this%PFLLSAN(i,j,k)           ! kg/m²/s
                  flux_ice = this%PFILSAN(i,j,k)           ! kg/m²/s
                  flux_tot = flux_liq + flux_ice           ! kg/m²/s
                  ! Skip if no precipitation
                  if(flux_tot .le. 0.) cycle
                  ! Skip if already saturated (no evaporation possible)
                  rh = relative_humidity(this%T(i, j, k), this%QV(i, j, k), this%PMID(i, j, k))
                  if(rh .ge. b0) cycle
                  ! 2. LAYER THICKNESS AND AIR MASS
                  ! lower edge - upper edge
                  air_mass = (this%PEDGE(i, j, k) - this%PEDGE(i, j, k+1)) / g0
                  if(air_mass .le. 0.0_fp) cycle
                  ! 3. TEMPERATURE-DEPENDENT EVAPORATION COEFFICIENT Abel & Boutle (2012)
                  if(this%T(i,j,k) .gt. T_liq) then
                     ! Pure liquid
                     C_evap = C_evap_liq

                  else if(this%T(i,j,k) .gt. T_ice) then
                     ! Mixed phase - linear interpolation
                     frac_liq = (this%T(i,j,k) - T_ice) / (T_liq - T_ice)
                     C_evap   = frac_liq * C_evap_liq + &
                        (1. - frac_liq) * C_evap_ice
                  else
                     ! Pure ice
                     C_evap = C_evap_ice
                  endif
                  ! 4. SUNDQVIST (1988) RH TERM
                  RH_term = MAX(0., 1. - RH/b0)   ! dimensionless

                  ! 5. LIQUID EVAPORATION
                  !    kg/m²/s → kg/kg/s: divide by air_mass
                  !    note: air_mass cancels with numerator
                  if(flux_liq .gt. 0.) then
                     ! kg/m²/s version: C_evap * RH_term * sqrt(flux) * air_mass
                     ! kg/kg/s version: air_mass cancels → simpler
                     reevap_liq = C_evap          &
                        * RH_term                 &   ! dimensionless
                        * sqrt(flux_liq)              ! (kg/m²/s)^0.5

                     ! Physical constraint: cannot exceed available flux
                     ! Convert flux to kg/kg/s for comparison
                     reevap_liq = MIN(reevap_liq, flux_liq/air_mass)
                     reevap_liq = MAX(0., reevap_liq)
                  else
                     reevap_liq = 0.
                  endif

                  ! 6. ICE SUBLIMATION
                  !    Only above T_ice threshold
                  if(flux_ice .gt. 0. .and. this%T(i,j,k) .gt. T_ice) then
                     reevap_ice = C_evap          &
                        * RH_term                 &
                        * sqrt(flux_ice)

                     reevap_ice = MIN(reevap_ice, flux_ice/air_mass)
                     reevap_ice = MAX(0., reevap_ice)
                  else
                     reevap_ice = 0.
                  endif

                  ! 7. TOTAL REEVAPORATION IN kg/kg/s
                  this%REEVAPLS(i,j,k) = reevap_liq + reevap_ice

                  ! Final safety constraint
                  this%REEVAPLS(i,j,k) = MAX(0., this%REEVAPLS(i,j,k))
                  this%REEVAPLS(i,j,k) = MIN(this%REEVAPLS(i,j,k), flux_tot/air_mass)

               enddo
            enddo
         enddo

       case default
         call error_mgr%report_error(ERROR_NOT_FOUND, &
            'Unknown derived field: ' // trim(field_name), rc, &
            thisLoc, 'Supported fields: AIRDEN,  TV,  BXHEIGHT')
         rc = CC_FAILURE
      end select

      call error_mgr%pop_context()
   end subroutine metstate_derive_field

   !> @brief Private helper to convert a string to uppercase for case-insensitive matching
   function to_upper(str) result(upper_str)
      implicit none
      character(len=*), intent(in) :: str
      character(len=len(str)) :: upper_str
      integer :: i, char_code
      do i = 1, len(str)
         char_code = ichar(str(i:i))
         if (char_code >= 97 .and. char_code <= 122) then
            upper_str(i:i) = char(char_code - 32)
         else
            upper_str(i:i) = str(i:i)
         end if
      end do
   end function to_upper

END MODULE MetState_Mod
