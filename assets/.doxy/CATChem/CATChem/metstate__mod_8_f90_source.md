

# File metstate\_mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**metstate\_mod.F90**](metstate__mod_8_f90.md)

[Go to the documentation of this file](metstate__mod_8_f90.md)


```Fortran
! \file metstate_mod.F90
!! \brief Module for meteorology state variables
!!
!! This module contains subroutines and functions related to the MetStateType instance of CATChem.
!! It includes subroutines for initializing of the MetStateType.
!!
!! \ingroup core_modules
!!!>
MODULE metstate_mod
   !
   ! USES:
   !
   USE error_mod
   USE precision_mod
   USE gridgeometry_mod
   USE met_utilities_mod
   USE timestate_mod, only: timestatetype



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
   TYPE, PUBLIC :: metstatetype
      CHARACTER(LEN=3)             :: State     = 'MET'
      INTEGER                      :: NLEVS     = 127
      TYPE(GridGeometryType) :: geometry
      INTEGER                      :: NSURFTYPE = 20
      ! Grid flags (2D: nx, ny)
      LOGICAL, ALLOCATABLE         :: IsLand(:,:)
      LOGICAL, ALLOCATABLE         :: IsWater(:,:)
      LOGICAL, ALLOCATABLE         :: IsIce(:,:)
      LOGICAL, ALLOCATABLE         :: IsSnow(:,:)
      ! Vertical flags and arrays (3D: nx, ny, nz)
      LOGICAL,  ALLOCATABLE        :: InStratMeso(:,:,:)
      LOGICAL,  ALLOCATABLE        :: InStratosphere(:,:,:)
      LOGICAL,  ALLOCATABLE        :: InTroposphere(:,:,:)
      LOGICAL,  ALLOCATABLE        :: InPbl(:,:,:)
      LOGICAL,  ALLOCATABLE        :: IsLocalNoon(:,:)
      ! Surface properties (2D: nx, ny)
      REAL(fp), ALLOCATABLE        :: AREA_M2(:,:)
      INTEGER,  ALLOCATABLE        :: LWI(:,:)
      INTEGER,  ALLOCATABLE        :: DLUSE(:,:)
      REAL(fp), ALLOCATABLE        :: FRVEG(:,:)
      REAL(fp), ALLOCATABLE        :: FRLAKE(:,:)
      REAL(fp), ALLOCATABLE        :: FRLAND(:,:)
      REAL(fp), ALLOCATABLE        :: FRLANDIC(:,:)
      REAL(fp), ALLOCATABLE        :: FROCEAN(:,:)
      REAL(fp), ALLOCATABLE        :: FRSEAICE(:,:)
      REAL(fp), ALLOCATABLE        :: FRSNO(:,:)
      REAL(fp), ALLOCATABLE        :: LAI(:,:)
      REAL(fp), ALLOCATABLE        :: GVF(:,:)
      ! Dust Only Variables
      REAL(fp), ALLOCATABLE        :: RDRAG(:,:)
      REAL(fp), ALLOCATABLE        :: USTAR_THRESHOLD(:,:)
      REAL(fp), ALLOCATABLE        :: SSM(:,:)
      ! Surface and ice properties (2D: nx, ny)
      REAL(fp), ALLOCATABLE        :: SEAICE00(:,:)
      REAL(fp), ALLOCATABLE        :: SEAICE10(:,:)
      REAL(fp), ALLOCATABLE        :: SEAICE20(:,:)
      REAL(fp), ALLOCATABLE        :: SEAICE30(:,:)
      REAL(fp), ALLOCATABLE        :: SEAICE40(:,:)
      REAL(fp), ALLOCATABLE        :: SEAICE50(:,:)
      REAL(fp), ALLOCATABLE        :: SEAICE60(:,:)
      REAL(fp), ALLOCATABLE        :: SEAICE70(:,:)
      REAL(fp), ALLOCATABLE        :: SEAICE80(:,:)
      REAL(fp), ALLOCATABLE        :: SEAICE90(:,:)
      REAL(fp), ALLOCATABLE        :: SNODP(:,:)
      REAL(fp), ALLOCATABLE        :: SNOMAS(:,:)

      ! Soil and land use arrays (2D for counts, 3D for fractions)
      INTEGER,  ALLOCATABLE        :: DSOILTYPE(:,:)
      REAL(fp), ALLOCATABLE        :: CLAYFRAC(:,:)
      REAL(fp), ALLOCATABLE        :: SANDFRAC(:,:)
      INTEGER,  ALLOCATABLE        :: nLNDTYPE(:,:)
      REAL(fp), ALLOCATABLE        :: GWETTOP(:,:)
      REAL(fp), ALLOCATABLE        :: GWETROOT(:,:)
      REAL(fp), ALLOCATABLE        :: WILT(:,:)
      INTEGER                      :: nSOIL
      INTEGER                      :: nSOILTYPE
      REAL(fp), ALLOCATABLE        :: SOILM(:,:,:)
      REAL(fp), ALLOCATABLE        :: SOILT(:,:,:)
      REAL(fp), ALLOCATABLE        :: FRLANDUSE(:,:,:)
      REAL(fp), ALLOCATABLE        :: FRSOIL(:,:,:)
      REAL(fp), ALLOCATABLE        :: FRLAI(:,:,:)
      INTEGER, ALLOCATABLE         :: ILAND(:,:,:)
      ! Location arrays (1D for single point, 2D for grid)
      real(fp), ALLOCATABLE        :: LAT(:,:)
      real(fp), ALLOCATABLE        :: LON(:,:)
      character(len=20)            :: LUCNAME
      ! Surface meteorological properties (2D: nx, ny)
      REAL(fp), ALLOCATABLE        :: ALBD_VIS(:,:)
      REAL(fp), ALLOCATABLE        :: ALBD_NIR(:,:)
      REAL(fp), ALLOCATABLE        :: ALBD_UV(:,:)
      REAL(fp), ALLOCATABLE        :: PARDR(:,:)
      REAL(fp), ALLOCATABLE        :: PARDF(:,:)
      REAL(fp), ALLOCATABLE        :: SUNCOS(:,:)
      REAL(fp), ALLOCATABLE        :: SUNCOSmid(:,:)
      REAL(fp), ALLOCATABLE        :: SUNCOSsum(:,:)
      REAL(fp), ALLOCATABLE        :: SZAFACT(:,:)
      REAL(fp), ALLOCATABLE        :: SWGDN(:,:)
      REAL(fp), ALLOCATABLE        :: EFLUX(:,:)
      REAL(fp), ALLOCATABLE        :: HFLUX(:,:)
      REAL(fp), ALLOCATABLE        :: U10M(:,:)
      REAL(fp), ALLOCATABLE        :: USTAR(:,:)
      REAL(fp), ALLOCATABLE        :: V10M(:,:)
      REAL(fp), ALLOCATABLE        :: Z0(:,:)
      REAL(fp), ALLOCATABLE        :: Z0H(:,:)
      REAL(fp), ALLOCATABLE        :: FRZ0(:,:,:)
      REAL(fp), ALLOCATABLE        :: PBLH(:,:)
      REAL(fp), ALLOCATABLE        :: SALINITY(:,:)
      REAL(fp), ALLOCATABLE        :: CMM(:,:)
      REAL(fp), ALLOCATABLE        :: ORO(:,:)
      REAL(fp), ALLOCATABLE        :: RCA(:,:)
      REAL(fp), ALLOCATABLE        :: WCA(:,:)          ! canopy water amount [kg/m2]
      ! 3D volumetric fields (3D: nx, ny, nz)
      REAL(fp), ALLOCATABLE        :: F_OF_PBL(:,:,:)
      REAL(fp), ALLOCATABLE        :: F_UNDER_PBLTOP(:,:,:)
      real(fp), ALLOCATABLE        :: OBK(:,:)
      ! Cloud and precipitation properties (2D for surface, 3D for volumetric)
      REAL(fp), ALLOCATABLE        :: CLDFRC(:,:)
      REAL(fp), ALLOCATABLE        :: CONV_DEPTH(:,:)
      REAL(fp), ALLOCATABLE        :: FLASH_DENS(:,:)
      REAL(fp), ALLOCATABLE        :: CNV_FRC(:,:)
      REAL(fp), ALLOCATABLE        :: CLDF(:,:,:)
      REAL(fp), ALLOCATABLE        :: CMFMC(:,:,:)
      REAL(fp), ALLOCATABLE        :: DQRCU(:,:,:)
      REAL(fp), ALLOCATABLE        :: DQRLSAN(:,:,:)
      REAL(fp), ALLOCATABLE        :: DTRAIN(:,:,:)
      REAL(fp), ALLOCATABLE        :: PRECANV(:,:)
      REAL(fp), ALLOCATABLE        :: PRECCON(:,:)
      REAL(fp), ALLOCATABLE        :: PRECLSC(:,:)
      real(fp), ALLOCATABLE        :: REEVAPLS(:,:,:)
      ! 3D cloud and precipitation arrays
      REAL(fp), ALLOCATABLE        :: QI(:,:,:)
      REAL(fp), ALLOCATABLE        :: QL(:,:,:)
      REAL(fp), ALLOCATABLE        :: PFICU(:,:,:)
      REAL(fp), ALLOCATABLE        :: PFILSAN(:,:,:)
      REAL(fp), ALLOCATABLE        :: PFLCU(:,:,:)
      REAL(fp), ALLOCATABLE        :: PFLLSAN(:,:,:)
      REAL(fp), ALLOCATABLE        :: TAUCLI(:,:,:)
      REAL(fp), ALLOCATABLE        :: TAUCLW(:,:,:)
      ! Surface scalars (now 2D: nx, ny)
      REAL(fp), ALLOCATABLE        :: PHIS(:,:)
      REAL(fp), ALLOCATABLE        :: PS_WET(:,:)
      REAL(fp), ALLOCATABLE        :: PS_DRY(:,:)
      REAL(fp), ALLOCATABLE        :: QV2M(:,:)
      REAL(fp), ALLOCATABLE        :: T2M(:,:)
      REAL(fp), ALLOCATABLE        :: TS(:,:)
      REAL(fp), ALLOCATABLE        :: TSKIN(:,:)
      REAL(fp), ALLOCATABLE        :: SST(:,:)
      REAL(fp), ALLOCATABLE        :: SLP(:,:)
      REAL(fp), ALLOCATABLE        :: PS(:,:)
      REAL(fp), ALLOCATABLE        :: TO3(:,:)
      REAL(fp), ALLOCATABLE        :: TROPP(:,:)
      INTEGER,  ALLOCATABLE        :: TropLev(:,:)
      REAL(fp), ALLOCATABLE        :: TropHt(:,:)
      ! 3D atmospheric variables (3D: nx, ny, nz)
      REAL(fp), ALLOCATABLE        :: Z(:,:,:)
      REAL(fp), ALLOCATABLE        :: ZMID(:,:,:)
      REAL(fp), ALLOCATABLE        :: BXHEIGHT(:,:,:)
      REAL(fp), ALLOCATABLE        :: QV(:,:,:)
      REAL(fp), ALLOCATABLE        :: T(:,:,:)
      REAL(fp), ALLOCATABLE        :: THETA(:,:,:)
      REAL(fp), ALLOCATABLE        :: TV(:,:,:)
      REAL(fp), ALLOCATABLE        :: V(:,:,:)
      REAL(fp), ALLOCATABLE        :: U(:,:,:)
      REAL(fp), ALLOCATABLE        :: OMEGA(:,:,:)
      REAL(fp), ALLOCATABLE        :: RH(:,:,:)
      REAL(fp), ALLOCATABLE        :: SPHU(:,:,:)
      REAL(fp), ALLOCATABLE        :: AIRDEN(:,:,:)
      REAL(fp), ALLOCATABLE        :: AIRDEN_DRY(:,:,:)
      REAL(fp), ALLOCATABLE        :: AIRNUMDEN(:,:,:)
      REAL(fp), ALLOCATABLE        :: MAIRDEN(:,:,:)
      REAL(fp), ALLOCATABLE        :: AVGW(:,:,:)
      REAL(fp), ALLOCATABLE        :: DELP(:,:,:)
      REAL(fp), ALLOCATABLE        :: DELP_DRY(:,:,:)
      REAL(fp), ALLOCATABLE        :: DAIRMASS(:,:,:)
      REAL(fp), ALLOCATABLE        :: AIRVOL(:,:,:)
      REAL(fp), ALLOCATABLE        :: PEDGE_DRY(:,:,:)
      REAL(fp), ALLOCATABLE        :: PEDGE(:,:,:)
      REAL(fp), ALLOCATABLE        :: PMID(:,:,:)
      REAL(fp), ALLOCATABLE        :: PMID_DRY(:,:,:)
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
      procedure, public :: get_2Dto0D_value => metstate_get_2dto0d_value
      procedure, public :: get_2Dto0D_value_int => metstate_get_2dto0d_value_int
      procedure, public :: get_2Dto0D_value_logical => metstate_get_2dto0d_value_logical
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
   end type metstatetype

CONTAINS

   subroutine metstate_init(this, nx, ny, nlevs, nsoil, nsoiltype, nsurftype, error_mgr, rc)
      use error_mod, only: errormanagertype, cc_success, error_memory_allocation

      implicit none
      class(MetStateType), intent(inout) :: this
      integer, intent(in) :: nx, ny, nlevs
      integer, intent(in), optional :: nsoil, nsoiltype, nsurftype
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc
      character(len=256) :: thisLoc

      thisloc = 'metstate_init (in core/metstate_mod.F90)'
      call error_mgr%push_context('metstate_init', 'initializing meteorological state')

      rc = cc_success

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

   subroutine allocate_metstate_arrays(this, field_name, error_mgr, rc)
      use error_mod, only: errormanagertype, cc_success, error_memory_allocation

      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc
      character(len=256) :: thisLoc
      integer :: nx, ny, nz, nsoil, nsoiltype, nSURFTYPE

      thisloc = 'allocate_metstate_arrays (in core/metstate_mod.F90)'
      rc = cc_success

      call this%geometry%get_dimensions(nx, ny, nz)

      ! Use the properly initialized values (no more defaults needed)
      nsoil = this%nSOIL
      nsoiltype = this%nSOILTYPE
      nsurftype = this%NSURFTYPE

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

   subroutine metstate_cleanup(this, field_name, rc)
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(out) :: rc

      rc = cc_success

      ! Deallocate only the requested field (or all if 'ALL')
      select case (trim(field_name))
#include "metstate_deallocate_fields.inc"
      end select

      this%State = ''
      this%NLEVS = 72  ! Reset to default
      this%NSURFTYPE = 1  ! Reset to default

   end subroutine metstate_cleanup

   subroutine metstate_validate(this, error_mgr, rc)
      use error_mod, only: errormanagertype, cc_success, error_invalid_input
      use utilities_mod, only: is_valid_temperature, is_valid_pressure

      implicit none
      class(MetStateType), intent(in) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      character(len=256) :: thisLoc

      thisloc = 'metstate_validate (in core/metstate_mod.F90)'
      call error_mgr%push_context('metstate_validate', 'validating meteorological state')

      rc = cc_success

      ! Check basic state
      if (this%NLEVS <= 0) then
         call error_mgr%report_error(error_invalid_input, &
            'Number of levels must be positive', rc, &
            thisloc, 'Set NLEVS to a positive integer')
         call error_mgr%pop_context()
         return
      endif

      ! Validate temperatures (use maxval/minval for array validation)
      if (allocated(this%T2M)) then
         if (maxval(this%T2M) > 400.0_fp .or. minval(this%T2M) < 100.0_fp) then
            call error_mgr%report_error(error_invalid_input, &
               '2m temperature out of physical range', rc, &
               thisloc, 'Check temperature units and values')
            call error_mgr%pop_context()
            return
         endif
      endif

      if (allocated(this%TS)) then
         if (maxval(this%TS) > 400.0_fp .or. minval(this%TS) < 100.0_fp) then
            call error_mgr%report_error(error_invalid_input, &
               'Surface temperature out of physical range', rc, &
               thisloc, 'Check temperature units and values')
            call error_mgr%pop_context()
            return
         endif
      endif

      ! Validate pressures (use maxval/minval for array validation)
      if (allocated(this%PS)) then
         if (maxval(this%PS) > 120000.0_fp .or. minval(this%PS) < 1000.0_fp) then
            call error_mgr%report_error(error_invalid_input, &
               'Surface pressure out of physical range', rc, &
               thisloc, 'Check pressure units and values')
            call error_mgr%pop_context()
            return
         endif
      endif

      if (allocated(this%SLP)) then
         if (maxval(this%SLP) > 120000.0_fp .or. minval(this%SLP) < 50000.0_fp) then
            call error_mgr%report_error(error_invalid_input, &
               'Sea level pressure out of physical range', rc, &
               thisloc, 'Check pressure units and values')
            call error_mgr%pop_context()
            return
         endif
      endif

      ! Check array allocation
      if (.not. this%is_allocated()) then
         call error_mgr%report_error(error_invalid_input, &
            'Required arrays not allocated', rc, &
            thisloc, 'Call init() before using MetState')
         call error_mgr%pop_context()
         return
      endif

      call error_mgr%pop_context()
   end subroutine metstate_validate

   subroutine metstate_reset(this, rc)
      implicit none
      class(MetStateType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success

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

   function metstate_is_allocated(this) result(is_alloc)
      implicit none
      class(MetStateType), intent(in) :: this
      logical :: is_alloc

      is_alloc = allocated(this%T) .and. allocated(this%U) .and. allocated(this%V) .and. &
         allocated(this%QV) .and. allocated(this%PMID) .and. allocated(this%DELP)
   end function metstate_is_allocated

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

   subroutine metstate_allocate_field(this, field_name, rc)
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(out) :: rc
      integer :: nx, ny, nz, nsoil, nsoiltype, nSURFTYPE
      rc = cc_success
      call this%geometry%get_dimensions(nx, ny, nz)
      ! Use the properly initialized values (no more defaults needed)
      nsoil = this%nSOIL
      nsoiltype = this%nSOILTYPE
      nsurftype = this%NSURFTYPE
      ! Only allocate the requested field
      select case (trim(field_name))
#include "metstate_allocate_fields.inc"
      end select
   end subroutine metstate_allocate_field

   subroutine metstate_deallocate_field(this, field_name, rc)
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(out) :: rc
      rc = cc_success
      ! Only deallocate the requested field
      select case (trim(field_name))
#include "metstate_deallocate_fields.inc"
      end select
   end subroutine metstate_deallocate_field

   subroutine metstate_get_3dto1d_ptr(this, field_name, i, j, col_ptr, rc)
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
   end subroutine metstate_get_3dto1d_ptr

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

   function metstate_get_2dto0d_value(this, field_name, i, j) result(scalar_val)
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
   end function metstate_get_2dto0d_value

   function metstate_get_scalar_value(this, field_name) result(scalar_val)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      real(fp) :: scalar_val
      select case (trim(field_name))
#include "metstate_scalar_accessor.inc"
      end select
   end function metstate_get_scalar_value

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

   function metstate_get_2dto0d_value_int(this, field_name, i, j) result(scalar_val)
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
   end function metstate_get_2dto0d_value_int

   function metstate_get_scalar_value_int(this, field_name) result(scalar_val)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      integer :: scalar_val
      select case (trim(field_name))
#include "metstate_scalar_accessor_int.inc"
      end select
   end function metstate_get_scalar_value_int

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

   function metstate_get_2dto0d_value_logical(this, field_name, i, j) result(scalar_val)
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
   end function metstate_get_2dto0d_value_logical

   function metstate_get_scalar_value_logical(this, field_name) result(scalar_val)
      class(MetStateType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      logical :: scalar_val
      select case (trim(field_name))
#include "metstate_scalar_accessor_logical.inc"
      end select
   end function metstate_get_scalar_value_logical

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

   subroutine metstate_get_column_ptr_subroutine(this, field_name, i, j, col_ptr, rc)
      use error_mod, only: cc_success, cc_failure

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

      rc = cc_failure
      nullify(col_ptr)

      call this%get_dimensions(nx, ny, nlev)

      ! First try to get as a 3D field (vertical column)
      temp_col_ptr => this%get_column_ptr_func(field_name, i, j)
      if (associated(temp_col_ptr)) then
         col_ptr => temp_col_ptr
         rc = cc_success
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
         rc = cc_failure
         return

       case default
         ! Try as scalar field - similar limitation
         scalar_val = this%get_scalar_value(field_name)
         rc = cc_failure
         return
      end select

   end subroutine metstate_get_column_ptr_subroutine


   !---------------------------------------------------------------------------
   !                 Dimensional MetState Set Field Subroutines
   !---------------------------------------------------------------------------

   subroutine metstate_set_field_scalar_real(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: errormanagertype, cc_success, cc_failure
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: field_data
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = cc_success

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
               call error_mgr%report_error(error_not_found, &
                  'Unknown REAL field name: ' // trim(field_name), rc)
               rc = cc_failure
            end select
         end select
      end select
   end subroutine metstate_set_field_scalar_real

   subroutine metstate_set_field_scalar_int(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: errormanagertype, cc_success, cc_failure
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: field_data
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = cc_success

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
               call error_mgr%report_error(error_not_found, &
                  'Unknown INTEGER field name: ' // trim(field_name), rc)
               rc = cc_failure
            end select
         end select
      end select
   end subroutine metstate_set_field_scalar_int

   subroutine metstate_set_field_scalar_logical(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: errormanagertype, cc_success, cc_failure
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      logical, intent(in) :: field_data
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = cc_success

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
               call error_mgr%report_error(error_not_found, &
                  'Unknown LOGICAL field name: ' // trim(field_name), rc)
               rc = cc_failure
            end select
         end select
      end select
   end subroutine metstate_set_field_scalar_logical

   subroutine metstate_set_field_2d_real(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: errormanagertype, cc_success, cc_failure
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: field_data(:,:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Generated include file for 2D REAL field assignment
      rc = cc_success
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_2d_real.inc"
       case default
         call error_mgr%report_error(error_not_found, &
            "Unknown field name: " // trim(field_name), rc)
         rc = cc_failure
      end select
   end subroutine metstate_set_field_2d_real

   subroutine metstate_set_field_2d_int(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: errormanagertype, cc_success, cc_failure
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: field_data(:,:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Generated include file for 2D INTEGER field assignment
      rc = cc_success
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_2d_int.inc"
       case default
         call error_mgr%report_error(error_not_found, &
            "Unknown field name: " // trim(field_name), rc)
         rc = cc_failure
      end select
   end subroutine metstate_set_field_2d_int

   subroutine metstate_set_field_2d_logical(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: errormanagertype, cc_success, cc_failure
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      logical, intent(in) :: field_data(:,:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Generated include file for 2D LOGICAL field assignment
      rc = cc_success
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_2d_logical.inc"
       case default
         call error_mgr%report_error(error_not_found, &
            "Unknown field name: " // trim(field_name), rc)
         rc = cc_failure
      end select
   end subroutine metstate_set_field_2d_logical

   subroutine metstate_set_field_3d_real(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: errormanagertype, cc_success, cc_failure
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: field_data(:,:,:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Generated include file for 3D REAL field assignment
      rc = cc_success
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_3d_real.inc"
       case default
         call error_mgr%report_error(error_not_found, &
            "Unknown field name: " // trim(field_name), rc)
         rc = cc_failure
      end select
   end subroutine metstate_set_field_3d_real

   subroutine metstate_set_field_3d_int(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: errormanagertype, cc_success, cc_failure
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: field_data(:,:,:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Generated include file for 3D INTEGER field assignment
      rc = cc_success
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_3d_int.inc"
       case default
         call error_mgr%report_error(error_not_found, &
            "Unknown field name: " // trim(field_name), rc)
         rc = cc_failure
      end select
   end subroutine metstate_set_field_3d_int

   subroutine metstate_set_field_3d_logical(this, field_name, field_data, error_mgr, rc)
      use error_mod, only: errormanagertype, cc_success, cc_failure
      implicit none
      class(MetStateType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      logical, intent(in) :: field_data(:,:,:)
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Generated include file for 3D LOGICAL field assignment
      rc = cc_success
      select case (trim(adjustl(field_name)))
#include "metstate_set_field_3d_logical.inc"
       case default
         call error_mgr%report_error(error_not_found, &
            "Unknown field name: " // trim(field_name), rc)
         rc = cc_failure
      end select
   end subroutine metstate_set_field_3d_logical

! Include the auto-generated multiple fields interface
#include "metstate_multiple_fields_interface.inc"

   subroutine metstate_derive_field(this, field_name, error_mgr, time_state, rc)
      use error_mod, only: errormanagertype, cc_success, cc_failure, error_invalid_input, error_not_found
      use constants, only: g0, rd, rdg0, airmw, h2omw

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

      thisloc = 'metstate_derive_field (in core/metstate_mod.F90)'
      call error_mgr%push_context('metstate_derive_field', 'deriving field: ' // trim(field_name))

      rc = cc_success
      call this%get_dimensions(nx, ny, nz)

      select case (trim(adjustl(field_name)))

       case ('MAIRDEN', 'mairden', 'AIRDEN', 'airden')
         ! Calculate dry air density from pressure and temperature
         ! ρ = P / (R_specific * T) where R_specific = R / MW
         if (.not. allocated(this%PMID) .or. .not. allocated(this%T)) then
            call error_mgr%report_error(error_invalid_input, &
               'PMID and T fields required for MAIRDEN/AIRDEN calculation', rc, &
               thisloc, 'Ensure pressure and temperature are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate MAIRDEN if not already allocated
         if (.not. allocated(this%MAIRDEN) .or. .not. allocated(this%AIRDEN)) then
            call error_mgr%report_error(rc, 'MAIRDEN/AIRDEN fields need to be allocated first!', rc, thisloc)
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
            call error_mgr%report_error(error_invalid_input, &
               'PMID and T fields required for AIRDEN_DRY calculation', rc, &
               thisloc, 'Ensure pressure and temperature are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate AIRDEN_DRY if not already allocated
         if (.not. allocated(this%AIRDEN_DRY) .or. .not. allocated(this%PMID_DRY) .or. &
            .not. allocated(this%PEDGE_DRY) .or. .not. allocated(this%DELP_DRY)) then
            call error_mgr%report_error(rc, 'AIRDEN_DRY/PMID_DRY/PEDGE_DRY/DELP_DRY fields need to be allocated first!', rc, thisloc)
            call error_mgr%pop_context()
            return
         endif

         ! Calculate dry air density: ρ = P / (R_dry * T)
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  avgw = airmw * this%QV(i,j,k) / ( h2omw * (1.0e+0_fp - this%QV(i,j,k)) )
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
            call error_mgr%report_error(error_invalid_input, &
               'T, PMID and QV fields required for RH calculation', rc, &
               thisloc, 'Ensure temperature, pressure and humidity are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate RH if not already allocated
         if (.not. allocated(this%RH)) then
            call error_mgr%report_error(rc, 'RH field needs to be allocated first!', rc, thisloc)
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
            call error_mgr%report_error(error_invalid_input, &
               'T and QV fields required for TV calculation', rc, &
               thisloc, 'Ensure temperature and humidity are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate TV if not already allocated
         if (.not. allocated(this%TV)) then
            call error_mgr%report_error(rc, 'TV field needs to be allocated first!', rc, thisloc)
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
            call error_mgr%report_error(error_invalid_input, &
               'TS, USTAR, AIRDEN and HFLUX fields required for OBK calculation', rc, &
               thisloc, 'Ensure temperature, ustar, air density, and sensible heat flux are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate OBK if not already allocated
         if (.not. allocated(this%OBK)) then
            call error_mgr%report_error(rc, 'OBK field needs to be allocated first!', rc, thisloc)
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
            call error_mgr%report_error(error_invalid_input, &
               'LAT and LON fields required for SUNCOS calculation', rc, &
               thisloc, 'Ensure latitude and longitude are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate OBK if not already allocated
         if (.not. allocated(this%SUNCOS)) then
            call error_mgr%report_error(rc, 'SUNCOS field needs to be allocated first!', rc, thisloc)
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
            call error_mgr%report_error(error_invalid_input, &
               'LAT and LON fields required for SUNCOSmid calculation', rc, &
               thisloc, 'Ensure latitude and longitude are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate OBK if not already allocated
         if (.not. allocated(this%SUNCOSmid)) then
            call error_mgr%report_error(rc, 'SUNCOSmid field needs to be allocated first!', rc, thisloc)
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
            call error_mgr%report_error(error_invalid_input, &
               'PEDGE field required for DELP calculation', rc, &
               thisloc, 'Ensure pressure edges are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate BXHEIGHT if not already allocated
         if (.not. allocated(this%DELP)) then
            call error_mgr%report_error(rc, 'BXHEIGHT field needs to be allocated first!', rc, thisloc)
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
            call error_mgr%report_error(error_invalid_input, &
               'PEDGE field required for BXHEIGHT calculation', rc, &
               thisloc, 'Ensure pressure edges are available')
            call error_mgr%pop_context()
            return
         endif

         ! Allocate BXHEIGHT if not already allocated
         if (.not. allocated(this%BXHEIGHT)) then
            call error_mgr%report_error(rc, 'BXHEIGHT field needs to be allocated first!', rc, thisloc)
            call error_mgr%pop_context()
            return
         endif

         ! Calculate box height as difference between edge heights
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  ! Refer to https://github.com/geoschem/geos-chem/GeosCore/calc_met_mod.F90
                  this%BXHEIGHT(i, j, k) = rdg0 * virtual_temperature(this%T(i, j, k), this%QV(i, j, k)) * &
                     log(this%PEDGE(i, j, k) / this%PEDGE(i, j, k+1))
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
         call error_mgr%report_error(error_not_found, &
            'Unknown derived field: ' // trim(field_name), rc, &
            thisloc, 'Supported fields: AIRDEN,  TV,  BXHEIGHT')
         rc = cc_failure
      end select

      call error_mgr%pop_context()
   end subroutine metstate_derive_field

END MODULE metstate_mod
```


