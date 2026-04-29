!> \file catchem_regrid_mod.F90
!! \brief ESMF regridding utilities for CATChem emission data
!!
!! Provides runtime regridding from global lat-lon grids to the model
!! grid (e.g. cubed-sphere). Caches ESMF RouteHandles so regrid weights
!! are computed only once per unique source grid shape.

module catchem_regrid_mod

   use ESMF
   use netcdf

   implicit none
   private

   public :: RegridCache
   public :: catchem_regrid_field
   public :: catchem_regrid_cleanup

   !> Maximum number of cached regrid route handles
   integer, parameter :: MAX_REGRID_CACHE = 32

   !> \brief Entry in the regrid route-handle cache
   type :: RegridCacheEntry
      integer :: nlon = 0         !< Source grid longitude count
      integer :: nlat = 0         !< Source grid latitude  count
      type(ESMF_Grid)        :: srcGrid
      type(ESMF_Field)       :: srcField
      type(ESMF_RouteHandle) :: routeHandle
      logical :: active = .false.
   end type RegridCacheEntry

   !> \brief Global cache of regrid route handles
   type :: RegridCache
      type(RegridCacheEntry) :: entries(MAX_REGRID_CACHE)
      integer :: count = 0
   contains
      procedure :: lookup  => regrid_cache_lookup
      procedure :: add     => regrid_cache_add
      procedure :: cleanup => regrid_cache_cleanup
   end type RegridCache

contains

   !--------------------------------------------------------------------------
   !> \brief Read a 2D/3D variable from a global lat-lon file and regrid
   !!        it onto an ESMF field defined on the model grid.
   !!
   !! If the file dimensions match the model grid tile, data is read
   !! directly (no regrid). Otherwise a rectilinear source grid is
   !! constructed from the file's latitude/longitude coordinate variables,
   !! the data is read onto that grid, and ESMF bilinear regridding is
   !! applied to fill the destination field.
   !!
   !! \param[inout] cache       Regrid cache (weights reused across calls)
   !! \param[in]    filename    Path to NetCDF file
   !! \param[in]    varname     Variable name in the file
   !! \param[in]    dstField    Destination ESMF field (on model grid)
   !! \param[in]    latname     Name of latitude  coordinate variable
   !! \param[in]    lonname     Name of longitude coordinate variable
   !! \param[in]    timeSlice   Optional time-slice index
   !! \param[out]   didRegrid   Optional flag indicating regrid was applied
   !! \param[out]   rc          Return code
   !--------------------------------------------------------------------------
   subroutine catchem_regrid_field(cache, filename, varname, dstField, &
                                   latname, lonname, regrid_method_name, &
                                   timeSlice, levelSlice, didRegrid, rc)
      type(RegridCache),    intent(inout)         :: cache
      character(len=*),     intent(in)            :: filename
      character(len=*),     intent(in)            :: varname
      type(ESMF_Field),     intent(inout)         :: dstField
      character(len=*),     intent(in)            :: latname
      character(len=*),     intent(in)            :: lonname
      character(len=*),     intent(in),  optional :: regrid_method_name
      integer,              intent(in),  optional :: timeSlice
      integer,              intent(in),  optional :: levelSlice
      logical,              intent(out), optional :: didRegrid
      integer,              intent(out), optional :: rc

      ! -- local variables
      integer :: localrc, ncid, ncStatus
      integer :: lonDimId, latDimId, lonVarId, latVarId, varId
      integer :: nlon, nlat, ndims, xtype, uid
      integer :: timeDimLen, idx
      integer, allocatable :: dimids(:), elemStart(:), elemCount(:)
      character(len=ESMF_MAXSTR) :: dimName
      real(ESMF_KIND_R8), allocatable :: lonCoord(:), latCoord(:)
      real(ESMF_KIND_R4), allocatable :: buf(:,:)
      real(ESMF_KIND_R4), pointer     :: srcPtr(:,:) => null()
      real(ESMF_KIND_R4), pointer     :: dstPtr(:,:) => null()
      type(ESMF_Grid) :: srcGrid
      type(ESMF_Field) :: srcField
      type(ESMF_RouteHandle) :: routeHandle
      logical :: cached
      integer :: exclusiveLBound(2)
      type(ESMF_RegridMethod_Flag) :: regridMethod

      if (present(rc)) rc = ESMF_SUCCESS
      if (present(didRegrid)) didRegrid = .false.

      ! Parse regrid method from name string (default: bilinear)
      regridMethod = ESMF_REGRIDMETHOD_BILINEAR
      if (present(regrid_method_name)) then
         select case (trim(regrid_method_name))
         case ('none', 'NONE')
            ! Should not reach here — caller should skip regridding
            call ESMF_LogWrite("catchem_regrid_field: regrid_method='none' but regrid was called", &
               ESMF_LOGMSG_WARNING)
            regridMethod = ESMF_REGRIDMETHOD_BILINEAR
         case ('bilinear', 'BILINEAR')
            regridMethod = ESMF_REGRIDMETHOD_BILINEAR
         case ('neareststod', 'NEARESTSTOD', 'nearest_stod')
            regridMethod = ESMF_REGRIDMETHOD_NEAREST_STOD
         case ('nearestdtos', 'NEARESTDTOS', 'nearest_dtos')
            regridMethod = ESMF_REGRIDMETHOD_NEAREST_DTOS
         case ('conserve', 'CONSERVE', 'conserve1')
            regridMethod = ESMF_REGRIDMETHOD_CONSERVE
         case ('patch', 'PATCH')
            regridMethod = ESMF_REGRIDMETHOD_PATCH
         case default
            call ESMF_LogWrite("catchem_regrid_field: Unknown regrid method '"// &
               trim(regrid_method_name)//"', using bilinear", ESMF_LOGMSG_WARNING)
            regridMethod = ESMF_REGRIDMETHOD_BILINEAR
         end select
      end if

      ! ---- Open the file (read-only) on every PET ----
      ncStatus = nf90_open(trim(filename), NF90_NOWRITE, ncid)
      if (ncStatus /= NF90_NOERR) then
         call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
            msg="catchem_regrid_field: Cannot open "//trim(filename), &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return
      end if

      ! ---- Read lon/lat coordinate dimensions ----
      call read_coord_dims(ncid, lonname, latname, nlon, nlat, localrc)
      if (localrc /= ESMF_SUCCESS) then
         ncStatus = nf90_close(ncid)
         call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
            msg="catchem_regrid_field: Cannot read coords from "//trim(filename), &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return
      end if

      ! ---- Check if we already have a cached route handle ----
      call cache%lookup(nlon, nlat, idx)
      cached = (idx > 0)

      if (.not. cached) then
         ! Read the actual coordinate values
         allocate(lonCoord(nlon), latCoord(nlat))
         call read_coord_values(ncid, lonname, latname, nlon, nlat, lonCoord, latCoord, localrc)
         if (localrc /= ESMF_SUCCESS) then
            ncStatus = nf90_close(ncid)
            deallocate(lonCoord, latCoord)
            call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
               msg="catchem_regrid_field: Cannot read coord values from "//trim(filename), &
               line=__LINE__, file=__FILE__, rcToReturn=rc)
            return
         end if

         ! Create rectilinear source grid
         call create_src_grid(nlon, nlat, lonCoord, latCoord, srcGrid, localrc)
         deallocate(lonCoord, latCoord)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) then
            ncStatus = nf90_close(ncid)
            return
         end if

         ! Create source field
         srcField = ESMF_FieldCreate(srcGrid, typekind=ESMF_TYPEKIND_R4, &
            staggerloc=ESMF_STAGGERLOC_CENTER, name="regrid_src", rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) then
            ncStatus = nf90_close(ncid)
            return
         end if

         ! Compute regrid weights using the selected method
         call ESMF_FieldRegridStore(srcField, dstField, &
            routehandle=routeHandle, &
            regridmethod=regridMethod, &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
            rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) then
            ncStatus = nf90_close(ncid)
            return
         end if

         ! Store in cache
         call cache%add(nlon, nlat, srcGrid, srcField, routeHandle, idx)

         call ESMF_LogWrite("catchem_regrid_field: Computed regrid weights for "// &
            trim(filename), ESMF_LOGMSG_INFO, rc=localrc)
      else
         srcField    = cache%entries(idx)%srcField
         routeHandle = cache%entries(idx)%routeHandle
      end if

      ! ---- Read data from file into source field ----
      ! Each PET reads its local chunk. Get farrayPtr and exclusiveLBound
      ! together (the ESMF_FieldGet overload requires farrayPtr to expose
      ! exclusiveLBound).
      call ESMF_FieldGet(srcField, localDe=0, farrayPtr=srcPtr, &
         exclusiveLBound=exclusiveLBound, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) then
         ncStatus = nf90_close(ncid)
         return
      end if

      call read_var_to_ptr(ncid, varname, srcPtr, exclusiveLBound, timeSlice, levelSlice, localrc)
      if (localrc /= ESMF_SUCCESS) then
         ncStatus = nf90_close(ncid)
         call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
            msg="catchem_regrid_field: Cannot read var "//trim(varname)// &
                " from "//trim(filename), &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return
      end if

      ncStatus = nf90_close(ncid)

      ! ---- Zero out destination, then apply regridding ----
      call ESMF_FieldGet(dstField, farrayPtr=dstPtr, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      dstPtr = 0.0_ESMF_KIND_R4

      call ESMF_FieldRegrid(srcField, dstField, routeHandle=routeHandle, &
         zeroregion=ESMF_REGION_SELECT, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      if (present(didRegrid)) didRegrid = .true.

   end subroutine catchem_regrid_field

   !--------------------------------------------------------------------------
   !> \brief Look up a cached entry by source grid shape
   !--------------------------------------------------------------------------
   subroutine regrid_cache_lookup(self, nlon, nlat, idx)
      class(RegridCache), intent(in)  :: self
      integer,            intent(in)  :: nlon, nlat
      integer,            intent(out) :: idx
      integer :: i

      idx = 0
      do i = 1, self%count
         if (self%entries(i)%active .and. &
             self%entries(i)%nlon == nlon .and. &
             self%entries(i)%nlat == nlat) then
            idx = i
            return
         end if
      end do
   end subroutine regrid_cache_lookup

   !--------------------------------------------------------------------------
   !> \brief Add a new entry to the regrid cache
   !--------------------------------------------------------------------------
   subroutine regrid_cache_add(self, nlon, nlat, srcGrid, srcField, routeHandle, idx)
      class(RegridCache),      intent(inout) :: self
      integer,                 intent(in)    :: nlon, nlat
      type(ESMF_Grid),        intent(in)    :: srcGrid
      type(ESMF_Field),       intent(in)    :: srcField
      type(ESMF_RouteHandle), intent(in)    :: routeHandle
      integer,                 intent(out)   :: idx

      if (self%count < MAX_REGRID_CACHE) then
         self%count = self%count + 1
         idx = self%count
      else
         ! Overwrite the oldest entry (slot 1) — simple eviction
         call ESMF_RouteHandleDestroy(self%entries(1)%routeHandle)
         call ESMF_FieldDestroy(self%entries(1)%srcField)
         call ESMF_GridDestroy(self%entries(1)%srcGrid)
         idx = 1
      end if

      self%entries(idx)%nlon        = nlon
      self%entries(idx)%nlat        = nlat
      self%entries(idx)%srcGrid     = srcGrid
      self%entries(idx)%srcField    = srcField
      self%entries(idx)%routeHandle = routeHandle
      self%entries(idx)%active      = .true.
   end subroutine regrid_cache_add

   !--------------------------------------------------------------------------
   !> \brief Destroy all cached regrid resources
   !--------------------------------------------------------------------------
   subroutine regrid_cache_cleanup(self, rc)
      class(RegridCache), intent(inout)        :: self
      integer,            intent(out), optional :: rc
      integer :: i, localrc

      if (present(rc)) rc = ESMF_SUCCESS

      do i = 1, self%count
         if (self%entries(i)%active) then
            call ESMF_RouteHandleDestroy(self%entries(i)%routeHandle, rc=localrc)
            call ESMF_FieldDestroy(self%entries(i)%srcField, rc=localrc)
            call ESMF_GridDestroy(self%entries(i)%srcGrid, rc=localrc)
            self%entries(i)%active = .false.
         end if
      end do
      self%count = 0
   end subroutine regrid_cache_cleanup

   !--------------------------------------------------------------------------
   !> \brief Public cleanup wrapper
   !--------------------------------------------------------------------------
   subroutine catchem_regrid_cleanup(cache, rc)
      type(RegridCache), intent(inout)        :: cache
      integer,           intent(out), optional :: rc

      call cache%cleanup(rc)
   end subroutine catchem_regrid_cleanup

   !--------------------------------------------------------------------------
   ! Private helpers
   !--------------------------------------------------------------------------

   !> Read dimension sizes of lon/lat coordinate variables
   subroutine read_coord_dims(ncid, lonname, latname, nlon, nlat, rc)
      integer,          intent(in)  :: ncid
      character(len=*), intent(in)  :: lonname, latname
      integer,          intent(out) :: nlon, nlat, rc

      integer :: ncStatus, varId, ndims
      integer, allocatable :: dimids(:)
      integer :: dimId

      rc = ESMF_SUCCESS

      ! Longitude
      ncStatus = nf90_inq_varid(ncid, trim(lonname), varId)
      if (ncStatus /= NF90_NOERR) then
         rc = ESMF_RC_NOT_FOUND; return
      end if
      ncStatus = nf90_inquire_variable(ncid, varId, ndims=ndims)
      if (ncStatus /= NF90_NOERR .or. ndims < 1) then
         rc = ESMF_RC_NOT_FOUND; return
      end if
      allocate(dimids(ndims))
      ncStatus = nf90_inquire_variable(ncid, varId, dimIds=dimids)
      ncStatus = nf90_inquire_dimension(ncid, dimids(1), len=nlon)
      deallocate(dimids)
      if (ncStatus /= NF90_NOERR) then
         rc = ESMF_RC_NOT_FOUND; return
      end if

      ! Latitude
      ncStatus = nf90_inq_varid(ncid, trim(latname), varId)
      if (ncStatus /= NF90_NOERR) then
         rc = ESMF_RC_NOT_FOUND; return
      end if
      ncStatus = nf90_inquire_variable(ncid, varId, ndims=ndims)
      if (ncStatus /= NF90_NOERR .or. ndims < 1) then
         rc = ESMF_RC_NOT_FOUND; return
      end if
      allocate(dimids(ndims))
      ncStatus = nf90_inquire_variable(ncid, varId, dimIds=dimids)
      ! For a 1-D lat variable, the dim is dimids(1).
      ! For a 2-D lat(lat,lon), the lat dim is also dimids(1) if lat varies along first dim.
      ! We handle the 1-D case here (most common for rectilinear grids).
      ncStatus = nf90_inquire_dimension(ncid, dimids(1), len=nlat)
      deallocate(dimids)
      if (ncStatus /= NF90_NOERR) then
         rc = ESMF_RC_NOT_FOUND; return
      end if

   end subroutine read_coord_dims

   !> Read coordinate values
   subroutine read_coord_values(ncid, lonname, latname, nlon, nlat, &
                                lonCoord, latCoord, rc)
      integer,             intent(in)  :: ncid, nlon, nlat
      character(len=*),    intent(in)  :: lonname, latname
      real(ESMF_KIND_R8),  intent(out) :: lonCoord(nlon), latCoord(nlat)
      integer,             intent(out) :: rc

      integer :: ncStatus, varId

      rc = ESMF_SUCCESS

      ncStatus = nf90_inq_varid(ncid, trim(lonname), varId)
      if (ncStatus /= NF90_NOERR) then; rc = ESMF_RC_NOT_FOUND; return; end if
      ncStatus = nf90_get_var(ncid, varId, lonCoord)
      if (ncStatus /= NF90_NOERR) then; rc = ESMF_RC_NOT_FOUND; return; end if

      ncStatus = nf90_inq_varid(ncid, trim(latname), varId)
      if (ncStatus /= NF90_NOERR) then; rc = ESMF_RC_NOT_FOUND; return; end if
      ncStatus = nf90_get_var(ncid, varId, latCoord)
      if (ncStatus /= NF90_NOERR) then; rc = ESMF_RC_NOT_FOUND; return; end if

   end subroutine read_coord_values

   !> Create a rectilinear ESMF Grid from 1-D lon/lat arrays
   subroutine create_src_grid(nlon, nlat, lonCoord, latCoord, srcGrid, rc)
      integer,            intent(in)  :: nlon, nlat
      real(ESMF_KIND_R8), intent(in)  :: lonCoord(nlon), latCoord(nlat)
      type(ESMF_Grid),    intent(out) :: srcGrid
      integer,            intent(out) :: rc

      integer :: localrc, i, j
      real(ESMF_KIND_R8) :: dlon_half, dlat_half
      real(ESMF_KIND_R8), allocatable :: lonCorner(:), latCorner(:)
      real(ESMF_KIND_R8) :: minCoord(2), maxCoord(2)
      real(ESMF_KIND_R8) :: lonSpan
      logical :: isPeriodic
      integer :: maxIndex(2)
      type(ESMF_VM) :: vm
      integer :: nPets, regDecomp(2), a, b

      rc = ESMF_SUCCESS

      ! Compute corner coordinates (edges of cells)
      allocate(lonCorner(nlon+1), latCorner(nlat+1))

      ! Interior corners at midpoints between centers
      do i = 2, nlon
         lonCorner(i) = 0.5_ESMF_KIND_R8 * (lonCoord(i-1) + lonCoord(i))
      end do
      dlon_half = 0.5_ESMF_KIND_R8 * (lonCoord(2) - lonCoord(1))
      lonCorner(1)      = lonCoord(1)    - dlon_half
      lonCorner(nlon+1) = lonCoord(nlon) + dlon_half

      do j = 2, nlat
         latCorner(j) = 0.5_ESMF_KIND_R8 * (latCoord(j-1) + latCoord(j))
      end do
      dlat_half = 0.5_ESMF_KIND_R8 * (latCoord(2) - latCoord(1))
      latCorner(1)      = latCoord(1)    - dlat_half
      latCorner(nlat+1) = latCoord(nlat) + dlat_half

      minCoord(1) = lonCorner(1)
      minCoord(2) = latCorner(1)
      maxCoord(1) = lonCorner(nlon+1)
      maxCoord(2) = latCorner(nlat+1)
      maxIndex(1) = nlon
      maxIndex(2) = nlat

      deallocate(lonCorner, latCorner)

      ! Check if longitude span is ~360 (global periodic grid)
      ! Float32 coords in NetCDF can cause precision loss, so use tolerance
      lonSpan = maxCoord(1) - minCoord(1)
      isPeriodic = (abs(lonSpan - 360.0_ESMF_KIND_R8) < 1.0_ESMF_KIND_R8)

      ! Compute safe decomposition: each DE must have at least 2 cells
      ! in each dimension (ESMF requirement for bilinear/conservative regridding).
      ! Default ESMF decomposition is (nPets, 1) which can create width-1 DEs
      ! when nlon < 2*nPets.
      call ESMF_VMGetGlobal(vm, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      call ESMF_VMGet(vm, petCount=nPets, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ! Find (a,b) with a*b=nPets, nlon/a >= 2, nlat/b >= 2
      ! Prefer more DEs along lon (a as large as possible)
      regDecomp(1) = 1
      regDecomp(2) = nPets
      do a = 1, nPets
         if (mod(nPets, a) /= 0) cycle
         b = nPets / a
         if (nlon >= 2*a .and. nlat >= 2*b) then
            regDecomp(1) = a
            regDecomp(2) = b
         end if
      end do

      if (isPeriodic) then
         ! Snap to exactly 360 to satisfy ESMF periodicity check
         maxCoord(1) = minCoord(1) + 360.0_ESMF_KIND_R8

         srcGrid = ESMF_GridCreate1PeriDimUfrm( &
            maxIndex=maxIndex, &
            minCornerCoord=minCoord, &
            maxCornerCoord=maxCoord, &
            regDecomp=regDecomp, &
            staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
            rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      else
         ! Regional (non-periodic) grid
         srcGrid = ESMF_GridCreateNoPeriDimUfrm( &
            maxIndex=maxIndex, &
            minCornerCoord=minCoord, &
            maxCornerCoord=maxCoord, &
            regDecomp=regDecomp, &
            staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
            rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      end if

   end subroutine create_src_grid

   !> Read a 2D slice of a variable into a Fortran pointer (local DE portion).
   !! Each PET reads its local chunk based on its exclusive bounds in the
   !! ESMF Field.
   subroutine read_var_to_ptr(ncid, varname, ptr, globalStart, timeSlice, levelSlice, rc)
      integer,             intent(in)    :: ncid
      character(len=*),    intent(in)    :: varname
      real(ESMF_KIND_R4),  intent(inout), pointer :: ptr(:,:)
      integer,             intent(in)    :: globalStart(2)
      integer, optional,   intent(in)    :: timeSlice
      integer, optional,   intent(in)    :: levelSlice
      integer,             intent(out)   :: rc

      integer :: ncStatus, varId, ndims, uid, xtype
      integer :: i1, i2, j1, j2
      integer, allocatable :: dimids(:), start(:), cnt(:)
      character(len=ESMF_MAXSTR) :: dimName
      integer :: timeDimLen
      real(ESMF_KIND_R4), allocatable :: buf(:,:)

      rc = ESMF_SUCCESS

      ncStatus = nf90_inq_varid(ncid, trim(varname), varId)
      if (ncStatus /= NF90_NOERR) then; rc = ESMF_RC_NOT_FOUND; return; end if

      ncStatus = nf90_inquire_variable(ncid, varId, ndims=ndims, xtype=xtype)
      if (ncStatus /= NF90_NOERR) then; rc = ESMF_RC_NOT_FOUND; return; end if

      ncStatus = nf90_inquire(ncid, unlimitedDimId=uid)

      allocate(dimids(ndims))
      ncStatus = nf90_inquire_variable(ncid, varId, dimIds=dimids)

      allocate(start(ndims), cnt(ndims))
      start = 1
      cnt = 1

      ! Identify spatial dims (first two that are not time)
      ! and handle time dimension
      ! Use global exclusive bounds for file read start position
      ! (farrayPtr bounds are 1-based locally, but file positions are global)
      i1 = lbound(ptr, 1)
      i2 = ubound(ptr, 1)
      j1 = lbound(ptr, 2)
      j2 = ubound(ptr, 2)

      start(1) = globalStart(1)
      cnt(1) = i2 - i1 + 1
      start(2) = globalStart(2)
      cnt(2) = j2 - j1 + 1

      ! Handle time dimension (last dim)
      if (ndims >= 3) then
         if (uid /= -1 .and. dimids(ndims) == uid) then
            if (present(timeSlice)) start(ndims) = timeSlice
            cnt(ndims) = 1
         else
            ! Check for named time dimension
            dimName = ''
            ncStatus = nf90_inquire_dimension(ncid, dimids(ndims), name=dimName, len=timeDimLen)
            if (ncStatus == NF90_NOERR .and. &
               (index(dimName,'time') > 0 .or. index(dimName,'Time') > 0 .or. &
                index(dimName,'TIME') > 0 .or. index(dimName,'month') > 0 .or. &
                index(dimName,'Month') > 0)) then
               if (present(timeSlice)) start(ndims) = timeSlice
               cnt(ndims) = 1
            end if
         end if
      end if

      ! Handle vertical level dimension (dim 3 for lon,lat,lev or lon,lat,lev,time)
      ! If levelSlice is provided, select a single level from the 3rd dimension
      if (present(levelSlice) .and. ndims >= 3) then
         ! Determine which dim is the level dim.
         ! For ndims=3 with no time: dim3 is level
         ! For ndims=4 with time as last dim: dim3 is level
         ! For ndims=3 with time as last dim: no level dim available, skip
         if (ndims == 4) then
            ! (lon, lat, lev, time) — dim 3 is level
            start(3) = levelSlice
            cnt(3) = 1
         else if (ndims == 3) then
            ! Check if dim 3 is time; if not, it's the level dim
            dimName = ''
            ncStatus = nf90_inquire_dimension(ncid, dimids(3), name=dimName)
            if (.not. (index(dimName,'time') > 0 .or. index(dimName,'Time') > 0 .or. &
                       index(dimName,'TIME') > 0 .or. index(dimName,'month') > 0 .or. &
                       index(dimName,'Month') > 0 .or. &
                       (uid /= -1 .and. dimids(3) == uid))) then
               ! dim 3 is the level dimension
               start(3) = levelSlice
               cnt(3) = 1
            end if
         end if
      end if

      ! Read into a temporary buffer, then copy to ptr
      allocate(buf(cnt(1), cnt(2)))
      ncStatus = nf90_get_var(ncid, varId, buf, start=start, count=cnt)
      if (ncStatus /= NF90_NOERR) then
         deallocate(buf, start, cnt, dimids)
         rc = ESMF_RC_NOT_FOUND
         return
      end if

      ptr(i1:i2, j1:j2) = buf(:,:)

      deallocate(buf, start, cnt, dimids)

   end subroutine read_var_to_ptr

end module catchem_regrid_mod
