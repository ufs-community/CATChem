!> \file catchem_latlon_output_mod.F90
!! \brief Regrid cubed-sphere diagnostic fields to a single lat/lon file
!!
!! Provides routines to regrid model fields from the cubed-sphere grid
!! to a regular lat/lon grid and write them to a single NetCDF file
!! (stitching all tiles together). Uses ESMF bilinear regridding with
!! a cached RouteHandle so weights are computed only once.
!!
!! The lat/lon resolution is matched to the cubed-sphere resolution:
!!   C{N} -> nlon = 4*N, nlat = 2*N
!!
!! Usage:
!!   call latlon_diag_init(model_grid, rc)         ! once at startup
!!   call latlon_diag_write_2d(data, varname, ...) ! for each 2D field
!!   call latlon_diag_write_3d(data, varname, ...) ! for each 3D field
!!   call latlon_diag_cleanup(rc)                  ! at finalize

module catchem_latlon_output_mod

   use ESMF
   use netcdf
   use, intrinsic :: ieee_exceptions, only: ieee_set_halting_mode, &
      ieee_get_halting_mode, ieee_invalid, ieee_divide_by_zero, &
      ieee_overflow

   implicit none
   private

   public :: latlon_diag_init
   public :: latlon_diag_write_2d
   public :: latlon_diag_write_3d
   public :: latlon_diag_cleanup
   public :: latlon_diag_is_init
   public :: latlon_diag_set_time

   ! Module state — persists across calls
   type(ESMF_Grid),        save :: ll_grid         !< Global lat/lon output grid
   type(ESMF_Field),       save :: ll_src_2d       !< Temp 2D field on model grid
   type(ESMF_Field),       save :: ll_dst_2d       !< Temp 2D field on lat/lon grid
   type(ESMF_RouteHandle), save :: ll_rh           !< Model -> lat/lon regrid handle
   logical,                save :: ll_initialized = .false.
   integer,                save :: ll_nlon = 0     !< Output longitude count
   integer,                save :: ll_nlat = 0     !< Output latitude count
   character(len=512),     save :: ll_current_file = '' !< Track current output file
   integer(ESMF_KIND_I4),  save :: ll_time_val = 0 !< Current time value (epoch seconds)
   integer,                save :: ll_localPet = -1 !< Cached local PET ID

contains

   !> Check if the lat/lon output system has been initialized
   logical function latlon_diag_is_init()
      latlon_diag_is_init = ll_initialized
   end function latlon_diag_is_init

   !--------------------------------------------------------------------------
   !> \brief Initialize the lat/lon diagnostic output system
   !!
   !! Creates a global regular lat/lon grid matching the cubed-sphere
   !! resolution and computes ESMF bilinear regrid weights from the
   !! model grid to the lat/lon grid. The RouteHandle is cached for reuse.
   !!
   !! \param[in] model_grid  The cubed-sphere model grid (from FV3)
   !! \param[out] rc         Return code
   !--------------------------------------------------------------------------
   subroutine latlon_diag_init(model_grid, rc)
      type(ESMF_Grid), intent(inout) :: model_grid
      integer,         intent(out)   :: rc

      integer :: localrc, tileCount, ntile
      integer :: minIndex(2), maxIndex(2)
      real(ESMF_KIND_R8) :: minCoord(2), maxCoord(2)
      integer :: maxIdx(2)
      character(len=128) :: logmsg
      type(ESMF_VM) :: vm
      logical :: halting_inv, halting_dzero, halting_ovf

      rc = ESMF_SUCCESS
      if (ll_initialized) return

      ! --- Determine cubed-sphere tile size ---
      call ESMF_GridGet(model_grid, tileCount=tileCount, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ! Skip lat/lon stitching for single-tile grids (no stitching needed)
      if (tileCount <= 1) then
         call ESMF_LogWrite('latlon_diag_init: single tile, skipping lat/lon output', &
            ESMF_LOGMSG_INFO, rc=localrc)
         return
      end if

      ! Get tile dimensions from tile 1
      call ESMF_GridGet(model_grid, tile=1, staggerloc=ESMF_STAGGERLOC_CENTER, &
         minIndex=minIndex, maxIndex=maxIndex, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ! Tile side length (e.g. 96 for C96)
      ntile = maxIndex(1) - minIndex(1) + 1

      ! Lat/lon resolution matching cubed-sphere: C{N} -> 4N x 2N
      ll_nlon = 4 * ntile
      ll_nlat = 2 * ntile

      write(logmsg, '(A,I0,A,I0,A,I0)') &
         'latlon_diag_init: C', ntile, ' -> ', ll_nlon, 'x', ll_nlat
      call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_INFO, rc=localrc)

      ! --- Create global regular lat/lon grid (periodic in longitude) ---
      maxIdx(1) = ll_nlon
      maxIdx(2) = ll_nlat
      minCoord(1) = 0.0_ESMF_KIND_R8
      minCoord(2) = -90.0_ESMF_KIND_R8
      maxCoord(1) = 360.0_ESMF_KIND_R8
      maxCoord(2) = 90.0_ESMF_KIND_R8

      ll_grid = ESMF_GridCreate1PeriDimUfrm( &
         maxIndex=maxIdx, &
         minCornerCoord=minCoord, &
         maxCornerCoord=maxCoord, &
         staggerLocList=(/ESMF_STAGGERLOC_CENTER/), &
         rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ! --- Create temporary 2D fields for regridding ---
      ll_src_2d = ESMF_FieldCreate(model_grid, typekind=ESMF_TYPEKIND_R4, &
         staggerloc=ESMF_STAGGERLOC_CENTER, name="ll_diag_src", rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ll_dst_2d = ESMF_FieldCreate(ll_grid, typekind=ESMF_TYPEKIND_R4, &
         staggerloc=ESMF_STAGGERLOC_CENTER, name="ll_diag_dst", rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ! --- Compute bilinear regrid weights (cubed-sphere -> lat/lon) ---
      ! Temporarily disable FPE trapping for ESMF weight computation
      call ieee_get_halting_mode(ieee_invalid, halting_inv)
      call ieee_get_halting_mode(ieee_divide_by_zero, halting_dzero)
      call ieee_get_halting_mode(ieee_overflow, halting_ovf)
      call ieee_set_halting_mode(ieee_invalid, .false.)
      call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
      call ieee_set_halting_mode(ieee_overflow, .false.)

      call ESMF_FieldRegridStore(ll_src_2d, ll_dst_2d, &
         routehandle=ll_rh, &
         regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
         unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
         ignoreDegenerate=.true., &
         rc=localrc)

      call ieee_set_halting_mode(ieee_invalid, halting_inv)
      call ieee_set_halting_mode(ieee_divide_by_zero, halting_dzero)
      call ieee_set_halting_mode(ieee_overflow, halting_ovf)

      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ll_initialized = .true.
      ll_current_file = ''

      ! Cache local PET ID for NetCDF I/O gating
      call ESMF_VMGetCurrent(vm, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      call ESMF_VMGet(vm, localPet=ll_localPet, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call ESMF_LogWrite("latlon_diag_init: Lat/lon diagnostic output ready", &
         ESMF_LOGMSG_INFO, rc=localrc)

   end subroutine latlon_diag_init

   !--------------------------------------------------------------------------
   !> \brief Regrid a 2D ESMF field to lat/lon and write to NetCDF
   !!
   !! Copies data into the cached source field, regrids to lat/lon,
   !! gathers the result to PET 0, and writes via direct NetCDF calls
   !! with shared lon/lat/time dimensions. Units and description are
   !! read from the source field's ESMF Info metadata.
   !!
   !! \param[in] srcField   2D ESMF field on the model grid
   !! \param[in] varname    Variable name in output file
   !! \param[in] filename   Output file path
   !! \param[in] timeslice  Time record index
   !! \param[out] rc        Return code
   !--------------------------------------------------------------------------
   subroutine latlon_diag_write_2d(srcField, varname, filename, timeslice, rc)
      type(ESMF_Field), intent(in)  :: srcField
      character(len=*), intent(in)  :: varname
      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: timeslice
      integer,          intent(out) :: rc

      integer :: localrc, localrc2
      real(ESMF_KIND_R4), pointer :: dataPtr(:,:), srcPtr(:,:), dstPtr(:,:)
      real(ESMF_KIND_R4), allocatable :: globalData(:,:)
      character(len=256) :: units, description
      type(ESMF_Info) :: info

      rc = ESMF_SUCCESS
      if (.not. ll_initialized) return

      ! Ensure output file exists (PET 0 creates on first access)
      if (trim(filename) /= trim(ll_current_file)) then
         if (ll_localPet == 0) call nc_create_file(trim(filename))
         ll_current_file = trim(filename)
      end if

      ! Copy source field data into internal src field
      call ESMF_FieldGet(srcField, farrayPtr=dataPtr, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      call ESMF_FieldGet(ll_src_2d, farrayPtr=srcPtr, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      srcPtr(:,:) = dataPtr(:,:)

      ! Regrid to lat/lon
      call ESMF_FieldGet(ll_dst_2d, farrayPtr=dstPtr, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      dstPtr = 0.0_ESMF_KIND_R4

      call ESMF_FieldRegrid(ll_src_2d, ll_dst_2d, routeHandle=ll_rh, &
         zeroregion=ESMF_REGION_SELECT, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ! Gather regridded data to PET 0
      allocate(globalData(ll_nlon, ll_nlat))
      globalData = 0.0_ESMF_KIND_R4
      call ESMF_FieldGather(ll_dst_2d, globalData, rootPet=0, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) then
         deallocate(globalData); return
      end if

      ! Extract field metadata (units, description)
      units = ''; description = ''
      call ESMF_InfoGetFromHost(srcField, info, rc=localrc)
      if (localrc == ESMF_SUCCESS) then
         localrc2 = ESMF_SUCCESS
         call ESMF_InfoGet(info, key="units", value=units, rc=localrc2)
         if (localrc2 /= ESMF_SUCCESS) units = ''
         localrc2 = ESMF_SUCCESS
         call ESMF_InfoGet(info, key="description", value=description, rc=localrc2)
         if (localrc2 /= ESMF_SUCCESS) description = ''
      end if

      ! PET 0 writes to NetCDF with shared dimensions
      if (ll_localPet == 0) then
         call nc_write_2d(trim(filename), trim(varname), globalData, &
            timeslice, trim(units), trim(description))
      end if

      deallocate(globalData)

   end subroutine latlon_diag_write_2d

   !--------------------------------------------------------------------------
   !> \brief Regrid a 3D ESMF field to lat/lon and write to NetCDF
   !!
   !! Regrids each vertical level independently using the cached 2D
   !! RouteHandle, gathers each level to PET 0, and writes the full
   !! 3D array via direct NetCDF calls with shared dimensions.
   !!
   !! \param[in] srcField   3D ESMF field on the model grid
   !! \param[in] varname    Variable name in output file
   !! \param[in] filename   Output file path
   !! \param[in] timeslice  Time record index
   !! \param[out] rc        Return code
   !--------------------------------------------------------------------------
   subroutine latlon_diag_write_3d(srcField, varname, filename, timeslice, rc)
      type(ESMF_Field), intent(in)  :: srcField
      character(len=*), intent(in)  :: varname
      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: timeslice
      integer,          intent(out) :: rc

      integer :: localrc, localrc2, k, nlev
      real(ESMF_KIND_R4), pointer :: dataPtr3d(:,:,:), srcPtr(:,:), dstPtr(:,:)
      real(ESMF_KIND_R4), allocatable :: globalData3d(:,:,:), gathered2d(:,:)
      character(len=256) :: units, description
      type(ESMF_Info) :: info

      rc = ESMF_SUCCESS
      if (.not. ll_initialized) return

      ! Ensure output file exists
      if (trim(filename) /= trim(ll_current_file)) then
         if (ll_localPet == 0) call nc_create_file(trim(filename))
         ll_current_file = trim(filename)
      end if

      ! Get 3D data from source field
      call ESMF_FieldGet(srcField, farrayPtr=dataPtr3d, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      nlev = size(dataPtr3d, 3)

      ! Get 2D work pointers
      call ESMF_FieldGet(ll_src_2d, farrayPtr=srcPtr, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      call ESMF_FieldGet(ll_dst_2d, farrayPtr=dstPtr, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      allocate(globalData3d(ll_nlon, ll_nlat, nlev))
      allocate(gathered2d(ll_nlon, ll_nlat))
      globalData3d = 0.0_ESMF_KIND_R4

      ! Regrid each vertical level and gather to PET 0
      do k = 1, nlev
         srcPtr(:,:) = dataPtr3d(:,:,k)
         dstPtr = 0.0_ESMF_KIND_R4
         call ESMF_FieldRegrid(ll_src_2d, ll_dst_2d, routeHandle=ll_rh, &
            zeroregion=ESMF_REGION_SELECT, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) then
            deallocate(globalData3d, gathered2d); return
         end if

         gathered2d = 0.0_ESMF_KIND_R4
         call ESMF_FieldGather(ll_dst_2d, gathered2d, rootPet=0, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) then
            deallocate(globalData3d, gathered2d); return
         end if
         if (ll_localPet == 0) globalData3d(:,:,k) = gathered2d(:,:)
      end do
      deallocate(gathered2d)

      ! Extract field metadata
      units = ''; description = ''
      call ESMF_InfoGetFromHost(srcField, info, rc=localrc)
      if (localrc == ESMF_SUCCESS) then
         localrc2 = ESMF_SUCCESS
         call ESMF_InfoGet(info, key="units", value=units, rc=localrc2)
         if (localrc2 /= ESMF_SUCCESS) units = ''
         localrc2 = ESMF_SUCCESS
         call ESMF_InfoGet(info, key="description", value=description, rc=localrc2)
         if (localrc2 /= ESMF_SUCCESS) description = ''
      end if

      ! PET 0 writes to NetCDF
      if (ll_localPet == 0) then
         call nc_write_3d(trim(filename), trim(varname), globalData3d, nlev, &
            timeslice, trim(units), trim(description))
      end if

      deallocate(globalData3d)

   end subroutine latlon_diag_write_3d

   !--------------------------------------------------------------------------
   !> \brief Clean up all lat/lon output resources
   !--------------------------------------------------------------------------
   subroutine latlon_diag_cleanup(rc)
      integer, intent(out), optional :: rc
      integer :: localrc

      if (present(rc)) rc = ESMF_SUCCESS
      if (.not. ll_initialized) return

      call ESMF_RouteHandleDestroy(ll_rh, rc=localrc)
      call ESMF_FieldDestroy(ll_src_2d, rc=localrc)
      call ESMF_FieldDestroy(ll_dst_2d, rc=localrc)
      call ESMF_GridDestroy(ll_grid, rc=localrc)
      ll_initialized = .false.
      ll_current_file = ''

   end subroutine latlon_diag_cleanup

   !--------------------------------------------------------------------------
   !> \brief Store current time value for coordinate output
   !--------------------------------------------------------------------------
   subroutine latlon_diag_set_time(time_seconds)
      integer(ESMF_KIND_I4), intent(in) :: time_seconds
      ll_time_val = time_seconds
   end subroutine latlon_diag_set_time

   !--------------------------------------------------------------------------
   ! Private NetCDF helpers (called only on PET 0)
   !--------------------------------------------------------------------------

   !> Create the lat/lon output file with shared dimensions and coordinate
   !! variables. Called once per output file.
   subroutine nc_create_file(filename)
      character(len=*), intent(in) :: filename

      integer :: ncid, ncStatus, lonDimId, latDimId, timeDimId
      integer :: lonVarId, latVarId, timeVarId, i
      real(ESMF_KIND_R8) :: dlon, dlat
      real(ESMF_KIND_R8), allocatable :: vals(:)

      ncStatus = nf90_create(trim(filename), ior(NF90_CLOBBER, NF90_NETCDF4), ncid)
      if (ncStatus /= NF90_NOERR) return

      ! Shared dimensions
      ncStatus = nf90_def_dim(ncid, 'lon', ll_nlon, lonDimId)
      ncStatus = nf90_def_dim(ncid, 'lat', ll_nlat, latDimId)
      ncStatus = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, timeDimId)

      ! Coordinate variables with CF attributes
      ncStatus = nf90_def_var(ncid, 'lon', NF90_DOUBLE, (/lonDimId/), lonVarId)
      ncStatus = nf90_put_att(ncid, lonVarId, 'units', 'degrees_east')
      ncStatus = nf90_put_att(ncid, lonVarId, 'long_name', 'longitude')
      ncStatus = nf90_put_att(ncid, lonVarId, 'standard_name', 'longitude')

      ncStatus = nf90_def_var(ncid, 'lat', NF90_DOUBLE, (/latDimId/), latVarId)
      ncStatus = nf90_put_att(ncid, latVarId, 'units', 'degrees_north')
      ncStatus = nf90_put_att(ncid, latVarId, 'long_name', 'latitude')
      ncStatus = nf90_put_att(ncid, latVarId, 'standard_name', 'latitude')

      ncStatus = nf90_def_var(ncid, 'time', NF90_INT, (/timeDimId/), timeVarId)
      ncStatus = nf90_put_att(ncid, timeVarId, 'units', &
         'seconds since 1970-01-01 00:00:00')
      ncStatus = nf90_put_att(ncid, timeVarId, 'long_name', 'time')
      ncStatus = nf90_put_att(ncid, timeVarId, 'calendar', 'standard')

      ncStatus = nf90_enddef(ncid)

      ! Write longitude values: centers at (i-0.5)*dlon
      dlon = 360.0_ESMF_KIND_R8 / real(ll_nlon, ESMF_KIND_R8)
      allocate(vals(ll_nlon))
      do i = 1, ll_nlon
         vals(i) = (real(i, ESMF_KIND_R8) - 0.5_ESMF_KIND_R8) * dlon
      end do
      ncStatus = nf90_put_var(ncid, lonVarId, vals)
      deallocate(vals)

      ! Write latitude values: centers at -90 + (j-0.5)*dlat
      dlat = 180.0_ESMF_KIND_R8 / real(ll_nlat, ESMF_KIND_R8)
      allocate(vals(ll_nlat))
      do i = 1, ll_nlat
         vals(i) = -90.0_ESMF_KIND_R8 + &
            (real(i, ESMF_KIND_R8) - 0.5_ESMF_KIND_R8) * dlat
      end do
      ncStatus = nf90_put_var(ncid, latVarId, vals)
      deallocate(vals)

      ncStatus = nf90_close(ncid)

   end subroutine nc_create_file

   !> Define (if needed) and write a 2D variable with shared lon/lat/time dims.
   subroutine nc_write_2d(filename, varname, data, timeslice, units, description)
      character(len=*),     intent(in) :: filename, varname
      real(ESMF_KIND_R4),   intent(in) :: data(:,:)
      integer,              intent(in) :: timeslice
      character(len=*),     intent(in) :: units, description

      integer :: ncid, ncStatus, varId
      integer :: lonDimId, latDimId, timeDimId, timeVarId

      ncStatus = nf90_open(trim(filename), NF90_WRITE, ncid)
      if (ncStatus /= NF90_NOERR) return

      ncStatus = nf90_inq_dimid(ncid, 'lon', lonDimId)
      ncStatus = nf90_inq_dimid(ncid, 'lat', latDimId)
      ncStatus = nf90_inq_dimid(ncid, 'time', timeDimId)

      ! Define variable if it doesn't exist yet
      ncStatus = nf90_inq_varid(ncid, trim(varname), varId)
      if (ncStatus /= NF90_NOERR) then
         ncStatus = nf90_redef(ncid)
         ncStatus = nf90_def_var(ncid, trim(varname), NF90_FLOAT, &
            (/lonDimId, latDimId, timeDimId/), varId)
         if (len_trim(units) > 0) &
            ncStatus = nf90_put_att(ncid, varId, 'units', trim(units))
         if (len_trim(description) > 0) &
            ncStatus = nf90_put_att(ncid, varId, 'long_name', trim(description))
         ncStatus = nf90_enddef(ncid)
      end if

      ! Write data slab
      ncStatus = nf90_put_var(ncid, varId, data, &
         start=(/1, 1, timeslice/), count=(/ll_nlon, ll_nlat, 1/))

      ! Write time value
      ncStatus = nf90_inq_varid(ncid, 'time', timeVarId)
      if (ncStatus == NF90_NOERR) then
         ncStatus = nf90_put_var(ncid, timeVarId, ll_time_val, &
            start=(/timeslice/))
      end if

      ncStatus = nf90_close(ncid)

   end subroutine nc_write_2d

   !> Define (if needed) and write a 3D variable with shared lon/lat/lev/time dims.
   subroutine nc_write_3d(filename, varname, data, nlev, timeslice, units, description)
      character(len=*),     intent(in) :: filename, varname
      real(ESMF_KIND_R4),   intent(in) :: data(:,:,:)
      integer,              intent(in) :: nlev, timeslice
      character(len=*),     intent(in) :: units, description

      integer :: ncid, ncStatus, varId
      integer :: lonDimId, latDimId, levDimId, timeDimId, timeVarId
      logical :: need_redef

      ncStatus = nf90_open(trim(filename), NF90_WRITE, ncid)
      if (ncStatus /= NF90_NOERR) return

      ncStatus = nf90_inq_dimid(ncid, 'lon', lonDimId)
      ncStatus = nf90_inq_dimid(ncid, 'lat', latDimId)
      ncStatus = nf90_inq_dimid(ncid, 'time', timeDimId)

      ! Check if lev dim and/or variable need to be created
      need_redef = .false.
      if (nf90_inq_dimid(ncid, 'lev', levDimId) /= NF90_NOERR) need_redef = .true.
      if (nf90_inq_varid(ncid, trim(varname), varId) /= NF90_NOERR) need_redef = .true.

      if (need_redef) then
         ncStatus = nf90_redef(ncid)
         if (nf90_inq_dimid(ncid, 'lev', levDimId) /= NF90_NOERR) then
            ncStatus = nf90_def_dim(ncid, 'lev', nlev, levDimId)
         end if
         if (nf90_inq_varid(ncid, trim(varname), varId) /= NF90_NOERR) then
            ncStatus = nf90_def_var(ncid, trim(varname), NF90_FLOAT, &
               (/lonDimId, latDimId, levDimId, timeDimId/), varId)
            if (len_trim(units) > 0) &
               ncStatus = nf90_put_att(ncid, varId, 'units', trim(units))
            if (len_trim(description) > 0) &
               ncStatus = nf90_put_att(ncid, varId, 'long_name', trim(description))
         end if
         ncStatus = nf90_enddef(ncid)
      end if

      ! Write data slab
      ncStatus = nf90_put_var(ncid, varId, data, &
         start=(/1, 1, 1, timeslice/), count=(/ll_nlon, ll_nlat, nlev, 1/))

      ! Write time value
      ncStatus = nf90_inq_varid(ncid, 'time', timeVarId)
      if (ncStatus == NF90_NOERR) then
         ncStatus = nf90_put_var(ncid, timeVarId, ll_time_val, &
            start=(/timeslice/))
      end if

      ncStatus = nf90_close(ncid)

   end subroutine nc_write_3d

end module catchem_latlon_output_mod
