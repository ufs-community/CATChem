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

   implicit none
   private

   public :: latlon_diag_init
   public :: latlon_diag_write_2d
   public :: latlon_diag_write_3d
   public :: latlon_diag_cleanup
   public :: latlon_diag_is_init

   ! Module state — persists across calls
   type(ESMF_Grid),        save :: ll_grid         !< Global lat/lon output grid
   type(ESMF_Field),       save :: ll_src_2d       !< Temp 2D field on model grid
   type(ESMF_Field),       save :: ll_dst_2d       !< Temp 2D field on lat/lon grid
   type(ESMF_RouteHandle), save :: ll_rh           !< Model -> lat/lon regrid handle
   logical,                save :: ll_initialized = .false.
   integer,                save :: ll_nlon = 0     !< Output longitude count
   integer,                save :: ll_nlat = 0     !< Output latitude count
   character(len=512),     save :: ll_current_file = '' !< Track current output file

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
      call ESMF_FieldRegridStore(ll_src_2d, ll_dst_2d, &
         routehandle=ll_rh, &
         regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
         unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
         rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ll_initialized = .true.
      ll_current_file = ''

      call ESMF_LogWrite("latlon_diag_init: Lat/lon diagnostic output ready", &
         ESMF_LOGMSG_INFO, rc=localrc)

   end subroutine latlon_diag_init

   !--------------------------------------------------------------------------
   !> \brief Regrid a 2D field to lat/lon and write to NetCDF
   !!
   !! Copies local data into the cached source field, applies the regrid,
   !! and writes the result using ESMF_FieldWrite.
   !!
   !! \param[in] data_2d    Local 2D data array (R4, on model grid DE)
   !! \param[in] varname    Variable name in output file
   !! \param[in] filename   Output file path
   !! \param[in] timeslice  Time record index
   !! \param[out] rc        Return code
   !--------------------------------------------------------------------------
   subroutine latlon_diag_write_2d(data_2d, varname, filename, timeslice, rc)
      real(ESMF_KIND_R4), intent(in) :: data_2d(:,:)
      character(len=*),   intent(in) :: varname
      character(len=*),   intent(in) :: filename
      integer,            intent(in) :: timeslice
      integer,            intent(out) :: rc

      integer :: localrc
      real(ESMF_KIND_R4), pointer :: srcPtr(:,:), dstPtr(:,:)
      type(ESMF_FileStatus_Flag) :: fstatus

      rc = ESMF_SUCCESS
      if (.not. ll_initialized) then
         call ESMF_LogSetError(ESMF_RC_NOT_SET, &
            msg="latlon_diag_write_2d: Not initialized", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return
      end if

      ! Copy local data into source field
      call ESMF_FieldGet(ll_src_2d, farrayPtr=srcPtr, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      srcPtr(:,:) = data_2d(:,:)

      ! Zero destination and apply regrid
      call ESMF_FieldGet(ll_dst_2d, farrayPtr=dstPtr, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      dstPtr = 0.0_ESMF_KIND_R4

      call ESMF_FieldRegrid(ll_src_2d, ll_dst_2d, routeHandle=ll_rh, &
         zeroregion=ESMF_REGION_SELECT, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ! Determine file status: REPLACE for first write to a new file, OLD after
      if (trim(filename) /= trim(ll_current_file)) then
         fstatus = ESMF_FILESTATUS_REPLACE
         ll_current_file = trim(filename)
      else
         fstatus = ESMF_FILESTATUS_OLD
      end if

      ! Write regridded field to single lat/lon NetCDF file
      call ESMF_FieldWrite(ll_dst_2d, fileName=trim(filename), &
         variableName=trim(varname), &
         overwrite=.true., &
         status=fstatus, &
         timeslice=timeslice, &
         rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

   end subroutine latlon_diag_write_2d

   !--------------------------------------------------------------------------
   !> \brief Regrid a 3D field to lat/lon and write to NetCDF
   !!
   !! Regrids each vertical level independently using the cached 2D
   !! RouteHandle, assembles into a 3D ESMF field on the lat/lon grid,
   !! and writes using ESMF_FieldWrite.
   !!
   !! \param[in] data_3d    Local 3D data array (R4, nx,ny,nlev)
   !! \param[in] varname    Variable name in output file
   !! \param[in] filename   Output file path
   !! \param[in] timeslice  Time record index
   !! \param[out] rc        Return code
   !--------------------------------------------------------------------------
   subroutine latlon_diag_write_3d(data_3d, varname, filename, timeslice, rc)
      real(ESMF_KIND_R4), intent(in) :: data_3d(:,:,:)
      character(len=*),   intent(in) :: varname
      character(len=*),   intent(in) :: filename
      integer,            intent(in) :: timeslice
      integer,            intent(out) :: rc

      integer :: localrc, k, nlev
      real(ESMF_KIND_R4), pointer :: srcPtr(:,:), dstPtr(:,:)
      type(ESMF_Field) :: dst_3d_field
      real(ESMF_KIND_R4), pointer :: dstPtr3d(:,:,:)
      type(ESMF_FileStatus_Flag) :: fstatus

      rc = ESMF_SUCCESS
      if (.not. ll_initialized) then
         call ESMF_LogSetError(ESMF_RC_NOT_SET, &
            msg="latlon_diag_write_3d: Not initialized", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return
      end if

      nlev = size(data_3d, 3)

      ! Create temporary 3D field on lat/lon grid
      dst_3d_field = ESMF_FieldCreate(ll_grid, typekind=ESMF_TYPEKIND_R4, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         ungriddedLBound=(/1/), ungriddedUBound=(/nlev/), &
         name=trim(varname), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call ESMF_FieldGet(dst_3d_field, farrayPtr=dstPtr3d, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) then
         call ESMF_FieldDestroy(dst_3d_field, rc=localrc)
         return
      end if
      dstPtr3d = 0.0_ESMF_KIND_R4

      ! Get pointers to the reusable 2D src/dst fields
      call ESMF_FieldGet(ll_src_2d, farrayPtr=srcPtr, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) then
         call ESMF_FieldDestroy(dst_3d_field, rc=localrc)
         return
      end if
      call ESMF_FieldGet(ll_dst_2d, farrayPtr=dstPtr, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) then
         call ESMF_FieldDestroy(dst_3d_field, rc=localrc)
         return
      end if

      ! Regrid each vertical level as a 2D slab, accumulate into 3D field
      do k = 1, nlev
         srcPtr(:,:) = data_3d(:,:,k)
         dstPtr = 0.0_ESMF_KIND_R4
         call ESMF_FieldRegrid(ll_src_2d, ll_dst_2d, routeHandle=ll_rh, &
            zeroregion=ESMF_REGION_SELECT, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) then
            call ESMF_FieldDestroy(dst_3d_field, rc=localrc)
            return
         end if
         dstPtr3d(:,:,k) = dstPtr(:,:)
      end do

      ! Determine file status
      if (trim(filename) /= trim(ll_current_file)) then
         fstatus = ESMF_FILESTATUS_REPLACE
         ll_current_file = trim(filename)
      else
         fstatus = ESMF_FILESTATUS_OLD
      end if

      ! Write the 3D regridded field
      call ESMF_FieldWrite(dst_3d_field, fileName=trim(filename), &
         variableName=trim(varname), &
         overwrite=.true., &
         status=fstatus, &
         timeslice=timeslice, &
         rc=localrc)

      call ESMF_FieldDestroy(dst_3d_field, rc=localrc)

      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

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

end module catchem_latlon_output_mod
