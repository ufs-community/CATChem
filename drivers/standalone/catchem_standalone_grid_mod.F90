!> \file catchem_standalone_grid_mod.F90
!! \brief ESMF grid creation utilities for the CATChem standalone driver
!!
!! \details
!! In a coupled configuration the meteorological driver (e.g. FV3) provides
!! the ESMF grid to CATChem. When running standalone there is no parent
!! component, so the standalone driver must create the grid itself. This
!! module provides reusable routines to build either:
!!
!!  - a single-column grid (1 x 1), used for fast process testing, or
!!  - a regular lat-lon grid (nx x ny), periodic in longitude, used for
!!    gridded runs.
!!
!! Both routines return a fully decomposed ESMF_Grid with center coordinates
!! populated in spherical degrees, suitable for handing to a CATChem cap that
!! supports a standalone (self-providing) initialization path.
!!
!! \author CATChem standalone driver
!! \ingroup catchem_nuopc_group
module catchem_standalone_grid_mod

   use ESMF

   implicit none
   private

   public :: CATChemGridConfig
   public :: create_standalone_grid
   public :: get_grid_num_levels

   character(len=*), parameter :: GRID_MODE_COLUMN  = "column"
   character(len=*), parameter :: GRID_MODE_GRIDDED = "gridded"

   !> ESMF_Info key under which the number of vertical levels is stored on the
   !! grid. The vertical dimension itself is NOT part of the ESMF_Grid geometry
   !! (it is an ungridded dimension on the fields); this key simply lets the
   !! level count travel with the grid object so the cap can realize fields
   !! with the correct number of levels.
   character(len=*), parameter, public :: CATCHEM_GRID_NLEV_KEY = "/catchem/num_levels"

   !> \brief Configuration describing the standalone grid to build
   type :: CATChemGridConfig
      character(len=32) :: mode = GRID_MODE_COLUMN  !< "column" or "gridded"
      integer  :: nx = 1                            !< number of grid cells in x (lon)
      integer  :: ny = 1                            !< number of grid cells in y (lat)
      integer  :: nz = 72                           !< number of vertical levels
      real(ESMF_KIND_R8) :: lon_start = 0.0_ESMF_KIND_R8   !< western edge (deg)
      real(ESMF_KIND_R8) :: lon_end   = 360.0_ESMF_KIND_R8 !< eastern edge (deg)
      real(ESMF_KIND_R8) :: lat_start = -90.0_ESMF_KIND_R8 !< southern edge (deg)
      real(ESMF_KIND_R8) :: lat_end   =  90.0_ESMF_KIND_R8 !< northern edge (deg)
      real(ESMF_KIND_R8) :: column_lon = 0.0_ESMF_KIND_R8  !< single-column longitude (deg)
      real(ESMF_KIND_R8) :: column_lat = 0.0_ESMF_KIND_R8  !< single-column latitude (deg)
      real(ESMF_KIND_R8) :: column_dlon = 0.0_ESMF_KIND_R8 !< nominal column width in lon (deg); 0 disables corners
      real(ESMF_KIND_R8) :: column_dlat = 0.0_ESMF_KIND_R8 !< nominal column width in lat (deg); 0 disables corners
   end type CATChemGridConfig

contains

   !> \brief Create an ESMF grid for the standalone CATChem driver
   !!
   !! Dispatches to the column or gridded builder based on \c cfg%mode.
   !!
   !! @param[in]  cfg  Grid configuration
   !! @param[out] grid Resulting ESMF grid (center stagger coordinates filled)
   !! @param[out] rc   ESMF return code (ESMF_SUCCESS on success)
   subroutine create_standalone_grid(cfg, grid, rc)
      type(CATChemGridConfig), intent(in)  :: cfg
      type(ESMF_Grid),         intent(out) :: grid
      integer,                 intent(out) :: rc

      rc = ESMF_SUCCESS

      select case (trim(cfg%mode))
      case (GRID_MODE_COLUMN)
         call create_column_grid(cfg, grid, rc)
      case (GRID_MODE_GRIDDED)
         call create_gridded_grid(cfg, grid, rc)
      case default
         call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="Unknown grid mode '"//trim(cfg%mode)//"' (expected 'column' or 'gridded')", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return
      end select
      if (rc /= ESMF_SUCCESS) return

      ! Record the vertical level count as grid metadata. The vertical is an
      ! ungridded dimension applied to fields, so it is stored here rather than
      ! baked into the horizontal grid geometry.
      call set_grid_num_levels(grid, cfg%nz, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

   end subroutine create_standalone_grid

   !> \brief Store the number of vertical levels on the grid as ESMF_Info
   subroutine set_grid_num_levels(grid, nz, rc)
      type(ESMF_Grid), intent(inout) :: grid
      integer,         intent(in)    :: nz
      integer,         intent(out)   :: rc

      type(ESMF_Info) :: info

      rc = ESMF_SUCCESS

      call ESMF_InfoGetFromHost(grid, info, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      call ESMF_InfoSet(info, key=CATCHEM_GRID_NLEV_KEY, value=nz, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

   end subroutine set_grid_num_levels

   !> \brief Retrieve the number of vertical levels stored on a grid
   !!
   !! Companion accessor for consumers of the grid (e.g. the cap's standalone
   !! initialization path) to read back the level count that was stamped on by
   !! \c create_standalone_grid.
   !!
   !! @param[in]  grid Grid previously built by create_standalone_grid
   !! @param[out] nz   Number of vertical levels (0 if the key is absent)
   !! @param[out] rc   ESMF return code (ESMF_SUCCESS on success)
   subroutine get_grid_num_levels(grid, nz, rc)
      type(ESMF_Grid), intent(in)  :: grid
      integer,         intent(out) :: nz
      integer,         intent(out) :: rc

      type(ESMF_Info) :: info
      logical :: isPresent

      rc = ESMF_SUCCESS
      nz = 0

      call ESMF_InfoGetFromHost(grid, info, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      isPresent = ESMF_InfoIsPresent(info, key=CATCHEM_GRID_NLEV_KEY, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      if (isPresent) then
         call ESMF_InfoGet(info, key=CATCHEM_GRID_NLEV_KEY, value=nz, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
      end if

   end subroutine get_grid_num_levels

   !> \brief Create a single-column (1 x 1) ESMF grid
   subroutine create_column_grid(cfg, grid, rc)
      type(CATChemGridConfig), intent(in)  :: cfg
      type(ESMF_Grid),         intent(out) :: grid
      integer,                 intent(out) :: rc

      real(ESMF_KIND_R8), pointer :: lonPtr(:,:), latPtr(:,:)
      real(ESMF_KIND_R8), pointer :: lonCorner(:,:), latCorner(:,:)
      integer :: localDECount
      integer :: clbnd(2), cubnd(2)
      integer :: ci, cj
      logical :: have_corners

      rc = ESMF_SUCCESS

      ! A single column is conceptually a point, so it has no intrinsic cell
      ! width. Corner coordinates (needed for conservative regridding and
      ! ESMF_FieldRegridGetArea) can only be defined if the user supplies a
      ! nominal column width via column_dlon/column_dlat. When either is <= 0
      ! we skip corners and warn that conservative regridding is unsupported.
      have_corners = (cfg%column_dlon > 0.0_ESMF_KIND_R8 .and. &
                      cfg%column_dlat > 0.0_ESMF_KIND_R8)
      if (.not. have_corners) then
         call ESMF_LogWrite('create_column_grid: column_dlon/column_dlat not '// &
            'set (or <= 0); no CORNER coordinates added. Conservative regridding '// &
            'and grid-cell area (AREA_M2) are not supported in column mode -- use '// &
            'neareststod/bilinear regridding for emissions.', &
            ESMF_LOGMSG_WARNING, rc=rc)
      end if

      grid = ESMF_GridCreateNoPeriDim( &
         maxIndex=(/1, 1/), &
         coordSys=ESMF_COORDSYS_SPH_DEG, &
         regDecomp=(/1, 1/), &
         indexflag=ESMF_INDEX_GLOBAL, &
         name="catchem_column_grid", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      call ESMF_GridAddCoord(grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      if (have_corners) then
         call ESMF_GridAddCoord(grid, staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
      end if

      ! A 1x1 grid has a single decomposition element (DE) that lives on one
      ! PET only. When the job is launched with more PETs than DEs, every other
      ! PET owns no local DE (localDECount == 0) and must NOT try to retrieve a
      ! coordinate array pointer -- doing so raises "localDeCount <= 0 prohibits
      ! request" and deadlocks the run. Only the owning PET fills coordinates.
      call ESMF_GridGet(grid, localDECount=localDECount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      if (localDECount > 0) then
         call ESMF_GridGetCoord(grid, coordDim=1, staggerloc=ESMF_STAGGERLOC_CENTER, &
            farrayPtr=lonPtr, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         call ESMF_GridGetCoord(grid, coordDim=2, staggerloc=ESMF_STAGGERLOC_CENTER, &
            farrayPtr=latPtr, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         if (associated(lonPtr)) lonPtr = cfg%column_lon
         if (associated(latPtr)) latPtr = cfg%column_lat

         ! Fill the four cell corners around the center using the nominal
         ! column width (only present when have_corners is true).
         if (have_corners) then
            call ESMF_GridGetCoord(grid, coordDim=1, staggerloc=ESMF_STAGGERLOC_CORNER, &
               computationalLBound=clbnd, computationalUBound=cubnd, farrayPtr=lonCorner, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return

            call ESMF_GridGetCoord(grid, coordDim=2, staggerloc=ESMF_STAGGERLOC_CORNER, &
               farrayPtr=latCorner, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return

            ! Corner (1,1) is the SW edge, (2,2) the NE edge of the single cell.
            do cj = clbnd(2), cubnd(2)
               do ci = clbnd(1), cubnd(1)
                  lonCorner(ci, cj) = cfg%column_lon + &
                     (real(ci, ESMF_KIND_R8) - 1.5_ESMF_KIND_R8) * cfg%column_dlon
                  latCorner(ci, cj) = cfg%column_lat + &
                     (real(cj, ESMF_KIND_R8) - 1.5_ESMF_KIND_R8) * cfg%column_dlat
               end do
            end do
         end if
      end if

   end subroutine create_column_grid

   !> \brief Create a regular lat-lon ESMF grid (periodic in longitude)
   !!
   !! Cell centers are placed at the midpoint of each cell. The grid is
   !! decomposed across the PETs of the current VM using an automatic 2D
   !! factorization of the PET count (longitude x latitude), bounded by nx
   !! and ny. PETs left over when petCount cannot be fully used simply own no
   !! local DE and are skipped below.
   subroutine create_gridded_grid(cfg, grid, rc)
      type(CATChemGridConfig), intent(in)  :: cfg
      type(ESMF_Grid),         intent(out) :: grid
      integer,                 intent(out) :: rc

      type(ESMF_VM) :: vm
      real(ESMF_KIND_R8), pointer :: lonPtr(:,:), latPtr(:,:)
      real(ESMF_KIND_R8), pointer :: lonCorner(:,:), latCorner(:,:)
      real(ESMF_KIND_R8) :: dlon, dlat
      integer :: petCount, decompX, decompY
      integer :: i, j
      integer :: lbnd(2), ubnd(2)
      integer :: clbnd(2), cubnd(2)
      integer :: localDECount

      rc = ESMF_SUCCESS

      if (cfg%nx < 1 .or. cfg%ny < 1) then
         call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="Gridded mode requires nx >= 1 and ny >= 1", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return
      end if

      call ESMF_VMGetCurrent(vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      call ESMF_VMGet(vm, petCount=petCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Decompose the grid into a 2D block layout (decompX x decompY) using an
      ! automatic factorization of petCount, bounded by nx and ny. The layout
      ! never asks for more blocks than PETs, so decompX*decompY <= petCount.
      call factor_2d(petCount, cfg%nx, cfg%ny, decompX, decompY)

      grid = ESMF_GridCreate1PeriDim( &
         maxIndex=(/cfg%nx, cfg%ny/), &
         coordSys=ESMF_COORDSYS_SPH_DEG, &
         regDecomp=(/decompX, decompY/), &
         indexflag=ESMF_INDEX_GLOBAL, &
         name="catchem_latlon_grid", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      call ESMF_GridAddCoord(grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Conservative regridding and grid-cell area calculation
      ! (ESMF_FieldRegridGetArea) require coordinates at the CORNER stagger
      ! location, not just cell centers. Add and fill them so emission
      ! remapping and AREA_M2 succeed.
      call ESMF_GridAddCoord(grid, staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      dlon = (cfg%lon_end - cfg%lon_start) / real(cfg%nx, ESMF_KIND_R8)
      dlat = (cfg%lat_end - cfg%lat_start) / real(cfg%ny, ESMF_KIND_R8)

      ! The grid is decomposed into `decompY` blocks in latitude. When the job
      ! runs with more PETs than blocks (petCount > decompY), the surplus PETs
      ! own no local DE (localDECount == 0). Retrieving a coordinate pointer on
      ! such a PET raises "localDeCount <= 0 prohibits request" and hangs the
      ! run, so only PETs that hold a DE fill coordinates.
      call ESMF_GridGet(grid, localDECount=localDECount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      if (localDECount > 0) then
         call ESMF_GridGetCoord(grid, coordDim=1, staggerloc=ESMF_STAGGERLOC_CENTER, &
            computationalLBound=lbnd, computationalUBound=ubnd, farrayPtr=lonPtr, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         call ESMF_GridGetCoord(grid, coordDim=2, staggerloc=ESMF_STAGGERLOC_CENTER, &
            farrayPtr=latPtr, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         ! Fill cell-center coordinates using global indices (ESMF_INDEX_GLOBAL)
         do j = lbnd(2), ubnd(2)
            do i = lbnd(1), ubnd(1)
               lonPtr(i, j) = cfg%lon_start + (real(i, ESMF_KIND_R8) - 0.5_ESMF_KIND_R8) * dlon
               latPtr(i, j) = cfg%lat_start + (real(j, ESMF_KIND_R8) - 0.5_ESMF_KIND_R8) * dlat
            end do
         end do

         ! Fill cell-corner coordinates. Corner (i, j) sits at the lower-left
         ! edge of cell (i, j), so it is offset by a full index from the
         ! domain origin (no half-cell shift).
         call ESMF_GridGetCoord(grid, coordDim=1, staggerloc=ESMF_STAGGERLOC_CORNER, &
            computationalLBound=clbnd, computationalUBound=cubnd, farrayPtr=lonCorner, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         call ESMF_GridGetCoord(grid, coordDim=2, staggerloc=ESMF_STAGGERLOC_CORNER, &
            farrayPtr=latCorner, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         do j = clbnd(2), cubnd(2)
            do i = clbnd(1), cubnd(1)
               lonCorner(i, j) = cfg%lon_start + (real(i, ESMF_KIND_R8) - 1.0_ESMF_KIND_R8) * dlon
               latCorner(i, j) = cfg%lat_start + (real(j, ESMF_KIND_R8) - 1.0_ESMF_KIND_R8) * dlat
            end do
         end do
      end if

   end subroutine create_gridded_grid

   !> \brief Factor petCount into a 2D block layout (px x py) for an nx x ny grid
   !!
   !! Chooses px in [1, nx] and py in [1, ny] so that px*py <= petCount and the
   !! number of used PETs (px*py) is maximised; ties are broken toward a block
   !! aspect ratio close to nx/ny. Any PETs not covered by the layout own no
   !! local DE and are handled gracefully by the localDECount guard.
   subroutine factor_2d(petCount, nx, ny, px, py)
      integer, intent(in)  :: petCount, nx, ny
      integer, intent(out) :: px, py

      integer :: cx, cy, used, best_used
      real    :: target_ratio, ratio, best_score, score

      px = 1
      py = 1
      best_used = 0
      best_score = huge(1.0)
      target_ratio = real(max(nx, 1)) / real(max(ny, 1))

      ! Try every candidate number of longitude blocks and take the largest
      ! matching latitude block count that still fits within petCount and ny.
      do cx = 1, min(petCount, nx)
         cy = min(petCount / cx, ny)
         if (cy < 1) cycle
         used = cx * cy
         ratio = real(cx) / real(cy)
         score = abs(ratio - target_ratio)
         if (used > best_used .or. (used == best_used .and. score < best_score)) then
            best_used = used
            best_score = score
            px = cx
            py = cy
         end if
      end do

      if (px < 1) px = 1
      if (py < 1) py = 1
   end subroutine factor_2d

end module catchem_standalone_grid_mod
