#define HAVE_NETCDF 1
module AQMIO

   use ESMF
#if HAVE_NETCDF
   use netcdf
#endif

   implicit none

   type AQMIOLayout
      logical :: localIOflag
      integer :: tile
      integer :: ncid
      integer :: iounit
      type(ESMF_GridComp) :: taskComp
   end type AQMIOLayout

   type ioData
      type(AQMIOLayout), pointer :: IOLayout(:) => null()
   end type ioData

   type ioWrapper
      type(ioData), pointer :: IO => null()
   end type ioWrapper

   integer, parameter :: AQMIO_FMT_BIN    = 101, &
      AQMIO_FMT_NETCDF = 102

   private

   public :: AQMIO_FMT_BIN
   public :: AQMIO_FMT_NETCDF

   public :: AQMIO_Create
   public :: AQMIO_Destroy
   public :: AQMIO_FileCreate
   public :: AQMIO_Open
   public :: AQMIO_Close
   public :: AQMIO_Read
   public :: AQMIO_ReadTimes
   public :: AQMIO_DataRead
   public :: AQMIO_Sync
   public :: AQMIO_Write

   ! Enhanced direct data I/O functions (no ESMF fields required)
   public :: AQMIO_Write1D
   public :: AQMIO_Read1D

contains

!------------------------------------------------------------------------------

   function AQMIO_Create(grid, vm, allpes, rc)
      type(ESMF_Grid), intent(in)            :: grid
      type(ESMF_VM),   intent(in),  optional :: vm
      logical,         intent(in),  optional :: allpes
      integer,         intent(out), optional :: rc

      type(ESMF_GridComp) :: AQMIO_Create

      ! -- local variables
      integer             :: localrc
      integer             :: i, iope, localDe, localDeCount, localpe, peCount, npe
      integer             :: tile, tileCount
      integer, dimension(:), allocatable :: localTile, tileToPet, pes, recvpes
      type(ESMF_GridComp) :: IOComp, taskComp
      type(ESMF_VM)       :: localVM, tasksVM
      type(ioWrapper)     :: is
      type(ioData), pointer :: IO => null()

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

      nullify(IO)

      call ESMF_GridGet(grid, localDeCount=localDeCount, tileCount=tileCount, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out


      if (present(vm)) then
         localVM = vm
      else
         call ESMF_VMGetCurrent(localVM, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
      end if

      call ESMF_VMGet(localVM, localPet=localpe, petCount=peCount, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out


      allocate(recvpes(peCount), pes(peCount), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, &
         msg="Unable to allocate internal memory for AQMIO initialization", &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      pes = 0
      pes(localpe+1) = -localDeCount

      call ESMF_VMAllReduce(localVM, pes, recvpes, peCount, ESMF_REDUCE_SUM, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      pes = -1
      npe = 0
      do i = 1, peCount
         if (recvpes(i) < 0) then
            npe = npe + 1
            pes(npe) = i - 1
         end if
      end do

      ! -- create IO component on this PET
      IOComp = ESMF_GridCompCreate(name="io_comp", petList=pes(1:npe), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      deallocate(recvpes, pes, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      call ESMF_GridCompSetServices(IOComp, IOCompSetServices, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (ESMF_GridCompIsPetLocal(IOComp)) then

         allocate(IO, stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, &
            msg="Unable to allocate internal memory for AQMIO initialization", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         allocate(IO % IOLayout(0:localDeCount-1), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, &
            msg="Unable to allocate internal memory for AQMIO initialization", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         is % IO => IO

      else

         is % IO => null()

      end if

      ! -- set internal state for IO component
      call ESMF_GridCompSetInternalState(IOComp, is, localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      ! -- save grid object in IO component
      call ESMF_GridCompSet(IOComp, grid=grid, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      allocate(localTile(tileCount), tileToPet(tileCount*peCount), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, &
         msg="Unable to allocate internal memory for AQMIO initialization", &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      ! -- store which tiles are assigned to this PET
      localTile = -1
      do localDe = 0, localDeCount-1
         call ESMF_GridGet(grid, localDE=localDe, tile=tile, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         localTile(tile) = localpe
         is % IO % IOLayout(localDe) % tile   = tile
         is % IO % IOLayout(localDe) % ncid   = 0
         is % IO % IOLayout(localDe) % iounit = 0
      end do

      tileToPet = -1
      call ESMF_VMAllGather(localVM, localTile, tileToPet, tileCount, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      deallocate(localTile, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      ! -- extract the list of PETs assigned to each tile and create MPI groups
      allocate(pes(peCount), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, &
         msg="Unable to allocate internal memory for AQMIO initialization", &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      ! -- gather PET list for each tile and create tile-specific VMs
      pes = -1
      do tile = 1, tileCount
         npe = 0
         do i = tile, tileCount*peCount, tileCount
            if (tileToPet(i) > -1) then
               npe = npe + 1
               pes(npe) = tileToPet(i)
            end if
         end do

         taskComp = ESMF_GridCompCreate(petList=pes(1:npe), rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         call ESMF_GridCompSetServices(taskComp, IOCompSetServices, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         do localDe = 0, localDeCount-1
            call ESMF_GridGet(grid, localDE=localDe, tile=i, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
            if (tile == i) then
               ! -- create new VM for tile
               is % IO % IOLayout(localDe) % taskComp = taskComp
            end if
         end do
      end do

      deallocate(pes, tileToPet, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      iope = 0
      if (present(allpes)) then
         if (allpes) iope = localpe
      end if


      ! -- flag PET if local I/O must be performed
      do localDe = 0, localDeCount - 1
         call ESMF_GridCompGet(is % IO % IOLayout(localDe) % taskComp, vm=tasksVM, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         call ESMF_VMGet(tasksVM, localPet=localpe, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         is % IO % IOLayout(localDe) % localIOflag = (localpe == iope)
      end do

      AQMIO_Create = IOComp

   end function AQMIO_Create

!------------------------------------------------------------------------------

   subroutine AQMIO_Destroy(IOComp, rc)
      type(ESMF_GridComp)            :: IOComp
      integer, intent(out), optional :: rc

      ! -- local variables
      integer :: localrc
      integer :: localDe
      type(ioWrapper) :: is

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

      if (ESMF_GridCompIsCreated(IOComp)) then
         call ESMF_GridCompGetInternalState(IOComp, is, localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         call AQMIO_Close(IOComp, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         if (associated(is % IO)) then
            if (associated(is % IO % IOLayout)) then
               do localDe = 0, size(is % IO % IOLayout) - 1
                  if (ESMF_GridCompIsCreated(is % IO % IOLayout(localDe) % taskComp)) then
                     call ESMF_GridCompDestroy(is % IO % IOLayout(localDe) % taskComp, rc=localrc)
                     if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                        line=__LINE__, &
                        file=__FILE__, &
                        rcToReturn=rc)) return  ! bail out
                  end if
               end do

               deallocate(is % IO % IOLayout, stat=localrc)
               if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
               nullify(is % IO % IOLayout)

               call ESMF_GridCompDestroy(IOComp, rc=localrc)
               if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
            end if
            nullify(is % IO)
         end if
      end if

   end subroutine AQMIO_Destroy

!------------------------------------------------------------------------------

   subroutine AQMIO_Open(IOComp, fileName, filePath, iomode, iofmt, rc)
      type(ESMF_GridComp),   intent(inout)         :: IOComp
      character(len=*),      intent(in)            :: fileName
      character(len=*),      intent(in),  optional :: filePath
      character(len=*),      intent(in),  optional :: iomode
      integer,               intent(in),  optional :: iofmt
      integer,               intent(out), optional :: rc

      ! -- local variables
      integer :: localrc
      integer :: ncStatus
      integer :: item, localDe, localDeCount, tileCount
      integer :: liofmt
      integer :: cmode
      logical :: create
      character(len=ESMF_MAXPATHLEN) :: fullName
      character(len=6) :: liomode, fmode
      type(ioWrapper) :: is
      type(ESMF_Grid) :: grid

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS


      if (.not.ESMF_GridCompIsPetLocal(IOComp)) then
         return
      end if

      call ESMF_GridCompGetInternalState(IOComp, is, localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (.not.associated(is % IO)) return
      if (.not.associated(is % IO % IOLayout)) return

      localDeCount = size(is % IO % IOLayout)


      liofmt = AQMIO_FMT_NETCDF
      if (present(iofmt)) liofmt = iofmt

      liomode = "read"
      if (present(iomode)) liomode = iomode

      create = .false.
      select case (trim(liomode))
       case ("r", "read")
         cmode = NF90_NOWRITE
         fmode = "read"
       case ("w", "write")
         cmode = NF90_WRITE
         fmode = "write"
       case ("c", "create")
         cmode = NF90_CLOBBER
         create = .true.
       case default
         call ESMF_LogSetError(ESMF_RC_ARG_WRONG, &
            msg="- Unsupported open mode: "//trim(liomode), &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return
      end select

      call ESMF_GridCompGet(IOComp, grid=grid, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      call ESMF_GridGet(grid, tileCount=tileCount, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      do localDe = 0, localDeCount - 1
         if (is % IO % IOLayout(localDe) % localIOflag) then
            if (tileCount > 1) then
               call AQMIO_FileNameGet(fullName, fileName, filePath=filePath, &
                  tile=is % IO % IOLayout(localDe) % tile)
            else
               call AQMIO_FileNameGet(fullName, fileName, filePath=filePath)
            end if
            if      (liofmt == AQMIO_FMT_NETCDF) then
#if HAVE_NETCDF
               if (create) then
                  ncStatus = nf90_create(trim(fullName), cmode, &
                     is % IO % IOLayout(localDe) % ncid)
                  if (ncStatus /= NF90_NOERR) then
                     call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                        msg="Error creating NetCDF data set: "//trim(fullName)//": "//nf90_strerror(ncStatus), &
                        line=__LINE__, &
                        file=__FILE__, &
                        rcToReturn=rc)
                     return  ! bail out
                  end if
                  ncStatus = nf90_enddef(is % IO % IOLayout(localDe) % ncid)
                  if (ncStatus /= NF90_NOERR) then
                     call ESMF_LogSetError(ESMF_RC_FILE_WRITE, &
                        msg="Error ending define mode: "//nf90_strerror(ncStatus), &
                        line=__LINE__, &
                        file=__FILE__, &
                        rcToReturn=rc)
                     return  ! bail out
                  end if
               else
                  ncStatus = nf90_open(trim(fullName), cmode, &
                     is % IO % IOLayout(localDe) % ncid)
                  if (ncStatus /= NF90_NOERR) then
                     call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                        msg="Error opening NetCDF data set: "//trim(fullName)//": "//nf90_strerror(ncStatus), &
                        line=__LINE__, &
                        file=__FILE__, &
                        rcToReturn=rc)
                     return  ! bail out
                  end if
               end if
#else
               call ESMF_LogSetError(rcToCheck=ESMF_RC_LIB_NOT_PRESENT, &
                  msg="- AQMIO was not built with NetCDF support", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return
#endif
            else if (liofmt == AQMIO_FMT_BIN) then

               call ESMF_UtilIOUnitGet (unit=is % IO % IOLayout(localDe) % iounit, rc=localrc)
               if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out

               open(is % IO % IOLayout(localDe) % iounit, file=trim(fullName), &
                  form='unformatted', action=fmode, position='rewind', iostat=localrc)
               if (localrc /= 0) then
                  call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                     msg=" - file: "//trim(fullName), &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return
               end if

            else
               call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
                  msg="I/O format not implemented", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return  ! bail out
            end if
         end if
      end do

   end subroutine AQMIO_Open

!------------------------------------------------------------------------------

   subroutine AQMIO_Close(IOComp, rc)
      type(ESMF_GridComp),   intent(inout)         :: IOComp
      integer,               intent(out), optional :: rc

      ! -- local variables
      integer :: localrc
      integer :: ncStatus
      integer :: item, localDe, localDeCount
      type(ioWrapper) :: is

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

      if (.not.ESMF_GridCompIsPetLocal(IOComp)) return

      call ESMF_GridCompGetInternalState(IOComp, is, localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (.not.associated(is % IO)) return
      if (.not.associated(is % IO % IOLayout)) return

      localDeCount = size(is % IO % IOLayout)

      do localDe = 0, localDeCount - 1
         if (is % IO % IOLayout(localDe) % localIOflag) then
            if (is % IO % IOLayout(localDe) % ncid > 0) then
#if HAVE_NETCDF
               ncStatus = nf90_close(is % IO % IOLayout(localDe) % ncid)
               if (ncStatus /= NF90_NOERR) then
                  call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                     msg="NetCDF error", &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return
               end if
               is % IO % IOLayout(localDe) % ncid = 0
#else
               call ESMF_LogSetError(rcToCheck=ESMF_RC_LIB_NOT_PRESENT, &
                  msg="AQMIO was not built with NetCDF support", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return
#endif
            end if

            if (is % IO % IOLayout(localDe) % iounit > 0) then

               close(is % IO % IOLayout(localDe) % iounit, iostat=localrc)
               if (localrc /= 0) then
                  call ESMF_LogSetError(ESMF_RC_FILE_CLOSE, &
                     msg="Error closing binary data set", &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return
               end if
               is % IO % IOLayout(localDe) % iounit = 0
            end if
         end if
      end do

   end subroutine AQMIO_Close

!------------------------------------------------------------------------------

   logical function AQMIO_IsOpen(IOComp, fileName, filePath, iofmt, rc)
      type(ESMF_GridComp),   intent(inout)         :: IOComp
      character(len=*),      intent(in)            :: fileName
      character(len=*),      intent(in),  optional :: filePath
      integer,               intent(in),  optional :: iofmt
      integer,               intent(out), optional :: rc

      ! -- local variables
      integer :: localrc
      integer :: ncStatus
      integer :: item, localDe, localDeCount, tileCount, pathLen
      integer :: liofmt
      logical :: isFileOpen
      character(len=ESMF_MAXPATHLEN) :: fullName, pathIn
      type(ioWrapper) :: is
      type(ESMF_Grid) :: grid

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

      AQMIO_IsOpen = .false.

      if (.not.ESMF_GridCompIsPetLocal(IOComp)) return

      call ESMF_GridCompGetInternalState(IOComp, is, localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (.not.associated(is % IO)) return
      if (.not.associated(is % IO % IOLayout)) return

      localDeCount = size(is % IO % IOLayout)

      liofmt = AQMIO_FMT_NETCDF
      if (present(iofmt)) liofmt = iofmt

      call ESMF_GridCompGet(IOComp, grid=grid, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      call ESMF_GridGet(grid, tileCount=tileCount, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      isFileOpen = .true.
      do localDe = 0, localDeCount - 1
         if (is % IO % IOLayout(localDe) % localIOflag) then
            if (tileCount > 1) then
               call AQMIO_FileNameGet(fullName, fileName, filePath=filePath, &
                  tile=is % IO % IOLayout(localDe) % tile)
            else
               call AQMIO_FileNameGet(fullName, fileName, filePath=filePath)
            end if
            if      (liofmt == AQMIO_FMT_NETCDF) then
#if HAVE_NETCDF
               if (is % IO % IOLayout(localDe) % ncid > 0) then
                  ncStatus = nf90_inq_path(is % IO % IOLayout(localDe) % ncid, &
                     pathLen, pathIn)
                  if (ncStatus /= NF90_NOERR) then
                     call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                        msg="NetCDF error", &
                        line=__LINE__, &
                        file=__FILE__, &
                        rcToReturn=rc)
                     return
                  end if
                  isFileOpen = isFileOpen .and. (trim(pathIn) == trim(fullName))
               else
                  exit
               end if
#else
               call ESMF_LogSetError(rcToCheck=ESMF_RC_LIB_NOT_PRESENT, &
                  msg="- AQMIO was not built with NetCDF support", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return
#endif
            end if
         end if
      end do

      AQMIO_IsOpen = isFileOpen

   end function AQMIO_IsOpen

!------------------------------------------------------------------------------

   subroutine AQMIO_Sync(IOComp, rc)
      type(ESMF_GridComp),   intent(inout)         :: IOComp
      integer,               intent(out), optional :: rc

      ! -- local variables
      integer :: localrc
      integer :: ncStatus
      integer :: item, localDe, localDeCount
      type(ioWrapper) :: is

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

#if HAVE_NETCDF
      if (.not.ESMF_GridCompIsPetLocal(IOComp)) return

      call ESMF_GridCompGetInternalState(IOComp, is, localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (.not.associated(is % IO)) return
      if (.not.associated(is % IO % IOLayout)) return

      localDeCount = size(is % IO % IOLayout)

      do localDe = 0, localDeCount - 1
         if (is % IO % IOLayout(localDe) % localIOflag) then
            if (is % IO % IOLayout(localDe) % ncid > 0) then
               ncStatus = nf90_sync(is % IO % IOLayout(localDe) % ncid)
               if (ncStatus /= NF90_NOERR) then
                  call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                     msg="NetCDF error", &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return
               end if
            end if
         end if
      end do
#endif

   end subroutine AQMIO_Sync

!------------------------------------------------------------------------------

   subroutine AQMIO_Write(IOComp, fieldList, fieldNameList, timeSlice, &
      fileName, filePath, iofmt, rc)
      type(ESMF_GridComp),   intent(inout)         :: IOComp
      type(ESMF_Field),      intent(in)            :: fieldList(:)
      character(len=*),      intent(in),  optional :: fieldNameList(:)
      integer,               intent(in),  optional :: timeSlice
      character(len=*),      intent(in),  optional :: fileName
      character(len=*),      intent(in),  optional :: filePath
      integer,               intent(in),  optional :: iofmt
      integer,               intent(out), optional :: rc

      ! -- local variables
      integer :: localrc
      integer :: item, localDe, localDeCount
      integer :: liofmt
      type(ioWrapper) :: is

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS


      if (.not.ESMF_GridCompIsPetLocal(IOComp)) then
         return
      end if

      liofmt = AQMIO_FMT_NETCDF
      if (present(iofmt)) liofmt = iofmt

      if (present(fieldNameList)) then
         if (size(fieldNameList) < size(fieldList)) then
            call ESMF_LogSetError(ESMF_RC_ARG_SIZE, &
               msg="size of fieldNameList must equal or larger than "// &
               "size of fieldList", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
         end if
      end if

      call ESMF_GridCompGetInternalState(IOComp, is, localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (.not.associated(is % IO)) return
      if (.not.associated(is % IO % IOLayout)) return

      localDeCount = size(is % IO % IOLayout)

      if (present(fileName)) then
         call AQMIO_Close(IOComp, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         ! Try to open in write mode first, fallback to create mode if file doesn't exist
         call AQMIO_Open(IOComp, fileName, filePath=filePath, iomode="write",&
            iofmt=iofmt, rc=localrc)
         if (localrc /= ESMF_SUCCESS) then
            ! If write mode fails, try create mode for new files
            call AQMIO_Open(IOComp, fileName, filePath=filePath, iomode="create",&
               iofmt=iofmt, rc=localrc)
         end if
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
      end if

      if (present(fieldNameList)) then
         do item = 1, size(fieldList)
            call AQMIO_FieldAccess(IOComp, fieldList(item), "write", &
               variableName=fieldNameList(item), timeSlice=timeSlice, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
         end do
      else
         do item = 1, size(fieldList)
            call AQMIO_FieldAccess(IOComp, fieldList(item), "write", &
               timeSlice=timeSlice, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
         end do
      end if

      if (present(fileName)) then
         call AQMIO_Close(IOComp, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
      end if

   end subroutine AQMIO_Write

!------------------------------------------------------------------------------

   subroutine AQMIO_Read(IOComp, fieldList, fieldNameList, timeSlice, &
      fileName, filePath, iofmt, rc)
      type(ESMF_GridComp),   intent(inout)         :: IOComp
      type(ESMF_Field),      intent(in)            :: fieldList(:)
      character(len=*),      intent(in),  optional :: fieldNameList(:)
      integer,               intent(in),  optional :: timeSlice
      character(len=*),      intent(in),  optional :: fileName
      character(len=*),      intent(in),  optional :: filePath
      integer,               intent(in),  optional :: iofmt
      integer,               intent(out), optional :: rc

      ! -- local variables
      integer :: localrc
      integer :: item, localDe, localDeCount
      logical :: isOpen
      type(ioWrapper) :: is

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

      if (.not.ESMF_GridCompIsPetLocal(IOComp)) return

      if (present(fieldNameList)) then
         if (size(fieldNameList) < size(fieldList)) then
            call ESMF_LogSetError(ESMF_RC_ARG_SIZE, &
               msg="size of fieldNameList must equal or larger than "// &
               "size of fieldList", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
         end if
      end if

      call ESMF_GridCompGetInternalState(IOComp, is, localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (.not.associated(is % IO)) return
      if (.not.associated(is % IO % IOLayout)) return

      localDeCount = size(is % IO % IOLayout)

      if (present(fileName)) then
         isOpen = AQMIO_IsOpen(IOComp, fileName, filePath=filePath, &
            iofmt=iofmt, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         if (.not.isOpen) then
            call AQMIO_Close(IOComp, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
            call AQMIO_Open(IOComp, fileName, filePath=filePath, &
               iomode="read", iofmt=iofmt, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
         end if
      end if

      if (present(fieldNameList)) then
         do item = 1, size(fieldList)
            call AQMIO_FieldAccess(IOComp, fieldList(item), "read", &
               variableName=fieldNameList(item), timeSlice=timeSlice, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
         end do
      else
         do item = 1, size(fieldList)
            call AQMIO_FieldAccess(IOComp, fieldList(item), "read", &
               timeSlice=timeSlice, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
         end do
      end if

      if (present(fileName)) then
         call AQMIO_Close(IOComp, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
      end if

   end subroutine AQMIO_Read

!------------------------------------------------------------------------------

   subroutine AQMIO_ReadTimes(IOComp, variableName, timesList, fileName, filePath, iofmt, rc)
      type(ESMF_GridComp),   intent(inout)           :: IOComp
      character(len=*),      intent(in)              :: variableName
      type(ESMF_Time),       intent(inout), pointer  :: timesList(:)
      character(len=*),      intent(in),    optional :: fileName
      character(len=*),      intent(in),    optional :: filePath
      integer,               intent(in),    optional :: iofmt
      integer,               intent(out),   optional :: rc

      ! -- local variables
      integer :: localrc
      integer :: item, localDe, localDeCount
      logical :: isOpen
      type(ioWrapper) :: is

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

      if (.not.ESMF_GridCompIsPetLocal(IOComp)) return

      call ESMF_GridCompGetInternalState(IOComp, is, localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (.not.associated(is % IO)) return
      if (.not.associated(is % IO % IOLayout)) return

      if (present(iofmt)) then
         if (.not.(iofmt == AQMIO_FMT_NETCDF)) then
            call ESMF_LogSetError(ESMF_RC_ARG_INCOMP, &
               msg="This function only supports NetCDF I/O", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return  ! bail out
         end if
      end if

      localDeCount = size(is % IO % IOLayout)

      if (present(fileName)) then
         isOpen = AQMIO_IsOpen(IOComp, fileName, filePath=filePath, &
            iofmt=iofmt, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         if (.not.isOpen) then
            call AQMIO_Close(IOComp, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
            call AQMIO_Open(IOComp, fileName, filePath=filePath, &
               iomode="read", iofmt=iofmt, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
         end if
      end if

      ! -- read times on localDe = 0 only, assuming tile-specific files are
      ! -- consistent
      call AQMIO_TimesRead(is % IO, variableName, timesList, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (present(fileName)) then
         call AQMIO_Close(IOComp, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
      end if

   end subroutine AQMIO_ReadTimes

!------------------------------------------------------------------------------
! Private methods below
!------------------------------------------------------------------------------

   subroutine AQMIO_FieldAccess(IOComp, field, action, variableName, timeSlice, rc)
      type(ESMF_GridComp),   intent(in)            :: IOComp
      type(ESMF_Field),      intent(in)            :: field
      character(len=*),      intent(in)            :: action
      character(len=*),      intent(in),  optional :: variableName
      integer,               intent(in),  optional :: timeSlice
      integer,               intent(out), optional :: rc

      ! -- local variables
      integer :: localrc
      integer :: localDe, localDeCount, rank
      integer :: de, deCount, dimCount, tile, tileCount, ungriddedCount
      integer :: iofmt
      integer, dimension(:),   pointer     :: ungriddedLBound, ungriddedUBound
      integer, dimension(:),   allocatable :: deToTileMap, localDeToDeMap
      integer, dimension(:,:), allocatable :: minIndexPDe, maxIndexPDe
      integer, dimension(:,:), allocatable :: minIndexPTile, maxIndexPTile
      type(ioWrapper) :: is
      type(ESMF_Grid) :: grid, iogrid
      type(ESMF_DistGrid) :: distgrid
      type(ESMF_VM) :: vm
      type(ESMF_GeomType_flag)      :: geomtype
      type(ESMF_StaggerLoc)         :: staggerloc
      type(ESMF_TypeKind_Flag)      :: typekind

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS


      if (.not.ESMF_GridCompIsPetLocal(IOComp)) then
         return
      end if

      call ESMF_GridCompGet(IOComp, grid=iogrid, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      call ESMF_FieldGet(field, geomtype=geomtype, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (geomtype == ESMF_GEOMTYPE_GRID) then
         call ESMF_FieldGet(field, grid=grid, rank=rank, &
            staggerloc=staggerloc, typekind=typekind, localDeCount=localDeCount, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         ! -- check if field is built on I/O component grid
         if (grid /= iogrid) then
            call ESMF_LogWrite("field and I/O component may not be on same grid", &
               ESMF_LOGMSG_WARNING, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
         else
         end if

         ! -- check
         if (rank < 2 .or. rank > 3) then
            call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
               msg="Only 2D and 3D fields are supported.", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return  ! bail out
         end if

         ! -- get domain decomposition
         call ESMF_GridGet(grid, staggerloc, distgrid=distgrid, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         call ESMF_DistGridGet(distgrid, deCount=deCount, dimCount=dimCount, &
            tileCount=tileCount, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out


         ungriddedCount = rank - dimCount

         if (ungriddedCount > 1) then
            call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
               msg="Only fields with one ungridded dimensions are supported", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return  ! bail out
         end if

         nullify(ungriddedLBound, ungriddedUbound)
         if (ungriddedCount > 0) then

            ! Allocate arrays for bounds and query them directly
            allocate(ungriddedLBound(ungriddedCount), ungriddedUBound(ungriddedCount), stat=localrc)
            if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            call ESMF_FieldGet(field, ungriddedLBound=ungriddedLBound, &
               ungriddedUBound=ungriddedUBound, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
         else
         end if

         allocate(minIndexPDe(dimCount, deCount), maxIndexPDe(dimCount, deCount),  &
            minIndexPTile(dimCount, tileCount), maxIndexPTile(dimCount, tileCount), &
            deToTileMap(deCount), localDeToDeMap(localDeCount), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         call ESMF_DistGridGet(distgrid, &
            minIndexPDe=minIndexPDe, maxIndexPDe=maxIndexPDe, &
            minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, &
            deToTileMap=deToTileMap, localDeToDeMap=localDeToDeMap, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         call ESMF_GridCompGetInternalState(IOComp, is, localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         do localDe = 0, localDeCount-1
            de   = localDeToDeMap(localDe+1) + 1
            tile = deToTileMap(de)

            select case (trim(action))
             case('r','read')
               call AQMIO_FieldRead(is % IO, field, &
                  minIndexPDe(:,de), maxIndexPDe(:,de), &
                  minIndexPTile(:,tile), maxIndexPTile(:,tile), &
                  ungriddedLBound=ungriddedLBound, ungriddedUBound=ungriddedUBound, &
                  variableName=variableName, timeSlice=timeSlice, localDe=localDe, &
                  rc=localrc)
               if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
             case('w','write')
               call AQMIO_FieldWrite(is % IO, field, &
                  minIndexPDe(:,de), maxIndexPDe(:,de), &
                  minIndexPTile(:,tile), maxIndexPTile(:,tile), &
                  ungriddedLBound=ungriddedLBound, ungriddedUBound=ungriddedUBound, &
                  variableName=variableName, timeSlice=timeSlice, localDe=localDe, &
                  rc=localrc)
               if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
             case default
               ! -- do nothing
            end select
         end do

         if (associated(ungriddedLBound)) then
            deallocate(ungriddedLBound, stat=localrc)
            if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
         end if

         if (associated(ungriddedUBound)) then
            deallocate(ungriddedUBound, stat=localrc)
            if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
         end if

         deallocate(minIndexPDe, maxIndexPDe, minIndexPTile, maxIndexPTile, &
            deToTileMap, localDeToDeMap, stat=localrc)
         if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

      else
         call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="I/O fields can only be defined on Grid objects.", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return  ! bail ou
      end if

   end subroutine AQMIO_FieldAccess

!------------------------------------------------------------------------------

   subroutine AQMIO_FieldRead(IO, field, &
      minIndexPDe, maxIndexPDe, minIndexPTile, maxIndexPTile, &
      ungriddedLBound, ungriddedUBound, variableName, timeSlice, localDe, rc)
      type(ioData),          intent(in)            :: IO
      type(ESMF_Field),      intent(in)            :: field
      integer, dimension(:), intent(in)            :: minIndexPDe
      integer, dimension(:), intent(in)            :: maxIndexPDe
      integer, dimension(:), intent(in)            :: minIndexPTile
      integer, dimension(:), intent(in)            :: maxIndexPTile
      integer,               intent(in),  optional :: ungriddedLBound(:)
      integer,               intent(in),  optional :: ungriddedUBound(:)
      character(len=*),      intent(in),  optional :: variableName
      integer,               intent(in),  optional :: timeSlice
      integer,               intent(in),  optional :: localDe
      integer,               intent(out), optional :: rc

      ! -- local variables
      integer :: localrc
      integer :: ilen, jlen, lbuf, lde, rank
      integer :: varId, ncStatus, ndims, xtype
      integer :: kmin, kmax, klen, uid
      integer, dimension(3) :: elb, eub
      integer,               dimension(:),     allocatable :: dimids
      integer,               dimension(:),     allocatable :: elemCount
      integer,               dimension(:),     allocatable :: elemStart
      integer(ESMF_KIND_I4), dimension(:),     allocatable :: bcstbuf_i4
      integer(ESMF_KIND_I4), dimension(:,:,:), allocatable :: buf_i4
      integer(ESMF_KIND_I4), dimension(:,:),   pointer     :: fp2d_i4 => null()
      integer(ESMF_KIND_I4), dimension(:,:,:), pointer     :: fp3d_i4 => null()
      real(ESMF_KIND_R4),    dimension(:),     allocatable :: bcstbuf_r4
      real(ESMF_KIND_R4),    dimension(:,:,:), allocatable :: buf_r4
      real(ESMF_KIND_R4),    dimension(:,:),   pointer     :: fp2d_r4 => null()
      real(ESMF_KIND_R4),    dimension(:,:,:), pointer     :: fp3d_r4 => null()
      real(ESMF_KIND_R8),    dimension(:),     allocatable :: bcstbuf_r8
      real(ESMF_KIND_R8),    dimension(:,:,:), allocatable :: buf_r8
      real(ESMF_KIND_R8),    dimension(:,:),   pointer     :: fp2d_r8 => null()
      real(ESMF_KIND_R8),    dimension(:,:,:), pointer     :: fp3d_r8 => null()
      character(len=ESMF_MAXSTR) :: fieldName, dataSetName
      type(ESMF_TypeKind_Flag) :: typekind
      type(ESMF_VM) :: vm

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

      lde = 0
      if (present(localDe)) lde = localDe

      call ESMF_FieldGet(field, name=fieldName, rank=rank, &
         typekind=typekind, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (present(variableName)) fieldName = variableName

      kmin = 1
      if (present(ungriddedLBound)) kmin = ungriddedLBound(1)

      kmax = 1
      if (present(ungriddedUBound)) kmax = ungriddedUBound(1)

      ilen = maxIndexPTile(1)-minIndexPTile(1)+1
      jlen = maxIndexPTile(2)-minIndexPTile(2)+1
      klen = kmax - kmin + 1
      lbuf = ilen * jlen * klen

      if      (typekind == ESMF_TYPEKIND_I4) then
         allocate(buf_i4(minIndexPTile(1):maxIndexPTile(1), &
            minIndexPTile(2):maxIndexPTile(2),kmin:kmax), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         buf_i4 = 0_ESMF_KIND_I4
      else if (typekind == ESMF_TYPEKIND_R4) then
         allocate(buf_r4(minIndexPTile(1):maxIndexPTile(1), &
            minIndexPTile(2):maxIndexPTile(2),kmin:kmax), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         buf_r4 = 0._ESMF_KIND_R4
      else if (typekind == ESMF_TYPEKIND_R8) then
         allocate(buf_r8(minIndexPTile(1):maxIndexPTile(1), &
            minIndexPTile(2):maxIndexPTile(2),kmin:kmax), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         buf_r8 = 0._ESMF_KIND_R8
      else
         call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="Field: "//trim(fieldName)//" - typekind not supported", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return  ! bail out
      end if

      call ESMF_GridCompGet(IO % IOLayout(lde) % taskComp, vm=vm, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (IO % IOLayout(lde) % localIOflag) then
         if (IO % IOLayout(lde) % ncid > 0) then
#if HAVE_NETCDF
            dataSetName = "NetCDF data set"

            ncStatus = nf90_inquire(IO % IOLayout(lde) % ncid, unlimitedDimId=uid)
            if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                  msg="Field "//trim(fieldName)//" not defined in "//trim(dataSetName), &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return
            end if

            ncStatus = nf90_inq_varid(IO % IOLayout(lde) % ncid, trim(fieldName), varId)
            if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                  msg="Field "//trim(fieldName)//" not defined in "//trim(dataSetName), &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return
            end if

            ncStatus = nf90_inquire_variable(IO % IOLayout(lde) % ncid, varId, &
               xtype=xtype, ndims=ndims)
            if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                  msg="Error inquiring variable "//trim(fieldName)//" in "//trim(dataSetName), &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return
            end if

            call AQMIO_VariableCheckType(fieldName, xtype, typekind, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            allocate(elemStart(ndims), elemCount(ndims), stat=localrc)
            if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
            elemStart = 1
            elemCount = 1

            if (uid == -1) then
               if (present(timeSlice)) then
                  if (timeSlice == 1) then
                     call ESMF_LogWrite("No time record found in "//trim(dataSetName) &
                        // " - proceed only for first time step", &
                        ESMF_LOGMSG_WARNING, rc=localrc)
                     if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                        line=__LINE__, &
                        file=__FILE__, &
                        rcToReturn=rc)) return  ! bail out
                  else
                     call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
                        msg="No time record found in "//dataSetName, &
                        line=__LINE__, &
                        file=__FILE__, &
                        rcToReturn=rc)
                     return  ! bail out
                  end if
               end if
            else
               allocate(dimids(ndims), stat=localrc)
               if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
               ncStatus = nf90_inquire_variable(IO % IOLayout(lde) % ncid, varId, dimIds=dimids)
               if (ncStatus /= NF90_NOERR) then
                  call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                     msg="NetCDF error", &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return
               end if
               if (dimids(ndims) == uid) then
                  if (present(timeSlice)) elemStart(ndims) = timeSlice
                  ndims = ndims - 1
               else
                  if (present(timeSlice)) then
                     call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
                        msg="No time record found for variable "//fieldName, &
                        line=__LINE__, &
                        file=__FILE__, &
                        rcToReturn=rc)
                     return  ! bail out
                  end if
               end if
               deallocate(dimids, stat=localrc)
               if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
            end if

            if (klen > 1) then
               if (rank /= ndims) localrc = ESMF_RC_ARG_INCOMP
            else
               if (rank > ndims .or. rank < ndims-1) localrc = ESMF_RC_ARG_INCOMP
            end if
            if (ESMF_LogFoundError(rcToCheck=localrc, &
               msg="Field rank incompatible with netCDF variable "//trim(fieldName), &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            elemCount(1) = ilen
            elemCount(2) = jlen
            if (ndims > 2) elemCount(3) = klen

            if      (typekind == ESMF_TYPEKIND_I4) then

               ncStatus = nf90_get_var(IO % IOLayout(lde) % ncid, varId, buf_i4, &
                  start=elemStart, count=elemCount)
               if (ncStatus /= NF90_NOERR) then
                  call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                     msg="NetCDF error", &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return
               end if

            else if (typekind == ESMF_TYPEKIND_R4) then

               ncStatus = nf90_get_var(IO % IOLayout(lde) % ncid, varId, buf_r4, &
                  start=elemStart, count=elemCount)
               if (ncStatus /= NF90_NOERR) then
                  call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                     msg="NetCDF error", &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return
               end if

            else if (typekind == ESMF_TYPEKIND_R8) then

               ncStatus = nf90_get_var(IO % IOLayout(lde) % ncid, varId, buf_r8, &
                  start=elemStart, count=elemCount)
               if (ncStatus /= NF90_NOERR) then
                  call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                     msg="NetCDF error", &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return
               end if

            end if

            deallocate(elemStart, elemCount, stat=localrc)
            if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

#else
            call ESMF_LogSetError(rcToCheck=ESMF_RC_LIB_NOT_PRESENT, &
               msg="- netCDF support is unavailable", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
#endif
         else if (IO % IOLayout(lde) % iounit > 0) then

            if      (typekind == ESMF_TYPEKIND_I4) then
               read(unit=IO % IOLayout(lde) % iounit, iostat=localrc) buf_i4
               if (localrc /= 0) then
                  call ESMF_LogSetError(ESMF_RC_FILE_READ, &
                     msg="Error reading field "//trim(fieldName), &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return  ! bail out
               end if
            else if (typekind == ESMF_TYPEKIND_R4) then
               read(unit=IO % IOLayout(lde) % iounit, iostat=localrc) buf_r4
               if (localrc /= 0) then
                  call ESMF_LogSetError(ESMF_RC_FILE_READ, &
                     msg="Error reading field "//trim(fieldName), &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return  ! bail out
               end if
            else if (typekind == ESMF_TYPEKIND_R8) then
               read(unit=IO % IOLayout(lde) % iounit, iostat=localrc) buf_r8
               if (localrc /= 0) then
                  call ESMF_LogSetError(ESMF_RC_FILE_READ, &
                     msg="Error reading field "//trim(fieldName), &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return  ! bail out
               end if
            end if
         end if
      end if

      if      (typekind == ESMF_TYPEKIND_I4) then

         allocate(bcstbuf_i4(lbuf), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         bcstbuf_i4 = reshape(buf_i4, (/lbuf/))

         call ESMF_VMBroadcast(vm, bcstbuf_i4, lbuf, 0, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         buf_i4 = reshape(bcstbuf_i4, (/ilen,jlen,klen/))
         deallocate(bcstbuf_i4, stat=localrc)
         if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         select case (rank)
          case (2)
            nullify(fp2d_i4)
            call ESMF_FieldGet(field, localDe=lde, farrayPtr=fp2d_i4, &
               exclusiveLBound=elb(1:2), exclusiveUBound=eub(1:2), rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            fp2d_i4(elb(1):eub(1),elb(2):eub(2)) = &
               buf_i4(minIndexPDe(1):maxIndexPDe(1), &
               minIndexPDe(2):maxIndexPDe(2), kmin)

          case (3)
            nullify(fp3d_i4)
            call ESMF_FieldGet(field, localDe=lde, farrayPtr=fp3d_i4, &
               exclusiveLBound=elb, exclusiveUBound=eub, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            fp3d_i4(elb(1):eub(1),elb(2):eub(2),elb(3):eub(3)) = &
               buf_i4(minIndexPDe(1):maxIndexPDe(1), &
               minIndexPDe(2):maxIndexPDe(2), &
               kmin:kmax)
         end select

         deallocate(buf_i4, bcstbuf_i4, stat=localrc)
         if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

      else if (typekind == ESMF_TYPEKIND_R4) then

         allocate(bcstbuf_r4(lbuf), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         bcstbuf_r4 = reshape(buf_r4, (/lbuf/))

         call ESMF_VMBroadcast(vm, bcstbuf_r4, lbuf, 0, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         buf_r4 = reshape(bcstbuf_r4, (/ilen,jlen,klen/))

         deallocate(bcstbuf_r4, stat=localrc)
         if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         select case (rank)
          case (2)
            call ESMF_FieldGet(field, localDe=lde, farrayPtr=fp2d_r4, &
               exclusiveLBound=elb(1:2), exclusiveUBound=eub(1:2), rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            fp2d_r4(elb(1):eub(1),elb(2):eub(2)) = &
               buf_r4(minIndexPDe(1):maxIndexPDe(1), &
               minIndexPDe(2):maxIndexPDe(2), kmin)
          case (3)
            call ESMF_FieldGet(field, localDe=lde, farrayPtr=fp3d_r4, &
               exclusiveLBound=elb, exclusiveUBound=eub, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            fp3d_r4(elb(1):eub(1),elb(2):eub(2),elb(3):eub(3)) = &
               buf_r4(minIndexPDe(1):maxIndexPDe(1), &
               minIndexPDe(2):maxIndexPDe(2), &
               kmin:kmax)
         end select

         deallocate(buf_r4, stat=localrc)
         if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

      else if (typekind == ESMF_TYPEKIND_R8) then

         allocate(bcstbuf_r8(lbuf), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         bcstbuf_r8 = reshape(buf_r8, (/lbuf/))

         call ESMF_VMBroadcast(vm, bcstbuf_r8, lbuf, 0, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         buf_r8 = reshape(bcstbuf_r8, (/ilen,jlen,klen/))

         deallocate(bcstbuf_r8, stat=localrc)
         if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         select case (rank)
          case (2)
            call ESMF_FieldGet(field, localDe=lde, farrayPtr=fp2d_r8, &
               exclusiveLBound=elb(1:2), exclusiveUBound=eub(1:2), rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            fp2d_r8(elb(1):eub(1),elb(2):eub(2)) = &
               buf_r8(minIndexPDe(1):maxIndexPDe(1), &
               minIndexPDe(2):maxIndexPDe(2), kmin)
          case (3)
            call ESMF_FieldGet(field, localDe=lde, farrayPtr=fp3d_r8, &
               exclusiveLBound=elb, exclusiveUBound=eub, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            fp3d_r8(elb(1):eub(1),elb(2):eub(2),elb(3):eub(3)) = &
               buf_r8(minIndexPDe(1):maxIndexPDe(1), &
               minIndexPDe(2):maxIndexPDe(2), &
               kmin:kmax)
         end select

         deallocate(buf_r8, stat=localrc)
         if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

      end if

   end subroutine AQMIO_FieldRead

!------------------------------------------------------------------------------

   subroutine AQMIO_DataRead(IOComp, fArray, variableName, timeSlice, localDe, rc)
      type(ESMF_GridComp),   intent(inout)         :: IOComp
      real(ESMF_KIND_R4),    pointer               :: fArray(:)
      character(len=*),      intent(in)            :: variableName
      integer,               intent(in),  optional :: timeSlice
      integer,               intent(in),  optional :: localDe
      integer,               intent(out), optional :: rc

      ! -- local variables
      integer :: localrc
      integer :: ilen, lde, rank
      integer :: varId, ncStatus, ndims, xtype
      integer :: uid
      integer :: ibuf(1)
      integer,               dimension(:),     allocatable :: dimids
      integer,               dimension(:),     allocatable :: elemCount
      integer,               dimension(:),     allocatable :: elemStart
      real(ESMF_KIND_R4),    dimension(:),     allocatable :: buf
      real(ESMF_KIND_R4),    dimension(:),     pointer     :: fp
      character(len=ESMF_MAXSTR) :: dataSetName
      type(ESMF_VM)         :: vm
      type(ioWrapper)       :: is
      type(ioData), pointer :: IO

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

      if (.not.ESMF_GridCompIsPetLocal(IOComp)) return

      call ESMF_GridCompGetInternalState(IOComp, is, localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (.not.associated(is % IO)) return
      if (.not.associated(is % IO % IOLayout)) return

      IO => is % IO

      lde = 0
      if (present(localDe)) lde = localDe

      call ESMF_GridCompGet(IO % IOLayout(lde) % taskComp, vm=vm, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (IO % IOLayout(lde) % localIOflag) then
         if (IO % IOLayout(lde) % ncid > 0) then
#if HAVE_NETCDF
            dataSetName = "NetCDF data set"

            ncStatus = nf90_inquire(IO % IOLayout(lde) % ncid, unlimitedDimId=uid)
            if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                  msg="NetCDF error", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return
            end if

            ncStatus = nf90_inq_varid(IO % IOLayout(lde) % ncid, trim(variableName), varId)
            if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                  msg="NetCDF error", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return
            end if

            ncStatus = nf90_inquire_variable(IO % IOLayout(lde) % ncid, varId, &
               xtype=xtype, ndims=ndims)
            if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                  msg="NetCDF error", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return
            end if

            call AQMIO_VariableCheckType(variableName, xtype, ESMF_TYPEKIND_R4, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            allocate(elemStart(ndims), elemCount(ndims), stat=localrc)
            if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
            elemStart = 1
            elemCount = 1

            allocate(dimids(ndims), stat=localrc)
            if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            ncStatus = nf90_inquire_variable(IO % IOLayout(lde) % ncid, varId, dimIds=dimids)
            if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                  msg="NetCDF error", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return
            end if

            if (uid == -1) then
               if (present(timeSlice)) then
                  if (timeSlice == 1) then
                     call ESMF_LogWrite("No time record found in "//trim(dataSetName) &
                        // " - proceed only for first time step", &
                        ESMF_LOGMSG_WARNING, rc=localrc)
                     if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                        line=__LINE__, &
                        file=__FILE__, &
                        rcToReturn=rc)) return  ! bail out
                  else
                     call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
                        msg="No time record found in "//dataSetName, &
                        line=__LINE__, &
                        file=__FILE__, &
                        rcToReturn=rc)
                     return  ! bail out
                  end if
               end if
            else
               if (dimids(ndims) == uid) then
                  if (present(timeSlice)) elemStart(ndims) = timeSlice
                  ndims = ndims - 1
               else
                  if (present(timeSlice)) then
                     call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
                        msg="No time record found for variable "// variableName, &
                        line=__LINE__, &
                        file=__FILE__, &
                        rcToReturn=rc)
                     return  ! bail out
                  end if
               end if
            end if

            rank = 1

            if (rank /= ndims) localrc = ESMF_RC_ARG_INCOMP
            if (ESMF_LogFoundError(rcToCheck=localrc, &
               msg="Variable rank incompatible with netCDF variable "//trim(variableName), &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            ! -- allocate array according to dimension on file
            ncStatus = nf90_inquire_dimension(IO % IOLayout(lde) % ncid, dimids(ndims), len=ibuf(1))
            if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                  msg="NetCDF error", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return
            end if
         end if
      end if

      call ESMF_VMBroadcast(vm, ibuf, 1, 0, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      ilen = ibuf(1)

      if (associated(fArray)) then
         ilen = min(ilen,size(fArray))
      else
         allocate(fp(ilen), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         fArray => fp
      end if

      fArray = 0._ESMF_KIND_R4

      if (IO % IOLayout(lde) % localIOflag) then
         if (IO % IOLayout(lde) % ncid > 0) then

            deallocate(dimids, stat=localrc)
            if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            elemCount(1) = ilen

            allocate(buf(ilen), stat=localrc)
            if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            ncStatus = nf90_get_var(IO % IOLayout(lde) % ncid, varId, fArray, &
               start=elemStart, count=elemCount)
            if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                  msg="NetCDF error", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return
            end if

            deallocate(elemStart, elemCount, stat=localrc)
            if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

         end if
      end if

      call ESMF_VMBroadcast(vm, fArray, ilen, 0, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (IO % IOLayout(lde) % localIOflag) then
         if (IO % IOLayout(lde) % ncid > 0) then
            ! -- nothing else to do
#else
            call ESMF_LogSetError(rcToCheck=ESMF_RC_LIB_NOT_PRESENT, &
               msg="- netCDF support is unavailable", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
#endif
         else if (IO % IOLayout(lde) % iounit > 0) then

            call ESMF_LogSetError(rcToCheck=ESMF_RC_LIB_NOT_PRESENT, &
               msg="- binary format is not supported", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return

         end if
      end if

   end subroutine AQMIO_DataRead

!------------------------------------------------------------------------------

   subroutine AQMIO_TimesRead(IO, variableName, timesList, localDe, rc)
      type(ioData),     intent(in)            :: IO
      character(len=*), intent(in)            :: variableName
      type(ESMF_Time),  pointer               :: timesList(:)
      integer,          intent(in),  optional :: localDe
      integer,          intent(out), optional :: rc

      ! -- local variables
      integer :: localrc
      integer :: ncStatus
      integer :: item, uid, varId, xtype, ndims, lde
      integer :: yy, mm, dd, h, m, s
      integer, dimension(:), allocatable :: dimIds, dimLen
      character(len=19), dimension(:), allocatable :: timeStrings
      type(ESMF_VM) :: vm

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

      if (associated(timesList)) then
         call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="timesList pointer must not be associated",&
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return
      end if

      lde = 0
      if (present(localDe)) lde = localDe

      if (IO % IOLayout(lde) % ncid > 0) then

#if HAVE_NETCDF
         ncStatus = nf90_inquire(IO % IOLayout(lde) % ncid, unlimitedDimId=uid)
         if (ncStatus /= NF90_NOERR) then
            call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
               msg="NetCDF error", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
         end if

         if (uid == -1) then
            call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
               msg="Time record not found", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
         end if

         ncStatus = nf90_inq_varid(IO % IOLayout(lde) % ncid, trim(variableName), varId)
         if (ncStatus /= NF90_NOERR) then
            call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
               msg="NetCDF error", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
         end if

         ncStatus = nf90_inquire_variable(IO % IOLayout(lde) % ncid, varId, &
            xtype=xtype, ndims=ndims)
         if (ncStatus /= NF90_NOERR) then
            call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
               msg="NetCDF error", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
         end if

         ! -- only time strings are supported
         if (xtype /= NF90_CHAR) then
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
               msg="Variable "//trim(variableName)//" must be string", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
         end if

         if (ndims /= 2) then
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
               msg="Variable "//trim(variableName)//" must have 2 dimensions", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
         end if

         allocate(dimIds(ndims), dimLen(ndims), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         ncStatus = nf90_inquire_variable(IO % IOLayout(lde) % ncid, varId, dimIds=dimIds)
         if (ncStatus /= NF90_NOERR) then
            call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
               msg="NetCDF error", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
         end if

         if (dimIds(ndims) /= uid) then
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
               msg="Variable "//trim(variableName)//" does not have unlimited dimensions", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
         end if

         do item = 1, ndims
            ncStatus = nf90_inquire_dimension(IO % IOLayout(lde) % ncid, dimIds(item), len=dimLen(item))
            if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                  msg="NetCDF error", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return
            end if
         end do

         if (dimLen(1) /= 19) then !YYYY-MM-DD HH:MM:SS has 19 length
            call ESMF_LogSetError(ESMF_RC_FILE_UNEXPECTED, &
               msg="String length must be 19 for variable "//trim(variableName), &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return  ! bail out
         end if

         allocate(timesList(dimLen(2)), timeStrings(dimlen(2)), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         timeStrings = ""
         ncStatus = nf90_get_var(IO % IOLayout(lde) % ncid, varId, timeStrings)
         if (ncStatus /= NF90_NOERR) then
            call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
               msg="NetCDF error", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
         end if

         do item = 1, dimlen(2)
            read(timeStrings(item), '(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)', iostat=localrc) yy, mm, dd, h, m, s
            if (localrc /= 0) then
               call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
                  msg="Unable to read timestamp", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return  ! bail out
            end if
            call ESMF_TimeSet(timesList(item), yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
         end do

         deallocate(timeStrings, stat=localrc)
         if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

#else
         call ESMF_LogSetError(rcToCheck=ESMF_RC_LIB_NOT_PRESENT, &
            msg="- AQMIO was not built with NetCDF support", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
#endif
      end if

   end subroutine AQMIO_TimesRead

!------------------------------------------------------------------------------

   subroutine AQMIO_FieldWrite(IO, field, &
      minIndexPDe, maxIndexPDe, minIndexPTile, maxIndexPTile, &
      ungriddedLBound, ungriddedUBound, variableName, timeSlice, localDe, rc)
      type(ioData),          intent(in)            :: IO
      type(ESMF_Field),      intent(in)            :: field
      integer, dimension(:), intent(in)            :: minIndexPDe
      integer, dimension(:), intent(in)            :: maxIndexPDe
      integer, dimension(:), intent(in)            :: minIndexPTile
      integer, dimension(:), intent(in)            :: maxIndexPTile
      integer,               intent(in),  optional :: ungriddedLBound(:)
      integer,               intent(in),  optional :: ungriddedUBound(:)
      character(len=*),      intent(in),  optional :: variableName
      integer,               intent(in),  optional :: timeSlice
      integer,               intent(in),  optional :: localDe
      integer,               intent(out), optional :: rc

      ! -- local variables
      integer :: localrc
      integer :: ilen, jlen, lbuf, lde
      integer :: varId, ncStatus
      integer :: ndims, rank
      integer :: kmin, kmax, klen
      integer, dimension(3) :: elb, eub
      integer, dimension(:), allocatable :: start
      integer(ESMF_KIND_I4), dimension(:),     allocatable :: recvbuf_i4
      integer(ESMF_KIND_I4), dimension(:,:,:), allocatable :: buf_i4
      integer(ESMF_KIND_I4), dimension(:,:),   pointer     :: fp2d_i4 => null()
      integer(ESMF_KIND_I4), dimension(:,:,:), pointer     :: fp3d_i4 => null()
      real(ESMF_KIND_R4),    dimension(:),     allocatable :: recvbuf_r4
      real(ESMF_KIND_R4),    dimension(:,:,:), allocatable :: buf_r4
      real(ESMF_KIND_R4),    dimension(:,:),   pointer     :: fp2d_r4 => null()
      real(ESMF_KIND_R4),    dimension(:,:,:), pointer     :: fp3d_r4 => null()
      real(ESMF_KIND_R8),    dimension(:),     allocatable :: recvbuf_r8
      real(ESMF_KIND_R8),    dimension(:,:,:), allocatable :: buf_r8
      real(ESMF_KIND_R8),    dimension(:,:),   pointer     :: fp2d_r8 => null()
      real(ESMF_KIND_R8),    dimension(:,:,:), pointer     :: fp3d_r8 => null()
      character(len=ESMF_MAXSTR) :: fieldName, dataSetName
      type(ESMF_VM) :: vm
      type(ESMF_TypeKind_Flag) :: typekind

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS


      lde = 0
      if (present(localDe)) lde = localDe

      call ESMF_FieldGet(field, name=fieldName, rank=rank, &
         typekind=typekind, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (present(variableName)) then
         fieldName = variableName
      else
      end if

      kmin = 1
      if (present(ungriddedLBound)) kmin = ungriddedLBound(1)

      kmax = 1
      if (present(ungriddedUBound)) kmax = ungriddedUBound(1)

      if (present(ungriddedLBound) .and. present(ungriddedUBound)) then
      else
      end if

      ilen = maxIndexPTile(1)-minIndexPTile(1)+1
      jlen = maxIndexPTile(2)-minIndexPTile(2)+1
      klen = kmax - kmin + 1
      lbuf = ilen * jlen * klen

      call ESMF_GridCompGet(IO % IOLayout(lde) % taskComp, vm=vm, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if      (typekind == ESMF_TYPEKIND_I4) then
         allocate(buf_i4(minIndexPTile(1):maxIndexPTile(1), &
            minIndexPTile(2):maxIndexPTile(2),kmin:kmax), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         buf_i4 = 0_ESMF_KIND_I4

         select case (rank)
          case(2)
            nullify(fp2d_i4)
            call ESMF_FieldGet(field, localDe=lde, farrayPtr=fp2d_i4, &
               exclusiveLBound=elb(1:2), exclusiveUBound=eub(1:2), rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            buf_i4(minIndexPDe(1):maxIndexPDe(1), &
               minIndexPDe(2):maxIndexPDe(2), &
               kmin) = fp2d_i4(elb(1):eub(1),elb(2):eub(2))
          case(3)
            nullify(fp3d_i4)
            call ESMF_FieldGet(field, localDe=lde, farrayPtr=fp3d_i4, &
               exclusiveLBound=elb, exclusiveUBound=eub, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            buf_i4(minIndexPDe(1):maxIndexPDe(1), &
               minIndexPDe(2):maxIndexPDe(2), &
               kmin:kmax) = fp3d_i4(elb(1):eub(1),elb(2):eub(2),elb(3):eub(3))
         end select

         allocate(recvbuf_i4(lbuf), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         call ESMF_VMReduce(vm, reshape(buf_i4, (/lbuf/)), recvbuf_i4, lbuf, &
            ESMF_REDUCE_SUM, 0, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         buf_i4 = reshape(recvbuf_i4, (/ilen, jlen, klen/))

         deallocate(recvbuf_i4, stat=localrc)
         if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

      else if (typekind == ESMF_TYPEKIND_R4) then
         allocate(buf_r4(minIndexPTile(1):maxIndexPTile(1), &
            minIndexPTile(2):maxIndexPTile(2),kmin:kmax), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         buf_r4 = 0._ESMF_KIND_R4

         select case (rank)
          case(2)
            nullify(fp2d_r4)
            call ESMF_FieldGet(field, localDe=lde, farrayPtr=fp2d_r4, &
               exclusiveLBound=elb(1:2), exclusiveUBound=eub(1:2), rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            buf_r4(minIndexPDe(1):maxIndexPDe(1), &
               minIndexPDe(2):maxIndexPDe(2), &
               kmin) = fp2d_r4(elb(1):eub(1),elb(2):eub(2))
          case(3)
            nullify(fp3d_r4)
            call ESMF_FieldGet(field, localDe=lde, farrayPtr=fp3d_r4, &
               exclusiveLBound=elb, exclusiveUBound=eub, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            buf_r4(minIndexPDe(1):maxIndexPDe(1), &
               minIndexPDe(2):maxIndexPDe(2), &
               kmin:kmax) = fp3d_r4(elb(1):eub(1),elb(2):eub(2),elb(3):eub(3))
         end select

         allocate(recvbuf_r4(lbuf), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         call ESMF_VMReduce(vm, reshape(buf_r4, (/lbuf/)), recvbuf_r4, lbuf, &
            ESMF_REDUCE_SUM, 0, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         buf_r4 = reshape(recvbuf_r4, (/ilen, jlen, klen/))

         deallocate(recvbuf_r4, stat=localrc)
         if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

      else if (typekind == ESMF_TYPEKIND_R8) then

         allocate(buf_r8(minIndexPTile(1):maxIndexPTile(1), &
            minIndexPTile(2):maxIndexPTile(2),kmin:kmax), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         buf_r8 = 0._ESMF_KIND_R8
         select case (rank)
          case(2)
            nullify(fp2d_r8)
            call ESMF_FieldGet(field, localDe=lde, farrayPtr=fp2d_r8, &
               exclusiveLBound=elb(1:2), exclusiveUBound=eub(1:2), rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            buf_r8(minIndexPDe(1):maxIndexPDe(1), &
               minIndexPDe(2):maxIndexPDe(2), &
               kmin) = fp2d_r8(elb(1):eub(1),elb(2):eub(2))
          case(3)
            nullify(fp3d_r8)
            call ESMF_FieldGet(field, localDe=lde, farrayPtr=fp3d_r8, &
               exclusiveLBound=elb, exclusiveUBound=eub, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            buf_r8(minIndexPDe(1):maxIndexPDe(1), &
               minIndexPDe(2):maxIndexPDe(2), &
               kmin:kmax) = fp3d_r8(elb(1):eub(1),elb(2):eub(2),elb(3):eub(3))
         end select

         allocate(recvbuf_r8(lbuf), stat=localrc)
         if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         call ESMF_VMReduce(vm, reshape(buf_r8, (/lbuf/)), recvbuf_r8, lbuf, &
            ESMF_REDUCE_SUM, 0, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

         buf_r8 = reshape(recvbuf_r8, (/ilen, jlen, klen/))

         deallocate(recvbuf_r8, stat=localrc)
         if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out

      else

         call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="Field: "//trim(fieldName)//" - typekind not supported", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return  ! bail out
      end if

      if (IO % IOLayout(lde) % localIOflag) then
         if (IO % IOLayout(lde) % ncid > 0) then
#if HAVE_NETCDF
            dataSetName = "NetCDF data set"

            ncStatus = nf90_inq_varid(IO % IOLayout(lde) % ncid, trim(fieldName), varId)
            if (ncStatus == NF90_ENOTVAR) then
               call AQMIO_VariableCreate(IO % IOLayout(lde), field, present(timeSlice), &
                  varId=varId, rc=localrc)
               if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
            else if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_READ, &
                  msg="Field "//trim(fieldName)//" query failed in NetCDF data set", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return  ! bail out
            else
            end if

            ncStatus = nf90_inquire_variable(IO % IOLayout(lde) % ncid, varId, ndims=ndims)
            if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_READ, &
                  msg="Error inquiring variable "//trim(fieldName)//" in NetCDF data set", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return  ! bail out
            end if

            allocate(start(ndims), stat=localrc)
            if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            start = 1
            if (present(timeSlice)) then
               start(ndims) = timeSlice
               ndims = ndims - 1
            end if

            if (ndims /= rank) then
               call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
                  msg="Variable "//trim(fieldName)//" has different rank than Field", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return  ! bail out
            end if

            if (typekind == ESMF_TYPEKIND_I4) then
               ncStatus = nf90_put_var(IO % IOLayout(lde) % ncid, varId, buf_i4, start=start)
               if (ncStatus /= NF90_NOERR) then
                  call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                     msg="NetCDF error", &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return
               end if
            else if (typekind == ESMF_TYPEKIND_R4) then
               ncStatus = nf90_put_var(IO % IOLayout(lde) % ncid, varId, buf_r4, start=start)
               if (ncStatus /= NF90_NOERR) then
                  call ESMF_LogSetError(ESMF_RC_FILE_WRITE, &
                     msg="Error writing field "//trim(fieldName)//" to NetCDF data set", &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return  ! bail out
               end if
            else if (typekind == ESMF_TYPEKIND_R8) then
               ncStatus = nf90_put_var(IO % IOLayout(lde) % ncid, varId, buf_r8, start=start)
               if (ncStatus /= NF90_NOERR) then
                  call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                     msg="NetCDF error", &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return
               end if
            end if

            deallocate(start, stat=localrc)
            if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
#else
            call ESMF_LogSetError(rcToCheck=ESMF_RC_LIB_NOT_PRESENT, &
               msg="- netCDF support is unavailable", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
#endif
         else if (IO % IOLayout(lde) % iounit > 0) then

            if      (typekind == ESMF_TYPEKIND_I4) then
               write(unit=IO % IOLayout(lde) % iounit, iostat=localrc) buf_i4
               if (localrc /= 0) then
                  call ESMF_LogSetError(ESMF_RC_FILE_READ, &
                     msg="Error writing field "//trim(fieldName), &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return  ! bail out
               end if
            else if (typekind == ESMF_TYPEKIND_R4) then
               write(unit=IO % IOLayout(lde) % iounit, iostat=localrc) buf_r4
               if (localrc /= 0) then
                  call ESMF_LogSetError(ESMF_RC_FILE_READ, &
                     msg="Error writing field "//trim(fieldName), &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return  ! bail out
               end if
            else if (typekind == ESMF_TYPEKIND_R8) then
               write(unit=IO % IOLayout(lde) % iounit, iostat=localrc) buf_r8
               if (localrc /= 0) then
                  call ESMF_LogSetError(ESMF_RC_FILE_READ, &
                     msg="Error writing field "//trim(fieldName), &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return  ! bail out
               end if
            else
               call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
                  msg="Field: "//trim(fieldName)//" - typekind not supported", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return  ! bail out

            end if
         end if
      end if

      if (typekind == ESMF_TYPEKIND_I4) then
         deallocate(buf_i4, stat=localrc)
         if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
      else if (typekind == ESMF_TYPEKIND_R4) then
         deallocate(buf_r4, stat=localrc)
         if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
      else if (typekind == ESMF_TYPEKIND_R8) then
         deallocate(buf_r8, stat=localrc)
         if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
      end if

   end subroutine AQMIO_FieldWrite

!------------------------------------------------------------------------------

   subroutine AQMIO_FileCreate(IOComp, fileName, filePath, &
      fieldList, fieldNameList, localDe, rc)
      type(ESMF_GridComp), intent(inout)         :: IOComp
      character(len=*),    intent(in)            :: fileName
      character(len=*),    intent(in),  optional :: filePath
      type(ESMF_Field),    intent(in),  optional :: fieldList(:)
      character(len=*),    intent(in),  optional :: fieldNameList(:)
      integer,             intent(in),  optional :: localDe
      integer,             intent(out), optional :: rc

      ! -- local variables
      integer :: localrc
      integer :: dimCount, item, sloc
      integer :: ncid, ncStatus, timeId, varId, xtype
      integer :: de, dimLen, tile, staggerlocCount, tileCount
      character(len=ESMF_MAXSTR) :: dimName, fieldName
      character(len=ESMF_MAXPATHLEN) :: fullName
      logical, dimension(:),   allocatable :: staggerlocList
      integer, dimension(:,:), allocatable :: dimIds
      integer, dimension(:,:), allocatable :: minIndexPTile, maxIndexPTile
      type(ioWrapper)          :: is
      type(ESMF_Grid)          :: grid
      type(ESMF_DistGrid)      :: distgrid
      type(ESMF_StaggerLoc)    :: staggerloc
      type(ESMF_TypeKind_Flag) :: typekind

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

      if (.not.ESMF_GridCompIsPetLocal(IOComp)) return

#if HAVE_NETCDF
      de = 0
      if (present(localDe)) de = localDe

      if (present(fieldList) .and. present(fieldNameList)) then
         if (size(fieldNameList) < size(fieldList)) then
            call ESMF_LogSetError(ESMF_RC_ARG_SIZE, &
               msg="size of fieldNameList must equal or larger than "// &
               "size of fieldList", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return
         end if
      end if

      call ESMF_GridCompGetInternalState(IOComp, is, localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (.not.is % IO % IOLayout(de) % localIOflag) return

      call ESMF_GridCompGet(IOComp, grid=grid, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      call ESMF_GridGet(grid, dimCount=dimCount, &
         staggerlocCount=staggerlocCount, tileCount=tileCount, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (tileCount > 1) then
         call AQMIO_FileNameGet(fullName, fileName, filePath=filePath, &
            tile=is % IO % IOLayout(de) % tile)
      else
         call AQMIO_FileNameGet(fullName, fileName, filePath=filePath)
      end if

      ! -- collect staggerloc values
      allocate(staggerlocList(0:staggerlocCount-1), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      staggerlocList = .false.
      if (present(fieldList)) then
         do item = 1, size(fieldList)
            call ESMF_FieldGet(fieldList(item), staggerloc=staggerloc, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out
            staggerlocList(staggerloc % staggerloc) = .true.
         end do
      else
         ! -- set default staggerloc as ESMF_STAGGERLOC_CENTER
         staggerlocList(ESMF_STAGGERLOC_CENTER % staggerloc) = .true.
      end if

      ncStatus = nf90_create(trim(fullName), NF90_CLOBBER, ncid)
      if (ncStatus /= NF90_NOERR) then
         call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
            msg="NetCDF error", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return
      end if

      allocate(dimIds(dimCount + 1, 0:staggerlocCount-1), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, &
         msg="Unable to allocate internal memory", &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      ! -- define dimensions
      dimIds = 0
      do sloc = 0, staggerlocCount-1

         if (staggerlocList(sloc)) then

            call ESMF_GridGet(grid, staggerloc=ESMF_StaggerLoc(sloc), &
               distgrid=distgrid, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            call ESMF_DistgridGet(distgrid, dimCount=dimCount, &
               tileCount=tileCount, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            allocate(minIndexPTile(dimCount, tileCount), &
               maxIndexPTile(dimCount, tileCount), stat=localrc)
            if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
               maxIndexPTile=maxIndexPTile, rc=localrc)
            if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

            do item = 1, dimCount
               tile = is % IO % IOLayout(de) % tile
               dimLen = maxIndexPTile(item, tile) - minIndexPTile(item, tile) + 1
               dimName = ""
               write(dimName, '("x",2i0)') sloc, item
               ncStatus = nf90_def_dim(ncid, trim(dimName), dimLen, dimIds(item,sloc))
               if (ncStatus /= NF90_NOERR) then
                  call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                     msg="NetCDF error", &
                     line=__LINE__, &
                     file=__FILE__, &
                     rcToReturn=rc)
                  return
               end if
            end do

            deallocate(minIndexPTile, maxIndexPTile, stat=localrc)
            if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)) return  ! bail out

         end if

      end do

      deallocate(staggerlocList, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      ! -- define unlimited dimension
      ncStatus = nf90_def_dim(ncid, "time", NF90_UNLIMITED, timeId)
      if (ncStatus /= NF90_NOERR) then
         call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
            msg="NetCDF error", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return
      end if

      dimIds(dimCount + 1, :) = timeId

      ! -- define Field variables
      if (present(fieldList)) then
         do item = 1, size(fieldList)

            if (present(fieldNameList)) then
               fieldName = fieldNameList(item)
            else
               call ESMF_FieldGet(fieldList(item), name=fieldName, &
                  staggerloc=staggerloc, typekind=typekind, rc=localrc)
               if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
            end if

            if      (typekind == ESMF_TYPEKIND_I4) then
               xtype = NF90_INT
            else if (typekind == ESMF_TYPEKIND_R4) then
               xtype = NF90_FLOAT
            else if (typekind == ESMF_TYPEKIND_R8) then
               xtype = NF90_DOUBLE
            else
               call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
                  msg="Field: "//trim(fieldName)//" - typekind not supported", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return  ! bail out
            end if

            ncStatus = nf90_def_var(ncid, trim(fieldName), xtype, &
               dimIds(:, staggerloc % staggerloc), varId)
            if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                  msg="NetCDF error", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return
            end if
         end do
      end if

      ncStatus = nf90_enddef(ncid)
      if (ncStatus /= NF90_NOERR) then
         call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
            msg="NetCDF error", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return
      end if

      deallocate(dimIds, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, &
         msg="Unable to deallocate internal memory for IONCCreate", &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      is % IO % IOLayout(de) % ncid = ncid
#else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_LIB_NOT_PRESENT, &
         msg="- HAVE_NETCDF not defined when lib was compiled", &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)
#endif

   end subroutine AQMIO_FileCreate

!------------------------------------------------------------------------------
!  Auxiliary methods
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  I/O component
!------------------------------------------------------------------------------

   subroutine IOCompNoOp(gcomp, importState, exportState, clock, rc)
      type(ESMF_GridComp)  :: gcomp
      type(ESMF_State)     :: importState, exportState
      type(ESMF_Clock)     :: clock
      integer, intent(out) :: rc
      rc = ESMF_SUCCESS
   end subroutine IOCompNoOp

!------------------------------------------------------------------------------

   subroutine IOCompSetServices(IOComp, rc)
      type(ESMF_GridComp)  :: IOComp
      integer, intent(out) :: rc

      rc = ESMF_SUCCESS

      call ESMF_GridCompSetEntryPoint(IOComp, ESMF_METHOD_INITIALIZE, IOCompNoOp, &
         rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) return  ! bail out

      call ESMF_GridCompSetEntryPoint(IOComp, ESMF_METHOD_RUN, IOCompNoOp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) return  ! bail out

      call ESMF_GridCompSetEntryPoint(IOComp, ESMF_METHOD_FINALIZE, IOCompNoOp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) return  ! bail out

   end subroutine IOCompSetServices

!------------------------------------------------------------------------------
!  I/O Utilities
!------------------------------------------------------------------------------

   recursive function AQMIO_StringReplaceWithInt(string, subString, intValue) &
      result (newString)
      character(len=*), intent(in) :: string
      character(len=*), intent(in) :: subString
      integer,          intent(in) :: intValue

      ! -- local variables
      integer :: idx
      character(len=ESMF_MAXPATHLEN) :: newString

      idx = index(string, subString)
      if (idx > 0) then
         write(newString, '(a,i0,a)') &
            string(1:idx-1), intValue, string(idx+len(subString):)
         newString = &
            AQMIO_StringReplaceWithInt(newString, subString, intValue)
      else
         newString = string
      end if

   end function AQMIO_StringReplaceWithInt

!------------------------------------------------------------------------------

   recursive function AQMIO_StringReplaceWithString(string, subString, &
      replaceString) result (newString)
      character(len=*), intent(in) :: string
      character(len=*), intent(in) :: subString
      character(len=*), intent(in) :: replaceString

      ! -- local variables
      integer :: idx
      character(len=ESMF_MAXPATHLEN) :: newString

      idx = index(string, subString)
      if (idx > 0) then
         newString = string(1:idx-1) // &
            replaceString // string(idx+len(subString):)
         newString = &
            AQMIO_StringReplaceWithString(newString, subString, replaceString)
      else
         newString = string
      end if

   end function AQMIO_StringReplaceWithString

!------------------------------------------------------------------------------

   subroutine AQMIO_FileNameGet(fullName, fileName, tile, filePath)
      character(len=*), intent(out)          :: fullName
      character(len=*), intent(in)           :: fileName
      integer,          intent(in), optional :: tile
      character(len=*), intent(in), optional :: filePath

      ! -- local variables
      integer :: lstr

      ! -- begin
      fullName = fileName

      if (present(filePath)) then
         if (len_trim(filePath) > 0) then
            lstr = len_trim(filePath)
            if (filePath(lstr:lstr) == "/") then
               fullName = trim(filePath) // fileName
            else
               fullName = trim(filePath) // "/" // fileName
            end if
         end if
      end if

      if (present(tile)) then
         fullName = AQMIO_StringReplaceWithInt(fullName, "<tile>", tile)
      else
         fullName = AQMIO_StringReplaceWithString(fullName, "/<tile>/", "/")
         fullName = AQMIO_StringReplaceWithString(fullName, ".<tile>.", ".")
         fullName = AQMIO_StringReplaceWithString(fullName, "<tile>", "")
      end if

   end subroutine AQMIO_FileNameGet

!------------------------------------------------------------------------------

#if HAVE_NETCDF
   subroutine AQMIO_VariableCreate(IOLayout, field, unlimited, varId, rc)
      type(AQMIOLayout), intent(in)            :: IOLayout
      type(ESMF_Field),  intent(in)            :: field
      logical,           intent(in)            :: unlimited
      integer,           intent(out), optional :: varId
      integer,           intent(out), optional :: rc

      ! -- local variables
      integer :: localrc, stat
      integer :: ncStatus
      integer :: rank, lrank
      integer :: dimCount, tileCount, tile
      integer :: item, length, dimId, lvarId, uid, ndims, xtype
      integer, dimension(:),   allocatable :: dimIds, dimLen
      integer, dimension(:),   allocatable :: ungriddedLBound, ungriddedUBound
      integer, dimension(:,:), allocatable :: minIndexPTile, maxIndexPTile
      character(len=ESMF_MAXSTR) :: fieldName
      character(len=ESMF_MAXSTR) :: dimName
      character(len=ESMF_MAXSTR) :: units
      character(len=ESMF_MAXSTR) :: description
      type(ESMF_DistGrid)      :: distgrid
      type(ESMF_Grid)          :: grid
      type(ESMF_Info)          :: info
      type(ESMF_StaggerLoc)    :: staggerloc
      type(ESMF_TypeKind_Flag) :: typekind

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS


      ncStatus = nf90_inquire(IOLayout % ncid, nDimensions=ndims, &
         unlimitedDimId=uid)
      if (ncStatus /= NF90_NOERR) then
         call ESMF_LogSetError(ESMF_RC_FILE_READ, &
            msg="Error inquiring NetCDF dataset", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return  ! bail out
      end if

      call ESMF_FieldGet(field, name=fieldName, rank=rank, grid=grid, &
         staggerloc=staggerloc, typekind=typekind, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      call ESMF_GridGet(grid, distgrid=distgrid, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      call ESMF_DistgridGet(distgrid, dimCount=dimCount, &
         tileCount=tileCount, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      allocate(minIndexPTile(dimCount, tileCount), &
         maxIndexPTile(dimCount, tileCount), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
         maxIndexPTile=maxIndexPTile, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      ncStatus = nf90_redef(IOLayout % ncid)
      if (ncStatus /= NF90_NOERR) then
         call ESMF_LogSetError(ESMF_RC_FILE_WRITE, &
            msg="Error switching to redef mode", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return  ! bail out
      end if

      lrank = rank
      if (unlimited) lrank = lrank + 1

      allocate(dimLen(lrank), dimIds(lrank), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      dimIds = -1
      dimLen = 0

      if (unlimited) then
         if (uid == -1) then
            ncStatus = nf90_def_dim(IOLayout % ncid, "Time", NF90_UNLIMITED, dimIds(lrank))
            if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_WRITE, &
                  msg="Error adding unlimited dimension", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return  ! bail out
            end if
         else
            dimIds(lrank) = uid
         end if
      else
      end if

      do item = 1, dimCount
         tile = IOLayout % tile
         dimLen(item) = maxIndexPTile(item, tile) - minIndexPTile(item, tile) + 1
      end do

      deallocate(minIndexPTile, maxIndexPTile, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (rank - dimCount > 0) then
         allocate(ungriddedLBound(rank-dimCount), ungriddedUBound(rank-dimCount), stat=stat)
         if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         call ESMF_FieldGet(field, ungriddedLBound=ungriddedLBound, &
            ungriddedUBound=ungriddedUBound, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
         dimLen(dimCount+1:rank) = ungriddedUBound - ungriddedLBound + 1
      else
      end if


      do dimId = 1, ndims
         ncStatus = nf90_inquire_dimension(IOLayout % ncid, dimId, len=length)
         if (ncStatus /= NF90_NOERR) then
            call ESMF_LogSetError(ESMF_RC_FILE_READ, &
               msg="Error inquiring existing dimension", &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return  ! bail out
         end if
         do item = 1, rank
            if (length == dimLen(item)) then
               dimIds(item) = dimId
               exit
            end if
         end do
      end do

      dimId = ndims
      do item = 1, rank
         if (dimIds(item) < 0) then
            dimid = dimid + 1
            dimName = ""
            write(dimName, '("x",i0.2)') dimid
            ncStatus = nf90_def_dim(IOLayout % ncid, trim(dimName), dimLen(item), dimIds(item))
            if (ncStatus /= NF90_NOERR) then
               call ESMF_LogSetError(ESMF_RC_FILE_WRITE, &
                  msg="Error defining dimension "//trim(dimName), &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
               return  ! bail out
            end if
         end if
      end do


      deallocate(dimLen, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if      (typekind == ESMF_TYPEKIND_I4) then
         xtype = NF90_INT
      else if (typekind == ESMF_TYPEKIND_R4) then
         xtype = NF90_FLOAT
      else if (typekind == ESMF_TYPEKIND_R8) then
         xtype = NF90_DOUBLE
      else
         call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="Field: "//trim(fieldName)//" - typekind not supported", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return  ! bail out
      end if

      ncStatus = nf90_def_var(IOLayout % ncid, trim(fieldName), xtype, dimIds, lvarId)
      if (ncStatus /= NF90_NOERR) then
         call ESMF_LogSetError(ESMF_RC_FILE_WRITE, &
            msg="Error defining NetCDF variable: "//trim(fieldName), &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return  ! bail out
      end if

      deallocate(dimIds, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      ! -- add units if available
      call ESMF_InfoGetFromHost(field, info, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      call ESMF_InfoGet(info, "units", units, default="", rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (len_trim(units) > 0) then
         ncStatus = nf90_put_att(IOLayout % ncid, lvarId, "units", trim(units))
         if (ncStatus /= NF90_NOERR) then
            call ESMF_LogSetError(ESMF_RC_FILE_WRITE, &
               msg="Error adding units to NetCDF variable: "//trim(fieldName), &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return  ! bail out
         end if
      else
      end if

      ! -- add description if available
      call ESMF_InfoGet(info, "description", description, default="", rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return  ! bail out

      if (len_trim(description) > 0) then
         ncStatus = nf90_put_att(IOLayout % ncid, lvarId, "description", trim(description))
         if (ncStatus /= NF90_NOERR) then
            call ESMF_LogSetError(ESMF_RC_FILE_WRITE, &
               msg="Error adding description to NetCDF variable: "//trim(fieldName), &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
            return  ! bail out
         end if
      else
      end if

      ncStatus = nf90_enddef(IOLayout % ncid)
      if (ncStatus /= NF90_NOERR) then
         call ESMF_LogSetError(ESMF_RC_FILE_WRITE, &
            msg="Error ending NetCDF define mode", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return  ! bail out
      end if

      if (present(varId)) varId = lvarId

   end subroutine AQMIO_VariableCreate

   subroutine AQMIO_VariableCheckType(name, xtype, typekind, rc)
      character(len=*),         intent(in)  :: name
      integer,                  intent(in)  :: xtype
      type(ESMF_TypeKind_Flag), intent(in)  :: typekind
      integer, optional,        intent(out) :: rc

      ! -- local variables
      integer          :: localrc
      logical          :: supported
      character(len=7) :: xtype_name, typekind_name

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

      ! -- identify NetCDF data type
      supported = .false.
      select case (xtype)
       case (NF90_BYTE)
         xtype_name = "byte"
       case (NF90_CHAR)
         xtype_name = "char"
       case (NF90_SHORT)
         xtype_name = "short"
         supported = .true.
       case (NF90_INT)
         xtype_name = "int"
         supported = .true.
       case (NF90_FLOAT)
         xtype_name = "float"
         supported = .true.
       case (NF90_DOUBLE)
         xtype_name = "double"
         supported = .true.
       case default
         xtype_name = "unknown"
      end select

      if (.not.supported) then
         call ESMF_LogSetError(ESMF_RC_NOT_IMPL, msg="Unsupported NetCDF data type", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return  ! bail out
      end if

      ! -- identify ESMF typekind
      supported = .true.
      if      (typekind == ESMF_TYPEKIND_I4) then
         typekind_name = "int"
      else if (typekind == ESMF_TYPEKIND_R4) then
         typekind_name = "float"
      else if (typekind == ESMF_TYPEKIND_R8) then
         typekind_name = "double"
      else
         typekind_name = "unknown"
         supported = .false.
      end if

      if (.not.supported) then
         call ESMF_LogSetError(ESMF_RC_NOT_IMPL, msg="Unsupported ESMF typekind", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
         return  ! bail out
      end if

      ! -- if not matching, issue warning
      if (xtype_name /= typekind_name) then
         call ESMF_LogWrite("Type mismatch for variable "//trim(name) &
            //" - found: "//trim(xtype_name)//", expected: "//trim(typekind_name) &
            //". Attempting automatic conversion ...", &
            logmsgFlag=ESMF_LOGMSG_WARNING, &
            line=__LINE__, &
            file=__FILE__, &
            rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return  ! bail out
      end if

   end subroutine AQMIO_VariableCheckType

   !> \brief Write 1D data directly to NetCDF file (consolidated function)
   !!
   !! This subroutine provides NetCDF I/O for writing 1D data arrays
   !! without requiring ESMF fields. Handles I4, R4, and R8 data types
   !! using case statements. Supports both creating new files and appending.
   !!
   !! \param filename NetCDF filename to write to
   !! \param data_i4 Integer data array (I4) - optional
   !! \param data_r4 Real data array (R4) - optional
   !! \param data_r8 Real data array (R8) - optional
   !! \param varname Variable name in NetCDF file
   !! \param append If true, read existing data and append new data
   !! \param del_old_file If true, delete old file before writing for the first time
   !! \param rc Return code
   !! \param current_size Current dimension size after writing - optional
   !! \param iocomp ESMF GridComp for MPI coordination - optional
   !!
   subroutine AQMIO_Write1D(filename, varname, append, del_old_file, rc, &
      data_i4, data_r4, data_r8, current_size, iocomp)
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      logical, intent(in), optional :: append
      logical, intent(in), optional :: del_old_file
      integer, intent(out), optional :: rc
      integer(ESMF_KIND_I4), intent(in), optional :: data_i4(:)
      real(ESMF_KIND_R4), intent(in), optional :: data_r4(:)
      real(ESMF_KIND_R8), intent(in), optional :: data_r8(:)
      integer, intent(out), optional :: current_size
      type(ESMF_GridComp), intent(inout), optional :: iocomp

#if HAVE_NETCDF
      ! Local variables
      integer :: localrc, ncid, varid, dimid
      logical :: file_exists, var_exists, append_mode
      character(len=ESMF_MAXPATHLEN), save :: current_filename = ""
      integer :: new_size, existing_size, total_size, data_type, netcdf_type
      type(ESMF_VM) :: vm
      integer :: localPet

      ! Data arrays for different types
      integer(ESMF_KIND_I4), allocatable :: combined_data_i4(:), existing_data_i4(:)
      real(ESMF_KIND_R4), allocatable :: combined_data_r4(:), existing_data_r4(:)
      real(ESMF_KIND_R8), allocatable :: combined_data_r8(:), existing_data_r8(:)

      if (present(rc)) rc = ESMF_SUCCESS

      ! MPI coordination: only PET 0 should perform file operations
      if (present(iocomp)) then
         call ESMF_GridCompGet(iocomp, vm=vm, rc=localrc)
         if (localrc == ESMF_SUCCESS) then
            call ESMF_VMGet(vm, localPet=localPet, rc=localrc)
            if (localrc == ESMF_SUCCESS .and. localPet /= 0) then
               ! Non-root PETs: return early (current_size will be set via broadcast)
               return
            end if
         end if
      end if

      ! Determine data type and size
      if (present(data_i4)) then
         data_type = 1
         new_size = size(data_i4)
         netcdf_type = NF90_INT
      else if (present(data_r4)) then
         data_type = 2
         new_size = size(data_r4)
         netcdf_type = NF90_REAL
      else if (present(data_r8)) then
         data_type = 3
         new_size = size(data_r8)
         netcdf_type = NF90_DOUBLE
      else
         if (present(rc)) rc = ESMF_FAILURE
         call ESMF_LogSetError(ESMF_RC_ARG_WRONG, &
            msg="No data array provided to AQMIO_Write1D", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return
      end if

      ! Set append mode (default false)
      append_mode = .false.
      if (present(append)) append_mode = append

      ! Check if file exists
      inquire(file=trim(filename), exist=file_exists)

      !delete file if we have a new file to write into for the first time
      if (present(del_old_file)) then
         if (filename /= current_filename .and. file_exists .and. del_old_file) then
            ! Use Fortran intrinsic instead of system call for safer operation
            open(unit=999, file=trim(filename), status='old', iostat=localrc)
            if (localrc == 0) then
               close(unit=999, status='delete', iostat=localrc)
               if (localrc == 0) then
                  current_filename = filename
                  file_exists = .false.  ! Update status since file was deleted
               end if
            end if
         end if
      end if

      ! open and read existing data if appending
      if (append_mode .and. file_exists) then
         ! Read existing data first based on data type
         select case (data_type)
          case (1) ! I4
            call AQMIO_Read1D(filename, varname, localrc, data_i4=existing_data_i4, iocomp=iocomp)
            if (localrc == ESMF_SUCCESS .and. allocated(existing_data_i4)) then
               existing_size = size(existing_data_i4)
               total_size = existing_size + new_size
               allocate(combined_data_i4(total_size))
               combined_data_i4(1:existing_size) = existing_data_i4(1:existing_size)
               combined_data_i4(existing_size+1:total_size) = data_i4(1:new_size)
               deallocate(existing_data_i4)
            else
               total_size = new_size
               allocate(combined_data_i4(total_size))
               combined_data_i4(1:new_size) = data_i4(1:new_size)
            end if
          case (2) ! R4
            call AQMIO_Read1D(filename, varname, localrc, data_r4=existing_data_r4, iocomp=iocomp)
            if (localrc == ESMF_SUCCESS .and. allocated(existing_data_r4)) then
               existing_size = size(existing_data_r4)
               total_size = existing_size + new_size
               allocate(combined_data_r4(total_size))
               combined_data_r4(1:existing_size) = existing_data_r4(1:existing_size)
               combined_data_r4(existing_size+1:total_size) = data_r4(1:new_size)
               deallocate(existing_data_r4)
            else
               total_size = new_size
               allocate(combined_data_r4(total_size))
               combined_data_r4(1:new_size) = data_r4(1:new_size)
            end if
          case (3) ! R8
            call AQMIO_Read1D(filename, varname, localrc, data_r8=existing_data_r8, iocomp=iocomp)
            if (localrc == ESMF_SUCCESS .and. allocated(existing_data_r8)) then
               existing_size = size(existing_data_r8)
               total_size = existing_size + new_size
               allocate(combined_data_r8(total_size))
               combined_data_r8(1:existing_size) = existing_data_r8(1:existing_size)
               combined_data_r8(existing_size+1:total_size) = data_r8(1:new_size)
               deallocate(existing_data_r8)
            else
               total_size = new_size
               allocate(combined_data_r8(total_size))
               combined_data_r8(1:new_size) = data_r8(1:new_size)
            end if
         end select
      else
         ! No append or file doesn't exist, use new data only
         total_size = new_size
         select case (data_type)
          case (1) ! I4
            allocate(combined_data_i4(total_size))
            combined_data_i4(1:new_size) = data_i4(1:new_size)
          case (2) ! R4
            allocate(combined_data_r4(total_size))
            combined_data_r4(1:new_size) = data_r4(1:new_size)
          case (3) ! R8
            allocate(combined_data_r8(total_size))
            combined_data_r8(1:new_size) = data_r8(1:new_size)
         end select
      end if

      ! Create or open file
      if (.not. file_exists) then
         localrc = nf90_create(trim(filename), NF90_CLOBBER, ncid)
      else
         localrc = nf90_open(trim(filename), NF90_WRITE, ncid)
      end if

      if (localrc /= NF90_NOERR) then
         if (present(rc)) rc = ESMF_FAILURE
         goto 999 ! Cleanup and return
      end if

      ! Check if variable exists
      localrc = nf90_inq_varid(ncid, trim(varname), varid)
      var_exists = (localrc == NF90_NOERR)

      if (.not. var_exists) then
         ! Define dimension
         if (trim(varname) == "time") then
            localrc = nf90_def_dim(ncid, "Time", NF90_UNLIMITED, dimid)
         else
            localrc = nf90_def_dim(ncid, trim(varname)//"_dim", total_size, dimid)
         end if

         if (localrc /= NF90_NOERR) then
            localrc = nf90_close(ncid)
            if (present(rc)) rc = ESMF_FAILURE
            goto 999 ! Cleanup and return
         end if

         ! Define variable
         localrc = nf90_def_var(ncid, trim(varname), netcdf_type, dimid, varid)
         if (localrc /= NF90_NOERR) then
            localrc = nf90_close(ncid)
            if (present(rc)) rc = ESMF_FAILURE
            goto 999 ! Cleanup and return
         end if

         ! Add attributes for time variable
         if (trim(varname) == "time") then
            localrc = nf90_put_att(ncid, varid, "units", "seconds since 1970-01-01 00:00:00")
            if (localrc /= NF90_NOERR) then
               call ESMF_LogWrite("Warning: Failed to add time units attribute", &
                  ESMF_LOGMSG_WARNING, rc=localrc)
            end if
            localrc = nf90_put_att(ncid, varid, "calendar", "gregorian")
            if (localrc /= NF90_NOERR) then
               call ESMF_LogWrite("Warning: Failed to add time calendar attribute", &
                  ESMF_LOGMSG_WARNING, rc=localrc)
            end if
         end if

         ! End define mode
         localrc = nf90_enddef(ncid)
         if (localrc /= NF90_NOERR) then
            localrc = nf90_close(ncid)
            if (present(rc)) rc = ESMF_FAILURE
            goto 999 ! Cleanup and return
         end if
      end if

      ! Write data based on type
      select case (data_type)
       case (1) ! I4
         localrc = nf90_put_var(ncid, varid, combined_data_i4)
       case (2) ! R4
         localrc = nf90_put_var(ncid, varid, combined_data_r4)
       case (3) ! R8
         localrc = nf90_put_var(ncid, varid, combined_data_r8)
      end select

      if (localrc /= NF90_NOERR) then
         localrc = nf90_close(ncid)
         if (present(rc)) rc = ESMF_FAILURE
         goto 999 ! Cleanup and return
      end if

      ! Set current_size output parameter if requested
      if (present(current_size)) then
         current_size = total_size
      end if

      ! Close file
      localrc = nf90_close(ncid)
      if (localrc /= NF90_NOERR) then
         if (present(rc)) rc = ESMF_FAILURE
      end if

999   continue
      ! Cleanup
      if (allocated(combined_data_i4)) deallocate(combined_data_i4)
      if (allocated(combined_data_r4)) deallocate(combined_data_r4)
      if (allocated(combined_data_r8)) deallocate(combined_data_r8)

#else
      if (present(rc)) rc = ESMF_FAILURE
      call ESMF_LogSetError(ESMF_RC_LIB_NOT_PRESENT, &
         msg="NetCDF not available", &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
#endif

   end subroutine AQMIO_Write1D

   !> \brief Read 1D data directly from NetCDF file (consolidated function)
   !!
   !! This subroutine provides NetCDF I/O for reading 1D data arrays
   !! without requiring ESMF fields. Handles I4, R4, and R8 data types
   !! using case statements. Returns allocated data array of requested type.
   !!
   !! \param filename NetCDF filename to read from
   !! \param varname Variable name in NetCDF file
   !! \param rc Return code
   !! \param data_i4 Integer data array (I4) - optional output
   !! \param data_r4 Real data array (R4) - optional output
   !! \param data_r8 Real data array (R8) - optional output
   !!
   subroutine AQMIO_Read1D(filename, varname, rc, &
      data_i4, data_r4, data_r8, iocomp)
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      integer, intent(out), optional :: rc
      integer(ESMF_KIND_I4), allocatable, intent(out), optional :: data_i4(:)
      real(ESMF_KIND_R4), allocatable, intent(out), optional :: data_r4(:)
      real(ESMF_KIND_R8), allocatable, intent(out), optional :: data_r8(:)
      type(ESMF_GridComp), intent(inout), optional :: iocomp

#if HAVE_NETCDF
      ! Local variables
      integer :: localrc, ncid, varid, var_size, data_type
      integer :: dimids(1)  ! Array to hold dimension IDs for 1D variable
      type(ESMF_VM) :: vm
      integer :: localPet

      if (present(rc)) rc = ESMF_SUCCESS

      ! MPI coordination: only PET 0 should perform file operations
      if (present(iocomp)) then
         call ESMF_GridCompGet(iocomp, vm=vm, rc=localrc)
         if (localrc == ESMF_SUCCESS) then
            call ESMF_VMGet(vm, localPet=localPet, rc=localrc)
            if (localrc == ESMF_SUCCESS .and. localPet /= 0) then
               ! Non-root PETs: return early
               return
            end if
         end if
      end if

      ! Determine requested data type
      if (present(data_i4)) then
         data_type = 1
      else if (present(data_r4)) then
         data_type = 2
      else if (present(data_r8)) then
         data_type = 3
      else
         if (present(rc)) rc = ESMF_FAILURE
         call ESMF_LogSetError(ESMF_RC_ARG_WRONG, &
            msg="No output data array provided to AQMIO_Read1D", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return
      end if

      ! Open file for reading
      localrc = nf90_open(trim(filename), NF90_NOWRITE, ncid)
      if (localrc /= NF90_NOERR) then
         if (present(rc)) rc = ESMF_FAILURE
         return
      end if

      ! Get variable ID
      localrc = nf90_inq_varid(ncid, trim(varname), varid)
      if (localrc /= NF90_NOERR) then
         localrc = nf90_close(ncid)
         if (present(rc)) rc = ESMF_FAILURE
         return
      end if

      ! Get dimension size
      localrc = nf90_inquire_variable(ncid, varid, dimids=dimids)
      if (localrc /= NF90_NOERR) then
         localrc = nf90_close(ncid)
         if (present(rc)) rc = ESMF_FAILURE
         return
      end if

      localrc = nf90_inquire_dimension(ncid, dimids(1), len=var_size)
      if (localrc /= NF90_NOERR) then
         localrc = nf90_close(ncid)
         if (present(rc)) rc = ESMF_FAILURE
         return
      end if

      ! Allocate data array based on type and read data
      select case (data_type)
       case (1) ! I4
         allocate(data_i4(var_size))
         localrc = nf90_get_var(ncid, varid, data_i4)
         if (localrc /= NF90_NOERR) then
            deallocate(data_i4)
            localrc = nf90_close(ncid)
            if (present(rc)) rc = ESMF_FAILURE
            return
         end if
       case (2) ! R4
         allocate(data_r4(var_size))
         localrc = nf90_get_var(ncid, varid, data_r4)
         if (localrc /= NF90_NOERR) then
            deallocate(data_r4)
            localrc = nf90_close(ncid)
            if (present(rc)) rc = ESMF_FAILURE
            return
         end if
       case (3) ! R8
         allocate(data_r8(var_size))
         localrc = nf90_get_var(ncid, varid, data_r8)
         if (localrc /= NF90_NOERR) then
            deallocate(data_r8)
            localrc = nf90_close(ncid)
            if (present(rc)) rc = ESMF_FAILURE
            return
         end if
      end select

      ! Close file
      localrc = nf90_close(ncid)
      if (localrc /= NF90_NOERR) then
         if (present(rc)) rc = ESMF_FAILURE
      end if

#else
      if (present(rc)) rc = ESMF_FAILURE
      call ESMF_LogSetError(ESMF_RC_LIB_NOT_PRESENT, &
         msg="NetCDF not available", &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
#endif

   end subroutine AQMIO_Read1D

#endif

!------------------------------------------------------------------------------

end module AQMIO
