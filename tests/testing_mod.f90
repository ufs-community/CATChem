module testing_mod
   use precision_mod, only: fp

   implicit none

   private
   public :: assert
   public :: assert_close
   public :: load_column_data

contains

   !> Error out if `cond` is false, with optional `msg`.
   subroutine assert(cond, msg)
      logical, intent(in) :: cond  !< A conditional expression or variable
      character(len=*), intent(in), optional :: msg  !< Brief description of the assertion

      character(len=:), allocatable :: msg_

      if (.not. present(msg)) then
         msg_ = "Assertion failed"
      else
         msg_ = "Assertion '"//trim(adjustl(msg))//"' failed"
      end if

      if (.not. cond) then
         print "(a)", msg_
         stop 1
      end if
   end subroutine assert

   !> Error out if `a` and `b` are not within `tol` of each other.
   subroutine assert_close(a, b, tol, msg)
      real(fp), intent(in) :: a, b
      real(fp), intent(in), optional :: tol  !< Absolute tolerance, defaults to TINY
      character(len=*), intent(in), optional :: msg  !< Brief description of the assertion

      real(fp) :: diff, tol_
      character(len=:), allocatable :: msg_

      if (.not. present(msg)) then
         msg_ = "Closeness assertion failed"
      else
         msg_ = "Closeness assertion '"//trim(adjustl(msg))//"' failed"
      end if

      if (.not. present(tol)) then
         tol_ = tiny(1.0_fp)
      else
         tol_ = tol
      end if

      diff = abs(a - b)
      if (diff > tol_) then
         print '(a, ": ", g11.4, " != ", g11.4)', msg_, a, b
         stop 1
      end if
   end subroutine assert_close

   subroutine load_column_data(filename, MetState, rc, verbose)
      use, intrinsic :: iso_fortran_env, only: iostat_end
      use MetState_Mod, only: MetStateType

      character(len=*), intent(in) :: filename
      type(MetStateType), intent(inout) :: MetState
      integer, intent(inout) :: rc
      logical, intent(in), optional :: verbose  !< Print info (defaults to false)

      ! Local Variables
      logical :: verbose_
      integer, parameter :: unum = 99  !< Unit number for loading column data
      integer :: i
      integer :: ios  !< File I/O status
      integer :: n  !< Size of data to load
      character(len=255) :: vn  !< Variable name
      character(len=5) :: tmp

      if (.not. present(verbose)) then
         verbose_ = .false.
      else
         verbose_ = verbose
      end if

      i = 0
      vn = ""
      n = 0

      ! Load the data
      open(unit=unum, file=TRIM(filename), iostat=ios)
      if (ios /= 0) error stop
      do
         if (mod(i, 2) == 0) then
            ! Variable info
            read(unum, *, iostat=ios) vn, n
            if (ios == iostat_end) exit
            if (ios /= 0) then
               print *, "Error reading header line", i
               rc = 1
               return
            end if
            if (verbose_) then
               write(tmp, "(i0)") n
               print*, "Reading: " // trim(vn) // " | size: " // trim(tmp)
            end if
         else
            select case (vn)
               ! >>> read column data >>>
             case ("bxheight")
               read(unum, *, iostat=ios) MetState%BXHEIGHT
               if (ios /= 0) then
                  print *, "Error reading BXHEIGHT:", ios
                  rc = 1
                  return
               end if
             case ("cldf")
               read(unum, *, iostat=ios) MetState%CLDF
               if (ios /= 0) then
                  print *, "Error reading CLDF:", ios
                  rc = 1
                  return
               end if
             case ("delp")
               read(unum, *, iostat=ios) MetState%DELP
               if (ios /= 0) then
                  print *, "Error reading DELP:", ios
                  rc = 1
                  return
               end if
             case ("dluse")
               read(unum, *, iostat=ios) MetState%DLUSE
               if (ios /= 0) then
                  print *, "Error reading DLUSE:", ios
                  rc = 1
                  return
               end if
             case ("dsoiltype")
               read(unum, *, iostat=ios) MetState%DSOILTYPE
               if (ios /= 0) then
                  print *, "Error reading DSOILTYPE:", ios
                  rc = 1
                  return
               end if
             case ("eflux")
               read(unum, *, iostat=ios) MetState%EFLUX
               if (ios /= 0) then
                  print *, "Error reading EFLUX:", ios
                  rc = 1
                  return
               end if
             case ("frsno")
               read(unum, *, iostat=ios) MetState%FRSNO
               if (ios /= 0) then
                  print *, "Error reading FRSNO:", ios
                  rc = 1
                  return
               end if
             case ("frveg")
               read(unum, *, iostat=ios) MetState%FRVEG
               if (ios /= 0) then
                  print *, "Error reading FRVEG:", ios
                  rc = 1
                  return
               end if
             case ("gwetroot")
               read(unum, *, iostat=ios) MetState%GWETROOT
               if (ios /= 0) then
                  print *, "Error reading GWETROOT:", ios
                  rc = 1
                  return
               end if
             case ("gwettop")
               read(unum, *, iostat=ios) MetState%GWETTOP
               if (ios /= 0) then
                  print *, "Error reading GWETTOP:", ios
                  rc = 1
                  return
               end if
             case ("hflux")
               read(unum, *, iostat=ios) MetState%HFLUX
               if (ios /= 0) then
                  print *, "Error reading HFLUX:", ios
                  rc = 1
                  return
               end if
             case ("lwi")
               read(unum, *, iostat=ios) MetState%LWI
               if (ios /= 0) then
                  print *, "Error reading LWI:", ios
                  rc = 1
                  return
               end if
             case ("pblh")
               read(unum, *, iostat=ios) MetState%PBLH
               if (ios /= 0) then
                  print *, "Error reading PBLH:", ios
                  rc = 1
                  return
               end if
             case ("pmid")
               read(unum, *, iostat=ios) MetState%PMID
               if (ios /= 0) then
                  print *, "Error reading PMID:", ios
                  rc = 1
                  return
               end if
             case ("ps")
               read(unum, *, iostat=ios) MetState%PS
               if (ios /= 0) then
                  print *, "Error reading PS:", ios
                  rc = 1
                  return
               end if
             case ("ql")
               read(unum, *, iostat=ios) MetState%QL
               if (ios /= 0) then
                  print *, "Error reading QL:", ios
                  rc = 1
                  return
               end if
             case ("qv")
               read(unum, *, iostat=ios) MetState%QV
               if (ios /= 0) then
                  print *, "Error reading QV:", ios
                  rc = 1
                  return
               end if
             case ("snodp")
               read(unum, *, iostat=ios) MetState%SNODP
               if (ios /= 0) then
                  print *, "Error reading SNODP:", ios
                  rc = 1
                  return
               end if
             case ("snomas")
               read(unum, *, iostat=ios) MetState%SNOMAS
               if (ios /= 0) then
                  print *, "Error reading SNOMAS:", ios
                  rc = 1
                  return
               end if
             case ("soilm")
               read(unum, *, iostat=ios) MetState%SOILM
               if (ios /= 0) then
                  print *, "Error reading SOILM:", ios
                  rc = 1
                  return
               end if
             case ("sphu")
               read(unum, *, iostat=ios) MetState%SPHU
               if (ios /= 0) then
                  print *, "Error reading SPHU:", ios
                  rc = 1
                  return
               end if
             case ("swgdn")
               read(unum, *, iostat=ios) MetState%SWGDN
               if (ios /= 0) then
                  print *, "Error reading SWGDN:", ios
                  rc = 1
                  return
               end if
             case ("t")
               read(unum, *, iostat=ios) MetState%T
               if (ios /= 0) then
                  print *, "Error reading T:", ios
                  rc = 1
                  return
               end if
             case ("t2m")
               read(unum, *, iostat=ios) MetState%T2M
               if (ios /= 0) then
                  print *, "Error reading T2M:", ios
                  rc = 1
                  return
               end if
             case ("ts")
               read(unum, *, iostat=ios) MetState%TS
               if (ios /= 0) then
                  print *, "Error reading TS:", ios
                  rc = 1
                  return
               end if
             case ("u")
               read(unum, *, iostat=ios) MetState%U
               if (ios /= 0) then
                  print *, "Error reading U:", ios
                  rc = 1
                  return
               end if
             case ("u10m")
               read(unum, *, iostat=ios) MetState%U10M
               if (ios /= 0) then
                  print *, "Error reading U10M:", ios
                  rc = 1
                  return
               end if
             case ("ustar")
               read(unum, *, iostat=ios) MetState%USTAR
               if (ios /= 0) then
                  print *, "Error reading USTAR:", ios
                  rc = 1
                  return
               end if
             case ("v")
               read(unum, *, iostat=ios) MetState%V
               if (ios /= 0) then
                  print *, "Error reading V:", ios
                  rc = 1
                  return
               end if
             case ("v10m")
               read(unum, *, iostat=ios) MetState%V10M
               if (ios /= 0) then
                  print *, "Error reading V10M:", ios
                  rc = 1
                  return
               end if
             case ("wilt")
               read(unum, *, iostat=ios) MetState%WILT
               if (ios /= 0) then
                  print *, "Error reading WILT:", ios
                  rc = 1
                  return
               end if
             case ("z")
               read(unum, *, iostat=ios) MetState%Z
               if (ios /= 0) then
                  print *, "Error reading Z:", ios
                  rc = 1
                  return
               end if
             case ("z0")
               read(unum, *, iostat=ios) MetState%Z0
               if (ios /= 0) then
                  print *, "Error reading Z0:", ios
                  rc = 1
                  return
               end if
             case ("zmid")
               read(unum, *, iostat=ios) MetState%ZMID
               if (ios /= 0) then
                  print *, "Error reading ZMID:", ios
                  rc = 1
                  return
               end if
               ! <<< read column data <<<
             case default
               print *, "Variable from file not used: " // trim(vn)
               read(unum, *, iostat=ios)
               if (ios /= 0) then
                  print *, "Error advancing to next line in file:", ios
                  rc = 1
                  return
               end if
            end select
         end if
         i = i + 1
      end do
      if (verbose_) then
         ! >>> print column data >>>
         print *, "BXHEIGHT:", size(MetState%BXHEIGHT), MetState%BXHEIGHT
         print *, "CLDF:", size(MetState%CLDF), MetState%CLDF
         print *, "DELP:", size(MetState%DELP), MetState%DELP
         print *, "DLUSE:", MetState%DLUSE
         print *, "DSOILTYPE:", MetState%DSOILTYPE
         print *, "EFLUX:", MetState%EFLUX
         print *, "FRSNO:", MetState%FRSNO
         print *, "FRVEG:", MetState%FRVEG
         print *, "GWETROOT:", MetState%GWETROOT
         print *, "GWETTOP:", MetState%GWETTOP
         print *, "HFLUX:", MetState%HFLUX
         print *, "LWI:", MetState%LWI
         print *, "PBLH:", MetState%PBLH
         print *, "PMID:", size(MetState%PMID), MetState%PMID
         print *, "PS:", MetState%PS
         print *, "QL:", size(MetState%QL), MetState%QL
         print *, "QV:", size(MetState%QV), MetState%QV
         print *, "SNODP:", MetState%SNODP
         print *, "SNOMAS:", MetState%SNOMAS
         print *, "SOILM:", size(MetState%SOILM), MetState%SOILM
         print *, "SPHU:", size(MetState%SPHU), MetState%SPHU
         print *, "SWGDN:", MetState%SWGDN
         print *, "T:", size(MetState%T), MetState%T
         print *, "T2M:", MetState%T2M
         print *, "TS:", MetState%TS
         print *, "U:", size(MetState%U), MetState%U
         print *, "U10M:", MetState%U10M
         print *, "USTAR:", MetState%USTAR
         print *, "V:", size(MetState%V), MetState%V
         print *, "V10M:", MetState%V10M
         print *, "WILT:", MetState%WILT
         print *, "Z:", size(MetState%Z), MetState%Z
         print *, "Z0:", MetState%Z0
         print *, "ZMID:", size(MetState%ZMID), MetState%ZMID
         ! <<< print column data <<<
      end if

   end subroutine load_column_data

end module testing_mod
