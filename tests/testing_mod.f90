module testing_mod
   use catchem_bridge_precision, only: fp

   implicit none

   private
   public :: assert
   public :: assert_close

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


end module testing_mod
