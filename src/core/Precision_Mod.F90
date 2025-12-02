!>
!! \file precision_mod.F90
!! \brief Module PRECISION\_MOD is used to change the precision of
!!  many variables throughout catchem at compile-time.  Also contains
!!  parameters that can be used to represent missing values.
!!
!! \ingroup core_modules
!!!>
module precision_mod
   implicit none

   ! KIND parameter for 4-byte precision
   INTEGER, PARAMETER, PUBLIC :: f4 = KIND( 0.0_4 ) !< KIND parameter for 4-byte precision

   ! KIND parameter for 8-byte precision
   INTEGER, PARAMETER, PUBLIC :: f8 = KIND( 0.0_8 ) !< KIND parameter for 8-byte precision

#ifdef USE_REAL8
   ! Use 8-byte floating point precision when asked.
   INTEGER, PARAMETER, PUBLIC :: fp = f8 !< KIND parameter for 8-byte precision
#else
   ! Use 4-byte floating point by default.
   INTEGER, PARAMETER, PUBLIC :: fp = f4 !< KIND parameter for 4-byte precision
#endif

   !=========================================================================
   ! Parameters for missing values
   !=========================================================================
   LOGICAL,          PARAMETER, PUBLIC :: MISSING_BOOL = .FALSE.   !< Missing boolean value
   INTEGER,          PARAMETER, PUBLIC :: MISSING_INT  = -999      !< Missing integer value
   REAL(fp),         PARAMETER, PUBLIC :: MISSING      = -999.0_fp !< Missing real value (kind=fp)
   REAL(f4),         PARAMETER, PUBLIC :: MISSING_REAL = -999.0_f4 !< Missing real value (kind=f4)
   REAL(f8),         PARAMETER, PUBLIC :: MISSING_DBLE = -999.0_f8 !< Missing real value (kind=f8)
   CHARACTER(LEN=7), PARAMETER, PUBLIC :: MISSING_STR  = "UNKNOWN" !< Missing string

   !=========================================================================
   ! Parameters for zero
   !=========================================================================
   REAL(fp),         PARAMETER, PUBLIC :: ZERO         =  0.0_fp   !< Zero value (kind=fp)
   REAL(f4),         PARAMETER, PUBLIC :: ZERO_REAL    =  0.0_f4   !< Zero value (kind=f4)
   REAL(f8),         PARAMETER, PUBLIC :: ZERO_DBLE    =  0.0_f8   !< Zero value (kind=f8)

   !=========================================================================
   ! Parameters for very tiny numbers
   !=========================================================================
   REAL(f4),         PARAMETER, PUBLIC :: TINY_REAL    =  1.0e-16_f4 !< A small value (kind=f4)
   REAL(f8),         PARAMETER, PUBLIC :: TINY_DBLE    =  1.0e-31_f8 !< A small value (kind=f8)
#ifdef USE_REAL8
   REAL(fp),         PARAMETER, PUBLIC :: TINY_        = TINY_DBLE
#else
   REAL(fp),         PARAMETER, PUBLIC :: TINY_        = TINY_REAL
#endif

   !=========================================================================
   ! Parameters for one
   !=========================================================================
   REAL(fp),         PARAMETER, PUBLIC :: ONE          =  1.0_fp !< One value (kind=fp)
   REAL(f4),         PARAMETER, PUBLIC :: ONE_REAL     =  1.0_f4 !< One value (kind=f4)
   REAL(f8),         PARAMETER, PUBLIC :: ONE_DBLE     =  1.0_f8 !< One value (kind=f8)

   interface rae
      module procedure rae_f4, rae_f8
   end interface rae

contains

   !> Real approximately equal: `abs(a - b) < tiny(a)`
   logical function rae_f4(a, b) result(res)
      real(f4), intent(in) :: a, b
      real(f4) :: diff

      diff = abs(a - b)
      res = diff < tiny(a)
   end function rae_f4

   !> Real approximately equal: `abs(a - b) < tiny(a)`
   logical function rae_f8(a, b) result(res)
      real(f8), intent(in) :: a, b
      real(f8) :: diff

      diff = abs(a - b)
      res = diff < tiny(a)
   end function rae_f8

end module precision_mod
