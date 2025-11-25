

# File Precision\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**Precision\_Mod.F90**](_precision___mod_8_f90.md)

[Go to the documentation of this file](_precision___mod_8_f90.md)


```Fortran

module precision_mod
   implicit none

   ! KIND parameter for 4-byte precision
   INTEGER, PARAMETER, PUBLIC :: f4 = kind( 0.0_4 ) 

   ! KIND parameter for 8-byte precision
   INTEGER, PARAMETER, PUBLIC :: f8 = kind( 0.0_8 ) 

#ifdef USE_REAL8
   ! Use 8-byte floating point precision when asked.
   INTEGER, PARAMETER, PUBLIC :: fp = f8 
#else
   ! Use 4-byte floating point by default.
   INTEGER, PARAMETER, PUBLIC :: fp = f4 
#endif

   !=========================================================================
   ! Parameters for missing values
   !=========================================================================
   LOGICAL,          PARAMETER, PUBLIC :: MISSING_BOOL = .false.   
   INTEGER,          PARAMETER, PUBLIC :: MISSING_INT  = -999
   REAL(fp),         PARAMETER, PUBLIC :: MISSING      = -999.0_fp 
   REAL(f4),         PARAMETER, PUBLIC :: MISSING_REAL = -999.0_f4 
   REAL(f8),         PARAMETER, PUBLIC :: MISSING_DBLE = -999.0_f8 
   CHARACTER(LEN=7), PARAMETER, PUBLIC :: MISSING_STR  = "UNKNOWN"

   !=========================================================================
   ! Parameters for zero
   !=========================================================================
   REAL(fp),         PARAMETER, PUBLIC :: ZERO         =  0.0_fp   
   REAL(f4),         PARAMETER, PUBLIC :: ZERO_REAL    =  0.0_f4   
   REAL(f8),         PARAMETER, PUBLIC :: ZERO_DBLE    =  0.0_f8   

   !=========================================================================
   ! Parameters for very tiny numbers
   !=========================================================================
   REAL(f4),         PARAMETER, PUBLIC :: TINY_REAL    =  1.0e-16_f4 
   REAL(f8),         PARAMETER, PUBLIC :: TINY_DBLE    =  1.0e-31_f8 
#ifdef USE_REAL8
   REAL(fp),         PARAMETER, PUBLIC :: TINY_        = tiny_dble
#else
   REAL(fp),         PARAMETER, PUBLIC :: TINY_        = tiny_real
#endif

   !=========================================================================
   ! Parameters for one
   !=========================================================================
   REAL(fp),         PARAMETER, PUBLIC :: ONE          =  1.0_fp 
   REAL(f4),         PARAMETER, PUBLIC :: ONE_REAL     =  1.0_f4 
   REAL(f8),         PARAMETER, PUBLIC :: ONE_DBLE     =  1.0_f8 

   interface rae
      module procedure rae_f4, rae_f8
   end interface rae

contains

   logical function rae_f4(a, b) result(res)
      real(f4), intent(in) :: a, b
      real(f4) :: diff

      diff = abs(a - b)
      res = diff < tiny(a)
   end function rae_f4

   logical function rae_f8(a, b) result(res)
      real(f8), intent(in) :: a, b
      real(f8) :: diff

      diff = abs(a - b)
      res = diff < tiny(a)
   end function rae_f8

end module precision_mod
```


