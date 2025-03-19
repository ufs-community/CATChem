
module PrepMetVars_Mod

   implicit none

   private

   public :: PrepMetVarsForGOCARTSUV
   public :: INCR_REAL_RANK3, INCR_REAL_RANK2
   public :: INCR_INT_RANK3, INCR_INT_RANK2


CONTAINS


   !>
   !! \brief PrepMetVarsForGOCARTSUV - Prep the meteorological variables for GOCART SUVolcanicEmissions scheme
   !!
   !! \param [INOUT] zmid
   !! \param [INOUT] pmid
   !!
   !! \ingroup core_modules
   !!!>

   ! Need to fix below subroutine to convert one variable at a time.

   subroutine PrepMetVarsForGOCARTSUV(km,        &
      delp,            &
      zbox,           &
      GOCART_DELP,     &
      GOCART_ZBOX)


      IMPLICIT NONE

      ! INPUTS
      INTEGER, intent(in)                     :: km     ! number of vertical levels
      REAL,  intent(in), DIMENSION(:), target :: delp   ! Temperature [K]
      REAL,  intent(in), DIMENSION(:), target :: zbox  ! Height [m]

      ! INPUT/OUTPUTS
      REAL, intent(inout), pointer :: GOCART_DELP(:,:,:)   !< temperature [K]
      REAL, intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_ZBOX  !< geometric height [m]

      ! OUTPUTS - Add error handling back in late
      !INTEGER :: rc !< Return code

      ! Error handling
      !character(len=255) :: thisloc

      allocate(GOCART_DELP(1, 1, km))
      allocate(GOCART_ZBOX(1, 1, km))

      GOCART_DELP(1,1,:) = delp !  pressure  in middle of layer
      GOCART_ZBOX(1,1,:) = zbox    ! mid layer geopotential height [m]

   end subroutine PrepMetVarsForGOCARTSUV


   SUBROUTINE INCR_REAL_RANK2(ARR, RESULT)
      REAL, INTENT(IN), TARGET :: ARR
      REAL, INTENT(INOUT), POINTER :: RESULT(:,:)

      ALLOCATE(RESULT(1, 1))
      RESULT(1,1)=ARR

   END SUBROUTINE INCR_REAL_RANK2

   SUBROUTINE INCR_REAL_RANK3(ARR, RESULT)
      REAL, INTENT(IN), TARGET :: ARR(:)
      REAL, INTENT(INOUT), POINTER :: RESULT(:,:,:)

      ALLOCATE(RESULT(1, 1, SIZE(ARR, 1)))
      RESULT(1,1,:)=ARR

   END SUBROUTINE INCR_REAL_RANK3

   SUBROUTINE INCR_INT_RANK2(ARR, RESULT)
      INTEGER, INTENT(IN), TARGET :: ARR
      INTEGER, INTENT(INOUT), POINTER :: RESULT(:,:)

      ALLOCATE(RESULT(1, 1))
      RESULT(1,1)=ARR

   END SUBROUTINE INCR_INT_RANK2

   SUBROUTINE INCR_INT_RANK3(ARR, RESULT)
      INTEGER, INTENT(IN), TARGET :: ARR(:)
      INTEGER, INTENT(INOUT), POINTER :: RESULT(:,:,:)

      ALLOCATE(RESULT(1, 1, SIZE(ARR, 1)))
      RESULT(1,1,:)=ARR

   END SUBROUTINE INCR_INT_RANK3


end module PrepMetVars_Mod
