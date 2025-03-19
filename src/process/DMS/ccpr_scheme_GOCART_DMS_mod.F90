!>
!! \file
!! \brief CCPr Scheme for DMS
!!
!!
!! Reference: Benchmarking GOCART-2G in the Goddard Earth Observing System (GEOS)
!! Allison B. Collow, Peter R. Colarco, Arlindo M. da Silva, Virginie Buchard,
!! Huisheng Bian, M Chin, Sampa Das, Ravi Govindaraju, Dongchul Kim, and Valentina Aquila,
!! Geosci. Model Development, 17, 14431468, 2024
!! https://doi.org/10.5194/gmd-17-1443-2024
!!
!! \author Lacey Holland
!! \date 01/2025
!!!>
module CCPr_Scheme_GOCART_DMS_Mod

   implicit none

   private

   public :: CCPr_Scheme_GOCART_DMS

contains

   !> \brief Brief description of the subroutine
   !!
   !! \param MetState     Meteorological Variables
   !! \param DiagState    Diagnostic Variables
   !! \param DMSState   DMS Variables
   !! \param RC           Success or Failure
   !!
   !! Note that other state types may be required, e.g. one specific to the process group.
   !!!>

   subroutine CCPr_Scheme_GOCART_DMS(km, cdt, g0, tmpu, u10m, v10m, lwi, delp, &
      dmso_conc, dms, SU_emis, ndms, RC)

      ! Uses
      USE GOCART2G_process, only: DMSemission
      USE PrepMetVars_Mod

      IMPLICIT NONE

      ! Arguments
      INTEGER, intent(in) :: km            ! number of vertical levels
      integer, intent(in) :: ndms      ! index of DMS relative to other sulfate tracers

      REAL, intent(in)    :: g0
      REAL, intent(in)    :: cdt               ! model timestep [sec]
      REAL, intent(in)    :: u10m                   ! 10-m u-wind component [m/sec]
      REAL, intent(in)    :: v10m                   ! 10-m v-wind component [m/sec]

      REAL, dimension(:,:),pointer  :: DMSO_CONC      ! DMS source concentration [units??]
      REAL, allocatable, DIMENSION(:) :: tmpu   ! Temperature [K]
      REAL, allocatable, DIMENSION(:) :: delp   ! Pressure Thickness for layer [Pa]

      INTEGER, intent(in)       :: lwi                   ! orography flag; Land, ocean, ice mask

      REAL, intent(inout),dimension(:,:,:),pointer  :: DMS      ! DMS [kg kg-1]
      REAL, intent(inout),dimension(:,:,:),pointer  :: SU_emis   ! SU emissions, kg/m2/s
      REAL, parameter :: fMassDMS=62.   ! g mol-1  -  should this go somewhere else in the future??

      integer, intent(out) :: RC                      ! Success or Failure

      ! Local Variables
      character(len=256) :: errMsg
      character(len=256) :: thisLoc

      real, pointer :: GOCART_TMPU(:,:,:)
      real, pointer :: GOCART_DELP(:,:,:)
      real, pointer :: GOCART_LWI(:,:)
      real, pointer :: GOCART_U10(:,:)
      real, pointer :: GOCART_V10(:,:)


      ! Initialize
      errMsg = ''
      thisLoc = ' -> at CCPr_Scheme_GOCART_DMS (in CCPr_Scheme_GOCART_DMS_mod.F90)'

      call INCR_REAL_RANK3(delp, GOCART_DELP)
      call INCR_REAL_RANK3(tmpu, GOCART_TMPU)
      call INCR_REAL_RANK2(u10m, GOCART_U10)
      call INCR_REAL_RANK2(v10m, GOCART_V10)
      call INCR_REAL_RANK2(real(LWI), GOCART_LWI)

      call DMSemission (km, cdt, g0, &
         GOCART_TMPU, GOCART_U10, &
         GOCART_V10, GOCART_LWI, &
         GOCART_DELP, fMassDMS, dmso_conc, &
         dms, SU_emis, ndms, rc)

      if (associated(GOCART_TMPU)) nullify(GOCART_TMPU)
      if (associated(GOCART_DELP)) nullify(GOCART_DELP)
      if (associated(GOCART_U10)) nullify(GOCART_U10)
      if (associated(GOCART_V10)) nullify(GOCART_V10)
      if (associated(GOCART_LWI)) nullify(GOCART_LWI)

   end subroutine CCPr_Scheme_GOCART_DMS


end module CCPr_Scheme_GOCART_DMS_Mod
