!>
!! \file
!! \brief CCPr Scheme for Volcanic Emissions
!!
!!
!! Reference: Benchmarking GOCART-2G in the Goddard Earth Observing System (GEOS)
!! Allison B. Collow, Peter R. Colarco, Arlindo M. da Silva, Virginie Buchard,
!! Huisheng Bian, M Chin, Sampa Das, Ravi Govindaraju, Dongchul Kim, and Valentina Aquila,
!! Geosci. Model Development, 17, 14431468, 2024
!! https://doi.org/10.5194/gmd-17-1443-2024
!!
!! \author Lacey Holland and Wei Li
!! \date 07/2024
!!!>
module CCPr_Scheme_GOCART_SUVolcanicEmissions_Mod

   implicit none

   private

   public :: CCPr_Scheme_GOCART_SUVolcanicEmissions

contains

   !> \brief Brief description of the subroutine
   !!
   !! \param km          Number of vertical levels
   !! \param cdt         Model timestep [sec]
   !! \param VStart      Emissions Start time [sec]
   !! \param VEnd        Emissions end time [sec]
   !! \param nVolc       Number of volcanic sources
   !! \param iPoint      Grid cell index i of each volcanic source
   !! \param jPoint      Grid cell index j of each volcanic source
   !! \param hms         Current model time [sec]
   !! \param g0          Gravity [m/s^2]
   !! \param zbox        Geopotential Height difference [m] for layer
   !! \param delp        Pressure Thickness for layer [Pa]
   !! \param area        Area of grid cell [m^2]
   !! \param vSO2        Volcanic emissions  [kg S/s]
   !! \param nSO2        Index of SO2 relative to other sulfate tracers
   !! \param SO2         SO2 emissions [kg kg-1]
   !! \param SU_emis     SU emissions, kg/m2/s
   !! \param vCloud      Top elevation of emissions [m]
   !! \param vElev       Bottom elevation of emissions [m]
   !! \param vLat        Latitude specified in file [degree]
   !! \param VLon        Longitude specified in file [degree]
   !! \param RC           Success or Failure
   !!
   !!!>
   subroutine CCPr_Scheme_GOCART_SUVolcanicEmissions(km, &
      cdt, &
      VStart, &
      VEnd, &
      nVolc, &
      iPoint, &
      jPoint, &
      hms, &
      g0, &
      zbox, &
      delp, &
      area, &
      vSO2, &
      nSO2, &
      SO2, &
      vCloud, &
      vElev, &
      vLat, &
      VLon, &
      RC )

      USE GOCART2G_Process, only: SUVolcanicEmissions
      USE PrepMetVars, only:  PrepMetVarsForGOCARTSUV

      IMPLICIT NONE

      ! Arguments
      INTEGER, intent(in)                   :: km          ! number of vertical levels
      REAL, intent(in)                      :: cdt         ! model timestep [sec]
      INTEGER, intent(inout),dimension(:)   :: vStart      ! Emissions Start time [sec]
      INTEGER, intent(inout),dimension(:)   :: vEnd        ! Emissions end time [sec]
      INTEGER, intent(inout)                :: nVolc       ! number of volcanic sources
      INTEGER, intent(inout),dimension(:)   :: iPoint, jPoint ! grid cell index of each volcanic source
      !INTEGER, intent(in)                  :: YMD
      INTEGER, intent(in)                   :: hms    ! current model time [sec]
      REAL, intent(in)                      :: g0
      REAL, allocatable, DIMENSION(:) :: zbox  ! geopotential Height difference [m] for layer
      REAL, allocatable, DIMENSION(:) :: delp   ! Pressure Thickness for layer [Pa]
      REAL, intent(inout),dimension(:,:)    :: area     ! area of grid cell [m^2]
      REAL, intent(inout),dimension(:)      :: vSO2   ! volcanic emissions  [kg S/s]
      INTEGER, intent(in)                   :: nSO2     ! index of SO2 relative to other sulfate tracers
      REAL, intent(inout),dimension(:,:,:),pointer  :: SO2       ! SO2 emissions [kg kg-1]
      !REAL, intent(inout),dimension(:,:,:),pointer  :: SU_emis   ! SU emissions, kg/m2/s
      REAL, intent(inout),dimension(:)        :: vCloud    ! top elevation of emissions [m]
      REAL, intent(inout),dimension(:)        :: vElev     ! bottom elevation of emissions [m]
      REAL, intent(inout),dimension(:)        :: vLat     ! latitude specified in file [degree]
      REAL, intent(inout),dimension(:)        :: VLon     ! longitude specified in file [degree]
      INTEGER, intent(inout)                :: rc          ! error code

      !local variables
      REAL, dimension(:,:,:),pointer  :: SU_emis   ! SU emissions [kg/m2/s; not really allocated]
      REAL, DIMENSION(:,:),pointer    :: SO2EMVN   ! non-explosive volcanic emissions [kg m-2 s-1; not really allocated]
      REAL, DIMENSION(:,:),pointer    :: SO2EMVE   ! explosive volcanic emissions [kg m-2 s-1; not really allocated]
      REAL, parameter :: fMassSulfur = 32.  !  gram molecular weights of species
      REAL, parameter :: fMassSO2 = 64.     !  gram molecular weights of species
      real, pointer :: GOCART_ZBOX(:,:,:)
      real, pointer :: GOCART_DELP(:,:,:)
      character(len=256) :: errMsg
      character(len=256) :: thisLoc

      ! Initialize
      errMsg = ''
      thisLoc = ' -> at CCPr_Scheme_GOCART_SUVolcanicEmissions &
      & (in CCPr_Scheme_GOCART_SUVolcanicEmissions_mod.F90)'
      RC = 0

      ! transform data for GOCART SUVolcanicEmissions call
      call PrepMetVarsForGOCARTSUV(km,  &
         delp,            &
         zbox,           &
         GOCART_DELP,     &
         GOCART_ZBOX)

      !convert SO2 unit from Kg S to Kg SO2
      vSO2 = vSO2 * fMassSO2 / fMassSulfur

      !call gocart emission function
      if (nVolc > 0) then

         !iPoint(1) = 0
         !jPoint(1) = 0

         allocate(SU_emis(1,1,nSO2)) !TODO: nSO2 =1 for now
         allocate(SO2EMVN, SO2EMVE, mold=area)

         call SUvolcanicEmissions (nVolc, vStart, vEnd, vSO2, &
            vElev, vCloud, &
            iPoint, jPoint, &
            hms, SO2EMVN, SO2EMVE, SO2, nSO2, &
            SU_emis, km, cdt, g0, gocart_ZBOX, gocart_DELP, area, &
            vLat, vLon, rc)

      end if

      if (associated(GOCART_DELP)) nullify(GOCART_DELP)
      if (associated(GOCART_zbox)) nullify(GOCART_zbox)
      if (associated(SU_emis)) nullify(SU_emis)
      if (associated(SO2EMVN)) nullify(SO2EMVN)
      if (associated(SO2EMVE)) nullify(SO2EMVE)


   end subroutine CCPr_Scheme_GOCART_SUVolcanicEmissions

end module CCPr_Scheme_GOCART_SUVolcanicEmissions_Mod
