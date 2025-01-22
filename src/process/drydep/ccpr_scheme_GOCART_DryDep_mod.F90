!>
!! \file
!! \brief CCPr Scheme for dry deposition
!!
!!
!! Reference: Benchmarking GOCART-2G in the Goddard Earth Observing System (GEOS)
!! Allison B. Collow, Peter R. Colarco, Arlindo M. da Silva, Virginie Buchard,
!! Huisheng Bian, M Chin, Sampa Das, Ravi Govindaraju, Dongchul Kim, and Valentina Aquila,
!! Geosci. Model Development, 17, 14431468, 2024
!! https://doi.org/10.5194/gmd-17-1443-2024
!!
!! \author Lacey Holland
!! \date 07/2024
!!!>
module CCPr_Scheme_GOCART_DryDep_Mod

   implicit none

   private

   public :: CCPr_Scheme_GOCART_DryDep

contains

   !> \brief Brief description of the subroutine
   !!
   !! \param MetState     Meteorological Variables
   !! \param DiagState    Diagnostic Variables
   !! \param DryDepState  DryDeposition Variables
   !! \param RC           Success or Failure
   !!
   !! Note that other state types may be required, e.g. one specific to the process group.
   !!!>

   subroutine CCPr_Scheme_GOCART_DryDep(km,              &
      tmpu,            &
      rhoa,            &
      hghte,           &
      lwi,             &
      ustar,           &
      pblh,            &
      hflux,           &
      von_karman,      &
      cp,              &
      g0,              &
      z0h,             &
      drydepf,         &
      ResuspensionOpt, &
      radius,          &
      rhop,            &
      u10m,             &
      v10m,             &
      fraclake,        &
      gwettop,         &
      RC)

      ! Uses
      USE GOCART2G_Process, only: DryDeposition

      IMPLICIT NONE

      ! Arguments
      INTEGER, intent(in)                     :: km            ! number of vertical levels
      INTEGER, intent(in)             :: lwi                   ! orography flag; Land, ocean, ice mask
      REAL, allocatable, intent(in), DIMENSION(:) :: tmpu   ! Temperature [K]
      REAL, allocatable, intent(in), DIMENSION(:) :: rhoa   ! Air density [kg/m^3]
      REAL, allocatable, intent(in), DIMENSION(:) :: hghte  ! Height [m]
      REAL,  intent(in)     :: radius                                ! particle radius [m]
      REAL,  intent(in)    :: rhop                                  ! particle density [kg/m^3]
      REAL,  intent(in)               :: ustar                                 ! friction speed [m/sec]
      REAL,  intent(in)               :: pblh                                  ! PBL height [m]
      REAL,  intent(in)               :: hflux                                 ! sfc. sens. heat flux [W m-2]
      REAL,  intent(in)               :: z0h                                   ! rough height, sens. heat [m]
      REAL,  intent(in)               :: u10m                   ! 10-m u-wind component [m/sec]
      REAL,  intent(in)               :: v10m                   ! 10-m v-wind component [m/sec]
      REAL,  intent(in)               :: fraclake               ! fraction covered by water [1]
      REAL,  intent(in)               :: gwettop                ! fraction soil moisture [1]
      real, intent(in)                        :: cp
      REAL, intent(in)                        :: g0
      real, intent(in)                        :: von_karman
      logical, intent(in)                     ::  ResuspensionOpt
      ! Output
      REAL, intent(out)  :: drydepf(1,1)
      integer, intent(out) :: RC                      ! Success or Failure

      ! Local Variables
      real, pointer :: GOCART_tmpu(:,:,:)
      real, pointer :: GOCART_rhoa(:,:,:)
      real, pointer :: GOCART_HGHTE(:,:,:)
      real, pointer :: GOCART_LWI(:,:)
      real, pointer :: GOCART_USTAR(:,:)
      real, pointer :: GOCART_PBLH(:,:)
      real, pointer :: GOCART_HFLUX(:,:)
      real, pointer :: GOCART_Z0H(:,:)
      real, pointer :: GOCART_U10(:,:)
      real, pointer :: GOCART_V10(:,:)
      real, pointer :: GOCART_FRACLAKE(:,:)
      real, pointer :: GOCART_GWETTOP(:,:)

      character(len=256) :: errMsg
      character(len=256) :: thisLoc

      ! Initialize
      errMsg = ''
      thisLoc = ' -> at CCPr_Scheme_GOCART_DryDep (in CCPr_Scheme_GOCART_mod.F90)'
      RC = 0

      ! transform data for GOCART DryDeposition call

      call PrepMetVarsForGOCART(km,     &
         tmpu,            &
         rhoa,            &
         hghte,           &
         u10m,             &
         v10m,             &
         fraclake,        &
         gwettop,         &
         lwi,             &
         ustar,           &
         pblh,            &
         hflux,           &
         z0h,             &
         GOCART_tmpu,     &
         GOCART_RHOA,     &
         GOCART_HGHTE,    &
         GOCART_U10,      &
         GOCART_V10,      &
         GOCART_FRACLAKE, &
         GOCART_GWETTOP,  &
         GOCART_LWI,      &
         GOCART_USTAR,    &
         GOCART_PBLH,     &
         GOCART_HFLUX,    &
         GOCART_Z0H)

      !------------------
      ! Begin Scheme Code
      !------------------

      if (ResuspensionOpt) then
         !if (.not.present(radius)) radius=DryDepState%particleradius
         !if (.not.present(rhop)) rhop=DryDepState%particledensity
         call DryDeposition(km, GOCART_TMPU, GOCART_RHOA, GOCART_HGHTE, GOCART_LWI, GOCART_USTAR, &
            GOCART_PBLH, GOCART_HFLUX, von_karman, cp, g0, GOCART_Z0H, DRYDEPF, RC, &
            radius, rhop, GOCART_U10, GOCART_V10, GOCART_FRACLAKE, GOCART_GWETTOP)
      else
         nullify(GOCART_U10, GOCART_V10, GOCART_FRACLAKE, GOCART_GWETTOP)
         call DryDeposition(km, GOCART_TMPU, GOCART_RHOA, GOCART_HGHTE, GOCART_LWI, GOCART_USTAR, &
            GOCART_PBLH, GOCART_HFLUX, von_karman, cp, g0, GOCART_Z0H, DRYDEPF, RC)
      endif

      if (associated(GOCART_TMPU)) nullify(GOCART_TMPU)
      if (associated(GOCART_RHOA)) nullify(GOCART_RHOA)
      if (associated(GOCART_HGHTE)) nullify(GOCART_HGHTE)
      if (associated(GOCART_U10)) nullify(GOCART_U10)
      if (associated(GOCART_FRACLAKE)) nullify(GOCART_FRACLAKE)
      if (associated(GOCART_GWETTOP)) nullify(GOCART_GWETTOP)
      if (associated(GOCART_LWI)) nullify(GOCART_LWI)
      if (associated(GOCART_USTAR)) nullify(GOCART_USTAR)
      if (associated(GOCART_LWI)) nullify(GOCART_LWI)
      if (associated(GOCART_HFLUX)) nullify(GOCART_HFLUX)
      if (associated(GOCART_Z0H)) nullify(GOCART_Z0H)

! End GOCART Code


   end subroutine CCPr_Scheme_GOCART_DryDep

   !>
   !! \brief PrepMetVarsForGOCART - Prep the meteorological variables for GOCART DryDeposition scheme
   !!
   !! \param [INOUT] metstate
   !! \param [INOUT] tmpu
   !! \param [INOUT] rhoa
   !! \param [INOUT] hghte
   !! \param [INOUT] oro
   !! \param [INOUT] ustar
   !! \param [INOUT] pblh
   !! \param [INOUT] shflux
   !! \param [INOUT] z0h
   !! \param [INOUT] u10m
   !! \param [INOUT] v10m
   !! \param [INOUT] fraclake
   !! \param [INOUT] gwettop
   !! \param [OUT] rc
   !!
   !! \ingroup core_modules
   !!!>
   subroutine PrepMetVarsForGOCART(km,              &
      tmpu,            &
      rhoa,            &
      hghte,           &
      u10m,             &
      v10m,             &
      fraclake,        &
      gwettop,         &
      lwi,             &
      ustar,           &
      pblh,            &
      hflux,           &
      z0h,             &
      GOCART_tmpu,     &
      GOCART_RHOA,     &
      GOCART_HGHTE,    &
      GOCART_U10,      &
      GOCART_V10,      &
      GOCART_FRACLAKE, &
      GOCART_GWETTOP,  &
      GOCART_LWI,      &
      GOCART_USTAR,    &
      GOCART_PBLH,     &
      GOCART_HFLUX,    &
      GOCART_Z0H)



      IMPLICIT NONE

      ! INPUTS
      INTEGER, intent(in)                     :: km     ! number of vertical levels
      INTEGER,  intent(in)                    :: lwi                                    ! orography flag; Land, ocean, ice mask
      REAL,  intent(in), DIMENSION(:), target :: tmpu   ! Temperature [K]
      REAL,  intent(in), DIMENSION(:), target :: rhoa   ! Air density [kg/m^3]
      REAL,  intent(in), DIMENSION(:), target :: hghte  ! Height [m]
      REAL,  intent(in), target               :: ustar                                 ! friction speed [m/sec]
      REAL,  intent(in), target              :: pblh                                  ! PBL height [m]
      REAL,  intent(in), target              :: hflux                                 ! sfc. sens. heat flux [W m-2]
      REAL,  intent(in), target              :: z0h                                   ! rough height, sens. heat [m]
      REAL,  intent(in), target :: u10m                   ! 10-m u-wind component [m/sec]
      REAL,  intent(in), target :: v10m                   ! 10-m v-wind component [m/sec]
      REAL,  intent(in), target :: fraclake               ! fraction covered by water [1]
      REAL,  intent(in), target :: gwettop                ! fraction soil moisture [1]

      ! INPUT/OUTPUTS
      REAL, intent(inout), pointer :: GOCART_TMPU(:,:,:)   !< temperature [K]
      REAL, intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_RHOA   !< air density [kg/m^3]
      REAL, intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_HGHTE  !< geometric height [m]
      REAL, intent(inout), pointer :: GOCART_U10(:,:)                 !< 10-m u-wind component [m/sec]
      REAL, intent(inout), pointer :: GOCART_V10 (:,:)                !< 10-m v-wind component [m/sec]
      REAL, intent(inout), pointer :: GOCART_FRACLAKE(:,:)            !< fraction covered by water [1]
      REAL, intent(inout), pointer :: GOCART_GWETTOP(:,:)             !< fraction soil moisture [1]
      real, intent(inout), pointer :: GOCART_LWI(:,:)                 !< orography flag; Land, ocean, ice mask
      REAL, intent(inout), pointer :: GOCART_USTAR(:,:)               !< friction speed [m/sec]
      REAL, intent(inout), pointer :: GOCART_PBLH(:,:)                !< PBL height [m]
      REAL, intent(inout), pointer :: GOCART_HFLUX(:,:)               !< sfc. sens. heat flux [W m-2]
      REAL, intent(inout), pointer :: GOCART_Z0H(:,:)                 !< rough height, sens. heat [m]

      ! OUTPUTS - Add error handling back in late
      !INTEGER :: rc !< Return code

      ! Error handling
      !character(len=255) :: thisloc

      allocate(GOCART_TMPU(1, 1, km))
      allocate(GOCART_RHOA(1, 1, km))
      allocate(GOCART_HGHTE(1, 1, km))
      allocate(GOCART_U10(1, 1))
      allocate(GOCART_V10(1, 1))
      allocate(GOCART_FRACLAKE(1, 1))
      allocate(GOCART_GWETTOP(1, 1))
      allocate(GOCART_LWI(1, 1))
      allocate(GOCART_USTAR(1, 1))
      allocate(GOCART_PBLH(1, 1))
      allocate(GOCART_HFLUX(1, 1))
      allocate(GOCART_Z0H(1, 1))

      GOCART_TMPU(1,1,:) = tmpu ! temperature [K]
      GOCART_RHOA = reshape(rhoa, (/1, 1, km/)) ! air density [kg/m^3]
      GOCART_HGHTE = reshape(hghte, (/1, 1, km/))    ! top of layer geopotential height [m]
      GOCART_LWI = real(LWI)     ! orography flag; Land, ocean, ice mask
      GOCART_USTAR  = ustar     ! friction speed [m/sec]
      GOCART_PBLH   = pblh      ! PBL height [m]
      GOCART_HFLUX = hflux     ! sfc. sens. heat flux [W m-2]
      GOCART_Z0H    = z0h       ! rough height, sens. heat [m]
      GOCART_U10 = u10m         ! zonal wind component (E/W) [m/s]
      GOCART_V10 = v10m         ! meridional wind component (N/S) [m/s]
      GOCART_FRACLAKE = fraclake   ! unitless, lake fraction (0-1)
      GOCART_GWETTOP = gwettop     ! unitless, soil moisture fraction (0-1)


   end subroutine PrepMetVarsForGOCART

end module CCPr_Scheme_GOCART_DryDep_Mod
