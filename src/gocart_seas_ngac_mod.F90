!! Revision History:
!! 03/2023, Updated it as NASA 2G GOCART, Kate.Zhang@noaa.gov
!! 06/2023, Restructure for CATChem, Jian.He@noaa.gov

module gocart_seas_ngac_mod

  use catchem_constants, only : kind_chem, g=>con_g, pi=>con_pi
  use catchem_config, only : num_emis_seas,num_chem, &
                             p_seas_1,p_seas_2,p_seas_3,p_seas_4,p_seas_5, &
                             p_eseas1,p_eseas2,p_eseas3,p_eseas4,p_eseas5
  use seas_data_mod

  implicit none

  ! -- NGAC parameters
  integer, parameter :: emission_scheme = 3    ! GEOSS 2012

  private

  public :: gocart_seas_ngac

CONTAINS

  subroutine gocart_seas_ngac(ktau,dt,u_phy,              &
          v_phy,chem_arr,dz8w,u10,         &
          v10,ustar,delp,tsk,        &
          frocean,fraci, &
          xlat,xlong,area,emis_seas,      &
          sstemisFlag,emission_scale,&
          random_factor)

     INTEGER,      INTENT(IN   ) :: ktau,sstemisFlag
     REAL(kind=kind_chem), INTENT(IN   ) :: dt,u_phy,v_phy, &
                                            dz8w,u10,v10,  &
                                            ustar,delp,           &
                                            tsk,            &
                                            frocean,fraci,        &
                                            xlat,xlong,           &
                                            area,random_factor
                                             
     REAL(kind=kind_chem), DIMENSION( num_chem ),                 &
           INTENT(INOUT ) ::                                   chem_arr
     REAL(kind=kind_chem), DIMENSION( num_emis_seas),                    &
           INTENT(OUT   ) ::                                   emis_seas
     REAL(kind=kind_chem),  DIMENSION( 5 ),                        &
            INTENT(IN   ) :: emission_scale
!
!
! local variables
!
    integer :: ipr,i,j,n,rc,ilwi
    real(kind=kind_chem) :: fsstemis, memissions, nemissions, tskin_c, ws10m
    real(kind=kind_chem) :: dummylon, fgridefficiency,deep_lakes_mask
    real(kind=kind_chem), DIMENSION (number_ss_bins) :: tc,bems
    real(kind=kind_chem) :: w10m,airmas,tskin
    real(kind=kind_chem) :: dxy

    real(kind=kind_chem) :: airmas1
    real(kind=kind_chem), dimension(number_ss_bins) :: tc1
    real(kind=kind_chem), dimension(number_ss_bins) :: bems1
    real(kind=kind_chem) :: one
!
! local parameters
!
    real(kind=kind_chem), parameter :: conver  = 1.e-9_kind_chem
    real(kind=kind_chem), parameter :: converi = 1.e+9_kind_chem
!
    one = 1.0
    emis_seas = 0.

! -- NGAC sea salt scheme
!Grid box efficiency to emission (fraction of sea water)
              
    deep_lakes_mask=1.0
    dummylon = xlong

    if( dummylon < 0.0 ) dummylon = dummylon + 360.0
    ! The Great Lakes: lon = [91W,75W], lat = [40.5N, 50N]
      if ((dummylon > 267.0) .and. &
          (dummylon < 285.0) .and. &
          (xlat >  40.5) .and. &
          (xlat <  50.0)) deep_lakes_mask = 0.0

    ! The Caspian Sea: lon = [45.0, 56], lat = 35, 48]
      if ((dummylon >  45.0) .and. &
          (dummylon <  56.0) .and. &
          (xlat >  35.0) .and. &
          (xlat <  48.0)) deep_lakes_mask = 0.0
      fgridefficiency = min(max(0.,(frocean-fraci)*deep_lakes_mask),1.)

      ! -- compute auxiliary variables
      !if (dz8w < 12.) then
      !  ws10m = sqrt(u_phy*u_phy+v_phy*v_phy)
      !else
        ws10m = sqrt(u10*u10+v10*v10)
      !end if

      ! -- compute NGAC SST correction
      if (sstemisFlag == 1) then          ! SST correction folowing Jaegle et al. 2011
        fsstemis = 0.0
        tskin_c  = tsk - 273.15
        fsstemis = (0.3 + 0.1*tskin_c - 0.0076*tskin_c**2 + 0.00021*tskin_c**3)
        fsstemis = max(fsstemis, 0.0)
      else if (sstemisFlag == 2) then
        fsstemis = 0.0
        tskin_c  = tsk - 273.15
        tskin_c  = min(max(tskin_c, -0.1), 36.0)    ! temperature range (0, 36) C

        fsstemis = (-1.107211 -0.010681*tskin_c -0.002276*tskin_c**2 &
                             + 60.288927*1.0/(40.0 - tskin_c))
        fsstemis = min(max(fsstemis, 0.0), 7.0)
      endif

      do n = 1, number_ss_bins
        memissions = 0.
        nemissions = 0.
        call SeasaltEmission( ra(n), rb(n), emission_scheme, &
                              ws10m, ustar, memissions, nemissions, rc )

        bems(n) = emission_scale(n)*fgridefficiency * fsstemis * memissions * random_factor
        tc(n) = bems(n) * dt * g / delp
      end do

      ! -- add sea salt emission increments to existing airborne concentrations
      chem_arr(p_seas_1) = chem_arr(p_seas_1) + tc(1)*converi
      chem_arr(p_seas_2) = chem_arr(p_seas_2) + tc(2)*converi
      chem_arr(p_seas_3) = chem_arr(p_seas_3) + tc(3)*converi
      chem_arr(p_seas_4) = chem_arr(p_seas_4) + tc(4)*converi
      chem_arr(p_seas_5) = chem_arr(p_seas_5) + tc(5)*converi

      ! for output diagnostics kg/m2/s
      emis_seas(p_eseas1) = bems(1)
      emis_seas(p_eseas2) = bems(2)
      emis_seas(p_eseas3) = bems(3)
      emis_seas(p_eseas4) = bems(4)
      emis_seas(p_eseas5) = bems(5)

  end subroutine gocart_seas_ngac


  subroutine SeasaltEmission ( rLow, rUp, method, w10m, ustar, &
                                memissions, nemissions, rc )

! !DESCRIPTION: Calculates the seasalt mass emission flux every timestep.
!  The particular method (algorithm) used for the calculation is based
!  on the value of "method" passed on input.  Mostly these algorithms are
!  a function of wind speed and particle size (nominally at 80% RH).
!  Routine is called once for each size bin, passing in the edge radii
!  "rLow" and "rUp" (in dry radius, units of um).  Returned in the emission
!  mass flux [kg m-2 s-1].  A sub-bin assumption is made to break (possibly)
!  large size bins into a smaller space.
!
! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   real(kind=kind_chem),    intent(in)           :: rLow, rUp   ! Dry particle bin edge radii [um]
   real(kind=kind_chem),    intent(in)           :: w10m        ! 10-m wind speed [m s-1]
   real(kind=kind_chem),    intent(in)           :: ustar       ! friction velocity [m s-1]
   integer, intent(in)           :: method      ! Algorithm to use

! !OUTPUT PARAMETERS:

   real(kind=kind_chem),    intent(inout)        :: memissions      ! Mass Emissions Flux [kg m-2 s-1]
   real(kind=kind_chem),    intent(inout)        :: nemissions      ! Number Emissions Flux [# m-2 s-1]
   integer, intent(out)          :: rc              ! Error return code:
                                                    !  0 - all is well
                                                    !  1 - 
! !Local Variables
   integer       :: ir
   real(kind=kind_chem)          :: w                               ! Intermediary wind speed [m s-1]
   real(kind=kind_chem)          :: r, dr                           ! sub-bin radius spacing (dry, um)
   real(kind=kind_chem)          :: rwet, drwet                     ! sub-bin radius spacing (rh=80%, um)
   real(kind=kind_chem)          :: aFac, bFac, scalefac, rpow, exppow, wpow

   integer, parameter :: nr = 10                    ! Number of (linear) sub-size bins

   character(len=*), parameter :: myname = 'SeasaltEmission'
   real(kind=kind_chem), parameter    :: r80fac = 1.65     ! ratio of radius(RH=0.8)/radius(RH=0.) [Gerber]
   real(kind=kind_chem), parameter    :: rhop = 2200.      ! dry seasalt density [kg m-3]

!  Define the sub-bins (still in dry radius)
   dr = (rUp - rLow)/nr
   r  = rLow + 0.5*dr

!  Loop over size bins
   nemissions = 0.
   memissions = 0.

   do ir = 1, nr

    rwet  = r80fac * r
    drwet = r80fac * dr

    select case(method)

     case(1)  ! Gong 2003
      aFac     = 4.7*(1.+30.*rwet)**(-0.017*rwet**(-1.44))
      bFac     = (0.433-log10(rwet))/0.433
      scalefac = 1.
      rpow     = 3.45
      exppow   = 1.607
      wpow     = 3.41
      w        =  w10m

     case(2)  ! Gong 1997
      aFac     = 3.
      bFac     = (0.380-log10(rwet))/0.650
      scalefac = 1.
      rpow     = 1.05
      exppow   = 1.19
      wpow     = 3.41
      w        =  w10m

     case(3)  ! GEOS5 2012
      aFac     = 4.7*(1.+30.*rwet)**(-0.017*rwet**(-1.44))
      bFac     = (0.433-log10(rwet))/0.433
      scalefac = 33.0e3
      rpow     = 3.45
      exppow   = 1.607
      wpow     = 3.41 - 1.
      w        =  ustar

     case default
!     if(chem_comm_isroot()) print *, 'SeasaltEmission missing algorithm method'
      rc = 1
      return

    end select

!   Number emissions flux (# m-2 s-1)
    nemissions = nemissions + SeasaltEmissionGong( rwet, drwet, w, scalefac, aFac, bFac, rpow, exppow, wpow )
!   Mass emissions flux (kg m-2 s-1)
    scalefac = scalefac * 4./3.*pi*rhop*r**3.*1.e-18
    memissions = memissions + SeasaltEmissionGong( rwet, drwet, w, scalefac, aFac, bFac, rpow, exppow, wpow )

    r = r + dr

   end do

   rc = 0

  end subroutine SeasaltEmission

! Function to compute sea salt emissions following the Gong style
! parameterization.  Functional form is from Gong 2003:
!  dN/dr = scalefac * 1.373 * (w^wpow) * (r^-aFac) * (1+0.057*r^rpow) *
!  10^(exppow*exp(-bFac^2))
! where r is the particle radius at 80% RH, dr is the size bin width at 80% RH,
! and w is the wind speed

  function SeasaltEmissionGong ( r, dr, w, scalefac, aFac, bFac, rpow, exppow, wpow )

   real(kind=kind_chem), intent(in)    :: r, dr     ! Wet particle radius, bin width [um]
   real(kind=kind_chem), intent(in)    :: w         ! Grid box mean wind speed [m s-1] (10-m or ustar wind)
   real(kind=kind_chem), intent(in)    :: scalefac, aFac, bFac, rpow, exppow, wpow
   real(kind=kind_chem)                :: SeasaltEmissionGong

!  Initialize
   SeasaltEmissionGong = 0.

!  Particle size distribution function
   SeasaltEmissionGong = scalefac * 1.373*r**(-aFac)*(1.+0.057*r**rpow) &
                         *10**(exppow*exp(-bFac**2.))*dr
!  Apply wind speed function
   SeasaltEmissionGong = w**wpow * SeasaltEmissionGong

  end function SeasaltEmissionGong


end module gocart_seas_ngac_mod
