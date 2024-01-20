module dust_fengsha_mod
!
!  This module developed by Barry Baker (NOAA ARL)
!  For serious questions contact barry.baker@noaa.gov
!
!  07/16/2019 - Adapted for NUOPC/GOCART, R. Montuoro
!  02/01/2020 - Adapted for FV3/CCPP, Haiqin Li
!  03/2023 - Updated by Kate.Zhang@noaa.gov to use consitent version as RRFS-SMOKE/DUST
!  06/2023, Restructure for CATChem, Jian.He@noaa.gov

  use catchem_constants, only : kind_chem, g=>con_g, pi=>con_pi
  use catchem_config, only: num_chem,num_emis_dust,&
                            p_dust_1,p_dust_2,p_dust_3,p_dust_4,p_dust_5, &
                            p_edust1,p_edust2,p_edust3,p_edust4,p_edust5
  use dust_data_mod

  implicit none

  private

  public :: gocart_dust_fengsha_driver

contains

  subroutine gocart_dust_fengsha_driver(dt,           &
          chem_arr,rho_phy,smois,delp,ssm,            &
          isltyp,vegfra,snowh,xland,flake,fsnow,area, &
          emis_dust,ust,znt,clay,sand,rdrag,uthr,     &
          num_soil_layers,random_factor)

  IMPLICIT NONE

     INTEGER,      INTENT(IN   ) :: isltyp, &
                                    num_soil_layers
     REAL(kind=kind_chem), INTENT(IN   ) :: dt,&
                                            rho_phy, &
                                            delp,&
                                            vegfra,snowh,&
                                            xland,flake,fsnow,&     
                                            area, &
                                            ust,znt, &
                                            uthr,&
                                            random_factor
     REAL(kind=kind_chem), INTENT(INOUT ) :: ssm,clay,sand,rdrag

     REAL(kind=kind_chem), DIMENSION( num_chem ),                 &
           INTENT(INOUT ) ::                                   chem_arr
     REAL(kind=kind_chem), DIMENSION(num_emis_dust),&
           OPTIONAL, INTENT(INOUT ) ::                  emis_dust
     REAL(kind=kind_chem), DIMENSION( num_soil_layers ) , &
        INTENT(INOUT) ::                               smois
     
    ! Local variables
    integer :: nmx,i,j,k,imx,jmx,lmx
    integer :: ilwi
    real(kind_chem) :: erodtot, gravsm, drylimit
    real(kind_chem), DIMENSION (5)   :: tc,bems
    real(kind_chem) :: airden,airmas,ustar,dxy
    real(kind_chem), dimension (3) :: massfrac
    real(kind_chem) :: conver,converi
    real(kind_chem) :: R
    real(kind_chem), parameter       :: ssm_thresh = 0.01    ! emit above this erodibility threshold [1]
    ! threshold values
    conver=1.e-9
    converi=1.e9

    ! Number of dust bins

    nmx=ndust

   if (isnan(ssm )) ssm=0.
   if (isnan(clay)) clay=-1.
   if (isnan(sand)) sand=-1.
   if (isnan(rdrag)) rdrag=-1.

    ! Don't do dust over water!!!

    ilwi=0
    if( xland .lt. 1.5 ) then
       ilwi=1

    ! Total concentration at lowest model level. This is still hardcoded for 5 bins.
    tc(1)=chem_arr(p_dust_1)*conver
    tc(2)=chem_arr(p_dust_2)*conver
    tc(3)=chem_arr(p_dust_3)*conver
    tc(4)=chem_arr(p_dust_4)*conver
    tc(5)=chem_arr(p_dust_5)*conver

    ! Air mass and density at lowest model level.
    airmas=area * delp / g
    airden=rho_phy
    ustar=ust
    dxy=area

    ! Mass fractions of clay, silt, and sand.
    massfrac(1)=clay
    massfrac(2)=1-(clay+sand)
    massfrac(3)=sand


             ! Total erodibility.
             
             if (isnan(ssm)) ssm=0.

             erodtot = ssm ! SUM(erod(i,j,:))
             
             ! Don't allow roughness lengths greater than 20 cm to be lofted.
             ! This kludge accounts for land use types like urban areas and
             ! forests which would otherwise show up as high dust emitters.
             ! This is a placeholder for a more widely accepted kludge
             ! factor in the literature, which reduces lofting for rough areas.
             ! Forthcoming...

             IF (znt .gt. 0.2) then
                ilwi=0
             endif

             ! limit where there is lots of vegetation
             if (vegfra .gt. .17) then
                ilwi = 0
             endif

             ! limit where there is snow on the ground
             if (snowh .gt. 0) then
                ilwi = 0
             endif

             ! Do not allow areas with bedrock, lava, or land-ice to loft

             IF (isltyp .eq. 15 .or. isltyp .eq. 16. .or. &
                  isltyp .eq. 18) then
                ilwi=0
             ENDIF
             IF (isltyp .eq. 0) then
                ilwi=0
             endif


             if ((ssm < 0.01) .or. (clay < 0.) &
              .or. (sand < 0.) .or. (rdrag < 0.)) then
                ilwi=0
             endif

             if(ilwi .ne. 0 ) then

             ! get drag partition
             ! FENGSHA uses the drag partition correction of MacKinnon et al 2004
             !     doi:10.1016/j.geomorph.2004.03.009
             if (dust_calcdrag .ne. 1) then
                call fengsha_drag(znt,R)
             else
                ! use the precalculated version derived from ASCAT; Prigent et al. (2012,2015)
                ! doi:10.1109/TGRS.2014.2338913 & doi:10.5194/amt-5-2703-2012
                ! pick only valid values
                if (rdrag > 0.) then
                  R = real(rdrag, kind=kind_chem)
                endif
             endif

             ! Call dust emission routine.
             
             call source_dust(imx,jmx, lmx, nmx, dt, tc, ustar, massfrac, & 
                  erodtot,flake,fsnow, dxy, smois(1), airden,&
                  airmas, bems, g, dust_alpha, dust_gamma, &
                  R, uthr,random_factor)

             ! convert back to concentration

             chem_arr(p_dust_1)=tc(1)*converi
             chem_arr(p_dust_2)=tc(2)*converi
             chem_arr(p_dust_3)=tc(3)*converi
             chem_arr(p_dust_4)=tc(4)*converi
             chem_arr(p_dust_5)=tc(5)*converi

             ! For output diagnostics

             emis_dust(p_edust1)=bems(1)
             emis_dust(p_edust2)=bems(2)
             emis_dust(p_edust3)=bems(3)
             emis_dust(p_edust4)=bems(4)
             emis_dust(p_edust5)=bems(5)
            endif
           endif
    !

  end subroutine gocart_dust_fengsha_driver


  subroutine source_dust(imx, jmx, lmx, nmx, dt1, tc, ustar, massfrac, &
                  erod, flake,fsnow,dxy, smois, airden, airmas, bems, g0, alpha, gamma, &
                  R, uthres,random_factor)

    ! ****************************************************************************
    ! *  Evaluate the source of each dust particles size bin by soil emission
    ! *
    ! *  Input:
    ! *         EROD      Fraction of erodible grid cell                (-)
    ! *         smois     Volumetric  soil moisture                     (m3/m3)
    ! *         ALPHA     Constant to fudge the total emission of dust  (1/m)
    ! *         GAMMA     Tuning constant for erodibility               (-)
    ! *         DXY       Surface of each grid cell                     (m2)
    ! *         AIRMAS    Mass of air for each grid box                 (kg)
    ! *         AIRDEN    Density of air for each grid box              (kg/m3)
    ! *         USTAR     Friction velocity                             (m/s)
    ! *         DT1       Time step                                     (s)
    ! *         NMX       Number of dust bins                           (-)
    ! *         IMX       Number of I points                            (-)
    ! *         JMX       Number of J points                            (-)
    ! *         LMX       Number of L points                            (-)
    ! *         R         Drag Partition                                (-)
    ! *         UTHRES    FENGSHA Dry Threshold Velocities              (m/s)
    ! *
    ! *  Data:
    ! *         MASSFRAC  Fraction of mass in each of 3 soil classes    (-) (clay silt sand) 
    ! *         DEN_DUST  Dust density                                  (kg/m3)
    ! *         DEN_SALT  Saltation particle density                    (kg/m3)
    ! *         REFF_SALT Reference saltation particle diameter         (m)
    ! *         REFF_DUST Reference dust particle diameter              (m)
    ! *         LO_DUST   Lower diameter limits for dust bins           (m)
    ! *         UP_DUST   Upper diameter limits for dust bins           (m)
    ! *         FRAC_SALT Soil class mass fraction for saltation bins   (-)
    ! *
    ! *  Parameters:
    ! *         CMB       Constant of proportionality                   (-)
    ! *         MMD_DUST  Mass median diameter of dust                  (m)
    ! *         GSD_DUST  Geometric standard deviation of dust          (-)
    ! *         LAMBDA    Side crack propagation length                 (m)
    ! *         CV        Normalization constant                        (-)
    ! *         G0        Gravitational acceleration                    (m/s2)
    ! *
    ! *  Working:
    ! *         RHOA      Density of air in cgs                         (g/cm3)
    ! *         DS_REL    Saltation surface area distribution           (-)
    ! *         DLNDP     Dust bin width                                (-)
    ! *         EMIT      Total vertical mass flux                      (kg/m2/s)
    ! *         EMIT_VOL  Total vertical volume flux                    (m/s)
    ! *         DSRC      Mass of emitted dust               (kg/timestep/cell)
    ! *
    ! *  Output:
    ! *         TC        Total concentration of dust        (kg/kg/timestep/cell)
    ! *         BEMS      Source of each dust type           (kg/timestep/cell)
    ! *
    ! ****************************************************************************
    implicit none

    ! Input
    INTEGER,            INTENT(IN)    :: imx,jmx,lmx,nmx
    REAL(kind_chem), INTENT(IN)    :: dt1
    REAL(kind_chem), INTENT(IN)    :: ustar
    REAL(kind_chem), INTENT(IN)    :: massfrac(3)
    REAL(kind_chem), INTENT(IN)    :: erod
    REAL(kind_chem), INTENT(IN)    :: flake
    REAL(kind_chem), INTENT(IN)    :: fsnow
    REAL(kind_chem), INTENT(IN)    :: dxy
    REAL(kind_chem), INTENT(IN)    :: smois
    REAL(kind_chem), INTENT(IN)    :: airden
    REAL(kind_chem), INTENT(IN)    :: airmas
    REAL(kind_chem), INTENT(IN)    :: g0
    REAL(kind_chem), INTENT(IN)    :: alpha
    REAL(kind_chem), INTENT(IN)    :: gamma
    REAL(kind_chem), INTENT(IN)    :: R
    REAL(kind_chem), INTENT(IN)    :: uthres
    REAL(kind_chem), INTENT(IN)    :: random_factor

    ! Output
    REAL(kind_chem), INTENT(INOUT) :: tc(nmx)

    ! Local Variables
    REAL(kind_chem), INTENT(OUT)   :: bems(nmx)
    
    REAL(kind_chem) :: dvol(nmx)
    REAL(kind_chem) :: distr_dust(nmx)
    REAL(kind_chem) :: dlndp(nmx)
    REAL(kind_chem) :: dsrc
    REAL(kind_chem) :: dvol_tot
    REAL(kind_chem) :: emit
    REAL(kind_chem) :: emit_vol
    REAL(kind_chem) :: rhoa
    INTEGER   :: i, j, n

    ! Constant of proportionality from Marticorena et al, 1997 (unitless)
    ! Arguably more ~consistent~ fudge than alpha, which has many walnuts
    ! sprinkled throughout the literature. - GC

    REAL(kind_chem), PARAMETER :: cmb=1.0
    REAL(kind_chem), PARAMETER :: kvhmax=2.0e-4

    ! Parameters used in Kok distribution function. Advise not to play with
    ! these without the expressed written consent of someone who knows what
    ! they're doing. - GC

    REAL(kind_chem), PARAMETER :: mmd_dust=3.4D-6  ! median mass diameter (m)
    REAL(kind_chem), PARAMETER :: gsd_dust=3.0     ! geom. std deviation
    REAL(kind_chem), PARAMETER :: lambda=12.0D-6   ! crack propagation length (m)
    REAL(kind_chem), PARAMETER :: cv=12.62D-6      ! normalization constant
!    REAL(kind_chem), PARAMETER :: lambda=8.0D-6   ! crack propogation length(m)
!    REAL(kind_chem), PARAMETER :: cv=6.83D-6      ! normalization constant
    REAL(kind_chem), PARAMETER :: RHOSOIL=2650.


    ! calculate the total vertical dust flux 

    emit = 0.0

    call DustEmissionFENGSHA(flake,fsnow,smois,massfrac(1),massfrac(3), massfrac(2), &
                                erod, R, airden, ustar, uthres, alpha, gamma, kvhmax, &
                                g0, RHOSOIL, emit)

    ! Now that we have the total dust emission, distribute into dust bins using
    ! lognormal distribution (Dr. Jasper Kok, in press), and
    ! calculate total mass emitted over the grid box over the timestep.
    !
    ! In calculating the Kok distribution, we assume upper and lower limits to each bin.
    ! For reff_dust=(/0.73D-6,1.4D-6,2.4D-6,4.5D-6,8.0D-6/) (default),
    ! lower limits were ASSUMED at lo_dust=(/0.1D-6,1.0D-6,1.8D-6,3.0D-6,6.0D-6/)
    ! upper limits were ASSUMED at up_dust=(/1.0D-6,1.8D-6,3.0D-6,6.0D-6,10.0D-6/)
    ! These may be changed within module_data_gocart_dust.F, but make sure it is
    ! consistent with reff_dust values.  These values were taken from the original
    ! GOCART bin configuration. We use them here to calculate dust bin width, dlndp.
    ! dVol is the volume distribution. You know...if you were wondering. GC

    dvol_tot=0.
    DO n=1,nmx
       dlndp(n)=LOG(up_dust(n)/lo_dust(n))
       dvol(n)=(2.0*reff_dust(n)/cv)*(1.+ERF(LOG(2.0*reff_dust(n)/mmd_dust)/(SQRT(2.)*LOG(gsd_dust))))*&
            EXP(-(2.0*reff_dust(n)/lambda)**3.0)*dlndp(n)
       dvol_tot=dvol_tot+dvol(n)
       ! Convert mass flux to volume flux
       !emit_vol=emit/den_dust(n) ! (m s^-1)
    END DO
    DO n=1,nmx
       distr_dust(n)=dvol(n)/dvol_tot
       !print *,"distr_dust(",n,")=",distr_dust(n)
    END DO

    ! Now distribute total vertical emission into dust bins and update concentration.

    DO n=1,nmx
       ! Calculate total mass emitted
       dsrc = emit*distr_dust(n)*dxy*dt1*random_factor  ! (kg)
       IF (dsrc < 0.0) dsrc = 0.0
       
       ! Update dust mixing ratio at first model level.
       tc(n) = tc(n) + dsrc / airmas ! (kg/kg)
       !   bems(i,j,n) = dsrc  ! diagnostic
       !bems(i,j,n) = 1000.*dsrc/(dxy(j)*dt1) ! diagnostic (g/m2/s)
       bems(n) = 1.e+9*dsrc/(dxy*dt1) ! diagnostic (ug/m2/s) !lzhang
       
    END DO

  END SUBROUTINE source_dust


  subroutine fengsha_drag(z0,R)
    implicit none

    real(kind_chem), intent(in) :: z0
    real(kind_chem), intent(out) :: R
    real(kind_chem), parameter :: z0s = 1.0e-04 !Surface roughness for ideal bare surface [m]
    ! ------------------------------------------------------------------------
    ! Function: Calculates the MacKinnon et al. 2004 Drag Partition Correction
    !
    !   R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)
    !
    !--------------------------------------------------------------------------
    ! Drag partition correction. See MacKinnon et al. (2004),
    !     doi:10.1016/j.geomorph.2004.03.009
    R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)

    ! Drag partition correction. See Marticorena et al. (1997),
    !     doi:10.1029/96JD02964
    !R = 1.0 - log(z0 / z0s) / log( 0.7 * (10./z0s) ** 0.8)

    return
  end subroutine fengsha_drag

  subroutine DustEmissionFENGSHA(flake,fsnow,slc, clay, sand, silt,  &
                                  ssm, rdrag, airdens, ustar, uthrs, alpha, gamma, &
                                  kvhmax, grav, rhop, emissions)
    
    ! !USES:
    implicit NONE
    
! !INPUT PARAMETERS:
    REAL(kind_chem), intent(in) :: flake      ! liquid water content of soil layer, volumetric fraction [1]
    REAL(kind_chem), intent(in) :: fsnow      ! liquid water content of soil layer, volumetric fraction [1]
    REAL(kind_chem), intent(in) :: slc      ! liquid water content of soil layer, volumetric fraction [1]
    REAL(kind_chem), intent(in) :: clay     ! fractional clay content [1]
    REAL(kind_chem), intent(in) :: sand     ! fractional sand content [1]
    REAL(kind_chem), intent(in) :: silt     ! fractional silt content [1]
    REAL(kind_chem), intent(in) :: ssm      ! erosion map [1]
    REAL(kind_chem), intent(in) :: rdrag    ! drag partition [1/m]
    REAL(kind_chem), intent(in) :: airdens  ! air density at lowest level [kg/m^3]
    REAL(kind_chem), intent(in) :: ustar    ! friction velocity [m/sec]
    REAL(kind_chem), intent(in) :: uthrs    ! threshold velocity [m/2]
    REAL(kind_chem), intent(in) :: alpha    ! scaling factor [1]
    REAL(kind_chem), intent(in) :: gamma    ! scaling factor [1]
    REAL(kind_chem), intent(in) :: kvhmax   ! max. vertical to horizontal mass flux ratio [1]
    REAL(kind_chem), intent(in) :: grav     ! gravity [m/sec^2]
    REAL(kind_chem), intent(in) :: rhop     ! soil class density [kg/m^3]
    
    ! !OUTPUT PARAMETERS:
    REAL(kind_chem), intent(inout) :: emissions ! binned surface emissions [kg/(m^2 sec)]
    
    ! !DESCRIPTION: Compute dust emissions using NOAA/ARL FENGSHA model
    !
    ! !REVISION HISTORY:
    !
    ! 22Feb2020 B.Baker/NOAA    - Original implementation
    ! 29Mar2021 R.Montuoro/NOAA - Refactored for process library
    ! 09Aug2022 B.Baker/NOAA    - Adapted for CCPP-Physics
    
    ! !Local Variables
    real(kind_chem)                  :: alpha_grav
    real(kind_chem)                  :: h
    real(kind_chem)                  :: kvh
    real(kind_chem)                  :: q
    real(kind_chem)                  :: rustar
    real(kind_chem)                  :: total_emissions
    real(kind_chem)                  :: u_sum, u_thresh
    real(kind_chem)                  :: fracland
    
!EOP
!-------------------------------------------------------------------------
!  Begin

!  Initialize emissions
!  --------------------
   emissions = 0.

!  Prepare scaling factor
!  ----------------------
   alpha_grav = alpha / grav

   fracland = max(0., min(1., 1.-flake)) &
                  * max(0., min(1., 1.-fsnow))

   ! Compute vertical-to-horizontal mass flux ratio
   ! ----------------------------------------------
   kvh = DustFluxV2HRatioMB95(clay, kvhmax)

   ! Compute total emissions
   ! -----------------------
   emissions = alpha_grav * (ssm ** gamma) * airdens * kvh * fracland

   !  Compute threshold wind friction velocity using drag partition
   !  -------------------------------------------------------------
   rustar = rdrag * ustar

   !  Now compute size-dependent total emission flux
   !  ----------------------------------------------
   ! Fecan moisture correction
   ! -------------------------
   h = moistureCorrectionFecan(slc, sand, clay, rhop)
   
   ! Adjust threshold
   ! ----------------
   u_thresh = uthrs * h
   
   u_sum = rustar + u_thresh
   
   ! Compute Horizontal Saltation Flux according to Eq (9) in Webb et al. (2020)
   ! ---------------------------------------------------------------------------
   q = max(0., rustar - u_thresh) * u_sum * u_sum
   
   ! Distribute emissions to bins and convert to mass flux (kg s-1)
   ! --------------------------------------------------------------
   emissions = emissions * q


 end subroutine DustEmissionFENGSHA
!-----------------------------------------------------------------
  real function soilMoistureConvertVol2Grav(vsoil, sandfrac, rhop)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(kind_chem), intent(in) :: vsoil       ! volumetric soil moisture fraction [1]
    REAL(kind_chem), intent(in) :: sandfrac    ! fractional sand content [1]
    REAL(kind_chem), intent(in) :: rhop        ! dry dust density [kg m-3]

! !DESCRIPTION: Convert soil moisture fraction from volumetric to gravimetric.
!
! !REVISION HISTORY:
!
!  02Apr2020, B.Baker/NOAA    - Original implementation
!  01Apr2020, R.Montuoro/NOAA - Adapted for GOCART process library

!  !Local Variables
    real :: vsat

!  !CONSTANTS:
    REAL(kind_chem), parameter :: rhow = 1000.    ! density of water [kg m-3]

!EOP
!-------------------------------------------------------------------------
!  Begin...

!  Saturated volumetric water content (sand-dependent) ! [m3 m-3]
    vsat = 0.489 - 0.00126 * ( 100. * sandfrac )

!  Gravimetric soil content
    soilMoistureConvertVol2Grav = vsoil * rhow / (rhop * (1. - vsat))

  end function soilMoistureConvertVol2Grav
!----------------------------------------------------------------
  real function moistureCorrectionFecan(slc, sand, clay, rhop)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(kind_chem), intent(in) :: slc     ! liquid water content of top soil layer, volumetric fraction [1]
    REAL(kind_chem), intent(in) :: sand    ! fractional sand content [1]
    REAL(kind_chem), intent(in) :: clay    ! fractional clay content [1]
    REAL(kind_chem), intent(in) :: rhop    ! dry dust density [kg m-3]

! !DESCRIPTION: Compute correction factor to account for Fecal soil moisture
!
! !REVISION HISTORY:
!
!  02Apr2020, B.Baker/NOAA    - Original implementation
!  01Apr2020, R.Montuoro/NOAA - Adapted for GOCART process library

!  !Local Variables
    real :: grvsoilm
    real :: drylimit

!EOP
!---------------------------------------------------------------
!  Begin...

!  Convert soil moisture from volumetric to gravimetric
    grvsoilm = soilMoistureConvertVol2Grav(slc, sand, 2650.)

!  Compute fecan dry limit
    drylimit = clay * (14.0 * clay + 17.0)

!  Compute soil moisture correction
    moistureCorrectionFecan = sqrt(1.0 + 1.21 * max(0., grvsoilm - drylimit)**0.68)

  end function moistureCorrectionFecan
!---------------------------------------------------------------
  real function DustFluxV2HRatioMB95(clay, kvhmax)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(kind_chem), intent(in) :: clay      ! fractional clay content [1]
    REAL(kind_chem), intent(in) :: kvhmax    ! maximum flux ratio [1]

!  !CONSTANTS:
    REAL(kind_chem), parameter :: clay_thresh = 0.2    ! clay fraction above which the maximum flux ratio is returned

! !DESCRIPTION: Computes the vertical-to-horizontal dust flux ratio according to
!               B.Marticorena, G.Bergametti, J.Geophys.Res., 100(D8), 164!               doi:10.1029/95JD00690
!
! !REVISION HISTORY:
!
! 22Feb2020 B.Baker/NOAA    - Original implementation
! 01Apr2021 R.Montuoro/NOAA - Adapted for GOCART process library
!
!EOP
!-------------------------------------------------------------------------
!  Begin...

    if (clay > clay_thresh) then
       DustFluxV2HRatioMB95 = kvhmax
    else
       DustFluxV2HRatioMB95 = 10.0**(13.4*clay-6.0)
    end if

  end function DustFluxV2HRatioMB95
  
end module dust_fengsha_mod
