module dust_fengsha_mod
!
!  This module developed by Barry Baker (NOAA ARL)
!  For serious questions contact barry.baker@noaa.gov
!
!  07/16/2019 - Adapted for NUOPC/GOCART, R. Montuoro
!  02/01/2020 - Adapted for FV3/CCPP, Haiqin Li
!  06/2023, Restructure for CATChem, Jian.He@noaa.gov
!  08/2023, Update to latest version Barry.Baker@noaa.gov
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
          chem_arr,rho_phy,dz8w,smois,         &
          delp,ssm,isltyp,vegfra,snowh,area,  &
          emis_dust,ust,znt,clay,sand, &
          rdrag,uthr,num_soil_layers,random_factor)

    IMPLICIT NONE

    ! Input Variables
    ! ---------------
    INTEGER, INTENT(IN) :: isltyp
    INTEGER, INTENT(IN) :: num_soil_layers ! number of soil layers 
    REAL(kind=kind_chem), INTENT(IN) :: dt ! timestep 
    REAL(kind=kind_chem), INTENT(IN) :: rho_phy ! air density
    REAL(kind=kind_chem), INTENT(IN) :: dz8w 
    REAL(kind=kind_chem), INTENT(IN) :: delp    ! del Pressure in first layer
    REAL(kind=kind_chem), INTENT(IN) :: ssm     ! sediment supply map 
    REAL(kind=kind_chem), INTENT(IN) :: vegfra  ! vegetation fraction 
    REAL(kind=kind_chem), INTENT(IN) :: snowh   ! snow height 
    REAL(kind=kind_chem), INTENT(IN) :: area    ! area of grid cell 
    REAL(kind=kind_chem), INTENT(IN) :: ust     ! ustar 
    REAL(kind=kind_chem), INTENT(IN) :: znt     ! surface roughness length 
    REAL(kind=kind_chem), INTENT(IN) :: clay    ! fractional clay content 
    REAL(kind=kind_chem), INTENT(IN) :: sand    ! fractional sand content 
    REAL(kind=kind_chem), INTENT(IN) :: rdrag   ! drag partition 
    REAL(kind=kind_chem), INTENT(IN) :: uthr    ! dry threshold friction velocity
    REAL(kind=kind_chem), INTENT(IN) :: random_factor
  
    REAL(kind=kind_chem), DIMENSION( num_chem ), INTENT(INOUT)                 :: chem_arr  ! chemical array 
    REAL(kind=kind_chem), DIMENSION( num_emis_dust ), OPTIONAL, INTENT(INOUT ) :: emis_dust ! final dust emission 
    REAL(kind=kind_chem), DIMENSION( num_soil_layers ), INTENT(IN)             :: smois     ! volumetric soil moisture at 0-5cm 

    ! Local variables
    ! ---------------
    integer :: do_dust ! 0 - no dust emission 1 - dust emission 
    integer :: n       ! looping variable 
    
    real(kind=kind_chem), DIMENSION (num_emis_dust) :: tc ! tracer concentration 
    REAL(kind=kind_chem), DIMENSION (num_emis_dust) :: bin_emis ! bin emisions 
    real(kind=kind_chem)                            :: emis ! total emission
    real(kind=kind_chem), DIMENSION (num_emis_dust) :: distribution ! fractional distribution for size bins 
    
    real(kind_chem), dimension (3)  :: massfrac      ! fractional mass of sand silt and clay 
    
    real(kind=kind_chem), PARAMETER :: conver=1.e-9 ! parameter to convert units   
    real(kind=kind_chem), PARAMETER :: converi=1.e9 ! parameter to convert units back 
    ! Total concentration at lowest model level. 
    do n=0,num_emis_dust-1 
       tc(n+1)=chem_arr(p_dust_1+n)*conver
    end do

    ! Air mass at lowest model level.
    airmas=area * delp / g

    ! ====================================
    ! Don't do dust over certain criteria 
    ! ====================================
    do_dust = 1
    
    ! limit where there is lots of vegetation
    if (vegfra .gt. .17) then
       do_dust = 0
    endif

    ! limit where there is snow on the ground
    if (snowh .gt. 0) then
       do_dust = 0
    endif

    ! Do not allow areas with bedrock, lava, or land-ice to loft
    IF (isltyp.eq. 15 .or. isltyp .eq. 16. .or. &
          isltyp .eq. 18) then
      do_dust=0
    ENDIF
    
    ! do not allow dust over the ocean 
    IF (isltyp .eq. 0)then
      do_dust=0
    endif
    
    ! check ssm input valid range
    if (ssm .lt. 0.05 .and. ssm .gt. 1.0)  then ! ensure values are realistic
       do_dust = 0
    end if
    
    ! check drag partition valid range 
    if (rdrag .lt. 0.05 .and. rdrag .gt. 1.0) then 
       do_dust = 0
    endif
    
    if(do_dust == 1 ) return
    ! ====================================

    ! Call fengsha dust emission for total emission 
    call DustEmissionFENGSHA(slc, clay, sand, ssm, rdrag, airdens, ustar, uthrs, area, dust_alpha, dust_gamma, emis)

    ! call dust distribution 
    call DustAerosolDistributionKok(reff_dust, lo_dust, up_dust, distribution)

    ! Distribute emissions to bins and convert to mass flux (kg s-1)
    ! --------------------------------------------------------------
    bin_emis = distribution * emis * dt * random_factor

    ! now convert to tracer concentration 
    ! -----------------------------------
    DO n=1,nmx
       dsrc = bin_emis(n)
       IF (dsrc < 0.0) dsrc = 0.0

       ! Update dust mixing ratio at first model level.
       tc(n) = tc(n) + dsrc / airmas ! (kg/kg)
       emis_dust(n)= 1.e+9*dsrc/(dxy*dt1) ! diagnostic (ug/m2/s) !lzhang
    END DO
    
    do n = 0, nmx-1
       chem_arr(p_dust_1+n))=tc(n + 1)*converi
    end do

  end subroutine gocart_dust_fengsha_driver


  SUBROUTINE source_dust(nmx, smx, dt1, tc, ustar, massfrac, &
       erod, dxy, gravsm, airden, airmas, bems, g0, drylimit, alpha,  &
       gamma, R, uthres, random_factor)

    ! ****************************************************************************
    ! *  Evaluate the source of each dust particles size bin by soil emission
    ! *
    ! *  Input:
    ! *         EROD      Fraction of erodible grid cell                (-)
    ! *         GRAVSM    Gravimetric soil moisture                     (g/g)
    ! *         DRYLIMIT  Upper GRAVSM limit for air-dry soil           (g/g)
    ! *         ALPHA     Constant to fudge the total emission of dust  (1/m)
    ! *         GAMMA     Tuning constant for erodibility               (-)
    ! *         DXY       Surface of each grid cell                     (m2)
    ! *         AIRMAS    Mass of air for each grid box                 (kg)
    ! *         AIRDEN    Density of air for each grid box              (kg/m3)
    ! *         USTAR     Friction velocity                             (m/s)
    ! *         DT1       Time step                                     (s)
    ! *         NMX       Number of dust bins                           (-)
    ! *         SMX       Number of saltation bins                      (-)
    ! *         IMX       Number of I points                            (-)
    ! *         JMX       Number of J points                            (-)
    ! *         LMX       Number of L points                            (-)
    ! *         R         Drag Partition                                (-)
    ! *         UTHRES    FENGSHA Dry Threshold Velocities              (m/s)
    ! *
    ! *  Data:
    ! *         MASSFRAC  Fraction of mass in each of 3 soil classes    (-)
    ! *         SPOINT    Pointer to 3 soil classes                     (-)
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
    ! *         G         Gravitational acceleration in cgs             (cm/s2)
    ! *
    ! *  Working:
    ! *         U_TS0     "Dry" threshold friction velocity             (m/s)
    ! *         U_TS      Moisture-adjusted threshold friction velocity (m/s)
    ! *         RHOA      Density of air in cgs                         (g/cm3)
    ! *         DEN       Dust density in cgs                           (g/cm3)
    ! *         DIAM      Dust diameter in cgs                          (cm)
    ! *         DMASS     Saltation mass distribution                   (-)
    ! *         DSURFACE  Saltation surface area per unit mass          (m2/kg)
    ! *         DS_REL    Saltation surface area distribution           (-)
    ! *         SALT      Saltation flux                                (kg/m/s)
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

    INTEGER,            INTENT(IN)    :: nmx,smx
    REAL(kind_chem), INTENT(IN)    :: dt1
    REAL(kind_chem), INTENT(INOUT) :: tc(nmx)
    REAL(kind_chem), INTENT(IN)    :: ustar
    REAL(kind_chem), INTENT(IN)    :: massfrac(3)
    REAL(kind_chem), INTENT(IN)    :: erod
    REAL(kind_chem), INTENT(IN)    :: dxy
    REAL(kind_chem), INTENT(IN)    :: gravsm
    REAL(kind_chem), INTENT(IN)    :: random_factor
    REAL(kind_chem), INTENT(IN)    :: airden
    REAL(kind_chem), INTENT(IN)    :: airmas
    REAL(kind_chem), INTENT(OUT)   :: bems(nmx)
    REAL(kind_chem), INTENT(IN)    :: g0
    REAL(kind_chem), INTENT(IN)    :: drylimit
    !! Sandblasting mass efficiency, aka "fudge factor" (based on Tegen et al,
    !! 2006 and Hemold et al, 2007)
    !
    !  REAL, PARAMETER :: alpha=1.8E-8  ! (m^-1)
    REAL(kind_chem), INTENT(IN)    :: alpha
    ! Experimental optional exponential tuning constant for erodibility.
    ! 0 < gamma < 1 -> more relative impact by low erodibility regions.
    REAL(kind_chem), INTENT(IN)    :: gamma
    REAL(kind_chem), INTENT(IN)    :: R
    REAL(kind_chem), INTENT(IN)    :: uthres

    REAL(kind_chem)    :: den(smx), diam(smx)
    REAL(kind_chem)    :: dvol(nmx), distr_dust(nmx), dlndp(nmx)
    REAL(kind_chem)    :: dsurface(smx), ds_rel(smx)
    REAL(kind_chem)    :: u_ts0, u_ts, dsrc, dmass, dvol_tot
    REAL(kind_chem)    :: salt,emit, emit_vol, stotal
    REAL(kind_chem)    :: rhoa, g
    INTEGER   :: i, j, n

    ! Sandblasting mass efficiency, beta.
    ! Beta maxes out for clay fractions above 0.2 = betamax.

    REAL(kind_chem), PARAMETER :: betamax=5.25E-4
    REAL(kind_chem) :: beta
    integer :: styp

    ! Constant of proportionality from Marticorena et al, 1997 (unitless)
    ! Arguably more ~consistent~ fudge than alpha, which has many walnuts
    ! sprinkled throughout the literature. - GC

    REAL(kind_chem), PARAMETER :: cmb=1.0
    ! REAL, PARAMETER :: cmb=2.61   ! from White,1979

    ! Parameters used in Kok distribution function. Advise not to play with
    ! these without the expressed written consent of someone who knows what
    ! they're doing. - GC

    REAL(kind_chem), PARAMETER :: mmd_dust=3.4D-6  ! median mass diameter (m)
    REAL(kind_chem), PARAMETER :: gsd_dust=3.0     ! geom. std deviation
    REAL(kind_chem), PARAMETER :: lambda=12.0D-6   ! crack propagation length (m)
    REAL(kind_chem), PARAMETER :: cv=12.62D-6      ! normalization constant

    ! Calculate saltation surface area distribution from sand, silt, and clay
    ! mass fractions and saltation bin fraction. This will later become a
    ! modifier to the total saltation flux.  The reasoning here is that the
    ! size and availability of saltators affects saltation efficiency. Based
    ! on Eqn. (32) in Marticorena & Bergametti, 1995 (hereon, MB95).

    DO n=1,smx
       dmass=massfrac(spoint(n))*frac_salt(n)
       dsurface(n)=0.75*dmass/(den_salt(n)*reff_salt(n))
    ENDDO

    ! The following equation yields relative surface area fraction.  It will only
    ! work if you are representing the "full range" of all three soil classes.
    ! For this reason alone, we have incorporated particle sizes that encompass
    ! the clay class, to account for the its relative area over the basal
    ! surface, even though these smaller bins would be unlikely to play any large
    ! role in the actual saltation process. - GC

    stotal=SUM(dsurface(:))
    DO n=1,smx
       ds_rel(n)=dsurface(n)/stotal
    ENDDO

    ! Calculate total dust emission due to saltation of sand sized particles.
    ! Begin by calculating DRY threshold friction velocity (u_ts0).  Next adjust
    ! u_ts0 for moisture to get threshold friction velocity (u_ts). Then
    ! calculate saltation flux (salt) where ustar has exceeded u_ts.  Finally,
    ! calculate total dust emission (tot_emit), taking into account erodibility.

    ! Set DRY threshold friction velocity to input value
    u_ts0 = uthres

    g = g0*1.0E2
    emit=0.0

    DO n = 1, smx
       den(n) = den_salt(n)*1.0D-3         ! (g cm^-3)
       diam(n) = 2.0*reff_salt(n)*1.0D2    ! (cm)
       rhoa = airden*1.0D-3                ! (g cm^-3)

             ! FENGSHA uses the 13 category soil type from the USDA
             ! call calc_fengsha_styp(massfrac(1),massfrac(3),massfrac(2),styp)
             ! Fengsha uses threshold velocities based on dale gilletes data
             ! call fengsha_utst(styp,uthres,u_ts0)

             ! Friction velocity threshold correction function based on physical
             ! properties related to moisture tension. Soil moisture greater than
             ! dry limit serves to increase threshold friction velocity (making
             ! it more difficult to loft dust). When soil moisture has not reached
             ! dry limit, treat as dry

             IF (gravsm > drylimit) THEN
                u_ts = MAX(0.0D+0,u_ts0*(sqrt(1.0+1.21*(gravsm-drylimit)**0.68)) / R)
             ELSE
                u_ts = u_ts0 / R
             END IF

             ! Calculate total vertical mass flux (note beta has units of m^-1)
             ! Beta acts to tone down dust in areas with so few dust-sized particles that the
             ! lofting efficiency decreases.  Otherwise, super sandy zones would be huge dust
             ! producers, which is generally not the case.  Equation derived from wind-tunnel
             ! experiments (see MB95).

             beta=10**(13.6*massfrac(1)-6.0)  ! (unitless)
             if (massfrac(1) <= 0.2) then
                beta=10**(13.4*massfrac(1)-6.0)
             else
                beta = 2.E-4
             endif

             !---------------------------------------------------------------------
             ! formula of Draxler & Gillette (2001) Atmos. Environ.
             ! F   =  K A (r/g) U* ( U*^2 - Ut*^2 )
             !
             ! where:
             !     F   = vertical emission flux  [g/m**2-s]
             !     K   = constant 2.0E-04                      [1/m]
             !     A   = 0~3.5  mean = 2.8  (fudge factor)
             !     U*  = friction velocity                     [m/s]
             !     Ut* = threshold friction velocity           [m/s]
             !
             !--------------------------------------------------------------------

             IF (ustar .gt. u_ts) then
                call fengsha_hflux(ustar,u_ts,beta, salt)
                salt = alpha * cmb * ds_rel(n) * airden / g0 * salt * (erod**gamma) * beta
             else
                salt = 0.
             endif
             ! EROD is taken into account above
             emit = emit + salt 
    END DO

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
             dsrc = emit*distr_dust(n)*dxy*dt1 *random_factor  ! (kg)
             IF (dsrc < 0.0) dsrc = 0.0

             ! Update dust mixing ratio at first model level.
             tc(n) = tc(n) + dsrc / airmas ! (kg/kg)
             !   bems(i,j,n) = dsrc  ! diagnostic
             !bems(i,j,n) = 1000.*dsrc/(dxy(j)*dt1) ! diagnostic (g/m2/s)
             bems(n) = 1.e+9*dsrc/(dxy*dt1) ! diagnostic (ug/m2/s) !lzhang
    END DO

  END SUBROUTINE source_dust

  subroutine DustEmissionFENGSHA(slc, clay, ssm, rdrag, airdens, ustar, uthrs, area, ,alpha, gamma, emissions)

      ! !USES:
      implicit NONE

      ! !INPUT PARAMETERS:
      real(kind=kind_chem), intent(in) :: slc      ! liquid water content of soil layer, volumetric fraction [1]
      real(kind=kind_chem), intent(in) :: clay     ! fractional clay content [1] - range: [0 1]
      real(kind=kind_chem), intent(in) :: sand     ! fractional clay content [1] - range: [0 1]
      real(kind=kind_chem), intent(in) :: ssm      ! erosion map [1] - range: [0 1]
      real(kind=kind_chem), intent(in) :: rdrag    ! drag partition [1/m] - range: [0 1]
      real(kind=kind_chem), intent(in) :: airdens  ! air density at lowest level [kg/m^3]
      real(kind=kind_chem), intent(in) :: ustar    ! friction velocity [m/sec]
      real(kind=kind_chem), intent(in) :: uthrs    ! threshold velocity [m/2]
      real(kind=kind_chem), intent(in) :: kvhmax   ! max. vertical to horizontal mass flux ratio [1]
      real(kind=kind_chem), intent(in) :: area     ! area of dust emission (can be fractional area of grid cell)

      ! this will need to be read in by a configuration file
      ! right now I just hard coded it in dust_data_mod.F90
      !=====================================================
      real(kind=kind_chem), intent(in) :: alpha    ! scaling factor [1]
      real(kind=kind_chem), intent(in) :: gamma    ! scaling factor [1]
      !=====================================================

      ! !OUTPUT PARAMETERS:
      REAL(kind=kind_chem), dimension(:), intent(inout) :: total_emissions ! binned surface emissions [kg/(m^2 sec)]

      ! !DESCRIPTION: Compute dust emissions using NOAA/ARL FENGSHA model
      !
      ! !REVISION HISTORY:
      !
      ! 22Feb2020 B.Baker/NOAA    - Original implementation
      ! 29Mar2021 R.Montuoro/NOAA - Refactored for process library
      ! 09Aug2022 B.Baker/NOAA    - Adapted for CCPP-Physics

      ! !Local Variables
      real(kind=kind_chem) :: alpha_grav
      real(kind=kind_chem) :: h
      real(kind=kind_chem) :: kvh
      real(kind=kind_chem) :: q
      real(kind=kind_chem) :: rustar
      real(kind=kind_chem) :: u_sum
      real(kind=kind_chem) :: u_thresh
      real(kind=kind_chem) :: emission
      real(kind=kind_chem) :: stotal
      real(kind=kind_chem) :: dmass
      real(kind=kind_chem), dimension(nsalt) :: dsurface
      real(kind=kind_chem), dimension(3) :: massfrac ! fractional soil content of sand(1) silt(2) and clay(3)

      real(kind=kind_chem), parameter:: clay_thresh = 0.2
      real(kind=kind_chem), parameter :: rhow = 1000.

      !EOP
      !-------------------------------------------------------------------------
      !  Begin

      !  Initialize emissions
      !  --------------------
      total_emissions = 0.

      !  Prepare scaling factor
      !  ----------------------
      alpha_grav = dust_alpha / con_g

      ! Compute vertical-to-horizontal mass flux ratio
      ! ----------------------------------------------
      if (clay > clay_thresh) then
         kvh = kvhmax
      else
         kvh = 10.0**(13.4*clay-6.0)
      end if

      ! Compute total emissions
      ! -----------------------
      emission = alpha_grav * (ssm ** dust_gamma) * airdens * kvh

      !  Compute threshold wind friction velocity using drag partition
      !  -------------------------------------------------------------
      rustar = rdrag * ustar

      !  Now compute size-dependent total emission flux
      !  ----------------------------------------------
      ! Fecan moisture correction
      ! -------------------------
      if (DUST_OPT_FENGSHA_FECAN .eq. .true.) then
         !  Saturated volumetric water content (sand-dependent) ! [m3 m-3]
         vsat = 0.489 - 0.00126 * ( 100. * sand )

         !  Gravimetric soil content
         grvsoilm = vsoil * rhow / (dust_den * (1. - vsat))

         !   Compute fecan dry limit
         drylimit = clay * (14.0 * clay + 17.0)

         !  Compute soil moisture correction
         h = sqrt(1.0 + 1.21 * max(0., grvsoilm - drylimit)**0.68)
      else
         ! Shao soil mositure
         !--------------------
         if (slc <= 0.03) then
            h = exp(22.7 * slc)
         else
            h = exp(93.5 * slc - 2.029)
         end if
      end if

      ! Adjust threshold
      ! ----------------
      u_thresh = uthrs * h

      u_sum = rustar + u_thresh

      ! Compute Horizontal Saltation Flux according to Eq (9) in Webb et al. (2020)
      ! ---------------------------------------------------------------------------
      q = max(0., rustar - u_thresh) * u_sum * u_sum

      ! Calculate total dust using the dust potential (q) 
      ! -------------------------------------------------
      emission = emission * q 
      
      ! Calculate saltation surface area distribution from sand, silt, and clay
      ! mass fractions and saltation bin fraction. Based on Eqn. (32) in 
      ! Marticorena & Bergametti, 1995 (hereon, MB95).
      ! ----------------------------------------------
      DO n=1,nsalt
         dmass=massfrac(spoint(n))*frac_salt(n)
         dsurface(n)=0.75*dmass/(den_salt(n)*reff_salt(n))
      ENDDO

      ! The following equation yields relative surface area fraction.  
      ! ------------------------------------------------------------
      stotal=SUM(dsurface(:))
      DO n=1,nsalt
         total_emissions = total_emissions + emission * dsurface(n)/stotal
      ENDDO
      
      ! Distribute emissions to bins and convert to mass flux (kg s-1)
      ! --------------------------------------------------------------
      !emissions = distribution * total_emissions * q


   end subroutine DustEmissionFENGSHA
   !-----------------------------------------------------------------
   subroutine DustAerosolDistributionKok ( radius, rLow, rUp, distribution )
     
     ! !USES:
     implicit NONE
     
     ! !INPUT PARAMETERS:
     real, dimension(:), intent(in)  :: radius      ! Dry particle bin effective radius [um]
     real, dimension(:), intent(in)  :: rLow, rUp   ! Dry particle bin edge radii [um]
     
     ! !OUTPUT PARAMETERS:
     real, dimension(:), intent(out) :: distribution    ! Normalized dust aerosol distribution [1]
     
     ! !DESCRIPTION: Computes lognormal aerosol size distribution for dust bins according to
     !               J.F.Kok, PNAS, Jan 2011, 108 (3) 1016-1021; doi:10.1073/pnas.1014798108
     !
     ! !REVISION HISTORY:
     !
     ! 22Feb2020 B.Baker/NOAA    - Original implementation
     ! 01Apr2021 R.Montuoro/NOAA - Refactored for GOCART process library
     !
     
     ! !Local Variables
     integer :: n, nbins
     real    :: diameter, dlam, dvol
     
     !   !CONSTANTS
     real, parameter    :: mmd    = 3.4          ! median mass diameter [um]
     real, parameter    :: stddev = 3.0          ! geometric standard deviation [1]
     real, parameter    :: lambda = 12.0         ! crack propagation length [um]
     real, parameter    :: factor = 1.e0 / (sqrt(2.e0) * log(stddev))  ! auxiliary constant
     
     character(len=*), parameter :: myname = 'DustAerosolDistributionKok'
     
     !EOP
     !-------------------------------------------------------------------------
     !  Begin...
     
     distribution = 0.
     
     !  Assume all arrays are dimensioned consistently
     nbins = size(radius)
     
     dvol = 0.
     do n = 1, nbins
        diameter = 2 * radius(n)
        dlam = diameter/lambda
        distribution(n) = diameter * (1. + erf(factor * log(diameter/mmd))) * exp(-dlam * dlam * dlam) * log(rUp(n)/rLow(n))
        dvol = dvol + distribution(n)
     end do
     
     !  Normalize distribution
     do n = 1, nbins
        distribution(n) = distribution(n) / dvol
     end do
     
   end subroutine DustAerosolDistributionKok
   
end module dust_fengsha_mod
