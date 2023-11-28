module dust_afwa_mod
!
!  this module developed by Sandra Jones (AFWA and AER) 
!  and Glenn Creighton (AFWA). For serious questions contact
!
!  08/31/2017 - Adapted for NUOPC/GOCART, R. Montuoro
!  06/2023, Restructure for CATChem, Jian.He@noaa.gov

  use catchem_constants, only : kind_chem, g=>con_g, pi=>con_pi
  use catchem_config, only: num_chem,num_emis_dust,p_dust_1, p_edust1
  use dust_data_mod

  implicit none

  private

  public :: gocart_dust_afwa_driver

contains

  subroutine gocart_dust_afwa_driver(ktau,dt,u_phy,              &
          v_phy,chem_arr,rho_phy,dz8w,smois,u10,         &
          v10,delp,erod,isltyp,area,       &
          emis_dust,ust,znt,clay,sand, &
          num_soil_layers)
    
    IMPLICIT NONE
    
    
    ! Input Variables 
    ! ---------------
    INTEGER, INTENT( IN ) :: ktau
    INTEGER, INTENT( IN ) :: isltyp
    INTEGER, INTENT( IN ) :: num_soil_layers
    
    REAL(kind=kind_chem), INTENT(IN) :: dt
    REAL(kind=kind_chem), INTENT(IN) :: u_phy
    REAL(kind=kind_chem), INTENT(IN) :: v_phy
    REAL(kind=kind_chem), INTENT(IN) :: rho_phy
    REAL(kind=kind_chem), INTENT(IN) :: dz8w
    REAL(kind=kind_chem), INTENT(IN) :: u10
    REAL(kind=kind_chem), INTENT(IN) :: v10
    REAL(kind=kind_chem), INTENT(IN) :: delp
    REAL(kind=kind_chem), INTENT(IN) :: area
    REAL(kind=kind_chem), INTENT(IN) :: ust
    REAL(kind=kind_chem), INTENT(IN) :: znt
    REAL(kind=kind_chem), INTENT(IN) :: clay
    REAL(kind=kind_chem), INTENT(IN) :: sand

    REAL(kind=kind_chem), DIMENSION(num_chem), INTENT(INOUT) ::  chem_arr
    REAL(kind=kind_chem), DIMENSION(ndust), INTENT(INOUT) :: emis_dust
    REAL(kind=kind_chem), DIMENSION(num_soil_layers), INTENT(INOUT) :: smois
    REAL(kind=kind_chem), DIMENSION( 3 ), INTENT(IN) :: erod

    ! Local variables

    integer :: n, nmx,smx,ilwi
    
    real(kind_chem) :: erodtot, 
    real(kind_chem) :: gravsm
    real(kind_chem) :: drylimit
    real(kind_chem) :: airden
    real(kind_chem) :: airmas
    real(kind_chem) :: ustar
    real(kind_chem) :: dxy
    real(kind_chem) :: conver
    real(kind_chem) :: converi
    
    real(kind_chem), dimension (3) :: massfrac
    
    real(kind_chem), DIMENSION (ndust) :: tc
    real(kind_chem), DIMENSION (ndust) :: bems
    
    

    conver=1.e-9
    converi=1.e9
    
    ! Number of dust bins
    
    nmx=ndust
    smx=nsalt
    
    ! Don't do dust over water!!!
    
    ilwi=1
    
    ! set initial tendency to zero
    tc = 0.
        
    ! Air mass and density at lowest model level.
    
    airmas=area * delp / g
    airden=rho_phy
    ustar=ust
    dxy=area
    
    ! Total erodibility.
    
    erodtot=SUM(erod(:))
    
    ! Mass fractions of clay, silt, and sand.
    
    massfrac(1)=clay
    massfrac(2)=1-(clay+sand)
    massfrac(3)=sand
    
    ! Don't allow roughness lengths greater than 20 cm to be lofted.
    ! This kludge accounts for land use types like urban areas and
    ! forests which would otherwise show up as high dust emitters.
    ! This is a placeholder for a more widely accepted kludge
    ! factor in the literature, which reduces lofting for rough areas.
    ! Forthcoming...
    
    IF (znt .gt. 0.2) then
       ilwi=0
    ENDIF
    
    ! Do not allow areas with bedrock, lava, or land-ice to loft
    ! ----------------------------------------------------------
    IF (isltyp.eq. 15 .or. isltyp .eq. 16. .or. isltyp .eq. 18) then
       ilwi=0
    ENDIF
    IF (isltyp .eq. 0)then
       ilwi=0
    endif

    if(ilwi == 0 ) return
    
    ! Calculate gravimetric soil moisture and drylimit.
    
    gravsm=100.*smois(1)/((1.-maxsmc(isltyp))*(2.65*(1.-clay)+2.50*clay))
    drylimit=14.0*clay*clay+17.0*clay
    
    ! Call dust emission routine.
    call source_dust(nmx, smx, dt, tc, ustar, massfrac, erodtot, ilwi, dxy, gravsm, airden, airmas, bems, g, drylimit, dust_alpha, dust_gamma)
    
    do n = 0, ndust-1
       ! for output diagnostics 
       ! ----------------------
       emis_dust(p_edust1 + n) = bems(n+1)
    end do
  
end subroutine gocart_dust_afwa_driver

  
  SUBROUTINE source_dust(nmx, smx, dt1, tc, ustar, massfrac, erod, ilwi, dxy, gravsm, airden, airmas, bems, g0, drylimit, alpha, gamma)
    
    ! ****************************************************************************
    ! *  Evaluate the source of each dust particles size bin by soil emission  
    ! *
    ! *  Input:
    ! *         EROD      Fraction of erodible grid cell                (-)
    ! *         ILWI      Land/water flag                               (-)
    ! *         GRAVSM    Gravimetric soil moisture                     (g/g)
    ! *         DRYLIMIT  Upper GRAVSM limit for air-dry soil           (g/g)
    ! *         ALPHA     Constant to fudge the total emission of dust  (1/m)
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
    ! *         LAMBDA    Side crack propogation length                 (m)
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

    INTEGER, INTENT(IN)   :: nmx,smx
    INTEGER, INTENT(IN)   :: ilwi
    REAL(kind_chem), INTENT(IN)    :: erod
    REAL(kind_chem), INTENT(IN)    :: ustar
    REAL(kind_chem), INTENT(IN)    :: gravsm
    REAL(kind_chem), INTENT(IN)    :: drylimit
    REAL(kind_chem), INTENT(IN)    :: dxy
    REAL(kind_chem), INTENT(IN)    :: airden
    REAL(kind_chem), INTENT(IN)    :: airmas
    REAL(kind_chem), INTENT(IN)    :: g0
    REAL(kind_chem), INTENT(IN)    :: dt1
    REAL(kind_chem), INTENT(INOUT) :: tc(nmx)
    REAL(kind_chem), INTENT(OUT)   :: bems(nmx) 
    
    ! Local Variables 
    ! ---------------
    REAL(kind_chem), dimension(smx) :: den
    REAL(kind_chem), dimension(smx) :: diam
    REAL(kind_chem), dimension(smx) :: dsurface
    REAL(kind_chem), dimension(smx) :: ds_rel
    REAL(kind_chem), dimension(nmx) :: dvol
    REAL(kind_chem), dimension(nmx) :: distr_dust
    REAL(kind_chem), dimension(nmx)  :: dlndp
  
    REAL(kind_chem), dimension(3) :: massfrac

    REAL(kind_chem) :: u_ts0
    REAL(kind_chem) :: u_ts
    REAL(kind_chem) :: dsrc
    REAL(kind_chem) :: srce
    REAL(kind_chem) :: dmass
    REAL(kind_chem) :: dvol_tot
    REAL(kind_chem) :: salt
    REAL(kind_chem) :: emit
    REAL(kind_chem) :: emit_vol
    REAL(kind_chem) :: stotal
    REAL(kind_chem) :: rhoa 
    REAL(kind_chem) :: g 
    
    INTEGER :: i, j, m, s, n

    ! Global tuning constant, alpha.  Sandblasting mass efficiency, beta.
    ! Beta maxes out for clay fractions above 0.2 = betamax.

    REAL(kind_chem), INTENT(IN)  :: alpha
    REAL(kind_chem), PARAMETER :: betamax=5.25E-4
    REAL(kind_chem) :: beta

    ! Experimental optional exponential tuning constant for erodibility.
    ! 0 < gamma < 1 -> more relative impact by low erodibility regions.
    
    REAL(kind_chem), INTENT(IN) :: gamma
    
    ! Constant of proportionality from Marticorena et al, 1997 (unitless)
    ! Arguably more ~consistent~ fudge than alpha, which has many walnuts
    ! sprinkled throughout the literature. - GC
    
    REAL(kind_chem), PARAMETER :: cmb=1.0    

    ! Parameters used in Kok distribution function. Advise not to play with 
    ! these without the expressed written consent of someone who knows what
    ! they're doing. - GC
    
    REAL(kind_chem), PARAMETER :: mmd_dust=3.4D-6  ! median mass diameter (m)
    REAL(kind_chem), PARAMETER :: gsd_dust=3.0     ! geom. std deviation
    REAL(kind_chem), PARAMETER :: lambda=12.0D-6   ! crack propogation length (m)
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
    
    g = g0*1.0E2
    emit=0.0
    
    DO n = 1, smx
       den(n) = den_salt(n)*1.0D-3         ! (g cm^-3)
       diam(n) = 2.0*reff_salt(n)*1.0D2    ! (cm)
       rhoa = airden*1.0D-3       ! (g cm^-3)
       
       ! Threshold friction velocity as a function of the dust density and
       ! diameter from Bagnold (1941) (m s^-1).
       
       u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
            SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
            SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0) 
       
       ! Friction velocity threshold correction function based on physical
       ! properties related to moisture tension. Soil moisture greater than
       ! dry limit serves to increase threshold friction velocity (making
       ! it more difficult to loft dust). When soil moisture has not reached
       ! dry limit, treat as dry (no correction to threshold friction
       ! velocity). GC
       
       IF (gravsm > drylimit) THEN
          u_ts = MAX(0.0D+0,u_ts0*(sqrt(1.0+1.21*(gravsm-drylimit)**0.68)))
       ELSE
          u_ts = u_ts0
       END IF
       
       ! Saltation flux from Marticorena & Bergametti 1995 (MB95). ds_rel is
       ! the relative surface area distribution
       
       IF (ustar .gt. u_ts .and. erod .gt. 0.0 .and. ilwi == 1) THEN
          salt = cmb*ds_rel(n)*(airden/g0)*(ustar**3)* &
               (1. + u_ts/ustar)*(1. - (u_ts**2)/(ustar**2))  ! (kg m^-1 s^-1)
       ELSE 
          salt = 0.0
       ENDIF

       ! Calculate total vertical mass flux (note beta has units of m^-1)
       ! Beta acts to tone down dust in areas with so few dust-sized particles that the
       ! lofting efficiency decreases.  Otherwise, super sandy zones would be huge dust
       ! producers, which is generally not the case.  Equation derived from wind-tunnel 
       ! experiments (see MB95).
       
       beta=10**(13.6*massfrac(1)-6.0)  ! (unitless)
       if (beta .gt. betamax) then
          beta=betamax
       endif
       emit=emit+salt*(erod**gamma)*alpha*beta    ! (kg m^-2 s^-1)
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
       ! emit_vol=emit/den_dust(n) ! (m s^-1)
    END DO
    DO n=1,nmx
       distr_dust(n)=dvol(n)/dvol_tot
    END DO
    
    ! Now distribute total vertical emission into dust bins and update concentration.
    
    DO n=1,nmx 
       ! Calculate total mass emitted
       dsrc = emit*distr_dust(n)*dxy*dt1  ! (kg)
       IF (dsrc < 0.0) dsrc = 0.0
       
       ! Update dust mixing ratio at first model level.
       tc(n) = tc(n) + dsrc / airmas ! (kg/kg)
       !bems(i,j,n) = dsrc/(dxy(j)*dt1) ! diagnostic (kg/m2/s)
       bems(n) = 1.e+9*dsrc/(dxy*dt1) ! diagnostic (ug/m2/s)
    END DO
    
  END SUBROUTINE source_dust
  
end module dust_afwa_mod
