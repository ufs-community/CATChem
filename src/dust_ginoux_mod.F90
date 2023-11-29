!! Revision History:
!! 06/2023, Restructure for CATChem, Jian.He@noaa.gov

module dust_ginoux_mod
  
  use catchem_constants, only : kind_chem, g=>con_g, pi=>con_pi
  use catchem_config, only: num_chem,num_emis_dust,p_dust_1,p_edust1
  use dust_data_mod
  
  implicit none
  
  private
  
  public :: gocart_dust_default
  
CONTAINS
  
  subroutine gocart_dust_ginoux(ktau,dt,u_phy, v_phy,rho_phy,dz8w,smois,u10,v10,delp,erod,isltyp,area,emis_dust,srce_dust,num_soil_layers,start_month)

    ! Input Variables 
    ! ---------------
    INTEGER, INTENT(IN) :: ktau
    INTEGER, INTENT(IN) :: isltyp
    INTEGER, INTENT(IN) :: start_month
    INTEGER, INTENT(IN) :: num_soil_layers
    
    
    REAL(kind=kind_chem), INTENT(IN) :: dt
    REAL(kind=kind_chem), INTENT(IN) :: u_phy
    REAL(kind=kind_chem), INTENT(IN) :: v_phy
    REAL(kind=kind_chem), INTENT(IN) :: rho_phy
    REAL(kind=kind_chem), INTENT(IN) :: dz8w
    REAL(kind=kind_chem), INTENT(IN) :: u10
    REAL(kind=kind_chem), INTENT(IN) :: v10
    REAL(kind=kind_chem), INTENT(IN) :: delp
    REAL(kind=kind_chem), INTENT(IN) :: area
    
    REAL(kind=kind_chem), DIMENSION(3), INTENT(IN) :: erod
    
    ! Output Variables 
    ! ----------------
    REAL(kind=kind_chem), DIMENSION(ndust), INTENT(INOUT ) :: emis_dust
    REAL(kind=kind_chem), DIMENSION(ndust), INTENT(INOUT ) :: srce_dust
    REAL(kind=kind_chem), DIMENSION(num_soil_layers), INTENT(INOUT) :: smois
    
    ! Local Variables 
    ! ---------------
    integer :: ipr,ilwi, n
    real(kind=kind_chem), DIMENSION (3) :: erodin  ! (ndcls,ndsrc)
    real(kind=kind_chem), DIMENSION (ndust) :: bems
    real(kind=kind_chem), DIMENSION (ndust) :: srce_out
    real(kind=kind_chem) :: w10m
    real(kind=kind_chem) :: gwet
    real(kind=kind_chem) :: airden
    real(kind=kind_chem) :: airmas
    real(kind=kind_chem) :: dxy
    real(kind=kind_chem), parameter :: min_default=0.
    real(kind=kind_chem), parameter :: conver  = 1.e-9
    real(kind=kind_chem), parameter :: converi = 1.e+9
    
    ilwi=1
    
    do n = 0, ndust -1 
       tc(n+1)=chem_arr(p_dust_1+ n)*conver
    end do
    
    ! Calculate 10m wind speed from u and v
    ! -------------------------------------
    w10m=sqrt(u10*u10+v10*v10)
    
    ! calculate air mass in first layer
    ! --------------------------------
    airmas=area * delp / g
    
    ! Calcualte 10m wind speed
    ! ------------------------
    if(dz8w.lt.12.) then
       ! Calculate 10m wind speed from first model layer if 
       ! first layer is thin
       ! --------------------------------------------------
       w10m=sqrt(u_phy*u_phy+v_phy*v_phy)
    else
       ! Calculate 10m wind speed from u10 and v10 
       ! -------------------------------------
       w10m=sqrt(u10*u10+v10*v10)
    end if
    
    ! get erosion potentials 
    ! ----------------------
    erodin(1)=erod(1)
    erodin(2)=erod(2)
    erodin(3)=erod(3)
    
    ! volumetric soil moisture over porosity
    ! --------------------------------------
    if(isltyp.eq.0)then
       ilwi=0
       gwet=1.
    else
       gwet=smois(1)/maxsmc(isltyp)
    endif

    ipr=0

    call source_du( ndust, dt, erodin, ilwi, area, w10m, gwet, rho_phy, airmas, bems,srce_out, start_month, g, ipr)

    do n = 0, ndust-1 
      !  ! Update Tracer concentrations 
      !  ! ----------------------------
      !  chem_arr(p_dust_1 + n)=max(min_default,tc(n + 1)*converi)

       ! update diagnostic emission rate
       ! -------------------------------
       emis_dust(p_edust1 + n) = max(0, bems(n+1) * converi)

       ! for output diagnostics of dust source
       ! -------------------------------------
       srce_dust(p_edust1 + n)=srce_out(n+1)
    end do

  end subroutine gocart_dust_default


  SUBROUTINE source_du( nmx, dt1,  & 
       erod, ilwi, dxy, w10m, gwet, airden, airmas, &
       bems,srce_out,month,g0,ipr) 

    ! ****************************************************************************
    ! *  Evaluate the source of each dust particles size classes  (kg/m3)        
    ! *  by soil emission.
    ! *  Input:
    ! *         EROD      Fraction of erodible grid cell                (-)
    ! *                   for 1: Sand, 2: Silt, 3: Clay
    ! *         DUSTDEN   Dust density                                  (kg/m3)
    ! *         DXY       Surface of each grid cell                     (m2)
    ! *         AIRVOL    Volume occupy by each grid boxes              (m3)
    ! *         NDT1      Time step                                     (s)
    ! *         W10m      Velocity at the anemometer level (10meters)   (m/s)
    ! *         u_tresh   Threshold velocity for particule uplifting    (m/s)
    ! *         CH_dust   Constant to fudge the total emission of dust  (s2/m2)
    ! *      
    ! *  Output:
    ! *         DSRC      Source of each dust type           (kg/timestep/cell) 
    ! *
    ! *  Working:
    ! *         SRC       Potential source                   (kg/m/timestep/cell)
    ! *
    ! ****************************************************************************

    INTEGER,            INTENT(IN)    :: nmx,ilwi,month

    REAL(kind=kind_chem), INTENT(IN)    :: dt1
    REAL(kind=kind_chem), INTENT(IN)    :: g0 
    REAL(kind=kind_chem), INTENT(IN)    :: erod(ndcls,ndsrc)
    REAL(kind=kind_chem), INTENT(IN)    :: w10m
    REAL(kind=kind_chem), INTENT(IN)    :: gwet
    REAL(kind=kind_chem), INTENT(IN)    :: dxy
    REAL(kind=kind_chem), INTENT(IN)    :: airden
    REAL(kind=kind_chem), INTENT(IN)    :: airmas
    
    !----------------------------------------------------------------------
    ! Output 
    !----------------------------------------------------------------------
   !  REAL(kind=kind_chem), INTENT(INOUT) :: tc(nmx)
    REAL(kind=kind_chem), INTENT(OUT)   :: bems(nmx)
    REAL(kind=kind_chem), INTENT(OUT)   :: srce_out(nmx) !dust source
    INTEGER,            INTENT(OUT)   :: ipr

    !-----------------------------------------------------------------------  
    ! local variables
    !-----------------------------------------------------------------------  
    INTEGER            :: i, j, n, m, k
   !  REAL(kind=kind_chem) :: g ! redundant variable 
    REAL(kind=kind_chem) :: den(nmx), diam(nmx)
    REAL(kind=kind_chem) :: rhoa, tsrc, u_ts0, u_ts, dsrc, srce
    LOGICAL            :: aerosol(nmx)

    REAL(kind=kind_chem), PARAMETER :: gthresh = 0.5_kind_chem

    ! executable statemenst
    ipr = 0

    DO n = 1, nmx
       ! Threshold velocity as a function of the dust density and the diameter
       ! from Bagnold (1941)
       den(n) = den_dust(n)*1.0D-3
       diam(n) = 2.0*reff_dust(n)*1.0D2
       g = g0*1.0E2
       ! Pointer to the 3 classes considered in the source data files
       m = ipoint(n)
       tsrc = 0.0_kind_chem
       
       ! No flux if wet soil 
       rhoa = airden*1.0D-3
       u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
               SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0) 

       ! Case of surface dry enough to erode
       IF (gwet < gthresh) THEN
          u_ts = MAX(0.0D+0,u_ts0*(1.2D+0+2.0D-1*LOG10(MAX(1.0D-3, gwet))))
       ELSE
          ! Case of wet surface, no erosion
          u_ts = 100.0_kind_chem
       END IF
       DO k =1, ndsrc
          srce = frac_s(n)*erod(m,k)*dxy  ! (m2)
          srce_out (k)=srce  ! output dust source
          IF (ilwi == 1 ) THEN
             dsrc = ch_dust(n,month)*srce*w10m**2 &
                  * (w10m - u_ts)*dt1  ! (kg)
          ELSE 
             dsrc = 0.0_kind_chem
          END IF
          IF (dsrc < 0.0_kind_chem) dsrc = 0.0_kind_chem

          ! Update dust mixing ratio at first model level.
          ! scale down dust by .7
         !  tc(n) = tc(n) + .7*dsrc / airmas
          bems(n) = .7*dsrc/(dxy*dt1) ! diagnostic (kg/m2/s)
       END DO
    END DO

  END SUBROUTINE source_du

end module gocart_dust_default_mod
