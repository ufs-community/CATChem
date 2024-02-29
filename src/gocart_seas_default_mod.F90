!! Revision History:
!! 06/2023, Restructure for CATChem, Jian.He@noaa.gov

module gocart_seas_default_mod

  use catchem_constants, only : kind_chem, g=>con_g, pi=>con_pi
  use catchem_config, only : num_emis_seas,num_chem, &
                             p_seas_1,p_seas_2,p_seas_3,p_seas_4,p_seas_5, &
                             p_eseas1,p_eseas2,p_eseas3,p_eseas4,p_eseas5
  use seas_data_mod

  implicit none

  private

  public :: gocart_seas_default

CONTAINS

  subroutine gocart_seas_default(ktau,dt,u_phy,v_phy,chem_arr,dz8w,u10,         &
       v10,delp,tsk,area,emis_seas)
    
    ! Input Variables 
    ! ---------------
    INTEGER, INTENT(IN) :: ktau
    REAL(kind=kind_chem), INTENT(IN) :: dt
    REAL(kind=kind_chem), INTENT(IN) :: u_phy
    REAL(kind=kind_chem), INTENT(IN) :: v_phy
    REAL(kind=kind_chem), INTENT(IN) :: dz8w
    REAL(kind=kind_chem), INTENT(IN) :: u10
    REAL(kind=kind_chem), INTENT(IN) :: v10
    REAL(kind=kind_chem), INTENT(IN) :: delp
    REAL(kind=kind_chem), INTENT(IN) :: tsk ! this isn't used anywhere
    REAL(kind=kind_chem), INTENT(IN) :: area
    
    ! Output Variables 
    ! ----------------
    REAL(kind=kind_chem), DIMENSION( num_chem ), INTENT(INOUT ) :: chem_arr
    REAL(kind=kind_chem), DIMENSION( num_emis_seas), INTENT(OUT) :: emis_seas
    
    ! local variables
    ! ---------------
    integer :: ipr,n,rc,ilwi
    real(kind=kind_chem) :: fsstemis
    real(kind=kind_chem) :: memissions
    real(kind=kind_chem) :: nemissions
    real(kind=kind_chem) :: ws10m
    real(kind=kind_chem) :: dummylon
    real(kind=kind_chem) :: fgridefficiency
    real(kind=kind_chem) :: deep_lakes_mask
    real(kind=kind_chem) :: w10m 
    real(kind=kind_chem) :: airmas
    real(kind=kind_chem), DIMENSION (number_ss_bins) :: tc
    real(kind=kind_chem), DIMENSION (number_ss_bins) :: bems
    
    real(kind=kind_chem) :: airmas1

    ! local parameters
    ! ----------------
    real(kind=kind_chem), parameter :: conver  = 1.e-9_kind_chem
    real(kind=kind_chem), parameter :: converi = 1.e+9_kind_chem


    ! -- original GOCART sea salt scheme
    ! -- compute auxiliary variables

    ! Calcualte 10m from first model layer if first layer is low
    if (dz8w < 12.) then
      ws10m = sqrt(u_phy*u_phy+v_phy*v_phy)
    else
      ws10m = sqrt(u10*u10+v10*v10)
    end if

    ilwi=0 ! this is always zero is it needed? 
    tc = 0.
        
    ! calcualte airmass 
    ! -----------------
    airmas=area * delp / g
    ipr=0

    call source_ss( number_ss_bins, dt, tc,ilwi, area, ws10m, airmas, bems,ipr)
    
    ! -- add sea salt emission increments to existing airborne concentrations
    do n = 0, number_ss_bins-1
       chem_arr(p_seas_1 + n)=chem_arr(p_seas_1 +n)+tc(n+1)*converi
    
       ! for output diagnostics
       emis_seas(p_eseas1 + n) = bems(n+1)
    end do
       
  end subroutine gocart_seas_default

  SUBROUTINE source_ss(nmx, dt1, tc, &
                       ilwi, dxy, w10m, airmas, &
                       bems,ipr)

! ****************************************************************************
! *  Evaluate the source of each seasalt particles size classes  (kg/m3) 
! *  by soil emission.
! *  Input:
! *         SSALTDEN  Sea salt density                               (kg/m3)
! *         DXY       Surface of each grid cell                     (m2)
! *         NDT1      Time step                                     (s)
! *         W10m      Velocity at the anemometer level (10meters)   (m/s)
! *      
! *  Output:
! *         DSRC      Source of each sea salt bins       (kg/timestep/cell) 
! *
! *
! * Number flux density: Original formula by Monahan et al. (1986) adapted
! * by Sunling Gong (JGR 1997 (old) and GBC 2003 (new)).  The new version is
! * to better represent emission of sub-micron sea salt particles.
!
! * dFn/dr = c1*u10**c2/(r**A) * (1+c3*r**c4)*10**(c5*exp(-B**2))
! * where B = (b1 -log(r))/b2
! * see c_old, c_new, b_old, b_new below for the constants.
! * number fluxes are at 80% RH.
! *
! * To calculate the flux:
! * 1) Calculate dFn based on Monahan et al. (1986) and Gong (2003)
! * 2) Assume that wet radius r at 80% RH = dry radius r_d *frh
! * 3) Convert particles flux to mass flux :
! *    dFM/dr_d = 4/3*pi*rho_d*r_d^3 *(dr/dr_d) * dFn/dr
! *             = 4/3*pi*rho_d*r_d^3 * frh * dFn/dr
! *               where rho_p is particle density [kg/m3]
! *    The factor 1.e-18 is to convert in micro-meter r_d^3
! ****************************************************************************
   

    IMPLICIT NONE

    INTEGER, INTENT(IN)    :: nmx,ipr
    INTEGER, INTENT(IN)    :: ilwi
    REAL(kind=kind_chem),    INTENT(IN)    :: dxy, w10m
    REAL(kind=kind_chem),    INTENT(IN)    :: airmas
    REAL(kind=kind_chem),    INTENT(INOUT) :: tc(nmx)
    REAL(kind=kind_chem),    INTENT(OUT)   :: bems(nmx)

    REAL(kind=kind_chem) :: c0(5), b0(2)
!  REAL(kind=kind_chem), PARAMETER :: c_old(5)=(/1.373, 3.41, 0.057, 1.05, 1.190/) 
!  REAL(kind=kind_chem), PARAMETER :: c_new(5)=(/1.373, 3.41, 0.057, 3.45, 1.607/)
    ! Change suggested by MC
    REAL(kind=kind_chem), PARAMETER :: c_old(5)=(/1.373, 3.2, 0.057, 1.05, 1.190/) 
    REAL(kind=kind_chem), PARAMETER :: c_new(5)=(/1.373, 3.2, 0.057, 3.45, 1.607/)
    REAL(kind=kind_chem), PARAMETER :: b_old(2)=(/0.380, 0.650/)
    REAL(kind=kind_chem), PARAMETER :: b_new(2)=(/0.433, 0.433/)
    REAL(kind=kind_chem), PARAMETER :: dr=5.0D-2 ! um   
    REAL(kind=kind_chem), PARAMETER :: theta=30.0
    ! Swelling coefficient frh (d rwet / d rd)
!!!  REAL(kind=kind_chem),    PARAMETER :: frh = 1.65
    REAL(kind=kind_chem),    PARAMETER :: frh = 2.d0
    LOGICAL, PARAMETER :: old=.TRUE., new=.FALSE.
    REAL(kind=kind_chem) :: rho_d, r0, r1, r, r_w, a, b, dfn, r_d, dfm, src
    INTEGER :: i, j, n, nr, ir
    REAL(kind=kind_chem) :: dt1,fudge_fac


    REAL(kind=kind_chem)    :: tcmw(nmx), ar(nmx), tcvv(nmx)
    REAL(kind=kind_chem)    :: ar_wetdep(nmx), kc(nmx)
    CHARACTER(LEN=20)     :: tcname(nmx), tcunits(nmx)
    LOGICAL               :: aerosol(nmx)


    REAL(kind=kind_chem) :: tc1(nmx)
    REAL(kind=kind_chem), TARGET :: tcms(nmx) ! tracer mass (kg; kgS for sulfur case)
    REAL(kind=kind_chem), TARGET :: tcgm(nmx) ! g/m3

    !-----------------------------------------------------------------------  
    ! emissions (input)
    !-----------------------------------------------------------------------  
    REAL(kind=kind_chem) :: e_an(2,nmx), e_bb(nmx), e_ac(nmx)

    !-----------------------------------------------------------------------  
    ! diagnostics (budget)
    !-----------------------------------------------------------------------
    !  ! tendencies per time step and process
    !  REAL(kind=kind_chem), TARGET :: bems(nmx), bdry(nmx), bstl(nmx)
    !  REAL(kind=kind_chem), TARGET :: bwet(nmx), bcnv(nmx)!
    
    !  ! integrated tendencies per process
    !  REAL(kind=kind_chem), TARGET :: tems(nmx), tstl(nmx)
    !  REAL(kind=kind_chem), TARGET :: tdry(nmx), twet(nmx), tcnv(nmx)
    
    ! global mass balance per time step 
    REAL(kind=kind_chem) :: tmas0(nmx), tmas1(nmx)
    REAL(kind=kind_chem) :: dtems(nmx), dttrp(nmx), dtdif(nmx), dtcnv(nmx)
    REAL(kind=kind_chem) :: dtwet(nmx), dtdry(nmx), dtstl(nmx)
    REAL(kind=kind_chem) :: dtems2(nmx), dttrp2(nmx), dtdif2(nmx), dtcnv2(nmx)
    REAL(kind=kind_chem) :: dtwet2(nmx), dtdry2(nmx), dtstl2(nmx)

    ! detailed integrated budgets for individual emissions
    REAL(kind=kind_chem), TARGET :: ems_an(nmx),    ems_bb(nmx), ems_tp
    REAL(kind=kind_chem), TARGET :: ems_ac(nmx)
    REAL(kind=kind_chem), TARGET :: ems_co(nmx)

    ! executable statements
    ! decrease seasalt emissions (Colarco et al. 2010)
    fudge_fac= .25 !lzhang

    DO n = 1,nmx
       bems(n) = 0.0
       rho_d = den_seas(n)
       r0 = ra(n)*frh
       r1 = rb(n)*frh
       r = r0
       nr = INT((r1-r0)/dr+.001)
       DO ir = 1,nr
          r_w = r + dr*0.5
          r = r + dr
          IF (new) THEN
             a = 4.7*(1.0 + theta*r_w)**(-0.017*r_w**(-1.44))
             c0 = c_new
             b0 = b_new
          ELSE
             a = 3.0
             c0 = c_old
             b0 = b_old
          END IF
          !
          b = (b0(1) - LOG10(r_w))/b0(2)
          dfn = (c0(1)/r_w**a)*(1.0 + c0(3)*r_w**c0(4))* &
               10**(c0(5)*EXP(-(b**2)))
          
          r_d = r_w/frh*1.0D-6  ! um -> m
          dfm = 4.0/3.0*pi*r_d**3*rho_d*frh*dfn*dr*dt1 ! 3600 !dt1

          IF (ilwi == 0) THEN
             src = dfm*dxy*w10m**c0(2)
             ! src = ch_ss(n,dt(1)%mn)*dfm*dxy(j)*w10m(i,j)**c0(2)
             tc(n) = tc(n) + fudge_fac*src/airmas
          ELSE
             src = 0.0
          END IF
          
          bems(n) = bems(n) + src*fudge_fac/(dxy*dt1) !kg/m2/s

       END DO ! ir
    END DO ! n

  END SUBROUTINE source_ss

end module gocart_seas_default_mod
