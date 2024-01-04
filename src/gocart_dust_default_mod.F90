!! Revision History:
!! 06/2023, Restructure for CATChem, Jian.He@noaa.gov

module gocart_dust_default_mod

   use catchem_constants, only : kind_chem, g=>con_g, pi=>con_pi
   use catchem_config, only: num_chem,num_emis_dust,&
      p_dust_1,p_dust_2,p_dust_3,p_dust_4,p_dust_5, &
      p_edust1,p_edust2,p_edust3,p_edust4,p_edust5
   use dust_data_mod

   implicit none

   private

   public :: gocart_dust_default

CONTAINS

   subroutine gocart_dust_default(ktau,dt,u_phy,              &
      v_phy,chem_arr,rho_phy,dz8w,smois,u10,         &
      v10,delp,erod,isltyp,area,       &
      emis_dust,srce_dust,num_soil_layers,start_month)

      INTEGER,      INTENT(IN   ) :: ktau, isltyp, start_month, &
         num_soil_layers
      REAL(kind=kind_chem), INTENT(IN   ) :: dt,u_phy,v_phy, &
         rho_phy, &
         dz8w,u10,v10,  &
         delp,           &
         area

      REAL(kind=kind_chem), DIMENSION( num_chem ),                 &
         INTENT(INOUT ) ::                                   chem_arr
      REAL(kind=kind_chem), DIMENSION(num_emis_dust),&
         INTENT(INOUT ) ::                  emis_dust, srce_dust
      REAL(kind=kind_chem), DIMENSION( num_soil_layers ) , &
         INTENT(INOUT) ::                               smois
      REAL(kind=kind_chem),  DIMENSION( 3 ) ,               &
         INTENT(IN   ) ::    erod

!
!
! local variables
!
      integer :: ipr,ilwi
      real(kind=kind_chem), DIMENSION (3,1) :: erodin  ! (ndcls,ndsrc)
      real(kind=kind_chem), DIMENSION (ndust) :: tc,bems,srce_out
      real(kind=kind_chem) :: w10m,gwet,airden,airmas
      real(kind=kind_chem) :: dxy
      real(kind=kind_chem)  tcs
      real(kind=kind_chem)  dttt
      real(kind=kind_chem), parameter::max_default=0.
      real(kind=kind_chem), parameter :: conver  = 1.e-9
      real(kind=kind_chem), parameter :: converi = 1.e+9

      ilwi=1
      tc(1)=chem_arr(p_dust_1)*conver
      tc(2)=chem_arr(p_dust_2)*conver
      tc(3)=chem_arr(p_dust_3)*conver
      tc(4)=chem_arr(p_dust_4)*conver
      tc(5)=chem_arr(p_dust_5)*conver
      w10m=sqrt(u10*u10+v10*v10)
      airmas=area * delp / g
!
! we don't trust the u10,v10 values, is model layers are very thin near surface
!
      if(dz8w.lt.12.)w10m=sqrt(u_phy*u_phy+v_phy*v_phy)
!
      erodin(1,1)=erod(1)!/area(i,j)
      erodin(2,1)=erod(2)!/area(i,j)
      erodin(3,1)=erod(3)!/area(i,j)

! -- volumetric soil moisture over porosity
      if(isltyp.eq.0)then
         ilwi=0
         gwet=1.
      else
         gwet=smois(1)/maxsmc(isltyp)
      endif

      airden=rho_phy
      dxy=area
      ipr=0

      call source_du( ndust, dt, tc, &
         erodin, ilwi, dxy, w10m, gwet, airden, airmas, &
         bems,srce_out, start_month,g,ipr)

      chem_arr(p_dust_1)=max(max_default,tc(1)*converi)
      chem_arr(p_dust_2)=max(max_default,tc(2)*converi)
      chem_arr(p_dust_3)=max(max_default,tc(3)*converi)
      chem_arr(p_dust_4)=max(max_default,tc(4)*converi)
      chem_arr(p_dust_5)=max(max_default,tc(5)*converi)

      ! -- for output diagnostics
      emis_dust(p_edust1)=bems(1)
      emis_dust(p_edust2)=bems(2)
      emis_dust(p_edust3)=bems(3)
      emis_dust(p_edust4)=bems(4)
      emis_dust(p_edust5)=bems(5)

      ! -- for output diagnostics of dust source
      srce_dust(p_edust1)=srce_out(1)
      srce_dust(p_edust2)=srce_out(2)
      srce_dust(p_edust3)=srce_out(3)
      srce_dust(p_edust4)=srce_out(4)
      srce_dust(p_edust5)=srce_out(5)

   end subroutine gocart_dust_default


   SUBROUTINE source_du( nmx, dt1, tc, &
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

      REAL(kind=kind_chem),               INTENT(IN)    :: dt1, g0
      REAL(kind=kind_chem), INTENT(IN)    :: erod(ndcls,ndsrc)
      REAL(kind=kind_chem), INTENT(IN)    :: w10m, gwet
      REAL(kind=kind_chem), INTENT(IN)    :: dxy
      REAL(kind=kind_chem), INTENT(IN)    :: airden, airmas
      REAL(kind=kind_chem), INTENT(INOUT) :: tc(nmx)
      REAL(kind=kind_chem), INTENT(OUT)   :: bems(nmx)
      REAL(kind=kind_chem), INTENT(OUT)   :: srce_out(nmx) !dust source
      INTEGER,            INTENT(OUT)   :: ipr

      !-----------------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------------
      INTEGER            :: i, j, n, m, k
      REAL(kind=kind_chem) :: g
      REAL(kind=kind_chem) :: den(nmx), diam(nmx)
      REAL(kind=kind_chem) :: rhoa, tsrc, u_ts0, u_ts, dsrc, srce
      REAL(kind=kind_chem) :: tcmw(nmx), ar(nmx), tcvv(nmx)
      REAL(kind=kind_chem) :: ar_wetdep(nmx), kc(nmx)
      CHARACTER(LEN=20)  :: tcname(nmx), tcunits(nmx)
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
         DO k = 1, ndsrc
            ! No flux if wet soil
            rhoa = airden*1.0D-3
            u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
               SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
               SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0)

            ! Case of surface dry enough to erode
            IF (gwet < gthresh) THEN
               u_ts = MAX(0.0D+0,u_ts0*(1.2D+0+2.0D-1*LOG10(MAX(1.0D-3, gwet))))
            ELSE
               ! Case of wet surface, no erosion
               u_ts = 100.0_kind_chem
            END IF
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
            tc(n) = tc(n) + .7*dsrc / airmas
            bems(n) = .7*dsrc/(dxy*dt1) ! diagnostic (kg/m2/s)
         END DO
      END DO

   END SUBROUTINE source_du

end module gocart_dust_default_mod
