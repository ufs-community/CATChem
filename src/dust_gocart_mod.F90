module dust_gocart_mod

  use catchem_constants ,         only : kind_chem
  use dust_data_mod

  implicit none

  private

  public :: gocart_dust_driver

contains

  subroutine gocart_dust_driver(chem_opt,ktau,dt,alt,t_phy,moist,u_phy,  &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,                  &
         ivgtyp,isltyp,vegfra,xland,xlat,xlong,gsw,area,g,emis_dust,srce_dust,  & 
         dusthelp,num_emis_dust,num_moist,num_chem,num_soil_layers,           &
         start_month,                                               &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )

    IMPLICIT NONE

     INTEGER,      INTENT(IN   ) :: chem_opt,ktau,start_month,               &
           num_emis_dust,num_moist,num_chem,num_soil_layers,                 &
                               ids,ide, jds,jde, kds,kde,               &
                                    ims,ime, jms,jme, kms,kme,               &
                                    its,ite, jts,jte, kts,kte
     INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
            INTENT(IN   ) ::                                                 &
                                                       ivgtyp,               &
                                                       isltyp
     REAL(kind=kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
           INTENT(IN ) ::                                   moist
     REAL(kind=kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
           INTENT(INOUT ) ::                                   chem
     REAL(kind=kind_chem), DIMENSION( ims:ime, 1, jms:jme,num_emis_dust),                    &
           INTENT(INOUT ) ::                                                 &
           emis_dust
     REAL(kind=kind_chem), DIMENSION( ims:ime, 1, jms:jme,num_emis_dust),                    &
           INTENT(INOUT ) ::                                                 &
           srce_dust  
     REAL(kind=kind_chem), DIMENSION( ims:ime, num_soil_layers, jms:jme ) ,      &
        INTENT(INOUT) ::                               smois
     REAL(kind=kind_chem),  DIMENSION( ims:ime , jms:jme, 3 )                   ,               &
            INTENT(IN   ) ::    erod
     REAL(kind=kind_chem),  DIMENSION( ims:ime , jms:jme )                   ,               &
            INTENT(IN   ) ::                                                 &
                                                       u10,                  &
                                                       v10,                  &
                                                       gsw,                  &
                                                    vegfra,                  &
                                                       xland,                &
                                                       xlat,                 &
                                                       xlong,area
     REAL(kind=kind_chem),  DIMENSION( ims:ime , jms:jme ),                        &
            INTENT(OUT  ) :: dusthelp
     REAL(kind=kind_chem),  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
            INTENT(IN   ) ::                                                 &
                                                          alt,               &
                                                        t_phy,               &
                                                       dz8w,p8w,             &
                                                u_phy,v_phy,rho_phy
   
    REAL(kind=kind_chem), INTENT(IN   ) :: dt,g
!
! local variables
!
    integer :: nmx,i,j,k,imx,jmx,lmx,ipr
    integer,dimension (1,1) :: ilwi
    real(kind=kind_chem), DIMENSION (1,1,3,1) :: erodin
    real(kind=kind_chem), DIMENSION (5) :: tc,bems,srce_out 
    real(kind=kind_chem), dimension (1,1) :: w10m,gwet,airden,airmas
    real(kind=kind_chem), dimension (1) :: dxy
    real(kind=kind_chem)  tcs,conver,converi
    real(kind=kind_chem) dttt
    real(kind=kind_chem),parameter::max_default=0.
! write(6,*)'in dust driver ',ktau,dt,start_month
! conver=1.e-9*mwdry
! converi=1.e9/mwdry
    conver=1.e-9
    converi=1.e9
!
! number of dust bins
!
    imx=1
    jmx=1
    lmx=1
    nmx=5 
    k=kts

    select case (chem_opt)
    case (304, 316, 317)

      dusthelp(:,:)=0.
      do j=jts,jte
        do i=its,ite
         if(xland(i,j).lt.1.5 .and. xland(i,j).gt.0.5)then
           ilwi(1,1)=1
           tc(1)=chem(i,kts,j,p_dust_1)*conver
           tcs=tc(1)
           tc(2)=1.d-30
           tc(3)=chem(i,kts,j,p_dust_2)*conver
           tc(4)=1.d-30
           tc(5)=1.d-30
           w10m(1,1)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
           airmas(1,1)=-(p8w(i,kts+1,j)-p8w(i,kts,j))*area(i,j)/g

           ! -- we don't trust the u10,v10 values, is model layers are very thin near surface
           if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
           erodin(1,1,1,1)=erod(i,j,1)!/area(i,j)
           erodin(1,1,2,1)=erod(i,j,2)!/area(i,j)
           erodin(1,1,3,1)=erod(i,j,3)!/area(i,j)

           ! -- volumetric soil moisture over porosity
           if(isltyp(i,j).eq.0)then
             ilwi(1,1)=0
             gwet(1,1)=1.
           else
             gwet(1,1)=smois(i,1,j)/maxsmc(isltyp(i,j))
           endif

           airden(1,1)=rho_phy(i,kts,j)
           dxy(1)=area(i,j)
           ipr=0
           call source_du( imx,jmx,lmx,nmx, dt, tc, & 
                            erodin, ilwi, dxy, w10m, gwet, airden, airmas, &
                            bems,srce_out, start_month,g,ipr) 
           chem(i,kts,j,p_dust_1)=max(max_default,(tc(1)+.3125*tc(2))*converi)
           chem(i,kts,j,p_dust_2)=max(max_default,(.67*tc(2)+tc(3))*converi)
           dusthelp(i,j)=max(max_default,tc(2)*converi)

           ! -- for output diagnostics
           emis_dust(i,1,j,p_edust1)=bems(1)
           emis_dust(i,1,j,p_edust2)=bems(2)
           emis_dust(i,1,j,p_edust3)=bems(3)
         endif
       enddo
     enddo

    case default
   
    do j=jts,jte
      do i=its,ite
       if(xland(i,j).lt.1.5 .and. xland(i,j).gt.0.5)then
!         print*,'default dust'
         ilwi(1,1)=1
         tc(1)=chem(i,kts,j,p_dust_1)*conver
         tc(2)=chem(i,kts,j,p_dust_2)*conver
         tc(3)=chem(i,kts,j,p_dust_3)*conver
         tc(4)=chem(i,kts,j,p_dust_4)*conver
         tc(5)=chem(i,kts,j,p_dust_5)*conver
         w10m(1,1)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
         airmas(1,1)=-(p8w(i,kts+1,j)-p8w(i,kts,j))*area(i,j)/g

         ! -- we don't trust the u10,v10 values, is model layers are very thin near surface
         if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
           erodin(1,1,1,1)=erod(i,j,1)!/area(i,j)
           erodin(1,1,2,1)=erod(i,j,2)!/area(i,j)
           erodin(1,1,3,1)=erod(i,j,3)!/area(i,j)

           ! -- volumetric soil moisture over porosity
           if(isltyp(i,j).eq.0)then
             ilwi(1,1)=0
             gwet(1,1)=1.
           else
             gwet(1,1)=smois(i,1,j)/maxsmc(isltyp(i,j))
           endif

           airden(1,1)=rho_phy(i,kts,j)
           dxy(1)=area(i,j)
           ipr=0

         call source_du( imx,jmx,lmx,nmx, dt, tc, & 
                          erodin, ilwi, dxy, w10m, gwet, airden, airmas, &
                          bems,srce_out, start_month,g,ipr) 

         chem(i,kts,j,p_dust_1)=max(max_default,tc(1)*converi)
         chem(i,kts,j,p_dust_2)=max(max_default,tc(2)*converi)
         chem(i,kts,j,p_dust_3)=max(max_default,tc(3)*converi)
         chem(i,kts,j,p_dust_4)=max(max_default,tc(4)*converi)
         chem(i,kts,j,p_dust_5)=max(max_default,tc(5)*converi)

         ! -- for output diagnostics
         emis_dust(i,1,j,p_edust1)=bems(1)
         emis_dust(i,1,j,p_edust2)=bems(2)
         emis_dust(i,1,j,p_edust3)=bems(3)
         emis_dust(i,1,j,p_edust4)=bems(4)
         emis_dust(i,1,j,p_edust5)=bems(5)

         ! -- for output diagnostics of dust source 
         srce_dust(i,1,j,p_edust1)=srce_out(1)
         srce_dust(i,1,j,p_edust2)=srce_out(2)
         srce_dust(i,1,j,p_edust3)=srce_out(3)
         srce_dust(i,1,j,p_edust4)=srce_out(4)
         srce_dust(i,1,j,p_edust5)=srce_out(5)
       endif
      enddo
    enddo

    end select

  end subroutine gocart_dust_driver

  
  SUBROUTINE source_du( imx,jmx,lmx,nmx, dt1, tc, & 
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

  INTEGER,            INTENT(IN)    :: imx,jmx,lmx,nmx
  INTEGER,            INTENT(IN)    :: ilwi(imx,jmx),month

  REAL(kind=kind_chem),               INTENT(IN)    :: dt1, g0
  REAL(kind=kind_chem), INTENT(IN)    :: erod(imx,jmx,ndcls,ndsrc)
  REAL(kind=kind_chem), INTENT(IN)    :: w10m(imx,jmx), gwet(imx,jmx)
  REAL(kind=kind_chem), INTENT(IN)    :: dxy(jmx)
  REAL(kind=kind_chem), INTENT(IN)    :: airden(imx,jmx,lmx), airmas(imx,jmx,lmx)
  REAL(kind=kind_chem), INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL(kind=kind_chem), INTENT(OUT)   :: bems(imx,jmx,nmx)
  REAL(kind=kind_chem), INTENT(OUT)   :: srce_out(imx,jmx,nmx) !dust source
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
        DO i = 1,imx
           DO j = 1,jmx
              rhoa = airden(i,j,1)*1.0D-3
              u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
                   SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
                   SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0) 
              
              ! Case of surface dry enough to erode
              IF (gwet(i,j) < gthresh) THEN
                 u_ts = MAX(0.0D+0,u_ts0*(1.2D+0+2.0D-1*LOG10(MAX(1.0D-3, gwet(i,j)))))
              ELSE
                 ! Case of wet surface, no erosion
                 u_ts = 100.0_kind_chem
              END IF
              srce = frac_s(n)*erod(i,j,m,k)*dxy(j)  ! (m2)
              srce_out (i,j,k)=srce  ! output dust source
              IF (ilwi(i,j) == 1 ) THEN
                 dsrc = ch_dust(n,month)*srce*w10m(i,j)**2 &
                      * (w10m(i,j) - u_ts)*dt1  ! (kg)
              ELSE 
                 dsrc = 0.0_kind_chem
              END IF
              IF (dsrc < 0.0_kind_chem) dsrc = 0.0_kind_chem
              
              ! Update dust mixing ratio at first model level.
              ! scale down dust by .7
              tc(i,j,1,n) = tc(i,j,1,n) + .7*dsrc / airmas(i,j,1)
              bems(i,j,n) = .7*dsrc/(dxy(j)*dt1) ! diagnostic (kg/m2/s)
           END DO
        END DO
     END DO
  END DO

END SUBROUTINE source_du

end module dust_gocart_mod
