module vash_settling_mod

  use catchem_constants ,         only : kind_chem
  use catchem_config

  implicit none
 
  private
  public :: vash_settling_driver, &
            vashshort_settling_driver

CONTAINS

SUBROUTINE vash_settling_driver(dt,t_phy,moist,                            &
         chem,rho_phy,dz8w,p8w,p_phy,area,                                 &
         ash_fall,g,num_moist,num_chem,                                 &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )

  IMPLICIT NONE
   INTEGER,      INTENT(IN   ) ::                                          &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte,               &
                                  num_chem,num_moist
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                                   moist
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
   REAL(kind_chem),  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::  t_phy,p_phy,dz8w,p8w,rho_phy
   REAL(kind_chem),  DIMENSION( ims:ime , jms:jme ),                        &
          INTENT(INOUT   ) ::  ash_fall
   REAL(kind_chem),  DIMENSION( ims:ime , jms:jme ),                        &
          INTENT(IN   ) ::  area

  REAL(kind_chem), INTENT(IN   ) :: dt,g
  integer :: nmx,i,j,k,kk,lmx,iseas,idust
  real(kind_chem), DIMENSION (1,1,kte-kts+1) :: tmp,airden,airmas,p_mid,delz,rh
  real(kind_chem), DIMENSION (1,1,kte-kts+1,5) :: dust
  real(kind_chem), DIMENSION (1,1,kte-kts+1,4) :: sea_salt
!srf
  real(kind_chem), DIMENSION (1,1,kte-kts+1,10) :: ash
  real(kind_chem), DIMENSION (10), PARAMETER :: den_ash(10)=(/2500.,2500.,2500.,2500.,2500., &
                                                     2500.,2500.,2500.,2500.,2500. /)
  real(kind_chem), DIMENSION (10), PARAMETER :: reff_ash(10)=(/0.5000D-3,&! 1.00 mm diameter 
                                                      0.3750D-3,&! 0.75 mm
                              0.1875D-3,&!
                              93.750D-6,&!
                              46.875D-6,&!
                              23.437D-6,&!
                              11.719D-6,&!
                              05.859D-6,&!
                              02.930D-6,&!
                              00.975D-6 /)! 3.9 um
  real(kind_chem), DIMENSION (10) :: bstl_ash
  real(kind_chem) :: maxash(10)
  real(kind_chem) :: are
  integer nv,iprt,iash
!srf

!
! bstl is for budgets
!

  real(kind_chem) conver,converi
       conver=1.e-9
       converi=1.e9
       lmx=kte-kts+1
       do j=jts,jte
       do i=its,ite
          kk=0
          are=area(i,j)
      bstl_ash(:)=0.
          do k=kts,kte 
          kk=kk+1
          p_mid(1,1,kk)=.01*p_phy(i,kte-k+kts,j) 
          delz(1,1,kk)=dz8w(i,kte-k+kts,j) 
          airmas(1,1,kk)=-(p8w(i,k+1,j)-p8w(i,k,j))/g
          airden(1,1,kk)=rho_phy(i,k,j)
          tmp(1,1,kk)=t_phy(i,k,j)
          rh(1,1,kk) = .95
          rh(1,1,kk) = MIN( .95, moist(i,k,j,p_qv) / &
               (3.80*exp(17.27*(t_phy(i,k,j)-273.)/ &
               (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
          rh(1,1,kk)=max(1.0D-1,rh(1,1,kk))
          enddo

!ash settling
          maxash(:)=0.
          kk=0
          do nv=p_vash_1,p_vash_10
          kk=kk+1
          do k=kts,kte
              if(chem(i,k,j,nv).gt.maxash(kk)) maxash(kk)=chem(i,k,j,nv)
            enddo
          enddo

          iseas=0
          idust=0
      iash =1

          kk=0
          do k=kts,kte
          kk=kk+1
          if(chem(i,k,j,p_vash_1).le.1.e-10)chem(i,k,j,p_vash_1)=0.
          if(chem(i,k,j,p_vash_2).le.1.e-10)chem(i,k,j,p_vash_2)=0.
          if(chem(i,k,j,p_vash_3).le.1.e-10)chem(i,k,j,p_vash_3)=0.
          if(chem(i,k,j,p_vash_4).le.1.e-10)chem(i,k,j,p_vash_4)=0.
          if(chem(i,k,j,p_vash_5).le.1.e-10)chem(i,k,j,p_vash_5)=0.
          if(chem(i,k,j,p_vash_6).le.1.e-10)chem(i,k,j,p_vash_6)=0.
          if(chem(i,k,j,p_vash_7).le.1.e-10)chem(i,k,j,p_vash_7)=0.
          if(chem(i,k,j,p_vash_8).le.1.e-10)chem(i,k,j,p_vash_8)=0.
          if(chem(i,k,j,p_vash_9).le.1.e-10)chem(i,k,j,p_vash_9)=0.
          if(chem(i,k,j,p_vash_10).le.1.e-10)chem(i,k,j,p_vash_10)=0.
          ash(1,1,kk,1)=chem(i,k,j,p_vash_1)*conver
          ash(1,1,kk,2)=chem(i,k,j,p_vash_2)*conver
          ash(1,1,kk,3)=chem(i,k,j,p_vash_3)*conver
          ash(1,1,kk,4)=chem(i,k,j,p_vash_4)*conver
          ash(1,1,kk,5)=chem(i,k,j,p_vash_5)*conver
          ash(1,1,kk,6)=chem(i,k,j,p_vash_6)*conver
          ash(1,1,kk,7)=chem(i,k,j,p_vash_7)*conver
          ash(1,1,kk,8)=chem(i,k,j,p_vash_8)*conver
          ash(1,1,kk,9)=chem(i,k,j,p_vash_9)*conver
          ash(1,1,kk,10)=chem(i,k,j,p_vash_10)*conver
          enddo
          iprt=0
          call vsettling(iprt,1, 1, lmx, 10, g,are,&
                    ash, tmp, p_mid, delz, airmas, &
                    den_ash, reff_ash, dt, bstl_ash, rh, idust, iseas,iash)
          kk=0
          ash_fall(i,j)=ash_fall(i,j)+sum(bstl_ash(1:10))
          do k=kts,kte
          kk=kk+1
            chem(i,k,j,p_vash_1)=min(maxash(1),ash(1,1,kk,1)*converi)
            chem(i,k,j,p_vash_2)=min(maxash(2),ash(1,1,kk,2)*converi)
            chem(i,k,j,p_vash_3)=min(maxash(3),ash(1,1,kk,3)*converi)
            chem(i,k,j,p_vash_4)=min(maxash(4),ash(1,1,kk,4)*converi)
            chem(i,k,j,p_vash_5)=min(maxash(5),ash(1,1,kk,5)*converi)
            chem(i,k,j,p_vash_6)=min(maxash(6),ash(1,1,kk,6)*converi)
            chem(i,k,j,p_vash_7)=min(maxash(7),ash(1,1,kk,7)*converi)
            chem(i,k,j,p_vash_8)=min(maxash(8),ash(1,1,kk,8)*converi)
            chem(i,k,j,p_vash_9)=min(maxash(9),ash(1,1,kk,9)*converi)
            chem(i,k,j,p_vash_10)=min(maxash(10),ash(1,1,kk,10)*converi)
          if(chem(i,k,j,p_vash_1).le.1.e-10)chem(i,k,j,p_vash_1)=0.
          if(chem(i,k,j,p_vash_2).le.1.e-10)chem(i,k,j,p_vash_2)=0.
          if(chem(i,k,j,p_vash_3).le.1.e-10)chem(i,k,j,p_vash_3)=0.
          if(chem(i,k,j,p_vash_4).le.1.e-10)chem(i,k,j,p_vash_4)=0.
          if(chem(i,k,j,p_vash_5).le.1.e-10)chem(i,k,j,p_vash_5)=0.
          if(chem(i,k,j,p_vash_6).le.1.e-10)chem(i,k,j,p_vash_6)=0.
          if(chem(i,k,j,p_vash_7).le.1.e-10)chem(i,k,j,p_vash_7)=0.
          if(chem(i,k,j,p_vash_8).le.1.e-10)chem(i,k,j,p_vash_8)=0.
          if(chem(i,k,j,p_vash_9).le.1.e-10)chem(i,k,j,p_vash_9)=0.
          if(chem(i,k,j,p_vash_10).le.1.e-10)chem(i,k,j,p_vash_10)=0.
          enddo

!ash settling end




       enddo
       enddo
END SUBROUTINE vash_settling_driver

SUBROUTINE vashshort_settling_driver(dt,t_phy,moist,                       &
         chem,rho_phy,dz8w,p8w,p_phy,area,                                 &
         ash_fall,g,num_moist,num_chem,                                    &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )

  IMPLICIT NONE
   INTEGER,      INTENT(IN   ) ::                                          &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte,               &
                                  num_chem,num_moist
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),               &
         INTENT(IN ) ::                                   moist
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
   REAL(kind_chem),  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::  t_phy,p_phy,dz8w,p8w,rho_phy
   REAL(kind_chem),  DIMENSION( ims:ime , jms:jme ),                        &
          INTENT(INOUT   ) ::  ash_fall
   REAL(kind_chem),  DIMENSION( ims:ime , jms:jme ),                        &
          INTENT(IN   ) ::  area

  REAL(kind_chem), INTENT(IN   ) :: dt,g
  integer :: nmx,i,j,k,kk,lmx,iseas,idust
  real(kind_chem), DIMENSION (1,1,kte-kts+1) :: tmp,airden,airmas,p_mid,delz,rh
  real(kind_chem), DIMENSION (1,1,kte-kts+1,5) :: dust
  real(kind_chem), DIMENSION (1,1,kte-kts+1,4) :: sea_salt
!srf
  real(kind_chem), DIMENSION (1,1,kte-kts+1,10) :: ash 
  real(kind_chem), DIMENSION (4), PARAMETER :: den_ash(4)=(/2500.,2500.,2500.,2500. /)
  real(kind_chem), DIMENSION (4), PARAMETER :: reff_ash(4)=(/ 11.719D-6,&!
                              05.859D-6,&!
                              02.930D-6,&!
                              00.975D-6 /)! 3.9 um
  real(kind_chem), DIMENSION (4) :: bstl_ash
  real(kind_chem) :: maxash(4)
  real(kind_chem) :: are
  integer nv,iprt,iash
!srf

!
! bstl is for budgets
!

  real(kind_chem) conver,converi

       conver=1.e-9
       converi=1.e9
       lmx=kte-kts+1
       do j=jts,jte
       do i=its,ite
          kk=0
          are=area(i,j)
      bstl_ash(:)=0.
          do k=kts,kte 
          kk=kk+1
          p_mid(1,1,kk)=.01*p_phy(i,kte-k+kts,j) 
          delz(1,1,kk)=dz8w(i,kte-k+kts,j) 
          airmas(1,1,kk)=-(p8w(i,k+1,j)-p8w(i,k,j))/g
          airden(1,1,kk)=rho_phy(i,k,j)
          tmp(1,1,kk)=t_phy(i,k,j)
          rh(1,1,kk) = .95
          rh(1,1,kk) = MIN( .95, moist(i,k,j,p_qv) / &
               (3.80*exp(17.27*(t_phy(i,k,j)-273.)/ &
               (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
          rh(1,1,kk)=max(1.0D-1,rh(1,1,kk))
          enddo
     

!ash settling
          kk=0
          maxash(:)=0.
          if(p_vash_4.gt.1)then
          do nv=p_vash_1,p_vash_4
          kk=kk+1
          do k=kts,kte
              if(chem(i,k,j,nv).gt.maxash(kk)) maxash(kk)=chem(i,k,j,nv)
            enddo
          enddo
! GOCART
          else if(p_bc2.gt.1)then
          nv=p_p25
          kk=kk+1
          do k=kts,kte
              if(chem(i,k,j,nv).gt.maxash(kk)) maxash(kk)=chem(i,k,j,nv)
          enddo
          nv=p_p10
          kk=kk+1
          do k=kts,kte
              if(chem(i,k,j,nv).gt.maxash(kk)) maxash(kk)=chem(i,k,j,nv)
          enddo
          endif


          iseas=0
          idust=0
      iash =1

          if(p_vash_4.gt.1)then
          kk=0
          do k=kts,kte 
          kk=kk+1
          if(chem(i,k,j,p_vash_1).le.1.e-10)chem(i,k,j,p_vash_1)=0.
          if(chem(i,k,j,p_vash_2).le.1.e-10)chem(i,k,j,p_vash_2)=0.
          if(chem(i,k,j,p_vash_3).le.1.e-10)chem(i,k,j,p_vash_3)=0.
          if(chem(i,k,j,p_vash_4).le.1.e-10)chem(i,k,j,p_vash_4)=0.
          ash(1,1,kk,1)=chem(i,k,j,p_vash_1)*conver
          ash(1,1,kk,2)=chem(i,k,j,p_vash_2)*conver
          ash(1,1,kk,3)=chem(i,k,j,p_vash_3)*conver
          ash(1,1,kk,4)=chem(i,k,j,p_vash_4)*conver
          enddo
!
! volc ash for gocart, this is crude
!
          else if(p_bc2.gt.1)then
             kk=0
             do k=kts,kte 
                kk=kk+1
                ash(1,1,kk,1)=0.
                ash(1,1,kk,4)=chem(i,k,j,p_p25)*conver
                ash(1,1,kk,3)=.67*chem(i,k,j,p_p10)*conver
                ash(1,1,kk,2)=(1.-.67)*chem(i,k,j,p_p10)*conver
                if(ash(1,1,kk,2).le.1.e-10)ash(1,1,kk,2)=0.
                if(ash(1,1,kk,3).le.1.e-10)ash(1,1,kk,3)=0.
                if(ash(1,1,kk,4).le.1.e-10)ash(1,1,kk,4)=0.
             enddo
          endif
             iprt=0
          call vsettling(iprt,1, 1, lmx, 4, g, are,&
                    ash, tmp, p_mid, delz, airmas, &
                    den_ash, reff_ash, dt, bstl_ash, rh, idust, iseas,iash)
          ash_fall(i,j)=ash_fall(i,j)+sum(bstl_ash(1:4))
          if(p_vash_4.gt.1)then
          kk=0
          do k=kts,kte
          kk=kk+1
            chem(i,k,j,p_vash_1)=min(maxash(1),ash(1,1,kk,1)*converi)
            chem(i,k,j,p_vash_2)=min(maxash(2),ash(1,1,kk,2)*converi)
            chem(i,k,j,p_vash_3)=min(maxash(3),ash(1,1,kk,3)*converi)
            chem(i,k,j,p_vash_4)=min(maxash(4),ash(1,1,kk,4)*converi)
          if(chem(i,k,j,p_vash_1).le.1.e-10)chem(i,k,j,p_vash_1)=0.
          if(chem(i,k,j,p_vash_2).le.1.e-10)chem(i,k,j,p_vash_2)=0.
          if(chem(i,k,j,p_vash_3).le.1.e-10)chem(i,k,j,p_vash_3)=0.
          if(chem(i,k,j,p_vash_4).le.1.e-10)chem(i,k,j,p_vash_4)=0.
          enddo
          else if(p_bc2.gt.1)then
          kk=0
          do k=kts,kte
          kk=kk+1
!           chem(i,k,j,p_p25)=min(maxash(1),ash(1,1,kk,4)*converi)
            chem(i,k,j,p_p10)=min(maxash(2),(ash(1,1,kk,2)+ash(1,1,kk,3))*converi)
!         if(chem(i,k,j,p_p25).le.1.e-16)chem(i,k,j,p_p25)=1.e-16
          if(chem(i,k,j,p_p10).le.1.e-16)chem(i,k,j,p_p10)=1.e-16
          enddo
          endif

!ash settling end

       enddo
       enddo

END SUBROUTINE vashshort_settling_driver


          subroutine vsettling(iprt,imx,jmx, lmx, nmx,g0,are, &
                    tc, tmp, p_mid, delz, airmas, &
                    den, reff, dt, bstl, rh, idust, iseas,iash)
! ****************************************************************************
! *                                                                          *
! *  Calculate the loss by settling, using an implicit method                *
! *                                                                          *
! *  Input variables:                                                        *
! *    SIGE(k)         - sigma coordinate of the vertical edges              *
! *    PS(i,j)         - Surface pressure (mb)                               *
! *    TMP(i,j,k)      - Air temperature  (K)                                *
! *    CT(i,j)         - Surface exchange coeff for moisture
! *                                                                          *
! **************************************************************************** 


  IMPLICIT  NONE

  INTEGER, INTENT(IN) :: imx, jmx, lmx, nmx,iseas,idust,iash
  INTEGER :: ntdt
  REAL(kind_chem), INTENT(IN) :: dt,g0,are ! ,dyn_visc
  REAL(kind_chem),    INTENT(IN) :: tmp(imx,jmx,lmx), delz(imx,jmx,lmx),  &
                         airmas(imx,jmx,lmx), rh(imx,jmx,lmx), &
                         den(nmx), reff(nmx), p_mid(imx,jmx,lmx)
  REAL(kind_chem), INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL(kind_chem), INTENT(OUT)   :: bstl(imx,jmx,nmx)

  REAL(kind_chem)    :: tc1(imx,jmx,lmx,nmx), dt_settl(nmx), rcm(nmx), rho(nmx),addmass(lmx,nmx)
  INTEGER :: ndt_settl(nmx)
  REAL(kind_chem)    :: dzmin, vsettl, dtmax, pres, rhb, rwet(nmx), ratio_r(nmx)
  REAL(kind_chem)    :: addmassf,c_stokes, free_path, c_cun, viscosity, vd_cor, growth_fac,mass_above
  REAL(kind_chem),    PARAMETER :: dyn_visc = 1.5E-5
  INTEGER :: iprt,k, n, i, j, l, l2
  ! for sea-salt:
  REAL(kind_chem), PARAMETER :: c1=0.7674, c2=3.079, c3=2.573E-11, c4=-1.424 

  ! for OMP:
  REAL(kind_chem) :: rwet_priv(nmx), rho_priv(nmx),vsettl_max(nmx)

  ! executable statements

! IF (type) /= 'dust' .AND. TRIM(aero_type) /= 'sea_salt') RETURN
  if(idust.ne.1.and.iseas.ne.1.and.iash.ne.1)return

  WHERE (tc(:,:,:,:) < 0.0) tc(:,:,:,:) = 1.0d-32

  dzmin = MINVAL(delz(:,:,:))
  IF (idust == 1)     growth_fac = 1.0
  IF (iseas == 1)     growth_fac = 3.0
  IF (iash  == 1)     growth_fac = 1.0

  DO k = 1,nmx

     ! Settling velocity (m/s) for each tracer (Stokes Law)
     ! DEN         density                        (kg/m3)
     ! REFF        effective radius               (m)
     ! dyn_visc    dynamic viscosity              (kg/m/s)
     ! g0          gravity                        (m/s2)
     ! 3.0         corresponds to a growth of a factor 3 of radius with 100% RH
     ! 0.5         upper limit with temp correction

     tc1(:,:,:,k) = tc(:,:,:,k)
     vsettl = 2.0/9.0 * g0 * den(k) * (growth_fac*reff(k))**2 / &
              (0.5*dyn_visc)
     vsettl_max(k)=vsettl
     ! Determine the maximum time-step satisying the CFL condition:
     ! dt <= (dz)_min / v_settl
     ntdt=INT(dt)
     dtmax = dzmin / vsettl
     ndt_settl(k) = MAX( 1, INT( ntdt /dtmax) )
     ! limit maximum number of iterations
!    IF (ndt_settl(k) > 12) ndt_settl(k) = 12
     IF (ndt_settl(k) > 12) then
         ndt_settl(k) = 12
         vsettl_max(k)=dzmin*ndt_settl(k)/dt
     endif
     dt_settl(k) = REAL(ntdt) / REAL(ndt_settl(k))

     ! Particles radius in centimeters
     IF (iseas.eq.1)rcm(k) = reff(k)*100.0
!srf     IF (idust.eq.1)then
     IF (idust.eq.1   .or. iash==1)then
          rwet(k) = reff(k)
          ratio_r(k) = 1.0
          rho(k) = den(k)
      endif
  END DO

  ! Solve the bidiagonal matrix (l,l)

!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( i,   j,   l,   l2, n,   k,   rhb, rwet_priv, ratio_r, c_stokes)&
!$OMP PRIVATE( free_path, c_cun, viscosity, rho_priv, vd_cor )

  ! Loop over latitudes
  DO j = 1,jmx
 
     DO k = 1,nmx
        IF (idust.eq.1 .or. iash==1) THEN
           rwet_priv(k) = rwet(k)
           rho_priv(k)  = rho(k)
        END IF

        DO n = 1,ndt_settl(k)

           ! Solve each vertical layer successively (layer l)
      
           DO l = lmx,1,-1
              l2 = lmx - l + 1

!           DO j = 1,jmx
              DO i = 1,imx

                 ! Dynamic viscosity
                 c_stokes = 1.458E-6 * tmp(i,j,l)**1.5/(tmp(i,j,l) + 110.4) 

                 ! Mean free path as a function of pressure (mb) and 
                 ! temperature (K)
                 ! order of p_mid is top->sfc
                 free_path = 1.1E-3/p_mid(i,j,l2)/SQRT(tmp(i,j,l))
!!!                 free_path = 1.1E-3/p_edge(i,j,l2)/SQRT(tmp(i,j,l))

                 ! Slip Correction Factor
                 c_cun = 1.0+ free_path/rwet_priv(k)* &
                      (1.257 + 0.4*EXP(-1.1*rwet_priv(k)/free_path))

                 ! Corrected dynamic viscosity (kg/m/s)
                 viscosity = c_stokes / c_cun

                 ! Settling velocity
!                IF (iseas.eq.1) THEN
!                   rho_priv(k) = ratio_r(k)*den(k) + (1.0 - ratio_r(k))*1000.0
!                END IF

                 vd_cor = min(vsettl_max(k),2.0/9.0*g0*rho_priv(k)*rwet_priv(k)**2/viscosity)

                 ! Update mixing ratio
                 ! Order of delz is top->sfc
                 IF (l == lmx) THEN
                    tc(i,j,l,k) = tc(i,j,l,k) / &
                         (1.0 + dt_settl(k)*vd_cor/delz(i,j,l2))
                 ELSE
                    tc(i,j,l,k) = 1.0/(1.0+dt_settl(k)*vd_cor/delz(i,j,l2))&
                         *(tc(i,j,l,k) + dt_settl(k)*vd_cor /delz(i,j,l2-1) &
                         * tc(i,j,l+1,k))
                 END IF
              END DO   !i
!           END DO   !j
        END DO  !l

     END DO  !n
  END DO  !k

  END DO   !j
!$OMP END PARALLEL DO

  DO n = 1,nmx
     DO i = 1,imx
        DO j = 1,jmx
           DO l = 1,lmx
              IF (tc(i,j,l,n) < 0.0) tc(i,j,l,n) = 1.0D-32
              addmass(l,n)=(tc(i,j,l,n) - tc1(i,j,l,n)) * airmas(i,j,l)
           END DO
! make sure this is not more mass then what  there was in the layer above
           DO l = lmx-1,1
              mass_above=tc1(i,j,l+1,n)*airmas(i,j,l+1)
              IF (addmass(l,n).gt.mass_above)then
               tc(i,j,l,n)=mass_above/airmas(i,j,l) + tc1(i,j,l,n)
               IF (tc(i,j,l,n) < 0.0) tc(i,j,l,n) = 1.0D-32
               addmass(l,n)=mass_above
              endif
           END DO
        END DO
     END DO
  END DO
  DO n = 1,nmx
     DO i = 1,imx
        DO j = 1,jmx
           bstl(i,j,n) = 0.0
           addmassf=0.
           DO l = 1,lmx
              addmassf=addmassf+(tc(i,j,l,n) - tc1(i,j,l,n)) * airmas(i,j,l)
!             IF (tc(i,j,l,n) < 0.0) tc(i,j,l,n) = 1.0D-32
           END DO
           if(addmassf.gt.0.)addmassf=0
           bstl(i,j,n) = bstl(i,j,n) - addmassf
        END DO
     END DO
  END DO
  
END SUBROUTINE vsettling

end module vash_settling_mod
