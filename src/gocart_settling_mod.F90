! 11/2023, Updated by Kate.Zhang@noaa.gov
module gocart_settling_mod

  use catchem_constants ,        only : kind_chem, grav => con_g
  use catchem_config, only : p_seas_1, p_seas_2, p_seas_3, p_seas_4, p_seas_5, &
                             p_dust_1, p_dust_2, p_dust_3, p_dust_4, p_dust_5, &
                             p_qv, seas_opt, dust_opt, chem_opt,               &
                             DUST_OPT_AFWA, DUST_OPT_FENGSHA, DUST_OPT_GOCART, SEAS_OPT_NONE, &
                             num_chem, num_moist   


  use dust_data_mod, only : den_dust, reff_dust, dyn_visc
  use seas_data_mod, only : den_seas, reff_seas

  implicit none

CONTAINS


SUBROUTINE settling_gocart_driver(dt,t_phy,moist,            &
                                  chem_arr,rho_phy,dz8w,     &
                                  p8w,p_phy,sedim,           &
                                  area,                      &
                                  kms,kme,kts,kte)

  IMPLICIT NONE

  INTEGER,      INTENT(IN   ) ::  kms,kme,kts,kte

  REAL(kind_chem), DIMENSION( kms:kme,num_moist ),                &
         INTENT(IN ) ::                                   moist
  REAL(kind_chem), DIMENSION( kms:kme, num_chem ),                 &
         INTENT(INOUT ) ::                                chem_arr
  REAL(kind_chem), INTENT(IN ) :: area
  REAL(kind_chem), DIMENSION( kms:kme ),                        &
          INTENT(IN   ) ::  t_phy,p_phy,dz8w,p8w,rho_phy
  REAL(kind_chem), INTENT(IN   ) :: dt

  REAL(kind_chem), DIMENSION( num_chem ), INTENT(OUT  ) :: sedim

  integer :: nv,i,j,k,kk,lmx,iseas,idust
  real(kind_chem), DIMENSION (1,1,kte-kts+1) :: tmp,airden,airmas,p_mid,delz,rh
  real(kind_chem), DIMENSION (1,1,kte-kts+1,5) :: dust
  real(kind_chem), DIMENSION (1,1,kte-kts+1,5) :: sea_salt
  real(kind_chem), DIMENSION (kme,num_chem) :: chem_before
  real(kind_chem), dimension (1:5) :: maxdust,maxseas
!
! bstl is for budgets
!
! real(kind_chem), DIMENSION (5), PARAMETER ::
! den_dust(5)=(/2500.,2650.,2650.,2650.,2650./)
! real(kind_chem), DIMENSION (5), PARAMETER ::
! reff_dust(5)=(/0.73D-6,1.4D-6,2.4D-6,4.5D-6,8.0D-6/)
! real(kind_chem), DIMENSION (4), PARAMETER ::
! den_seas(4)=(/2200.,2200.,2200.,2290./)
! real(kind_chem), DIMENSION (4), PARAMETER ::
! reff_seas(4)=(/0.30D-6,1.00D-6,3.25D-6,7.50D-6/)
  real(kind_chem), DIMENSION (5) :: bstl_dust
  real(kind_chem), DIMENSION (5) :: bstl_seas
  real(kind_chem) conver,converi
  real(kind_chem),parameter::max_default=0.

  sedim = 0.
!      conver=1.e-9*mwdry
!      converi=1.e9/mwdry
       conver=1.e-9
       converi=1.e9
       lmx=kte-kts+1

!
! run with all GOCART variables, GOCART sort of HEAVY!
!
! 
! initialize some met stuff
!
          kk=0
          bstl_dust(:)=0.
          bstl_seas(:)=0.
          do k=kts,kte
          kk=kk+1
          p_mid(1,1,kk)=.01*p_phy(kte-k+kts)
          delz(1,1,kk)=dz8w(kte-k+kts)
          airmas(1,1,kk)=-(p8w(k+1)-p8w(k))*area/grav
          airden(1,1,kk)=rho_phy(k)
          tmp(1,1,kk)=t_phy(k)
          rh(1,1,kk) = .95
          rh(1,1,kk) = MIN( .95, moist(k,p_qv) / &
               (3.80*exp(17.27*(t_phy(k)-273.)/ &
               (t_phy(k)-36.))/(.01*p_phy(k))))
          rh(1,1,kk)=max(1.0D-1,rh(1,1,kk))
          do nv = 1, num_chem
            chem_before(k,nv) =  chem_arr(k,nv)
          enddo
          enddo
!
! max dust in column
!
          if((dust_opt == DUST_OPT_GOCART) .or. &
               (dust_opt == DUST_OPT_AFWA) .or. &
               (dust_opt == DUST_OPT_FENGSHA)) then
          iseas=0
          idust=1
          maxdust(:)=0.
          kk=0
          do nv = p_dust_1,p_dust_5
            kk=kk+1
            do k=kts,kte
              if(chem_arr(k,nv).gt.maxdust(kk)) maxdust(kk)=chem_arr(k,nv)
            enddo
          enddo
          kk=0
          do k=kts,kte
            kk=kk+1
            dust(1,1,kk,1)=chem_arr(k,p_dust_1)*conver
            dust(1,1,kk,2)=chem_arr(k,p_dust_2)*conver
            dust(1,1,kk,3)=chem_arr(k,p_dust_3)*conver
            dust(1,1,kk,4)=chem_arr(k,p_dust_4)*conver
            dust(1,1,kk,5)=chem_arr(k,p_dust_5)*conver
          enddo


          call settling(1, 1, lmx, 5,grav,dyn_visc, &
                    dust, tmp, p_mid, delz, airmas, &
                    den_dust, reff_dust, dt, bstl_dust, rh, idust, iseas,airden)

           kk = 0
          do k = kts,kte
             kk = kk+1
             chem_arr(k,p_dust_1)=dust(1,1,kk,1)*converi          ! dust for size bin 1 [ug/kg]
             chem_arr(k,p_dust_2)=dust(1,1,kk,2)*converi          ! ...
             chem_arr(k,p_dust_3)=dust(1,1,kk,3)*converi          ! ...
             chem_arr(k,p_dust_4)=dust(1,1,kk,4)*converi          ! ...
             chem_arr(k,p_dust_5)=dust(1,1,kk,5)*converi          ! dust for size bin 5 (dust_opt 3: for all size bins) [ug/kg]
          enddo
          endif ! dust_opt 
!
!
!
          if(seas_opt /= SEAS_OPT_NONE) then
          iseas=1
          idust=0
          kk=0
          do k=kts,kte
          kk=kk+1
             sea_salt(1,1,kk,1)=chem_arr(k,p_seas_1)*conver
             sea_salt(1,1,kk,2)=chem_arr(k,p_seas_2)*conver
             sea_salt(1,1,kk,3)=chem_arr(k,p_seas_3)*conver
             sea_salt(1,1,kk,4)=chem_arr(k,p_seas_4)*conver
             sea_salt(1,1,kk,5)=chem_arr(k,p_seas_5)*conver
          enddo
!
! max seasalt in column
!
          maxseas(:)=0.
          kk=0
          do nv = p_seas_1,p_seas_5
             kk=kk+1
             do k=kts,kte
                 if(chem_arr(k,nv).gt.maxseas(kk)) maxseas(kk)=chem_arr(k,nv)
             enddo
          enddo
             call settling(1, 1, lmx, 5, grav,dyn_visc,&
                    sea_salt, tmp, p_mid, delz, airmas, &
                    den_seas, reff_seas, dt, bstl_seas, rh, iseas, iseas,airden)
          kk=0
          do k=kts,kte
             kk=kk+1
             chem_arr(k,p_seas_1)=sea_salt(1,1,kk,1)*converi
             chem_arr(k,p_seas_2)=sea_salt(1,1,kk,2)*converi
             chem_arr(k,p_seas_3)=sea_salt(1,1,kk,3)*converi
             chem_arr(k,p_seas_4)=sea_salt(1,1,kk,4)*converi
             chem_arr(k,p_seas_5)=sea_salt(1,1,kk,5)*converi
          enddo


          endif   ! end seasopt==1


          do nv = 1, num_chem
            do k = kts,kte
              sedim(nv) = sedim(nv)+(chem_before(k,nv) - chem_arr(k,nv))*p8w(k)/grav
            enddo
            sedim(nv) = sedim(nv) / dt  !ug/m2/s
          enddo
!
!
!
!
END SUBROUTINE settling_gocart_driver


          subroutine settling(imx,jmx, lmx, nmx,g0,dyn_visc, &
                    tc, tmp, p_mid, delz, airmas, &
                    den, reff, dt, bstl, rh, idust, iseas,airden)
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

  INTEGER, INTENT(IN) :: imx, jmx, lmx, nmx,iseas,idust
  INTEGER :: ntdt
  REAL(kind_chem), INTENT(IN) :: dt,g0,dyn_visc
  REAL(kind_chem),    INTENT(IN) :: tmp(imx,jmx,lmx), delz(imx,jmx,lmx),  &
                         airmas(imx,jmx,lmx), rh(imx,jmx,lmx), &
                         den(nmx), reff(nmx),p_mid(imx,jmx,lmx),&
                         airden(imx,jmx,lmx)
  REAL(kind_chem), INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL(kind_chem), INTENT(OUT)   :: bstl(imx,jmx,nmx)

  REAL(kind_chem)    :: tc1(imx,jmx,lmx,nmx), dt_settl(nmx), rcm(nmx), rho(nmx)
  INTEGER :: ndt_settl(nmx)
  REAL(kind_chem)    :: dzmin, vsettl, dtmax, rhb, rwet(nmx), ratio_r(nmx)
  REAL(kind_chem)    :: c_stokes, free_path, c_cun, viscosity,  growth_fac
  REAL(kind_chem)    :: vd_cor(lmx),vd_wk1 
  INTEGER :: k, n, i, j, l, l2
  REAL(kind_chem)    :: transfer_to_below_level,temp_tc
  ! for sea-salt:
  REAL(kind_chem), PARAMETER :: c1=0.7674, c2=3.079, c3=2.573E-11, c4=-1.424 

  ! for OMP:
  REAL(kind_chem) :: rwet_priv(nmx), rho_priv(nmx)

  ! executable statements

  bstl = 0._kind_chem

! IF (type) /= 'dust' .AND. TRIM(aero_type) /= 'sea_salt') RETURN
  if(idust.ne.1.and.iseas.ne.1)return

!!!  WHERE (tc(:,:,:,:) < 0.0) tc(:,:,:,:) = 1.0E-32

  dzmin = MINVAL(delz(:,:,:))
  IF (idust == 1)     growth_fac = 1.0
  IF (iseas == 1)     growth_fac = 3.0

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

     ! Determine the maximum time-step satisying the CFL condition:
     ! dt <= (dz)_min / v_settl
     ntdt=INT(dt)
     dtmax = dzmin / vsettl
     ndt_settl(k) = MAX( 1, INT( ntdt /dtmax) )
     ! limit maximum number of iterations
     IF (ndt_settl(k) > 12) ndt_settl(k) = 12
     dt_settl(k) = REAL(ntdt) / REAL(ndt_settl(k))

     ! Particles radius in centimeters
     IF (iseas.eq.1)rcm(k) = reff(k)*100.0
     IF (idust.eq.1)then
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
        IF (idust.eq.1) THEN
           rwet_priv(k) = rwet(k)
           rho_priv(k)  = rho(k)
        END IF

        DO n = 1,ndt_settl(k)

           ! Solve each vertical layer successively (layer l)
        transfer_to_below_level=0
 
           DO l = lmx,1,-1
              l2 = lmx - l + 1

!           DO j = 1,jmx
              DO i = 1,imx

                 IF (iseas.eq.1) THEN
                    rhb = MIN(9.9D-1, rh(i,j,l))  
                    ! Aerosol growth with relative humidity (Gerber, 1985)
! td 
! changed to LOG10
                    rwet_priv(k) = 0.01*(c1*rcm(k)**c2/(c3*rcm(k)**c4 - &
                         LOG10(rhb)) + rcm(k)**3)**0.33
                    ratio_r(k) = (reff(k)/rwet_priv(k))**3.0
                 END IF

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
                 IF (iseas.eq.1) THEN
                    rho_priv(k) = ratio_r(k)*den(k) + (1.0 - ratio_r(k))*1000.0
                 END IF

                 vd_cor(l) = 2.0/9.0*g0*rho_priv(k)*rwet_priv(k)**2/viscosity

            ! Update mixing ratio; order of delz: top->sfc
            temp_tc=tc(i,j,l,k)      !temp_tc - for temporal storage [ug/kg]            
            vd_wk1 = dt_settl(k)*vd_cor(l)/delz(i,j,l2)   !fraction to leave level

            tc(i,j,l,k)   =  tc(i,j,l,k)*(1.- vd_wk1)+transfer_to_below_level ! [ug/kg]
            
           if (l.gt.1) transfer_to_below_level =(temp_tc*vd_wk1)*((delz(i,j,l2) &
                   *airden(i,j,l))/(delz(i,j,l2+1)*airden(i,j,l-1)))          ! [ug/kg]

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
           bstl(i,j,n) = 0._kind_chem
           DO l = 1,lmx
              IF (tc(i,j,l,n) < 0.0) tc(i,j,l,n) = 1.0D-32
              bstl(i,j,n) = bstl(i,j,n) + &
                   (tc(i,j,l,n) - tc1(i,j,l,n)) * airmas(i,j,l)
           END DO
        END DO
     END DO
  END DO
  
END SUBROUTINE settling

end module gocart_settling_mod
