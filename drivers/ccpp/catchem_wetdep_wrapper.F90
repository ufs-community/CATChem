!>\file catchem_wetdep_wrapper.F90
!! This file is GSDChem large-scale wet deposition wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 06/2020
!! Kate.Zhang@noaa.gov 02/2023
!! Revision History:
!! 05/2023, Restructure for CATChem, Jian.He@noaa.gov
!! Kate.Zhang@noaa.gov 11/2023

 module catchem_wetdep_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config
   use wetdep_ls_mod
   use gocart_aerosols_mod
   use dust_data_mod
   use gocart_diag_mod

   implicit none

   private

   public :: catchem_wetdep_wrapper_init, catchem_wetdep_wrapper_run, catchem_wetdep_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine catchem_wetdep_wrapper_init()
      
      end subroutine catchem_wetdep_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_wetdep_wrapper_finalize Argument Table
!!
      subroutine catchem_wetdep_wrapper_finalize()
      end subroutine catchem_wetdep_wrapper_finalize

!> \defgroup catchem_group CATChem wetdep wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_wetdep_wrapper CATChem wetdep wrapper Module  
!> \ingroup catchem_wetdep_group
!! This is the CATChem wetdep wrapper Module
!! \section arg_table_catchem_wetdep_wrapper_run Argument Table
!! \htmlinclude catchem_wetdep_wrapper_run.html
!!
!>\section catchem_wetdep_wrapper CATChem Scheme General Algorithm
!> @{
    subroutine catchem_wetdep_wrapper_run(im, kte, kme, ktau, dt,       &
                   imp_physics, imp_physics_gfdl, imp_physics_thompson, &
                   rain_cplchm, rainc_cpl,rlat,                         &
                   pr3d, ph3d,phl3d, prl3d, tk3d, us3d, vs3d, spechum,  &
                   w, dqdt, ntrac,ntchmdiag,                            &
                   ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                     &
                   ntbc1,ntbc2,ntoc1,ntoc2,                             &
                   ntss1,ntss2,ntss3,ntss4,ntss5,                       &
                   ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,      &
                   gq0,qgrs,wetdpl,wetdep_ls_opt_in,                    &
                   errmsg,errflg)

    implicit none


    integer,        intent(in) :: im,kte,kme,ktau,ntchmdiag
    integer,        intent(in) :: ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
    integer,        intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer,        intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    real(kind_phys),intent(in) :: dt

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    real(kind_phys), dimension(im),     intent(in) :: rain_cplchm,rainc_cpl,rlat
    real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte), intent(in) :: phl3d, prl3d, tk3d,        &
                us3d, vs3d, spechum, w, dqdt
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
    real(kind_phys), dimension(im,ntchmdiag), intent(inout) :: wetdpl
    integer,           intent(in) :: wetdep_ls_opt_in
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson
    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy, u_phy, v_phy,       &
                     p_phy, z_at_w, dz8w, p8w, t8w, rho_phy, vvel, dqdti

    real(kind_phys), dimension(ims:im, jms:jme) :: rcav, rnav,xlat

!>- vapor & chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist 
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
    real(kind_phys), dimension(ims:im, jms:jme, 1:num_chem )  :: var_rmv

    integer :: ide, ime, ite, kde

    real(kind_phys) :: dtstep
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg

    ! -- output tracers
    real(kind_phys), dimension(im, 1, ntchmdiag, 4) :: trdf


!>-- local variables
    integer :: i, j, jp, k, kp, n
  

    errmsg = ''
    errflg = 0

    wetdep_ls_opt     = wetdep_ls_opt_in

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
   !ppm2ugkg(p_so2 ) = 1.e+03_kind_phys * mw_so2_aer / mwdry
    ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry

    ! -- initialize large-sacle wet depostion
!    if (ktau==1) then
!     call dep_wet_ls_init()
!    endif

    ! -- set control flags

    ! -- compute accumulated large-scale and convective rainfall since last call
    if (ktau > 1) then
      dtstep = call_chemistry * dt
    else
      dtstep = dt
    end if

    ! -- compute incremental convective and large-scale rainfall
    do i=its,ite
     rcav(i,1)=max(rainc_cpl(i)*1000.              , 0.) ! meter to mm
     rnav(i,1)=max((rain_cplchm(i)-rainc_cpl(i))*1000., 0.) ! meter to mm
    enddo

!!!

!>- get ready for chemistry run
    call catchem_prep_wetdep(ktau,dtstep,                               &
        imp_physics, imp_physics_gfdl, imp_physics_thompson,            &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w, dqdt,           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                   &
        t8w,dqdti,z_at_w,vvel,rlat,xlat,                                &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                                &
        ntbc1,ntbc2,ntoc1,ntoc2,                                        &
        ntss1,ntss2,ntss3,ntss4,ntss5,                                  &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                 &
        ntrac,gq0,num_chem, num_moist,                                  &
        ppm2ugkg,moist,chem,                                            &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)

    ! -- ls wet deposition
    var_rmv(:,:,:)=0.

    select case (wetdep_ls_opt)
      case (WDLS_OPT_GSD)
        do j=jts,jte
          do i=its,ite
            call wetdep_ls(dt,chem(i,:,j,:),rnav(i,j),moist(i,:,j,:),   &
                        rho_phy(i,:,j),var_rmv(i,j,:),xlat(i,j),        &
                        p_qc,p_qi,dz8w(i,:,j),vvel(i,:,j),              &
                        kms,kme,kts,kte)
          enddo
        end do

      case (WDLS_OPT_NGAC)
        do j=jts,jte
          do i=its,ite
            call WetRemovalGOCART(kts,kte, 1,1, dt,   &
                               var_rmv(i,j,:),chem(i,:,j,:),    &
                               p_phy(i,:,j),t_phy(i,:,j),       &
                               rho_phy(i,:,j),dqdti(i,:,j),     &
                               rcav(i,j),rnav(i,j),             &
                               kms,kme)
          enddo
        enddo

         !if (chem_rc_check(localrc, msg="Failure in NGAC wet removal scheme", &
         !  file=__FILE__, line=__LINE__, rc=rc)) return
      case default
        ! -- no further option implemented
        errmsg = 'Logic error in catchem_wetdep_wrapper_run: invalid wdls_opt'
        errflg = 1
        return
    end select


    ! -- put chem stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       gq0(i,k,ntso2  )=ppm2ugkg(p_so2   ) * max(epsilc,chem(i,k,1,p_so2))
       gq0(i,k,ntsulf )=ppm2ugkg(p_sulf  ) * max(epsilc,chem(i,k,1,p_sulf))
       gq0(i,k,ntdms  )=ppm2ugkg(p_dms   ) * max(epsilc,chem(i,k,1,p_dms)) 
       gq0(i,k,ntmsa  )=ppm2ugkg(p_msa   ) * max(epsilc,chem(i,k,1,p_msa))
       gq0(i,k,ntpp25 )=ppm2ugkg(p_p25   ) * max(epsilc,chem(i,k,1,p_p25))
       gq0(i,k,ntbc1  )=ppm2ugkg(p_bc1   ) * max(epsilc,chem(i,k,1,p_bc1))
       gq0(i,k,ntbc2  )=ppm2ugkg(p_bc2   ) * max(epsilc,chem(i,k,1,p_bc2))
       gq0(i,k,ntoc1  )=ppm2ugkg(p_oc1   ) * max(epsilc,chem(i,k,1,p_oc1))
       gq0(i,k,ntoc2  )=ppm2ugkg(p_oc2   ) * max(epsilc,chem(i,k,1,p_oc2))
       gq0(i,k,ntdust1)=ppm2ugkg(p_dust_1) * max(epsilc,chem(i,k,1,p_dust_1))
       gq0(i,k,ntdust2)=ppm2ugkg(p_dust_2) * max(epsilc,chem(i,k,1,p_dust_2))
       gq0(i,k,ntdust3)=ppm2ugkg(p_dust_3) * max(epsilc,chem(i,k,1,p_dust_3))
       gq0(i,k,ntdust4)=ppm2ugkg(p_dust_4) * max(epsilc,chem(i,k,1,p_dust_4))
       gq0(i,k,ntdust5)=ppm2ugkg(p_dust_5) * max(epsilc,chem(i,k,1,p_dust_5))
       gq0(i,k,ntss1  )=ppm2ugkg(p_seas_1) * max(epsilc,chem(i,k,1,p_seas_1))
       gq0(i,k,ntss2  )=ppm2ugkg(p_seas_2) * max(epsilc,chem(i,k,1,p_seas_2))
       gq0(i,k,ntss3  )=ppm2ugkg(p_seas_3) * max(epsilc,chem(i,k,1,p_seas_3))
       gq0(i,k,ntss4  )=ppm2ugkg(p_seas_4) * max(epsilc,chem(i,k,1,p_seas_4))
       gq0(i,k,ntss5  )=ppm2ugkg(p_seas_5) * max(epsilc,chem(i,k,1,p_seas_5))
       gq0(i,k,ntpp10 )=ppm2ugkg(p_p10   ) * max(epsilc,chem(i,k,1,p_p10))
     enddo
    enddo

    do k=kts,kte
     do i=its,ite
       qgrs(i,k,ntso2  )=gq0(i,k,ntso2  )
       qgrs(i,k,ntsulf )=gq0(i,k,ntsulf )
       qgrs(i,k,ntdms  )=gq0(i,k,ntdms  )
       qgrs(i,k,ntmsa  )=gq0(i,k,ntmsa  )
       qgrs(i,k,ntpp25 )=gq0(i,k,ntpp25 )
       qgrs(i,k,ntbc1  )=gq0(i,k,ntbc1  )
       qgrs(i,k,ntbc2  )=gq0(i,k,ntbc2  )
       qgrs(i,k,ntoc1  )=gq0(i,k,ntoc1  )
       qgrs(i,k,ntoc2  )=gq0(i,k,ntoc2  )
       qgrs(i,k,ntdust1)=gq0(i,k,ntdust1)
       qgrs(i,k,ntdust2)=gq0(i,k,ntdust2)
       qgrs(i,k,ntdust3)=gq0(i,k,ntdust3)
       qgrs(i,k,ntdust4)=gq0(i,k,ntdust4)
       qgrs(i,k,ntdust5)=gq0(i,k,ntdust5)
       qgrs(i,k,ntss1  )=gq0(i,k,ntss1  )
       qgrs(i,k,ntss2  )=gq0(i,k,ntss2  )
       qgrs(i,k,ntss3  )=gq0(i,k,ntss3  )
       qgrs(i,k,ntss4  )=gq0(i,k,ntss4  )
       qgrs(i,k,ntss5  )=gq0(i,k,ntss5  )
       qgrs(i,k,ntpp10 )=gq0(i,k,ntpp10 )
     enddo
    enddo

    ! -- output large-scale wet deposition
    call gocart_diag_store(3, var_rmv, trdf)

     wetdpl (:,:)=trdf(:,1,:,3)

!
   end subroutine catchem_wetdep_wrapper_run
!> @}

  subroutine catchem_prep_wetdep(ktau,dtstep,                          &
        imp_physics, imp_physics_gfdl, imp_physics_thompson,           &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,dqdt,           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                  &
        t8w,dqdti,z_at_w,vvel,rlat,xlat,                               &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                               &
        ntbc1,ntbc2,ntoc1,ntoc2,                                       &
        ntss1,ntss2,ntss3,ntss4,ntss5,                                 &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                &
        ntrac,gq0,num_chem, num_moist,                                 &
        ppm2ugkg,moist,chem,                                           &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)

    !Chem input configuration
    integer, intent(in) :: ktau
    real(kind=kind_phys), intent(in) :: dtstep

    !FV3 input variables
    integer, intent(in) :: ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
    integer, intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer, intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson
    integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    real(kind=kind_phys), dimension(ims:ime), intent(in) ::rlat
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,dqdt
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !GSD Chem variables
    integer,intent(in) ::  num_chem, num_moist
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg
    
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, u_phy, v_phy, p_phy, rho_phy, dz8w, p8w, t8w, vvel, dqdti
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, 1:num_moist), intent(out) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w
    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) :: xlat
    ! -- local variables
!   real(kind=kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: p_phy
    real(kind_phys) ::  factor,factor2,pu,pl,aln,pwant
    real(kind_phys) ::  xhour,xmin,xlonn,xtime,real_time
    real(kind_phys), DIMENSION (1,1) :: sza,cosszax
    integer i,ip,j,jp,k,kp,kk,kkp,nv,jmax,jmaxi,l,ll,n,ndystep,ixhour

    ! -- initialize output arrays
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    u_phy          = 0._kind_phys
    v_phy          = 0._kind_phys
    p_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    t8w            = 0._kind_phys
    vvel           = 0._kind_phys
    dqdti          = 0._kind_phys
    moist          = 0._kind_phys  
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys
    xlat           = 0._kind_phys


    do i=its,ite
    xlat (i,1)=rlat(i)*180./pi
    enddo


    do j=jts,jte
      jp = j - jts + 1
      do i=its,ite
         ip = i - its + 1
         z_at_w(i,kts,j)=max(0.,ph3d(ip,1)/g)
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=abs(ph3d(ip,kp+1)-ph3d(ip,kp))/g
          z_at_w(i,k+1,j)=z_at_w(i,k,j)+dz8w(i,k,j)
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          p8w(i,k,j)=pr3d(ip,kp)
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
        kk=min(k,kte)
        kkp = kk - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=z_at_w(i,kk+1,j)-z_at_w(i,kk,j)
          t_phy(i,k,j)=tk3d(ip,kkp)
          p_phy(i,k,j)=prl3d(ip,kkp)
          u_phy(i,k,j)=us3d(ip,kkp)
          dqdti(i,k,j)=dqdt(ip,kkp)
          v_phy(i,k,j)=vs3d(ip,kkp)
          rho_phy(i,k,j)=p_phy(i,k,j)/(287.04*t_phy(i,k,j)*(1.+.608*spechum(ip,kkp)))
          rri(i,k,j)=1./rho_phy(i,k,j)
          vvel(i,k,j)=-w(ip,kkp)*rri(i,k,j)/g 
          moist(i,k,j,1)=gq0(ip,kkp,p_atm_shum)
          moist(i,k,j,2)=gq0(ip,kkp,p_atm_cldq)
          if (moist(i,k,j,2) < 1.e-30) moist(i,k,j,2)=0.
            if (imp_physics==imp_physics_thompson) then
            moist(i,k,j,3)=gq0(ip,kkp,p_atm_cldi)
            else if (imp_physics==imp_physics_gfdl) then
            moist(i,k,j,3)=gq0(ip,kkp,p_atm_cldii)
            endif
          if(moist(i,k,j,3) < 1.e-30) moist(i,k,j,3)=0.
          !--
        enddo
      enddo
    enddo

    do j=jts,jte
      do k=2,kte
        do i=its,ite
          t8w(i,k,j)=.5*(t_phy(i,k,j)+t_phy(i,k-1,j))
        enddo
      enddo
    enddo

    ! -- only used in phtolysis....
    do j=jts,jte
      do i=its,ite
        t8w(i,1,j)=t_phy(i,1,j)
        t8w(i,kte+1,j)=t_phy(i,kte,j)
      enddo
    enddo

 
    do k=kms,kte
     do i=ims,ime
       chem(i,k,jts,p_so2   )=max(epsilc,gq0(i,k,ntso2  )/ppm2ugkg(p_so2))
       chem(i,k,jts,p_sulf  )=max(epsilc,gq0(i,k,ntsulf )/ppm2ugkg(p_sulf))
       chem(i,k,jts,p_dms   )=max(epsilc,gq0(i,k,ntdms  )/ppm2ugkg(p_dms))
       chem(i,k,jts,p_msa   )=max(epsilc,gq0(i,k,ntmsa  )/ppm2ugkg(p_msa))
       chem(i,k,jts,p_p25   )=max(epsilc,gq0(i,k,ntpp25 )/ppm2ugkg(p_p25))
       chem(i,k,jts,p_bc1   )=max(epsilc,gq0(i,k,ntbc1  )/ppm2ugkg(p_bc1))
       chem(i,k,jts,p_bc2   )=max(epsilc,gq0(i,k,ntbc2  )/ppm2ugkg(p_bc2))
       chem(i,k,jts,p_oc1   )=max(epsilc,gq0(i,k,ntoc1  )/ppm2ugkg(p_oc1))
       chem(i,k,jts,p_oc2   )=max(epsilc,gq0(i,k,ntoc2  )/ppm2ugkg(p_oc2))
       chem(i,k,jts,p_dust_1)=max(epsilc,gq0(i,k,ntdust1)/ppm2ugkg(p_dust_1))
       chem(i,k,jts,p_dust_2)=max(epsilc,gq0(i,k,ntdust2)/ppm2ugkg(p_dust_2))
       chem(i,k,jts,p_dust_3)=max(epsilc,gq0(i,k,ntdust3)/ppm2ugkg(p_dust_3))
       chem(i,k,jts,p_dust_4)=max(epsilc,gq0(i,k,ntdust4)/ppm2ugkg(p_dust_4))
       chem(i,k,jts,p_dust_5)=max(epsilc,gq0(i,k,ntdust5)/ppm2ugkg(p_dust_5))
       chem(i,k,jts,p_seas_1)=max(epsilc,gq0(i,k,ntss1  )/ppm2ugkg(p_seas_1))
       chem(i,k,jts,p_seas_2)=max(epsilc,gq0(i,k,ntss2  )/ppm2ugkg(p_seas_2))
       chem(i,k,jts,p_seas_3)=max(epsilc,gq0(i,k,ntss3  )/ppm2ugkg(p_seas_3))
       chem(i,k,jts,p_seas_4)=max(epsilc,gq0(i,k,ntss4  )/ppm2ugkg(p_seas_4))
       chem(i,k,jts,p_seas_5)=max(epsilc,gq0(i,k,ntss5  )/ppm2ugkg(p_seas_5))
       chem(i,k,jts,p_p10   )=max(epsilc,gq0(i,k,ntpp10 )/ppm2ugkg(p_p10))
     enddo
    enddo


  end subroutine catchem_prep_wetdep
!> @}
  end module catchem_wetdep_wrapper
