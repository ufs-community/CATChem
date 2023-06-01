!>\file catchem_drydep_wrapper.F90
!! This file is GSDChem dry deposition wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 06/2020
!! Revision History:
!! 05/2023, Restructure for CATChem, Jian.He@noaa.gov

 module catchem_drydep_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config
   use dep_dry_mod
   use gocart_diag_mod

   implicit none

   private

   public :: catchem_drydep_wrapper_init, catchem_drydep_wrapper_run, catchem_drydep_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine catchem_drydep_wrapper_init()
      end subroutine catchem_drydep_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_drydep_wrapper_finalize Argument Table
!!
      subroutine catchem_drydep_wrapper_finalize()
      end subroutine catchem_drydep_wrapper_finalize

!> \defgroup catchem_group CATChem drydep wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_drydep_wrapper CATChem drydep wrapper Module  
!> \ingroup catchem_drydep_group
!! This is the CATChem drydep wrapper Module
!! \section arg_table_catchem_drydep_wrapper_run Argument Table
!! \htmlinclude catchem_drydep_wrapper_run.html
!!
!> @{
    subroutine catchem_drydep_wrapper_run(im, kte, kme, ktau, dt, land,      &
                   ustar, rlat, rlon, tskin, julian, rainc_cpl, hf2d, pb2d,   &
                   pr3d, ph3d, phl3d, prl3d, tk3d, spechum, exch,             &
                   vegtype, sigmaf, jdate, idat, dswsfc, zorl, snow_cplchm,      & 
                   ntrac,ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                     &
                   ntbc1,ntbc2,ntoc1,ntoc2,ntss1,ntss2,ntss3,ntss4,ntss5,     &
                   ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,            &
                   ntchmdiag,gq0,qgrs,drydep,chem_conv_tr_in,                 &
                   errmsg,errflg)

    implicit none


    integer,        intent(in) :: im,kte,kme,ktau,jdate(8),idat(8)
    integer,        intent(in) :: ntrac,ntchmdiag,ntss1,ntss2,ntss3,ntss4,ntss5
    integer,        intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer,        intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    real(kind_phys),intent(in) :: dt,julian

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    integer, dimension(im), intent(in) :: land, vegtype     
    real(kind_phys), dimension(im), intent(in) :: ustar,                  &
                rlat,rlon, tskin, rainc_cpl,                              &
                hf2d, pb2d, sigmaf, dswsfc, zorl, snow_cplchm
    real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte), intent(in) :: phl3d, prl3d, tk3d, spechum, exch
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
    real(kind_phys), dimension(im,ntchmdiag), intent(inout) :: drydep
    integer,        intent(in) :: chem_conv_tr_in
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy,        &
                     p_phy, z_at_w, dz8w, p8w, t8w, rho_phy, zmid, exch_h 

    real(kind_phys), dimension(ims:im, jms:jme) :: ust, tsk,              &
                     xland, xlat, xlong, rcav, hfx, pbl

!>- chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist 
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem

    real(kind_phys), dimension(ims:im, jms:jme, 1:num_chem )  :: dry_fall
    real(kind_phys), dimension(im, 1, ntchmdiag, 4) :: trdf

    integer :: ide, ime, ite, kde, julday

    real(kind_phys), dimension(ims:im, jms:jme) :: vegfrac, rmol, gsw, znt
    real(kind_phys), dimension(ims:im, jms:jme) :: snowh  
    integer,         dimension(ims:im, jms:jme) :: ivgtyp

    integer :: current_month

    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: ac3, ahno3, anh3, asulf, cor3, h2oai, h2oaj, nu3
    real(kind_phys), dimension(ims:im, jms:jme) :: dep_vel_o3, e_co

    real(kind_phys) :: gmt
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg

!>-- local variables
    integer :: i, j, jp, k, kp, n
  

    errmsg = ''
    errflg = 0

    chem_conv_tr      = chem_conv_tr_in

    h2oai = 0.
    h2oaj = 0.
    nu3   = 0.
    ac3   = 0.
    cor3  = 0.
    asulf = 0.
    ahno3 = 0.
    anh3  = 0.
    e_co  = 0.
    dep_vel_o3 = 0.


    gmt = real(idat(5))
    julday = real(julian)                                       

    current_month=jdate(2)

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
   !ppm2ugkg(p_so2 ) = 1.e+03_kind_phys * mw_so2_aer / mwdry
    ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry

    ! -- compute incremental convective and large-scale rainfall
    do i=its,ite
     rcav(i,1)=max(rainc_cpl(i)*1000.              , 0.) ! meter to mm
    enddo

!!!

!>- get ready for chemistry run
    call catchem_prep_drydep(                                          &
        ustar,land,rlat,rlon,tskin,                                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,spechum,exch,                        &
        vegtype,sigmaf,dswsfc,zorl,snow_cplchm,hf2d,pb2d,                  &
        ust,tsk,xland,xlat,xlong,                                       &
        rri,t_phy,p_phy,rho_phy,dz8w,p8w,                               &
        t8w,exch_h,z_at_w,zmid,                                         &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                                &
        ntbc1,ntbc2,ntoc1,ntoc2,ntss1,ntss2,ntss3,ntss4,ntss5,          &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                 &
        ntrac,gq0,num_chem, num_moist,                                  &
        ppm2ugkg,moist,chem,                                            &
        ivgtyp,vegfrac,rmol,gsw,znt,hfx,pbl,snowh,                      &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)


     !>-- compute dry deposition
     call dry_dep_driver(ktau,dt,julday,current_month,t_phy,p_phy,&
       moist,p8w,rmol,rri,gmt,t8w,rcav,                           &
       chem,rho_phy,dz8w,exch_h,hfx,                              &
       ivgtyp,tsk,gsw,vegfrac,pbl,ust,znt,zmid,z_at_w,            &
       xland,xlat,xlong,h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,     &
       anh3,dry_fall,dep_vel_o3,g,                             &
       e_co,kemit,snowh,numgas,                                   &
       num_chem,num_moist,                                        &
       ids,ide, jds,jde, kds,kde,                                 &
       ims,ime, jms,jme, kms,kme,                                 &
       its,ite, jts,jte, kts,kte)


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

    ! -- output dry deposition
    call gocart_diag_store(2, dry_fall, trdf)

    drydep (:,:)=trdf(:,1,:,2)

!
   end subroutine catchem_drydep_wrapper_run
!> @}
  subroutine catchem_prep_drydep(                                     &
        ustar,land,rlat,rlon,ts2d,                                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,spechum,exch,                       &
        vegtype,sigmaf,dswsfc,zorl,snow_cplchm,hf2d,pb2d,                 &
        ust,tsk,xland,xlat,xlong,                                      &
        rri,t_phy,p_phy,rho_phy,dz8w,p8w,                              &
        t8w,exch_h,z_at_w,zmid,                                        &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                               &
        ntbc1,ntbc2,ntoc1,ntoc2,ntss1,ntss2,ntss3,ntss4,ntss5,         &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                &
        ntrac,gq0,num_chem, num_moist,                                 &
        ppm2ugkg,moist,chem,                                           &
        ivgtyp,vegfrac,rmol,gsw,znt,hfx,pbl,snowh,                     &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)

    !FV3 input variables
    integer, dimension(ims:ime), intent(in) :: land, vegtype
    integer, intent(in) :: ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
    integer, intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer, intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer, intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    real(kind=kind_phys), dimension(ims:ime), intent(in) ::                & 
         ustar, rlat, rlon, ts2d, sigmaf, dswsfc, zorl, snow_cplchm, hf2d, pb2d
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) :: phl3d,tk3d,prl3d,spechum,exch
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !GSD Chem variables
    integer,intent(in) ::  num_chem, num_moist
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg

    
    integer,dimension(ims:ime, jms:jme), intent(out) :: ivgtyp
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, p_phy, rho_phy, dz8w, p8w, t8w, zmid,         &
         exch_h
    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) ::              &
         ust, tsk, xland, xlat, xlong, vegfrac, rmol, gsw, znt, hfx,    &
         pbl, snowh
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w

    ! -- local variables
    integer i,ip,j,jp,k,kp,kk,kkp,nv,l,ll,n


    ! -- initialize output arrays
    ivgtyp         = 0._kind_phys
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    p_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    t8w            = 0._kind_phys
    zmid           = 0._kind_phys
    exch_h         = 0._kind_phys
    ust            = 0._kind_phys
    tsk            = 0._kind_phys
    xland          = 0._kind_phys
    xlat           = 0._kind_phys
    xlong          = 0._kind_phys
    vegfrac        = 0._kind_phys
    rmol           = 0._kind_phys
    gsw            = 0._kind_phys
    znt            = 0._kind_phys
    hfx            = 0._kind_phys
    pbl            = 0._kind_phys
    snowh          = 0._kind_phys
    moist          = 0._kind_phys  
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys


    do i=its,ite
     tsk  (i,1)=ts2d (i)
     ust  (i,1)=ustar(i)
     xland(i,1)=real(land(i))
     xlat (i,1)=rlat(i)*180./pi
     xlong(i,1)=rlon(i)*180./pi
     gsw  (i,1)=dswsfc(i)
     znt  (i,1)=zorl(i)*0.01
     hfx  (i,1)=hf2d(i)
     pbl  (i,1)=pb2d(i)
     snowh(i,1)=snow_cplchm(i)*0.001
     ivgtyp (i,1)=vegtype(i)
     vegfrac(i,1)=sigmaf (i)
    enddo
   
    rmol=0.

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
          rho_phy(i,k,j)=p_phy(i,k,j)/(287.04*t_phy(i,k,j)*(1.+.608*spechum(ip,kkp)))
          rri(i,k,j)=1./rho_phy(i,k,j)
          moist(i,k,j,:)=0.
          moist(i,k,j,1)=gq0(ip,kkp,p_atm_shum)
          if (t_phy(i,k,j) > 265.) then
            moist(i,k,j,2)=gq0(ip,kkp,p_atm_cldq)
            moist(i,k,j,3)=0.
            if (moist(i,k,j,2) < 1.e-8) moist(i,k,j,2)=0.
          else
            moist(i,k,j,2)=0.
            moist(i,k,j,3)=gq0(ip,kkp,p_atm_cldq)
            if(moist(i,k,j,3) < 1.e-8)moist(i,k,j,3)=0.
          endif
          !--
          zmid(i,k,j)=phl3d(ip,kkp)/g
        enddo
      enddo
    enddo

    ! -- the imported atmospheric heat diffusivity is only available up to kte-1
    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte-1
        kkp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          exch_h(i,k,j)=exch(ip,kkp)
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


  end subroutine catchem_prep_drydep

!> @}
  end module catchem_drydep_wrapper
