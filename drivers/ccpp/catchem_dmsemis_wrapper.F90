!>\file catchem_dmsemis_wrapper.F90
!! This file is GSDChem dms emission wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 06/2020
!! Revision History:
!! 05/2023, Restructure for CATChem, Jian.He@noaa.gov

 module catchem_dmsemis_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config
   use gocart_dmsemis_mod
   use plume_rise_mod

   implicit none

   private

   public :: catchem_dmsemis_wrapper_init, catchem_dmsemis_wrapper_run, catchem_dmsemis_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine catchem_dmsemis_wrapper_init()
      end subroutine catchem_dmsemis_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_dmsemis_wrapper_finalize Argument Table
!!
      subroutine catchem_dmsemis_wrapper_finalize()
      end subroutine catchem_dmsemis_wrapper_finalize

!> \defgroup catchem_group CATChem dmsemis wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_dmsemis_wrapper CATChem dmsemis wrapper Module  
!> \ingroup catchem_dmsemis_group
!! This is the CATChem dmsemis wrapper Module
!! \section arg_table_catchem_dmsemis_wrapper_run Argument Table
!! \htmlinclude catchem_dmsemis_wrapper_run.html
!!
!>\section catchem_dmsemis_wrapper CATChem Scheme General Algorithm
!> @{
    subroutine catchem_dmsemis_wrapper_run(im, kte, kme, dt, garea,       &
                   land, u10m, v10m, tskin,                                &
                   pr3d, ph3d,phl3d, prl3d, tk3d, us3d, vs3d, spechum,     &
                   vegtype, soiltyp,                                       & 
                   emi_in, ntrac, ntdms, gq0, qgrs, dmsemis_opt_in,        &
                   errmsg,errflg)

    implicit none


    integer,        intent(in) :: im,kte,kme
    integer,        intent(in) :: ntrac
    integer,        intent(in) :: ntdms
    real(kind_phys),intent(in) :: dt

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    integer, dimension(im), intent(in) :: land, vegtype, soiltyp        
    real(kind_phys), dimension(im,10), intent(in) :: emi_in
    real(kind_phys), dimension(im),    intent(in) :: u10m, v10m, garea, tskin
    real(kind_phys), dimension(im,kme),intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte),intent(in) :: phl3d, prl3d, tk3d, us3d, vs3d, spechum
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
    integer,           intent(in) :: dmsemis_opt_in
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy, u_phy, v_phy,       &
                     dz8w, p8w, rho_phy

    real(kind_phys), dimension(ims:im, jms:jme) :: u10, v10, tsk, xland, dxy

    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
    real(kind_phys), dimension(ims:im, jms:jme) :: dms_0
    integer,         dimension(ims:im, jms:jme) :: isltyp, ivgtyp
    integer :: ide, ime, ite, kde

    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg

!>-- local variables
    integer :: i, j, jp, k, kp, n
  

    errmsg = ''
    errflg = 0

    dmsemis_opt       = dmsemis_opt_in

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
   !ppm2ugkg(p_so2 ) = 1.e+03_kind_phys * mw_so2_aer / mwdry
    ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry


!>- get ready for chemistry run
    call catchem_dmsemis_prep(                                         &
        u10m,v10m,land,garea,tskin,                                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,                   &
        vegtype,soiltyp,emi_in,u10,v10,tsk,xland,dxy,                   &
        rri,t_phy,u_phy,v_phy,rho_phy,dz8w,p8w,                         &
        ntdms,ntrac,gq0,num_chem,ppm2ugkg,                              &
        chem,ivgtyp,isltyp,dms_0,                                       &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)


    if (dmsemis_opt == DMSE_OPT_ENABLE) then
      call gocart_dmsemis(dt,rri,t_phy,u_phy,v_phy,                     &
         chem,rho_phy,dz8w,u10,v10,p8w,dms_0,tsk,                       &
         ivgtyp,isltyp,xland,dxy,g,mwdry,                               &
         num_chem,p_dms,                                                &
         ids,ide, jds,jde, kds,kde,                                     &
         ims,ime, jms,jme, kms,kme,                                     &
         its,ite, jts,jte, kts,kte)
    endif


    ! -- put chem stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       gq0(i,k,ntdms  )=ppm2ugkg(p_dms   ) * max(epsilc,chem(i,k,1,p_dms)) 
     enddo
    enddo

    do k=kts,kte
     do i=its,ite
       qgrs(i,k,ntdms )=gq0(i,k,ntdms  )
     enddo
    enddo

!
   end subroutine catchem_dmsemis_wrapper_run
!> @}
  subroutine catchem_dmsemis_prep(                                    &
        u10m,v10m,land,garea,ts2d,                                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,                  &
        vegtype,soiltyp,emi_in,u10,v10,tsk,xland,dxy,                  &
        rri,t_phy,u_phy,v_phy,rho_phy,dz8w,p8w,                        &
        ntdms,ntrac,gq0,num_chem,ppm2ugkg,                             &
        chem,ivgtyp,isltyp,dms_0,                                      &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)


    !FV3 input variables
    integer, dimension(ims:ime), intent(in) :: land, vegtype, soiltyp
    integer, intent(in) :: ntrac
    integer,        intent(in) :: ntdms
    real(kind=kind_phys), dimension(ims:ime), intent(in) ::                & 
         u10m, v10m, garea, ts2d
    real(kind=kind_phys), dimension(ims:ime,    10),   intent(in) :: emi_in
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) ::     &
         pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,us3d,vs3d,spechum
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !GSD Chem variables
    integer,intent(in) ::  num_chem
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg

    
    integer,dimension(ims:ime, jms:jme), intent(out) :: isltyp, ivgtyp
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, u_phy, v_phy, rho_phy, dz8w, p8w
    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) ::              &
         u10, v10, tsk, xland, dxy,     &
         dms_0
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem


    ! -- local variables
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: z_at_w, p_phy
    integer i,ip,j,jp,k,kp,kk,kkp,l,ll,n

    ! -- initialize output arrays
    isltyp         = 0._kind_phys
    ivgtyp         = 0._kind_phys
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    u_phy          = 0._kind_phys
    v_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    u10            = 0._kind_phys
    v10            = 0._kind_phys
    tsk            = 0._kind_phys
    xland          = 0._kind_phys
    dxy            = 0._kind_phys
    dms_0          = 0._kind_phys
    chem           = 0._kind_phys


    do i=its,ite
     u10  (i,1)=u10m (i)
     v10  (i,1)=v10m (i)
     tsk  (i,1)=ts2d (i)
     dxy  (i,1)=garea(i)
     xland(i,1)=real(land(i))
     ivgtyp (i,1)=vegtype(i)
     isltyp (i,1)=soiltyp(i)
     dms_0(i,1  )=emi_in(i,7) ! --dm0
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
          v_phy(i,k,j)=vs3d(ip,kkp)
          rho_phy(i,k,j)=p_phy(i,k,j)/(287.04*t_phy(i,k,j)*(1.+.608*spechum(ip,kkp)))
          rri(i,k,j)=1./rho_phy(i,k,j)
          !--
        enddo
      enddo
    enddo

    do k=kms,kte
     do i=ims,ime
       chem(i,k,jts,p_dms   )=max(epsilc,gq0(i,k,ntdms  )/ppm2ugkg(p_dms))
     enddo
    enddo


  end subroutine catchem_dmsemis_prep
!> @}
  end module catchem_dmsemis_wrapper
