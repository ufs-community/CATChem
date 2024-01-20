!>\file catchem_rad_wrapper.F90
!! This file is GSDChem radiation wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 06/2020
!! Kate.Zhang@noaa.gov 04/2023
!! Revision History:
!! 05/2023, Restructure for CATChem, Jian.He@noaa.gov

 module catchem_rad_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config
   use opt_mod
   use dust_data_mod

   implicit none

   private

   public :: catchem_rad_wrapper_init, catchem_rad_wrapper_run, catchem_rad_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine catchem_rad_wrapper_init()
      end subroutine catchem_rad_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_rad_wrapper_finalize Argument Table
!!
      subroutine catchem_rad_wrapper_finalize()
      end subroutine catchem_rad_wrapper_finalize

!> \defgroup catchem_group CATChem rad wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_rad_wrapper CATChem rad wrapper Module  
!> \ingroup catchem_rad_group
!! This is the CATChem rad wrapper Module
!! \section arg_table_catchem_rad_wrapper_run Argument Table
!! \htmlinclude catchem_rad_wrapper_run.html
!!
!>\section catchem_rad_wrapper CATChem Scheme General Algorithm
!> @{
    subroutine catchem_rad_wrapper_run(im, kte, kme, ktau, dt,         &
                   ph3d,prl3d, tk3d, spechum,                           &
                   ntrac,ntso2,ntsulf,ntDMS,ntmsa,ntpp25,               &
                   ntbc1,ntbc2,ntoc1,ntoc2,                             &
                   ntss1,ntss2,ntss3,ntss4,ntss5,                       &
                   ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,      &
                   gq0,abem,                                            &
!                   cplchp_rad_opt,lmk,faersw_cpl,                       &
                   chem_opt_in,aer_ra_feedback_in,aer_ra_frq_in,        &
                   errmsg,errflg)

    implicit none


    integer,        intent(in) :: im,kte,kme,ktau
    integer,        intent(in) :: ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
    integer,        intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer,        intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    real(kind_phys),intent(in) :: dt

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    real(kind_phys), dimension(im,kme), intent(in) :: ph3d
    real(kind_phys), dimension(im,kte), intent(in) :: prl3d, tk3d, spechum
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0
    real(kind_phys), dimension(im,12        ), intent(inout) :: abem
!    integer,         intent(in) :: lmk
!    real(kind_phys), dimension(im, lmk, 14, 3),intent(inout) :: faersw_cpl
!    logical, intent(in) :: cplchp_rad_opt
    integer,        intent(in) :: chem_opt_in
    integer,        intent(in) :: aer_ra_feedback_in,aer_ra_frq_in
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, dz8w

    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem

    integer :: ide, ime, ite, kde

    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: h2oai, h2oaj

!>- optical variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: relhum
    real(kind_phys), dimension(ims:im,         jms:jme) :: aod
    real(kind_phys), dimension(ims:im,         jms:jme,5) :: aerodp ! dust, soot, waso, suso, ssam
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:nbands) :: extt
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:nbands) :: ssca
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:nbands) :: asympar
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:4) ::                   &
        tauaersw, gaersw, waersw, bscoefsw,                                        &
        l2aer,  l3aer, l4aer, l5aer, l6aer, l7aer           
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:16) :: tauaerlw
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_ext_coef) :: ext_coeff
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_bscat_coef) :: bscat_coeff
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_asym_par)   :: asym_par
    real(kind_phys), dimension(im) :: aod2d
    real(kind_phys), dimension(im, kte, 1:nbands) :: ext_cof, sscal, asymp

!>-- local variables
    real(kind_phys) :: dtstep
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg
    real(kind_phys) :: curr_secs
    logical :: call_radiation
    logical :: store_arrays
    integer :: nbegin, nv, nvv
    integer :: i, j, jp, k, kp, n
  

    errmsg = ''
    errflg = 0

    chem_opt          = chem_opt_in
    aer_ra_feedback   = aer_ra_feedback_in
    aer_ra_frq        = aer_ra_frq_in

    h2oai = 0.
    h2oaj = 0.
    extt =0.
    ssca   =0.
    asympar=0.
    aod=0.
    aod2d=0.
    aerodp=0.

    curr_secs = ktau * dt

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
   !ppm2ugkg(p_so2 ) = 1.e+03_kind_phys * mw_so2_aer / mwdry
    ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry

    ! -- set control flags
    call_radiation   = (mod(int(curr_secs), max(1, 60*aer_ra_frq)) == 0) .or. (ktau == 1)

    if (ktau > 1) then
      dtstep = call_chemistry * dt
    else
      dtstep = dt
    end if

!!!

!>- get ready for chemistry run
    call catchem_prep_rad(                                             &
        ktau,dtstep,ph3d,tk3d,prl3d,spechum,rri,dz8w,                   &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                                &
        ntbc1,ntbc2,ntoc1,ntoc2,                                        &
        ntss1,ntss2,ntss3,ntss4,ntss5,                                  &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                 &
        ntrac,gq0,num_chem,ppm2ugkg,chem,relhum,                        &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)


    if (call_radiation .and. aer_ra_feedback > 0 ) then
      store_arrays = .false.
      select case (aer_ra_feedback)
        case (1)
          call optical_driver(curr_secs,dtstep,             &
               chem,dz8w,rri,relhum,                        &
               h2oai,h2oaj,                                 &
               tauaersw,gaersw,waersw,bscoefsw,tauaerlw,    &
               l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,         &
               num_chem,chem_opt,ids,ide, jds,jde, kds,kde, &
               ims,ime, jms,jme, kms,kme,                   &
               its,ite, jts,jte, kts,kte)
          call aer_opt_out(aod,dz8w,                        &
               ext_coeff,bscat_coeff,asym_par,              &
               tauaersw,gaersw,waersw,tauaerlw,             &
               num_ext_coef,num_bscat_coef,num_asym_par,    &
               ids,ide, jds,jde, kds,kde,                   &
               ims,ime, jms,jme, kms,kme,                   &
               its,ite, jts,jte, kts,kte)
          call aer_ra(dz8w                                  &
               ,extt,ssca,asympar,nbands                    &
               ,tauaersw,gaersw,waersw,tauaerlw             &
               ,ids,ide, jds,jde, kds,kde                   &
               ,ims,ime, jms,jme, kms,kme                   &
               ,its,ite, jts,jte, kts,kte)
          store_arrays = .true.
        case (2)
          call aero_opt('sw',dz8w,chem                                  &
                   ,rri,relhum,aod                                      &
                   ,extt,ssca,asympar,num_chem                          &
                   ,ids,ide, jds,jde, kds,kde                           &
                   ,ims,ime, jms,jme, kms,kme                           &
                   ,its,ite, jts,jte, kts,kte)
          store_arrays = .true.
        case (3)
          call aero_opt_new('sw',dz8w,chem                              &
                   ,rri,relhum                                          &
                   ,extt,ssca,asympar,num_chem                          &
                   ,ids,ide, jds,jde, kds,kde                           &
                   ,ims,ime, jms,jme, kms,kme                           &
                   ,its,ite, jts,jte, kts,kte                           &
                   ,aod,aerodp)
          store_arrays = .true.
        case default
          ! -- no feedback
      end select
      if (store_arrays) then
        do nv = 1, nbands
          do k = kts, kte
            do i = its, ite
              ext_cof(i,k,nv) = extt   (i,k,1,nv)
              sscal  (i,k,nv) = ssca   (i,k,1,nv)
              asymp  (i,k,nv) = asympar(i,k,1,nv)
            end do
          end do
        end do
        aod2d(its:ite) = aod(its:ite,1)
      end if

    abem(its:ite,7)=aod2d(its:ite)
    abem(:,8)=aerodp(:,1,1)
    abem(:,9)=aerodp(:,1,2)
    abem(:,10)=aerodp(:,1,3)
    abem(:,11)=aerodp(:,1,4)
    abem(:,12)=aerodp(:,1,5)
    endif
!>---- feedback to radiation
!    if (cplchp_rad_opt) then
!     do nv = 1, nbands
!      do k = kts, kte
!       do i = its, ite
!        faersw_cpl(i,k,nv,1) =  ext_cof(i,k,nv)
!        faersw_cpl(i,k,nv,2) =  sscal  (i,k,nv)
!        faersw_cpl(i,k,nv,3) =  asymp  (i,k,nv)
!       end do
!      end do
!     end do
!    endif

!
   end subroutine catchem_rad_wrapper_run
!> @}
   subroutine catchem_prep_rad(                                       &
        ktau,dtstep,ph3d,tk3d,prl3d,spechum,rri,dz8w,                  &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                               &
        ntbc1,ntbc2,ntoc1,ntoc2,                                       &
        ntss1,ntss2,ntss3,ntss4,ntss5,                                 &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                &
        ntrac,gq0,num_chem,ppm2ugkg,chem,relhum,                       &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)

    integer, intent(in) :: ktau
    real(kind=kind_phys), intent(in) :: dtstep

    !FV3 input variables
    integer, intent(in) :: ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
    integer, intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer, intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) :: tk3d,prl3d,spechum
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !GSD Chem variables
    integer,intent(in) ::  num_chem
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg
    
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: rri, dz8w, relhum
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

    ! -- local variables
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_moist) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: z_at_w, p_phy, rho_phy, t_phy
    integer i,ip,j,jp,k,kp,kk,kkp,nv,l,ll,n

    ! -- initialize output arrays
    rri            = 0._kind_phys
    dz8w           = 0._kind_phys
    relhum         = 0._kind_phys
    chem           = 0._kind_phys


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
          relhum(i,k,j) = .99
          relhum(i,k,j) = MIN( .99, moist(i,k,j,1) / &
            (3.80*exp(17.27*(t_phy(i,k,j)-273.)/ &
            (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
          relhum(i,k,j)=max(0.1,relhum(i,k,j))
          !--
        enddo
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


  end subroutine catchem_prep_rad
!> @}
  end module catchem_rad_wrapper
