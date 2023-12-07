!>\file catchem_plume_wrapper.F90
!! This file is GSDChem plume wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 07/2020
!! Kate.Zhang@noaa.gov 10/2023
!! Revision History:
!! 05/2023, Restructure for CATChem, Jian.He@noaa.gov

 module catchem_plume_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config
   use plume_data_mod
   use plume_rise_mod

   implicit none

   private

   public :: catchem_plume_wrapper_init, catchem_plume_wrapper_run, catchem_plume_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine catchem_plume_wrapper_init()
      end subroutine catchem_plume_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_plume_wrapper_finalize Argument Table
!!
      subroutine catchem_plume_wrapper_finalize()
      end subroutine catchem_plume_wrapper_finalize

!> \defgroup catchem_plume_group CATChem plume wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_plume_wrapper CATChem plume wrapper Module  
!> \ingroup catchem_plume_group
!! This is the CATChem plume wrapper Module
!! \section arg_table_catchem_plume_wrapper_run Argument Table
!! \htmlinclude catchem_plume_wrapper_run.html
!!
!>\section catchem_plume_wrapper CATChem Scheme General Algorithm
!> @{
    subroutine catchem_plume_wrapper_run(im, kte, kme, ktau, dt, julian,         &
                   pr3d, ph3d,phl3d, prl3d, tk3d, us3d, vs3d, spechum,           &
                   w,vegtype,fire_GBBEPx,fire_MODIS,                             &
                   ntrac,ntso2,ntpp25,ntbc1,ntoc1,ntpp10,igb,                    &
                   gq0,qgrs,ebu,abem,biomass_burn_opt_in,plumerise_flag_in,      &
                   plumerisefire_frq_in,pert_scale_plume,ca_emis_plume,          &
                   do_ca,ca_sgs,ca_sgs_emis,vegtype_cpl,ca_sgs_gbbepx_frp,       &
                   emis_amp_plume, do_sppt_emis, sppt_wts, errmsg,errflg)

    implicit none


    integer,        intent(in) :: im,kte,kme,ktau,igb
    integer,        intent(in) :: ntrac,ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    real(kind_phys),intent(in) :: dt, emis_amp_plume, pert_scale_plume,julian

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    logical,        intent(in) :: do_sppt_emis, do_ca, ca_sgs_emis, ca_sgs
    real(kind_phys), intent(in) :: sppt_wts(:,:), ca_emis_plume(:)
    integer, dimension(im), intent(in) :: vegtype    
    integer, dimension(im), intent(out) :: vegtype_cpl
    real(kind_phys), dimension(im,    5, igb), intent(in) :: fire_GBBEPx
    real(kind_phys), dimension(im,   13), intent(in) :: fire_MODIS
    real(kind_phys), intent(out) :: ca_sgs_gbbepx_frp(:)
    real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte), intent(in) :: phl3d, prl3d, tk3d,        &
                us3d, vs3d, spechum, w
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
    real(kind_phys), dimension(im,12        ), intent(inout) :: abem
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_ebu), intent(inout) :: ebu
    integer,        intent(in) :: biomass_burn_opt_in, plumerise_flag_in, plumerisefire_frq_in
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy, u_phy, v_phy,       &
                     p_phy, z_at_w, dz8w, p8w, rho_phy, vvel

!>- sea salt & chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist 
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem

    integer :: ide, ime, ite, kde
    integer,         dimension(ims:im, jms:jme) :: ivgtyp

!>- plume variables
    ! -- buffers
    real(kind_phys), dimension(ims:im, jms:jme, num_ebu_in) :: ebu_in
    real(kind_phys), dimension(ims:im, jms:jme) ::                              &
         mean_fct_agef, mean_fct_aggr, mean_fct_agsv, mean_fct_agtf,            &
         firesize_agef, firesize_aggr, firesize_agsv, firesize_agtf,            &
         ca_sgs_gbbepx_frp_with_j
    real(kind_phys), dimension(ims:im, jms:jme, num_frp_plume ) :: plume_frp
    real(kind_phys) :: dtstep
   !integer,parameter :: plumerise_flag = 2  ! 1=MODIS, 2=GBBEPx
    logical :: call_plume, scale_fire_emiss, doing_sgs_emis
    logical, save :: firstfire = .true.
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg
    real(kind_phys), parameter :: ugkg = 1.e-09_kind_phys !lzhang

!>-- local variables
    real(kind_phys) :: curr_secs
    real(kind_phys) :: factor, factor2, factor3, random_factor(ims:im)
    integer :: nbegin
    integer :: i, j, jp, k, kp, n
  
    errmsg = ''
    errflg = 0

    biomass_burn_opt  = biomass_burn_opt_in
    plumerise_flag    = plumerise_flag_in
    plumerisefire_frq = plumerisefire_frq_in
    random_factor = 1.0
    curr_secs = ktau * dt
    doing_sgs_emis = do_ca .and. ca_sgs_emis .and. .not. ca_sgs

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
    call_plume       = (biomass_burn_opt == BURN_OPT_ENABLE) .and. (plumerisefire_frq > 0)
    if (call_plume) &
       call_plume    = (mod(int(curr_secs), max(1, 60*plumerisefire_frq)) == 0)      &
                        .or. (ktau == 1)
    scale_fire_emiss = .false.

    ! -- compute accumulated large-scale and convective rainfall since last call
    if (ktau > 1) then
      dtstep = call_chemistry * dt
    else
      dtstep = dt
    end if

!>- get ready for chemistry run
    call catchem_prep_plume(ktau,dtstep,julian,                         &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,                 &
        vegtype,fire_GBBEPx,fire_MODIS,                                 &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                   &
        z_at_w,vvel,                                                    &
        ntso2,ntpp25,ntbc1,ntoc1,ntpp10,igb,ntrac,gq0,                  &
        num_chem, num_moist,num_ebu_in,ca_sgs_gbbepx_frp_with_j,        &
        plumerise_flag,num_plume_data,ppm2ugkg,                         &
        mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,        &
        firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr,        &
        moist,chem,plume_frp,ebu_in,ivgtyp,ca_emis_plume,doing_sgs_emis,&
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)

    ! Input to cellular automata
    if(doing_sgs_emis) then
      !do i=ids,ide
      do i=its,ite
        ca_sgs_gbbepx_frp(i) = ca_sgs_gbbepx_frp_with_j(i,jds)
        vegtype_cpl(i) = vegtype(i)
      enddo
    endif


    ! compute wild-fire plumes
    if (call_plume) then
      call plumerise_driver (ktau,dtstep,num_chem,num_ebu,num_ebu_in,   &
        ebu,ebu_in,                                                     &
        mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,        &
        firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr,        &
        'GOCART','BIOMASSB', t_phy,moist(:,:,:,p_qv),                   &
        rho_phy,vvel,u_phy,v_phy,p_phy,                                 &
        z_at_w,scale_fire_emiss,plume_frp,plumerise_flag,               &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte                                )
    end if
    ! -- add biomass burning emissions at every timestep
    if (biomass_burn_opt == BURN_OPT_ENABLE) then
      jp = jte

      if(plumerise_flag == FIRE_OPT_GBBEPx .and. do_sppt_emis) then
        random_factor(:) = pert_scale_plume*max(min(1+(sppt_wts(:,kme/2)-1)*emis_amp_plume,2.0),0.0)
      endif

      factor3 = 0._kind_phys
      select case (plumerise_flag)
        case (FIRE_OPT_MODIS)
          factor3 = 4.828e-04_kind_phys/60.
          kp = kte    ! full column
        case (FIRE_OPT_GBBEPx)
          factor3 = 1.e-03_kind_phys * mwdry / mw_so2_aer
          if (plumerisefire_frq > 0) then
            kp = kte  ! full column
          else
            kp = kts  ! surface only
          end if
        case default
          ! -- no further options available, skip this step
          jp = jts - 1
      end select

      if (kp == kts) then
        ! -- only include surface emissions
        k = kts
        do j = jts, jp
          do i = its, ite
            ! -- factor for pm emissions, factor2 for burn emissions
            factor  = dt*rri(i,k,j)/dz8w(i,k,j)*random_factor(i)
            factor2 = factor * factor3
            chem(i,k,j,p_oc1) = chem(i,k,j,p_oc1) + factor  * ebu_in(i,j,p_ebu_in_oc  )
            chem(i,k,j,p_bc1) = chem(i,k,j,p_bc1) + factor  * ebu_in(i,j,p_ebu_in_bc  )
            chem(i,k,j,p_p25) = chem(i,k,j,p_p25) + factor  * ebu_in(i,j,p_ebu_in_pm25)
            chem(i,k,j,p_p10) = chem(i,k,j,p_p10) + factor  * ebu_in(i,j,p_ebu_in_pm10)
            chem(i,k,j,p_so2) = chem(i,k,j,p_so2) + factor2 * ebu_in(i,j,p_ebu_in_so2 )
          end do
        end do

      else
        ! -- use full-column emissions
        do j = jts, jp
          do k = kts, kp
            do i = its, ite
              ! -- factor for pm emissions, factor2 for burn emissions
              factor  = dt*rri(i,k,j)/dz8w(i,k,j)*random_factor(i)
              factor2 = factor * factor3
              chem(i,k,j,p_oc1) = chem(i,k,j,p_oc1) + factor  * ebu(i,k,j,p_ebu_oc  )
              chem(i,k,j,p_bc1) = chem(i,k,j,p_bc1) + factor  * ebu(i,k,j,p_ebu_bc  )
              chem(i,k,j,p_p25) = chem(i,k,j,p_p25) + factor  * ebu(i,k,j,p_ebu_pm25)
              chem(i,k,j,p_p10) = chem(i,k,j,p_p10) + factor  * ebu(i,k,j,p_ebu_pm10)
              chem(i,k,j,p_so2) = chem(i,k,j,p_so2) + factor2 * ebu(i,k,j,p_ebu_so2 )
            end do
          end do
        end do
      end if

    end if


    ! -- put chem stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       gq0(i,k,ntso2  )=ppm2ugkg(p_so2   ) * max(epsilc,chem(i,k,1,p_so2))
       gq0(i,k,ntpp25 )=ppm2ugkg(p_p25   ) * max(epsilc,chem(i,k,1,p_p25))
       gq0(i,k,ntbc1  )=ppm2ugkg(p_bc1   ) * max(epsilc,chem(i,k,1,p_bc1))
       gq0(i,k,ntoc1  )=ppm2ugkg(p_oc1   ) * max(epsilc,chem(i,k,1,p_oc1))
       gq0(i,k,ntpp10 )=ppm2ugkg(p_p10   ) * max(epsilc,chem(i,k,1,p_p10))
     enddo
    enddo

    do k=kts,kte
     do i=its,ite
       qgrs(i,k,ntso2 )=gq0(i,k,ntso2  )
       qgrs(i,k,ntpp25)=gq0(i,k,ntpp25 )
       qgrs(i,k,ntbc1 )=gq0(i,k,ntbc1  )
       qgrs(i,k,ntoc1 )=gq0(i,k,ntoc1  )
       qgrs(i,k,ntpp10)=gq0(i,k,ntpp10 )
     enddo
    enddo

    abem(:,4)=ugkg*ebu_in  (:,kts,p_ebu_in_bc )*random_factor(:)
    abem(:,5)=ugkg*ebu_in  (:,kts,p_ebu_in_oc )*random_factor(:)
    abem(:,6)=ugkg*ebu_in  (:,kts,p_ebu_in_so2)*random_factor(:)


   end subroutine catchem_plume_wrapper_run
!> @}

   subroutine catchem_prep_plume(                                        &
        ktau,dtstep, julian,                    &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,                &
        vegtype,                  &
        fire_GBBEPx,fire_MODIS,                              &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                  &
        z_at_w,vvel,                                              &
        ntso2,ntpp25,                               &
        ntbc1,ntoc1,                                       &
        ntpp10,igb,                &
        ntrac,gq0,                                                     &
        num_chem, num_moist,num_ebu_in,ca_sgs_gbbepx_frp_with_j,       &
        plumerise_flag,num_plume_data,                    &
        ppm2ugkg,                                             &
        mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,       &
        firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr,       &
        moist,chem,plumedist,ebu_in,                                   &
        ivgtyp,ca_emis_plume,doing_sgs_emis,              &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)

    !Chem input configuration
    integer, intent(in) :: ktau
    real(kind=kind_phys), intent(in) :: dtstep,julian

    !FV3 input variables
    integer, dimension(ims:ime), intent(in) :: vegtype
    integer, intent(in) :: ntrac,igb
    integer, intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    logical, intent(in) :: doing_sgs_emis
    real(kind=kind_phys), intent(in) :: ca_emis_plume(:)
    real(kind=kind_phys), dimension(ims:ime,     5,igb),   intent(in) :: fire_GBBEPx
    real(kind=kind_phys), dimension(ims:ime,    13),   intent(in) :: fire_MODIS
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) ::     &
         pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,us3d,vs3d,spechum,w
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !GSD Chem variables
    integer,intent(in) ::  num_chem, num_moist, num_ebu_in,                &
                           plumerise_flag, num_plume_data
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg
    real(kind_phys), dimension(:, :), intent(out) :: ca_sgs_gbbepx_frp_with_j
    real(kind_phys), dimension(ims:ime, jms:jme, num_ebu_in),intent(out) :: ebu_in
    
    integer,dimension(ims:ime, jms:jme), intent(out) :: ivgtyp
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, u_phy, v_phy, p_phy, rho_phy, dz8w, p8w, vvel
         
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w
    real(kind_phys), dimension(ims:ime, jms:jme, num_frp_plume), intent(out) :: plumedist
    real(kind_phys), dimension(ims:ime, jms:jme   ), intent(out) ::                    &
                   mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,            &
                   firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr       
    real(kind_phys), dimension(ims:ime, jms:jme, num_ebu_in) :: emiss_abu
    real(kind_phys), dimension(ims:ime, jms:jme, num_plume_data) :: plume
    real(kind_phys), parameter :: frp2plume = 1.e+06_kind_phys  ! FRP-to-plume conversion factor
    real(kind_phys), parameter :: frpc  = 1.e+09_kind_phys      ! FRP conversion factor

    ! -- local variables
!   real(kind=kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: p_phy
    real(kind_phys) ::  factor,factor2
    integer i,ip,j,jp,k,kp,kk,kkp,l,ll,n,ii
    integer, save :: curr_day
    ! -- initialize output arrays
    ebu_in         = 0._kind_phys
    ivgtyp         = 0._kind_phys
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    u_phy          = 0._kind_phys
    v_phy          = 0._kind_phys
    p_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    vvel           = 0._kind_phys
    moist          = 0._kind_phys  
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys

    ! -- initialize fire emissions
    plume          = 0._kind_phys
    plumedist      = 0._kind_phys
    mean_fct_agtf  = 0._kind_phys
    mean_fct_agef  = 0._kind_phys
    mean_fct_agsv  = 0._kind_phys
    mean_fct_aggr  = 0._kind_phys
    firesize_agtf  = 0._kind_phys
    firesize_agef  = 0._kind_phys
    firesize_agsv  = 0._kind_phys
    firesize_aggr  = 0._kind_phys


    do i=its,ite
     ivgtyp (i,1)=vegtype(i)
    enddo
   

    if (ktau <= 1) then
    curr_day = int(julian)
    end if
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
          vvel(i,k,j)=-w(ip,kkp)*rri(i,k,j)/g 
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
        enddo
      enddo
    enddo

    ! -- fire
    emiss_abu = 0.   ! fire emission

    !print*,'hli ',plumerise_flag,FIRE_OPT_MODIS,FIRE_OPT_GBBEPx
    select case (plumerise_flag)
      case (FIRE_OPT_MODIS)
        do j=jts,jte
         do i=its,ite
          emiss_abu(i,j,p_e_bc)   =fire_MODIS(i,1)
          emiss_abu(i,j,p_e_oc)   =fire_MODIS(i,2)
          emiss_abu(i,j,p_e_pm_25)=fire_MODIS(i,3)
          emiss_abu(i,j,p_e_so2)  =fire_MODIS(i,4)
          emiss_abu(i,j,p_e_pm_10)=fire_MODIS(i,5)
          plume(i,j,1)            =fire_MODIS(i,6)
          plume(i,j,2)            =fire_MODIS(i,7)
          plume(i,j,3)            =fire_MODIS(i,8)
          plume(i,j,4)            =fire_MODIS(i,9)
          plume(i,j,5)            =fire_MODIS(i,10)
          plume(i,j,6)            =fire_MODIS(i,11)
          plume(i,j,7)            =fire_MODIS(i,12)
          plume(i,j,8)            =fire_MODIS(i,13)
         enddo
        enddo
      case (FIRE_OPT_GBBEPx)
          if (igb == 1) then
            ii=1
          else
            ii=int(julian)-curr_day+1
          endif
        do j=jts,jte
         do i=its,ite
          emiss_abu(i,j,p_e_bc)   =fire_GBBEPx(i,1,ii)
          emiss_abu(i,j,p_e_oc)   =fire_GBBEPx(i,2,ii)
          emiss_abu(i,j,p_e_pm_25)=fire_GBBEPx(i,3,ii)
          emiss_abu(i,j,p_e_so2)  =fire_GBBEPx(i,4,ii)
         enddo
        enddo
        if(doing_sgs_emis) then
          do j=jts,jte
           do i=its,ite
             ca_sgs_gbbepx_frp_with_j(i,j) = fire_GBBEPx(i,5,ii)
             plume(i,j,1)            =ca_emis_plume(i)! *0.5 + fire_GBBEPx(i,5)*0.5
           enddo
          enddo
        else
          do j=jts,jte
           do i=its,ite
             plume(i,j,1)            =fire_GBBEPx(i,5,ii)
           enddo
          enddo
        endif
!        print*,'hli GBBEPx plume',maxval(plume(:,:,1))
      case default
          ! -- no further option available
    end select


    factor=0.
    k=kts
    if (p_bc2 > 1) then
      do j=jts,jte
        do i=its,ite


          ebu_in(i,j,p_ebu_in_pm10)=emiss_abu(i,j,p_e_pm_10)
          ebu_in(i,j,p_ebu_in_dms)= 0._kind_phys

          select case (plumerise_flag)
            case (FIRE_OPT_MODIS)
              ebu_in(i,j,p_ebu_in_oc)   = emiss_abu(i,j,p_e_oc)
              ebu_in(i,j,p_ebu_in_bc)   = emiss_abu(i,j,p_e_bc)
              ebu_in(i,j,p_ebu_in_pm25) = emiss_abu(i,j,p_e_pm_25)
              ebu_in(i,j,p_ebu_in_so2)  = emiss_abu(i,j,p_e_so2)
              mean_fct_agtf(i,j)=plume(i,j,1)
              mean_fct_agef(i,j)=plume(i,j,2)
              mean_fct_agsv(i,j)=plume(i,j,3)
              mean_fct_aggr(i,j)=plume(i,j,4)
              firesize_agtf(i,j)=plume(i,j,5)
              firesize_agef(i,j)=plume(i,j,6)
              firesize_agsv(i,j)=plume(i,j,7)
              firesize_aggr(i,j)=plume(i,j,8)
            case (FIRE_OPT_GBBEPx)
              ebu_in(i,j,p_ebu_in_oc)   = frpc * emiss_abu(i,j,p_e_oc)
              ebu_in(i,j,p_ebu_in_bc)   = frpc * emiss_abu(i,j,p_e_bc)
              ebu_in(i,j,p_ebu_in_pm25) = frpc * (emiss_abu(i,j,p_e_pm_25) - emiss_abu(i,j,p_e_bc) - emiss_abu(i,j,p_e_oc))
              ebu_in(i,j,p_ebu_in_so2)  = frpc * emiss_abu(i,j,p_e_so2)
              plumedist(i,j,p_frp_flam_frac) = flaming(catb(ivgtyp(i,j)))
              plumedist(i,j,p_frp_mean     ) = frp2plume * plume(i,j,1)
              plumedist(i,j,p_frp_std      ) = 0.3_kind_phys   * frp2plume * plume(i,j,1)
              plumedist(i,j,p_frp_mean_size) = msize(ivgtyp(i,j)) * frp2plume * plume(i,j,1)
              plumedist(i,j,p_frp_std_size ) = 0.5_kind_phys * plumedist(i,j,p_frp_mean_size)
            case default
              ! -- no further option available
          end select
        enddo
      enddo
    endif
 
    do k=kms,kte
     do i=ims,ime
       chem(i,k,jts,p_so2   )=max(epsilc,gq0(i,k,ntso2  )/ppm2ugkg(p_so2))
       chem(i,k,jts,p_p25   )=max(epsilc,gq0(i,k,ntpp25 )/ppm2ugkg(p_p25))
       chem(i,k,jts,p_bc1   )=max(epsilc,gq0(i,k,ntbc1  )/ppm2ugkg(p_bc1))
       chem(i,k,jts,p_oc1   )=max(epsilc,gq0(i,k,ntoc1  )/ppm2ugkg(p_oc1))
       chem(i,k,jts,p_p10   )=max(epsilc,gq0(i,k,ntpp10 )/ppm2ugkg(p_p10))
     enddo
    enddo


  end subroutine catchem_prep_plume

!> @}
  end module catchem_plume_wrapper
