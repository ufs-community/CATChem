!>\file catchem_gocart_wrapper.F90
!! This file is GSDChem gocart wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 06/2020
!! Revision History:
!! 05/2023, Restructure for CATChem, Jian.He@noaa.gov

module catchem_gocart_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config
   use gocart_aerosols_mod
   use gocart_chem_mod
   use dust_data_mod

   implicit none

   private

   public :: catchem_gocart_wrapper_init, catchem_gocart_wrapper_run, catchem_gocart_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
   subroutine catchem_gocart_wrapper_init()
   end subroutine catchem_gocart_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_gocart_wrapper_finalize Argument Table
!!
   subroutine catchem_gocart_wrapper_finalize()
   end subroutine catchem_gocart_wrapper_finalize

!> \defgroup catchem_group CATChem gocart wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_gocart_wrapper CATChem gocart wrapper Module
!> \ingroup catchem_gocart_group
!! This is the CATChem gocart wrapper Module
!! \section arg_table_catchem_gocart_wrapper_run Argument Table
!! \htmlinclude catchem_gocart_wrapper_run.html
!!
!>\section catchem_gocart_wrapper CATChem Scheme General Algorithm
!> @{
   subroutine catchem_gocart_wrapper_run(im, kte, kme, ktau, dt, garea,   &
      rlat, rlon, julian, xcosz,                               &
      pr3d, ph3d, prl3d, tk3d, spechum, idat, emi2_in, ntrac,  &
      ntso2, ntsulf, ntDMS, ntmsa, ntpp25,                     &
      ntbc1, ntbc2, ntoc1, ntoc2, ntpp10,                      &
      chem_in_opt,chem_opt_in,                                 &
      gq0, qgrs, tile_num, errmsg, errflg)

      implicit none


      integer,        intent(in) :: im,kte,kme,ktau,idat(8),tile_num
      integer,        intent(in) :: ntrac
      integer,        intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
      integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
      real(kind_phys),intent(in) :: dt,julian

      integer, parameter :: ids=1,jds=1,jde=1, kds=1
      integer, parameter :: ims=1,jms=1,jme=1, kms=1
      integer, parameter :: its=1,jts=1,jte=1, kts=1

      real(kind_phys), dimension(im,64, 3), intent(in) :: emi2_in
      real(kind_phys), dimension(im), intent(in) :: garea, rlat,rlon, xcosz
      real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
      real(kind_phys), dimension(im,kte), intent(in) :: prl3d, tk3d, spechum
      real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
      integer,           intent(in) :: chem_in_opt
      integer,           intent(in) :: chem_opt_in
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy,        &
         p_phy, dz8w, p8w, t8w, rho_phy

      real(kind_phys), dimension(ims:im, jms:jme) :: xlat, xlong, dxy

!>- chemistry variables
      real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist
      real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem

      integer :: ide, ime, ite, kde, julday

!   integer, parameter :: SEAS_OPT_DEFAULT = 1
!   integer, parameter :: chem_in_opt = 0  ! 0 for coldstart, 1 for restart
      logical, parameter :: readrestart = .false.
      integer, parameter :: nvl_gocart  = 64  ! number of input levels from gocart file

      real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: pm10, pm2_5_dry, pm2_5_dry_ec

!>- chemical background variables
      real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: backg_oh,backg_h2o2,backg_no3

      real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: oh_t, h2o2_t, no3_t
      real(kind_phys), dimension(ims:im, jms:jme) :: ttday, tcosz

      real(kind_phys) :: dtstep, gmt
      real(kind_phys), dimension(1:num_chem) :: ppm2ugkg

      ! -- output tracers
      real(kind_phys), dimension(ims:im, jms:jme, 1:kme) :: p10, pm25!, ebu_oc
      real(kind_phys), dimension(ims:im, jms:jme, 1:kme) :: oh_bg, h2o2_bg, no3_bg


!>-- local variables
      logical :: call_gocart
      integer :: i, j, jp, k, kp, n


      errmsg = ''
      errflg = 0

      chem_opt          = chem_opt_in

      gmt = real(idat(5))
      julday = real(julian)

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
      call_gocart      = (mod(ktau, call_chemistry) == 0) .or. (ktau == 1)

      ! -- compute accumulated large-scale and convective rainfall since last call
      if (ktau > 1) then
         dtstep = call_chemistry * dt
      else
         dtstep = dt
      end if

!!!

!>- get ready for chemistry run
      call catchem_gocart_prep(                                             &
         readrestart,chem_in_opt,ktau,dtstep,xcosz,                      &
         garea,rlat,rlon,                     &
         pr3d,ph3d,tk3d,prl3d,spechum,                 &
         emi2_in,                                &
         xlat,xlong,dxy,                           &
         rri,t_phy,p_phy,rho_phy,dz8w,p8w,                   &
         t8w,                                               &
         ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                                &
         ntbc1,ntbc2,ntoc1,ntoc2,ntpp10,                                        &
         ntrac,gq0,                                                      &
         num_chem, num_moist,                                 &
         call_gocart,nvl_gocart,                                         &
         ttday,tcosz,gmt,julday,                                         &
         backg_oh,backg_h2o2,backg_no3,                                  &
         ppm2ugkg,                                              &
         moist,chem,                                   &
         ids,ide, jds,jde, kds,kde,                                      &
         ims,ime, jms,jme, kms,kme,                                      &
         its,ite, jts,jte, kts,kte)
      !write (*,*) 'hli test2 ktau',call_gocart

      if (call_gocart) then
         call gocart_chem_driver(ktau,dt,dtstep,gmt,julday,xcosz,          &
            t_phy,moist,chem,rho_phy,dz8w,p8w,backg_oh,oh_t,             &
            backg_h2o2,h2o2_t,backg_no3,no3_t,                           &
            dxy,g,xlat,xlong,ttday,tcosz,                             &
            chem_opt,num_chem,num_moist,                                 &
            ids,ide, jds,jde, kds,kde,                                   &
            ims,ime, jms,jme, kms,kme,                                   &
            its,ite, jts,jte, kts,kte                        )
         call gocart_aerosols_driver(ktau,dtstep,t_phy,moist,              &
            chem,rho_phy,dz8w,p8w,dxy,g,                              &
            chem_opt,num_chem,num_moist,                                 &
            ids,ide, jds,jde, kds,kde,                                   &
            ims,ime, jms,jme, kms,kme,                                   &
            its,ite, jts,jte, kts,kte                        )
      endif

      call sum_pm_gocart (                                                &
         rri, chem,pm2_5_dry, pm2_5_dry_ec, pm10,                       &
         num_chem,chem_opt,                                             &
         ids,ide, jds,jde, kds,kde,                                     &
         ims,ime, jms,jme, kms,kme,                                     &
         its,ite, jts,jte, kts,kte)

      ! -- pm25 and pm10 for output , not for tracer options
      do j = jts, jte
         do k = kts, kte
            do i = its, ite
               pm25  (i,j,k) = pm2_5_dry(i,k,j)
               p10   (i,j,k) = pm10     (i,k,j)
               !ebu_oc(i,j,k) = ebu      (i,k,j,p_ebu_oc)
            end do
         end do
      end do

      if (call_gocart) then
         do j = jts, jte
            do k = kts, kte
               do i = its, ite
                  oh_bg  (i,j,k) = max(0., oh_t  (i,k,j))
                  h2o2_bg(i,j,k) = max(0., h2o2_t(i,k,j))
                  no3_bg (i,j,k) = max(0., no3_t (i,k,j))
               end do
            end do
         end do
      end if


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
            qgrs(i,k,ntpp10 )=gq0(i,k,ntpp10 )
         enddo
      enddo

!
   end subroutine catchem_gocart_wrapper_run
!> @}

   subroutine catchem_gocart_prep(                                        &
      readrestart,chem_in_opt,ktau,dtstep,xcosz,                     &
      garea,rlat,rlon,                     &
      pr3d,ph3d,tk3d,prl3d,spechum,                &
      emi2_in,                               &
      xlat,xlong,dxy,                          &
      rri,t_phy,p_phy,rho_phy,dz8w,p8w,                  &
      t8w,                                              &
      ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                               &
      ntbc1,ntbc2,ntoc1,ntoc2,ntpp10,                                       &
      ntrac,gq0,                                                     &
      num_chem, num_moist,                                &
      call_gocart,nvl_gocart,                                        &
      ttday,tcosz,gmt,julday,                                        &
      backg_oh,backg_h2o2,backg_no3,                                 &
      ppm2ugkg,                                             &
      moist,chem,                                   &
      ids,ide, jds,jde, kds,kde,                                     &
      ims,ime, jms,jme, kms,kme,                                     &
      its,ite, jts,jte, kts,kte)

      !Chem input configuration
      logical, intent(in) :: readrestart
      integer, intent(in) :: chem_in_opt, ktau, julday
      real(kind=kind_phys), intent(in) :: dtstep, gmt

      !FV3 input variables
      integer, intent(in) :: ntrac
      integer, intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
      integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
      real(kind=kind_phys), dimension(ims:ime), intent(in) :: garea, rlat, rlon, xcosz
      real(kind=kind_phys), dimension(ims:ime, 64, 3),   intent(in) :: emi2_in
      real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d,ph3d
      real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         tk3d,prl3d,spechum
      real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


      !GSD Chem variables
      integer,intent(in) ::  num_chem, num_moist,                 &
         nvl_gocart
      logical,intent(in) ::  call_gocart
      integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
         ims,ime, jms,jme, kms,kme,                      &
         its,ite, jts,jte, kts,kte

      real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg
      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme),    intent(out) ::          &
         backg_oh,backg_h2o2,backg_no3


      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              &
         rri, t_phy, p_phy, rho_phy, dz8w, p8w, t8w
      real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) ::              &
         xlat, xlong, dxy,       &
         ttday, tcosz
      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: z_at_w
      real(kind_phys), dimension(nvl_gocart) :: p_gocart

      ! -- local variables
      real(kind_phys), dimension(ims:ime, jms:jme, nvl_gocart) :: oh_backgd,h2o2_backgd,no3_backgd
      real(kind_phys) ::  pu,pl,aln,pwant
      real(kind_phys) ::  xhour,xmin,gmtp,xlonn,xtime,real_time
      real(kind_phys), DIMENSION (1,1) :: sza,cosszax
      integer i,ip,j,jp,k,kp,kk,kkp,nv,jmax,jmaxi,l,ll,n,ndystep,ixhour

      p_gocart = (/ 1000., 992.5, 985., 977.5, 970., 955., 940., 925., 910.,               &
         895., 880., 865., 850., 825., 800., 775., 750., 712.5,  675., 637.5, 600.,     &
         562.5, 525., 487.5, 450., 412.5, 375., 337.5, 288.08, 244.88, 208.15, 176.93,  &
         150.39, 127.84, 108.66, 92.37, 78.51, 66.6, 56.39, 47.64, 40.18, 33.81, 28.37, &
         23.73, 19.79,  16.46, 13.64, 11.28, 9.29, 7.62, 6.22, 5.05, 4.08, 3.28, 2.62,  &
         2.08, 1.65, 1.3, 1.02, 0.8, 0.62, 0.48, 0.37, 0.28 /)

      ! -- initialize output arrays
      backg_oh       = 0._kind_phys
      backg_h2o2     = 0._kind_phys
      backg_no3      = 0._kind_phys
      rri            = 0._kind_phys
      t_phy          = 0._kind_phys
      p_phy          = 0._kind_phys
      rho_phy        = 0._kind_phys
      dz8w           = 0._kind_phys
      p8w            = 0._kind_phys
      t8w            = 0._kind_phys
      xlat           = 0._kind_phys
      xlong          = 0._kind_phys
      dxy            = 0._kind_phys
      ttday          = 0._kind_phys
      tcosz          = 0._kind_phys
      moist          = 0._kind_phys
      chem           = 0._kind_phys
      z_at_w         = 0._kind_phys


      do i=its,ite
         dxy  (i,1)=garea(i)
         xlat (i,1)=rlat(i)*180./pi
         xlong(i,1)=rlon(i)*180./pi
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
            chem(i,k,jts,p_p10   )=max(epsilc,gq0(i,k,ntpp10 )/ppm2ugkg(p_p10))
         enddo
      enddo

      if (.NOT. readrestart) then
         if (chem_in_opt == 0 ) then
            if(ktau.le.1)then
!           if(chem_opt > 0 ) then
               do j=jts,jte
                  jp = j - jts + 1
                  do k=kts,kte
                     do i=its,ite
                        ip = i - its + 1
                        if (chem_opt == CHEM_OPT_GOCART) then
                           do n=1,num_chem
                              chem(i,k,j,n)=1.e-30
                           enddo
                        endif  ! chem_opt==300
                        if ((chem_opt > CHEM_OPT_GOCART) .and. (chem_opt < CHEM_OPT_MAX)) then
                           chem(i,k,j,p_so2)=5.e-10
                           chem(i,k,j,p_sulf)=3.e-10
                           chem(i,k,j,p_msa)=1.e-10
                           chem(i,k,j,p_dms)=1.e-10
                        endif !chem_opt >= 300 .and. chem_opt <  500

!                if ((chem_opt == CHEM_OPT_GOCART_RACM) .or. (chem_opt == CHEM_OPT_RACM_SOA_VBS)) then  !added o3 background !lzhang
!                  kk=min(k,kte)
!                  kkp = kk - kts + 1
!                  ! -- add initial constant into O3,CH4 and CO ect.
!                  chem(i,k,j,p_o3)=epsilc
!                  ! -- this section needs to be revisited before enabling the
!                  ! corresponding chem_opt options
!                  ! maxth=min(400.,th_pvsrf(i,j))
!                  ! if (tr3d(ip,jp,kkp,p_atm_ptem) > maxth) then
!                  !   chem(i,k,j,p_o3)=(airmw/48.)*tr3d(ip,jp,kkp,p_atm_o3mr)*1e6
!                  !   !convert kg/kg to ppm
!                  ! else
!                  !   chem(i,k,j,p_o3)=0.03 !ppm
!                  ! endif
!                  chem(i,k,j,p_ch4)=1.85 !ppm
!                  chem(i,k,j,p_co)=0.06 !ppm
!                  chem(i,k,j,p_co2)=380.
!                  chem(i,k,j,p_ete)=epsilc
!                  chem(i,k,j,p_udd)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_hket)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_api)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_lim)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_dien)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_macr)=chem(i,k,j,p_ete)
!                endif !( (chem_opt == 301.or.chem_opt==108))
                     enddo
                  enddo
               enddo
            endif !(ktau<=1)

         else !(chem_in_opt == 0 )

            if ((ktau<=1).and.((chem_opt == CHEM_OPT_GOCART_RACM).or.(chem_opt == CHEM_OPT_RACM_SOA_VBS))) then  !added GFS o3 background above 380K!lzhang
               do j=jts,jte
                  jp = j - jts + 1
                  do k=kts,kte+1
                     kk=min(k,kte)
                     kkp = kk - kts + 1
                     do i=its,ite
                        ip = i - its + 1
                        ! -- this section needs to be revisited before enabling the
                        ! corresponding chem_opt options
                        ! maxth=min(400.,th_pvsrf(i,j))
                        ! if (tr3d(ip,jp,kkp,p_atm_ptem) >= maxth) then
                        !   chem(i,k,j,p_o3)=(airmw/48.)*tr3d(ip,jp,kkp,p_atm_o3mr)*1e6 !convert kg/kg to ppm
                        ! endif !380K
                     enddo
                  enddo
               enddo
            endif ! chem_opt == 301.or.chem_opt==108

         endif !(chem_in_opt == 1 )
      endif ! readrestart

      !-- assgin read in 3D background chemical species
      do i=its,ite
         do k=1,nvl_gocart
            h2o2_backgd(i,1,k)=emi2_in(i,k,1)
            no3_backgd (i,1,k)=emi2_in(i,k,2)
            oh_backgd  (i,1,k)=emi2_in(i,k,3)
         enddo
      enddo

      !
      ! -- gocart background fields only if gocart is called
      !
      !if (.NOT. readrestart) then
      if (call_gocart .and. (chem_opt == CHEM_OPT_GOCART))then
         do j=jts,jte
            do i=its,ite
               do k=kts,kte
                  do ll=2,nvl_gocart
                     l=ll
                     if (p_gocart(l) < .01*p_phy(i,k,j)) exit
                  enddo
                  pu=alog(p_gocart(l))
                  pl=alog(p_gocart(l-1))
                  pwant=alog(.01*p_phy(i,k,j))
                  if (pwant > pl)then
                     backg_oh(i,k,j)=oh_backgd(i,j,l)
                     backg_h2o2(i,k,j)=h2o2_backgd(i,j,l)
                     backg_no3(i,k,j)=no3_backgd(i,j,l)
                  else
                     aln=(oh_backgd(i,j,l)*(pwant-pl)+            &
                        oh_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
                     backg_oh(i,k,j)=aln
                     aln=(h2o2_backgd(i,j,l)*(pwant-pl)+            &
                        h2o2_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
                     backg_h2o2(i,k,j)=aln
                     aln=(no3_backgd(i,j,l)*(pwant-pl)+            &
                        no3_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
                     backg_no3(i,k,j)=aln
                  endif
               enddo
            enddo
         enddo
      endif   ! end gocart stuff
      !endif !restart


      if ((chem_opt == CHEM_OPT_RACM_SOA_VBS) .or. (chem_opt >= CHEM_OPT_GOCART .and. chem_opt < CHEM_OPT_MAX)) then
         !ndystep=86400/ifix(dtstepc)
         ndystep=86400/ifix(dtstep)
         do j=jts,jte
            do i=its,ite
               tcosz(i,j)=0.
               ttday(i,j)=0.
!         rlat=xlat(i,j)*3.1415926535590/180.
               xlonn=xlong(i,j)
               do n=1,ndystep
                  xtime=n*dtstep/60.
                  ixhour=ifix(gmt+.01)+ifix(xtime/60.)
                  xhour=float(ixhour)
                  xmin=60.*gmt+(xtime-xhour*60.)
                  gmtp=mod(xhour,24.)
                  gmtp=gmtp+xmin/60.
                  CALL szangle(1, 1, julday, gmtp, sza, cosszax,xlonn,rlat(i))
                  TCOSZ(i,j)=TCOSZ(I,J)+cosszax(1,1)
                  if (cosszax(1,1) > 0.) ttday(i,j)=ttday(i,j)+dtstep
                  !--use physics inst cosine zenith -- hli 03/06/2020
!            TCOSZ(i,j)=TCOSZ(I,J)+xcosz(i)
!            if (xcosz(i) > 0.) ttday(i,j)=ttday(i,j)+dtstep
               enddo
            enddo
         enddo
      endif !chem_opt >= 300 .and. chem_opt <  500


   end subroutine catchem_gocart_prep

!> @}
end module catchem_gocart_wrapper
