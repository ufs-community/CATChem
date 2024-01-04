!>\file catchem_settling_wrapper.F90
!! This file is GSDChem settling wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 06/2020
!! Revision History:
!! 05/2023, Restructure for CATChem, Jian.He@noaa.gov

module catchem_settling_wrapper

   use physcons, g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config
   use gocart_settling_mod
   use vash_settling_mod
   use gocart_diag_mod

   implicit none

   private

   public :: catchem_settling_wrapper_init, catchem_settling_wrapper_run, catchem_settling_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
   subroutine catchem_settling_wrapper_init()
   end subroutine catchem_settling_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_settling_wrapper_finalize Argument Table
!!
   subroutine catchem_settling_wrapper_finalize()
   end subroutine catchem_settling_wrapper_finalize

!> \defgroup catchem_group CATChem settling wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_settling_wrapper CATChem settling wrapper Module
!> \ingroup catchem_settling_group
!! This is the CATChem settling wrapper Module
!! \section arg_table_catchem_settling_wrapper_run Argument Table
!! \htmlinclude catchem_settling_wrapper_run.html
!!
!>\section catchem_settling_wrapper CATChem Scheme General Algorithm
!> @{
   subroutine catchem_settling_wrapper_run(im, kte, kme, ktau, dt, garea,      &
      pr3d, ph3d, prl3d, tk3d, spechum,                             &
      ntrac,ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                        &
      ntbc1,ntbc2,ntoc1,ntoc2,ntss1,ntss2,ntss3,ntss4,ntss5,        &
      ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,ntchmdiag,     &
      gq0,qgrs,sedimio, dust_opt_in, seas_opt_in,                   &
      errmsg,errflg)

      implicit none


      integer,        intent(in) :: im,kte,kme,ktau
      integer,        intent(in) :: ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
      integer,        intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
      integer,        intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10,ntchmdiag
      integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
      real(kind_phys),intent(in) :: dt

      integer, parameter :: ids=1,jds=1,jde=1, kds=1
      integer, parameter :: ims=1,jms=1,jme=1, kms=1
      integer, parameter :: its=1,jts=1,jte=1, kts=1

      real(kind_phys), dimension(im), intent(in) :: garea
      real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
      real(kind_phys), dimension(im,kte), intent(in) :: prl3d, tk3d, spechum
      real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
      real(kind_phys), dimension(im,ntchmdiag), intent(inout) :: sedimio
      integer,        intent(in) :: dust_opt_in, seas_opt_in
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: t_phy, p_phy, dz8w, p8w, rho_phy

      real(kind_phys), dimension(ims:im, jms:jme) :: dxy

!>- sea salt & chemistry variables
      real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist
      real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
      real(kind_phys), dimension(ims:im, jms:jme, 1:num_chem )  :: sedim
      real(kind_phys), dimension(ims:im, jms:jme) :: seashelp, dusthelp

      integer :: ide, ime, ite, kde

      real(kind_phys), dimension(ims:im, jms:jme) :: ash_fall
      real(kind_phys), dimension(im, 1, ntchmdiag, 4) :: trdf

      real(kind_phys), dimension(1:num_chem) :: ppm2ugkg

!>-- local variables
      integer :: i, j, jp, k, kp, n


      errmsg = ''
      errflg = 0

      dust_opt          = dust_opt_in
      seas_opt          = seas_opt_in

      ash_fall   = 0.
      trdf       = 0.

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
      call catchem_prep_settling(ktau,garea,                             &
         pr3d,ph3d,tk3d,prl3d,spechum,dxy,                               &
         t_phy,p_phy,rho_phy,dz8w,p8w,                                   &
         ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                                &
         ntbc1,ntbc2,ntoc1,ntoc2,                                        &
         ntss1,ntss2,ntss3,ntss4,ntss5,                                  &
         ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                 &
         ntrac,gq0,num_chem, num_moist,ppm2ugkg,moist,chem,              &
         ids,ide, jds,jde, kds,kde,                                      &
         ims,ime, jms,jme, kms,kme,                                      &
         its,ite, jts,jte, kts,kte)

      if ((dust_opt /= DUST_OPT_NONE) .or.                                &
         (seas_opt /= SEAS_OPT_NONE)) then

         select case (chem_opt)
          case (304, 316, 317)
            !
            !  GOCART "very" light
            !
            do j=jts,jte
               do i=its,ite
                  call settling_simple_driver(dt,t_phy(i,:,j),              &
                     moist(i,:,j,:),chem(i,:,j,:),rho_phy(i,:,j),            &
                     dz8w(i,:,j),p8w(i,:,j),p_phy(i,:,j),                  &
                     sedim(i,j,:),dusthelp(i,j),seashelp(i,j),             &
                     dxy(i,j),kms,kme,kts,kte)
               enddo
            enddo

          case default
            !
            ! run with all GOCART variables, GOCART sort of HEAVY!
            !
            do j=jts,jte
               do i=its,ite
                  call settling_gocart_driver(dt,t_phy(i,:,j),              &
                     moist(i,:,j,:),chem(i,:,j,:),rho_phy(i,:,j),            &
                     dz8w(i,:,j),p8w(i,:,j),p_phy(i,:,j),                  &
                     sedim(i,j,:),dxy(i,j),kms,kme,kts,kte)
               enddo
            enddo

         end select
      end if

      ! -- 4 volcanic size bins
      do j=jts,jte
         do i=its,ite
            call vashshort_settling_driver(dt,t_phy(i,:,j),              &
               moist(i,:,j,:),chem(i,:,j,:),rho_phy(i,:,j),            &
               dz8w(i,:,j),p8w(i,:,j),p_phy(i,:,j),                  &
               dxy(i,j),ash_fall(i,j),kms,kme,kts,kte)
         enddo
      enddo

      ! -- put chem stuff back into tracer array
      do k=kts,kte
         do i=its,ite
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

      ! -- output sedimentation
      call gocart_diag_store(1, sedim, trdf)

      sedimio(:,:)=trdf(:,1,:,1)

!
   end subroutine catchem_settling_wrapper_run
!> @}
   subroutine catchem_prep_settling(ktau,garea,                       &
      pr3d,ph3d,tk3d,prl3d,spechum,dxy,                              &
      t_phy,p_phy,rho_phy,dz8w,p8w,                                  &
      ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                               &
      ntbc1,ntbc2,ntoc1,ntoc2,                                       &
      ntss1,ntss2,ntss3,ntss4,ntss5,                                 &
      ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                &
      ntrac,gq0,num_chem, num_moist,ppm2ugkg,moist,chem,             &
      ids,ide, jds,jde, kds,kde,                                     &
      ims,ime, jms,jme, kms,kme,                                     &
      its,ite, jts,jte, kts,kte)

      !Chem input configuration
      integer, intent(in) :: ktau

      !FV3 input variables
      integer, intent(in) :: ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
      integer, intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
      integer, intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
      integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
      real(kind=kind_phys), dimension(ims:ime), intent(in) :: garea
      real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d,ph3d
      real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) :: tk3d,prl3d,spechum
      real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


      !GSD Chem variables
      integer,intent(in) ::  num_chem, num_moist
      integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
         ims,ime, jms,jme, kms,kme,                      &
         its,ite, jts,jte, kts,kte

      real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg

      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: t_phy, p_phy, rho_phy, dz8w, p8w
      real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) :: dxy
      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

      ! -- local variables
      integer i,ip,j,jp,k,kp,kk,kkp,l,ll,n
      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: z_at_w

      ! -- initialize output arrays
      t_phy          = 0._kind_phys
      p_phy          = 0._kind_phys
      rho_phy        = 0._kind_phys
      dz8w           = 0._kind_phys
      p8w            = 0._kind_phys
      dxy            = 0._kind_phys
      moist          = 0._kind_phys
      chem           = 0._kind_phys

      do i=its,ite
         dxy  (i,1)=garea(i)
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

      do k=kms,kte
         do i=ims,ime
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
         enddo
      enddo



   end subroutine catchem_prep_settling
!> @}
end module catchem_settling_wrapper
