!11/2023, Optimized and Updated by Kate.Zhang@noaa.gov for 2 threads job
module plume_rise_mod

  use catchem_constants, only : kind_chem,g => con_g, cp => con_cp, &
                                 r_d => con_rd, r_v => con_rv 

  use catchem_config

  use plume_data_mod,  only : nveg_agreg
  use plume_zero_mod
  use plume_scalar_mod

  implicit none

  real(kind=kind_chem),parameter :: p1000mb = 100000.  ! p at 1000mb (pascals)

  private

  public :: num_frp_plume
  public :: plumerise_driver

contains

  subroutine plumerise_driver (ktau,dtstep,num_chem,num_ebu,num_ebu_in,          &
             ebu,ebu_in,                                                         &
             mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,            &
             firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr,            &
             chem_opt,burn_opt,t_phy,q_vap,                                      &
             rho_phy,vvel,u_phy,v_phy,p_phy,                                     &
             z_at_w,scale_fire_emiss,plume_frp,plumerise_flag,                   &
             ids,ide, jds,jde, kds,kde,                                          &
             ims,ime, jms,jme, kms,kme,                                          &
             its,ite, jts,jte, kts,kte                                         )

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ktau,num_chem,num_ebu,               &
                                  num_ebu_in,plumerise_flag,           &
                                  ids,ide, jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte
   REAL(kind=kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme, num_ebu ),              &
         INTENT(INOUT ) ::                                   ebu
   REAL(kind=kind_chem), DIMENSION( ims:ime, jms:jme, num_ebu_in ),                    &
         INTENT(INOUT ) ::                                   ebu_in
   REAL(kind=kind_chem), DIMENSION( ims:ime, jms:jme, num_frp_plume ),                 &
         INTENT(INOUT ) ::                                   plume_frp
   REAL(kind=kind_chem), DIMENSION( ims:ime, jms:jme ),INTENT(IN ) ::     &
           mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,       &
           firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr

!
!
!
   REAL(kind=kind_chem), DIMENSION( ims:ime , kms:kme , jms:jme ),                     &
         INTENT(IN   ) ::                                              &
                                                      t_phy,           &
                 z_at_w,vvel,u_phy,v_phy,rho_phy,p_phy,q_vap
   REAL(kind=kind_chem), INTENT(IN   ) ::                             dtstep

   LOGICAL,      INTENT(IN   ) :: scale_fire_emiss
   character (len=*), intent(in) :: chem_opt,burn_opt

!
! Local variables...
!
      INTEGER :: nv, i, j, k, ksub, nspecies

      character(len=100)  :: errmsg
      integer           :: errflg
!     integer, parameter :: nspecies=num_ebu
      real(kind=kind_chem), dimension (num_ebu) :: eburn_in 
      real(kind=kind_chem), dimension (kte,num_ebu) :: eburn_out
      real(kind=kind_chem), dimension (kte) :: u_in ,v_in ,w_in ,theta_in ,pi_in  &
                              ,rho_phyin ,qv_in ,zmid    &
                              ,z_lev
      real(kind=kind_chem), dimension(nveg_agreg) :: firesize,mean_fct
      real(kind=kind_chem) :: sum, ffirs, rcp,ratio,zk,cpor
!     real(kind=kind_chem),save,dimension(its:ite,jts:jte) :: ffirs
      ffirs=0.
      rcp=r_d/cp
      cpor=cp/r_d

      nspecies=num_ebu

      if( scale_fire_emiss ) then
        if( chem_opt /= 'MOZCART_KPP' .and. &
            chem_opt /= 'MOZART_KPP' .and. &
            chem_opt /= 'MOZART_MOSAIC_4BIN_VBS0_KPP' ) then
          print*,"Fire emission scaling only supported for MOZART_KPP, MOZCART_KPP chem chem_opts"
          !call chem_comm_abort(msg="Fire emission scaling only supported for MOZART_KPP, MOZCART_KPP chem chem_opts")
        endif
      endif

       if ( burn_opt == 'BIOMASSB' ) then
         do j=jts,jte
            do i=its,ite
            if ( chem_opt == 'GOCARTRACM'.or.chem_opt == 'RACMSOAVBS' ) then
!               ebu(i,kts,j,p_ebu_no)=ebu_in(i,j,p_ebu_in_no)
!               ebu(i,kts,j,p_ebu_no2)=ebu_in(i,j,p_ebu_in_no2)
!               ebu(i,kts,j,p_ebu_co)=ebu_in(i,j,p_ebu_in_co)
!               ebu(i,kts,j,p_ebu_co2)=ebu_in(i,j,p_ebu_in_co2)
!               ebu(i,kts,j,p_ebu_eth)=ebu_in(i,j,p_ebu_in_eth)
!               ebu(i,kts,j,p_ebu_hc3)=ebu_in(i,j,p_ebu_in_hc3)
!               ebu(i,kts,j,p_ebu_hc5)=ebu_in(i,j,p_ebu_in_hc5)
!               ebu(i,kts,j,p_ebu_hc8)=ebu_in(i,j,p_ebu_in_hc8)
!              ! ebu(i,kts,j,p_ebu_ete)=ebu_in(i,j,p_ebu_in_ete)
!               ebu(i,kts,j,p_ebu_olt)=ebu_in(i,j,p_ebu_in_olt)
!               ebu(i,kts,j,p_ebu_oli)=ebu_in(i,j,p_ebu_in_oli)
!              ! ebu(i,kts,j,p_ebu_dien)=ebu_in(i,j,p_ebu_in_dien)
!               ebu(i,kts,j,p_ebu_iso)=ebu_in(i,j,p_ebu_in_iso)
!              ! ebu(i,kts,j,p_ebu_api)=ebu_in(i,j,p_ebu_in_api)
!              ! ebu(i,kts,j,p_ebu_lim)=ebu_in(i,j,p_ebu_in_lim)
!               ebu(i,kts,j,p_ebu_tol)=ebu_in(i,j,p_ebu_in_tol)
!               ebu(i,kts,j,p_ebu_xyl)=ebu_in(i,j,p_ebu_in_xyl)
!               ebu(i,kts,j,p_ebu_csl)=ebu_in(i,j,p_ebu_in_csl)
!               ebu(i,kts,j,p_ebu_hcho)=ebu_in(i,j,p_ebu_in_hcho)
!               ebu(i,kts,j,p_ebu_ald)=ebu_in(i,j,p_ebu_in_ald)
!               ebu(i,kts,j,p_ebu_ket)=ebu_in(i,j,p_ebu_in_ket)
!              ! ebu(i,kts,j,p_ebu_macr)=ebu_in(i,j,p_ebu_in_macr)
!              !ebu(i,kts,j,p_ebu_ora1)=ebu_in(i,j,p_ebu_in_ora1)
!               ebu(i,kts,j,p_ebu_ora2)=ebu_in(i,j,p_ebu_in_ora2)
               endif
               ebu(i,kts,j,p_ebu_dms)=ebu_in(i,j,p_ebu_in_dms)
!               ebu(i,kts,j,p_ebu_sulf)=ebu_in(i,j,p_ebu_in_sulf)
               ebu(i,kts,j,p_ebu_bc)=ebu_in(i,j,p_ebu_in_bc)
               ebu(i,kts,j,p_ebu_oc)=ebu_in(i,j,p_ebu_in_oc)
               ebu(i,kts,j,p_ebu_so2)=ebu_in(i,j,p_ebu_in_so2)
               ebu(i,kts,j,p_ebu_pm25)=ebu_in(i,j,p_ebu_in_pm25)
               ebu(i,kts,j,p_ebu_pm10)=ebu_in(i,j,p_ebu_in_pm10)
            enddo
         enddo
       elseif ( burn_opt == 'BIOMASSB_MOZC' .or. &
                burn_opt == 'BIOMASSB_MOZ' ) then
         do j=jts,jte
            do i=its,ite
!               ebu(i,kts,j,p_ebu_no)=ebu_in(i,j,p_ebu_in_no)
!               ebu(i,kts,j,p_ebu_co)=ebu_in(i,j,p_ebu_in_co)
!               ebu(i,kts,j,p_ebu_bigalk) = ebu_in(i,j,p_ebu_in_bigalk)
!               ebu(i,kts,j,p_ebu_bigene) = ebu_in(i,j,p_ebu_in_bigene)
!               ebu(i,kts,j,p_ebu_c2h4) = ebu_in(i,j,p_ebu_in_c2h4)
!               ebu(i,kts,j,p_ebu_c2h5oh) = ebu_in(i,j,p_ebu_in_c2h5oh)
!               ebu(i,kts,j,p_ebu_c2h6) = ebu_in(i,j,p_ebu_in_c2h6)
!               ebu(i,kts,j,p_ebu_c3h6) = ebu_in(i,j,p_ebu_in_c3h6)
!               ebu(i,kts,j,p_ebu_c3h8) = ebu_in(i,j,p_ebu_in_c3h8)
!               ebu(i,kts,j,p_ebu_ch2o) = ebu_in(i,j,p_ebu_in_ch2o)
!               ebu(i,kts,j,p_ebu_ch3cho) = ebu_in(i,j,p_ebu_in_ch3cho)
!               ebu(i,kts,j,p_ebu_ch3coch3) = ebu_in(i,j,p_ebu_in_ch3coch3)
!               ebu(i,kts,j,p_ebu_ch3oh) = ebu_in(i,j,p_ebu_in_ch3oh)
!               ebu(i,kts,j,p_ebu_mek) = ebu_in(i,j,p_ebu_in_mek)
!               ebu(i,kts,j,p_ebu_so2) = ebu_in(i,j,p_ebu_in_so2)
!               ebu(i,kts,j,p_ebu_toluene) = ebu_in(i,j,p_ebu_in_toluene)
!               ebu(i,kts,j,p_ebu_nh3) = ebu_in(i,j,p_ebu_in_nh3)
!               ebu(i,kts,j,p_ebu_no2)=ebu_in(i,j,p_ebu_in_no2)
!               ebu(i,kts,j,p_ebu_open) = ebu_in(i,j,p_ebu_in_open)
!               ebu(i,kts,j,p_ebu_c10h16) = ebu_in(i,j,p_ebu_in_c10h16)
!               ebu(i,kts,j,p_ebu_mgly) = ebu_in(i,j,p_ebu_in_mgly)
!               ebu(i,kts,j,p_ebu_ch3cooh) = ebu_in(i,j,p_ebu_in_ch3cooh)
!               ebu(i,kts,j,p_ebu_cres) = ebu_in(i,j,p_ebu_in_cres)
!               ebu(i,kts,j,p_ebu_glyald) = ebu_in(i,j,p_ebu_in_glyald)
!               ebu(i,kts,j,p_ebu_gly)=ebu_in(i,j,p_ebu_in_gly)
!               ebu(i,kts,j,p_ebu_acetol) = ebu_in(i,j,p_ebu_in_acetol)
!               ebu(i,kts,j,p_ebu_isop) = ebu_in(i,j,p_ebu_in_isop)
!               ebu(i,kts,j,p_ebu_macr) = ebu_in(i,j,p_ebu_in_macr)
!               ebu(i,kts,j,p_ebu_mvk)=ebu_in(i,j,p_ebu_in_mvk)
            enddo
         enddo
         if( burn_opt == 'BIOMASSB_MOZC' ) then
           do j=jts,jte
!             ebu(its:ite,kts,j,p_ebu_pm10) = ebu_in(its:ite,j,p_ebu_in_pm10)
!             ebu(its:ite,kts,j,p_ebu_pm25) = ebu_in(its:ite,j,p_ebu_in_pm25)
!             ebu(its:ite,kts,j,p_ebu_oc) = ebu_in(its:ite,j,p_ebu_in_oc)
!             ebu(its:ite,kts,j,p_ebu_bc) = ebu_in(its:ite,j,p_ebu_in_bc)
           enddo
         endif
       elseif ( burn_opt == 'BIOMASSB_GHG' ) then
         do j=jts,jte
            do i=its,ite
!               ebu(i,kts,j,p_ebu_co)  = ebu_in(i,j,p_ebu_in_co)
!               ebu(i,kts,j,p_ebu_co2) = ebu_in(i,j,p_ebu_in_co2)
!               ebu(i,kts,j,p_ebu_ch4) = ebu_in(i,j,p_ebu_in_ch4)
            enddo
          enddo
       endif
!
       do nv=1,num_ebu
          do j=jts,jte
            do k=kts+1,kte
               do i=its,ite
                 ebu(i,k,j,nv)=0.
               enddo
            enddo
          enddo
       enddo
       
       do j=jts,jte
          do i=its,ite
            select case (plumerise_flag)
              case (FIRE_OPT_MODIS)
                sum=mean_fct_agtf(i,j)+mean_fct_agef(i,j)+mean_fct_agsv(i,j)    &
                        +mean_fct_aggr(i,j)
                if(sum.lt.1.e-6)Cycle
    !           ffirs=ffirs+1
                sum=firesize_agtf(i,j)+firesize_agef(i,j)+firesize_agsv(i,j)    &
                        +firesize_aggr(i,j)
                if(sum.lt.1.e-6)Cycle
                eburn_out=0.
                mean_fct(1)=mean_fct_agtf(i,j)
                mean_fct(2)=mean_fct_agef(i,j)
                mean_fct(3)=mean_fct_agsv(i,j)
                mean_fct(4)=mean_fct_aggr(i,j)
                firesize(1)=firesize_agtf(i,j)
                firesize(2)=firesize_agef(i,j)
                firesize(3)=firesize_agsv(i,j)
                firesize(4)=firesize_aggr(i,j)
              case (FIRE_OPT_GBBEPx)
                if (plume_frp(i,j,p_frp_mean) < 1.e-06) cycle
              case default
                ! -- no further option available
            end select
            do nv=1,num_ebu
            eburn_in(nv)=ebu(i,kts,j,nv)
            enddo
            if( maxval( eburn_in(:) ) == 0. ) cycle
            do k=kts,kte
              u_in(k)=u_phy(i,k,j)
              v_in(k)=v_phy(i,k,j)
              w_in(k)=vvel(i,k,j)
              qv_in(k)=q_vap(i,k,j)
              pi_in(k)=cp*(p_phy(i,k,j)/p1000mb)**rcp
              !zk=.5*(z_at_w(i,k+1,j)-z_at_w(i,k,j))
              zk=.5*(z_at_w(i,k+1,j)+z_at_w(i,k,j)) !lzhang
              zmid(k)=zk-z_at_w(i,kts,j)
              z_lev(k)=z_at_w(i,k,j)-z_at_w(i,kts,j)
              rho_phyin(k)=rho_phy(i,k,j)
              theta_in(k)=t_phy(i,k,j)/pi_in(k)*cp
            enddo
!!$              pi_in(kte)=pi_in(kte-1)  !wig: These are no longer needed after changing definition
!!$              u_in(kte)=u_in(kte-1)    !     of kte in chem_driver (12-Oct-2007)
!!$              v_in(kte)=v_in(kte-1)
!!$              w_in(kte)=w_in(kte-1)
!!$              qv_in(kte)=qv_in(kte-1)
!!$              zmid(kte)=z(i,kte,j)-z_at_w(i,kts,j)
!!$              z_lev(kte)=z_at_w(i,kte,j)-z_at_w(i,kts,j)
!!$              rho_phyin(kte)=rho_phyin(kte-1)
!!$              theta_in(kte)=theta_in(kte-1)
            call plumerise(kte,1,1,1,1,1,1,firesize,mean_fct  &
                    ,nspecies,eburn_in,eburn_out              &
                    ,u_in ,v_in ,w_in ,theta_in ,pi_in        &
                    ,rho_phyin ,qv_in ,zmid                   &
                    ,z_lev,plume_frp(i,j,:),plumerise_flag    &
                    ,g, cp,r_d, cpor,errmsg, errflg )

            do nv=1,num_ebu
              do k=kts+1,kte
                ebu(i,k,j,nv)=eburn_out(k,nv)*(z_at_w(i,k+1,j)-z_at_w(i,k,j))
              enddo
            enddo
!            print*,'hli plumerise',maxval(eburn_out),maxval(plume_frp(i,j,:))

has_total_emissions : &
            if( scale_fire_emiss ) then
is_mozcart : &
              if( (chem_opt == 'MOZCART_KPP' .and. &
                   burn_opt == 'BIOMASSB_MOZC') .or. &
                  (chem_opt == 'MOZART_KPP' .and. &
                   burn_opt == 'BIOMASSB_MOZ') .or. &
                  (chem_opt == 'MOZART_MOSAIC_4BIN_VBS0_KPP' .and. &
                   burn_opt == 'BIOMASSB_MOZC') ) then
!-------------------------------------------------------------------
! we input total emissions instead of smoldering emissions:
! ratio of smolderling to total
!-------------------------------------------------------------------
                sum = 0.
                do k = kts,kte
                  sum = sum + ebu(i,k,j,p_ebu_co)
                end do
                if( sum > 0. ) then             
                  ratio = ebu(i,kts,j,p_ebu_co)/sum
                else
                  ratio = 0.
                endif

                do k = kts,kte
!                  ebu(i,k,j,p_ebu_no) = ebu(i,k,j,p_ebu_no)*ratio
!                  ebu(i,k,j,p_ebu_co) = ebu(i,k,j,p_ebu_co)*ratio
!                  ebu(i,k,j,p_ebu_bigalk) = ebu(i,k,j,p_ebu_bigalk)*ratio
!                  ebu(i,k,j,p_ebu_bigene) = ebu(i,k,j,p_ebu_bigene)*ratio
!                  ebu(i,k,j,p_ebu_c2h4)   = ebu(i,k,j,p_ebu_c2h4)*ratio
!                  ebu(i,k,j,p_ebu_c2h5oh) = ebu(i,k,j,p_ebu_c2h5oh)*ratio
!                  ebu(i,k,j,p_ebu_c2h6) = ebu(i,k,j,p_ebu_c2h6)*ratio
!                  ebu(i,k,j,p_ebu_c3h6) = ebu(i,k,j,p_ebu_c3h6)*ratio
!                  ebu(i,k,j,p_ebu_c3h8) = ebu(i,k,j,p_ebu_c3h8)*ratio
!                  ebu(i,k,j,p_ebu_ch2o) = ebu(i,k,j,p_ebu_ch2o)*ratio
!                  ebu(i,k,j,p_ebu_ch3cho) = ebu(i,k,j,p_ebu_ch3cho)*ratio
!                  ebu(i,k,j,p_ebu_ch3coch3) = ebu(i,k,j,p_ebu_ch3coch3)*ratio
!                  ebu(i,k,j,p_ebu_ch3oh)    = ebu(i,k,j,p_ebu_ch3oh)*ratio
!                  ebu(i,k,j,p_ebu_mek) = ebu(i,k,j,p_ebu_mek)*ratio
!                  ebu(i,k,j,p_ebu_so2) = ebu(i,k,j,p_ebu_so2)*ratio
!                  ebu(i,k,j,p_ebu_toluene) = ebu(i,k,j,p_ebu_toluene)*ratio
!                  ebu(i,k,j,p_ebu_nh3) = ebu(i,k,j,p_ebu_nh3)*ratio
!                  ebu(i,k,j,p_ebu_no2)  = ebu(i,k,j,p_ebu_no2)*ratio
!                  ebu(i,k,j,p_ebu_open) = ebu(i,k,j,p_ebu_open)*ratio
!                  ebu(i,k,j,p_ebu_c10h16) = ebu(i,k,j,p_ebu_c10h16)*ratio
!                  ebu(i,k,j,p_ebu_mgly) = ebu(i,k,j,p_ebu_mgly)*ratio
!                  ebu(i,k,j,p_ebu_ch3cooh) = ebu(i,k,j,p_ebu_ch3cooh)*ratio
!                  ebu(i,k,j,p_ebu_cres) = ebu(i,k,j,p_ebu_cres)*ratio
!                  ebu(i,k,j,p_ebu_glyald) = ebu(i,k,j,p_ebu_glyald)*ratio
!                  ebu(i,k,j,p_ebu_gly) = ebu(i,k,j,p_ebu_gly)*ratio
!                  ebu(i,k,j,p_ebu_acetol) = ebu(i,k,j,p_ebu_acetol)*ratio
!                  ebu(i,k,j,p_ebu_isop) = ebu(i,k,j,p_ebu_isop)*ratio
!                  ebu(i,k,j,p_ebu_macr) = ebu(i,k,j,p_ebu_macr)*ratio
!                  ebu(i,k,j,p_ebu_mvk)  = ebu(i,k,j,p_ebu_mvk)*ratio
                end do
                if( chem_opt == 'MOZCART_KPP' .or. &
                    chem_opt == 'MOZART_MOSAIC_4BIN_VBS0_KPP' ) then
                  do k = kts,kte
!                    ebu(i,k,j,p_ebu_pm10) = ebu(i,k,j,p_ebu_pm10)*ratio
!                    ebu(i,k,j,p_ebu_pm25) = ebu(i,k,j,p_ebu_pm25)*ratio
!                    ebu(i,k,j,p_ebu_oc)  = ebu(i,k,j,p_ebu_oc)*ratio
!                    ebu(i,k,j,p_ebu_bc)  = ebu(i,k,j,p_ebu_bc)*ratio
                  end do
                endif
              end if is_mozcart
            end if has_total_emissions

          enddo
          enddo
  end subroutine plumerise_driver

end module plume_rise_mod
