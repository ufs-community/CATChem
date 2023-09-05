!>\file catchem_diag_wrapper.F90
!! This file is GSDChem diagnostic wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 09/2020
!! Revision History:
!! 05/2023, Restructure for CATChem, Jian.He@noaa.gov

module catchem_diag_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config
   use gocart_diag_mod

   implicit none

   private

   public :: catchem_diag_wrapper_init, catchem_diag_wrapper_run, catchem_diag_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
   subroutine catchem_diag_wrapper_init()
   end subroutine catchem_diag_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_diag_wrapper_finalize Argument Table
!!
   subroutine catchem_diag_wrapper_finalize()
   end subroutine catchem_diag_wrapper_finalize

!> \defgroup catchem_group CATChem diag wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_diag_wrapper CATChem diag wrapper Module
!> \ingroup catchem_diag_group
!! This is the CATChem diag wrapper Module
!! \section arg_table_catchem_diag_wrapper_run Argument Table
!! \htmlinclude catchem_diag_wrapper_run.html
!!
!>\section catchem_diag_wrapper CATChem Scheme General Algorithm
!> @{
   subroutine catchem_diag_wrapper_run(im, kte, kme, ktau,                  &
      pr3d, ntrac, ntso2, gq0, aecm, ntchmdiag, ntchm,           &
      wetdpc, wetdpc_deep, wetdpc_mid,  wetdpc_shal,             &
      imfdeepcnv, imfdeepcnv_samf, imfdeepcnv_gf, chem_opt_in,   &
      errmsg,errflg)

      implicit none


      integer,        intent(in) :: im,kte,kme,ktau
      integer,        intent(in) :: ntrac,ntso2,ntchmdiag,ntchm
      integer,        intent(in)  :: imfdeepcnv, imfdeepcnv_samf, imfdeepcnv_gf

      integer, parameter :: ids=1,jds=1,jde=1, kds=1
      integer, parameter :: ims=1,jms=1,jme=1, kms=1
      integer, parameter :: its=1,jts=1,jte=1, kts=1

      real(kind_phys), dimension(im,kme), intent(in) :: pr3d
      real(kind_phys), dimension(im,ntchmdiag), intent(inout) :: wetdpc
      real(kind_phys), dimension(im,ntchm), intent(in) :: wetdpc_deep, wetdpc_mid, wetdpc_shal
      real(kind_phys), dimension(im,kte,ntrac), intent(in) :: gq0
      real(kind_phys), dimension(im,6        ), intent(inout) :: aecm
      integer,        intent(in) :: chem_opt_in
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! -- for diagnostics
      real(kind_phys), dimension(ims:im, jms:jme, 6) :: trcm  ! inst tracer column mass density
      real(kind_phys), dimension(ims:im, jms:jme, ntchmdiag, 4) :: trdf
      real(kind_phys), dimension(im,jme,kte,ntrac) :: gq0j
      real(kind_phys), dimension(im,jme,kme) :: pr3dj
      real(kind_phys), dimension(ims:im, jms:jme, 1:ntchm) :: wet_dep

      integer :: ide, ime, ite, kde

      integer :: i, j, jp, k, kp, n, nbegin


      errmsg = ''
      errflg = 0

      chem_opt          = chem_opt_in

      ! -- set domain
      ide=im
      ime=im
      ite=im
      kde=kte

      nbegin = ntso2-1
      pr3dj(:,1,:  )=pr3d(:,:  )
      gq0j (:,1,:,:)=gq0 (:,:,:)
      ! -- calculate column mass density
      call gocart_diag_cmass(chem_opt, nbegin, g, pr3dj, gq0j, trcm)
      aecm(:,:)=trcm(:,1,:)

      ! -- calculate convective wet deposition
      if (imfdeepcnv == imfdeepcnv_samf) then
         do n=1,ntchm
            do i=1,im
               wet_dep(i,1,n) = (max(0.,wetdpc_deep(i,n)) +max(0.,wetdpc_shal(i,n)))
            enddo
         enddo
      elseif (imfdeepcnv == imfdeepcnv_gf) then
         do n=1,ntchm
            do i=1,im
               wet_dep(i,1,n) = (max(0.,wetdpc_deep(i,n))+max(0.,wetdpc_mid(i,n)) +max(0.,wetdpc_shal(i,n)))
            enddo
         enddo
      endif

      call gocart_diag_store(4, wet_dep, trdf)
      wetdpc (:,:)=trdf(:,1,:,4)

!
   end subroutine catchem_diag_wrapper_run
!> @}

end module catchem_diag_wrapper
