module gocart_diag_mod

   use catchem_constants ,       only : kind_chem
   use catchem_config

   implicit none

   private

   public :: gocart_diag_cmass
   public :: gocart_diag_store

contains

   subroutine gocart_diag_cmass(chem_opt, nbegin, g, pr, tr, trcm)

      integer,                                intent(in)  :: chem_opt
      integer,                                intent(in)  :: nbegin
      real,                                   intent(in)  :: g
      real(kind_chem), dimension(:,:,:),   intent(in)  :: pr
      real(kind_chem), dimension(:,:,:,:), intent(in)  :: tr
      real(kind_chem), dimension(:,:,:),   intent(out) :: trcm

      ! -- local variables
      integer            :: i, j, k, ni, nj, nk
      real(kind_chem) :: coef

      real(kind_chem), parameter :: fdust2 = 0.38_kind_chem
      real(kind_chem), parameter :: fseas2 = 0.83_kind_chem

      ! -- begin

      ni = size(pr,1)
      nj = size(pr,2)
      nk = size(pr,3) - 1

      trcm = 0._kind_chem

      select case (chem_opt)
       case (CHEM_OPT_GOCART, CHEM_OPT_GOCART_RACM)

         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  coef = 1.e-9_kind_chem * (pr(i,j,k)-pr(i,j,k+1)) / g  !kg/m2
                  ! -- black carbon
                  trcm(i,j,2) = trcm(i,j,2) + coef * (tr(i,j,k,nbegin + p_bc1) + tr(i,j,k,nbegin + p_bc2))
                  ! -- organic carbon
                  trcm(i,j,3) = trcm(i,j,3) + coef * (tr(i,j,k,nbegin + p_oc1) + tr(i,j,k,nbegin + p_oc2))
                  ! -- sulfate
                  trcm(i,j,4) = trcm(i,j,4) + coef * tr(i,j,k,nbegin + p_sulf)
                  ! -- dust
                  trcm(i,j,5) = trcm(i,j,5) + coef * (tr(i,j,k,nbegin + p_dust_1) + fdust2 * tr(i,j,k,nbegin + p_dust_2))
                  ! -- seas
                  trcm(i,j,6) = trcm(i,j,6) + coef * (tr(i,j,k,nbegin + p_seas_1) + tr(i,j,k,nbegin + p_seas_2 ) &
                     + fseas2* tr(i,j,k,nbegin + p_seas_3))
               end do
            end do
         end do
         ! -- pm2.5 aerosol includes all tracers above (note: p25 emissions are added to oc1)
         do k = 2, 6
            do j = 1, nj
               do i = 1, ni
                  trcm(i,j,1) = trcm(i,j,1) + trcm(i,j,k)
               end do
            end do
         end do

       case (CHEM_OPT_RACM_SOA_VBS)

!        do k = 1, nk
!          do j = 1, nj
!            do i = 1, ni
!              coef = 1.e-9_kind_chem * (pr(i,j,k)-pr(i,j,k+1)) / g !kg/m2
!              ! -- pm2.5 aerosol
!              trcm(i,j,1) = trcm(i,j,1) + coef * (tr(i,j,k,nbegin + p_p25i) + tr(i,j,k,nbegin + p_p25j))
!              ! -- black carbon
!              trcm(i,j,2) = trcm(i,j,2) + coef * (tr(i,j,k,nbegin + p_eci) + tr(i,j,k,nbegin + p_ecj))
!              ! -- organic carbon
!              trcm(i,j,3) = trcm(i,j,3) + coef * (tr(i,j,k,nbegin + p_orgpai) + tr(i,j,k,nbegin + p_orgpaj))
!              ! -- sulfate
!              trcm(i,j,4) = trcm(i,j,4) + coef * (tr(i,j,k,nbegin + p_so4ai) + tr(i,j,k,nbegin + p_so4aj))
!              ! -- dust
!              trcm(i,j,5) = trcm(i,j,5) + coef * tr(i,j,k,nbegin + p_soila)
!            end do
!          end do
!        end do
!        ! -- pm2.5 aerosol includes all tracers above
!        do k = 2, 5
!          do j = 1, nj
!            do i = 1, ni
!              trcm(i,j,1) = trcm(i,j,1) + trcm(i,j,k)
!            end do
!          end do
!        end do

       case default
         ! -- not yet implemented
      end select

   end subroutine gocart_diag_cmass


   subroutine gocart_diag_store(ipos, v, w)

      integer,                                intent(in)    :: ipos
      real(kind_chem), dimension(:,:,:),   intent(in)    :: v
      real(kind_chem), dimension(:,:,:,:), intent(inout) :: w
      real(kind_chem), parameter :: ugkg = 1.e-09_kind_chem !lzhang

      ! -- local variables
      integer :: m, n, nd, nt

      nd = size(w, dim=4)
      if (ipos > nd) return

      nt = size(v, dim=3)
      if (nt > size(w, dim=3) + 2) return

      m = 0
      do n = 1, nt
         if (n == p_so2) cycle
         if (n == p_msa) cycle
         m = m + 1
         w(:,:,m,ipos) = ugkg*v(:,:,n) !kg/m2/s
      end do

   end subroutine gocart_diag_store

end module gocart_diag_mod
