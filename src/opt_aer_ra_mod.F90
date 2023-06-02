MODULE opt_aer_ra_mod
   use catchem_constants ,        only : kind_chem
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: aer_ra
CONTAINS
   SUBROUTINE aer_ra(dz8w  &
                   ,extt,ssca,asympar,nbands                      &
                   ,tauaersw,gaersw,waersw,tauaerlw                &
                   ,ids,ide, jds,jde, kds,kde                     &
                   ,ims,ime, jms,jme, kms,kme                     &
                   ,its,ite, jts,jte, kts,kte )
   IMPLICIT NONE

   INTEGER,    INTENT(IN   ) ::        ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte,nbands
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme,nbands ), INTENT (OUT) :: extt
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme,nbands ), INTENT (OUT) :: ssca
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme, nbands), INTENT (OUT) :: asympar
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme, 4 ),        &
         INTENT(IN    ) :: tauaersw,gaersw,waersw
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme ),           &
         INTENT(IN    ) :: dz8w
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme, 16 ),        &
         INTENT(IN    ) :: tauaerlw
   real(kind_chem) :: ang,slope,slopeg,slopessa,onemang    
   integer :: i,j,k,ib
   real(kind_chem), dimension(NBANDS) :: midbands  ! jcb
   REAL(kind_chem),    PARAMETER ::   thresh=1.e-9
!  ---  band wavenumber intervals
   real(kind_chem), dimension(14):: wvnum1, wvnum2
      data wvnum1/                                        &
     &         2600.0, 3251.0, 4001.0, 4651.0, 5151.0, 6151.0, 7701.0,  &
     &         8051.0,12851.0,16001.0,22651.0,29001.0,38001.0,  820.0 /
      data wvnum2/                                        &
     &         3250.0, 4000.0, 4650.0, 5150.0, 6150.0, 7700.0, 8050.0,  &
     &        12850.0,16000.0,22650.0,29000.0,38000.0,50000.0, 2600.0 /

!  data midbands/.2,.235,.27,.2875,.3025,.305,.3625,.55,1.92,1.745,6.135/


! As in PNNL MOSAIC, extrapolate or interpolate based on 300-999 nm Angstrom coefficient,
! or linear interpolation/extrapolation between 300 and 999 nm for assymetry coefficient
      do ib=1,nbands
      midbands(ib)=(1./wvnum1(ib)+1./wvnum2(ib))*.5e4
      do j = jts,jte
      do k = kts,kte
      do i = its,ite
        if(tauaersw(i,k,j,1).gt.thresh .and. tauaersw(i,k,j,4).gt.thresh) then
           ang=log(tauaersw(i,k,j,1)/tauaersw(i,k,j,4))/log(999./300.)
           extt(i,k,j,ib)=tauaersw(i,k,j,2)*(0.4/midbands(ib))**ang 

! ssa - linear interpolation; extrapolation
           slope=(waersw(i,k,j,3)-waersw(i,k,j,2))/.2
           ssca(i,k,j,ib) = slope*(midbands(ib)-.6)+waersw(i,k,j,3) 
           if(ssca(i,k,j,ib).lt.0.4) ssca(i,k,j,ib)=0.4
           if(ssca(i,k,j,ib).ge.1.0) ssca(i,k,j,ib)=1.0

! g - linear interpolation;extrapolation
           slope=(gaersw(i,k,j,3)-gaersw(i,k,j,2))/.2
           asympar(i,k,j,ib) = slope*(midbands(ib)-.6)+gaersw(i,k,j,3) 
           if(asympar(i,k,j,ib).lt.0.5) asympar(i,k,j,ib)=0.5
           if(asympar(i,k,j,ib).ge.1.0) asympar(i,k,j,ib)=1.0
        else
           extt(i,k,j,ib)=0.
           ssca(i,k,j,ib)=1.
           asympar(i,k,j,ib)=0.
        endif

      end do
      end do
      end do
      end do
!
      do ib=1,nbands
      do j = jts,jte
      do i = its,ite
         slope = 0.  !use slope as a sum holder
         do k = kts,kte
            slope = slope + extt(i,k,j,ib)
         end do
         if( slope < 0. ) then
            write(0,*)'ERROR: Negative total optical depth',j,slope
         else if( slope > 6. ) then
            write(0,*)'adjusting extt ',ib,j,slope
            do k = kts,kte
              extt(i,k,j,ib)=extt(i,k,j,ib)*6./slope
            enddo
         endif
      end do
      end do
      end do

  END SUBROUTINE aer_ra
END MODULE opt_aer_ra_mod

