!! Revision History:
!! 06/2023, Restructure for CATChem, Jian.He@noaa.gov

module gocart_dmsemis_mod

   use catchem_constants ,         only : kind_chem, g=>con_g, pi=>con_pi
   use catchem_config,  only : smw, NDMS, mwdry, num_chem, p_dms

   implicit none

   public :: gocart_dmsemis

contains

   subroutine gocart_dmsemis(dt,u_phy,v_phy,          &
      chem_arr,dz8w,u10,v10,    &
      delp,dms_0,tsk,area)

      IMPLICIT NONE

      REAL(kind=kind_chem), DIMENSION( num_chem ),                 &
         INTENT(INOUT ) ::                                   chem_arr
      REAL(kind_chem), INTENT(IN ) :: dt, u_phy,v_phy, &
         dz8w,u10,v10,  &
         delp,dms_0,tsk, area

!
! local variables
!
      integer :: i,j,k,ndt,imx,jmx,lmx,nmx
      integer,dimension (1,1) :: ilwi
      real(kind_chem), DIMENSION (1,1,1,1) :: tc
      real(kind_chem), DIMENSION (1,1,1) :: bems,airmas
      real(kind_chem), DIMENSION (1,1) :: emsdms
      real(kind_chem), dimension (1,1) :: w10m,gwet,airden,tskin,dmso
      real(kind_chem), dimension (1) :: dxy
      real(kind_chem),parameter::max_default=1.e-30
!
! number of dust bins
!
      imx=1
      jmx=1
      lmx=1
      nmx=1
      ndt=int(dt)

      ilwi(1,1)=0

      tc(1,1,1,1)=chem_arr(p_dms)*1.d-6
      dmso(1,1)=dms_0
! w10m(1,1)=sqrt(u10*u10+v10*v10)
      airmas(1,1,1)=delp*area/g
      dxy(1)=area
      tskin(1,1)=tsk
      emsdms(1,1)=0.
!
! we don't trust the u10,v10 values, is model layers are very thin near surface
!
! if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
      w10m=sqrt(u_phy*u_phy+v_phy*v_phy)

      call  srcdms(imx, jmx, lmx, nmx, ndt, tc, mwdry,&
         tskin, ilwi, dmso, w10m, airmas, dxy, emsdms, bems)
!    chem(i,kts,j,p_dms)=max(1.e-30,tc(1,1,1,1)*1.e6)
      chem_arr(p_dms)=max(max_default,tc(1,1,1,1)*1.e6)
!
   end subroutine gocart_dmsemis

   SUBROUTINE srcdms(imx, jmx, lmx, nmx, ndt1, tc,airmw, &
      tskin, ilwi, dmso, w10m, airmas, dxy, emsdms, bems)

      ! **************************************************************************
      ! **                                                                      **
      ! **  This subroutine calculates DMS emissions from the ocean.            **
      ! **                                                                      **
      ! **************************************************************************


      IMPLICIT NONE
      REAL(kind_chem),    PARAMETER :: dms_mw = 62.00
      REAL(kind_chem),    PARAMETER :: tcmw(1)=dms_mw
      INTEGER, INTENT(IN)    :: imx, jmx, lmx, nmx, ndt1
      REAL(kind_chem),    INTENT(IN)    :: tskin(imx,jmx), dmso(imx,jmx)
      INTEGER, INTENT(IN)    :: ilwi(imx,jmx)
      REAL(kind_chem),    INTENT(IN)    :: dxy(jmx), w10m(imx,jmx)
      REAL(kind_chem),      INTENT(IN)    :: airmw
      REAL(kind_chem),    INTENT(IN)    :: airmas(imx,jmx,lmx)
      REAL(kind_chem),    INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
      REAL(kind_chem),    INTENT(INOUT) :: emsdms(imx,jmx)
      REAL(kind_chem),    INTENT(OUT)   :: bems(imx,jmx,nmx)

      INTEGER :: i,j
      REAL(kind_chem)    :: sst, sc, conc, w10, scco2, akw, erate, dmssrc, c

      ! **************************************************************************
      ! *  ilwi = 0: water                                                      **
      ! *  If ilwi = 0: DMSEMS = seawaterDMS * transfer velocity.               **
      ! *  Otherwise,  DMSEMS = 0.0                                             **
      ! **************************************************************************

      ! executable statements
! tcmw(NDMS) = dms_mw
      lat_loop: DO j = 1,jmx
         lon_loop: DO i = 1,imx
            ! convert tskin (=sst over water) from K to degC
            sst = tskin(i,j) - 273.15
!       if_water: IF (ilwi(i,j) == 0) THEN

            ! -- Schmidt number for DMS (Saltzman et al., 1993)
            sc = 2674.0 - 147.12*sst + 3.726*(sst**2) - 0.038*(sst**3)

! ****************************************************************************
! *  Calculate transfer velocity in cm/hr  (AKw)                             *
! *                                                                          *
! *  Tans et al. transfer velocity (1990) for CO2 at 25oC (Erickson, 1993)   *
! *                                                                          *
! *  Tans et al. assumed AKW=0 when W10<=3. I modified it to let             *
! *  DMS emit at low windseeds too. Chose 3.6m/s as the threshold.           *
! *                                                                          *
! *  Schmidt number for CO2:       Sc = 600  (20oC, fresh water)             *
! *                                Sc = 660  (20oC, seawater)                *
! *                                Sc = 428  (25oC, Erickson 93)             *
! ****************************************************************************

            conc = dmso(i,j)

            w10  = w10m(i,j)
!           ! --- GEOS-1 or GEOS-STRAT: using SSMI winds
!           IF (lmx <= 26) w10 = wssmi(i,j)

! ---  Tans et al. (1990) -----------------
!       ScCO2 = 428.
!       if (W10 .le. 3.6) then
!        AKw = 1.0667 * W10
!       else
!        AKw = 6.4 * (W10 - 3.)
!       end if

! ---  Wanninkhof (1992) ------------------
!       ScCO2 = 660.
!       AKw = 0.31 * W10**2

            ! ---  Liss and Merlivat (1986) -----------
            scco2 = 600.0
            IF (w10 <= 3.6) THEN
               akw = 0.17 * w10
            ELSE IF (w10 <= 13.0) THEN
               akw = 2.85 * w10 - 9.65
            ELSE
               akw = 5.90 * w10 - 49.3
            END IF
            ! ------------------------------------------
            IF (w10 <= 3.6) THEN
               if (sc .le. 0.) then
                  akw=1.0E-8
               else
                  akw = akw * ((scco2/sc) ** 0.667)
               endif
            ELSE
               if (sc .le. 0.) then
                  akw=1.0E-8
               else
                  akw = akw * SQRT(scco2/sc)
               endif
            END IF

! ****************************************************************************
! *  Calculate emission flux in kg/box/timestep                              *
! *                                                                          *
! *   AKw is in cm/hr:                 AKw/100/3600    -> m/sec              *
! *   CONC is in nmol/L (nmol/dm3):    CONC*1E-12*1000 -> kmol/m3            *
! *   TCMW(NDMS)       : kgDMS/kmol                                          *
! *   ERATE            : kgDMS/m2/timestep                                   *
! *   DMSSRC           : kgDMS/box/timestep                                  *
! ****************************************************************************

            erate  = akw/100.0/3600.0*conc*1.0E-12*1000.0*REAL(ndt1)*tcmw(NDMS)
            dmssrc = erate * dxy(j)

!       ELSE   ! ilwi /= 0 (water)

!          dmssrc = 0.0

!       END IF if_water

! ****************************************************************************
! *  Update DMS concentration in level 1 (where emission occurs)             *
! ****************************************************************************

            ! -- Convert emission from kg/box/timestep to mixing ratio/timestep:
            c = dmssrc / airmas(i,j,1) * airmw / tcmw(NDMS)
            tc(i,j,1,NDMS) = tc(i,j,1,NDMS) + c

            !    ---------------------------------------------------------------
            !     Diagnostics:      DMS surface emission in kgS/timestep
            !    ---------------------------------------------------------------
            emsdms(i,j) = emsdms(i,j) + dmssrc * smw / tcmw(NDMS) ! kgS
!        bems(i,j,NDMS) = c * airmas(i,j,1) / airmw * smw ! kgS
            bems(i,j,NDMS) = dmssrc  ! kgDMS

         END DO lon_loop
      END DO lat_loop

   END SUBROUTINE srcdms

end module gocart_dmsemis_mod
