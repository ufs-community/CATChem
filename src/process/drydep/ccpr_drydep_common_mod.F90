!>
!! \file ccpr_drydep_common_mod.F90
!! \brief Contains module ccpr_drydep_common_mod
!!
!! \ingroup catchem_drydep_process
!!
!! \author Wei Li
!! \date 02/2025
!!!>
module CCPr_drydep_Common_Mod
   use precision_mod, only: fp, ZERO, rae
   use Error_Mod
   use constants
   implicit none
   !private

   public  :: Wesely_Rc_Gas
   public  :: AERO_SFCRSII
   public  :: Wesely_Ra_Rb
   public  :: INIT_WEIGHTSS


   ! module variables (mainly some constants dependent on land use in the scheme)

   integer,  parameter :: NDRYDTYPE   = 11    !< # of drydep land types following GEOS-Chem
   real(fp), parameter :: TWO_THIRDS  = 2.0_fp / 3.0_fp
   !real(fp), parameter :: H2OMW = 18.0_fp !declared in constant module
   real(fp), parameter :: SMALL = 1.0e-10_fp !< Small number
   integer,  parameter :: IWATER = 1      !< Index for water in Olson land use
   ! Arrays that hold information for each of the 11 drydep land types
   integer :: IDRYDTYPE(NDRYDTYPE)
   real(fp) :: IRAC(NDRYDTYPE),  IRCLO(NDRYDTYPE), IRCLS(NDRYDTYPE)
   real(fp) :: IRGSS(NDRYDTYPE), IRGSO(NDRYDTYPE), IRLU(NDRYDTYPE)
   real(fp) :: IRI(NDRYDTYPE),   IVSMAX(NDRYDTYPE)
   !some Olson land use (74 types) related parameters
   real(fp):: DRYCOEFF(20) !< DRYCOEFF : Baldocchi polynomial coeffs
   integer :: IOLSON (74), IDEP(74)
   integer :: IZO(74) !< Roughness height for each Olson land types

   !assign some drydep values to arrays based on the 11 drydep land use in GEOS-Chem.
   !Wesely (1989) is separated into seasons, but not sure how GEOS-Chem gets its values (TODO).
   !You can find the values and references in https://wiki.seas.harvard.edu/geos-chem/index.php/Dry_deposition
   !***********************************************************************
   !* The land types within each grid square are defined using the Olson
   !* land-type database.  Each of the Olson land types is assigned a
   !* corresponding "deposition land type" with characteristic values of
   !* surface resistance components.  There are 74 Olson land-types but only
   !* 11 deposition land-types (i.e., many of the Olson land types share the
   !* same deposition characteristics).  Surface resistance components for
   !* the "deposition land types" are from Wesely [1989] except for tropical
   !* forests [Jacob and Wofsy, 1990] and for tundra [Jacob et al., 1992].
   !* All surface resistance components are normalized to a leaf area index
   !* of unity.
   !*
   !* Olson land types, deposition land types, and surface resistance
   !* components are read from file 'Olson_2001_Drydep_Inputs.nc'; check that file for
   !* further details.
   !***********************************************************************

   !           (1)        (2)        (3)         (4)         (5)        (6)      (7)     (8)     (9)     (10)   (11)
   !        snow/ice  deciduous  coniferous  agricultural  shub/     Amozaon   tundra  Desert  wetland  urban   water
   !                   forest     forest        land      grassland   forest
   ! Note IRI(3) = 200 is hardcoded and not the same as in the file.
   DATA IRI    /9999,    200,       200,          200,        200,       200,      200,   9999,   200,   9999,  9999/
   DATA IRLU   /9999,   9000,      9000,         9000,       9000,      1000,     4000,   9999,  9000,   9999,  9999/
   DATA IRAC   /   0,   2000,      2000,          200,        100,      2000,        0,      0,   300,    100,     0/
   DATA IRGSS  / 100,    500,       500,          150,        350,       200,      340,   1000,     0,    400,     0/
   DATA IRGSO  /3500,    200,       200,          150,        200,       200,      340,    400,  1000,    300,  2000/
   DATA IRCLS  /9999,   2000,      2000,         2000,       2000,      9999,     9999,   9999,  2500,   9999,  9999/
   DATA IRCLO  /1000,   1000,      1000,         1000,       1000,      9999,     9999,   9999,  1000,   9999,  9999/
   DATA IVSMAX /100,    100,       100,          100,        100,       100,      100,     10,   100,    100,    10/
   DATA IDRYDTYPE /1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11/
   ! Olson land use related parameters (https://wiki.seas.harvard.edu/geos-chem/index.php/Olson_land_map)
   DATA IOLSON /1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, &
      28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, &
      53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74/
   DATA IDEP  / 11,10, 5, 3, 3, 2, 2, 5, 8, 7,  5,  8,  1,  9,  11, 11, 5,   5,  5,  2, 6,  3,   3,  2,  2,  2,  2, &
      3, 6,   6,  4,  4,  2,  6,  2,  4,  9,  4,  4,  4,  5, 5,   5,  2,  5, 9,  5,   5,  2,  8,  8,  5, &
      5, 7,   2,  4,  2,  2,  2,  5,  2,  2,  3,  5,  5,  9, 9,   9,  9,  8, 8,  8,   9,  11/
   DATA IZO  /  10,  25000, 100,  10000, 10000, 10000, 10000, 100,  10, 2000, 100, 10,  1,   100,  1000,  1000,  1000, 100, 100, 2000, &
      10000,  10000,10000, 10000, 10000, 10000, 10000,10000,1000,10000,1000,1000,2000,10000,10000, 1000,  100, 1000, 1000,1000, &
      100,    100,   100, 2000,  100,    100,  1000, 1000, 1000, 1000,1000, 50,  50,  50,  2000,  2000, 2000, 2000, 1000, 100, &
      2000,   2000, 10000, 2000,  1000,  1000,  1000, 1000, 1000, 10,  1000,1000,500, 100 /
   ! Baldocchi polynomial coeffs
   DATA DRYCOEFF /-0.358, 3.02,  3.85, -0.0978,  -3.66,   12,   0.252,  -7.8,  0.226,  0.274,  &
      1.14,  -2.19,  0.261, -4.62,   0.685, -0.254, 4.37,  -0.266, -0.159, -0.206  /

   !There are 15 land types in Zhang et al., 2001 **aerosol deposition** scheme.
   !The land types in the model need to be mapped to these 15 land types.
   !=======================================================================
   !   #  LUC [Zhang et al., 2001]                GEOS-CHEM LUC (Corr. #)
   !-----------------------------------------------------------------------
   !   1 - Evergreen needleleaf trees             Snow/Ice          (12)
   !   2 - Evergreen broadleaf trees              Deciduous forest  ( 4)
   !   3 - Deciduous needleleaf trees             Coniferous forest ( 1)
   !   4 - Deciduous broadleaf trees              Agricultural land ( 7)
   !   5 - Mixed broadleaf and needleleaf trees   Shrub/grassland   (10)
   !   6 - Grass                                  Amazon forest     ( 2)
   !   7 - Crops and mixed farming                Tundra            ( 9)
   !   8 - Desert                                 Desert            ( 8)
   !   9 - Tundra                                 Wetland           (11)
   !  10 - Shrubs and interrupted woodlands       Urban             (15)
   !  11 - Wet land with plants                   Water             (14)
   !  12 - Ice cap and glacier
   !  13 - Inland water
   !  14 - Ocean
   !  15 - Urban
   !=======================================================================
   ! GEOS-CHEM LUC                 1, 2, 3, 4, 5, 6, 7  8, 9,10,11 (TODO:may add other land types later)
   INTEGER :: LUCINDEX_GC(11) = (/12, 4, 1, 7,10, 2, 9, 8,11,15,14/)

   !=======================================================================
   !   LUC       1,    2,    3,    4,    5,    6,    7,    8,
   !   alpha   1.0,  0.6,  1.1,  0.8,  0.8,  1.2,  1.2, 50.0,
   !   gamma  0.56, 0.58, 0.56, 0.56, 0.56, 0.54, 0.54, 0.54
   !
   !   LUC       9,   10,   11,   12,   13,   14,   15
   !   alpha  50.0,  1,3,  2.0, 50.0,100.0,100.0,  1.5
   !   gamma  0.54, 0.54, 0.54, 0.54, 0.50, 0.50, 0.56
   !=======================================================================
   REAL(fp)  :: ALPHA(15) = (/   1.0e+0_fp,   0.6e+0_fp,  1.1e+0_fp, &
      0.8e+0_fp,   0.8e+0_fp,  1.2e+0_fp, &
      1.2e+0_fp,  50.0e+0_fp, 50.0e+0_fp, &
      1.3e+0_fp,   2.0e+0_fp, 50.0e+0_fp, &
      100.0e+0_fp, 100.0e+0_fp,  1.5e+0_fp  /)

   REAL(fp)  :: GAMMA(15) = (/ 0.56e+0_fp, 0.58e+0_fp, 0.56e+0_fp, &
      0.56e+0_fp, 0.56e+0_fp, 0.54e+0_fp, &
      0.54e+0_fp, 0.54e+0_fp, 0.54e+0_fp, &
      0.54e+0_fp, 0.54e+0_fp, 0.54e+0_fp, &
      0.50e+0_fp, 0.50e+0_fp, 0.56e+0_fp  /)

   !...A unit is (mm) so multiply by 1.D-3 to (m)
   !   LUC       1,    2,    3,    4,    5,    6,    7,    8,
   !   SC1     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
   !   SC2     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
   ! A SC3     2.0,  5.0,  5.0, 10.0,  5.0,  5.0,  5.0,-999.,
   !   SC4     2.0,  5.0,  5.0, 10.0,  5.0,  5.0,  5.0,-999.,
   !   SC5     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
   !
   !   LUC       9,   10,   11,   12,   13,   14,   15
   !   SC1   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
   !   SC2   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
   ! A SC3   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
   !   SC4   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
   !   SC5   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
   REAL(fp)  :: A(15,5)

   DATA   A / 2.0e+0_fp,   5.0e+0_fp,   2.0e+0_fp,   5.0e+0_fp,  5.0e+0_fp, &
      2.0e+0_fp,   2.0e+0_fp, -999.e+0_fp, -999.e+0_fp, 10.0e+0_fp, &
      10.0e+0_fp, -999.e+0_fp, -999.e+0_fp, -999.e+0_fp, 10.0e+0_fp, &

      2.0e+0_fp,   5.0e+0_fp,   2.0e+0_fp,   5.0e+0_fp,  5.0e+0_fp, &
      2.0e+0_fp,   2.0e+0_fp, -999.e+0_fp, -999.e+0_fp, 10.0e+0_fp, &
      10.0e+0_fp, -999.e+0_fp, -999.e+0_fp, -999.e+0_fp, 10.0e+0_fp, &

      2.0e+0_fp,   5.0e+0_fp,   5.0e+0_fp,  10.0e+0_fp,  5.0e+0_fp, &
      5.0e+0_fp,   5.0e+0_fp, -999.e+0_fp, -999.e+0_fp, 10.0e+0_fp, &
      10.0e+0_fp, -999.e+0_fp, -999.e+0_fp, -999.e+0_fp, 10.0e+0_fp, &

      2.0e+0_fp,   5.0e+0_fp,   5.0e+0_fp,  10.0e+0_fp,  5.0e+0_fp, &
      5.0e+0_fp,   5.0e+0_fp, -999.e+0_fp, -999.e+0_fp, 10.0e+0_fp, &
      10.0e+0_fp, -999.e+0_fp, -999.e+0_fp, -999.e+0_fp, 10.0e+0_fp, &

      2.0e+0_fp,   5.0e+0_fp,   2.0e+0_fp,   5.0e+0_fp,  5.0e+0_fp, &
      2.0e+0_fp,   2.0e+0_fp, -999.e+0_fp, -999.e+0_fp, 10.0e+0_fp, &
      10.0e+0_fp, -999.e+0_fp, -999.e+0_fp, -999.e+0_fp, 10.0e+0_fp  /

   ! Annual average of A; put in the function now
   !REAL(fp)  :: Aavg(15)
   !Aavg(:) = (A(:,1)+A(:,2)+A(:,3)+A(:,4)+A(:,5))/5.

   ! Allocatable arrays for sea salt volume size bins
   REAL(fp),   ALLOCATABLE :: DMID    (:    )
   REAL(fp),   ALLOCATABLE :: SALT_V  (:    )
   !TODO:put sea salt size bins here for now; may be read from input file later
   real(fp), parameter :: SALA_REDGE_um(2)=(/0.01, 0.5/) !< accumulation mode Sea salt radius bin [um]
   real(fp), parameter :: SALC_REDGE_um(2)=(/0.5, 8.0/) !< coarse mode Sea salt radius bin [um]



contains
   !>
   !! \brief Computes the bulk surface resistance (Rc) for the gas species
   !!
   !!References:
   !! Wesely [1989]
   !!
   !! \param RADIAT      Solar radiation [W/m2]
   !! \param TEMP        Temperature [K]
   !! \param SUNCOS      Cosine of solar zenith angle
   !! \param F0          React. factor for oxidation depends on species
   !! \param HSTAR       Henry's law constant depends on species
   !! \param XMW         Molecular weight [kg/mol]
   !! \param USTAR       Friction velocity [m/s]
   !! \param CFRAC       Surface cloud fraction
   !! \param PRESSU      Surface pressure [Pa]
   !! \param XLAI        Leaf area index
   !! \param II          Index of the drydep land type
   !! \param SPC         Species name
   !! \param SALINITY    Salinity of the ocean
   !! \param TSKIN       Skin temperature
   !! \param IODIDE      Iodide concentration
   !! \param XLON        Longitude
   !! \param YLAT        Latitude
   !! \param CO2_EFFECT  CO2 effect on RS
   !! \param CO2_LEVEL   CO2 level
   !! \param CO2_REF     CO2 reference level
   !! \param RSURFC      Bulk Surface resistance [s/m]
   !! \param RC          Success or failure?
   !!
   !! \ingroup catchem_drydep_process
   !!!>
   subroutine Wesely_Rc_Gas(  RADIAT, TEMP, SUNCOS, F0, HSTAR, XMW, USTAR, CFRAC, PRESSU,     &
      XLAI,   II,   SPC, SALINITY, TSKIN, IODIDE, XLON, YLAT,         &
      CO2_EFFECT, CO2_LEVEL, CO2_REF, RSURFC,   RC)
      IMPLICIT NONE
      ! Parameters
      !-----------
      real(fp), intent(in)  :: RADIAT      !< Solar radiation [W/m2]
      real(fp), intent(in)  :: TEMP        !< Temperature [K]
      real(fp), intent(in)  :: SUNCOS      !< Cosine of solar zenith angle
      real(fp), intent(inout)  :: F0          !< React. factor for oxidation depends on species (inout because it is changed in the function)
      real(fp), intent(in)  :: HSTAR       !< Henry's law constant depends on species
      real(fp), intent(in)  :: XMW         !< Molecular weight [kg/mol]
      real(fp), intent(in)  :: USTAR       !< Friction velocity [m/s]
      real(fp), intent(in)  :: CFRAC       !< Surface cloud fraction
      real(fp), intent(in)  :: PRESSU      !< Surface pressure [Pa]
      real(fp), intent(in)  :: XLAI        !< Leaf area index
      integer,  intent(in)  :: II          !< Index of the drydep land type
      !integer,  intent(in)  :: N_SPC      !< Species ID (TODO: may be changed to species name)
      character(len=20), intent(in) :: SPC !< Species name
      !some inputs are for O3 over water and Hg over Amazon forest (not sure if we should include them for now)
      real(fp), intent(in)  :: SALINITY    !< Salinity of the ocean
      real(fp), intent(in)  :: TSKIN       !< Skin temperature
      real(fp), intent(in)  :: IODIDE      !< Iodide concentration
      real(fp), intent(in)  :: XLON        !< Longitude
      real(fp), intent(in)  :: YLAT        !< Latitude
      logical,  intent(in)  :: CO2_EFFECT  !< CO2 effect on RS
      real(fp), intent(in)  :: CO2_LEVEL   !< CO2 level
      real(fp), intent(in)  :: CO2_REF     !< CO2 reference level
      !output
      real(fp), intent(out) :: RSURFC      !< Bulk Surface resistance [s/m]
      integer, intent(out)  :: RC          !< Success or failure?

      ! Local Variables
      !----------------
      real(fp) :: RI, RLU, RAC, RGSS, RGSO, RCLS, RCLO
      real(fp) :: RT,RIX,GFACT,GFACI,RS_SCALE
      real(fp) :: RDC,RLUXX,RGSX,DTMP1,DTMP2,DTMP3,DTMP4
      real(fp) :: XMWH2O,TEMPK,TEMPC,DEPVw,alpha0
      real(fp) :: RCLX,RIXX    !,BIOFIT
      !string
      character(len=255)       :: thisLoc
      character(len=512)       :: ErrMsg

      !--------------------------------------------
      ! main function
      !--------------------------------------------

      ! Assume success
      RC      =  CC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> at Wesely_Rc_Gas (in process/drydep/CCPr_drydep_Commmon_Mod.F90)'

      ! Zero variables that aren't zeroed below
      RSURFC     = 0.0_fp
      RI         = 0.0_fp
      RLU        = 0.0_fp
      RAC        = 0.0_fp
      RGSS       = 0.0_fp
      RGSO       = 0.0_fp
      RCLS       = 0.0_fp
      RCLO       = 0.0_fp
      RIX        = 0.0_fp
      GFACT      = 0.0_fp
      GFACI      = 0.0_fp
      RDC        = 0.0_fp
      XMWH2O     = 0.0_fp
      RIXX       = 0.0_fp
      RLUXX      = 0.0_fp
      RGSX       = 0.0_fp
      RCLX       = 0.0_fp
      DTMP1      = 0.0_fp
      DTMP2      = 0.0_fp
      DTMP3      = 0.0_fp
      DTMP4      = 0.0_fp
      !N_SPC      = 0
      alpha0     = 0.0_fp
      DEPVw      = 0.0_fp

      !** TEMPK and TEMPC are surface air temperatures in K and in C
      TEMPK = TEMP
      TEMPC = TEMP - 273.15e+0_fp

      !* Adjust external surface resistances for temperature;
      !* from Wesely [1989], expression given in text on p. 1296.
      !*
      !* BUG FIX!  Wesely [1989] gives RT = 1000.0*EXP(-TEMPC-4.0)
      !*        RT = 1000.0*EXP(-(TEMPC-4.0))
      RT = 1000.0e+0_fp*EXP(-TEMPC-4.0e+0_fp)

      !If the surface to be snow or ice, set II to 1 instead.
      !IF((State_Met%isSnow(I,J)).OR.(State_Met%isIce(I,J))) II=1

      !************************************************************************
      !* Read the internal resistance RI (minimum stomatal resistance for
      !* water vapor,per unit area of leaf) from the IRI array; a '9999'
      !* value means no deposition to stomata so we impose a very large
      !* value for RI.
      !
      !*    Adjust stomatal resistances for insolation and temperature:
      !*     Temperature adjustment is from Wesely [1989], equation (3).
      !*
      !*     Light adjustment by the function BIOFIT is described by Wang
      !*     [1996]. It combines
      !*       - Local dependence of stomal resistance on the intensity I
      !*         of light impinging the leaf; this is expressed as a
      !*         multiplicative factor I/(I+b) to the stomatal resistance
      !*         where b = 50 W m-2 (equation (7) of Baldocchi et al.[1987])
      !*       - radiative transfer of direct and diffuse radiation in the
      !*         canopy using equations (12)-(16) from Guenther et al.[1995]
      !*       - separate accounting of sunlit and shaded leaves using
      !*         equation (12) of Guenther et al. [1995]
      !*       - partitioning of the radiation at the top of the canopy into
      !*         direct and diffuse components using a parameterization to
      !*         results from an atmospheric radiative transfer model
      !*         [Wang, 1996]
      !*     The dependent variables of the function BIOFIT are the leaf
      !*     area index (XYLAI), the cosine of zenith angle (SUNCOS) and
      !*     the fractional cloud cover (CFRAC).  The factor GFACI
      !*     integrates the light dependence over the canopy depth; sp even
      !*     though RI is input per unit area of leaf it need not be scaled
      !*     by LAI to yield a bulk canopy value because that's already
      !*     done in the GFACI formulation.
      !********************************************************************
      !RI = DBLE(IRI(II))
      RI = IRI(II)
      IF (RI   .GE. 9999.e+0_fp) THEN
         RI   = 1.e+12_fp
      ELSE
         GFACT = 100.0e+0_fp
         IF (TEMPC .GT. 0.e+0_fp .AND. TEMPC .LT. 40.e+0_fp) THEN
            GFACT = 400.e+0_fp/TEMPC/(40.0e+0_fp-TEMPC)
         ENDIF

         GFACI = 100.e+0_fp
         IF ( RADIAT > 0.e+0_fp .and. XLAI > 0.e+0_fp ) THEN
            GFACI = 1.e+0_fp / BIOFIT( DRYCOEFF,  XLAI, SUNCOS, CFRAC, SIZE(DRYCOEFF) )
         ENDIF
         RIX = RI*GFACT*GFACI
         ! Apply scaling factor to RIX when CO2 effect is turned
         ! on based on Franks et al. (2013)
         If (CO2_EFFECT) THEN
            RS_SCALE = CO2_LEVEL / CO2_REF *                   &
               (CO2_LEVEL + 80.0_fp) *          &
               (CO2_REF   - 40.0_fp) /          &
               (CO2_LEVEL - 40.0_fp) /          &
               (CO2_REF   + 80.0_fp)
            RIX = RIX * RS_SCALE
         ENDIF

      ENDIF

      !*Cuticular resistances IRLU array defined above are per unit area of leaf;
      !*divide them by the leaf area index to get a cuticular resistance for the bulk canopy.
      !*If IRLU is '9999' it means there are no cuticular surfaces on which to deposit so
      !*we impose a very large value for RLU.
      !TODO: not sure if XLAI is land type dependent or not.
      IF ( IRLU(II) >= 9999 .or. XLAI <= 0.e+0_fp ) THEN
         RLU = 1.e+6_fp
      ELSE
         !RLU = DBLE( IRLU(II) ) / XLAI
         RLU =  IRLU(II)  / XLAI
         ! Additional resistance at low temperatures.Limit increase to a factor of 2.
         ! Ref Jaegle et al. 2018
         RLU = MIN( RLU + RT, 2.e+0_fp * RLU )
      ENDIF

      !*The following are the remaining resistances for the Wesely model for a surface canopy
      !*(Wesely 1989, Fig.1).
      !RAC  = MAX(DBLE(IRAC(II)), 1.e+0_fp)
      RAC  = MAX(IRAC(II), 1.e+0_fp)
      IF (RAC  .GE. 9999.e+0_fp) RAC  = 1.e+12_fp
      !RGSS = MAX(DBLE(IRGSS(II)), 1.e+0_fp)
      RGSS = MAX(IRGSS(II), 1.e+0_fp)
      ! Additional resistance at low temperatures.Limit increase to a factor of 2.
      ! Ref Jaegle et al. 2018
      RGSS = MIN( RGSS + RT, 2.e+0_fp * RGSS )
      IF (RGSS .GE. 9999.e+0_fp) RGSS = 1.e12_fp
      !RGSO = MAX(DBLE(IRGSO(II)) ,1.e+0_fp)
      RGSO = MAX(IRGSO(II) ,1.e+0_fp)
      RGSO = MIN( RGSO + RT, 2.e+0_fp * RGSO )
      IF (RGSO .GE. 9999.e+0_fp) RGSO = 1.e+12_fp
      !RCLS = DBLE(IRCLS(II))
      RCLS = IRCLS(II)
      RCLS = MIN( RCLS + RT, 2.e+0_fp * RCLS )
      IF (RCLS .GE. 9999.e+0_fp) RCLS = 1.e+12_fp
      !RCLO = DBLE(IRCLO(II))
      RCLO = IRCLO(II)
      RCLO = MIN( RCLO + RT, 2.e+0_fp * RCLO )
      IF (RCLO .GE. 9999.e+0_fp) RCLO = 1.e+12_fp

      !* Compute aerodynamic resistance to lower elements in lower part
      !* of the canopy or structure, assuming level terrain -
      !* equation (5) of Wesely [1989].
      !* species-dependent corrections to resistances
      !* are from equations (6)-(9) of Wesely [1989].

      RDC = 100.e+0_fp*(1.0e+0_fp+1000.0e+0_fp/(RADIAT+10.e+0_fp))

      IF ( SPC .EQ. 'O3' ) THEN
         !O3 over water
         IF ((II .EQ. 11)) THEN
            IF (SALINITY .GT. 20.0_fp) THEN
               ! Now apply the Luhar et al. [2018] equations for the
               ! special treatment of O3 dry deposition to the ocean
               CALL OCEANO3(TSKIN,USTAR,IODIDE,DEPVw)
               ! Now convert to the new rc value(s)
               alpha0 = 10.0_fp**(-0.25-0.013 * (TSKIN-273.16_fp))
               RSURFC = 1.0_fp/(alpha0*DEPVw)
            ELSE
               ! It's not saline enough for 'ocean' so we instead don't change it from
               ! 'default' rc to water
               RSURFC = 2000.0_fp
            ENDIF
         ENDIF

         !O3 over snow/ice, the surface resistance is set to an observation derived value
         IF ((II .EQ. 1)) THEN
            RSURFC = 10000.0_fp
         ENDIF
      ELSE
         !set a different F0 for Hg0
         IF (SPC .EQ. 'Hg0') THEN
            ! Assume lower reactivity
            F0 = 3.0e-05_fp
            ! But if this is the rainforest land type and we fall
            ! within the bounding box of the Amazon rainforest,
            ! then increase reactivity as inferred from observations.
            IF ( II  ==  6          .AND.             &
               XLON >  -82.0_fp   .AND.             &
               XLON <  -33.0_fp   .AND.             &
               YLAT >  -34.0_fp   .AND.             &
               YLAT <   14.0_fp ) THEN
               F0 = 2.0e-01_fp
            ENDIF
         ENDIF

         XMWH2O = H2OMW * 1.e-3_fp
         RIXX = RIX*DIFFG(TEMPK,PRESSU,XMWH2O)/ DIFFG(TEMPK,PRESSU,XMW) &
            + 1.e+0_fp/(HSTAR/3000.e+0_fp+100.e+0_fp*F0)
         RLUXX = 1.e+12_fp
         IF (RLU .LT. 9999.e+0_fp) RLUXX = RLU/(HSTAR/1.0e+05_fp + F0)
         RGSX = 1.e+0_fp/(HSTAR/1.0e+05_fp/RGSS + F0/RGSO)
         RCLX = 1.e+0_fp/(HSTAR/1.0e+05_fp/RCLS + F0/RCLO)
         !** Get the bulk surface resistance of the canopy, RSURFC, from
         !** the network of resistances in parallel and in series (Fig.1 of Wesely [1989])
         DTMP1=1.e+0_fp/RIXX
         DTMP2=1.e+0_fp/RLUXX
         DTMP3=1.e+0_fp/(RAC+RGSX)
         DTMP4=1.e+0_fp/(RDC+RCLX)
         RSURFC = 1.e+0_fp/(DTMP1 + DTMP2 + DTMP3 + DTMP4)
      ENDIF

      !TODO: this should be put in the main scheme function since it is also applied to aerosols
      !*Set max and min values for bulk surface resistances
      !!RSURFC = MAX(1.e+0_fp, MIN(RSURFC,9999.e+0_fp))
      !*because of high resistance values, different rule applied for ocean ozone
      !!IF ((SPC .EQ. 'O3') .AND. (II .EQ. 11)) THEN
      !!   RSURFC = MAX(1.e+0_fp, MIN(RSURFC,999999.e+0_fp))
      !!ENDIF
      ! Set Rc for strong acids (HNO3,HCl,HBr) to 1 s/m
      ! Ref. Jaegle et al. 2018, cf. Erisman,van Pul,Ayers 1994
      !!IF ( HSTAR .gt. 1.e+10_fp ) RSURFC= 1.e+0_fp

      return
   end subroutine Wesely_Rc_Gas

   !>
   !! \brief calculates the dry deposition velocity of O3 to ocean
   !!
   !!References:
   !! Pound, R. J., Sherwen, T., Helmig, D., Carpenter, L. J., and Evans, M. J.:
   !! Influence of oceanic ozone deposition on tropospheric photochemistry,
   !! Atmos. Chem. Phys., https://doi.org/10.5194/acp-20-4227-2020, 2020.
   !!
   !! \param TEMPK      Temperatue [K]
   !! \param USTAR      Fictional Velocity [m/s]
   !! \param IODIDE_IN  Surface iodide concentration [nM]
   !! \param DEPV       output of the new deposition vel [cm/s]
   !!
   !! \ingroup catchem_drydep_process
   !!!>
   SUBROUTINE OCEANO3( TEMPK, USTAR, IODIDE_IN, DEPV )

      IMPLICIT NONE

      !INPUT PARAMETERS:
      REAL(fp), INTENT(IN)         :: TEMPK ! Temperature [K]
      REAL(fp), INTENT(IN)         :: USTAR ! Fictional Velocity [m/s]
      REAL(fp), INTENT(IN)         :: IODIDE_IN ! Surface iodide concentration [nM]
      REAL(fp), INTENT(OUT)        :: DEPV  ! the new deposition vel [cm/s]
      !LOCAL VARIABLES:
      REAL(fp) :: a0,D,DelM,b,PSI,LAM,EP,USTARWater,K0,K1,Iodide

      !=================================================================
      ! OCEANO3 begins here!
      !=================================================================

      USTARWater = 0.0345_fp * USTAR !waterside friction velocity

      Iodide = IODIDE_IN*1.0E-9_fp ! Convert from nM to M

      a0 = Iodide*EXP((-8772.2/TEMPK)+51.5) !chemical reactivity

      D = 1.1E-6*EXP(-1896.0/TEMPK) ! diffusivity

      DelM = SQRT(D/a0) ! reaction-diffusion length

      b = 2.0_fp/(0.4_fp*USTARWater)

      LAM = DelM*SQRT(a0/D) ! this cancels to 1 but here for completeness of equations

      EP = SQRT(2.0_fp*a0*b*(DelM+(b*D/2.0_fp)))

      PSI = EP/SQRT(a0*b**2*D)

      CALL K0K1_APROX(EP,K0,K1)

      DEPV = SQRT(a0*D)*((PSI*K1*COSH(LAM)+K0*SINH(LAM))/(PSI*K1* SINH(LAM)+K0*COSH(LAM)))

   END SUBROUTINE OCEANO3

   !>
   !! \brief estimate the modified Bessel functions of the second order zero (K0) and one (K1).
   !!
   !!References:
   !! Approach initially described in Numerical Recipes in Fortran 90 second edition
   !! (1996). This implementation is designed to be specific to the use
   !! case required for calculating oceanic deposition velocity. Uses a
   !! polynomial fit of each type of modified bessel function to
   !! estimate the value of the function for each input.
   !!
   !! \param input_arg    !the value we want the soln for
   !! \param K0, K1       output of the modified bessel functions
   !!
   !! \ingroup catchem_drydep_process
   !!!>
   SUBROUTINE K0K1_APROX( input_arg, K0, K1 )

      IMPLICIT NONE
      !INPUT PARAMETERS:
      REAL(fp), INTENT(IN)  :: input_arg !the value we want the soln for
      REAL(fp), INTENT(OUT) :: K0,K1     !the values of the modified bessel fncs
      !LOCAL VARIABLES:
      REAL(fp), DIMENSION(7) :: coeff !coefficients for polynomial fit
      ! of each bessel function
      REAL(fp)               :: I0,I1 !modified bessel functions of
      ! first kind order 0 and 1

      ! determine which fit method is best for the bessel functions
      IF (input_arg <= 2.0_fp) THEN
         ! begin the calculation of k0 by estimating i0
         coeff = (/1.0,3.5156229,3.0899424,1.2067492,0.2659732, &
            0.360768e-1,0.45813e-2/)
         I0 = poly_fit((input_arg/3.75_fp)**2,coeff)
         !now we can use this estimate of i0 to calculate k0
         coeff = (/-0.57721566,0.42278420,0.23069756,0.3488590e-1, &
            0.262698e-2,0.10750e-3,0.74e-5/)
         K0 = (-LOG(0.5_fp*input_arg)*I0)+ &
            poly_fit(0.25_fp*input_arg**2,coeff)

         !begin the calculation of k0 by estimating i1
         coeff = (/0.5,0.87890594,0.51498869,0.15084934,0.2658733e-1, &
            0.301532e-2,0.32411e-3/)
         I1 = input_arg*poly_fit((input_arg/3.75_fp)**2,coeff)
         ! now we can use this to estimate to get a value for k1
         coeff = (/1.0,0.15443144,-0.67278579,-0.18156897, &
            -0.1919402e-1,-0.110404e-2,-0.4686e-4/)
         K1 = (LOG(0.5_fp*input_arg)*I1)+(1.0_fp/input_arg)* &
            poly_fit(0.25_fp*input_arg**2,coeff)
      ELSE !use a different approximation that doesn't need I0/I1
         coeff = (/1.25331414,-0.7832358e-1,0.2189568e-1,-0.1062446e-1, &
            0.587872e-2,-0.251540e-2,0.53208e-3/)
         K0 = (EXP(-input_arg)/SQRT(input_arg))* &
            poly_fit((2.0_fp/input_arg),coeff)
         coeff = (/1.25331414,0.23498619,-0.3655620e-1,0.1504268e-1, &
            -0.780353e-2,0.325614e-2,-0.68245e-3/)
         K1 = (EXP(-input_arg)/SQRT(input_arg))* &
            poly_fit((2.0_fp/input_arg),coeff)
      ENDIF

   END SUBROUTINE K0K1_APROX

   !>
   !! \brief calculate the value of a polynomial fit used in
   !! the K0K1_APPROX function in estimating the values of a
   !! modified bessel function.
   !!
   !!References:
   !!
   !! \param input
   !! \param coeffs
   !!
   !! \ingroup catchem_drydep_process
   !!!>
   FUNCTION poly_fit ( input, coeffs )

      !INPUT PARAMETERS:
      REAL(fp), INTENT(IN)               :: input
      REAL(fp), DIMENSION(:), INTENT(IN) :: coeffs
      !LOCAL VARIABLES:
      REAL(fp)                           :: poly_fit
      INTEGER                            :: i

      poly_fit = 0

      DO i = 1,7,1
         poly_fit = poly_fit+coeffs(i)*input**i
      ENDDO

   END FUNCTION poly_fit

   !>
   !! \brief calculates the molecular diffusivity [m2/s] in air for a gas X
   !!  of molecular weight XM [kg] at temperature TK [K] and pressure PRESS [Pa].
   !!
   !!References:
   !!
   !! \param TK      Temperatue [K]
   !! \param PRESS   Pressure [Pa]
   !! \param XM      Molecular weight of gas [kg]
   !!
   !! \ingroup catchem_drydep_process
   !!!>

   FUNCTION DIFFG( TK, PRESS, XM ) RESULT( DIFF_G )

      !INPUT PARAMETERS:
      REAL(fp), INTENT(IN) :: TK     ! Temperature [K]
      REAL(fp), INTENT(IN) :: PRESS  ! Pressure [Pa]
      REAL(fp), INTENT(IN) :: XM     ! Molecular weight of gas [kg]
      !LOCAL VARIABLES:
      REAL(fp)             :: AIRDEN, Z, DIAM, FRPATH, SPEED, DIFF_G

      !REMARKS:
      !We specify the molecular weight of air (XMAIR) and the hard-sphere molecular
      !radii of air (RADAIR) and of the diffusing gas (RADX).  The molecular
      !radius of air is given in a Table on p. 479 of Levine [1988].  The Table
      !also gives radii for some other molecules.  Rather than requesting the user
      !to supply a molecular radius we specify here a generic value of 1.2E-10 m for
      !all molecules, which is good enough in terms of calculating the diffusivity
      !as long as molecule is not too big.

      !DEFINED PARAMETERS:
      REAL(fp), PARAMETER  :: XMAIR  = 28.8e-3_fp ! Moist air molec wt?
      REAL(fp), PARAMETER  :: RADAIR = 1.2e-10_fp
      REAL(fp), PARAMETER  :: RADX   = 1.5e-10_fp

      !=================================================================
      ! DIFFG begins here!
      !=================================================================

      ! Air density [molec/m3]
      AIRDEN = ( PRESS * AVO ) / ( RSTARG * TK )

      ! DIAM is the collision diameter for gas X with air.
      DIAM   = RADX + RADAIR

      ! Calculate the mean free path for gas X in air:
      ! eq. 8.5 of Seinfeld [1986];
      Z      = XM  / XMAIR
      FRPATH = 1e+0_fp /( PI * SQRT( 1e+0_fp + Z ) * AIRDEN * ( DIAM**2 ) )

      ! Calculate average speed of gas X; eq. 15.47 of Levine [1988]
      SPEED  = SQRT( 8e+0_fp * RSTARG * TK / ( PI * XM ) )

      ! Calculate diffusion coefficient of gas X in air;
      ! eq. 8.9 of Seinfeld [1986]
      DIFF_G = ( 3e+0_fp * PI / 32e+0_fp ) * ( 1e+0_fp + Z ) * FRPATH * SPEED

   END FUNCTION DIFFG

   !>
   !! \brief computes the light correction used in the dry deposition and canopy NOx modules.
   !! It was part of the old Harvard-GISS CTM and was ported into GEOS-Chem
   !!
   !!References:
   !! Wang, Y., D.J. Jacob, and J.A. Logan, "Global simulation of tropospheric
   !! O3-NOx-hydrocarbon chemistry, 1. Model formulation", J. Geophys. Res.,
   !! 103/D9, 10,713-10,726, 1998.
   !!
   !! \param COEFF1     Baldocchi drydep coefficients
   !! \param XLAI1      Leaf area index [cm2/cm2]
   !! \param SUNCOS1    Cosine( Solar Zenith Angle )
   !! \param CFRAC1     Cloud fraction [unitless]
   !! \param NPOLY      # of drydep coefficients
   !!
   !! \ingroup catchem_drydep_process
   !!!>
   FUNCTION BioFit( COEFF1, XLAI1, SUNCOS1, CFRAC1, NPOLY ) RESULT( BIO_FIT )

      !INPUT PARAMETERS:
      INTEGER,   INTENT(IN) :: NPOLY           ! # of drydep coefficients
      REAL(fp),  INTENT(IN) :: COEFF1(NPOLY)   ! Baldocchi drydep coefficients
      REAL(fp),  INTENT(IN) :: XLAI1           ! Leaf area index [cm2/cm2]
      REAL(fp),  INTENT(IN) :: SUNCOS1         ! Cosine( Solar Zenith Angle )
      REAL(fp),  INTENT(IN) :: CFRAC1          ! Cloud fraction [unitless]
      !RETURN VALUE:
      REAL(fp)              :: BIO_FIT         ! Resultant light correction
      !DEFINED PARAMETERS:
      INTEGER, PARAMETER    :: KK = 4
      INTEGER, PARAMETER    :: NN = 3  ! # of variables (LAI, SUNCOS, CLDFRC)
      REAL(fp)  :: ND(NN) = (/ 55.0e0_fp, 20.0e0_fp, 11.0e0_fp /) !scaling factor for each variable !codespell:ignore
      REAL(fp)  :: X0(NN) = (/ 11.0e0_fp, 1.0e0_fp,  1.0e0_fp /) !maximum for each variable
      !LOCAL VARIABLES:
      REAL(fp)              :: XLOW !minimum for each variable
      REAL(fp)              :: TERM(KK)
      REAL(fp)              :: REALTERM(NPOLY)
      INTEGER               :: K,K1,K2,K3,I,I2

      !=================================================================
      ! BIOFIT begins here!
      !=================================================================
      TERM(1) = 1.0e0_fp
      TERM(2) = XLAI1
      TERM(3) = SUNCOS1
      TERM(4) = CFRAC1
      !we replace SUNPARAM_R4( TERM(2:4) ) as below
      !outdate lai,suncos,cloud fraction
      DO I = 1, NN
         I2 = I + 1 !variable index in TERM is from 2
         TERM(I2) = MIN( TERM(I2), X0(I) )
         ! XLOW = minimum for each variable
         IF ( I .NE. 3 ) THEN
            XLOW = X0(I) / ND(I) !codespell:ignore
         ELSE
            XLOW = 0.0e0_fp
         ENDIF
         TERM(I2) = MAX( TERM(I2), XLOW )
         TERM(I2) = TERM(I2) / X0(I)
      ENDDO

      !get realterm
      K = 0
      DO K3 = 1, KK
         DO K2 = K3, KK
            DO K1 = K2, KK
               K = K + 1
               REALTERM(K)=TERM(K1)*TERM(K2)*TERM(K3)
            ENDDO
         ENDDO
      ENDDO

      BIO_FIT = 0e0_fp
      DO K = 1, NPOLY
         BIO_FIT = BIO_FIT + COEFF1(K)*REALTERM(K)
      END DO
      IF ( BIO_FIT .LT. 0.1e0_fp ) BIO_FIT = 0.1e0_fp

   END FUNCTION BioFit


   !>
   !! \brief computes the aerodynamic resistance of aerosols.
   !!
   !!References:
   !! Zhang, L., Gong, S., Padro, J., & Barrie, L. (2001). A size-segregated particle
   !! dry deposition scheme for an atmospheric aerosol module.
   !! Atmospheric environment., https://doi.org/10.1016/S1352-2310(00)00326-5
   !!
   !! Emerson, E. W., et al. (2020). Revisiting particle dry deposition and its role
   !! in radiative effect estimates. PNAS, 117(42), 26076-26082.
   !! https://doi.org/10.1073/pnas.2014761117
   !!
   !! \param SPC        Species name
   !! \param II         Surface type index
   !! \param IS_DUST    Is dust species?
   !! \param IS_SEASALT Is seasalt species?
   !! \param LUCINDEX   mapping above II to the 15 drydep land use categories
   !! \param A_RADI     Aerosol radius [m]
   !! \param A_DEN      Aerosol density [kg/m3]
   !! \param PRESS      Pressure [KPa]
   !! \param TEMP       Temperature [K]
   !! \param USTAR      Fictional Velocity [m/s]
   !! \param RHB        Relative humidity [fraction]
   !! \param W10        10m wind speed [m/s]
   !! \param VTSout     output of setttling velocity [m/s]
   !! \param RC         success flag
   !! \param RS         return value of surface resistance [s/m]
   !!
   !! \ingroup catchem_drydep_process
   !!!>

   FUNCTION AERO_SFCRSII(  SPC, II, IS_DUST, IS_SEASALT, LUCINDEX, A_RADI, A_DEN, &
      PRESS, TEMP, USTAR, RHB, W10, VTSout, RC) RESULT( RS )

      IMPLICIT NONE
      !INPUT PARAMETERS
      CHARACTER(len=20), INTENT(IN) :: SPC    ! Species name
      !TODO: not sure if SPC or index is better
      !INTEGER,  INTENT(IN) :: K    ! Drydep species index (range: 1-NUMDEP)
      INTEGER,  INTENT(IN) :: II    ! Surface type index of host model (e.g., GEOS-CHEM)
      LOGICAL,  INTENT(IN) :: IS_DUST, IS_SEASALT ! Is dust or seasalt species?
      INTEGER,  DIMENSION(:), INTENT(IN) :: LUCINDEX !mapping above II to the 15 drydep land use categories
      REAL(fp), INTENT(IN) :: A_RADI ! Aerosol radius [m]
      REAL(fp), INTENT(IN) :: A_DEN  ! Aerosol density [kg/m3]
      REAL(fp), INTENT(IN) :: PRESS ! Pressure [kPa] (1 mb = 100 Pa = 0.1 kPa)
      REAL(fp), INTENT(IN) :: TEMP  ! Temperature [K]
      REAL(fp), INTENT(IN) :: USTAR ! Friction velocity [m/s]
      REAL(fp), INTENT(IN) :: RHB   ! Relative humidity (fraction)
      REAL(fp), INTENT(IN) :: W10   ! 10m wind speed [m/s]; only need for SeaSalt over water
      !OUTPUT PARAMETERS
      REAL(fp), INTENT(OUT) :: VTSout ! Settling velocity [m/s]
      INTEGER,  INTENT(OUT) :: RC     ! success flag
      !RETURN VALUE
      REAL(fp)             :: RS    ! Surface resistance for particles [s/m]

      !define constants
      REAL(fp), PARAMETER   :: E0       =  3.0_fp
      ! Emerson et al. (2020) added parameters
      REAL(fp), PARAMETER   :: UPSILON  =  0.8_fp
      REAL(fp), PARAMETER   :: BETA     =  1.7_fp
      REAL(fp), PARAMETER   :: CB       =  0.2_fp
      REAL(fp), PARAMETER   :: CIN      =  2.5_fp
      REAL(fp), PARAMETER   :: CIM      =  0.4_fp
      !increment of radius for integration of settling velocity (um)
      REAL(fp), PARAMETER   :: DR       =  5.0e-2_fp
      !LOCAL VARIABLES
      INTEGER   :: LUC
      INTEGER   :: ID,NR
      REAL(fp)  :: AIRVS       ! kinematic viscosity of Air (m^2/s)
      REAL(fp)  :: DP          ! Diameter of aerosol [um]
      REAL(fp)  :: PDP         ! Press * Dp
      REAL(fp)  :: CONST       ! Constant for settling velocity calculations
      REAL(fp)  :: SLIP        ! Slip correction factor
      REAL(fp)  :: VISC        ! Viscosity of air (Pa s)
      REAL(fp)  :: SALT_MASS, SALT_MASS_TOTAL, VTS_WEIGHT, DMIDW
      real(fp)  :: D0, D1      !lower and upper bounds of sea-salt dry diameter bins
      REAL(fp)  :: DIFF        ! Brownian Diffusion constant for particles (m2/s)
      REAL(fp)  :: SC, ST      ! Schmidt and Stokes number (nondim)
      REAL(fp)  :: DIAM, RDRY, RWET, RUM, DEN
      REAL(fp)  :: EB, EIM, EIN, R1, AA, VTS
      REAL(fp)  :: RHBL        ! Relative humidity local
      REAL(fp)  :: Aavg(15)    ! annual average of A
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at AERO_SFCRSII (in process/drydep/ccpr_drydep_common_mod.F90)'

      !=================================================================
      ! ADUST_SFCRII begins here!
      !=================================================================

      ! Annual average of A
      Aavg(:) = (A(:,1)+A(:,2)+A(:,3)+A(:,4)+A(:,5))/5.

      LUC     = LUCINDEX(II)
      AA      = Aavg(LUC) * 1.e-3_fp
      RS = 0e+0_fp !initialize returned value first

      !=================================================================
      !...Ref. Zhang et al., AE 35(2001) 549-560
      !.
      !...Model theroy
      !    Vd = Vs + 1./(Ra+Rs)
      !      where Vs is the gravitational settling velocity,
      !      Ra is the aerodynamic resistance above the canopy
      !      Rs  is the surface resistance
      !    Here we calculate Rs only..
      !    Rs = 1 / (Eo*Ustar*(Eb+Eim+Ein)*R1)
      !      where Eo is an empirical constant ( = 3.)
      !      Ustar is the friction velocity
      !      Collection efficiency from
      !        Eb,  [Brownian diffusion]
      !        Eim, [Impaction]
      !        Ein, [Interception]
      !      R1 is the correction factor representing the fraction
      !         of particles that stick to the surface.
      !=======================================================================
      !      Eb is a funciont of Schmidt number, Eb = Sc^(-gamma)
      !         Sc = v/D, v (the kinematic viscosity of air)
      !                   D (particle brownian diffusivity)
      !         r usually lies between 1/2 and 2/3
      !      Eim is a function of Stokes number, St
      !          St = Vs * Ustar / (g0 * A)   for vegetated surfaces
      !          St = Vs * Ustar * Ustar / v  for smooth surface
      !          A is the characteristic radius of collectors.
      !
      !       1) Slinn (1982)
      !           Eim = 10^(-3/St)          for smooth surface
      !           Eim = St^2 / ( 1 + St^2 ) for vegetative canopies
      !       2) Peters and Eiden (1992)
      !           Eim = ( St / ( alpha + St ) )^(beta)
      !                alpha(=0.8) and beta(=2) are constants
      !       3) Giorgi (1986)
      !           Eim = St^2 / ( 400 + St^2 )     for smooth surface
      !           Eim = ( St / (0.6 + St) )^(3.2) for vegetative surface
      !       4) Davidson et al.(1982)
      !           Eim = St^3 / (St^3+0.753*St^2+2.796St-0.202) for grassland
      !       5) Zhang et al.(2001) used 2) method with alpha varying with
      !          vegetation type and beta equal to 2
      !
      !      Ein = 0.5 * ( Dp / A )^2
      !
      !      R1 (Particle rebound)  = exp(-St^0.5)
      !=================================================================

      ! Particle diameter [m]
      ! A_RADI & A_DEN are read from inputs.
      DIAM  = A_RADI * 2.e+0_fp

      ! Particle density [kg/m3]
      DEN   = A_DEN

      !update DIAM of dust species; no hygroscopic growth for dust
      !TODO: diameter for dust is hardcoded here; may need to change it to be flexible for different dust bins
      !IF ( K == idd_DST1 .or. K == idd_DSTAL1 .or. K == idd_NITD1 .or. K == idd_SO4D1 ) THEN
      IF ( SPC == 'DST1' .or. SPC == 'DSTAL1' .or. SPC == 'NITD1' .or. SPC == 'SO4D1' .or. SPC == 'dust1') THEN
         DIAM = 0.66895E-6
      ENDIF

      !IF ( K == idd_DST2 .or. K == idd_DSTAL2 .or. K == idd_NITD2 .or. K == idd_SO4D2 ) THEN
      IF ( SPC == 'DST2' .or. SPC == 'DSTAL2' .or. SPC == 'NITD2' .or. SPC == 'SO4D2' .or. SPC == 'dust2') THEN
         DIAM = 2.4907E-6
      ENDIF

      !IF ( K == idd_DST3  .or. K == idd_DSTAL3 .or. K == idd_NITD3 .or. K == idd_SO4D3 ) THEN
      IF ( SPC == 'DST3'  .or. SPC == 'DSTAL3' .or. SPC == 'NITD3' .or. SPC == 'SO4D3' ) THEN
         DIAM = 4.164E-6
      ENDIF

      !IF ( K == idd_DST4  .or. K == idd_DSTAL4 .or. K == idd_NITD4 .or. K == idd_SO4D4 ) THEN
      IF ( SPC == 'DST4'  .or. SPC == 'DSTAL4' .or. SPC == 'NITD4' .or. SPC == 'SO4D4' ) THEN
         DIAM = 6.677E-6
      ENDIF

      ! Hygroscopic growth following Latimer and Martin (2019) ACP
      RHBL    = MAX( TINY(RHB), RHB )

      ! Over oceans the RH in the viscous sublayer is set to 98%,
      ! following Lewis and Schwartz (2004)
      IF (LUC == 14) THEN
         RHBL = 0.98
      ENDIF

      IF (.NOT. IS_DUST) THEN
         !update DIAM and DEN after hygroscopic growth for non-dust species
         call New_DIAM_DEN( SPC, IS_SEASALT, RHBL, RDRY, RWET, DIAM, DEN, RC)
         if (RC /= CC_SUCCESS ) then
            errMsg = 'New_DIAM_DEN failed.'
            CALL CC_Error( errMsg, RC, thisLoc )
            RETURN
         endif
      ENDIF

      ! Dp [m] --> [um] = particle diameter
      DP    = DIAM * 1.e+6_fp

      ! Constant for settling velocity calculation
      CONST = DEN * DIAM**2 * g0 / 18.e+0_fp

      !=================================================================
      ! Slip correction factor calculations following Seinfeld,
      ! pp464 which is thought to be more accurate but more computation
      ! required.
      !   # air molecule number density
      !   num = P * 1d3 * 6.023d23 / (8.314 * Temp)
      !   # gas mean free path
      !   lambda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 )
      !   # Slip correction
      !   Slip = 1. + 2. * lambda * (1.257 + 0.4 * exp( -1.1 * Dp &
      !          / (2. * lambda))) / Dp
      !
      ! Note, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
      ! which produce slip correction factore with small error
      ! compared to the above with less computation.
      !=================================================================

      ! Slip correction factor as function of (P*dp)
      PDP  = PRESS * DP
      SLIP = 1e+0_fp + ( 15.60e+0_fp + 7.0e+0_fp * &
         EXP( -0.059e+0_fp * PDP) ) / PDP

      ! Viscosity [Pa s] of air as a function of temp (K)
      VISC = 1.458e-6_fp * (TEMP)**(1.5e+0_fp) / (TEMP + 110.4e+0_fp)

      ! Kinematic viscosity (Dynamic viscosity/Density)
      AIRVS= VISC / 1.2928e+0_fp

      ! Settling velocity [m/s]
      VTS  = CONST * SLIP / VISC
      !sea salt VTS update
      IF (IS_SEASALT) THEN
         ! This settling velocity is for the mid-point of the size bin.
         ! Need to integrate over the size bin, taking into account the
         ! mass distribution of sea-salt and the dependence of VTS on aerosol
         ! size. See WET_SETTLING in SEASALT_MOD.f for more details.

         !TODO: this may be used in initialization of the scheme
         !Number of bins for sea salt size distribution
         !SALA_radius_bin_in_um: [0.01, 0.5];  SALC_radius_bin_in_um: [0.5,  8.0]

         ! Make sure that SALA, SALC bins are contiguous
         !IF ( SALA_REDGE_um(2) /= SALC_REDGE_um(1) ) THEN
         !   MSG = 'SALA and SALC bin edges are not contiguous!'
         !   CALL ERROR_STOP( MSG, LOCATION )
         !ENDIF
         !TODO: need to figure out how to read in these values from the namelist
         NR = INT((( SALC_REDGE_um(2) - SALA_REDGE_um(1) ) &
            / DR ) + 0.5e+0_fp )

         SALT_MASS_TOTAL = 0e+0_fp
         VTS_WEIGHT      = 0e+0_fp

         ! Dry particle radius [m] --> [um]
         RUM  = RDRY * 1.e+6_fp

         ! Check what the min/max range of the SS size bins are
         IF ( RUM .le. SALA_REDGE_um(2) ) THEN
            D0 = SALA_REDGE_um(1)*2e+0_fp
            D1 = SALA_REDGE_um(2)*2e+0_fp
         ELSE
            D0 = SALC_REDGE_um(1)*2e+0_fp
            D1 = SALC_REDGE_um(2)*2e+0_fp
         ENDIF

         DO ID = 1, NR
            ! Calculate mass of wet aerosol (Dw = wet diameter, D = dry diameter):
            ! Overall = dM/dDw = dV/dlnD * Rwet/Rdry * DEN /Rw
            !TODO: DMID is not defined in this module. Need to define it.
            IF (DMID(ID) .ge. D0 .and. DMID(ID) .le. D1 ) THEN
               DMIDW = DMID(ID) * RWET/RDRY   ! wet radius [um]
               SALT_MASS   = SALT_V(ID) * RWET/RDRY * DEN / &
                  (DMIDW*0.5e+0_fp)
               VTS_WEIGHT  = VTS_WEIGHT + &
               !SALT_MASS * VTS * (DMIDW/(RWET*1d6*2e+0_fp) )** &
                  SALT_MASS * VTS * (DMIDW/(RWET*1e+6_fp*2e+0_fp) )** &
                  2e+0_fp * (2e+0_fp * DR *  RWET/RDRY)
               SALT_MASS_TOTAL = SALT_MASS_TOTAL+SALT_MASS * &
                  (2e+0_fp * DR *  RWET/RDRY)
            ENDIF
         ENDDO

         ! Final mass weighted setting velocity:
         VTS = VTS_WEIGHT/SALT_MASS_TOTAL
      END IF

      VTSout = VTS !need to save out for final Vd calculation

      ! Brownian diffusion constant for particle (m2/s)
      DIFF = BOLTZ * TEMP * SLIP / (3.e+0_fp * PI * VISC * DIAM)

      ! Schmidt number
      SC   = AIRVS / DIFF
      !EB   = 1.e+0_fp/SC**(gamma(LUC))

      !--------------------------------------------------------------
      ! NOTE: This loses precision, use TWO_THIRDS parameter instead
      !EB   = CB/SC**(0.6667e+0_fp) ! Emerson 2020 update JRP
      !--------------------------------------------------------------
      EB   = CB/SC**TWO_THIRDS ! Emerson 2020 update JRP

      ! Stokes number
      IF ( AA < 0e+0_fp ) then
         ST   = VTS * USTAR * USTAR / ( AIRVS * g0 ) ! for smooth surface
         EIN  = 0e+0_fp
      ELSE
         ST   = VTS   * USTAR / ( g0 * AA )          ! for vegetated surfaces
         !EIN  = 0.5e+0_fp * ( DIAM / AA )**2
         EIN  = CIN * ( DIAM / AA )**(UPSILON) ! Emerson 2020 update JRP
      ENDIF

      IF (LUC == 14 .and. IS_SEASALT) THEN
         EIM  = 10.e+0_fp**( -3.e+0_fp/ ST )         ! for water surface
         ! JRP: Emerson doesn't describe what to do here, so I'm leaving as is
      ELSE
         !EIM  = ( ST / ( ALPHA(LUC) + ST ) )**(BETA)
         EIM  = CIM * ( ST / ( ALPHA(LUC) + ST ) )**(BETA) ! Emerson 2020 update JRP
         EIM  = MIN( EIM, 0.6e+0_fp )
      ENDIF

      IF (LUC == 11 .OR. LUC == 13 .OR. LUC == 14) THEN
         R1 = 1.e+0_fp
      ELSE
         R1 = EXP( -1e+0_fp * SQRT( ST ) )
      ENDIF

      ! surface resistance for particle
      IF (LUC == 14 .and. IS_SEASALT) THEN
         ! Use the formulation of Slinn and Slinn (1980) for the impaction over
         ! water surfaces for sea salt
         RS   = 1.e+0_fp / (USTAR**2.e+0_fp/ (W10*VON_KARMAN) * &
            (EB + EIM ) + VTS)
      ELSE
         RS   = 1.e0_fp / (E0 * USTAR * (EB + EIM + EIN) * R1 )
      ENDIF

   END FUNCTION AERO_SFCRSII


   !>
   !! \brief updates the diameter and density of non-dust aerosols
   !!
   !!References:
   !! Adapted from GEOS-Chem source code (GeosCore/drydep_mod.F90)
   !! ADUST_SFCRSII and AERO_SFCRSII functions
   !!
   !! \param SPC        Species name
   !! \param IS_SEASALT Is seasalt species?
   !! \param RHBL       Relative humidity [unitless]
   !! \param RDRY       Dry radius of particle [m]
   !! \param RWET       Wet radius of particle [m]
   !! \param DIAM       diameter of wet particle [m]
   !! \param DEN        density of particle [kg/m3]
   !! \param RC         return code
   !!
   !! \ingroup catchem_drydep_process
   !!!>
   SUBROUTINE New_DIAM_DEN( SPC, IS_SEASALT, RHBL, RDRY, RWET, DIAM, DEN, RC)

      IMPLICIT NONE

      !input parameters
      character(len=20), INTENT(IN) :: SPC    ! Species name
      logical, INTENT(IN) :: IS_SEASALT  ! Is sea salt species?
      !INTEGER,  INTENT(IN) :: K    ! Drydep species index (range: 1-NUMDEP)
      real(fp), INTENT(IN)    :: RHBL    ! Relative humidity local
      real(fp), INTENT(OUT)   :: RDRY    ! dry radius of particle [m]
      real(fp), INTENT(OUT)   :: RWET    ! wet radius of particle [m]
      real(fp), INTENT(INOUT) :: DIAM    ! diameter of wet particle [m]
      real(fp), INTENT(INOUT) :: DEN     ! density of particle [kg/m3]
      integer, INTENT(OUT)    :: RC      ! success flag
      !defined parameters
      REAL(fp), PARAMETER   :: C1       =  0.7674_fp
      REAL(fp), PARAMETER   :: C2       =  3.079_fp
      REAL(fp), PARAMETER   :: C3       =  2.573e-11_fp
      REAL(fp), PARAMETER   :: C4       = -1.424_fp
      !REAL(fp), PARAMETER   :: E0       =  3.0_fp

      ! Parameters for polynomial coefficients to derive seawater
      ! density. From Tang et al. (1997)
      REAL(fp),  PARAMETER  :: A1       =  7.93e-3_fp
      REAL(fp),  PARAMETER  :: A2       = -4.28e-5_fp
      REAL(fp),  PARAMETER  :: A3       =  2.52e-6_fp
      REAL(fp),  PARAMETER  :: A4       = -2.35e-8_fp
      REAL(fp),  PARAMETER  :: EPSI     =  1.0e-4_fp

      ! parameters for assumed size distribution of accumulation and coarse
      ! mode sea salt aerosols, as described in Jaegle et al. (ACP, 11, 2011)
      ! 1) geometric dry mean diameters (microns)
      !REAL(fp),  PARAMETER  :: RG_A     =  0.085e+0_fp
      !REAL(fp),  PARAMETER  :: RG_C     =  0.4e+0_fp
      ! 2) sigma of the size distribution
      !REAL(fp),  PARAMETER  :: SIG_A    =  1.5e+0_fp
      !REAL(fp),  PARAMETER  :: SIG_C    =  1.8e+0_fp

      !increment of radius for integration of settling velocity (um)
      !REAL(fp), PARAMETER   :: DR       =  5.0e-2_fp
      !local variables
      real(fp)    :: FAC1, FAC2  !Exponential factors for hygroscopic growth
      real(fp)    :: RUM         !Radius of dry particle in micronmeters [um]
      REAL(fp)    :: RATIO_R     !Ratio dry over wet radii
      REAL(fp)    :: DEN0, DEN1, WTP
      integer     :: I          !Loop index
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at New_DIAM_DEN (in process/drydep/ccpr_drydep_common_mod.F90)'

      IF ( .NOT. IS_SEASALT ) THEN

         ! Particle diameter [m]
         DIAM  = 0.17378e-6_fp
         RDRY = DIAM / 2.0e+0_fp !Not needed for further calculations for dust species

         ! SIA
         !IF ( K == idd_NIT .or. K == idd_NH4 .or. K == idd_SO4 ) THEN
         IF ( SPC == 'NIT' .or. SPC == 'NH4' .or. SPC == 'SO4' ) THEN
            ! Efflorescence transitions
            IF (RHBL .LT. 0.35) THEN
               ! DIAM is not changed
            ELSE IF ((RHBL .GE. 0.35) .AND. (RHBL .LE. 0.40)) THEN
               ! Linear hygroscopic growth
               DIAM = DIAM + (DIAM * ((1.0_fp + 0.61_fp * 0.40_fp /             &
                  (1.0_fp - 0.40_fp)) ** (1.0_fp / 3.0_fp)) - DIAM) /        &
                  (0.40_fp - 0.35_fp) * (RHBL - 0.35_fp)
            ELSE
               ! Kohler hygroscopic growth
               DIAM = DIAM * ((1.0_fp + 0.61_fp * RHBL / (1.0_fp - RHBL))       &
                  ** (1.0_fp / 3.0_fp))
            ENDIF

            !BC
            !ELSE IF ( K == idd_BCPI .OR. K == idd_BCPO )  THEN
         ELSE IF ( SPC == 'BCPI' .OR. SPC == 'BCPO' )  THEN
            ! DIAM is not changed

            !OA
         ELSE
            DIAM = DIAM * ((1.0_fp + 0.1_fp * RHBL / (1.0_fp - RHBL))             &
               ** (1.0_fp / 3.0_fp))
         ENDIF

         !get RWET
         RWET = DIAM / 2.0e+0_fp
         ! Particle density [kg/m3]; same for all aerosols except sea salt and  dust
         DEN   = 1500

      ELSE !sea salt aerosol case

         !drydepRadius = A_RADI(K)
         RDRY = DIAM / 2.0e+0_fp

         ! Coarse seasalt
         !IF ( K == idd_NITS .or. K == idd_SALC .or. K == idd_SO4S .or. K == idd_BRSALC .or. K == idd_ISALC ) THEN
         IF ( SPC == 'NITS' .or. SPC == 'SALC' .or. SPC == 'SO4S' .or. SPC == 'BRSALC' .or. SPC == 'ISALC' ) THEN
            RDRY = 0.74025E-6
         ENDIF

         !IF ( K == idd_SALA .OR. K == idd_BRSALA .or. K == idd_ISALA ) THEN
         IF ( SPC == 'SALA' .OR. SPC == 'BRSALA' .or. SPC == 'ISALA' ) THEN
            RDRY = 0.114945E-6
         ENDIF

         ! Dry particle radius [um]
         RUM  = RDRY * 1.e+6_fp

         ! Exponential factors used for hygroscopic growth (not used now)
         FAC1 = C1 * ( RUM**C2 )
         FAC2 = C3 * ( RUM**C4 )

         ! Corrected bug in Gerber formulation: use of LOG10  (jaegle 5/11/11)
         !RWET    = 0.01e+0_fp*(FAC1/(FAC2-DLOG(RHBL))+RCM**3.e+0_fp)**0.33e+0_fp
         !RWET = 1.d-6*(FAC1/(FAC2-LOG10(RHBL))+RUM**3.e+0_fp)**0.33333e+0_fp

         ! Use equation 5 in Lewis and Schwartz (2006) for sea salt growth [m]
         ! (jaegle 5/11/11)
         RWET = RDRY * (4.e+0_fp / 3.7e+0_fp) * &
            ( (2.e+0_fp - RHBL)/(1.e+0_fp - RHBL) )**(1.e+0_fp/3.e+0_fp)

         ! Ratio dry over wet radii at the cubic power
         !RATIO_R = ( A_RADI(K) / RWET )**3.e+0_fp

         ! Diameter of the wet aerosol [m]
         DIAM  = RWET * 2.e+0_fp

         ! Density of the wet aerosol [kg/m3] (bec, 12/8/04)
         !DEN   = RATIO_R * A_DEN(K) + ( 1.e+0_fp - RATIO_R ) * 1000.e+0_fp

         ! Above density calculation is chemically unsound because it ignores chemical solvation.
         ! Iteratively solve Tang et al., 1997 equation 5 to calculate density of wet aerosol (kg/m3)
         ! Redefine RATIO_R
         RATIO_R = RDRY / RWET

         ! Assume an initial density of 1000 kg/m3
         DEN0 = DEN !assign initial DEN to DEN0
         DEN  = 1000.e+0_fp
         DEN1 = 0.e+0_fp !initialize
         i = 0 !initialize loop index
         !Note that if RH is too low, the loop will not converge and will run forever
         DO WHILE ( ABS( DEN1-DEN ) .gt. EPSI )
            ! First calculate weight percent of aerosol (kg_RH=0.8/kg_wet)
            WTP    = 100.e+0_fp * DEN0/DEN * RATIO_R**3.e+0_fp
            ! Then calculate density of wet aerosol using equation 5
            ! in Tang et al., 1997 [kg/m3]
            DEN1   = ( 0.9971e+0_fp + (A1 * WTP) + (A2 * WTP**2) + &
               (A3 * WTP**3) + (A4 * WTP**4) ) * 1000.e+0_fp

            ! Now calculate new weight percent using above density calculation
            WTP    = 100.e+0_fp * DEN0/DEN1 * RATIO_R**3.e+0_fp
            ! Now recalculate new wet density [kg/m3]
            DEN   = ( 0.9971e+0_fp + (A1 * WTP) + (A2 * WTP**2) + &
               (A3 * WTP**3) + (A4 * WTP**4) ) * 1000.e+0_fp

            ! add some protection against infinite loop
            i = i+1
            IF ( i .GT. 500 ) THEN
               errMsg = 'Error in calculating new density for sea salt aerosol due to very low RH input!'
               CALL CC_Error( errMsg, RC, thisLoc )
               RETURN
            ENDIF

         ENDDO
      ENDIF

   END SUBROUTINE New_DIAM_DEN

   !>
   !! \brief calculates the volume size distribution of sea-salt.
   !! This only has to be done once. We assume that sea-salt is the
   !! combination of a coarse mode and accumulation model log-normal
   !! distribution functions. The resulting arrays are: DMID = diameter
!  of bin and SALT_V = dV/dln(D) [in um3].
   !!
   !!References:
   !! Adapted from GEOS-Chem source code (GeosCore/drydep_mod.F90)
   !! INIT_WEIGHTSS function
   !!
   !!
   !! \param SALT_RLOW_um  lowest edge of sea salt radius [um]
   !! \param SALT_RUP_um   uppest edge of sea sakt radius [um]
   !!
   !! \ingroup catchem_drydep_process
   !!!>

   SUBROUTINE INIT_WEIGHTSS( SALT_RLOW_um, SALT_RUP_um, RC )

      IMPLICIT NONE
      !INPUT PARAMETERS:
      real(fp), INTENT(IN) :: SALT_RLOW_um ! lowest edge of sea salt radius [um]
      real(fp), INTENT(IN) :: SALT_RUP_um  ! uppest edge of sea sakt radius [um]
      INTEGER,         INTENT(INOUT) :: RC       ! Success or failure
      !LOCAL VARIABLES:
      !INTEGER             :: N
      REAL(fp)            :: DEDGE
      INTEGER             :: ID,NR
      !DEFINED PARAMETERS:
      ! increment of radius for integration of settling velocity (um)
      REAL(fp), PARAMETER :: DR    = 5.e-2_fp

      ! parameters for assumed size distribution of acc and coarse mode
      ! sea salt aerosols
      ! geometric dry mean diameters (microns)
      REAL(fp), PARAMETER :: RG_A  = 0.085e+0_fp
      REAL(fp), PARAMETER :: RG_C  = 0.4e+0_fp
      ! sigma of the size distribution
      REAL(fp), PARAMETER :: SIG_A = 1.5e+0_fp
      REAL(fp), PARAMETER :: SIG_C = 1.8e+0_fp
      ! Error handling
      !---------------
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc

      !=================================================================
      ! INIT_WEIGHTSS begins here!
      !=================================================================
      ErrMsg = ''
      ThisLoc = ' -> at INIT_WEIGHTSS (in process/drydep/ccpr_dryde_common_mod.F90)'

      ! Number of bins between the lowest bound of of the accumulation mode
      ! sea salt and the upper bound of the coarse mode sea salt.
      NR = INT((( SALT_RUP_um - SALT_RLOW_um )  / DR ) + 0.5e+0_fp )

      ALLOCATE( DMID( NR ), STAT=RC )
      IF ( RC /= CC_SUCCESS ) THEN
         ErrMsg = 'Could not allocate array DMID'
         CALL CC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      END IF
      DMID = 0e+0_fp

      ALLOCATE( SALT_V( NR ), STAT=RC )
      IF ( RC /= CC_SUCCESS ) THEN
         ErrMsg = 'Could not allocate array SALT_V'
         CALL CC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      END IF
      SALT_V = 0e+0_fp

      !=================================================================
      ! Define the volume size distribution of sea-salt. This only has
      ! to be done once. We assume that sea-salt is the combination of a
      ! coarse mode and accumulation model log-normal distribution functions
      !=================================================================

      ! Lower edge of 0th bin diameter [um]
      DEDGE=SALT_RLOW_um * 2e+0_fp

      ! Loop over diameters
      DO ID = 1, NR

         ! Diameter of mid-point in microns
         DMID(ID)  = DEDGE + ( DR )

         ! Calculate the dry volume size distribution as the sum of two
         ! log-normal size distributions. The parameters for the size
         ! distribution are based on Reid et al. and Quinn et al.
         ! The scaling factors 13. and 0.8 for acc and coarse mode aerosols
         ! are chosen to obtain a realistic distribution
         ! SALT_V (D) = dV/dln(D) [um3]
         SALT_V(ID) = PI / 6e+0_fp* (DMID(ID)**3) * (         &
            13e+0_fp*exp(-0.5*( LOG(DMID(ID))-       &
            LOG(RG_A*2e+0_fp) )**2e+0_fp/            &
            LOG(SIG_A)**2e+0_fp )           &
            /( sqrt(2e+0_fp * PI) * LOG(SIG_A) )  +  &
            0.8e+0_fp*exp(-0.5*( LOG(DMID(ID))-      &
            LOG(RG_C*2e+0_fp) )**2e+0_fp/            &
            LOG(SIG_C)**2e+0_fp)            &
            /( sqrt(2e+0_fp * PI) * LOG(SIG_C) )  )

         ! update the next edge
         DEDGE = DEDGE + DR*2e+0_fp
      ENDDO

   END SUBROUTINE INIT_WEIGHTSS

   !>
   !! \brief calculates the Ra and Rb term in the Wesely scheme
   !!
   !!References:
   !! Wesely et al., 1989
   !!
   !! \param TEMPK       Temperatue [K]
   !! \param PRESSU      Pressure [Pa]
   !! \param XMW         Molecular weight [kg/mol]
   !! \param USTAR       Fictional Velocity [m/s]
   !! \param OBK         Monin-Obhukov length [m]
   !! \param ZO          Roughness length [m]
   !! \param THIK        height of first model layer [m]
   !! \param LNLPBL      flag to use non-local mixing
   !! \param IS_GAS      flag for gas
   !! \param Ra          output of aerodynamic resistance [s/m]
   !! \param Rb          output of quasi-laminar boundary layer resistance [s/m]
   !! \param RC          Success or failure?
   !!
   !! \ingroup catchem_drydep_process
   !!!>
   subroutine Wesely_Ra_Rb(TEMP, PRESSU, XMW, USTAR, OBK, ZO, THIK, LNLPBL, IS_GAS, Ra, Rb, RC)
      IMPLICIT NONE
      ! Parameters
      !-----------
      real(fp), intent(in)  :: TEMP        !< Temperature [K]
      real(fp), intent(in)  :: PRESSU      !< Pressure [Pa]
      real(fp), intent(in)  :: XMW         !< Molecular weight [kg/mol]
      real(fp), intent(in)  :: USTAR       !< Friction velocity [m/s]
      real(fp), intent(in)  :: OBK         !< Monin-Obhukov length [m]
      real(fp), intent(in)  :: ZO          !< Roughness length [m]
      real(fp), intent(in)  :: THIK        !< height of first model layer [m]
      logical, intent(in)   :: LNLPBL      !< flag to use non-local mixing
      logical, intent(in)   :: IS_GAS      !< flag for gas
      !output
      real(fp), intent(out) :: Ra          !< aerodynamic resistance [s/m]
      real(fp), intent(out) :: Rb          !< quasi-laminar boundary layer resistance [s/m]
      integer, intent(out)  :: RC          !< Success or failure?

      ! Local Variables
      !----------------
      real(fp) :: C1,CZ,XNU
      real(fp) :: CKUSTR,REYNO,CORR1,CORR2,Z0OBK
      real(fp) :: DUMMY1,DUMMY2,DUMMY3,DUMMY4
      real(fp) :: DAIR,TEMPK,TEMPC
      logical  :: LRGERA !stable atmosphere; a high aerodynamic resistance (RA=1.E4 m s-1) is imposed; else RA is calculated
      !string
      character(len=255)       :: thisLoc
      character(len=512)       :: ErrMsg

      !--------------------------------------------
      ! main function
      !--------------------------------------------

      ! Assume success
      RC      =  CC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> at Wesely_Ra_Rb (in process/drydep/CCPr_drydep_Commmon_Mod.F90)'

      ! Zero variables that aren't zeroed below
      CZ         = 0.0_fp
      CKUSTR     = 0.0_fp
      REYNO      = 0.0_fp
      CORR1      = 0.0_fp
      CORR2      = 0.0_fp
      Z0OBK      = 0.0_fp
      DUMMY1     = 0.0_fp
      DUMMY2     = 0.0_fp
      DUMMY3     = 0.0_fp
      DUMMY4     = 0.0_fp
      DAIR       = 0.0_fp
      Ra         = 0.0_fp
      Rb         = 0.0_fp

      !CZ is Altitude (m) at which deposition velocity is computed
      !use Midpoint height of first model level [m]
      CZ = THIK / 2.0e+0_fp

      !** TEMPK and TEMPC are surface air temperatures in K and in C
      TEMPK = TEMP
      TEMPC = TEMP-273.15e+0_fp

      !** Calculate the kinematic viscosity XNU (m2 s-1) of air
      !** as a function of temperature.
      !** The kinematic viscosity is used to calculate the roughness heights
      !** over water surfaces and to diagnose whether such surfaces are
      !** aerodynamically rough or smooth using a Reynolds number criterion.
      !** The expression for the temperature dependence of XNU
      !** is from the FORTRAN code in Appendix II of Wesely [1988];
      !** I wasn't able to find an original reference but it seems benign enough.
      C1  = TEMPK/273.15e+0_fp
      XNU = 0.151e+0_fp*(C1**1.77e+0_fp)*1.0e-04_fp

      !***** Get aerodynamic resistances Ra and Rb. ***********
      !   The aerodynamic resistance Ra is integrated from altitude z0+d up
      !   to the altitude z1 at which the dry deposition velocity is to be
      !   referenced. The integration corrects for stability using Monin-
      !   Obukhov similarity formulas from Businger et al. [1971] which
      !   apply over the range -2.5 < z/zMO < 1.5 (see their Figure 2).
      !   Under very unstable conditions when z1 > -2.5 zMO, we assume that
      !   there is no resistance to transfer in the convective column
      !   between zMO and z1. Under very stable conditions when z1 > 1.5 zMO
      !   we assume that vertical transfer in the column between zMO and z1
      !   is strongly suppressed so that the deposition velocity at altitude
      !   z1 is very low.  Under these conditions we just specify a very
      !   large Ra=1.E4 s m-1 (LRGERA = T).
      !**
      !   The Reynolds number REYNO diagnoses whether a surface is
      !   aerodynamically rough (REYNO > 1) or smooth.
      !
      !   NOTE: The criterion "REYNO > 1" was originally "REYNO > 10". See
      !   below for an explanation of why it was changed (hyl, 10/15/99)
      !
      !   Surface is rough in all cases except over water with low wind
      !   speeds. In the smooth case, vertical transport IN THE SUBLAYER
      !   near the surface is limited by molecular diffusion and is
      !   therefore very slow; we assign a large value we assign a large
      !   value of Ra + Rb to account for this effect.  [In Versions 3.2
      !   and earlier we used the formulation for Ra + Rb given in Equation
      !   (12) of Walcek et al [1986] to calculate the aerodynamic
      !   resistance over smooth surfaces.  However, that expression fails
      !   when u* is very small, as it yields negative values of Ra + Rb].
      !   (djj, hyl, bmy, 5/8/00)
      !**
      !   In the aerodynamically rough case, the expression for Ra is as
      !   given in equation (5) of Jacob et al. [1992]:
      !
      !          Ra = (1/ku*)*int(from z0 to z1) (phi(x)/z)dz
      !
      !   where x = (z-D)/zMO, z is the height above ground, and D is the
      !   displacement height which is typically 70-80% of the canopy
      !   height [Brutsaert, 1982].  We change the vertical coordinate so
      !   that z=0 at the displacement height; that's OK since for all
      !   practical applications z1 >> D.  In this manner we don't need
      !   to assume any specific value for the displacement height.
      !   Applying the variable transformation z -> x = z/zMO, the equation
      !   above becomes
      !
      !          Ra = (1/ku*)*int(from x0 to x1) (phi(x)/x)dx with x=z/zMO
      !
      !   Here phi is a stability correction function originally formulated
      !   by Businger et al. [1971] and given in eqns 5a and 5b of Jacob et
      !   al. [1992]. For unstable conditions,
      !
      !          phi(x) = a/sqrt(1-bx)  where a=0.74, b = 9
      !
      !   The analytical solution to the integral is [Dwight, 1957,
      !   integral 192.11]:
      !
      !          int(dx/(x*sqrt(1-bx))) = log(abs((sqrt(1-bx)-1)
      !                                   /(sqrt(1-bx)+1)))
      !
      !   which yields the expression for Ra used in the code for
      !   unstable conditions.  For stable conditions,
      !
      !          phi(x) = a + bx        where a=0.74, b = 4.7
      !
      !   and the analytical solution to the integral is
      !
      !          int((a/x)+b)dx = a*ln(x) + bx
      !
      !   which yields the expression of Ra used in the code for stable
      !   conditions.
      !**
      !   The formulation of RB for gases is equation (12) of Walcek et al.
      !   [1986].  The parameterization for deposition of aerosols does not
      !   include an RB term so RB for aerosols is set to zero.
      !   Modify phi(x) according to the non-local mixing scheme
      !   by Holtslag and Boville [1993] ( Lin, 07/18/08 )
      !   For unstable conditions,
      !          phi(x) = a/sqrt(1-bx)  where a=1.0, b=15.0
      !
      !   For stable conditions,
      !          phi(x) = a + bx
      !              where a=1.0, b=5.0 for 0 <= x <= 1, and
      !                    a=5.0, b=1.0 for x > 1.0
      !********************************************************

      CKUSTR = VON_KARMAN * USTAR
      REYNO = USTAR*ZO/XNU
      CORR1 = CZ/OBK
      ! Define Z0OBK
      Z0OBK = ZO/OBK

      LRGERA = .FALSE.
      ! Add option for non-local PBL
      IF (.NOT. LNLPBL) THEN
         IF (CORR1 .GT. 0.e+0_fp) THEN
            IF (CORR1 .GT.  1.5e+0_fp) LRGERA = .TRUE.
         ELSEIF(CORR1 .LE. 0.e+0_fp) THEN
            IF (CORR1 .LE. -2.5e+0_fp) CORR1 = -2.5e+0_fp
            CORR2 = LOG(-CORR1)
         ENDIF
      ENDIF

      !use rae function from pecision_mod to avoid "equality comparison for real" warning
      IF ( rae(CKUSTR, 0.0e+0_fp) ) THEN
         ErrMsg = 'CKUSTR cannot be zero.'
         CALL CC_Error( ErrMsg, RC, ThisLoc )
         RETURN             ! debug
      ENDIF

      !...aerodynamically rough or smooth surface
      ! "In the classic study by Nikuradse (1933) the transition from
      ! smooth to rough was examined in pipe flow. He introduced a
      ! roughness Reynolds number Rr = U* Z0 / Nu and found the flow to
      ! be smooth for Rr < 0.13 and rough for Rr > 2.5 with a transition
      ! regime in between." (E.B. Kraus and J.A. Businger, Atmosphere-Ocean
      ! Interaction, second edition, P.144-145, 1994).
      ! Similar statements can be found in the books: Evaporation into the
      ! atmosphere, by Wilfried Brutsaert ,P.59,89, 1982; or Seinfeld &
      ! Pandis, P.858, 1998.
      ! Here we assume a sudden transition point Rr = 1 from smooth to
      ! rough, following L. Merlivat (1978, The dependence of bulk
      ! evaporation coefficients on air-water interfacial conditions as
      ! determined by the isotopic method, J. Geophys. Res., Oceans &
      ! Atmos., 83, C6, 2977-2980). Also refer to Brutsaert's book, P.125.
      ! We used to use the criterion "REYNO > 10" for aerodynamically rough
      ! surface and now change to "REYNO > 1". (hyl, 10/15/99)
      ! D. J. Jacob change the criterion for aerodynamically rough
      ! surface to REYNO > 0.1
      IF ( REYNO > 0.1e+0_fp ) THEN !rough surface
         ! Add option for non-local PBL
         !TODO: do we need to include both options?
         IF (.NOT. LNLPBL) THEN

            !...aerodynamically rough surface.
            !*
            IF (CORR1.LE.0.0e+0_fp .AND. Z0OBK .LT. -1.e+0_fp)THEN
               !*... unstable condition; set RA to zero.
               !*    (first implemented in V. 3.2)
               RA     = 0.e+0_fp
               !*... error trap: prevent CORR1 or Z0OBK from being
               !*... zero or close to zero (ckeller, 3/15/16)
            ELSEIF ( ABS(CORR1)<=SMALL .OR. ABS(Z0OBK)<=SMALL ) THEN
               RA = 0.e+0_fp
            ELSEIF (CORR1.LE.0.0e+0_fp .AND. Z0OBK .GE. -1.e+0_fp) THEN
               !*... unstable conditions;
               !*... compute Ra as described above
               DUMMY1 = (1.e+0_fp - 9e+0_fp*CORR1)**0.5e+0_fp
               DUMMY2 = (1.e+0_fp - 9e+0_fp*Z0OBK)**0.5e+0_fp
               DUMMY3 = ABS((DUMMY1 - 1.e+0_fp)/(DUMMY1 + 1.e+0_fp))
               DUMMY4 = ABS((DUMMY2 - 1.e+0_fp)/(DUMMY2 + 1.e+0_fp))
               RA = 0.74e+0_fp* (1.e+0_fp/CKUSTR) * LOG(DUMMY3/DUMMY4)

            ELSEIF((CORR1.GT.0.0e+0_fp).AND.(.NOT.LRGERA))  THEN
               !*... moderately stable conditions (z/zMO <1);
               !*... compute Ra as described above
               RA = (1e+0_fp/CKUSTR) * (.74e+0_fp*LOG(CORR1/Z0OBK) + &
                  4.7e+0_fp*(CORR1-Z0OBK))
            ELSEIF(LRGERA) THEN
               !*... very stable conditions
               RA     = 1.e+04_fp
            ENDIF
            !* check that RA is positive; if RA is negative (as occasionally
            !* happened in version 3.1) send a warning message.

         ELSE !not using non-local PBL

            IF (CORR1.LT.0.0e+0_fp) THEN
               !*... unstable conditions; compute Ra as described
               !*... above.
               !coef_a=1.e+0_fp
               !coef_b=15.e+0_fp
               DUMMY1 = (1.e+0_fp - 15.e+0_fp*CORR1)**0.5e+0_fp
               DUMMY2 = (1.e+0_fp - 15.e+0_fp*Z0OBK)**0.5e+0_fp
               DUMMY3 = ABS((DUMMY1 - 1.e+0_fp)/(DUMMY1 + 1.e+0_fp))
               DUMMY4 = ABS((DUMMY2 - 1.e+0_fp)/(DUMMY2 + 1.e+0_fp))
               RA = 1.e+0_fp * (1.e+0_fp/CKUSTR) * LOG(DUMMY3/DUMMY4)

            ELSEIF((CORR1.GE.0.0e+0_fp).AND.(CORR1.LE.1.0e+0_fp)) THEN
               !coef_a=1.e+0_fp
               !coef_b=5.e+0_fp
               RA = (1.e+0_fp/CKUSTR) * (1.e+0_fp*LOG(CORR1/Z0OBK) + &
                  5.e+0_fp*(CORR1-Z0OBK))

            ELSE ! CORR1 .GT. 1.0D0
               !coef_a=5e+0_fp
               !coef_b=1.e+0_fp
               RA = (1.e+0_fp/CKUSTR) * (5.e+0_fp*LOG(CORR1/Z0OBK) + &
                  1.e+0_fp*(CORR1-Z0OBK))
            ENDIF

            !* check that RA is positive and maximize at 1.E4 s m-1
            RA   = MIN(RA,1.e+4_fp)
            ! If RA is < 0, set RA = 0
            IF (RA .LT. 0.e+0_fp) RA = 0.0e+0_fp

         END IF !PBL or non-local PBL options

         !get Rb for a gas species; arosol Rb is set to zero
         !** DAIR is the thermal diffusivity of air; value 0.2*1.E-4 m2 s-1 cited on p. 16,476 of
         !** Jacob et al. [1992]
         DAIR = 0.2e0_fp*1.e-4_fp
         IF (IS_GAS) THEN
            RB = (2.e+0_fp/CKUSTR)* (DAIR/DIFFG(TEMPK,PRESSU,XMW)) &
               **0.667e+0_fp
         END IF

      ELSE  !smooth surface
         !** suppress drydep over smooth surfaces by setting Ra to
         !** a large value (1e4).  This prevents negative dry deposition
         !** velocities when u* is very small. Rb is not important in that case since
         !** the total resistentce is Ra + Rb. So we set Rb to zero.
         RA     = 1.e+4_fp

      END IF

   end subroutine Wesely_Ra_Rb


end module CCPr_drydep_Common_Mod
