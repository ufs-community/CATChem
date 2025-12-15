!> \file DryDepScheme_ZHANG_Mod.F90
!! \brief Zhang et al. [2001] scheme with Emerson et al. [2020] updates.
!! The Ra and Rb are still from Wesely (1989) for now.
!!
!! Pure science kernel for zhang scheme in drydep process.
!! This module contains ONLY the computational algorithm with NO infrastructure dependencies.
!! Uses only basic Fortran types for maximum portability and reusability.
!! Reference:
!! (1) Wesely, M. L. (1989). Parameterization of surface resistances to gaseous dry
!!     deposition in regional-scale numerical models. Atmospheric Environment.
!! (2) Zhang, L., Gong, S., Padro, J., & Barrie, L. (2001). A size-segregated particle
!!     dry deposition scheme for an atmospheric aerosol module. Atmospheric environment.
!! (3) Emerson, E. W., et al. (2020). Revisiting particle dry deposition and its role
!!     in radiative effect estimates. PNAS, 117(42), 26076-26082.
!! (4) Most of the codes are adopted from GEOS-Chem drydep_mod.F90 module.
!!     https://github.com/geoschem/geos-chem
!!
!! SCIENCE CUSTOMIZATION GUIDE:
!! 1. Modify the algorithm in compute_zhang (search for "TODO")
!! 2. Add scheme-specific helper subroutines as needed
!! 3. Update physical constants for your scheme
!! 4. Customize the environmental response functions
!!
!! INFRASTRUCTURE RESPONSIBILITIES (handled by host model):
!! - Parameter initialization and validation
!! - Input array validation and error handling
!! - Memory management and array allocation
!! - Integration with host model time stepping
!!
!! Generated on: 2025-11-13T17:12:59.281598
!! Author: Wei Li
!! Reference: Zhang et al., 2001; Emerson et al., 2020
module DryDepScheme_ZHANG_Mod

   use precision_mod, only: fp, rae, f8
   use error_mod, only: CC_SUCCESS, CC_Error
   use DryDepCommon_Mod, only: DryDepSchemeZHANGConfig
   use Constants, only: PI, AVO, VON_KARMAN, RSTARG, g0, BOLTZ  !load the constants needed for this scheme

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_zhang

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
   ! Noah-MP LUC mapping
   INTEGER :: LUCINDEX_NOAH(20) = (/ 1, 2, 3, 4, 5, 10, 10, 10, 10, 6, 11, 7, 15, 7, 12, 8, 14, 9, 9, 9 /)
   ! IGBP LUC mapping
   INTEGER :: LUCINDEX_IGBP(17) = (/ 1, 2, 3, 4, 5, 10, 10, 10, 10, 6, 11, 7, 15, 7, 12, 8, 14 /)

   !same as in the Wesely scheme
   integer :: IDEP_IOLSON(74)
   DATA IDEP_IOLSON  / 11,10, 5, 3, 3, 2, 2, 5, 8, 7,  5,  8,  1,  9,  11, 11, 5,   5,  5,  2, 6,  3,   3,  2,  2,  2,  2, &
      3, 6,   6,  4,  4,  2,  6,  2,  4,  9,  4,  4,  4,  5, 5,   5,  2,  5, 9,  5,   5,  2,  8,  8,  5, &
      5, 7,   2,  4,  2,  2,  2,  5,  2,  2,  3,  5,  5,  9, 9,   9,  9,  8, 8,  8,   9,  11/
   real(fp), parameter :: TWO_THIRDS  = 2.0_fp / 3.0_fp

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

   ! REAL(fp)  :: GAMMA(15) = (/ 0.56e+0_fp, 0.58e+0_fp, 0.56e+0_fp, &
   !    0.56e+0_fp, 0.56e+0_fp, 0.54e+0_fp, &
   !    0.54e+0_fp, 0.54e+0_fp, 0.54e+0_fp, &
   !    0.54e+0_fp, 0.54e+0_fp, 0.54e+0_fp, &
   !    0.50e+0_fp, 0.50e+0_fp, 0.56e+0_fp  /)

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

   ! SeaSalt bin boundaries
   real(fp), allocatable, save :: SeaSalt_Lower_Bin(:), SeaSalt_UPPER_Bin(:)
   ! Allocatable arrays for sea salt volume size bins (persistent across calls)
   REAL(fp),   ALLOCATABLE, SAVE :: DMID    (:    )
   REAL(fp),   ALLOCATABLE, SAVE :: SALT_V  (:    )
   !TODO:put sea salt size bins here for now; may be read from input file later
   !real(fp), parameter :: SALA_REDGE_um(2)=(/0.01_fp, 0.5_fp/) !< accumulation mode Sea salt radius bin [um]
   !real(fp), parameter :: SALC_REDGE_um(2)=(/0.5_fp, 8.0_fp/) !< coarse mode Sea salt radius bin [um]

contains

   !> Pure science computation for zhang scheme
   !!
   !! This is a pure computational kernel implementing Zhang et al. [2001] scheme with Emerson et al. [2020] updates.
   !! NO error checking, validation, or infrastructure concerns.
   !! Host model must ensure all inputs are valid before calling.
   !!
   !! @param[in]  num_layers     Number of vertical layers
   !! @param[in]  num_species    Number of chemical species
   !! @param[in]  params         Scheme parameters (pre-validated by host)
   !! @param[in]  bxheight    BXHEIGHT field [appropriate units]
   !! @param[in]  frlanduse    FRLANDUSE field [appropriate units]
   !! @param[in]  iland    ILAND field [appropriate units]
   !! @param[in]  isice    IsIce field [appropriate units]
   !! @param[in]  issnow    IsSnow field [appropriate units]
   !! @param[in]  lucname    LUCNAME field [appropriate units]
   !! @param[in]  obk    OBK field [appropriate units]
   !! @param[in]  ps    PS field [appropriate units]
   !! @param[in]  rh    RH field [appropriate units]
   !! @param[in]  ts    TS field [appropriate units]
   !! @param[in]  tstep    Time step [s] - retrieved from process interface
   !! @param[in]  u10m    U10M field [appropriate units]
   !! @param[in]  ustar    USTAR field [appropriate units]
   !! @param[in]  v10m    V10M field [appropriate units]
   !! @param[in]  z0    Z0 field [appropriate units]
   !! @param[in]  species_mw_g    Species mw_g property
   !! @param[in]  species_radius    Species radius property
   !! @param[in]  species_density    Species density property
   !! @param[in]  species_short_name    Species short_name property
   !! @param[in]  species_dd_hstar    Species dd_hstar property
   !! @param[in]  species_dd_DvzAerSnow    Species dd_DvzAerSnow property
   !! @param[in]  species_dd_DvzMinVal_snow    Species dd_DvzMinVal_snow property
   !! @param[in]  species_dd_DvzMinVal_land    Species dd_DvzMinVal_land property
   !! @param[in]  species_lower_radius    Species lower_radius property
   !! @param[in]  species_upper_radius    Species upper_radius property
   !! @param[in]  species_is_dust    Species is_dust property
   !! @param[in]  species_is_seasalt    Species is_seasalt property
   !! @param[in]  species_conc   Species concentrations [mol/mol] (num_layers, num_species)
   !! @param[inout] species_tendencies  Species tendency terms [mol/mol/s] (num_layers, num_species)
   !! @param[inout] drydep_con_per_species    Dry deposition concentration per species [ug/kg or ppm] (num_species)
   !! @param[inout] drydep_velocity_per_species    Dry deposition velocity [m/s] (num_species)
   !! @param[in] diagnostic_species_id Indices mapping diagnostic species to species array (optional, for per-species diagnostics)
   subroutine compute_zhang( &
      num_layers, &
      num_species, &
      params, &
      bxheight, &
      frlanduse, &
      iland, &
      isice, &
      issnow, &
      lucname, &
      obk, &
      ps, &
      rh, &
      ts, &
      tstep, &
      u10m, &
      ustar, &
      v10m, &
      z0, &
      species_mw_g, &
      species_radius, &
      species_density, &
      species_short_name, &
      species_dd_hstar, &
      species_dd_DvzAerSnow, &
      species_dd_DvzMinVal_snow, &
      species_dd_DvzMinVal_land, &
      species_lower_radius, &
      species_upper_radius, &
      species_is_dust, &
      species_is_seasalt, &
      species_conc, &
      species_tendencies, &
      is_gas, &
      drydep_con_per_species, &
      drydep_velocity_per_species, &
      diagnostic_species_id &
      )

      ! Arguments
      integer, intent(in) :: num_layers
      integer, intent(in) :: num_species
      type(DryDepSchemeZHANGConfig), intent(in) :: params
      real(fp), intent(in) :: bxheight(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: frlanduse(:)  ! Categorical field - variable dimension array
      integer, intent(in) :: iland(:)  ! Categorical field - variable dimension array
      logical, intent(in) :: isice  ! Surface field - scalar
      logical, intent(in) :: issnow  ! Surface field - scalar
      character(len=255), intent(in) :: lucname  ! Surface field - scalar
      real(fp), intent(in) :: obk  ! Surface field - scalar
      real(fp), intent(in) :: ps  ! Surface field - scalar
      real(fp), intent(in) :: rh(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: ts  ! Surface field - scalar
      real(fp), intent(in) :: tstep  ! Time step [s] - from process interface
      real(fp), intent(in) :: u10m  ! Surface field - scalar
      real(fp), intent(in) :: ustar  ! Surface field - scalar
      real(fp), intent(in) :: v10m  ! Surface field - scalar
      real(fp), intent(in) :: z0  ! Surface field - scalar
      real(fp), intent(in) :: species_mw_g(num_species)  ! Species mw_g property
      real(fp), intent(in) :: species_radius(num_species)  ! Species radius property
      real(fp), intent(in) :: species_density(num_species)  ! Species density property
      character(len=32), intent(in) :: species_short_name(num_species)  ! Species short_name property
      real(fp), intent(in) :: species_dd_hstar(num_species)  ! Species dd_hstar property
      real(fp), intent(in) :: species_dd_DvzAerSnow(num_species)  ! Species dd_DvzAerSnow property
      real(fp), intent(in) :: species_dd_DvzMinVal_snow(num_species)  ! Species dd_DvzMinVal_snow property
      real(fp), intent(in) :: species_dd_DvzMinVal_land(num_species)  ! Species dd_DvzMinVal_land property
      real(fp), intent(in) :: species_lower_radius(num_species)  ! Species lower_radius property
      real(fp), intent(in) :: species_upper_radius(num_species)  ! Species upper_radius property
      logical, intent(in) :: species_is_dust(num_species)  ! Species is_dust property
      logical, intent(in) :: species_is_seasalt(num_species)  ! Species is_seasalt property
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(inout) :: species_tendencies(num_layers, num_species)
      logical, intent(in) :: is_gas(num_species)  ! Species type flags (true=gas, false=aerosol)
      real(fp), intent(inout), optional :: drydep_con_per_species(:)
      real(fp), intent(inout), optional :: drydep_velocity_per_species(:)
      integer, intent(in), optional :: diagnostic_species_id(:)  ! Indices mapping diagnostic species to species array

      ! Local variables
      integer :: rc, k, species_idx
      integer :: diag_idx  ! For diagnostic species indexing
      real(fp) :: VD       ! Dry deposition velocity
      real(fp) :: DDFreq   ! Dry deposition frequency
      real(fp) :: C1X, RA, RB, RSURFC, VTSoutput, VK, DVZ
      real(fp) :: HSTAR, XMW, W10
      integer  :: II     !< Index of the drydep land type
      integer  :: ILDT   !< index of the land types in the grid box
      integer  :: LDT    !loop index of land types
      integer  :: LUCINDEX !mapping above II to Zhang's 15 land types for aerosols
      !logical
      logical, save :: firsttime = .true.
      !string
      character(255) :: SPC  !current species name
      character(len=255)       :: thisLoc
      character(len=512)       :: ErrMsg

      ! Note: species_tendencies and diagnostic arrays are already initialized
      ! by the host ProcessInterface before calling this subroutine.
      ! Do not re-initialize them here.

      RC = CC_SUCCESS
      ErrMsg = ''
      thisLoc = ' -> at compute_zhang (in src/process/drydep/schemes/DryDepScheme_ZHANG_Mod.F90)'

      !calculate the volume distribution of sea salt aerosols (only need to do this once)
      IF ( firsttime ) THEN
         ! Derive seasalt bin boundaries from species properties
         call get_seasalt_bin_boundaries(num_species, species_is_seasalt, &
            species_lower_radius, species_upper_radius, &
            SeaSalt_Lower_Bin, SeaSalt_UPPER_Bin)

         CALL INIT_WEIGHTSS(MINVAL(SeaSalt_Lower_Bin), MAXVAL(SeaSalt_UPPER_Bin), RC)
         IF ( RC /= CC_SUCCESS ) THEN
            ErrMsg = 'Could not Allocate arrays in INIT_WEIGHTSS'
            CALL CC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
         firsttime = .false.
      ENDIF

      ! Calculate 10m wind speed
      W10 = sqrt(U10M**2 + V10M**2)

      ! Main computation loop - CUSTOMIZE THIS SECTION FOR YOUR SCHEME
      do k = 1, num_layers

         ! Apply to each species
         do species_idx = 1, num_species
            ! Skip species that don't match scheme type (gas vs aerosol)
            if (is_gas(species_idx)) cycle

            ! Add option for non-local PBL mixing scheme: THIK must be the first box height.
            ! TODO: we only use non-local mixing here
            !IF (.NOT. LNLPBL) THIK = MAX( ZH, THIK )

            ! Zero variables that aren't zeroed below
            VD         = 0.0_fp
            DDFreq     = 0.0_fp
            DVZ        = 0.0_fp
            RSURFC     = 0.0_fp
            RA         = 0.0_fp
            RB         = 0.0_fp
            C1X        = 0.0_fp
            VK         = 0.0_fp
            VTSoutput  = 0.0_fp

            !property for current species
            HSTAR = species_dd_hstar(species_idx)
            XMW   = species_mw_g(species_idx)*1e-3_fp   !convert from g/mol to kg/mole
            SPC   = trim(species_short_name(species_idx))

            ! Better test for depositing species: We need both HSTAR and XMW
            ! to be nonzero, OR the value of AIROSOL to be true.  This should
            ! avoid any further floating point invalid issues caused by putting
            ! a zero value in a denominator.
            DO LDT =1 , SIZE(frlanduse)
               ! If the land type is not represented in grid
               ! box, then skip to the next land type
               IF ( frlanduse(LDT) <= 0 ) CYCLE

               ILDT = ILAND(LDT)
               IF ( LUCNAME == 'OLSON' ) THEN
                  ! Olson land type index + 1
                  ILDT = ILDT + 1
                  ! Dry deposition land type index
                  II   = IDEP_IOLSON(ILDT)
                  LUCINDEX = LUCINDEX_GC(II)
               ELSE IF ( LUCNAME == 'NOAH' ) THEN
                  ! it is possible that water is given as 0 not 17 in GFS CCPP
                  IF (ILDT == 0) ILDT = 17
                  !II   = IDEP_NOAH(ILDT)
                  !Note: we use ILDT, instead of II,  to get LUCINDEX here
                  LUCINDEX = LUCINDEX_NOAH(ILDT)
               ELSE IF ( LUCNAME == 'IGBP' ) THEN
                  ! it is possible that water is given as 0 not 17
                  IF (ILDT == 0) ILDT = 17
                  !II   = IDEP_IGBP(ILDT)
                  LUCINDEX = LUCINDEX_IGBP(ILDT)
               ENDIF

               !get bulk surface resistances (Rs)
               !Note to change pressure unit from Pa to kPa
               RSURFC = AERO_SFCRSII ( SPC, species_is_dust(species_idx), species_is_seasalt(species_idx), LUCINDEX, &
                  species_radius(species_idx)*1e-6_fp, species_density(species_idx), PS*1e-3_fp, & !um to m; Pa to kPa
                  TS, USTAR, RH(1), W10, SeaSalt_Lower_Bin, SeaSalt_UPPER_Bin,VTSoutput, RC)

               if (RC /= CC_SUCCESS ) then
                  errMsg = 'Error in getting bulk surface resistances (RSURFC)'
                  CALL CC_Error( errMsg, RC, thisLoc )
                  RETURN
               endif

               !*Set max and min values for bulk surface resistances
               RSURFC = MAX(1.e+0_fp, MIN(RSURFC,9999.e+0_fp))
               ! Set Rc for strong acids (HNO3,HCl,HBr) to 1 s/m
               ! Ref. Jaegle et al. 2018, cf. Erisman,van Pul,Ayers 1994
               IF ( HSTAR .gt. 1.e+10_fp ) RSURFC= 1.e+0_fp

               !get Ra and Rb
               call Wesely_Ra_Rb(TS, PS, XMW, USTAR, OBK, Z0, bxheight(1), .FALSE., Ra, Rb,  RC)

               !get VD (TODO: IUSE is decimal not percent or permille as in GEOS-Chem)
               C1X = RSURFC + Ra + Rb
               VK = VD
               !VD = VK + DBLE( IUSE(LDT) ) / C1X + DBLE( IUSE(LDT) ) * VTSoutput
               VD = VK +  frlanduse(LDT)  / C1X +  frlanduse(LDT) * VTSoutput
            END DO

            !apply spectial treatment or scaling factor to Vd
            DVZ = VD *100.e+0_fp !m/s -- > cm/s

            !-----------------------------------------------------------
            ! Special treatment for snow and ice
            !-----------------------------------------------------------
            IF ( (ISSNOW) .OR. (ISICE) ) THEN

               !-------------------------------------
               ! %%% SURFACE IS SNOW OR ICE %%%
               !-------------------------------------
               IF ( species_DD_DvzAerSnow(species_idx) > 0.0_fp ) THEN

                  ! For most aerosol species (basically everything
                  ! except sea salt and dust species), we just set
                  ! the deposition velocity over snow to a fixed value
                  !DVZ = DBLE( DD_DvzAerSnow )
                  DVZ = species_DD_DvzAerSnow(species_idx)

               ELSE

                  ! Otherwise, enforce a minimum drydep velocity over snow
                  ! (cf. the GOCART model).  NOTE: In practice this will
                  ! only apply to the species SO2, SO4, MSA, NH3, NH4, NIT.
                  !DVZ = MAX( DVZ, DBLE( DD_DvzMinVal_Snow ) )
                  DVZ = MAX( DVZ,  species_DD_DvzMinVal_Snow(species_idx) )

               ENDIF

            ELSE

               !-------------------------------------
               ! %%% SURFACE IS NOT SNOW OR ICE %%%
               !-------------------------------------

               ! Enforce a minimum drydep velocity over land (cf. the
               ! GOCART model).  NOTE: In practice this will only apply
               ! to the species SO2, SO4, MSA, NH3, NH4, NIT.
               !DVZ = MAX( DVZ, DBLE( DD_DvzMinVal_Land ) )
               DVZ = MAX( DVZ,  species_DD_DvzMinVal_Land(species_idx) )

            ENDIF

            !-----------------------------------------------------------
            ! Compute drydep velocity and frequency
            !-----------------------------------------------------------

            ! Dry deposition velocities [m/s]
            VD = DVZ / 100.e+0_fp * params%scale_factor

            ! Dry deposition frequency [1/s]
            DDFreq = VD / bxheight(1)

            ! Ensure non-negative emissions
            species_tendencies(k, species_idx) = max(0.0_fp, DDFreq)

            ! TODO: Update diagnostic fields here based on your scheme's requirements
            ! Each process should implement custom diagnostic calculations
            ! Example patterns:
            ! Per-species diagnostic: only update for diagnostic species
            if (present(drydep_con_per_species) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == species_idx) then
                     ! Add your custom dry deposition concentration per species calculation
                     drydep_con_per_species(diag_idx) =  &
                        MAX(0.0_fp, species_conc(k,species_idx) * (1.0_fp - exp(-1.0_fp * species_tendencies(k, species_idx) * tstep)))
                     exit
                  end if
               end do
            end if
            ! Per-species diagnostic: only update for diagnostic species
            if (present(drydep_velocity_per_species) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == species_idx) then
                     ! Add your custom dry deposition velocity calculation
                     drydep_velocity_per_species(diag_idx) = VD
                     exit
                  end if
               end do
            end if
         end do

      end do

   end subroutine compute_zhang

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines
   ! Examples: environmental response functions, species-specific calculations, etc.

   !> Extract and sort seasalt bin boundaries from species properties
   !!
   !! This subroutine identifies seasalt species and extracts their lower and upper
   !! radius boundaries, then sorts them in ascending order using the same approach
   !! as Find_SeaSalt_Bin in ChemState_Mod.
   !!
   !! @param[in]  num_species       Number of chemical species
   !! @param[in]  is_seasalt        Logical array indicating which species are seasalt
   !! @param[in]  lower_radius      Lower radius bounds for each species [μm]
   !! @param[in]  upper_radius      Upper radius bounds for each species [μm]
   !! @param[out] lower_bin         Sorted lower bin boundaries [μm]
   !! @param[out] upper_bin         Sorted upper bin boundaries [μm]
   subroutine get_seasalt_bin_boundaries(num_species, is_seasalt, lower_radius, upper_radius, &
      lower_bin, upper_bin)
      implicit none

      ! Arguments
      integer, intent(in) :: num_species
      logical, intent(in) :: is_seasalt(num_species)
      real(fp), intent(in) :: lower_radius(num_species), upper_radius(num_species)
      real(fp), allocatable, intent(out) :: lower_bin(:), upper_bin(:)

      ! Local variables
      integer :: n, n_bin, rc
      real(fp), allocatable :: temp_lower(:), temp_upper(:), temp_lower1(:), temp_upper1(:)
      logical, allocatable :: mask(:)
      ! Error handling
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: thisLoc

      RC = CC_SUCCESS
      ErrMsg = ''
      thisLoc = ' -> at get_seasalt_bin_boundaries (in src/process/drydep/schemes/DryDepScheme_ZHANG_Mod.F90)'

      ! Count unique seasalt bins (similar to Find_SeaSalt_Bin approach)
      n_bin = 0
      allocate(temp_lower1(num_species), temp_upper1(num_species))

      do n = 1, num_species
         if (.not. is_seasalt(n)) cycle

         if (n_bin == 0) then
            ! First seasalt species
            n_bin = 1
            temp_lower1(n_bin) = lower_radius(n)
            temp_upper1(n_bin) = upper_radius(n)
         else
            ! Check if this radius already exists
            if (ALL(ABS(temp_lower1(1:n_bin) - lower_radius(n)) > 0.0_fp)) then
               n_bin = n_bin + 1
               temp_lower1(n_bin) = lower_radius(n)
               temp_upper1(n_bin) = upper_radius(n)
            endif
         endif
      enddo

      if (n_bin == 0) then
         ! No seasalt species - allocate empty arrays
         allocate(lower_bin(0), upper_bin(0))
         deallocate(temp_lower1, temp_upper1)
         return
      end if

      ! Allocate output arrays
      allocate(lower_bin(n_bin), upper_bin(n_bin))
      allocate(temp_lower(n_bin), temp_upper(n_bin))
      allocate(mask(n_bin))

      !copy temp_lower1 and temp_upper1 to temp_lower and temp_upper for sorting
      temp_lower = temp_lower1(1:n_bin)
      temp_upper = temp_upper1(1:n_bin)

      !sort bins by radius from low to high for lower_radius
      mask(1:n_bin) = .TRUE.
      do n = 1, n_bin
         lower_bin(n) =  MINVAL(temp_lower,mask)
         mask(MINLOC(temp_lower,mask)) = .FALSE.
      enddo

      !sort bins by radius from low to high for upper_radius
      mask(1:n_bin) = .TRUE.
      do n = 1, n_bin
         upper_bin(n) =  MINVAL(temp_upper,mask)
         mask(MINLOC(temp_upper,mask)) = .FALSE.
      enddo

      !check if the bins are continuous
      do n = 1, n_bin-1
         if ( .not. rae(upper_bin(n), lower_bin(n+1)) ) then
            errMsg = 'Sea Salt Bins are not continuous'
            call CC_Error(errMsg, RC, thisLoc)
            RETURN
         endif
      enddo

      ! Clean up
      deallocate(temp_lower, temp_upper, temp_lower1, temp_upper1,mask)

   end subroutine get_seasalt_bin_boundaries

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
   !! \param LOWERBIN   lower bound of sea-salt bins
   !! \param UPPERBIN   upper bound of sea-salt bins
   !! \param VTSout     output of setttling velocity [m/s]
   !! \param RC         success flag
   !! \param RS         return value of surface resistance [s/m]
   !!
   !! \ingroup catchem_drydep_process
   !!!>

   FUNCTION AERO_SFCRSII(  SPC, IS_DUST, IS_SEASALT, LUCINDEX, A_RADI, A_DEN, &
      PRESS, TEMP, USTAR, RHB, W10, LOWERBIN, UPPERBIN, VTSout, RC) RESULT( RS )

      IMPLICIT NONE
      !INPUT PARAMETERS
      CHARACTER(len=255), INTENT(IN) :: SPC    ! Species name
      !TODO: not sure if SPC or index is better
      !INTEGER,  INTENT(IN) :: K    ! Drydep species index (range: 1-NUMDEP)
      !INTEGER,  INTENT(IN) :: II    ! Surface type index of host model (e.g., GEOS-CHEM)
      LOGICAL,  INTENT(IN) :: IS_DUST, IS_SEASALT ! Is dust or seasalt species?
      INTEGER,  INTENT(IN) :: LUCINDEX !mapping to the 15 drydep land use categories
      REAL(fp), INTENT(IN) :: A_RADI ! Aerosol radius [m]
      REAL(fp), INTENT(IN) :: A_DEN  ! Aerosol density [kg/m3]
      REAL(fp), INTENT(IN) :: PRESS ! Pressure [kPa] (1 mb = 100 Pa = 0.1 kPa)
      REAL(fp), INTENT(IN) :: TEMP  ! Temperature [K]
      REAL(fp), INTENT(IN) :: USTAR ! Friction velocity [m/s]
      REAL(fp), INTENT(IN) :: RHB   ! Relative humidity (fraction)
      REAL(fp), INTENT(IN) :: W10   ! 10m wind speed [m/s]; only need for SeaSalt over water
      REAL(fp), DIMENSION(:), INTENT(IN) :: LOWERBIN ! lower bound of sea-salt bins
      REAL(fp), DIMENSION(:), INTENT(IN) :: UPPERBIN ! upper bound of sea-salt bins
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
      INTEGER   :: ID,NR,n, n1
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
      thisLoc = ' -> at AERO_SFCRSII (in process/drydep/scheme/DryDepScheme_ZHANG_Mod.F90)'

      !=================================================================
      ! ADUST_SFCRII begins here!
      !=================================================================

      ! Annual average of A
      Aavg(:) = (A(:,1)+A(:,2)+A(:,3)+A(:,4)+A(:,5))/5.0_fp

      LUC     = LUCINDEX
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
      !TODO: here we comment out these hardcoded lines and use the radius from inputs directly
      !IF ( K == idd_DST1 .or. K == idd_DSTAL1 .or. K == idd_NITD1 .or. K == idd_SO4D1 ) THEN
      !IF ( SPC == 'DST1' .or. SPC == 'DSTAL1' .or. SPC == 'NITD1' .or. SPC == 'SO4D1' .or. SPC == 'dust1') THEN
      !   DIAM = 0.66895E-6
      !ENDIF

      !IF ( K == idd_DST2 .or. K == idd_DSTAL2 .or. K == idd_NITD2 .or. K == idd_SO4D2 ) THEN
      !IF ( SPC == 'DST2' .or. SPC == 'DSTAL2' .or. SPC == 'NITD2' .or. SPC == 'SO4D2' .or. SPC == 'dust2') THEN
      !   DIAM = 2.4907E-6
      !ENDIF

      !IF ( K == idd_DST3  .or. K == idd_DSTAL3 .or. K == idd_NITD3 .or. K == idd_SO4D3 ) THEN
      !IF ( SPC == 'DST3'  .or. SPC == 'DSTAL3' .or. SPC == 'NITD3' .or. SPC == 'SO4D3' ) THEN
      !   DIAM = 4.164E-6
      !ENDIF

      !IF ( K == idd_DST4  .or. K == idd_DSTAL4 .or. K == idd_NITD4 .or. K == idd_SO4D4 ) THEN
      !IF ( SPC == 'DST4'  .or. SPC == 'DSTAL4' .or. SPC == 'NITD4' .or. SPC == 'SO4D4' ) THEN
      !   DIAM = 6.677E-6
      !ENDIF

      ! Hygroscopic growth following Latimer and Martin (2019) ACP
      RHBL    = MAX( TINY(RHB), RHB )

      ! Over oceans the RH in the viscous sublayer is set to 98%,
      ! following Lewis and Schwartz (2004)
      !I added condition when RHBL=1 to avoid DIAM = infinity issue in the New_DIAM_DEN subroutine later for SO4 (Wei Li)
      IF (LUC == 14 .or. rae(RHBL, 1.0_fp)) THEN
         RHBL = 0.98_fp
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

      ! Dp [m] --> [um] = particle diameter if necessary
      IF (DIAM > 0.001) THEN !here use 0.001 to determine if the unit is in m or um
         DP    = DIAM
         DIAM  = DIAM * 1.e-6_fp
      ELSE
         DP    = DIAM * 1.e+6_fp
      ENDIF

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
         NR = INT((( MAXVAL(UPPERBIN) - MINVAL(LOWERBIN) ) / DR ) + 0.5e+0_fp )

         SALT_MASS_TOTAL = 0e+0_fp
         VTS_WEIGHT      = 0e+0_fp

         ! Dry particle radius [m] --> [um]
         RUM  = RDRY * 1.e+6_fp

         ! Check what the min/max range of the SS size bins are
         !IF ( RUM .le. SALA_REDGE_um(2) ) THEN
         !   D0 = SALA_REDGE_um(1)*2e+0_fp
         !   D1 = SALA_REDGE_um(2)*2e+0_fp
         !ELSE
         !   D0 = SALC_REDGE_um(1)*2e+0_fp
         !   D1 = SALC_REDGE_um(2)*2e+0_fp
         !ENDIF
         D0 = 0e+0_fp; D1 = 0e+0_fp; n1=1
         DO n =1, size(UPPERBIN)
            IF ( (RUM .ge. LOWERBIN(n)) .and. (RUM .le. UPPERBIN(n)) ) THEN
               D0 = LOWERBIN(n)*2e+0_fp
               D1 = UPPERBIN(n)*2e+0_fp
               n1=0
               EXIT
            ENDIF
         ENDDO
         IF (n1 > 0) THEN ! D0 and D1 may not be set properly
            errMsg = 'Sea salt radius is not in any bins. Check the species namelist.'
            CALL CC_Error( errMsg, RC, thisLoc )
            RETURN
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
         R1 = MAX( tiny(R1), R1 ) !avoid R1 = 0 when ST is large under very low TEMP and AA < 0 (Wei Li)
      ENDIF

      !add error check here to make sure RS below is not a infinite value
      IF (rae(R1, 0.0_fp) .or. rae(USTAR, 0.0_fp)) THEN
         !write(*,*) 'DEBUG INFO: SPC=', trim(SPC), LUC, USTAR, R1, ST, AA, VTS, CONST, DEN, DIAM, RHBL, RHB, AIRVS
         errMsg = 'USTAR or R1 is zero. Check met field or diameter (in m) of aerosol is too big.'
         CALL CC_Error( errMsg, RC, thisLoc )
         RETURN
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
      character(len=255), INTENT(IN) :: SPC    ! Species name
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
      REAL(f8),  PARAMETER  :: EPSI     =  1.0e-4_f8 !!!Note we changed it fp to f8 otherwise the do while loop may not converge

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
      REAL(f8)    :: RATIO_R     !Ratio dry over wet radii
      REAL(f8)    :: DEN0, DEN1, WTP, DEN_f8
      integer     :: I          !Loop index
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = ' -> at New_DIAM_DEN (in process/drydep/scheme/DryDepScheme_ZHANG_Mod.F90)'

      IF ( .NOT. IS_SEASALT ) THEN

         ! Particle diameter [m]
         DIAM  = 0.17378e-6_fp
         RDRY = DIAM / 2.0e+0_fp !Not needed for further calculations for dust species

         ! SIA (TODO: keep this for now and need to be consistent with real species names in the future)
         !IF ( K == idd_NIT .or. K == idd_NH4 .or. K == idd_SO4 ) THEN
         IF ( SPC == 'NIT' .or. SPC == 'NH4' .or. SPC == 'SO4' .or. SPC == 'ASO4J' .or. &
            SPC == 'nit' .or. SPC == 'nh4' .or. SPC == 'so4' .or. SPC == 'aso4j' ) THEN
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
         ELSE IF ( SPC == 'BCPI' .OR. SPC == 'BCPO' .OR. &
            SPC == 'bcpi' .OR. SPC == 'bcpo' )  THEN
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

         ! Coarse seasalt (TODO: we commented out these hardcoded lines and use the radius from inputs directly)
         !IF ( K == idd_NITS .or. K == idd_SALC .or. K == idd_SO4S .or. K == idd_BRSALC .or. K == idd_ISALC ) THEN
         !IF ( SPC == 'NITS' .or. SPC == 'SALC' .or. SPC == 'SO4S' .or. SPC == 'BRSALC' .or. SPC == 'ISALC' ) THEN
         !   RDRY = 0.74025E-6
         !ENDIF

         !IF ( K == idd_SALA .OR. K == idd_BRSALA .or. K == idd_ISALA ) THEN
         !IF ( SPC == 'SALA' .OR. SPC == 'BRSALA' .or. SPC == 'ISALA' ) THEN
         !   RDRY = 0.114945E-6
         !ENDIF

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
         RATIO_R = real(RDRY / RWET, f8)

         ! Assume an initial density of 1000 kg/m3
         DEN0 = real(DEN, f8) !assign initial DEN to DEN0
         DEN_f8  = 1000.e+0_f8
         DEN1 = 0.e+0_f8 !initialize
         i = 0 !initialize loop index
         !Note that if RH is too low, the loop will not converge and will run forever
         DO WHILE ( ABS( DEN1-DEN_f8 ) .gt. EPSI )
            ! First calculate weight percent of aerosol (kg_RH=0.8/kg_wet)
            WTP    = 100.e+0_f8 * DEN0/DEN_f8 * RATIO_R**3
            ! Then calculate density of wet aerosol using equation 5
            ! in Tang et al., 1997 [kg/m3]
            DEN1   = ( 0.9971e+0_f8 + (A1 * WTP) + (A2 * WTP**2) + &
               (A3 * WTP**3) + (A4 * WTP**4) ) * 1000.e+0_f8

            ! Now calculate new weight percent using above density calculation
            WTP    = 100.e+0_f8 * DEN0/DEN1 * RATIO_R**3
            ! Now recalculate new wet density [kg/m3]
            DEN_f8   = ( 0.9971e+0_f8 + (A1 * WTP) + (A2 * WTP**2) + &
               (A3 * WTP**3) + (A4 * WTP**4) ) * 1000.e+0_f8

            ! add some protection against infinite loop
            i = i+1
            IF ( i .GT. 500 ) THEN
               !write(*,*) 'Test NEW_DIAM_DEN output: ', trim(SPC), RHBL, RDRY, RWET, DIAM, DEN0, DEN,DEN1
               errMsg = 'Error in calculating new density for sea salt aerosol due to very low RH input!'
               CALL CC_Error( errMsg, RC, thisLoc )
               RETURN
            ENDIF

         ENDDO
         ! Convert back to fp
         DEN = REAL(DEN_f8, fp)
      ENDIF

   END SUBROUTINE New_DIAM_DEN

   !>
   !! \brief calculates the volume size distribution of sea-salt.
   !! This only has to be done once. We assume that sea-salt is the
   !! combination of a coarse mode and accumulation model log-normal
   !! distribution functions. The resulting arrays are: DMID = diameter
   !! of bin and SALT_V = dV/dln(D) [in um3].
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
      ThisLoc = ' -> at INIT_WEIGHTSS (in process/drydep/schemes/DryDepScheme_ZHANG_Mod.F90)'

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
            13e+0_fp*exp(-0.5_fp*( LOG(DMID(ID))-       &
            LOG(RG_A*2e+0_fp) )**2e+0_fp/            &
            LOG(SIG_A)**2e+0_fp )           &
            /( sqrt(2e+0_fp * PI) * LOG(SIG_A) )  +  &
            0.8e+0_fp*exp(-0.5_fp*( LOG(DMID(ID))-      &
            LOG(RG_C*2e+0_fp) )**2e+0_fp/            &
            LOG(SIG_C)**2e+0_fp)            &
            /( sqrt(2e+0_fp * PI) * LOG(SIG_C) )  )

         ! update the next edge
         DEDGE = DEDGE + DR*2e+0_fp
      ENDDO

   END SUBROUTINE INIT_WEIGHTSS

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
   !! \brief calculates the Ra and Rb term in the Wesely scheme
   !!
   !!References:
   !! Wesely, M. L. "Parameterization of surface resistances to gaseous dry deposition in
   !! regional-scale numerical models." Atmospheric environment 41 (2007): 52-63.
   !! https://doi.org/10.1016/0004-6981(89)90153-4
   !!
   !! \param TEMPK       Temperatue [K]
   !! \param PRESSU      Pressure [Pa]
   !! \param XMW         Molecular weight [kg/mol]
   !! \param USTAR       Fictional Velocity [m/s]
   !! \param OBK         Monin-Obhukov length [m]
   !! \param ZO          Roughness length [m]
   !! \param THIK        height of first model layer [m]
   !! \param IS_GAS      flag for gas
   !! \param Ra          output of aerodynamic resistance [s/m]
   !! \param Rb          output of quasi-laminar boundary layer resistance [s/m]
   !! \param RC          Success or failure?
   !!
   !! \ingroup catchem_drydep_process
   !!!>
   subroutine Wesely_Ra_Rb(TEMP, PRESSU, XMW, USTAR, OBK, ZO, THIK, IS_GAS, Ra, Rb, RC)
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
      ThisLoc = ' -> at Wesely_Ra_Rb (in process/drydep/scheme/DryDepScheme_ZHANG_Mod.F90)'

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
      !IF (.NOT. LNLPBL) THEN
      !   IF (CORR1 .GT. 0.e+0_fp) THEN
      !      IF (CORR1 .GT.  1.5e+0_fp) LRGERA = .TRUE.
      !   ELSEIF(CORR1 .LE. 0.e+0_fp) THEN
      !      IF (CORR1 .LE. -2.5e+0_fp) CORR1 = -2.5e+0_fp
      !      CORR2 = LOG(-CORR1)
      !   ENDIF
      !ENDIF

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
      IF ( REYNO >= 0.1e+0_fp ) THEN !rough surface
         ! Add option for non-local PBL
         !TODO: we only use non-local option
         !IF (.NOT. LNLPBL) THEN

         !...aerodynamically rough surface.
         !*
         !   IF (CORR1.LE.0.0e+0_fp .AND. Z0OBK .LT. -1.e+0_fp)THEN
         !*... unstable condition; set RA to zero.
         !*    (first implemented in V. 3.2)
         !      RA     = 0.e+0_fp
         !*... error trap: prevent CORR1 or Z0OBK from being
         !*... zero or close to zero (ckeller, 3/15/16)
         !   ELSEIF ( ABS(CORR1)<=SMALL .OR. ABS(Z0OBK)<=SMALL ) THEN
         !      RA = 0.e+0_fp
         !   ELSEIF (CORR1.LE.0.0e+0_fp .AND. Z0OBK .GE. -1.e+0_fp) THEN
         !*... unstable conditions;
         !*... compute Ra as described above
         !      DUMMY1 = (1.e+0_fp - 9e+0_fp*CORR1)**0.5e+0_fp
         !      DUMMY2 = (1.e+0_fp - 9e+0_fp*Z0OBK)**0.5e+0_fp
         !      DUMMY3 = ABS((DUMMY1 - 1.e+0_fp)/(DUMMY1 + 1.e+0_fp))
         !      DUMMY4 = ABS((DUMMY2 - 1.e+0_fp)/(DUMMY2 + 1.e+0_fp))
         !      RA = 0.74e+0_fp* (1.e+0_fp/CKUSTR) * LOG(DUMMY3/DUMMY4)

         !   ELSEIF((CORR1.GT.0.0e+0_fp).AND.(.NOT.LRGERA))  THEN
         !*... moderately stable conditions (z/zMO <1);
         !*... compute Ra as described above
         !      RA = (1e+0_fp/CKUSTR) * (.74e+0_fp*LOG(CORR1/Z0OBK) + &
         !         4.7e+0_fp*(CORR1-Z0OBK))
         !   ELSEIF(LRGERA) THEN
         !*... very stable conditions
         !      RA     = 1.e+04_fp
         !   ENDIF
         !* check that RA is positive; if RA is negative (as occasionally
         !* happened in version 3.1) send a warning message.

         !ELSE !not using non-local PBL

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

         !END IF !PBL or non-local PBL options

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



end module DryDepScheme_ZHANG_Mod
