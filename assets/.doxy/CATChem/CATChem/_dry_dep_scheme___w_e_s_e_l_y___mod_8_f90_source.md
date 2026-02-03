

# File DryDepScheme\_WESELY\_Mod.F90

[**File List**](files.md) **>** [**drydep**](dir_57fb5aa14ddb2cd518a6d90b65ffd000.md) **>** [**schemes**](dir_5a3c86e36f17958630366ebc2b7ca21b.md) **>** [**DryDepScheme\_WESELY\_Mod.F90**](_dry_dep_scheme___w_e_s_e_l_y___mod_8_f90.md)

[Go to the documentation of this file](_dry_dep_scheme___w_e_s_e_l_y___mod_8_f90.md)


```Fortran

module drydepscheme_wesely_mod

   use precision_mod, only: fp, rae
   use error_mod, only: cc_success, cc_error
   use drydepcommon_mod, only: drydepschemeweselyconfig
   use constants, only: pi, h2omw, avo, von_karman, rstarg  !load the constants needed for this scheme

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_wesely

   ! module variables (mainly some constants dependent on land use in the scheme)
   integer,  parameter :: NDRYDTYPE   = 11
   real(fp), parameter :: TWO_THIRDS  = 2.0_fp / 3.0_fp
   !real(fp), parameter :: H2OMW = 18.0_fp !declared in constant module
   real(fp), parameter :: SMALL = 1.0e-10_fp 
   integer,  parameter :: IWATER = 1
   ! Arrays that hold information for each of the 11 drydep land types
   integer :: IDRYDTYPE(NDRYDTYPE)
   real(fp) :: IRAC(NDRYDTYPE),  IRCLO(NDRYDTYPE), IRCLS(NDRYDTYPE)
   real(fp) :: IRGSS(NDRYDTYPE), IRGSO(NDRYDTYPE), IRLU(NDRYDTYPE)
   real(fp) :: IRI(NDRYDTYPE),   IVSMAX(NDRYDTYPE)
   !some Olson land use (74 types) related parameters
   real(fp):: DRYCOEFF(20)
   integer :: IOLSON (74), IDEP_IOLSON(74)
   integer :: IDEP_NOAH(20), IDEP_IGBP(17)
   !integer :: IZO(74) !< Roughness height for each Olson land types

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
   DATA iri    /9999,    200,       200,          200,        200,       200,      200,   9999,   200,   9999,  9999/
   DATA irlu   /9999,   9000,      9000,         9000,       9000,      1000,     4000,   9999,  9000,   9999,  9999/
   DATA irac   /   0,   2000,      2000,          200,        100,      2000,        0,      0,   300,    100,     0/
   DATA irgss  / 100,    500,       500,          150,        350,       200,      340,   1000,     0,    400,     0/
   DATA irgso  /3500,    200,       200,          150,        200,       200,      340,    400,  1000,    300,  2000/
   DATA ircls  /9999,   2000,      2000,         2000,       2000,      9999,     9999,   9999,  2500,   9999,  9999/
   DATA irclo  /1000,   1000,      1000,         1000,       1000,      9999,     9999,   9999,  1000,   9999,  9999/
   DATA ivsmax /100,    100,       100,          100,        100,       100,      100,     10,   100,    100,    10/
   DATA idrydtype /1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11/
   ! Olson land use related parameters (https://wiki.seas.harvard.edu/geos-chem/index.php/Olson_land_map)
   DATA iolson /1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, &
      28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, &
      53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74/
   DATA idep_iolson  / 11,10, 5, 3, 3, 2, 2, 5, 8, 7,  5,  8,  1,  9,  11, 11, 5,   5,  5,  2, 6,  3,   3,  2,  2,  2,  2, &
      3, 6,   6,  4,  4,  2,  6,  2,  4,  9,  4,  4,  4,  5, 5,   5,  2,  5, 9,  5,   5,  2,  8,  8,  5, &
      5, 7,   2,  4,  2,  2,  2,  5,  2,  2,  3,  5,  5,  9, 9,   9,  9,  8, 8,  8,   9,  11/
   DATA idep_noah / 3, 6, 3, 2, 2, 5, 5, 2, 5, 5, 9, 4, 10, 5, 1, 8, 11, 7, 7, 7 /
   DATA idep_igbp / 3, 6, 3, 2, 2, 5, 5, 2, 5, 5, 9, 4, 10, 5, 1, 8, 11 / !same as NOAH without the last three tundra types
   !roughness height is not used now and is read from MET directly; so comment out
   !DATA IZO  /  10,  25000, 100,  10000, 10000, 10000, 10000, 100,  10, 2000, 100, 10,  1,   100,  1000,  1000,  1000, 100, 100, 2000, &
   !   10000,  10000,10000, 10000, 10000, 10000, 10000,10000,1000,10000,1000,1000,2000,10000,10000, 1000,  100, 1000, 1000,1000, &
   !   100,    100,   100, 2000,  100,    100,  1000, 1000, 1000, 1000,1000, 50,  50,  50,  2000,  2000, 2000, 2000, 1000, 100, &
   !   2000,   2000, 10000, 2000,  1000,  1000,  1000, 1000, 1000, 10,  1000,1000,500, 100 /
   ! Baldocchi polynomial coeffs
   DATA drycoeff /-0.358_fp, 3.02_fp,  3.85_fp, -0.0978_fp,  -3.66_fp,   12_fp,   0.252_fp,  -7.8_fp,  0.226_fp,  0.274_fp,  &
      1.14_fp,  -2.19_fp,  0.261_fp, -4.62_fp,   0.685_fp, -0.254_fp, 4.37_fp,  -0.266_fp, -0.159_fp, -0.206_fp  /

contains

   subroutine compute_wesely( &
      num_layers, &
      num_species, &
      params, &
      bxheight, &
      cldfrc, &
      frlai, &
      frlanduse, &
      iland, &
      isice, &
      island, &
      issnow, &
      lat, &
      lon, &
      lucname, &
      obk, &
      ps, &
      salinity, &
      suncosmid, &
      swgdn, &
      ts, &
      tskin, &
      tstep, &
      ustar, &
      z0, &
      species_mw_g, &
      species_dd_f0, &
      species_short_name, &
      species_dd_hstar, &
      species_dd_DvzAerSnow, &
      species_dd_DvzMinVal_snow, &
      species_dd_DvzMinVal_land, &
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
      type(DryDepSchemeWESELYConfig), intent(in) :: params
      real(fp), intent(in) :: bxheight(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: cldfrc  ! Surface field - scalar
      real(fp), intent(in) :: frlai(:)  ! Categorical field - variable dimension array
      real(fp), intent(in) :: frlanduse(:)  ! Categorical field - variable dimension array
      integer, intent(in) :: iland(:)  ! Categorical field - variable dimension array
      logical, intent(in) :: isice  ! Surface field - scalar
      logical, intent(in) :: island  ! Surface field - scalar
      logical, intent(in) :: issnow  ! Surface field - scalar
      real(fp), intent(in) :: lat  ! Surface field - scalar
      real(fp), intent(in) :: lon  ! Surface field - scalar
      character(len=255), intent(in) :: lucname  ! Surface field - scalar
      real(fp), intent(in) :: obk  ! Surface field - scalar
      real(fp), intent(in) :: ps  ! Surface field - scalar
      real(fp), intent(in) :: salinity  ! Surface field - scalar
      real(fp), intent(in) :: suncosmid  ! Surface field - scalar
      real(fp), intent(in) :: swgdn  ! Surface field - scalar
      real(fp), intent(in) :: ts  ! Surface field - scalar
      real(fp), intent(in) :: tskin  ! Surface field - scalar
      real(fp), intent(in) :: tstep  ! Time step [s] - from process interface
      real(fp), intent(in) :: ustar  ! Surface field - scalar
      real(fp), intent(in) :: z0  ! Surface field - scalar
      real(fp), intent(in) :: species_mw_g(num_species)  ! Species mw_g property
      real(fp), intent(in) :: species_dd_f0(num_species)  ! Species dd_f0 property
      character(len=32), intent(in) :: species_short_name(num_species)  ! Species short_name property
      real(fp), intent(in) :: species_dd_hstar(num_species)  ! Species dd_hstar property
      real(fp), intent(in) :: species_dd_DvzAerSnow(num_species)  ! Species dd_DvzAerSnow property
      real(fp), intent(in) :: species_dd_DvzMinVal_snow(num_species)  ! Species dd_DvzMinVal_snow property
      real(fp), intent(in) :: species_dd_DvzMinVal_land(num_species)  ! Species dd_DvzMinVal_land property
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
      real(fp) :: XLAI_IN, C1X, RA, RB, RSURFC, VK, DVZ
      real(fp) :: HSTAR, XMW, IODIDE, F0
      integer  :: II
      integer  :: ILDT
      integer  :: LDT    !loop index of land types
      character(255) :: SPC  !current species name
      character(len=255)       :: thisLoc
      character(len=512)       :: ErrMsg

      ! Note: species_tendencies and diagnostic arrays are already initialized
      ! by the host ProcessInterface before calling this subroutine.
      ! Do not re-initialize them here.
      ! Assume success
      rc      =  cc_success
      errmsg  = ''
      thisloc = ' -> at compute_Wesely (in process/drydep/scheme/DryDepScheme_WESELY_Mod.F90)'

      !TODO : Placeholder for iodide concentration; zero means O3 deposition to ocean through halogen chemistry is not
      !doing anything although the codes are there. If we have iodide as a species in the future, we can get its
      !concentration from the species_conc array here.
      iodide = 0.0_fp

      ! Main computation loop - CUSTOMIZE THIS SECTION FOR YOUR SCHEME
      do k = 1, num_layers

         ! Apply to each species
         do species_idx = 1, num_species
            ! Skip species that don't match scheme type (gas vs aerosol)
            if (.not. is_gas(species_idx)) cycle
            ! Add option for non-local PBL mixing scheme: THIK must be the first box height.
            ! TODO: we only use non-local mixing here
            !IF (.NOT. LNLPBL) THIK = MAX( ZH, THIK )

            ! Zero variables that aren't zeroed below
            vd         = 0.0_fp
            ddfreq     = 0.0_fp
            dvz        = 0.0_fp
            rsurfc     = 0.0_fp
            ra         = 0.0_fp
            rb         = 0.0_fp
            c1x        = 0.0_fp
            vk         = 0.0_fp
            xlai_in    = 0.0_fp

            !property for current species
            hstar = species_dd_hstar(species_idx)
            xmw   = species_mw_g(species_idx)*1e-3_fp   !convert from g/mol to kg/mole
            spc   = species_short_name(species_idx)
            f0    = species_dd_f0(species_idx)

            ! Better test for depositing species: We need both HSTAR and XMW
            ! to be nonzero, OR the value of AIROSOL to be true.  This should
            ! avoid any further floating point invalid issues caused by putting
            ! a zero value in a denominator.
            IF ( ( hstar > 0e+0_fp .and. xmw > 0e+0_fp ) ) THEN
               DO ldt =1 , SIZE(frlanduse)
                  ! If the land type is not represented in grid
                  ! box, then skip to the next land type
                  IF ( frlanduse(ldt) <= 0 ) cycle

                  ildt = iland(ldt)
                  IF ( lucname == 'OLSON' ) THEN
                     ! Olson land type index + 1
                     ildt = ildt + 1
                     ! Dry deposition land type index
                     ii   = idep_iolson(ildt)
                  ELSE IF ( lucname == 'NOAH' ) THEN
                     ! it is possible that water is given as 0 not 17 in GFS CCPP
                     IF (ildt == 0) ildt = 17
                     ii   = idep_noah(ildt)
                  ELSE IF ( lucname == 'IGBP' ) THEN
                     ! it is possible that water is given as 0 not 17
                     IF (ildt == 0) ildt = 17
                     ii   = idep_igbp(ildt)
                  ENDIF

                  !LAI of the landtype in the subgrid
                  !XLAI_IN = XLAI * DBLE(IUSE(LDT)) !TODO: may be able to calculate online if fraction LAI is not provided
                  xlai_in = frlai(ldt)

                  !If the surface to be snow or ice;set II to 1 instead
                  !We do not use II index to specify directly since IS_SNOW and IS_ICE are given at each grid not subgrid as ILAND
                  IF( (issnow) .OR. (isice) ) ii=1

                  !get bulk surface resistances (Rs)
                  call wesely_rc_gas( swgdn, ts, suncosmid,  f0, hstar, xmw, ustar, cldfrc, &
                     ps, xlai_in, ii,  spc, salinity, tskin, iodide, lon, lat, &
                     params%co2_effect, params%co2_level, params%co2_reference, rsurfc,   rc)

                  if (rc /= cc_success ) then
                     errmsg = 'Error in getting bulk surface resistances (RSURFC)'
                     CALL cc_error( errmsg, rc, thisloc )
                     RETURN
                  endif

                  !*Set max and min values for bulk surface resistances
                  rsurfc = max(1.e+0_fp, min(rsurfc,9999.e+0_fp))
                  !*because of high resistance values, different rule applied for ocean ozone
                  IF ( (spc .EQ. 'O3' .OR. spc .EQ. 'o3') &
                     .AND. (ii .EQ. 11)) THEN
                     rsurfc = max(1.e+0_fp, min(rsurfc,999999.e+0_fp))
                  ENDIF
                  ! Set Rc for strong acids (HNO3,HCl,HBr) to 1 s/m
                  ! Ref. Jaegle et al. 2018, cf. Erisman,van Pul,Ayers 1994
                  IF ( hstar .gt. 1.e+10_fp ) rsurfc= 1.e+0_fp

                  !get Ra and Rb
                  call wesely_ra_rb(ts, ps, xmw, ustar, obk, z0, bxheight(1), .true., ra, rb,  rc)

                  !get VD (TODO: IUSE is decimal not percent or permille as in GEOS-Chem)
                  c1x = rsurfc + ra + rb
                  vk = vd
                  !VD = VK + DBLE( IUSE(LDT) ) / C1X !This seems to be useless in the original codes

                  !VD = VK + DBLE( IUSE(LDT) ) / C1X
                  vd = vk + frlanduse(ldt)  / c1x

               END DO
            ENDIF

            !apply spectial treatment or scaling factor to Vd
            dvz = vd *100.e+0_fp !m/s -- > cm/s

            ! Scale relative to specified species(Note:we do not use FLAG but match names instead)
            !TODO: We simply hardcode the scaling factor here

            !IF ( FLAG(D) .eq. 1 )  THEN
            IF ((spc .eq. 'N2O5') .or. (spc .eq. 'HC187') .or. &
               (spc .eq. 'n2o5') .or. (spc .eq. 'hc187') ) THEN

               ! Scale species to HNO3 (MW_g = 63.012 g/mol)
               dvz = dvz * sqrt(63.01_fp) / sqrt( xmw*1e3_fp )

               !ELSE IF ( FLAG(D) .eq. 2 ) THEN
            ELSE IF ((spc .eq. 'MPAN') .or. (spc .eq. 'PPN') .or. (spc .eq. 'R4N2') .or. &
               (spc .eq. 'mpan') .or. (spc .eq. 'ppn') .or. (spc .eq. 'r4n2') ) THEN

               ! Scale species to PAN (MW_g = 121.06 g/mol)
               dvz = dvz * sqrt(121.06_fp) / sqrt( xmw*1e3_fp )

               !ELSE IF ( FLAG(D) .eq. 3 ) THEN
            ELSE IF ((spc .eq. 'MONITS') .or. (spc .eq. 'MONITU') .or. (spc .eq. 'HONIT') .or. &
               (spc .eq. 'monits') .or. (spc .eq. 'monitu') .or. (spc .eq. 'honit') ) THEN

               ! Scale species to ISOPN (MW_g = 147.15 g/mol)
               dvz = dvz * sqrt(147.15_fp)  / sqrt(xmw*1e3_fp)

            ENDIF

            !-----------------------------------------------------------
            ! Special treatment for snow and ice
            !-----------------------------------------------------------
            IF ( (issnow) .OR. (isice) ) THEN

               !-------------------------------------
               ! %%% SURFACE IS SNOW OR ICE %%%
               !-------------------------------------
               IF ( species_dd_dvzaersnow(species_idx) > 0.0_fp ) THEN

                  ! For most aerosol species (basically everything
                  ! except sea salt and dust species), we just set
                  ! the deposition velocity over snow to a fixed value
                  !DVZ = DBLE( DD_DvzAerSnow )
                  dvz = species_dd_dvzaersnow(species_idx)

               ELSE

                  ! Otherwise, enforce a minimum drydep velocity over snow
                  ! (cf. the GOCART model).  NOTE: In practice this will
                  ! only apply to the species SO2, SO4, MSA, NH3, NH4, NIT.
                  !DVZ = MAX( DVZ, DBLE( DD_DvzMinVal_Snow ) )
                  dvz = max( dvz,  species_dd_dvzminval_snow(species_idx) )

               ENDIF

            ELSE

               !-------------------------------------
               ! %%% SURFACE IS NOT SNOW OR ICE %%%
               !-------------------------------------

               ! Enforce a minimum drydep velocity over land (cf. the
               ! GOCART model).  NOTE: In practice this will only apply
               ! to the species SO2, SO4, MSA, NH3, NH4, NIT.
               !DVZ = MAX( DVZ, DBLE( DD_DvzMinVal_Land ) )
               dvz = max( dvz,  species_dd_dvzminval_land(species_idx) )

            ENDIF

            !-----------------------------------------------------------
            ! Special treatment for ACETONE
            !-----------------------------------------------------------

            ! For ACET, we need to only do drydep over the land
            ! and not over the oceans.
            !IF ( N == id_ACET ) THEN
            IF ( (spc == 'ACET') .or. (spc == 'acet') ) THEN
               IF ( island ) THEN
                  dvz = 0.1e+0_fp
               ELSE
                  dvz = 0e+0_fp
               ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! Special treatment for ALD2,MENO3,ETNO3,MOH
            !-----------------------------------------------------------

            ! we need to only do drydep over the land
            ! and not over the oceans.
            !IF ( N == id_ALD2 ) THEN
            IF ( (spc == 'ALD2') .or. (spc == 'MENO3') .or. (spc == 'ETNO3') .or. (spc == 'MOH') .or. &
               (spc == 'ald2') .or. (spc == 'meno3') .or. (spc == 'etno3') .or. (spc == 'moh') ) THEN
               IF ( .not. island ) THEN
                  dvz = 0e+0_fp
               ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! Compute drydep velocity and frequency
            !-----------------------------------------------------------

            ! Dry deposition velocities [m/s]
            vd = dvz / 100.e+0_fp * params%scale_factor

            ! Dry deposition frequency [1/s]
            ddfreq = vd / bxheight(1)

            ! Ensure non-negative emissions
            species_tendencies(k, species_idx) = max(0.0_fp, ddfreq)

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
                        max(0.0_fp, species_conc(k,species_idx) * (1.0_fp - exp(-1.0_fp * species_tendencies(k, species_idx) * tstep)))
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
                     drydep_velocity_per_species(diag_idx) = vd
                     exit
                  end if
               end do
            end if
         end do

      end do

   end subroutine compute_wesely

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines

   subroutine wesely_rc_gas(  RADIAT, TEMP, SUNCOS, F0, HSTAR, XMW, USTAR, CFRAC, PRESSU,     &
      XLAI,   II,   SPC, SALINITY, TSKIN, IODIDE, XLON, YLAT,         &
      CO2_EFFECT, CO2_LEVEL, CO2_REF, RSURFC,   RC)
      IMPLICIT NONE
      ! Parameters
      !-----------
      real(fp), intent(in)  :: RADIAT
      real(fp), intent(in)  :: TEMP
      real(fp), intent(in)  :: SUNCOS
      real(fp), intent(inout)  :: F0
      real(fp), intent(in)  :: HSTAR
      real(fp), intent(in)  :: XMW
      real(fp), intent(in)  :: USTAR
      real(fp), intent(in)  :: CFRAC
      real(fp), intent(in)  :: PRESSU
      real(fp), intent(in)  :: XLAI
      integer,  intent(in)  :: II
      !integer,  intent(in)  :: N_SPC      !< Species ID (TODO: may be changed to species name)
      character(len=20), intent(in) :: SPC
      !some inputs are for O3 over water and Hg over Amazon forest (not sure if we should include them for now)
      real(fp), intent(in)  :: SALINITY
      real(fp), intent(in)  :: TSKIN
      real(fp), intent(in)  :: IODIDE
      real(fp), intent(in)  :: XLON
      real(fp), intent(in)  :: YLAT
      logical,  intent(in)  :: CO2_EFFECT
      real(fp), intent(in)  :: CO2_LEVEL
      real(fp), intent(in)  :: CO2_REF
      !output
      real(fp), intent(out) :: RSURFC
      integer, intent(out)  :: RC

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
      rc      =  cc_success
      errmsg  = ''
      thisloc = ' -> at Wesely_Rc_Gas (in process/drydep/scheme/DryDepScheme_WESELY_Mod.F90)'

      ! Zero variables that aren't zeroed below
      rsurfc     = 0.0_fp
      ri         = 0.0_fp
      rlu        = 0.0_fp
      rac        = 0.0_fp
      rgss       = 0.0_fp
      rgso       = 0.0_fp
      rcls       = 0.0_fp
      rclo       = 0.0_fp
      rix        = 0.0_fp
      gfact      = 0.0_fp
      gfaci      = 0.0_fp
      rdc        = 0.0_fp
      xmwh2o     = 0.0_fp
      rixx       = 0.0_fp
      rluxx      = 0.0_fp
      rgsx       = 0.0_fp
      rclx       = 0.0_fp
      dtmp1      = 0.0_fp
      dtmp2      = 0.0_fp
      dtmp3      = 0.0_fp
      dtmp4      = 0.0_fp
      !N_SPC      = 0
      alpha0     = 0.0_fp
      depvw      = 0.0_fp

      !** TEMPK and TEMPC are surface air temperatures in K and in C
      tempk = temp
      tempc = temp - 273.15e+0_fp

      !* Adjust external surface resistances for temperature;
      !* from Wesely [1989], expression given in text on p. 1296.
      !*
      !* BUG FIX!  Wesely [1989] gives RT = 1000.0*EXP(-TEMPC-4.0)
      !*        RT = 1000.0*EXP(-(TEMPC-4.0))
      rt = 1000.0e+0_fp*exp(-tempc-4.0e+0_fp)

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
      ri = iri(ii)
      IF (ri   .GE. 9999.e+0_fp) THEN
         ri   = 1.e+12_fp
         rix   = 1.e+12_fp
      ELSE
         gfact = 100.0e+0_fp
         IF (tempc .GT. 0.e+0_fp .AND. tempc .LT. 40.e+0_fp) THEN
            gfact = 400.e+0_fp/tempc/(40.0e+0_fp-tempc)
         ENDIF

         gfaci = 100.e+0_fp
         IF ( radiat > 0.e+0_fp .and. xlai > 0.e+0_fp ) THEN
            gfaci = 1.e+0_fp / biofit( drycoeff,  xlai, suncos, cfrac, SIZE(drycoeff) )
         ENDIF
         rix = ri*gfact*gfaci
         ! Apply scaling factor to RIX when CO2 effect is turned
         ! on based on Franks et al. (2013)
         If (co2_effect) THEN
            rs_scale = co2_level / co2_ref *                   &
               (co2_level + 80.0_fp) *          &
               (co2_ref   - 40.0_fp) /          &
               (co2_level - 40.0_fp) /          &
               (co2_ref   + 80.0_fp)
            rix = rix * rs_scale
         ENDIF

      ENDIF

      !*Cuticular resistances IRLU array defined above are per unit area of leaf;
      !*divide them by the leaf area index to get a cuticular resistance for the bulk canopy.
      !*If IRLU is '9999' it means there are no cuticular surfaces on which to deposit so
      !*we impose a very large value for RLU.
      !TODO: not sure if XLAI is land type dependent or not.
      IF ( irlu(ii) >= 9999 .or. xlai <= 0.e+0_fp ) THEN
         rlu = 1.e+6_fp
      ELSE
         !RLU = DBLE( IRLU(II) ) / XLAI
         rlu =  irlu(ii)  / xlai
         ! Additional resistance at low temperatures.Limit increase to a factor of 2.
         ! Ref Jaegle et al. 2018
         rlu = min( rlu + rt, 2.e+0_fp * rlu )
      ENDIF

      !*The following are the remaining resistances for the Wesely model for a surface canopy
      !*(Wesely 1989, Fig.1).
      !RAC  = MAX(DBLE(IRAC(II)), 1.e+0_fp)
      rac  = max(irac(ii), 1.e+0_fp)
      IF (rac  .GE. 9999.e+0_fp) rac  = 1.e+12_fp
      !RGSS = MAX(DBLE(IRGSS(II)), 1.e+0_fp)
      rgss = max(irgss(ii), 1.e+0_fp)
      ! Additional resistance at low temperatures.Limit increase to a factor of 2.
      ! Ref Jaegle et al. 2018
      rgss = min( rgss + rt, 2.e+0_fp * rgss )
      IF (rgss .GE. 9999.e+0_fp) rgss = 1.e12_fp
      !RGSO = MAX(DBLE(IRGSO(II)) ,1.e+0_fp)
      rgso = max(irgso(ii) ,1.e+0_fp)
      rgso = min( rgso + rt, 2.e+0_fp * rgso )
      IF (rgso .GE. 9999.e+0_fp) rgso = 1.e+12_fp
      !RCLS = DBLE(IRCLS(II))
      rcls = ircls(ii)
      rcls = min( rcls + rt, 2.e+0_fp * rcls )
      IF (rcls .GE. 9999.e+0_fp) rcls = 1.e+12_fp
      !RCLO = DBLE(IRCLO(II))
      rclo = irclo(ii)
      rclo = min( rclo + rt, 2.e+0_fp * rclo )
      IF (rclo .GE. 9999.e+0_fp) rclo = 1.e+12_fp

      !* Compute aerodynamic resistance to lower elements in lower part
      !* of the canopy or structure, assuming level terrain -
      !* equation (5) of Wesely [1989].
      !* species-dependent corrections to resistances
      !* are from equations (6)-(9) of Wesely [1989].

      rdc = 100.e+0_fp*(1.0e+0_fp+1000.0e+0_fp/(radiat+10.e+0_fp))

      IF ( spc .EQ. 'O3' .or. spc .EQ. 'o3') THEN
         !O3 over water
         IF ((ii .EQ. 11)) THEN
            IF (salinity .GT. 20.0_fp) THEN
               ! Now apply the Luhar et al. [2018] equations for the
               ! special treatment of O3 dry deposition to the ocean
               CALL oceano3(tskin,ustar,iodide,depvw)
               ! Now convert to the new rc value(s)
               alpha0 = 10.0_fp**(-0.25_fp-0.013_fp * (tskin-273.16_fp))
               rsurfc = 1.0_fp/(alpha0*depvw)
            ELSE
               ! It's not saline enough for 'ocean' so we instead don't change it from
               ! 'default' rc to water
               rsurfc = 2000.0_fp
            ENDIF
         ENDIF

         !O3 over snow/ice, the surface resistance is set to an observation derived value
         IF ((ii .EQ. 1)) THEN
            rsurfc = 10000.0_fp
         ENDIF
      ELSE
         !set a different F0 for Hg0
         IF (spc .EQ. 'Hg0' .or. spc .EQ. 'hg0') THEN
            ! Assume lower reactivity
            f0 = 3.0e-05_fp
            ! But if this is the rainforest land type and we fall
            ! within the bounding box of the Amazon rainforest,
            ! then increase reactivity as inferred from observations.
            IF ( ii  ==  6          .AND.             &
               xlon >  -82.0_fp   .AND.             &
               xlon <  -33.0_fp   .AND.             &
               ylat >  -34.0_fp   .AND.             &
               ylat <   14.0_fp ) THEN
               f0 = 2.0e-01_fp
            ENDIF
         ENDIF

         xmwh2o = h2omw * 1.e-3_fp
         rixx = rix*diffg(tempk,pressu,xmwh2o)/ diffg(tempk,pressu,xmw) &
            + 1.e+0_fp/(hstar/3000.e+0_fp+100.e+0_fp*f0)
         rluxx = 1.e+12_fp
         IF (rlu .LT. 9999.e+0_fp) rluxx = rlu/(hstar/1.0e+05_fp + f0)
         rgsx = 1.e+0_fp/(hstar/1.0e+05_fp/rgss + f0/rgso)
         rclx = 1.e+0_fp/(hstar/1.0e+05_fp/rcls + f0/rclo)
         !** Get the bulk surface resistance of the canopy, RSURFC, from
         !** the network of resistances in parallel and in series (Fig.1 of Wesely [1989])
         dtmp1=1.e+0_fp/rixx
         dtmp2=1.e+0_fp/rluxx
         dtmp3=1.e+0_fp/(rac+rgsx)
         dtmp4=1.e+0_fp/(rdc+rclx)
         rsurfc = 1.e+0_fp/(dtmp1 + dtmp2 + dtmp3 + dtmp4)
      ENDIF

      !TODO: this is put in the main scheme function
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
   end subroutine wesely_rc_gas


   SUBROUTINE oceano3( TEMPK, USTAR, IODIDE_IN, DEPV )

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

      ustarwater = 0.0345_fp * ustar !waterside friction velocity

      iodide = iodide_in*1.0e-9_fp ! Convert from nM to M

      a0 = iodide*exp((-8772.2_fp/tempk)+51.5_fp) !chemical reactivity

      d = 1.1e-6_fp*exp(-1896.0_fp/tempk) ! diffusivity

      delm = sqrt(d/a0) ! reaction-diffusion length

      b = 2.0_fp/(0.4_fp*ustarwater)

      lam = delm*sqrt(a0/d) ! this cancels to 1 but here for completeness of equations

      ep = sqrt(2.0_fp*a0*b*(delm+(b*d/2.0_fp)))

      psi = ep/sqrt(a0*b**2*d)

      CALL k0k1_aprox(ep,k0,k1)

      depv = sqrt(a0*d)*((psi*k1*cosh(lam)+k0*sinh(lam))/(psi*k1* sinh(lam)+k0*cosh(lam)))

   END SUBROUTINE oceano3

   SUBROUTINE k0k1_aprox( input_arg, K0, K1 )

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
         coeff = (/1.0_fp,3.5156229_fp,3.0899424_fp,1.2067492_fp,0.2659732_fp, &
            0.360768e-1_fp,0.45813e-2_fp/)
         i0 = poly_fit((input_arg/3.75_fp)**2,coeff)
         !now we can use this estimate of i0 to calculate k0
         coeff = (/-0.57721566_fp,0.42278420_fp,0.23069756_fp,0.3488590e-1_fp, &
            0.262698e-2_fp,0.10750e-3_fp,0.74e-5_fp/)
         k0 = (-log(0.5_fp*input_arg)*i0)+ &
            poly_fit(0.25_fp*input_arg**2,coeff)

         !begin the calculation of k0 by estimating i1
         coeff = (/0.5_fp,0.87890594_fp,0.51498869_fp,0.15084934_fp,0.2658733e-1_fp, &
            0.301532e-2_fp,0.32411e-3_fp/)
         i1 = input_arg*poly_fit((input_arg/3.75_fp)**2,coeff)
         ! now we can use this to estimate to get a value for k1
         coeff = (/1.0_fp,0.15443144_fp,-0.67278579_fp,-0.18156897_fp, &
            -0.1919402e-1_fp,-0.110404e-2_fp,-0.4686e-4_fp/)
         k1 = (log(0.5_fp*input_arg)*i1)+(1.0_fp/input_arg)* &
            poly_fit(0.25_fp*input_arg**2,coeff)
      ELSE !use a different approximation that doesn't need I0/I1
         coeff = (/1.25331414_fp,-0.7832358e-1_fp,0.2189568e-1_fp,-0.1062446e-1_fp, &
            0.587872e-2_fp,-0.251540e-2_fp,0.53208e-3_fp/)
         k0 = (exp(-input_arg)/sqrt(input_arg))* &
            poly_fit((2.0_fp/input_arg),coeff)
         coeff = (/1.25331414_fp,0.23498619_fp,-0.3655620e-1_fp,0.1504268e-1_fp, &
            -0.780353e-2_fp,0.325614e-2_fp,-0.68245e-3_fp/)
         k1 = (exp(-input_arg)/sqrt(input_arg))* &
            poly_fit((2.0_fp/input_arg),coeff)
      ENDIF

   END SUBROUTINE k0k1_aprox

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


   FUNCTION diffg( TK, PRESS, XM ) RESULT( DIFF_G )

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
      airden = ( press * avo ) / ( rstarg * tk )

      ! DIAM is the collision diameter for gas X with air.
      diam   = radx + radair

      ! Calculate the mean free path for gas X in air:
      ! eq. 8.5 of Seinfeld [1986];
      z      = xm  / xmair
      frpath = 1e+0_fp /( pi * sqrt( 1e+0_fp + z ) * airden * ( diam**2 ) )

      ! Calculate average speed of gas X; eq. 15.47 of Levine [1988]
      speed  = sqrt( 8e+0_fp * rstarg * tk / ( pi * xm ) )

      ! Calculate diffusion coefficient of gas X in air;
      ! eq. 8.9 of Seinfeld [1986]
      diff_g = ( 3e+0_fp * pi / 32e+0_fp ) * ( 1e+0_fp + z ) * frpath * speed

   END FUNCTION diffg

   FUNCTION biofit( COEFF1, XLAI1, SUNCOS1, CFRAC1, NPOLY ) RESULT( BIO_FIT )

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
      term(1) = 1.0e0_fp
      term(2) = xlai1
      term(3) = suncos1
      term(4) = cfrac1
      !we replace SUNPARAM_R4( TERM(2:4) ) as below
      !outdate lai,suncos,cloud fraction
      DO i = 1, nn
         i2 = i + 1 !variable index in TERM is from 2
         term(i2) = min( term(i2), x0(i) )
         ! XLOW = minimum for each variable
         IF ( i .NE. 3 ) THEN
            xlow = x0(i) / nd(i) !codespell:ignore
         ELSE
            xlow = 0.0e0_fp
         ENDIF
         term(i2) = max( term(i2), xlow )
         term(i2) = term(i2) / x0(i)
      ENDDO

      !get realterm
      k = 0
      DO k3 = 1, kk
         DO k2 = k3, kk
            DO k1 = k2, kk
               k = k + 1
               realterm(k)=term(k1)*term(k2)*term(k3)
            ENDDO
         ENDDO
      ENDDO

      bio_fit = 0e0_fp
      DO k = 1, npoly
         bio_fit = bio_fit + coeff1(k)*realterm(k)
      END DO
      IF ( bio_fit .LT. 0.1e0_fp ) bio_fit = 0.1e0_fp

   END FUNCTION biofit

   subroutine wesely_ra_rb(TEMP, PRESSU, XMW, USTAR, OBK, ZO, THIK, IS_GAS, Ra, Rb, RC)
      IMPLICIT NONE
      ! Parameters
      !-----------
      real(fp), intent(in)  :: TEMP
      real(fp), intent(in)  :: PRESSU
      real(fp), intent(in)  :: XMW
      real(fp), intent(in)  :: USTAR
      real(fp), intent(in)  :: OBK
      real(fp), intent(in)  :: ZO
      real(fp), intent(in)  :: THIK
      logical, intent(in)   :: IS_GAS
      !output
      real(fp), intent(out) :: Ra
      real(fp), intent(out) :: Rb
      integer, intent(out)  :: RC

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
      rc      =  cc_success
      errmsg  = ''
      thisloc = ' -> at Wesely_Ra_Rb (in process/drydep/scheme/DryDepScheme_WESELY_Mod.F90)'

      ! Zero variables that aren't zeroed below
      cz         = 0.0_fp
      ckustr     = 0.0_fp
      reyno      = 0.0_fp
      corr1      = 0.0_fp
      corr2      = 0.0_fp
      z0obk      = 0.0_fp
      dummy1     = 0.0_fp
      dummy2     = 0.0_fp
      dummy3     = 0.0_fp
      dummy4     = 0.0_fp
      dair       = 0.0_fp
      ra         = 0.0_fp
      rb         = 0.0_fp

      !CZ is Altitude (m) at which deposition velocity is computed
      !use Midpoint height of first model level [m]
      cz = thik / 2.0e+0_fp

      !** TEMPK and TEMPC are surface air temperatures in K and in C
      tempk = temp
      tempc = temp-273.15e+0_fp

      !** Calculate the kinematic viscosity XNU (m2 s-1) of air
      !** as a function of temperature.
      !** The kinematic viscosity is used to calculate the roughness heights
      !** over water surfaces and to diagnose whether such surfaces are
      !** aerodynamically rough or smooth using a Reynolds number criterion.
      !** The expression for the temperature dependence of XNU
      !** is from the FORTRAN code in Appendix II of Wesely [1988];
      !** I wasn't able to find an original reference but it seems benign enough.
      c1  = tempk/273.15e+0_fp
      xnu = 0.151e+0_fp*(c1**1.77e+0_fp)*1.0e-04_fp

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

      ckustr = von_karman * ustar
      reyno = ustar*zo/xnu
      corr1 = cz/obk
      ! Define Z0OBK
      z0obk = zo/obk

      lrgera = .false.
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
      IF ( rae(ckustr, 0.0e+0_fp) ) THEN
         errmsg = 'CKUSTR cannot be zero.'
         CALL cc_error( errmsg, rc, thisloc )
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
      IF ( reyno >= 0.1e+0_fp ) THEN !rough surface
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

         IF (corr1.LT.0.0e+0_fp) THEN
            !*... unstable conditions; compute Ra as described
            !*... above.
            !coef_a=1.e+0_fp
            !coef_b=15.e+0_fp
            dummy1 = (1.e+0_fp - 15.e+0_fp*corr1)**0.5e+0_fp
            dummy2 = (1.e+0_fp - 15.e+0_fp*z0obk)**0.5e+0_fp
            dummy3 = abs((dummy1 - 1.e+0_fp)/(dummy1 + 1.e+0_fp))
            dummy4 = abs((dummy2 - 1.e+0_fp)/(dummy2 + 1.e+0_fp))
            ra = 1.e+0_fp * (1.e+0_fp/ckustr) * log(dummy3/dummy4)

         ELSEIF((corr1.GE.0.0e+0_fp).AND.(corr1.LE.1.0e+0_fp)) THEN
            !coef_a=1.e+0_fp
            !coef_b=5.e+0_fp
            ra = (1.e+0_fp/ckustr) * (1.e+0_fp*log(corr1/z0obk) + &
               5.e+0_fp*(corr1-z0obk))

         ELSE ! CORR1 .GT. 1.0D0
            !coef_a=5e+0_fp
            !coef_b=1.e+0_fp
            ra = (1.e+0_fp/ckustr) * (5.e+0_fp*log(corr1/z0obk) + &
               1.e+0_fp*(corr1-z0obk))
         ENDIF

         !* check that RA is positive and maximize at 1.E4 s m-1
         ra   = min(ra,1.e+4_fp)
         ! If RA is < 0, set RA = 0
         IF (ra .LT. 0.e+0_fp) ra = 0.0e+0_fp

         !END IF !PBL or non-local PBL options

         !get Rb for a gas species; arosol Rb is set to zero
         !** DAIR is the thermal diffusivity of air; value 0.2*1.E-4 m2 s-1 cited on p. 16,476 of
         !** Jacob et al. [1992]
         dair = 0.2e0_fp*1.e-4_fp
         IF (is_gas) THEN
            rb = (2.e+0_fp/ckustr)* (dair/diffg(tempk,pressu,xmw)) &
               **0.667e+0_fp
         END IF

      ELSE  !smooth surface
         !** suppress drydep over smooth surfaces by setting Ra to
         !** a large value (1e4).  This prevents negative dry deposition
         !** velocities when u* is very small. Rb is not important in that case since
         !** the total resistentce is Ra + Rb. So we set Rb to zero.
         ra     = 1.e+4_fp

      END IF

   end subroutine wesely_ra_rb


end module drydepscheme_wesely_mod
```


