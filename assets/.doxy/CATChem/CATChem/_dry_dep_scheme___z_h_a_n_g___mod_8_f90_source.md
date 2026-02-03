

# File DryDepScheme\_ZHANG\_Mod.F90

[**File List**](files.md) **>** [**drydep**](dir_57fb5aa14ddb2cd518a6d90b65ffd000.md) **>** [**schemes**](dir_5a3c86e36f17958630366ebc2b7ca21b.md) **>** [**DryDepScheme\_ZHANG\_Mod.F90**](_dry_dep_scheme___z_h_a_n_g___mod_8_f90.md)

[Go to the documentation of this file](_dry_dep_scheme___z_h_a_n_g___mod_8_f90.md)


```Fortran

module drydepscheme_zhang_mod

   use precision_mod, only: fp, rae, f8
   use error_mod, only: cc_success, cc_error
   use drydepcommon_mod, only: drydepschemezhangconfig
   use constants, only: pi, avo, von_karman, rstarg, g0, boltz  !load the constants needed for this scheme

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
   DATA idep_iolson  / 11,10, 5, 3, 3, 2, 2, 5, 8, 7,  5,  8,  1,  9,  11, 11, 5,   5,  5,  2, 6,  3,   3,  2,  2,  2,  2, &
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

   DATA   a / 2.0e+0_fp,   5.0e+0_fp,   2.0e+0_fp,   5.0e+0_fp,  5.0e+0_fp, &
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
      integer  :: II
      integer  :: ILDT
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

      rc = cc_success
      errmsg = ''
      thisloc = ' -> at compute_zhang (in src/process/drydep/schemes/DryDepScheme_ZHANG_Mod.F90)'

      !calculate the volume distribution of sea salt aerosols (only need to do this once)
      IF ( firsttime ) THEN
         ! Derive seasalt bin boundaries from species properties
         call get_seasalt_bin_boundaries(num_species, species_is_seasalt, &
            species_lower_radius, species_upper_radius, &
            seasalt_lower_bin, seasalt_upper_bin)

         CALL init_weightss(minval(seasalt_lower_bin), maxval(seasalt_upper_bin), rc)
         IF ( rc /= cc_success ) THEN
            errmsg = 'Could not Allocate arrays in INIT_WEIGHTSS'
            CALL cc_error( errmsg, rc, thisloc )
            RETURN
         ENDIF
         firsttime = .false.
      ENDIF

      ! Calculate 10m wind speed
      w10 = sqrt(u10m**2 + v10m**2)

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
            vd         = 0.0_fp
            ddfreq     = 0.0_fp
            dvz        = 0.0_fp
            rsurfc     = 0.0_fp
            ra         = 0.0_fp
            rb         = 0.0_fp
            c1x        = 0.0_fp
            vk         = 0.0_fp
            vtsoutput  = 0.0_fp

            !property for current species
            hstar = species_dd_hstar(species_idx)
            xmw   = species_mw_g(species_idx)*1e-3_fp   !convert from g/mol to kg/mole
            spc   = trim(species_short_name(species_idx))

            ! Better test for depositing species: We need both HSTAR and XMW
            ! to be nonzero, OR the value of AIROSOL to be true.  This should
            ! avoid any further floating point invalid issues caused by putting
            ! a zero value in a denominator.
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
                  lucindex = lucindex_gc(ii)
               ELSE IF ( lucname == 'NOAH' ) THEN
                  ! it is possible that water is given as 0 not 17 in GFS CCPP
                  IF (ildt == 0) ildt = 17
                  !II   = IDEP_NOAH(ILDT)
                  !Note: we use ILDT, instead of II,  to get LUCINDEX here
                  lucindex = lucindex_noah(ildt)
               ELSE IF ( lucname == 'IGBP' ) THEN
                  ! it is possible that water is given as 0 not 17
                  IF (ildt == 0) ildt = 17
                  !II   = IDEP_IGBP(ILDT)
                  lucindex = lucindex_igbp(ildt)
               ENDIF

               !get bulk surface resistances (Rs)
               !Note to change pressure unit from Pa to kPa
               rsurfc = aero_sfcrsii( spc, species_is_dust(species_idx), species_is_seasalt(species_idx), lucindex, &
                  species_radius(species_idx)*1e-6_fp, species_density(species_idx), ps*1e-3_fp, & !um to m; Pa to kPa
                  ts, ustar, rh(1), w10, seasalt_lower_bin, seasalt_upper_bin,vtsoutput, rc)

               if (rc /= cc_success ) then
                  errmsg = 'Error in getting bulk surface resistances (RSURFC)'
                  CALL cc_error( errmsg, rc, thisloc )
                  RETURN
               endif

               !*Set max and min values for bulk surface resistances
               rsurfc = max(1.e+0_fp, min(rsurfc,9999.e+0_fp))
               ! Set Rc for strong acids (HNO3,HCl,HBr) to 1 s/m
               ! Ref. Jaegle et al. 2018, cf. Erisman,van Pul,Ayers 1994
               IF ( hstar .gt. 1.e+10_fp ) rsurfc= 1.e+0_fp

               !get Ra and Rb
               call wesely_ra_rb(ts, ps, xmw, ustar, obk, z0, bxheight(1), .false., ra, rb,  rc)

               !get VD (TODO: IUSE is decimal not percent or permille as in GEOS-Chem)
               c1x = rsurfc + ra + rb
               vk = vd
               !VD = VK + DBLE( IUSE(LDT) ) / C1X + DBLE( IUSE(LDT) ) * VTSoutput
               vd = vk +  frlanduse(ldt)  / c1x +  frlanduse(ldt) * vtsoutput
            END DO

            !apply spectial treatment or scaling factor to Vd
            dvz = vd *100.e+0_fp !m/s -- > cm/s

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

   end subroutine compute_zhang

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines
   ! Examples: environmental response functions, species-specific calculations, etc.

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

      rc = cc_success
      errmsg = ''
      thisloc = ' -> at get_seasalt_bin_boundaries (in src/process/drydep/schemes/DryDepScheme_ZHANG_Mod.F90)'

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
            if (all(abs(temp_lower1(1:n_bin) - lower_radius(n)) > 0.0_fp)) then
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
      mask(1:n_bin) = .true.
      do n = 1, n_bin
         lower_bin(n) =  minval(temp_lower,mask)
         mask(minloc(temp_lower,mask)) = .false.
      enddo

      !sort bins by radius from low to high for upper_radius
      mask(1:n_bin) = .true.
      do n = 1, n_bin
         upper_bin(n) =  minval(temp_upper,mask)
         mask(minloc(temp_upper,mask)) = .false.
      enddo

      !check if the bins are continuous
      do n = 1, n_bin-1
         if ( .not. rae(upper_bin(n), lower_bin(n+1)) ) then
            errmsg = 'Sea Salt Bins are not continuous'
            call cc_error(errmsg, rc, thisloc)
            RETURN
         endif
      enddo

      ! Clean up
      deallocate(temp_lower, temp_upper, temp_lower1, temp_upper1,mask)

   end subroutine get_seasalt_bin_boundaries


   FUNCTION aero_sfcrsii(  SPC, IS_DUST, IS_SEASALT, LUCINDEX, A_RADI, A_DEN, &
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
      rc = cc_success
      errmsg = ''
      thisloc = ' -> at AERO_SFCRSII (in process/drydep/scheme/DryDepScheme_ZHANG_Mod.F90)'

      !=================================================================
      ! ADUST_SFCRII begins here!
      !=================================================================

      ! Annual average of A
      aavg(:) = (a(:,1)+a(:,2)+a(:,3)+a(:,4)+a(:,5))/5.0_fp

      luc     = lucindex
      aa      = aavg(luc) * 1.e-3_fp
      rs = 0e+0_fp !initialize returned value first

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
      diam  = a_radi * 2.e+0_fp

      ! Particle density [kg/m3]
      den   = a_den

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
      rhbl    = max( tiny(rhb), rhb )

      ! Over oceans the RH in the viscous sublayer is set to 98%,
      ! following Lewis and Schwartz (2004)
      !I added condition when RHBL=1 to avoid DIAM = infinity issue in the New_DIAM_DEN subroutine later for SO4 (Wei Li)
      IF (luc == 14 .or. rae(rhbl, 1.0_fp)) THEN
         rhbl = 0.98_fp
      ENDIF

      IF (.NOT. is_dust) THEN
         !update DIAM and DEN after hygroscopic growth for non-dust species
         call new_diam_den( spc, is_seasalt, rhbl, rdry, rwet, diam, den, rc)
         if (rc /= cc_success ) then
            errmsg = 'New_DIAM_DEN failed.'
            CALL cc_error( errmsg, rc, thisloc )
            RETURN
         endif
      ENDIF

      ! Dp [m] --> [um] = particle diameter if necessary
      IF (diam > 0.001) THEN !here use 0.001 to determine if the unit is in m or um
         dp    = diam
         diam  = diam * 1.e-6_fp
      ELSE
         dp    = diam * 1.e+6_fp
      ENDIF

      ! Constant for settling velocity calculation
      const = den * diam**2 * g0 / 18.e+0_fp

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
      pdp  = press * dp
      slip = 1e+0_fp + ( 15.60e+0_fp + 7.0e+0_fp * &
         exp( -0.059e+0_fp * pdp) ) / pdp

      ! Viscosity [Pa s] of air as a function of temp (K)
      visc = 1.458e-6_fp * (temp)**(1.5e+0_fp) / (temp + 110.4e+0_fp)

      ! Kinematic viscosity (Dynamic viscosity/Density)
      airvs= visc / 1.2928e+0_fp

      ! Settling velocity [m/s]
      vts  = const * slip / visc
      !sea salt VTS update
      IF (is_seasalt) THEN
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
         nr = int((( maxval(upperbin) - minval(lowerbin) ) / dr ) + 0.5e+0_fp )

         salt_mass_total = 0e+0_fp
         vts_weight      = 0e+0_fp

         ! Dry particle radius [m] --> [um]
         rum  = rdry * 1.e+6_fp

         ! Check what the min/max range of the SS size bins are
         !IF ( RUM .le. SALA_REDGE_um(2) ) THEN
         !   D0 = SALA_REDGE_um(1)*2e+0_fp
         !   D1 = SALA_REDGE_um(2)*2e+0_fp
         !ELSE
         !   D0 = SALC_REDGE_um(1)*2e+0_fp
         !   D1 = SALC_REDGE_um(2)*2e+0_fp
         !ENDIF
         d0 = 0e+0_fp; d1 = 0e+0_fp; n1=1
         DO n =1, size(upperbin)
            IF ( (rum .ge. lowerbin(n)) .and. (rum .le. upperbin(n)) ) THEN
               d0 = lowerbin(n)*2e+0_fp
               d1 = upperbin(n)*2e+0_fp
               n1=0
               EXIT
            ENDIF
         ENDDO
         IF (n1 > 0) THEN ! D0 and D1 may not be set properly
            errmsg = 'Sea salt radius is not in any bins. Check the species namelist.'
            CALL cc_error( errmsg, rc, thisloc )
            RETURN
         ENDIF

         DO id = 1, nr
            ! Calculate mass of wet aerosol (Dw = wet diameter, D = dry diameter):
            ! Overall = dM/dDw = dV/dlnD * Rwet/Rdry * DEN /Rw
            !TODO: DMID is not defined in this module. Need to define it.
            IF (dmid(id) .ge. d0 .and. dmid(id) .le. d1 ) THEN
               dmidw = dmid(id) * rwet/rdry   ! wet radius [um]
               salt_mass   = salt_v(id) * rwet/rdry * den / &
                  (dmidw*0.5e+0_fp)
               vts_weight  = vts_weight + &
               !SALT_MASS * VTS * (DMIDW/(RWET*1d6*2e+0_fp) )** &
                  salt_mass * vts * (dmidw/(rwet*1e+6_fp*2e+0_fp) )** &
                  2e+0_fp * (2e+0_fp * dr *  rwet/rdry)
               salt_mass_total = salt_mass_total+salt_mass * &
                  (2e+0_fp * dr *  rwet/rdry)
            ENDIF
         ENDDO

         ! Final mass weighted setting velocity:
         vts = vts_weight/salt_mass_total
      END IF

      vtsout = vts !need to save out for final Vd calculation

      ! Brownian diffusion constant for particle (m2/s)
      diff = boltz * temp * slip / (3.e+0_fp * pi * visc * diam)

      ! Schmidt number
      sc   = airvs / diff
      !EB   = 1.e+0_fp/SC**(gamma(LUC))

      !--------------------------------------------------------------
      ! NOTE: This loses precision, use TWO_THIRDS parameter instead
      !EB   = CB/SC**(0.6667e+0_fp) ! Emerson 2020 update JRP
      !--------------------------------------------------------------
      eb   = cb/sc**two_thirds ! Emerson 2020 update JRP

      ! Stokes number
      IF ( aa < 0e+0_fp ) then
         st   = vts * ustar * ustar / ( airvs * g0 ) ! for smooth surface
         ein  = 0e+0_fp
      ELSE
         st   = vts   * ustar / ( g0 * aa )          ! for vegetated surfaces
         !EIN  = 0.5e+0_fp * ( DIAM / AA )**2
         ein  = cin * ( diam / aa )**(upsilon) ! Emerson 2020 update JRP
      ENDIF

      IF (luc == 14 .and. is_seasalt) THEN
         eim  = 10.e+0_fp**( -3.e+0_fp/ st )         ! for water surface
         ! JRP: Emerson doesn't describe what to do here, so I'm leaving as is
      ELSE
         !EIM  = ( ST / ( ALPHA(LUC) + ST ) )**(BETA)
         eim  = cim * ( st / ( alpha(luc) + st ) )**(beta) ! Emerson 2020 update JRP
         eim  = min( eim, 0.6e+0_fp )
      ENDIF

      IF (luc == 11 .OR. luc == 13 .OR. luc == 14) THEN
         r1 = 1.e+0_fp
      ELSE
         r1 = exp( -1e+0_fp * sqrt( st ) )
         r1 = max( tiny(r1), r1 ) !avoid R1 = 0 when ST is large under very low TEMP and AA < 0 (Wei Li)
      ENDIF

      !add error check here to make sure RS below is not a infinite value
      IF (rae(r1, 0.0_fp) .or. rae(ustar, 0.0_fp)) THEN
         !write(*,*) 'DEBUG INFO: SPC=', trim(SPC), LUC, USTAR, R1, ST, AA, VTS, CONST, DEN, DIAM, RHBL, RHB, AIRVS
         errmsg = 'USTAR or R1 is zero. Check met field or diameter (in m) of aerosol is too big.'
         CALL cc_error( errmsg, rc, thisloc )
         RETURN
      ENDIF

      ! surface resistance for particle
      IF (luc == 14 .and. is_seasalt) THEN
         ! Use the formulation of Slinn and Slinn (1980) for the impaction over
         ! water surfaces for sea salt
         rs   = 1.e+0_fp / (ustar**2.e+0_fp/ (w10*von_karman) * &
            (eb + eim ) + vts)
      ELSE
         rs   = 1.e0_fp / (e0 * ustar * (eb + eim + ein) * r1 )
      ENDIF

   END FUNCTION aero_sfcrsii

   SUBROUTINE new_diam_den( SPC, IS_SEASALT, RHBL, RDRY, RWET, DIAM, DEN, RC)

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
      rc = cc_success
      errmsg = ''
      thisloc = ' -> at New_DIAM_DEN (in process/drydep/scheme/DryDepScheme_ZHANG_Mod.F90)'

      IF ( .NOT. is_seasalt ) THEN

         ! Particle diameter [m]
         diam  = 0.17378e-6_fp
         rdry = diam / 2.0e+0_fp !Not needed for further calculations for dust species

         ! SIA (TODO: keep this for now and need to be consistent with real species names in the future)
         !IF ( K == idd_NIT .or. K == idd_NH4 .or. K == idd_SO4 ) THEN
         IF ( spc == 'NIT' .or. spc == 'NH4' .or. spc == 'SO4' .or. spc == 'ASO4J' .or. &
            spc == 'nit' .or. spc == 'nh4' .or. spc == 'so4' .or. spc == 'aso4j' ) THEN
            ! Efflorescence transitions
            IF (rhbl .LT. 0.35) THEN
               ! DIAM is not changed
            ELSE IF ((rhbl .GE. 0.35) .AND. (rhbl .LE. 0.40)) THEN
               ! Linear hygroscopic growth
               diam = diam + (diam * ((1.0_fp + 0.61_fp * 0.40_fp /             &
                  (1.0_fp - 0.40_fp)) ** (1.0_fp / 3.0_fp)) - diam) /        &
                  (0.40_fp - 0.35_fp) * (rhbl - 0.35_fp)
            ELSE
               ! Kohler hygroscopic growth
               diam = diam * ((1.0_fp + 0.61_fp * rhbl / (1.0_fp - rhbl))       &
                  ** (1.0_fp / 3.0_fp))
            ENDIF

            !BC
            !ELSE IF ( K == idd_BCPI .OR. K == idd_BCPO )  THEN
         ELSE IF ( spc == 'BCPI' .OR. spc == 'BCPO' .OR. &
            spc == 'bcpi' .OR. spc == 'bcpo' )  THEN
            ! DIAM is not changed

            !OA
         ELSE
            diam = diam * ((1.0_fp + 0.1_fp * rhbl / (1.0_fp - rhbl))             &
               ** (1.0_fp / 3.0_fp))
         ENDIF

         !get RWET
         rwet = diam / 2.0e+0_fp
         ! Particle density [kg/m3]; same for all aerosols except sea salt and  dust
         den   = 1500

      ELSE !sea salt aerosol case

         !drydepRadius = A_RADI(K)
         rdry = diam / 2.0e+0_fp

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
         rum  = rdry * 1.e+6_fp

         ! Exponential factors used for hygroscopic growth (not used now)
         fac1 = c1 * ( rum**c2 )
         fac2 = c3 * ( rum**c4 )

         ! Corrected bug in Gerber formulation: use of LOG10  (jaegle 5/11/11)
         !RWET    = 0.01e+0_fp*(FAC1/(FAC2-DLOG(RHBL))+RCM**3.e+0_fp)**0.33e+0_fp
         !RWET = 1.d-6*(FAC1/(FAC2-LOG10(RHBL))+RUM**3.e+0_fp)**0.33333e+0_fp

         ! Use equation 5 in Lewis and Schwartz (2006) for sea salt growth [m]
         ! (jaegle 5/11/11)
         rwet = rdry * (4.e+0_fp / 3.7e+0_fp) * &
            ( (2.e+0_fp - rhbl)/(1.e+0_fp - rhbl) )**(1.e+0_fp/3.e+0_fp)

         ! Ratio dry over wet radii at the cubic power
         !RATIO_R = ( A_RADI(K) / RWET )**3.e+0_fp

         ! Diameter of the wet aerosol [m]
         diam  = rwet * 2.e+0_fp

         ! Density of the wet aerosol [kg/m3] (bec, 12/8/04)
         !DEN   = RATIO_R * A_DEN(K) + ( 1.e+0_fp - RATIO_R ) * 1000.e+0_fp

         ! Above density calculation is chemically unsound because it ignores chemical solvation.
         ! Iteratively solve Tang et al., 1997 equation 5 to calculate density of wet aerosol (kg/m3)
         ! Redefine RATIO_R
         ratio_r = real(rdry / rwet, f8)

         ! Assume an initial density of 1000 kg/m3
         den0 = real(den, f8) !assign initial DEN to DEN0
         den_f8  = 1000.e+0_f8
         den1 = 0.e+0_f8 !initialize
         i = 0 !initialize loop index
         !Note that if RH is too low, the loop will not converge and will run forever
         DO WHILE ( abs( den1-den_f8 ) .gt. epsi )
            ! First calculate weight percent of aerosol (kg_RH=0.8/kg_wet)
            wtp    = 100.e+0_f8 * den0/den_f8 * ratio_r**3
            ! Then calculate density of wet aerosol using equation 5
            ! in Tang et al., 1997 [kg/m3]
            den1   = ( 0.9971e+0_f8 + (a1 * wtp) + (a2 * wtp**2) + &
               (a3 * wtp**3) + (a4 * wtp**4) ) * 1000.e+0_f8

            ! Now calculate new weight percent using above density calculation
            wtp    = 100.e+0_f8 * den0/den1 * ratio_r**3
            ! Now recalculate new wet density [kg/m3]
            den_f8   = ( 0.9971e+0_f8 + (a1 * wtp) + (a2 * wtp**2) + &
               (a3 * wtp**3) + (a4 * wtp**4) ) * 1000.e+0_f8

            ! add some protection against infinite loop
            i = i+1
            IF ( i .GT. 500 ) THEN
               !write(*,*) 'Test NEW_DIAM_DEN output: ', trim(SPC), RHBL, RDRY, RWET, DIAM, DEN0, DEN,DEN1
               errmsg = 'Error in calculating new density for sea salt aerosol due to very low RH input!'
               CALL cc_error( errmsg, rc, thisloc )
               RETURN
            ENDIF

         ENDDO
         ! Convert back to fp
         den = real(den_f8, fp)
      ENDIF

   END SUBROUTINE new_diam_den


   SUBROUTINE init_weightss( SALT_RLOW_um, SALT_RUP_um, RC )

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
      errmsg = ''
      thisloc = ' -> at INIT_WEIGHTSS (in process/drydep/schemes/DryDepScheme_ZHANG_Mod.F90)'

      ! Number of bins between the lowest bound of of the accumulation mode
      ! sea salt and the upper bound of the coarse mode sea salt.
      nr = int((( salt_rup_um - salt_rlow_um )  / dr ) + 0.5e+0_fp )

      ALLOCATE( dmid( nr ), stat=rc )
      IF ( rc /= cc_success ) THEN
         errmsg = 'Could not allocate array DMID'
         CALL cc_error( errmsg, rc, thisloc )
         RETURN
      END IF
      dmid = 0e+0_fp

      ALLOCATE( salt_v( nr ), stat=rc )
      IF ( rc /= cc_success ) THEN
         errmsg = 'Could not allocate array SALT_V'
         CALL cc_error( errmsg, rc, thisloc )
         RETURN
      END IF
      salt_v = 0e+0_fp

      !=================================================================
      ! Define the volume size distribution of sea-salt. This only has
      ! to be done once. We assume that sea-salt is the combination of a
      ! coarse mode and accumulation model log-normal distribution functions
      !=================================================================

      ! Lower edge of 0th bin diameter [um]
      dedge=salt_rlow_um * 2e+0_fp

      ! Loop over diameters
      DO id = 1, nr

         ! Diameter of mid-point in microns
         dmid(id)  = dedge + ( dr )

         ! Calculate the dry volume size distribution as the sum of two
         ! log-normal size distributions. The parameters for the size
         ! distribution are based on Reid et al. and Quinn et al.
         ! The scaling factors 13. and 0.8 for acc and coarse mode aerosols
         ! are chosen to obtain a realistic distribution
         ! SALT_V (D) = dV/dln(D) [um3]
         salt_v(id) = pi / 6e+0_fp* (dmid(id)**3) * (         &
            13e+0_fp*exp(-0.5_fp*( log(dmid(id))-       &
            log(rg_a*2e+0_fp) )**2e+0_fp/            &
            log(sig_a)**2e+0_fp )           &
            /( sqrt(2e+0_fp * pi) * log(sig_a) )  +  &
            0.8e+0_fp*exp(-0.5_fp*( log(dmid(id))-      &
            log(rg_c*2e+0_fp) )**2e+0_fp/            &
            log(sig_c)**2e+0_fp)            &
            /( sqrt(2e+0_fp * pi) * log(sig_c) )  )

         ! update the next edge
         dedge = dedge + dr*2e+0_fp
      ENDDO

   END SUBROUTINE init_weightss


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
      thisloc = ' -> at Wesely_Ra_Rb (in process/drydep/scheme/DryDepScheme_ZHANG_Mod.F90)'

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



end module drydepscheme_zhang_mod
```


