!>
!! \file
!! \brief CCPr Scheme for dry deposition
!!
!!
!! Reference:
!! (1) Wesely, M. L. (1989). Parameterization of surface resistances to gaseous dry
!!     deposition in regional-scale numerical models. Atmospheric Environment.
!! (2) Zhang, L., Gong, S., Padro, J., & Barrie, L. (2001). A size-segregated particle
!!     dry deposition scheme for an atmospheric aerosol module. Atmospheric environment.
!! (3) Most of the codes are adopted from GEOS-Chem drydep_mod.F90 module.
!!     https://github.com/geoschem/geos-chem
!!
!! \author Wei Li
!! \date 02/2025
!!!>
module CCPr_Scheme_Wesely_Mod

   implicit none

   private

   public :: CCPr_Scheme_Wesely

contains

   !>
   !! \brief Computes the dry deposition velocity using the Wesely scheme
   !!
   !!References: Wesely, M. L. (1989).
   !!
   !! \param RADIAT      Solar radiation [W/m2]
   !! \param TEMP        Surface Temperature [K]
   !! \param SUNCOS      Cosine of solar zenith angle at middle of current chem timestep
   !! \param F0          React. factor for oxidation depends on species
   !! \param HSTAR       Henry's law constant depends on species
   !! \param XMW         Molecular weight [kg/mol]
   !! \param A_RADI      Aerosol radius [m]
   !! \param A_DEN       Aerosol density [kg/m3]
   !! \param USTAR       Friction velocity [m/s]
   !! \param OBK         Monin-Obhukov length [m]
   !! \param CFRAC       Surface cloud fraction
   !! \param ZH          PBL height [m]
   !! \param THIK        height of first model layer [m]
   !! \param ZO          Roughness length [m]
   !! \param RHB         Relative humidity at surface [uniteless]
   !! \param PRESSU      Surface pressure [Pa]
   !! \param W10         Wind speed at 10m [m/s]
   !! \param SPC         Species name
   !! \param XLAI        Leaf area index (Note: change to fraction LAI of each land type)
   !! \param ILAND       Land type ID in current grid box (mapped to deposition surface types
   !! \param IUSE        Fraction of gridbox area occupied by each land type
   !! \param SALINITY    Salinity of the ocean
   !! \param TSKIN       Skin temperature
   !! \param IODIDE      Iodide concentration
   !! \param XLON        Longitude
   !! \param YLAT        Latitude
   !! \param IS_GAS      Flag for gas species
   !! \param IS_DUST     Flag for dust species
   !! \param IS_SEASALT  Flag for sea salt species
   !! \param IS_SNOW     Flag for snow surface
   !! \param IS_ICE      Flag for ice surface
   !! \param IS_LAND     Flag for land surface
   !! \param DD_DvzAerSnow  Fixed VD for some aerosols over snow and ice [cm/s]
   !! \param DD_DvzMinVal_SNOW  Minimum VD for some sulfate species over snow and ice [cm/s]
   !! \param DD_DvzMinVal_LAND  Minimum VD for some sulfate species over land [cm/s]
   !! \param VD          output of dry deposition velocity [m/s]
   !! \param DDFreq      output of dry deposition frequency [1/s]
   !! \param RC          Success or failure?
   !!
   !! \ingroup catchem_drydep_process
   !!!>
   subroutine CCPr_Scheme_Wesely(  RADIAT,    TEMP,       SUNCOS,                &
      F0,        HSTAR,     XMW,    A_RADI, A_DEN,        &
      USTAR,     OBK,        CFRAC,  ZH,  THIK,  &
      ZO,        RHB,        PRESSU,     &
      W10,       SPC, XLAI,  ILAND, IUSE, &
      SALINITY, TSKIN, IODIDE, XLON, YLAT, CO2_EFFECT, CO2_LEVEL, CO2_REF, LNLPBL, &
      IS_GAS, IS_DUST, IS_SEASALT, IS_SNOW, IS_ICE, IS_LAND, DD_DvzAerSnow, DD_DvzMinVal_SNOW, DD_DvzMinVal_LAND, VD, DDFreq, RC)
      ! Uses
      !USE Constants,     Only : PI_180      !pull in a constant from the CONSTANTS MODULE
      use precision_mod, only : fp           !pull in a precision from the PRECISION MODULE
      Use Error_Mod,     Only : CC_SUCCESS   ! Error Check Success
      USE CCPr_Drydep_Common_Mod

      IMPLICIT NONE
      ! Parameters
      !-----------
      real(fp), intent(in)  :: RADIAT      !< Solar radiation [W/m2]
      real(fp), intent(in)  :: TEMP        !< Temperature [K]
      real(fp), intent(in)  :: SUNCOS      !< Cosine of solar zenith angle at middle of current chem timestep
      real(fp), intent(inout)  :: F0          !< React. factor for oxidation depends on species
      real(fp), intent(in)  :: HSTAR       !< Henry's law constant depends on species
      real(fp), intent(in)  :: XMW         !< Molecular weight [kg/mol]
      real(fp), intent(in)  :: A_RADI      !< Aerosol radius [m]
      real(fp), intent(in)  :: A_DEN       !< Aerosol density [kg/m3]
      real(fp), intent(in)  :: USTAR       !< Friction velocity [m/s]
      real(fp), intent(in)  :: OBK         !< Monin-Obhukov length [m]
      real(fp), intent(in)  :: CFRAC       !< Surface cloud fraction [unitless]
      real(fp), intent(in)  :: ZH          !< PBL height [m]
      real(fp), intent(inout)  :: THIK        !< height of first model layer [m]
      real(fp), intent(in)  :: ZO          !< Roughness length [m]
      real(fp), intent(in)  :: RHB         !< Relative humidity at surface [uniteless]
      real(fp), intent(in)  :: PRESSU      !< Surface pressure [Pa]
      real(fp), intent(in)  :: W10         !< Wind speed at 10m [m/s]
      !integer,  intent(in)  :: N_SPC      !< Species ID (TODO: may be changed to species name)
      character(len=20), intent(in) :: SPC !< Species name
      real(fp), dimension(:), intent(in)  :: XLAI        !< Leaf area index (Note: change to fraction LAI of each land type)
      integer,  dimension(:), intent(in)  :: ILAND       !< Land type ID in current grid box (mapped to deposition surface types
      real(fp), dimension(:), intent(in)  :: IUSE        !< Fraction (per mille) of gridbox area occupied by each land type (TODO!!)
      !some inputs are for O3 over water and Hg over Amazon forest (not sure if we should include them for now)
      real(fp), intent(in)  :: SALINITY    !< Salinity of the ocean
      real(fp), intent(in)  :: TSKIN       !< Skin temperature
      real(fp), intent(in)  :: IODIDE      !< Iodide concentration
      real(fp), intent(in)  :: XLON        !< Longitude
      real(fp), intent(in)  :: YLAT        !< Latitude
      ! CO2 effect on Rs
      logical, intent(in)   :: CO2_EFFECT  !< Flag for CO2 effect on Rs
      real(fp), intent(in)  :: CO2_LEVEL   !< CO2 level
      real(fp), intent(in)  :: CO2_REF     !< Reference CO2 level
      logical, intent(in)   :: LNLPBL      !< Flag for non-local PBL mixing instead of full PBL mixing
      logical, intent(in)   :: IS_GAS, IS_DUST, IS_SEASALT
      logical, intent(in)   :: IS_SNOW, IS_ICE, IS_LAND !< Flags for snow, ice or land
      !set range of dry deposition velocities
      real(fp), intent(in)  :: DD_DvzAerSnow !< Fixed VD for some aerosols over snow and ice [cm/s]
      real(fp), intent(in)  :: DD_DvzMinVal_SNOW !< Minimum VD for some sulfate species over snow and ice [cm/s]
      real(fp), intent(in)  :: DD_DvzMinVal_LAND !< Minimum VD for some sulfate species over land [cm/s]
      !set scaling factor for dry deposition velocity if needed
      !real(fp), intent(in)  :: VD_scaling !< Scaling factor for dry deposition velocity
      !output
      real(fp), intent(out) :: VD          !< dry deposition velocity [m/s]
      real(fp), intent(out) :: DDFreq      !< dry deposition frequency [1/s]
      integer, intent(out)  :: RC          !< Success or failure?

      ! Local Variables
      !----------------
      real(fp) :: XLAI_IN, C1X, RA, RB, RSURFC, VTSoutput, VK, DVZ
      integer  :: II     !< Index of the drydep land type
      integer  :: ILDT   !< index of the land types in the grid box
      integer  :: LDT    !loop index of land types
      !string
      character(len=255)       :: thisLoc
      character(len=512)       :: ErrMsg

      !--------------------------------------------
      ! main function
      !--------------------------------------------

      ! Assume success
      RC      =  CC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> at CCPr_scheme_Wesely (in process/drydep/CCPr_Scheme_Wesely_Mod.F90)'

      ! Add option for non-local PBL mixing scheme: THIK must
      ! be the first box height.
      ! We do not use PBL_DRYDEP as the original codes
      IF (.NOT. LNLPBL) THIK = MAX( ZH, THIK )

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
      XLAI_IN    = 0.0_fp

      ! Better test for depositing species: We need both HSTAR and XMW
      ! to be nonzero, OR the value of AIROSOL to be true.  This should
      ! avoid any further floating point invalid issues caused by putting
      ! a zero value in a denominator.
      IF ( ( HSTAR > 0e+0_fp .and. XMW > 0e+0_fp ) .or. ( .NOT. IS_GAS ) ) THEN
         DO LDT =1 , SIZE(IUSE)
            ! If the land type is not represented in grid
            ! box, then skip to the next land type
            IF ( IUSE(LDT) <= 0 ) CYCLE
            ! Olson land type index + 1
            !TODO: not always +1???
            ILDT = ILAND(LDT)+1
            ! Dry deposition land type index
            !TODO: this IDEP is predefined based on Olson
            II   = IDEP(ILDT)

            !LAI of the landtype in the subgrid
            !XLAI_IN = XLAI * DBLE(IUSE(LDT)) !TODO: may be able to calculate online if fraction LAI is not provided
            XLAI_IN = XLAI(LDT)

            !If the surface to be snow or ice;set II to 1 instead
            !We do not use II index to specify directly since IS_SNOW and IS_ICE are given at each grid not subgrid as ILAND
            IF( (IS_SNOW) .OR. (IS_ICE) ) II=1

            !get bulk surface resistances (Rs)
            IF (IS_GAS) THEN
               call Wesely_Rc_Gas( RADIAT,    TEMP,       SUNCOS,                &
                  F0,        HSTAR,     XMW,            &
                  USTAR,     CFRAC,  PRESSU,     &
                  XLAI_IN, II,  SPC, SALINITY, TSKIN, IODIDE, XLON, YLAT, &
                  CO2_EFFECT, CO2_LEVEL, CO2_REF, &
                  RSURFC,   RC)
            ELSE
               !Note to change pressure unit from Pa to kPa
               RSURFC = AERO_SFCRSII ( SPC, II, IS_DUST, IS_SEASALT, LUCINDEX_GC, A_RADI, A_DEN, &
                  PRESSU*1e-3_fp, TEMP, USTAR, RHB, W10, VTSoutput, RC)

            ENDIF
            if (RC /= CC_SUCCESS ) then
               errMsg = 'Error in getting bulk surface resistances (RSURFC)'
               CALL CC_Error( errMsg, RC, thisLoc )
               RETURN
            endif

            !*Set max and min values for bulk surface resistances
            RSURFC = MAX(1.e+0_fp, MIN(RSURFC,9999.e+0_fp))
            !*because of high resistance values, different rule applied for ocean ozone
            IF ((SPC .EQ. 'O3') .AND. (II .EQ. 11)) THEN
               RSURFC = MAX(1.e+0_fp, MIN(RSURFC,999999.e+0_fp))
            ENDIF
            ! Set Rc for strong acids (HNO3,HCl,HBr) to 1 s/m
            ! Ref. Jaegle et al. 2018, cf. Erisman,van Pul,Ayers 1994
            IF ( HSTAR .gt. 1.e+10_fp ) RSURFC= 1.e+0_fp

            !get Ra and Rb
            call Wesely_Ra_Rb(TEMP, PRESSU, XMW, USTAR, OBK, ZO, THIK, LNLPBL, IS_GAS, Ra, Rb,  RC)

            !get VD (TODO: IUSE is decimal not percent or permille as in GEOS-Chem)
            C1X = RSURFC + Ra + Rb
            VK = VD
            !VD = VK + DBLE( IUSE(LDT) ) / C1X !This seems to be useless in the original codes
            IF (IS_GAS) THEN
               !VD = VK + DBLE( IUSE(LDT) ) / C1X
               VD = VK + IUSE(LDT)  / C1X
            ELSE
               !VD = VK + DBLE( IUSE(LDT) ) / C1X + DBLE( IUSE(LDT) ) * VTSoutput
               VD = VK +  IUSE(LDT)  / C1X +  IUSE(LDT) * VTSoutput
            ENDIF
         END DO
      ENDIF

      !apply spectial treatment or scaling factor to Vd
      DVZ = VD *100.e+0_fp !m/s -- > cm/s

      ! Scale relative to specified species(Note:we do not use FLAG but match names instead)
      !TODO: We simply hardcode the scaling factor here
      !DVZ = DVZ * VD_scaling

      !IF ( FLAG(D) .eq. 1 )  THEN
      IF ((SPC .eq. 'N2O5') .or. (SPC .eq. 'HC187') ) THEN

         ! Scale species to HNO3 (MW_g = 63.012 g/mol)
         DVZ = DVZ * sqrt(63.01) / sqrt( XMW*1e3_fp )

         !ELSE IF ( FLAG(D) .eq. 2 ) THEN
      ELSE IF ((SPC .eq. 'MPAN') .or. (SPC .eq. 'PPN') .or. (SPC .eq. 'R4N2')) THEN

         ! Scale species to PAN (MW_g = 121.06 g/mol)
         DVZ = DVZ * sqrt(121.06) / sqrt( XMW*1e3_fp )

         !ELSE IF ( FLAG(D) .eq. 3 ) THEN
      ELSE IF ((SPC .eq. 'MONITS') .or. (SPC .eq. 'MONITU') .or. (SPC .eq. 'HONIT')) THEN

         ! Scale species to ISOPN (MW_g = 147.15 g/mol)
         DVZ = DVZ * sqrt(147.15)  / sqrt(XMW*1e3_fp)

      ENDIF

      !-----------------------------------------------------------
      ! Special treatment for snow and ice
      !-----------------------------------------------------------
      IF ( (IS_SNOW) .OR. (IS_ICE) ) THEN

         !-------------------------------------
         ! %%% SURFACE IS SNOW OR ICE %%%
         !-------------------------------------
         IF ( DD_DvzAerSnow > 0.0_fp ) THEN

            ! For most aerosol species (basically everything
            ! except sea salt and dust species), we just set
            ! the deposition velocity over snow to a fixed value
            !DVZ = DBLE( DD_DvzAerSnow )
            DVZ = DD_DvzAerSnow

         ELSE

            ! Otherwise, enforce a minimum drydep velocity over snow
            ! (cf. the GOCART model).  NOTE: In practice this will
            ! only apply to the species SO2, SO4, MSA, NH3, NH4, NIT.
            !DVZ = MAX( DVZ, DBLE( DD_DvzMinVal_Snow ) )
            DVZ = MAX( DVZ,  DD_DvzMinVal_Snow )

         ENDIF

      ELSE

         !-------------------------------------
         ! %%% SURFACE IS NOT SNOW OR ICE %%%
         !-------------------------------------

         ! Enforce a minimum drydep velocity over land (cf. the
         ! GOCART model).  NOTE: In practice this will only apply
         ! to the species SO2, SO4, MSA, NH3, NH4, NIT.
         !DVZ = MAX( DVZ, DBLE( DD_DvzMinVal_Land ) )
         DVZ = MAX( DVZ,  DD_DvzMinVal_Land )

      ENDIF

      !-----------------------------------------------------------
      ! Special treatment for ACETONE
      !-----------------------------------------------------------

      ! For ACET, we need to only do drydep over the land
      ! and not over the oceans.
      !IF ( N == id_ACET ) THEN
      IF ( SPC == 'ACET' ) THEN
         IF ( Is_Land ) THEN
            DVZ = 0.1e+0_fp
         ELSE
            DVZ = 0e+0_fp
         ENDIF
      ENDIF

      !-----------------------------------------------------------
      ! Special treatment for ALD2,MENO3,ETNO3,MOH
      !-----------------------------------------------------------

      ! we need to only do drydep over the land
      ! and not over the oceans.
      !IF ( N == id_ALD2 ) THEN
      IF ( (SPC == 'ALD2') .or. (SPC == 'MENO3') .or. (SPC == 'ETNO3') .or. (SPC == 'MOH') ) THEN
         IF ( .not. Is_Land ) THEN
            DVZ = 0e+0_fp
         ENDIF
      ENDIF

      !-----------------------------------------------------------
      ! Compute drydep velocity and frequency
      !-----------------------------------------------------------

      ! Dry deposition velocities [m/s]
      VD = DVZ / 100.e+0_fp

      ! Dry deposition frequency [1/s]
      DDFreq = VD / THIK

      !test only
      !write(*,*) 'Test finish for species () with Vd (): ', SPC, VD


   end subroutine CCPr_Scheme_Wesely

end module CCPr_Scheme_Wesely_Mod
