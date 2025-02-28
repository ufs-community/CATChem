!>
!! \file
!! \brief CCPr Scheme for dry deposition of aeroosl species from Zhang et al., (2001) with
!! Emerson's updates. The Ra and Rb are still from Wesely (1989) for now.
!!
!!
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
!! \author Wei Li
!! \date 02/2025
!!!>
module CCPr_Scheme_Zhang_aerosol_Mod

   implicit none

   private

   public :: CCPr_Scheme_Zhang_Aero

contains

   !>
   !! \brief Computes the dry deposition velocity using the Wesely scheme
   !!
   !!References: Wesely, M. L. (1989).
   !!
   !! \param TEMP        Surface Temperature [K]
   !! \param HSTAR       Henry's law constant depends on species
   !! \param XMW         Molecular weight [kg/mol]
   !! \param A_RADI      Aerosol radius [m]
   !! \param A_DEN       Aerosol density [kg/m3]
   !! \param USTAR       Friction velocity [m/s]
   !! \param OBK         Monin-Obhukov length [m]
   !! \param THIK        height of first model layer [m]
   !! \param ZO          Roughness length [m]
   !! \param RHB         Relative humidity at surface [uniteless]
   !! \param PRESSU      Surface pressure [Pa]
   !! \param W10         Wind speed at 10m [m/s]
   !! \param SPC         Species name
   !! \param ILAND       Land type ID in current grid box (mapped to deposition surface types
   !! \param IUSE        Fraction of gridbox area occupied by each land type
   !! \param SeaSalt_Lower_Bin Lower bin boundary of sea salt radius [um]
   !! \param SeaSalt_UPPER_Bin Upper bin boundary of sea salt radius [um]
   !! \param LUC         name of land use category (one of OLSON, NOAH and IGBP for now)
   !! \param IS_DUST     Flag for dust species
   !! \param IS_SEASALT  Flag for sea salt species
   !! \param IS_SNOW     Flag for snow surface
   !! \param IS_ICE      Flag for ice surface
   !! \param DD_DvzAerSnow  Fixed VD for some aerosols over snow and ice [cm/s]
   !! \param DD_DvzMinVal_SNOW  Minimum VD for some sulfate species over snow and ice [cm/s]
   !! \param DD_DvzMinVal_LAND  Minimum VD for some sulfate species over land [cm/s]
   !! \param VD          output of dry deposition velocity [m/s]
   !! \param DDFreq      output of dry deposition frequency [1/s]
   !! \param RC          Success or failure?
   !!
   !! \ingroup catchem_drydep_process
   !!!>
   subroutine CCPr_Scheme_Zhang_Aero( TEMP, HSTAR, XMW, A_RADI, A_DEN,                   &
      USTAR,  OBK, THIK,  ZO,  RHB,    PRESSU,  W10,    SPC, ILAND, IUSE,                &
      SeaSalt_Lower_Bin, SeaSalt_UPPER_Bin, LUC, IS_DUST, IS_SEASALT, IS_SNOW, IS_ICE,   &
      DD_DvzAerSnow, DD_DvzMinVal_SNOW, DD_DvzMinVal_LAND, VD, DDFreq, RC)
      ! Uses
      !USE Constants,     Only : PI_180      !pull in a constant from the CONSTANTS MODULE
      use precision_mod, only : fp           !pull in a precision from the PRECISION MODULE
      Use Error_Mod,     Only : CC_SUCCESS   ! Error Check Success
      USE CCPr_Drydep_Common_Mod

      IMPLICIT NONE
      ! Parameters
      !-----------
      !real(fp), intent(in)  :: RADIAT      !< Solar radiation [W/m2]
      real(fp), intent(in)  :: TEMP        !< Temperature [K]
      !real(fp), intent(in)  :: SUNCOS      !< Cosine of solar zenith angle at middle of current chem timestep
      !real(fp), intent(inout)  :: F0          !< React. factor for oxidation depends on species
      real(fp), intent(in)  :: HSTAR       !< Henry's law constant depends on species
      real(fp), intent(in)  :: XMW         !< Molecular weight [kg/mol]
      real(fp), intent(in)  :: A_RADI      !< Aerosol radius [m]
      real(fp), intent(in)  :: A_DEN       !< Aerosol density [kg/m3]
      real(fp), intent(in)  :: USTAR       !< Friction velocity [m/s]
      real(fp), intent(in)  :: OBK         !< Monin-Obhukov length [m]
      !real(fp), intent(in)  :: CFRAC       !< Surface cloud fraction [unitless]
      !real(fp), intent(in)  :: ZH          !< PBL height [m]
      real(fp), intent(in)  :: THIK        !< height of first model layer [m]
      real(fp), intent(in)  :: ZO          !< Roughness length [m]
      real(fp), intent(in)  :: RHB         !< Relative humidity at surface [uniteless]
      real(fp), intent(in)  :: PRESSU      !< Surface pressure [Pa]
      real(fp), intent(in)  :: W10         !< Wind speed at 10m [m/s]
      !integer,  intent(in)  :: N_SPC      !< Species ID (TODO: may be changed to species name)
      character(len=20), intent(in) :: SPC !< Species name
      !real(fp), dimension(:), intent(in)  :: XLAI        !< Leaf area index (Note: change to fraction LAI of each land type)
      integer,  dimension(:), intent(in)  :: ILAND       !< Land type ID in current grid box (mapped to deposition surface types
      real(fp), dimension(:), intent(in)  :: IUSE        !< Fraction (per mille) of gridbox area occupied by each land type (TODO!!)
      real(fp), dimension(:), intent(in)  :: SeaSalt_Lower_Bin !< Lower bin boundary of sea salt radius [um]
      real(fp), dimension(:), intent(in)  :: SeaSalt_UPPER_Bin !< Upper bin boundary of sea salt radius [um]
      !some inputs are for O3 over water and Hg over Amazon forest (not sure if we should include them for now)
      !real(fp), intent(in)  :: SALINITY    !< Salinity of the ocean
      !real(fp), intent(in)  :: TSKIN       !< Skin temperature
      !real(fp), intent(in)  :: IODIDE      !< Iodide concentration
      !real(fp), intent(in)  :: XLON        !< Longitude
      !real(fp), intent(in)  :: YLAT        !< Latitude
      character(len=20), intent(in) :: LUC !< name of land use category (one of OLSON, NOAH and IGBP for now)
      ! CO2 effect on Rs
      !logical, intent(in)   :: CO2_EFFECT  !< Flag for CO2 effect on Rs
      !real(fp), intent(in)  :: CO2_LEVEL   !< CO2 level
      !real(fp), intent(in)  :: CO2_REF     !< Reference CO2 level
      logical, intent(in)   :: IS_DUST, IS_SEASALT
      logical, intent(in)   :: IS_SNOW, IS_ICE !< Flags for snow, ice
      !set range of dry deposition velocities
      real(fp), intent(in)  :: DD_DvzAerSnow !< Fixed VD for some aerosols over snow and ice [cm/s]
      real(fp), intent(in)  :: DD_DvzMinVal_SNOW !< Minimum VD for some sulfate species over snow and ice [cm/s]
      real(fp), intent(in)  :: DD_DvzMinVal_LAND !< Minimum VD for some sulfate species over land [cm/s]
      !output
      real(fp), intent(out) :: VD          !< dry deposition velocity [m/s]
      real(fp), intent(out) :: DDFreq      !< dry deposition frequency [1/s]
      integer, intent(out)  :: RC          !< Success or failure?

      ! Local Variables
      !----------------
      real(fp) :: C1X, RA, RB, RSURFC, VTSoutput, VK, DVZ
      integer  :: II     !< Index of the drydep land type
      integer  :: ILDT   !< index of the land types in the grid box
      integer  :: LDT    !loop index of land types
      integer  :: LUCINDEX !mapping above II to Zhang's 15 land types for aerosols
      !string
      character(len=255)       :: thisLoc
      character(len=512)       :: ErrMsg

      !--------------------------------------------
      ! main function
      !--------------------------------------------

      ! Assume success
      RC      =  CC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> at CCPr_scheme_Zhang_Aero (in process/drydep/CCPr_Scheme_Zhang_aerosol_Mod.F90)'

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

      ! Better test for depositing species: We need both HSTAR and XMW
      ! to be nonzero, OR the value of AIROSOL to be true.  This should
      ! avoid any further floating point invalid issues caused by putting
      ! a zero value in a denominator.
      DO LDT =1 , SIZE(IUSE)
         ! If the land type is not represented in grid
         ! box, then skip to the next land type
         IF ( IUSE(LDT) <= 0 ) CYCLE

         ILDT = ILAND(LDT)
         IF ( LUC == 'OLSON' ) THEN
            ! Olson land type index + 1
            ILDT = ILDT + 1
            ! Dry deposition land type index
            II   = IDEP_IOLSON(ILDT)
            LUCINDEX = LUCINDEX_GC(II)
         ELSE IF ( LUC == 'NOAH' ) THEN
            ! it is possible that water is given as 0 not 17 in GFS CCPP
            IF (ILDT == 0) ILDT = 17
            II   = IDEP_NOAH(ILDT)
            !Note: we use ILDT, instead of II,  to get LUCINDEX here
            LUCINDEX = LUCINDEX_NOAH(ILDT)
         ELSE IF ( LUC == 'IGBP' ) THEN
            ! it is possible that water is given as 0 not 17
            IF (ILDT == 0) ILDT = 17
            II   = IDEP_IGBP(ILDT)
            LUCINDEX = LUCINDEX_IGBP(ILDT)
         ENDIF

         !get bulk surface resistances (Rs)
         !Note to change pressure unit from Pa to kPa
         RSURFC = AERO_SFCRSII ( SPC, IS_DUST, IS_SEASALT, LUCINDEX, A_RADI, A_DEN, PRESSU*1e-3_fp, &
            TEMP, USTAR, RHB, W10, SeaSalt_Lower_Bin, SeaSalt_UPPER_Bin,VTSoutput, RC)

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
         call Wesely_Ra_Rb(TEMP, PRESSU, XMW, USTAR, OBK, ZO, THIK, .FALSE., Ra, Rb,  RC)

         !get VD (TODO: IUSE is decimal not percent or permille as in GEOS-Chem)
         C1X = RSURFC + Ra + Rb
         VK = VD
         !VD = VK + DBLE( IUSE(LDT) ) / C1X + DBLE( IUSE(LDT) ) * VTSoutput
         VD = VK +  IUSE(LDT)  / C1X +  IUSE(LDT) * VTSoutput
      END DO


      !apply spectial treatment or scaling factor to Vd
      DVZ = VD *100.e+0_fp !m/s -- > cm/s

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
      ! Compute drydep velocity and frequency
      !-----------------------------------------------------------

      ! Dry deposition velocities [m/s]
      VD = DVZ / 100.e+0_fp

      ! Dry deposition frequency [1/s]
      DDFreq = VD / THIK

      !test only
      !write(*,*) 'Test finish for species () with Vd (): ', SPC, VD


   end subroutine CCPr_Scheme_Zhang_Aero


end module CCPr_Scheme_Zhang_aerosol_Mod
