!>
!! \file ccpr_bvoc_common_mod.F90
!! \brief Contains module ccpr_bvoc_common_mod
!!
!! \ingroup catchem_bvoc_process
!!
!! \author Barry Baker
!! \date 05/2024
!!!>
module CCPr_BVOC_Common_Mod
   use precision_mod, only: fp, ZERO, rae
   use Error_Mod
   use constants
   implicit none
   private

   public  :: GET_MEGAN_PARAMS
   public  :: GET_GAMMA_PAR_PCEEA
   public  :: GET_GAMMA_PAR_C
   public  :: GET_GAMMA_T_LI
   public  :: GET_GAMMA_T_LD
   public  :: GET_GAMMA_T_LD_C
   public  :: GET_GAMMA_LAI
   public  :: GET_GAMMA_AGE
   public  :: GET_GAMMA_SM
   public  :: CALC_NORM_FAC
   public  :: Calc_Sun_Frac
   public  :: SOLAR_ANGLE
   public  :: GET_GAMMA_CO2
   public  :: CALC_AEF
   public  :: GET_CDEA
   public  :: BvocStateType

   !> \brief Type for CATCHem BVOC Process
   !!
   !! \details Contains all the information needed to run the CATChem BVOC Process
   !!
   !! This type contains the following variables:
   !! - Activate : Activate Process (True/False)
   !! - nBvocSpecies : Number of BVOC processes
   !! - BvocSpeciesIndex : Index of BVOC species
   !! - SpcIDs : CATChem species IDs
   !! - CO2Inhib      : CO2 inhibition for isoprene Option [true/false]
   !! - CO2conc       : CO2 concentration [ppmv]
   !! - ISOPscale     : factors to scale isoprene emissions
   !! - ISOPtoSOAP    : isoprene conversion factor to SOAP
   !! - ISOPtoSOAS    : isoprene conversion factor to SOAS
   !! - MONOtoSOAP    : monoterpene conversion factor to SOAP
   !! - MONOtoSOAS    : monoterpene conversion factor to SOAS
   !! - TERPtoSOAP    : other terpene conversion factor to SOAP
   !! - TERPtoSOAS    : other terpene conversion factor to SOAS
   !! - TotalEmission : Total Emission          [kg/m^2/s]
   !! - EmissionPerSpecies : Emission Rate per Bvoc species  [kg/m^2/s]
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   TYPE :: BvocStateType
      ! Generic Variables for Every Process
      Logical                         :: Activate                !< Activate Process (True/False)
      INTEGER                         :: SchemeOpt               !< Scheme Option
      integer                         :: nBvocSpecies           !< Number of BVOC processes
      integer, allocatable            :: BvocSpeciesIndex(:)    !< Index of BVOC species
      character(len=31), allocatable  :: BvocSpeciesName(:)     !< name of BVOC species
      integer, allocatable            :: SpcIDs(:)               !< CATChem species IDs
      integer                         :: CatIndex                !< Index of emission category in EmisState

      ! Process Specific Parameters
      real(fp)                        :: TotalEmission           !< Total emission          [kg/m^2/s]
      real(fp), allocatable           :: EmissionPerSpecies(:)   !< Emission per species    [kg/m^2/s]
      real(fp), allocatable           :: EmisNormFactor(:)       !< Emission normalized factor (onle one used for all now)

      ! Scheme Options (ISOP scaling is turned off at the moment)
      Logical                         :: CO2Inhib                !< CO2 inhibition for isoprene Option [True/False]
      real(fp)                        :: CO2conc                 !< CO2 concentration [ppmv]
      !real(fp)                        :: ISOPscale               !< factors to scale isoprene emissions
      !real(fp)                        :: ISOPtoSOAP              !< isoprene conversion factor to SOAP
      !real(fp)                        :: ISOPtoSOAS              !< isoprene conversion factor to SOAS
      !real(fp)                        :: MONOtoSOAP              !< monoterpene conversion factor to SOAP
      !real(fp)                        :: MONOtoSOAS              !< monoterpene conversion factor to SOAS
      !real(fp)                        :: TERPtoSOAP              !< other terpenes conversion factor to SOAP
      !real(fp)                        :: TERPtoSOAS              !< other terpenes conversion factor to SOAS


      !=================================================================
      ! Module specific variables/arrays/data pointers come below
      !=================================================================

   END TYPE BvocStateType

contains
   !>
   !! \brief Computes the light-independent fraction of emissions
   !!
   !!References:
   !! new function from below Sam et al, 2020
   !! Silva, S. J., Heald, C. L., and Guenther, A. B.: Development of a reduced-complexity
   !! plant canopy physics surrogate model for use in chemical transport models: a case study
   !! with GEOS-Chem v12.3.0, Geosci. Model Dev., 13, 2569–2585,
   !! https://doi.org/10.5194/gmd-13-2569-2020, 2020.
   !!
   !! \param LAI leaf area index
   !! \param Sinbeta
   !! \param Distgauss   Gauss distance
   !! \param SunFrac  output of light-dependent emission factor
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   subroutine Calc_Sun_Frac( LAI, Sinbeta, Distgauss, SunFrac)
      IMPLICIT NONE
      ! Parameters
      !-----------
      real(fp), intent(in)  :: LAI         !< leaf area index
      real(fp), intent(in)  :: Sinbeta     !<
      real(fp), intent(in)  :: Distgauss   !<
      real(fp), intent(out) :: SunFrac     !< Activity factor for the light-independent fraction of emissions

      ! Local Variables
      !----------------
      real(fp), parameter :: Cluster = 0.9 !< Standard reference temperature [K]
      real(fp), parameter :: CANTRAN = 0.2 !<
      real(fp) :: Kb, LAIadj, LAIdepth     !<

      !--------------------------------------------
      ! main function
      !--------------------------------------------
      Kb = Cluster * 0.5 / Sinbeta
      LAIadj = LAI / ( 1 - CANTRAN )
      LAIdepth   = LAIadj  * Distgauss

      if ((Sinbeta  > 0.002) .AND. (LAIadj  > 0.001)) then
         SunFrac = EXP(-Kb * LAIdepth)
      else
         SunFrac = 0.2
      endif
      return

   end subroutine Calc_Sun_Frac

   !>
   !! \returns the emission parameters for each MEGAN compound needed to compute emissions.
   !!
   !! Guenther et al, (GMD 2012) and associated MEGANv2.1 source code
   !!
   !! \param CPD,
   !! \param BTA,   LIDF,  C_T1,  C_EO, A_NEW, A_GRO, A_MAT, A_OLD, BI_DIR
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   subroutine GET_MEGAN_PARAMS( CPD,   BTA,   LIDF,  C_T1,  C_EO, A_NEW, A_GRO, A_MAT, A_OLD, BI_DIR, RC)

      ! Uses
      use precision_mod, only : fp
      use Error_Mod,     Only : CC_FAILURE, CC_Error
      IMPLICIT NONE
      ! input Parameters
      character(LEN=256), intent(in) :: CPD       ! Compound name
      ! input/output Parameters
      real(fp), intent(inout)  :: BTA      !< Beta coefficient for temperature activity factor for light-independent fraction

      real(fp), intent(inout)  :: LIDF     !< Light-dependent fraction of emissions
      real(fp), intent(inout)  :: C_T1     !< CT1 parameter for temperature activity factor for light-dependent fraction
      real(fp), intent(inout)  :: C_EO     !< Ceo parameter for temperature activity factor for light-dependent fraction
      real(fp), intent(inout)  :: A_NEW    !< Relative emission factor (new leaves)
      real(fp), intent(inout)  :: A_GRO    !< Relative emission factor (growing leaves)
      real(fp), intent(inout)  :: A_MAT    !< Relative emission factor (mature leaves)
      real(fp), intent(inout)  :: A_OLD    !< Relative emission factor (old leaves)
      logical,  intent(inout)  :: BI_DIR   !< Logical flag to indicate bidirectional exchange
      integer,  intent(out)    :: RC       !< Success or Failure

      !local variables
      character(len=255)       :: MSG, thisLoc

      !=================================================================
      ! GET_MEGAN_PARAMS begins here!
      !=================================================================

      ! Initialize values
      BTA    = 0.0_fp
      LIDF   = 0.0_fp
      C_T1   = 0.0_fp
      C_EO   = 0.0_fp
      A_NEW  = 0.0_fp
      A_GRO  = 0.0_fp
      A_MAT  = 0.0_fp
      A_OLD  = 0.0_fp
      BI_DIR = .FALSE.

      ! ----------------------------------------------------------------
      ! Note that not all the above compounds are used in standard chemistry
      ! simulations, but they are provided here for future incorporation or
      ! specialized applications. More compounds can be added as needed
      ! by adding the corresponding CPD name and the appropriate parameters.
      !
      ! Values are from Table 4 in Guenther et al., 2012
      ! ----------------------------------------------------------------

      ! Isoprene, MBO
      IF ( TRIM(CPD) == 'ISOP' .OR. &
         TRIM(CPD) == 'MBOX' ) THEN
         BTA    = 0.13_fp  ! Not actually used for ISOP, MBO
         LIDF   = 1.0_fp
         C_T1   = 95.0_fp
         C_EO   = 2.0_fp
         A_NEW  = 0.05_fp
         A_GRO  = 0.6_fp
         A_MAT  = 1.0_fp
         A_OLD  = 0.9_fp
         BI_DIR = .FALSE.

         ! Myrcene, sabinene, alpha-pinene
      ELSE IF ( TRIM(CPD) == 'MYRC' .OR. &
         TRIM(CPD) == 'SABI' .OR. &
         TRIM(CPD) == 'APIN' ) THEN
         BTA    = 0.10_fp
         LIDF   = 0.6_fp
         C_T1   = 80.0_fp
         C_EO   = 1.83_fp
         A_NEW  = 2.0_fp
         A_GRO  = 1.8_fp
         A_MAT  = 1.0_fp
         A_OLD  = 1.05_fp
         BI_DIR = .FALSE.

         ! Limonene, 3-carene, beta-pinene
      ELSE IF ( TRIM(CPD) == 'LIMO' .OR. &
         TRIM(CPD) == 'CARE' .OR. &
         TRIM(CPD) == 'BPIN' ) THEN
         BTA    = 0.10_fp
         LIDF   = 0.2_fp
         C_T1   = 80.0_fp
         C_EO   = 1.83_fp
         A_NEW  = 2.0_fp
         A_GRO  = 1.8_fp
         A_MAT  = 1.0_fp
         A_OLD  = 1.05_fp
         BI_DIR = .FALSE.

         ! t-beta-ocimene
      ELSE IF ( TRIM(CPD) == 'OCIM' ) THEN
         BTA    = 0.10_fp
         LIDF   = 0.8_fp
         C_T1   = 80.0_fp
         C_EO   = 1.83_fp
         A_NEW  = 2.0_fp
         A_GRO  = 1.8_fp
         A_MAT  = 1.0_fp
         A_OLD  = 1.05_fp
         BI_DIR = .FALSE.

         ! Other monoterpenes (lumped)
      ELSE IF ( TRIM(CPD) == 'OMON' ) THEN
         BTA    = 0.10_fp
         LIDF   = 0.4_fp
         C_T1   = 80.0_fp
         C_EO   = 1.83_fp
         A_NEW  = 2.0_fp
         A_GRO  = 1.8_fp
         A_MAT  = 1.0_fp
         A_OLD  = 1.05_fp
         BI_DIR = .FALSE.

         ! Methanol
      ELSE IF ( TRIM(CPD) == 'MOH' ) THEN
         BTA    = 0.08_fp
         LIDF   = 0.8_fp
         C_T1   = 60.0_fp
         C_EO   = 1.6_fp
         A_NEW  = 3.5_fp
         A_GRO  = 3.0_fp
         A_MAT  = 1.0_fp
         A_OLD  = 1.2_fp
         BI_DIR = .FALSE.

         ! Acetone
      ELSE IF ( TRIM(CPD) == 'ACET' ) THEN
         BTA    = 0.1_fp
         LIDF   = 0.2_fp
         C_T1   = 80.0_fp
         C_EO   = 1.83_fp
         A_NEW  = 1.0_fp
         A_GRO  = 1.0_fp
         A_MAT  = 1.0_fp
         A_OLD  = 1.0_fp
         BI_DIR = .FALSE.

         ! Bidirectional VOC: Ethanol, formaldehyde, acetaldehyde, formic acid,
         ! acetic acid
      ELSE IF ( TRIM(CPD) == 'EOH'  .OR. &
         TRIM(CPD) == 'CH2O' .OR. &
         TRIM(CPD) == 'ALD2' .OR. &
         TRIM(CPD) == 'FAXX' .OR. &
         TRIM(CPD) == 'AAXX' ) THEN
         BTA    = 0.13_fp
         LIDF   = 0.8_fp
         C_T1   = 95.0_fp
         C_EO   = 2.0_fp
         A_NEW  = 1.0_fp
         A_GRO  = 1.0_fp
         A_MAT  = 1.0_fp
         A_OLD  = 1.0_fp
         BI_DIR = .TRUE.

         ! Stress VOCs: ethene, toluene, HCN
         ! There are others species in this category but none are currently
         ! used in GEOS-Chem
      ELSE IF ( TRIM(CPD) == 'C2H4' .OR. &
         TRIM(CPD) == 'TOLU' .OR. &
         TRIM(CPD) == 'HCNX' ) THEN
         BTA    = 0.1_fp
         LIDF   = 0.8_fp
         C_T1   = 80.0_fp
         C_EO   = 1.83_fp
         A_NEW  = 1.0_fp
         A_GRO  = 1.0_fp
         A_MAT  = 1.0_fp
         A_OLD  = 1.0_fp
         BI_DIR = .FALSE.

         ! Other VOCs: >C2 alkenes
         ! This includes propene, butene and very minor contribution from
         ! larger alkenes
      ELSE IF ( TRIM(CPD) == 'PRPE' ) THEN
         BTA    = 0.1_fp
         LIDF   = 0.2_fp
         C_T1   = 80.0_fp
         C_EO   = 1.83_fp
         A_NEW  = 1.0_fp
         A_GRO  = 1.0_fp
         A_MAT  = 1.0_fp
         A_OLD  = 1.0_fp
         BI_DIR = .FALSE.

         ! SOAupdate: Sesquiterpenes hotp 3/2/10
         ! alpha-Farnesene, beta-Caryophyllene, other sesquiterpenes
      ELSE IF ( TRIM(CPD) == 'FARN' .OR. &
         TRIM(CPD) == 'BCAR' .OR. &
         TRIM(CPD) == 'OSQT' ) THEN
         BTA    = 0.17_fp
         LIDF   = 0.5_fp
         C_T1   = 130.0_fp
         C_EO   = 2.37_fp
         A_NEW  = 0.4_fp
         A_GRO  = 0.6_fp
         A_MAT  = 1.0_fp
         A_OLD  = 0.95_fp
         BI_DIR = .FALSE.

         ! Calls for any other MEGAN compounds (e.g. sesquiterpenes, etc.) can
         ! be added following the above format based on the parameters in
         ! Guenther 2012 or the MEGAN source code.
      ELSE
         RC = CC_FAILURE
         MSG = 'Invalid compound name'
         thisLoc = ' -> at CCPr_BVOC_Common (in process/bvoc/ccpr_bvoc_common_mod.F90)'
         call CC_Error( MSG, RC , thisLoc)
         return

      ENDIF
      return

   end subroutine GET_MEGAN_PARAMS

   !>
   !! \brief Computes the PCEEA gamma activity factor with sensitivity to light
   !!
   !! References:
   !! Guenther et al, 2006; Guenther et al, 2007, MEGAN v2.1 user guide
   !! Code was taken & adapted directly from the MEGAN v2.1 source code.
   !!
   !! \param Q_DIR_2, Q_DIFF_2
   !! \param PARDR_AVG_SIM, PARDF_AVG_SIM
   !! \param LAT, DOY, LocalHour
   !! \param D2RAD, RAD2D
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   subroutine GET_GAMMA_PAR_PCEEA(Q_DIR_2, Q_DIFF_2, PARDR_AVG_SIM, PARDF_AVG_SIM,  &
      LAT, DOY, LocalHour, D2RAD, RAD2D, GAMMA_P_PCEEA)

      IMPLICIT NONE
      ! Parameters
      real(fp),  intent(in)  :: Q_DIR_2           !< Direct PAR [W/m2]
      real(fp),  intent(in)  :: Q_DIFF_2          !< Diffuse PAR [W/m2]
      real(fp),  intent(in)  :: PARDR_AVG_SIM     !< Avg direct PAR [W/m2]
      real(fp),  intent(in)  :: PARDF_AVG_SIM     !< Avg diffuse PAR [W/m2]
      real(fp),  intent(in)  :: LAT               !< Note: this may need a new function and put in local variable
      integer,   intent(in)  :: DOY               !< Note: this may need a new function and put in local variable
      real(fp),  intent(in)  :: LocalHour         !< Note: this may need a new function and put in local variable
      real(fp),  intent(in)  :: D2RAD, RAD2D      !< Note: this could be put in constants
      real(fp),  intent(out) :: GAMMA_P_PCEEA     !< GAMMA factor for light

      ! Local Variables
      real(fp) :: mmPARDR_DAILY      !<
      real(fp) :: mmPARDF_DAILY      !<
      real(fp) :: PAC_DAILY, PAC_INSTANT, C_PPFD
      real(fp) :: PTOA, PHI
      real(fp) :: BETA,   SINbeta
      real(fp) :: AAA, BBB
      !real(fp) :: LocalHour, LAT
      !integer  :: DOY
      !integer  :: RC

      ! Constants
      !real(fp), parameter :: mmd = 3.4         !< median mass diameter [microns]
      ! W/m2 -> umol/m2/s
      REAL(fp), PARAMETER  :: WM2_TO_UMOLM2S = 4.766_fp

      !-----------------------------------------------------
      ! Compute GAMMA_PAR_PCEEA
      !-----------------------------------------------------
      ! Initialize
      C_PPFD = 0.0_fp
      PTOA   = 0.0_fp

      ! Convert past light conditions to micromol/m2/s
      mmPARDR_DAILY = PARDR_AVG_SIM  * WM2_TO_UMOLM2S
      mmPARDF_DAILY = PARDF_AVG_SIM  * WM2_TO_UMOLM2S

      ! Work out the light at the top of the canopy.
      PAC_DAILY    = mmPARDR_DAILY + mmPARDF_DAILY
      PAC_INSTANT  = (Q_DIR_2  +  Q_DIFF_2) * WM2_TO_UMOLM2S

      ! Get latitude
      !LAT = HcoState%Grid%YMID%Val(I,J)

      ! Get day of year, local-time and latitude
      !CALL HcoClock_Get( HcoState%Clock, cDOY = DOY, RC=RC )
      !CALL HcoClock_GetLocal( HcoState, I, J, cH = LocalHour, RC=RC )

      ! Get solar elevation angle
      CALL  SOLAR_ANGLE( DOY, LocalHour, LAT, D2RAD, SINbeta )
      BETA    =  ASIN( SINbeta ) * RAD2D

      IF ( SINbeta < 0.0_fp ) THEN

         GAMMA_P_PCEEA = 0.0_fp

      ELSEIF ( SINbeta > 0.0_fp ) THEN

         ! PPFD at top of atmosphere
         PTOA    = 3000.0_fp + 99.0_fp * &
            COS( 2._fp * 3.14159265358979323_fp * &
            ( DOY - 10.0_fp ) / 365.0_fp )

         ! Above canopy transmission
         PHI     = PAC_INSTANT / ( SINbeta * PTOA )

         ! Work out gamma P
         BBB     = 1.0_fp + 0.0005_fp *( PAC_DAILY - 400.0_fp )
         AAA     = ( 2.46_fp * BBB * PHI ) - ( 0.9_fp * PHI**2 )

         GAMMA_P_PCEEA = SINbeta * AAA

      ENDIF

      ! Screen unforced errors. IF solar elevation angle is
      ! less than 1 THEN gamma_p can not be greater than 0.1.
      IF ( BETA < 1.0_fp .AND. GAMMA_P_PCEEA > 0.1_fp ) THEN
         GAMMA_P_PCEEA  = 0.0_fp
      ENDIF

      ! Prevent negative values
      GAMMA_P_PCEEA = MAX( GAMMA_P_PCEEA , 0.0_fp )

      return
   end subroutine GET_GAMMA_PAR_PCEEA

   !>
   !! \brief Computes the local solar angle for a given day of year, latitude and longitude (or local time).
   !!
   !! References:
   !! (1 ) Guenther et al, 2006
   !! (2 ) Guenther et al, MEGAN v2.1 user manual 2007-09
   !! This code was taken directly from the MEGAN v2.1 source code
   !!
   !! \param DOY, SHOUR
   !! \param LAT
   !! \param D2RAD
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   subroutine SOLAR_ANGLE(DOY, SHOUR, LAT, D2RAD, SINbeta)

      IMPLICIT NONE
      ! Parameters
      integer,  intent(in)    :: DOY      !< Day of year
      real(fp), intent(in)    :: SHOUR    !< Local time
      real(fp), intent(in)    :: LAT      !< Latitude
      real(fp), intent(in)    :: D2RAD    !< Degree to radiance
      real(fp), intent(out)   :: SINbeta  !< Sin of the local solar angle


      ! local variable
      real(fp) :: sindelta, cosdelta, A, B

      ! Calculation of sin beta
      sindelta = -SIN( 0.40907_fp ) * COS( 6.28_fp * ( DOY + 10 ) / 365 )

      cosdelta = (1-sindelta**2.0_fp)**0.5_fp

      A = SIN( LAT * D2RAD ) * sindelta
      B = COS( LAT * D2RAD ) * cosdelta

      SINbeta = A + B * COS( 2.0_fp * PI * ( SHOUR-12 )/24 )

      return
   end subroutine SOLAR_ANGLE

   !>
   !! \brief Computes the temperature activity factor for the light-independent fraction of emissions
   !!
   !!References:
   !! (1 ) Guenther et al, 2006
   !! (2 ) Guenther et al, MEGAN user manual 2007-08
   !!(3 ) Guenther et al., GMD 2012 and MEGANv2.1 source code
   !!
   !! \param T
   !! \param T_Leaf_Int, T_Leaf_Temp
   !! \param BETA
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   subroutine GET_GAMMA_T_LI(T, BETA, T_Leaf_Int, T_Leaf_Temp, GAMMA_T_LI)
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)  :: T             !<
      real(fp), intent(in)  :: BETA          !< Temperature factor per species
      real(fp), intent(in)  :: T_Leaf_Int    !<
      real(fp), intent(in)  :: T_Leaf_Temp   !< Soil Moisture Attenuation Factor
      real(fp), intent(out) :: GAMMA_T_LI  !< Factor for the light-independent emissions

      ! Local Variables
      !----------------
      real(fp) :: L_T     !<
      !real(fp) :: L_PT_T  !<
      real(fp), parameter :: T_STANDARD = 303.0_fp

      !--------------------------------------------
      ! GET_GAMMAT_T_LI begins here!
      !--------------------------------------------
      L_T = T * T_Leaf_Temp + T_Leaf_Int
      GAMMA_T_LI = EXP( BETA * ( T - T_STANDARD ) )

      return
   end subroutine GET_GAMMA_T_LI

   !>
   !! \brief Computes the temperature sensitivity for the light-dependent fraction of emissions.
   !!
   !!  References:
   !!  (1 ) Guenther et al, 1995
   !!  (2 ) Guenther et al, 2006
   !!  (3 ) Guenther et al, MEGAN v2.1 user manual 2007-08
   ! ! (4 ) Guenther et al., GMD 2012 and MEGANv2.1 source code.
   !!
   !! \param T
   !! \param PT_15, PT_1
   !! \param CT1, CEO
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   !subroutine GET_GAMMA_T_LD(T, PT_15, PT_1, CT1, CEO, GAMMA_T_LD)
   subroutine GET_GAMMA_T_LD(T, PT_15, CT1, CEO, GAMMA_T_LD)
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)  :: T             !< Current leaf temperature [K]
      real(fp), intent(in)  :: PT_15         !< Average leaf temperature over the past 15 days
      !real(fp), intent(in)  :: PT_1          !< Average leaf temperature over the past arbitrary day(s). Not used at present
      real(fp), intent(in)  :: CT1, CEO      !< Compound-specific parameters for light-dependent temperature activity
      real(fp), intent(out) :: GAMMA_T_LD  !< Temperature activity factor for the light-dependent emissions

      ! Local Variables
      real(fp) :: C_T, CT2        !<
      real(fp) :: E_OPT, T_OPT, X !<
      ! Ideal gas constant [J/mol/K] (!!!! Note: the constant module is 8.314)
      real(fp), parameter :: R   = 8.3144598e-3_fp

      !--------------------------------------------
      ! GET_GAMMA_T_LD begins here!
      !--------------------------------------------
      E_OPT = CEO * EXP( 0.08_fp * ( PT_15  - 2.97e2_fp ) )
      T_OPT = 3.13e2_fp + ( 6.0e-1_fp * ( PT_15 - 2.97e2_fp ) )
      CT2   = 200.0_fp

      ! Variable related to temperature
      X     = ( 1.0_fp/T_OPT - 1.0_fp/T ) / R

      ! C_T: Effect of temperature on leaf BVOC emission, including
      ! effect of average temperature over previous 15 days, based on
      ! Eq 5a, 5b, 5c from Guenther et al, 1999.
      C_T   = E_OPT * CT2 * EXP( CT1 * X ) / &
         ( CT2 - CT1 * ( 1.0_fp - EXP( CT2 * X ) ) )

      ! Hourly emission activity = C_T
      ! Prevent negative values
      GAMMA_T_LD = MAX( C_T, 0.0_fp )

      return
   end subroutine GET_GAMMA_T_LD

   !>
   !! \brief Computes the gamma exchange activity factor which is sensitive to leaf area.
   !!
   !!  References:
   !!  (1 ) Guenther et al, 2006
   !!  (2 ) Guenther et al, MEGAN user manual 2007-08
   !!  (3 ) Guenther et al., GMD 2012 and MEGANv2.1 source code.
   !!
   !! \param CMLAI
   !! \param BIDIREXCH
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   subroutine GET_GAMMA_LAI(CMLAI, BIDIREXCH, GAMMA_LAI)
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)  :: CMLAI       !< Current month's LAI [cm2/cm2]
      logical,  intent(in)  :: BIDIREXCH   !< Logical flag for bidirectional exchange
      real(fp), intent(out) :: GAMMA_LAI   !< LAI factor

      ! Local Variables

      !--------------------------------------------
      ! GET_GAMMA_LAI begins here!
      !--------------------------------------------

      ! Formulation for birectional compounds is as described for
      ! ALD2 in Millet et al., ACP 2010
      IF ( BIDIREXCH ) THEN

         IF ( CMLAI <= 6.0_fp) THEN

            ! if lai less than 2:
            IF ( CMLAI <= 2.0_fp ) THEN
               GAMMA_LAI = 0.5_fp * CMLAI

               ! if between 2 and 6:
            ELSE
               GAMMA_LAI = 1.0_fp - 0.0625_fp * ( CMLAI - 2.0_fp )
            END IF

         ELSE
            ! keep at 0.75 for LAI > 6
            GAMMA_LAI = 0.75_fp
         END IF

         ! For all other compounds use the standard gamma_lai formulation
      ELSE
         !GAMMA_LAI = 0.49_fp * CMLAI / SQRT( 1.0_fp + 0.2_fp * CMLAI*CMLAI)
         GAMMA_LAI = 1.0_fp   !canopy add
      ENDIF

      return

   end subroutine GET_GAMMA_LAI

   !>
   !! \brief Computes the temperature sensitivity for the light-dependent
   !! fraction of emissions using the updated Canopy Model (Sam Silva's paper).
   !!
   !! References:
   !!
   !!
   !! \param T
   !! \param PT_15, PT_24
   !! \param CT1, CEO
   !! \param T_Leaf_Int, T_Leaf_Temp
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   !subroutine GET_GAMMA_T_LD_C(T, PT_15, PT_24, CT1, CEO, T_Leaf_Int, T_Leaf_Temp, GAMMA_T_LD_C )
   subroutine GET_GAMMA_T_LD_C(T, PT_24, CT1, CEO, T_Leaf_Int, T_Leaf_Temp, GAMMA_T_LD_C )

      IMPLICIT NONE

      ! Input Parameters
      !-----------------
      real(fp), intent(in) :: T             !<
      real(fp), intent(in) :: T_leaf_Int    !<
      real(fp), intent(in) :: T_Leaf_Temp   !<
      !!!!!!!TODO: Sam's version has no 15-day average temperature effects
      !real(fp), intent(in) :: PT_15         !< Average leaf temperature over the past 15 days. Not used at present
      real(fp), intent(in) :: PT_24         !< Average leaf temperature over the past day
      real(fp), intent(in) :: CT1, CEO      !< Compound-specific parameters

      ! Output Parameters
      !------------------
      real(fp), intent(out) :: GAMMA_T_LD_C !< threshold friction velocity


      ! Local Variables
      !-----------------
      real(fp) :: C_T, CT2         !<
      real(fp) :: E_OPT, T_OPT, X  !<
      real(fp) :: L_T, L_PT_T      !<
      ! Ideal gas constant [J/mol/K] (!!!! Note: the constant module is 8.314)
      real(fp), parameter :: R   = 8.3144598e-3_fp

      !=================================================================
      ! GET_GAMMA_T_LD begins here!
      !=================================================================

      L_T = T * T_Leaf_Temp + T_Leaf_Int
      L_PT_T = PT_24 * T_Leaf_Temp + T_Leaf_Int

      E_OPT = CEO * EXP( 0.1_fp * ( L_PT_T  - 2.97e2_fp ) )
      T_OPT = 3.125e2_fp + ( 6.0e-1_fp * ( L_PT_T - 2.97e2_fp ) )
      CT2   = 230.0_fp

      ! Variable related to temperature
      X     = ( 1.0_fp/T_OPT - 1.0_fp/L_T ) / R

      ! C_T: Effect of temperature on leaf BVOC emission, including
      ! effect of average temperature over previous 15 days, based on
      ! Eq 5a, 5b, 5c from Guenther et al, 1999.
      C_T   = E_OPT * CT2 * EXP( CT1 * X ) /       &
         ( CT2 - CT1 * ( 1.0_fp - EXP( CT2 * X ) ) )

      ! Hourly emission activity = C_T
      ! Prevent negative values
      IF (T < 260) THEN
         GAMMA_T_LD_C = 0.0_fp
      ELSE
         GAMMA_T_LD_C = MAX( C_T, 0.0_fp )
      ENDIF

      return

   end subroutine GET_GAMMA_T_LD_C


   !>
   !! \brief Computes the PCEEA gamma activity factor with sensitivity to light.
   !!
   !!  References:
   !!  (1 ) Guenther et al, 2006
   !!  (2 ) Guenther et al, 2007,   MEGAN v2.1 user guide
   !!
   !! \param Q_DIR_2,Q_DIFF_2
   !! \param PARDR_AVG_SIM, PARDF_AVG_SIM
   !! \param P_Leaf_LAI, P_Leaf_Int, LAI, PSTD
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   subroutine GET_GAMMA_PAR_C(Q_DIR_2, Q_DIFF_2, PARDR_AVG_SIM, PARDF_AVG_SIM, &
      P_Leaf_LAI, P_Leaf_Int, LAI, PSTD, GAMMA_P_C)
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)  :: Q_DIR_2          !< Direct PAR [W/m2]
      real(fp), intent(in)  :: Q_DIFF_2         !< Diffuse PAR [W/m2]
      real(fp), intent(in)  :: PARDR_AVG_SIM    !< Avg direct PAR [W/m2]
      real(fp), intent(in)  :: PARDF_AVG_SIM    !< Avg diffuse PAR [W/m2]
      real(fp), intent(in)  :: P_Leaf_LAI       !<
      real(fp), intent(in)  :: P_Leaf_Int       !<
      real(fp), intent(in)  :: LAI              !<
      real(fp), intent(in)  :: PSTD             !<
      real(fp), intent(out) :: GAMMA_P_C        !< GAMMA factor for light

      ! Local Variables
      real(fp) :: mmPARDR_DAILY    !<
      real(fp) :: mmPARDF_DAILY    !<
      real(fp) :: PAC_DAILY        !<
      real(fp) :: PAC_INSTANT      !<
      real(fp) :: C_PPFD           !<
      real(fp) :: PTOA             !<
      !real(fp) :: PHI              !<
      !real(fp) :: BETA,   SINbeta  !<
      real(fp) :: C1               !<
      real(fp) :: Alpha            !<

      !--------------------------------------------
      ! GET_GAMMA_PAR_C begins here!
      !--------------------------------------------
      ! Initialize
      C_PPFD   = 0.0_fp
      PTOA     = 0.0_fp

      ! Convert past light conditions to micromol/m2/s
      mmPARDR_DAILY   = PARDR_AVG_SIM
      mmPARDF_DAILY   = PARDF_AVG_SIM

      ! Work out the light at the top of the canopy.
      PAC_DAILY    = mmPARDR_DAILY + mmPARDF_DAILY
      PAC_INSTANT  = Q_DIR_2       +  Q_DIFF_2

      PAC_DAILY = PAC_DAILY * exp(P_Leaf_Int + P_Leaf_LAI * LAI)
      PAC_INSTANT = PAC_INSTANT * exp(P_Leaf_Int + P_Leaf_LAI * LAI)

      IF ( PAC_DAILY < 0.01_fp ) THEN

         GAMMA_P_C = 0.0_fp

      ELSE
         Alpha  = 0.004
         Alpha  = 0.004 - 0.0005*LOG(PAC_DAILY)
         C1 = 1.03
         C1 = 0.0468 * EXP(0.0005 * (PAC_DAILY - PSTD)) *     &
            (PAC_DAILY **  0.6)
         GAMMA_P_C = (Alpha * C1 * PAC_INSTANT) /             &
            ((1 + Alpha**2. * PAC_INSTANT**2.)**0.5)
      ENDIF
      ! Prevent negative values
      GAMMA_P_C = MAX( GAMMA_P_C , 0.0_fp )

      return
   end subroutine GET_GAMMA_PAR_C


   !>
   !! \brief Computes the
   !!
   !!References:
   !! (1 ) Probably from Sam Silva paper !!!Note
   !!
   !! \param CMLAI
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   subroutine GET_CDEA(CMLAI, CDEA)
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)  :: CMLAI         !< Current month's LAI [cm2/cm2]
      real(fp), intent(out) :: CDEA(5)       !<

      ! Local Variables
      !----------------
      real(fp) :: LAIdepth    !<
      REAL(fp) :: Cdepth(5)
      integer  :: K
      real(fp), parameter :: CCD1 = -0.2_fp
      real(fp), parameter :: CCD2 = 1.3_fp

      !--------------------------------------------
      ! GET_CDEA begins here!
      !--------------------------------------------
      Cdepth (1)   = 0.0469101
      Cdepth (2)   = 0.2307534
      Cdepth (3)   = 0.5
      Cdepth (4)   = 0.7692465
      Cdepth (5)   = 0.9530899
      DO K = 1, 5
         LAIdepth = CMLAI * Cdepth(K)
         IF ( LAIdepth .GT. 3 ) THEN
            LAIdepth = 3.0
         ENDIF
         CDEA(K) = CCD1 * LAIdepth + CCD2
      ENDDO

      return
   end subroutine GET_CDEA

   !>
   !! \brief Computes the gamma exchange activity factor which is sensitive to leaf age
   !!
   !!References:
   !! (1 ) Guenther et al, 2006
   !! (2 ) Guenther et al, MEGAN user manual 2007-08
   !! (3 ) Guenther et al., GMD 2012 and MEGANv2.1 source code
   !!
   !! \param CMLAI, PMLAI
   !! \param DBTWN, TT
   !! \param AN, AG, AM, AO
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   subroutine GET_GAMMA_AGE(CMLAI, PMLAI, DBTWN, TT, AN, AG, AM, AO, GAMMA_AGE)
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)  :: CMLAI      !< Current month's LAI [cm2/cm2]
      real(fp), intent(in)  :: PMLAI      !< Previous months LAI [cm2/cm2]
      real(fp), intent(in)  :: DBTWN      !< Number of days between
      real(fp), intent(in)  :: TT         !< Daily average temperature [K]
      real(fp), intent(in)  :: AN         !< Relative emission factor (new leaves)
      real(fp), intent(in)  :: AG         !< Relative emission factor (growing leaves)
      real(fp), intent(in)  :: AM         !< Relative emission factor (mature leaves)
      real(fp), intent(in)  :: AO         !< Relative emission factor (old leaves)
      real(fp), intent(out) :: GAMMA_AGE  !< leaf age activity factor

      ! Local Variables
      !----------------
      real(fp) :: FNEW  !< Fraction of new leaves in canopy
      real(fp) :: FGRO  !< Fraction of growing leaves
      real(fp) :: FMAT  !< Fraction of mature leaves
      real(fp) :: FOLD  !< Fraction of old leaves
      real(fp) :: TI    !< number of days after budbreak required to induce emissions
      real(fp) :: TM    !< number of days after budbreak required to reach peak emissions

      !--------------------------------------------
      ! GET_GAMMAT_AGE begins here!
      !--------------------------------------------

      !-----------------------
      ! Compute TI and TM
      ! (mpb,2009)
      !-----------------------
      TI = 0.0_fp !avoid uninitialized warning when compiling
      IF ( TT <= 303.0_fp ) THEN
         TI = 5.0_fp + 0.7_fp * ( 300.0_fp - TT )
      ELSEIF ( TT >  303.0_fp ) THEN
         TI = 2.9_fp
      ENDIF
      TM = 2.3_fp * TI

      !-----------------------
      ! Compute GAMMA_AGE
      !-----------------------

      !use rae function from pecision_mod to avoid "equality comparison for real" warning
      IF ( rae(CMLAI, PMLAI) ) THEN !(i.e. LAI stays the same)

         FNEW = 0.0_fp
         FGRO = 0.1_fp
         FMAT = 0.8_fp
         FOLD = 0.1_fp

      ELSE IF ( CMLAI > PMLAI ) THEN !(i.e. LAI has INcreased)

         ! Calculate Fnew
         IF ( DBTWN > TI ) THEN
            FNEW = ( TI / DBTWN ) * ( 1.0_fp -  PMLAI / CMLAI )
         ELSE
            FNEW = 1.0_fp - ( PMLAI / CMLAI )
         ENDIF

         ! Calculate FMAT
         IF ( DBTWN > TM ) THEN
            FMAT = ( PMLAI / CMLAI ) + &
               (( DBTWN - TM ) / DBTWN )*( 1.0_fp -  PMLAI / CMLAI )
         ELSE
            FMAT = ( PMLAI / CMLAI )
         ENDIF

         ! Calculate Fgro and Fold
         FGRO = 1.0_fp - FNEW - FMAT
         FOLD = 0.0_fp

      ELSE ! This is the case if  PMLAI > CMLAI (i.e. LAI has DEcreased)

         FNEW = 0.0_fp
         FGRO = 0.0_fp
         FOLD = ( PMLAI - CMLAI ) / PMLAI
         FMAT = 1.0_fp - FOLD

      ENDIF

      ! Age factor
      GAMMA_AGE = FNEW * AN + FGRO * AG + FMAT * AM + FOLD * AO

      ! Prevent negative values
      GAMMA_AGE = MAX( GAMMA_AGE , 0.0_fp )

      return
   end subroutine GET_GAMMA_AGE


   !>
   !! \brief Computes ctivity factor for soil moisture
   !!
   !!References:
   !! (1 ) Guenther et al, 2006
   !! (2 ) Guenther et al., GMD 2012 and MEGANv2.1 source code
   !!
   !! \param GWETROOT
   !! \param CMPD
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   subroutine GET_GAMMA_SM(GWETROOT, CMPD, GAMMA_SM)
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)            :: GWETROOT   !< Relative root zone wetness (unitless)
      character(len=256), intent(in)  :: CMPD       !< Compound name
      real(fp), intent(out)           :: GAMMA_SM   !< Activity factor

      ! Local Variables
      !----------------
      real(fp) :: GWETROOT2    !<

      !--------------------------------------------
      ! GET_GAMMAT_SM begins here!
      !--------------------------------------------
      ! By default gamma_sm is 1.0
      GAMMA_SM = 1.0_fp

      ! Error trap: GWETROOT must be between 0.0 and 1.0 (ckeller, 4/16/15)
      !GWETROOT = MIN(MAX(ExtState%GWETROOT%Arr%Val(I,J),0.0_fp),1.0_fp)
      GWETROOT2 = MIN(MAX(GWETROOT,0.0_fp),1.0_fp)

      IF ( TRIM( CMPD ) == 'ALD2' .OR. TRIM ( CMPD ) == 'EOH' ) THEN

         ! GWETROOT = degree of saturation or wetness in the root-zone
         ! (top meter of soil). This is defined as the ratio of the volumetric
         ! soil moisture to the porosity. We use a soil moisture activity factor
         ! for ALD2 to account for stimulation of emission by flooding.
         ! (Millet et al., ACP 2010)
         ! Constant value of 1.0 for GWETROOT = 0-0.9, increasing linearly to
         ! 3.0 at GWETROOT =1.
         GAMMA_SM = MAX( 20.0_fp * GWETROOT - 17.0_fp, 1.0_fp)

      ENDIF

      return
   end subroutine GET_GAMMA_SM

   !>
   !! \brief Computes the CO2 activity factor associated with CO2 inhibition of isoprene emission.
   !!
   !!References:
   !! (1 ) Heald, C. L., Wilkinson, M. J., Monson, R. K., Also, C. A.,
   !!      Wang, G. L., and Guenther, A.: Response of isoprene emission
   !!      to ambient co(2) changes and implications for global budgets,
   !!      Global Change Biology, 15, 1127-1140, 2009.
   !! (2 ) Wilkinson, M. J., Monson, R. K., Trahan, N., Lee, S., Brown, E.,
   !!      Jackson, R. B., Polley, H. W., Fay, P. A., and Fall, R.: Leaf
   !!      isoprene emission rate as a function of atmospheric CO2
   !!      concentration, Global Change Biology, 15, 1189-1200, 2009.
   !!(3 )  Possell, M., and Hewitt, C. N.: Isoprene emissions from plants
   !!      are mediated by atmospheric co2 concentrations, Global Change
   !!      Biology, 17, 1595-1610, 2011.
   !!
   !! \param CO2a
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   subroutine GET_GAMMA_CO2(CO2a, GAMMA_CO2)
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)  :: CO2a         !< Atmospheric CO2 conc [ppmv]
      real(fp), intent(out)  :: GAMMA_CO2   !< CO2 activity factor [unitless]

      ! Local Variables
      !----------------
      real(fp) :: CO2i        !< Intercellular CO2 conc [ppmv]
      real(fp) :: ISMAXi      !< Asymptote for intercellular CO2
      real(fp) :: HEXPi       !< Exponent for intercellular CO2
      real(fp) :: CSTARi      !< Scaling coef for intercellular CO2
      real(fp) :: ISMAXa      !< Asymptote for atmospheric CO2
      real(fp) :: HEXPa       !< Exponent for atmospheric CO2
      real(fp) :: CSTARa      !< Scaling coef for atmospheric CO2
      logical  :: LPOSSELL    !< Use Possell & Hewitt (2011)?
      logical  :: LWILKINSON  !< Use Wilkinson et al. (2009)?

      !--------------------------------------------
      ! GET_GAMMAT_CO2 begins here!
      !--------------------------------------------

      !----------------------------------------------------------
      ! Choose between two alternative CO2 inhibition schemes
      !----------------------------------------------------------

      ! Empirical relationship of Possell & Hewitt (2011) based on nine
      ! experimental studies including Wilkinson et al. (2009). This is
      ! especially recommended for sub-ambient CO2 concentrations:
      LPOSSELL    = .TRUE.   ! Default option

      ! Semi-process-based parameterization of Wilkinson et al. (2009),
      ! taking into account of sensitivity to intercellular CO2
      ! fluctuation, which is here set as a constant fraction of
      ! atmospheric CO2:
      LWILKINSON  = .FALSE.   ! Set .TRUE. only if LPOSSELL = .FALSE.

      !-----------------------
      ! Compute GAMMA_CO2
      !-----------------------

      IF ( LPOSSELL ) THEN

         ! Use empirical relationship of Possell & Hewitt (2011):
         GAMMA_CO2 = 8.9406_fp / ( 1.0_fp + 8.9406_fp * 0.0024_fp * CO2a )

      ELSEIF ( LWILKINSON ) THEN

         ! Use parameterization of Wilkinson et al. (2009):

         ! Parameters for intercellular CO2 using linear interpolation:
         IF ( CO2a <= 600.0_fp ) THEN
            ISMAXi = 1.036_fp  - (1.036_fp - 1.072_fp) / &
               (600.0_fp - 400.0_fp) * (600.0_fp - CO2a)
            HEXPi  = 2.0125_fp - (2.0125_fp - 1.7000_fp) / &
               (600.0_fp - 400.0_fp) * (600.0_fp - CO2a)
            CSTARi = 1150.0_fp - (1150.0_fp - 1218.0_fp) / &
               (600.0_fp - 400.0_fp) * (600.0_fp - CO2a)
         ELSEIF ( CO2a > 600.0_fp .AND. CO2a < 800.0_fp ) THEN
            ISMAXi = 1.046_fp  - (1.046_fp - 1.036_fp) / &
               (800.0_fp - 600.0_fp) * (800.0_fp - CO2a)
            HEXPi  = 1.5380_fp - (1.5380_fp - 2.0125_fp) / &
               (800.0_fp - 600.0_fp) * (800.0_fp - CO2a)
            CSTARi = 2025.0_fp - (2025.0_fp - 1150.0_fp) / &
               (800.0_fp - 600.0_fp) * (800.0_fp - CO2a)
         ELSE
            ISMAXi = 1.014_fp - (1.014_fp - 1.046_fp) / &
               (1200.0_fp - 800.0_fp) * (1200.0_fp - CO2a)
            HEXPi  = 2.8610_fp - (2.8610_fp - 1.5380_fp) / &
               (1200.0_fp - 800.0_fp) * (1200.0_fp - CO2a)
            CSTARi = 1525.0_fp - (1525.0_fp - 2025.0_fp) / &
               (1200.0_fp - 800.0_fp) * (1200.0_fp - CO2a)
         ENDIF

         ! Parameters for atmospheric CO2:
         ISMAXa    = 1.344_fp
         HEXPa     = 1.4614_fp
         CSTARa    = 585.0_fp

         ! For now, set CO2_Ci = 0.7d0 * CO2_Ca as recommended by Heald
         ! et al. (2009):
         CO2i      = 0.7_fp * CO2a

         ! Compute GAMMA_CO2:
         GAMMA_CO2 = ( ISMAXi -  ISMAXi * CO2i**HEXPi / &
            ( CSTARi**HEXPi + CO2i**HEXPi ) )  &
            * ( ISMAXa - ISMAXa * ( 0.7_fp * CO2a )**HEXPa / &
            ( CSTARa**HEXPa + ( 0.7_fp * CO2a )**HEXPa ) )

      ELSE

         ! No CO2 inhibition scheme is used; GAMMA_CO2 set to unity:
         GAMMA_CO2 = 1.0_fp

      ENDIF

      return
   end subroutine GET_GAMMA_CO2

   !>
   !! \brief Computes the normalization factor needed to compute emissions
   !!
   !!References:
   !!(1 ) Guenther et al., GMD 2012 and MEGANv2.1 source code
   !!(2 ) Created by dbm 11/2012. We calculate only 1 normalization factor for all
   !!     compounds based on the isoprene gamma values. Formally there should be a
   !!     different normalization factor for each compound, but we are following
   !!     Alex Guenther's approach here and the MEGAN source code.
   !!
   !! \param D2RAD_FAC
   !! \param NORM_FAC
   !! \param RC
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   subroutine CALC_NORM_FAC(D2RAD_FAC, NORM_FAC, RC)
      use Error_Mod,     Only : CC_SUCCESS !, CC_FAILURE, CC_Error
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)    :: D2RAD_FAC   !<
      real(fp), intent(out)    :: NORM_FAC    !<
      integer,  intent(out) :: RC          !<

      ! Local Variables
      !----------------
      real(fp) :: PAC_DAILY
      real(fp) :: GAMMA_T_LI_STANDARD
      real(fp) :: GAMMA_SM_STANDARD
      real(fp) :: CMLAI, GAMMA_LAI_STANDARD
      real(fp) :: GAMMA_AGE_STANDARD
      real(fp) :: PT_15, T, R, CEO, CT1, E_OPT, T_OPT, CT2, X
      real(fp) :: LDF, GAMMA_STANDARD
      !canopy add by Sam Silva (some variables are unused because of this)
      !real(fp) :: PHI, AAA, BBB, EA2L, GAMMA_P_STANDARD, GAMMA_T_LD_STANDARD
      real(fp) ::  SunF, GAMMA_TP_STANDARD
      real(fp) :: T_Leaf_Int_Sun(5)
      real(fp) :: T_Leaf_Int_Shade(5)
      real(fp) :: T_Leaf_Temp_Sun(5)
      real(fp) :: T_Leaf_Temp_Shade(5)
      !real(fp) :: T_Leaf_Wind_Sun(5)
      !real(fp) :: T_Leaf_Wind_Shade(5)
      real(fp) :: P_Leaf_Int_Sun(5)
      real(fp) :: P_Leaf_Int_Shade(5)
      real(fp) :: P_Leaf_LAI_Sun(5)
      real(fp) :: P_Leaf_LAI_Shade(5)
      real(fp) :: Distgauss(5), CDEA(5), VPGWT(5)
      real(fp) :: EA1L, GAMMA_T_LD_SUN, GAMMA_T_LD_SHADE
      real(fp) :: L_PT_T, L_T, C1, LAI, GAMMA_PAR_SUN
      real(fp) :: GAMMA_PAR_SHADE, ALPHA, PAC_INSTANT
      integer  :: Q

      !--------------------------------------------
      ! CALC_NORM_FAC begins here!
      !--------------------------------------------

      ! -----------------
      ! GAMMA_P for standard conditions
      ! -----------------
      ! Based on Eq. 11b from Guenther et al., 2006
      ! Using standard conditions of phi = 0.6, solar angle of 60 deg,
      ! and P_daily = 400
      ! Note corrigendum for Eq. 11b in that paper, should be a
      ! minus sign before the 0.9. canopy add
      !PAC_DAILY = 400.0_fp
      !PHI       = 0.6_fp
      !BBB       = 1.0_fp + 0.0005_fp *( PAC_DAILY - 400.0_fp )
      !AAA       = ( 2.46_fp * BBB * PHI ) - ( 0.9_fp * PHI**2 )
      ! sin(60) = 0.866
      !GAMMA_P_STANDARD = 0.866_fp * AAA

      ! -----------------
      ! GAMMA_T_LI for standard conditions
      ! -----------------
      ! gamma_t_li = EXP( Beta * ( T - T_Standard ) )
      ! This is 1.0 for T = T_Standard
      GAMMA_T_LI_STANDARD = 1.0_fp

      ! -----------------
      ! GAMMA_SM for standard conditions
      ! -----------------
      ! Standard condition is soil moisture = 0.3 m^3/m^3
      ! GAMMA_SM = 1.0 for all compounds under this condition
      GAMMA_SM_STANDARD = 1.0_fp
      ! -----------------
      ! GAMMA_TP for standard conditions canopy add
      ! -----------------
      ! gamma_t_li = EXP( Beta * ( T - T_Standard ) )
      ! This is 1.0 for T = T_Standard

      CALL GET_CDEA( 5.0_fp, CDEA )
      Distgauss = (/0.0469101, 0.2307534, 0.5, 0.7692465,               &
         0.9530899/)
      VPGWT = (/0.1184635, 0.2393144, 0.284444444,                      &
         0.2393144, 0.1184635/)
      P_Leaf_Int_Sun  = (/1.0831_fp, 1.0964_fp, 1.1036_fp,              &
         1.0985_fp, 1.0901_fp/)
      P_Leaf_Int_Shade = (/0.8706_fp, 0.8895_fp, 0.9160_fp,             &
         0.9407_fp, 0.9564_fp/)

      P_Leaf_LAI_Sun = (/0.0018_fp, -0.1281_fp, -0.2977_fp,             &
         -0.4448_fp, -0.5352_fp/)
      P_Leaf_LAI_Shade = (/0.0148_fp, -0.1414_fp, -0.3681_fp,           &
         -0.5918_fp, -0.7425_fp/)

      T_Leaf_Int_Sun  = (/-13.891_fp, -12.322_fp, -1.032_fp,            &
         -5.172_fp, -5.589_fp/)
      T_Leaf_Int_Shade = (/-12.846_fp, -11.343_fp, -1.068_fp,           &
         -5.551_fp, -5.955_fp/)

      T_Leaf_Temp_Sun = (/1.064_fp, 1.057_fp, 1.031_fp,                 &
         1.050_fp, 1.051_fp/)
      T_Leaf_Temp_Shade = (/1.060_fp, 1.053_fp, 1.031_fp,               &
         1.051_fp, 1.052_fp/)

      GAMMA_TP_STANDARD = 0.0_fp
      LDF = 1.0_fp
      LAI = 5.0_fp
      DO Q = 1, 5

         PAC_INSTANT  = 1500.0_fp/4.766_fp
         PAC_DAILY = 740.0_fp/4.766_fp

         PAC_INSTANT = PAC_INSTANT * exp(P_Leaf_Int_Sun(Q) +            &
            P_Leaf_LAI_Sun(Q) * LAI)
         PAC_DAILY = PAC_DAILY * exp(P_Leaf_Int_Sun(Q) +                &
            P_Leaf_LAI_Sun(Q) * LAI)
         Alpha  = 0.004 - 0.0005*LOG(PAC_DAILY)
         C1 = 0.0468 * EXP(0.0005 * (PAC_DAILY - 200.0_fp)) *           &
            (PAC_DAILY **  0.6)
         GAMMA_PAR_Sun = (Alpha * C1 * PAC_INSTANT) /                   &
            ((1 + Alpha**2. * PAC_INSTANT**2.)**0.5)

         PAC_INSTANT  = 1500.0_fp/4.766_fp
         PAC_DAILY = 740.0_fp/4.766_fp
         PAC_DAILY = PAC_DAILY * exp(P_Leaf_Int_Shade(Q) +              &
            P_Leaf_LAI_Shade(Q) * LAI)
         PAC_INSTANT = PAC_INSTANT * exp(P_Leaf_Int_Shade(Q) +          &
            P_Leaf_LAI_Shade(Q) * LAI)
         Alpha  = 0.004 - 0.0005*LOG(PAC_DAILY)
         C1 = 0.0468 * EXP(0.0005 * (PAC_DAILY - 50.0_fp)) *            &
            (PAC_DAILY **  0.6)
         GAMMA_PAR_Shade = (Alpha * C1 * PAC_INSTANT) /                 &
            ((1 + Alpha**2. * PAC_INSTANT**2.)**0.5)

         PT_15 = 298.5_fp
         T     = 303.0_fp
         R     = 8.3144598e-3_fp
         CEO = 2.0_fp
         CT1 = 95.0_fp
         CT2   = 230.0_fp

         L_T = T * T_Leaf_Temp_Sun(Q) + T_Leaf_Int_Sun(Q)
         L_PT_T = PT_15 * T_Leaf_Temp_Sun(Q) + T_Leaf_Int_Sun(Q)
         E_OPT = CEO * EXP( 0.1_fp * ( L_PT_T  - 2.97e2_fp ) )
         T_OPT = 3.125e2_fp + ( 6.0e-1_fp * ( L_PT_T - 2.97e2_fp ) )
         X     = ( 1.0_fp/T_OPT - 1.0_fp/L_T ) / R
         GAMMA_T_LD_Sun   = E_OPT * CT2 * EXP( CT1 * X ) /              &
            ( CT2 - CT1 * ( 1.0_fp - EXP( CT2 * X ) ) )

         L_T = T * T_Leaf_Temp_Shade(Q) + T_Leaf_Int_Shade(Q)
         L_PT_T = PT_15 * T_Leaf_Temp_Shade(Q) + T_Leaf_Int_Shade(Q)
         E_OPT = CEO * EXP( 0.08_fp * ( L_PT_T  - 2.97e2_fp ) )
         T_OPT = 3.125e2_fp + ( 6.0e-1_fp * ( L_PT_T - 2.97e2_fp ) )
         X     = ( 1.0_fp/T_OPT - 1.0_fp/L_T ) / R
         GAMMA_T_LD_Shade   = E_OPT * CT2 * EXP( CT1 * X ) /            &
            ( CT2 - CT1 * ( 1.0_fp - EXP( CT2 * X ) ) )

         call Calc_Sun_Frac(5.0_fp,SIN(60.0_fp*D2RAD_FAC),            &
            Distgauss(Q), SunF)
         Ea1L  =  CDEA(Q) * GAMMA_PAR_Sun * GAMMA_T_LD_Sun * SunF +     &
            GAMMA_PAR_Shade * GAMMA_T_LD_Shade * (1-SunF)

         !   write(*,*) ' '
         !   write(*,*) '--- GET_MEGAN_NormFrac --- '
         !   write(*,*) 'GAMMA_TP_STANDARD      : ', GAMMA_TP_STANDARD
         !   write(*,*) 'Ea1L      : ', Ea1L
         !   write(*,*) 'SunF      : ', SunF
         !   write(*,*) 'CDEA      : ', CDEA(Q)
         !   write(*,*) 'VPGWT      : ', VPGWT(Q)
         !   write(*,*) 'Q      : ', Q
         !   write(*,*) 'Q      : ', Q
         !   write(*,*) 'GAMMA_PAR_Sun      : ', GAMMA_PAR_Sun
         !   write(*,*) 'GAMMA_PAR_Shade      : ', GAMMA_PAR_Shade
         !   write(*,*) 'GAMMA_T_LD_Sun      : ', GAMMA_T_LD_Sun
         !   write(*,*) 'GAMMA_T_LD_Shade      : ', GAMMA_T_LD_Shade

         GAMMA_TP_STANDARD  = GAMMA_TP_STANDARD +                       &
            (Ea1L*LDF)* VPGWT(Q)
      ENDDO

      !  write(*,*) ' '
      !   write(*,*) '--- GET_MEGAN_NormFrac --- '
      !   write(*,*) 'GAMMA_TP_STANDARD      : ', GAMMA_TP_STANDARD
      !   write(*,*) 'Ea1L      : ', Ea1L
      !   write(*,*) 'SunF      : ', SunF
      !   write(*,*) 'CDEA      : ', CDEA(Q)
      !   write(*,*) 'VPGWT      : ', VPGWT(Q)
      !   write(*,*) 'Q      : ', Q

      ! -----------------
      ! GAMMA_LAI for standard conditions
      ! -----------------
      ! Standard condition is LAI = 5
      CMLAI = 5.0_fp
      !canopy add
      !GAMMA_LAI_STANDARD = 0.49_fp * CMLAI / SQRT( 1.0_fp + 0.2_fp * CMLAI*CMLAI )
      GAMMA_LAI_STANDARD = CMLAI
      ! -----------------
      ! GAMMA_AGE for standard conditions
      ! -----------------
      ! Standard condition is 0% new, 10% growing, 80% mature, 10% old foliage
      ! Isoprene uses A_NEW = 0.05d0, A_GRO = 0.6d0, A_MAT = 1.d0, A_OLD = 0.9d0
      GAMMA_AGE_STANDARD = 0.1_fp*0.6_fp + 0.8_fp*1.0_fp + 0.1_fp*0.9_fp

      ! -----------------
      ! GAMMA_T_LD for standard conditions
      ! -----------------
      ! Standard condition is
      ! PT_15 = average leaf temp over past 24-240 hours = 297K
      ! T = air temperature = 303K
      !PT_15 = 297.0_fp
      !T     = 303.0_fp
      !R     = 8.3144598e-3_fp
      ! parameters for isoprene
      !CEO = 2.0_fp
      !CT1 = 95.0_fp

      !E_OPT = CEO * EXP( 0.08_fp * ( PT_15  - 2.97e2_fp ) )
      !T_OPT = 3.13e2_fp + ( 6.0e-1_fp * ( PT_15 - 2.97e2_fp ) )
      !CT2   = 200.0_fp

      ! Variable related to temperature
      !X     = ( 1.0_fp/T_OPT - 1.0_fp/T ) / R

      !GAMMA_T_LD_STANDARD = E_OPT * CT2 * EXP( CT1 * X ) / &
      !                      ( CT2 - CT1 * ( 1.0_fp - EXP( CT2 * X ) ) )

      ! -----------------
      ! Overall GAMMA_STANDARD
      ! -----------------
      ! LDF = 1d0 for isoprene
      !LDF = 1.0_fp
      !GAMMA_STANDARD = &
      !     GAMMA_AGE_STANDARD * GAMMA_SM_STANDARD * GAMMA_LAI_STANDARD &
      !     * ((1.0_fp - LDF) * GAMMA_T_LI_STANDARD &
      !     + (LDF * GAMMA_P_STANDARD * GAMMA_T_LD_STANDARD))
      ! This ends up being 1.0101081.
      !canopy add
      GAMMA_STANDARD = GAMMA_AGE_STANDARD * GAMMA_SM_STANDARD *   &
         GAMMA_TP_STANDARD * GAMMA_LAI_STANDARD

      NORM_FAC = 1.0_fp / GAMMA_STANDARD

      !write(*,*) ' '
      !write(*,*) '--- GET_MEGAN_NormFrac --- '
      !write(*,*) 'GAMMA_STANDARD      : ', GAMMA_STANDARD
      !write(*,*) 'GAMMA_AGE_STANDARD      : ', GAMMA_AGE_STANDARD
      !write(*,*) 'GAMMA_SM_STANDARD      : ', GAMMA_SM_STANDARD
      !write(*,*) 'GAMMA_TP_STANDARD      : ', GAMMA_TP_STANDARD
      !write(*,*) 'GAMMA_LAI_STANDARD      : ', GAMMA_LAI_STANDARD

      RC = CC_SUCCESS
      return
   end subroutine CALC_NORM_FAC

   !>
   !! \brief Computes Emission Factors for all biogenic VOC species
   !!  Note; I separate the reading AE into the main GET_MEGAN_EMIS function
   !!        Here only the species need to be calculaed are included
   !!References:
   !! (1 ) Guenther et al, 2004
   !!
   !! \param PFT_16
   !! \param CMPD
   !! \param AE, RC
   !!
   !! \ingroup catchem_bvoc_process
   !!!>
   subroutine CALC_AEF(PFT_16, CMPD ,AE, RC)
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)            :: PFT_16(16)  !< 16 PFT array (TODO:read from MetState??)
      character(len=256), intent(in)  :: CMPD        !< compound name
      real(fp), intent(out)           :: AE          !< annual emission factor for the compound
      integer,  intent(out)    :: RC                 !< Success or Failure

      ! Local Variables
      !----------------
      integer                 :: P, ARR_IND
      real(fp)                :: FACTOR
      character(len=255)       :: MSG, thisLoc
      REAL(fp)                :: PFT_EF_OMON(15), PFT_EF_MOH(15)
      REAL(fp)                :: PFT_EF_ACET(15), PFT_EF_BIDR(15)
      REAL(fp)                :: PFT_EF_STRS(15), PFT_EF_OTHR(15)
      ! --->
      ! dbm, compute EF maps for a-pinene and myrcene as well since there seems to
      ! be an issue with the EF maps for these species provided on the MEGAN
      ! data portal
      REAL(fp)                :: PFT_EF_APIN(15), PFT_EF_MYRC(15)
      ! <---
      REAL(fp)                :: PFT_EF_FARN(15), PFT_EF_BCAR(15)
      REAL(fp)                :: PFT_EF_OSQT(15)
      REAL(fp)                :: EM_FRAC_ALD2(15), EM_FRAC_EOH(15)
      REAL(fp)                :: EM_FRAC_FAXX(15), EM_FRAC_AAXX(15)
      REAL(fp)                :: EM_FRAC_CH2O(15)
      !try to calculate the seven read-in species online
      REAL(fp)                :: PFT_EF_ISOP(15), PFT_EF_MBOX(15)
      REAL(fp)                :: PFT_EF_BPIN(15), PFT_EF_CARE(15)
      REAL(fp)                :: PFT_EF_LIMO(15), PFT_EF_OCIM(15)
      REAL(fp)                :: PFT_EF_SABI(15)
      !-----------------------------------------------------------------
      ! Point to PFT fractions (the array needs to be in that order)
      !-----------------------------------------------------------------
      ! CLM4 PFT coverage (unitless)
      ! From Table 3 in Guenther et al., 2012
      ! PFT_BARE                : Bare
      ! PFT_NDLF_EVGN_TMPT_TREE : Needleleaf evergreen temperate tree
      ! PFT_NDLF_EVGN_BORL_TREE : Needleleaf evergreen boreal tree
      ! PFT_NDLF_DECD_BORL_TREE : Needleleaf deciduous boreal tree
      ! PFT_BDLF_EVGN_TROP_TREE : Broadleaf evergreen tropical tree
      ! PFT_BDLF_EVGN_TMPT_TREE : Broadleaf evergreen temperate tree
      ! PFT_BDLF_DECD_TROP_TREE : Broadleaf deciduous tropical tree
      ! PFT_BDLF_DECD_TMPT_TREE : Broadleaf deciduous temperate tree
      ! PFT_BDLF_DECD_BORL_TREE : Broadleaf deciduous boreal tree
      ! PFT_BDLF_EVGN_SHRB      : Broadleaf evergreen temperate shrub
      ! PFT_BDLF_DECD_TMPT_SHRB : Broadleaf deciduous temperate shrub
      ! PFT_BDLF_DECD_BORL_SHRB : Broadleaf deciduous boreal shrub
      ! PFT_C3_ARCT_GRSS        : Arctic C3 grass
      ! PFT_C3_NARC_GRSS        : Cool C3 grass
      ! PFT_C4_GRSS             : Warm C4 grass
      ! PFT_CROP                : Crop

      ! --------------------------------------------------------------------------------
      ! PFT-specific EFs from Table 2 in Guenther et al., 2012
      ! in ug compound/m2/h
      ! PFTs 1-15 in the table correspond to #2-16
      ! (i.e., excluding bare ground #1) in the above array.
      ! --------------------------------------------------------------------------------
      ! Compound Class EF1 EF2 EF3 EF4 EF5 EF6 EF7 EF8 EF9 EF10 EF11 EF12 EF13 EF14 EF15
      ! --------------------------------------------------------------------------------
      ! Other Monoterp 180 180 170 150 150 150 150 150 110 200  110  5    5    5    5
      ! Methanol       900 900 900 500 900 500 900 900 900 900  900  500  500  500  900
      ! Acetone        240 240 240 240 240 240 240 240 240 240  240  80   80   80   80
      ! Bidirect VOC   500 500 500 500 500 500 500 500 500 500  500  80   80   80   80
      ! Stress VOC     300 300 300 300 300 300 300 300 300 300  300  300  300  300  300
      ! Other VOC      140 140 140 140 140 140 140 140 140 140  140  140  140  140  140
      ! a-Pinene       500 500 510 600 400 600 400 400 200 300  200    2    2    2    2
      ! Myrcene         70  70  60  80  30  80  30  30  30  50   30  0.3  0.3  0.3  0.3
      ! a-Farnesene     40  40  40  60  40  60  40  40  40  40   40    3    3    3    4
      ! b-Carophyllene  80  80  80  60  40  60  40  40  50  50   50    1    1    1    4
      ! Other sesqt.   120 120 120 120 100 120 100 100 100 100  100    2    2    2    2
      ! --------------------------------------------------------------------------------

      ! One thing to note is these are net emissions to the canopy atmosphere
      ! but not net emissions to the above canopy atmosphere since they don't
      !  account for within-canopy deposition. Only an issue for OVOCs.

      !try to calculate the seven read-in species online now
      !               EF1        EF2        EF3         EF4        EF5
      PFT_EF_ISOP = (/600.0_fp,  3000.0_fp, 1.0_fp,     7000.0_fp, 10000.0_fp, &
      !               EF6        EF7        EF8         EF9        EF10
         7000.0_fp, 10000.0_fp,11000.0_fp, 2000.0_fp, 4000.0_fp,  &
      !               EF11       EF12       EF13        EF14       EF15
         4000.0_fp, 1600.0_fp, 800.0_fp,   200.0_fp,  1.0_fp/)

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_MBOX = (/700.0_fp, 60.0_fp,  0.01_fp,  0.01_fp,  0.01_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         0.01_fp,  0.01_fp,  2.0_fp,   0.01_fp,  0.01_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         0.01_fp,  0.01_fp,  0.01_fp,  0.01_fp,  0.01_fp/)

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_BPIN = (/300.0_fp, 300.0_fp, 200.0_fp, 120.0_fp, 130.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         120.0_fp, 130.0_fp, 130.0_fp, 100.0_fp, 150.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         100.0_fp, 1.5_fp  , 1.5_fp  , 1.5_fp  , 1.5_fp/)

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_CARE = (/160.0_fp, 160.0_fp, 80.0_fp, 40.0_fp,   30.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         40.0_fp,  30.0_fp,  30.0_fp,  30.0_fp,  100.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         30.0_fp,  0.3_fp  , 0.3_fp  , 0.3_fp  , 0.3_fp/)

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_LIMO = (/100.0_fp, 100.0_fp, 130.0_fp, 80.0_fp,  80.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         80.0_fp,  80.0_fp,  80.0_fp,  60.0_fp,  100.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         60.0_fp,  0.7_fp  , 0.7_fp  , 0.7_fp  , 0.7_fp/)

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_OCIM = (/70.0_fp,  70.0_fp,  60.0_fp,  150.0_fp, 120.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         150.0_fp, 120.0_fp, 120.0_fp, 90.0_fp,  150.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         90.0_fp,  2.0_fp  , 2.0_fp  , 2.0_fp  , 2.0_fp/)

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_SABI = (/70.0_fp,  70.0_fp,  40.0_fp,  80.0_fp,  50.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         80.0_fp,  50.0_fp,  50.0_fp,  50.0_fp,  70.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         50.0_fp,  0.7_fp  , 0.7_fp  , 0.7_fp  , 0.7_fp/)

      ! default EF online calculation below
      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_OMON = (/180.0_fp, 180.0_fp, 170.0_fp, 150.0_fp, 150.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         150.0_fp, 150.0_fp, 150.0_fp, 110.0_fp, 200.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         110.0_fp, 5.0_fp  , 5.0_fp  , 5.0_fp  , 5.0_fp/)

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_MOH  = (/900.0_fp, 900.0_fp, 900.0_fp, 500.0_fp, 900.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         500.0_fp, 900.0_fp, 900.0_fp, 900.0_fp, 900.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         900.0_fp, 500.0_fp, 500.0_fp, 500.0_fp, 900.0_fp/)

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_ACET = (/240.0_fp, 240.0_fp, 240.0_fp, 240.0_fp, 240.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         240.0_fp, 240.0_fp, 240.0_fp, 240.0_fp, 240.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         240.0_fp, 80.0_fp , 80.0_fp , 80.0_fp , 80.0_fp/)

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_BIDR = (/500.0_fp, 500.0_fp, 500.0_fp, 500.0_fp, 500.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         500.0_fp, 500.0_fp, 500.0_fp, 500.0_fp, 500.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         500.0_fp, 80.0_fp , 80.0_fp , 80.0_fp , 80.0_fp/)

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_STRS = (/300.0_fp, 300.0_fp, 300.0_fp, 300.0_fp, 300.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         300.0_fp, 300.0_fp, 300.0_fp, 300.0_fp, 300.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         300.0_fp, 300.0_fp, 300.0_fp, 300.0_fp, 300.0_fp/)

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_OTHR = (/140.0_fp, 140.0_fp, 140.0_fp, 140.0_fp, 140.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         140.0_fp, 140.0_fp, 140.0_fp, 140.0_fp, 140.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         140.0_fp, 140.0_fp, 140.0_fp, 140.0_fp, 140.0_fp/)

      ! ---> Now compute EFs for a-pinene and myrcene as well (dbm, 12/2012)
      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_APIN = (/500.0_fp, 500.0_fp, 510.0_fp, 600.0_fp, 400.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         600.0_fp, 400.0_fp, 400.0_fp, 200.0_fp, 300.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         200.0_fp, 2.0_fp,   2.0_fp,   2.0_fp,   2.0_fp/)

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_MYRC = (/70.0_fp,  70.0_fp,  60.0_fp,  80.0_fp,  30.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         80.0_fp,  30.0_fp,  30.0_fp,  30.0_fp,  50.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         30.0_fp,  0.3_fp,   0.3_fp,   0.3_fp,   0.3_fp/)
      ! <---

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_FARN = (/40.0_fp,  40.0_fp,  40.0_fp,  60.0_fp,  40.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         60.0_fp,  40.0_fp,  40.0_fp,  40.0_fp,  40.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         40.0_fp,  3.0_fp,   3.0_fp,   3.0_fp,   4.0_fp/)

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_BCAR = (/80.0_fp,  80.0_fp,  80.0_fp,  60.0_fp,  40.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         60.0_fp,  40.0_fp,  40.0_fp,  50.0_fp,  50.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         50.0_fp,  1.0_fp,   1.0_fp,   1.0_fp,   4.0_fp/)

      !               EF1       EF2       EF3       EF4       EF5
      PFT_EF_OSQT = (/120.0_fp, 120.0_fp, 120.0_fp, 120.0_fp, 100.0_fp, &
      !               EF6       EF7       EF8       EF9       EF10
         120.0_fp, 100.0_fp, 100.0_fp, 100.0_fp, 100.0_fp, &
      !               EF11      EF12      EF13      EF14      EF15
         100.0_fp, 2.0_fp,   2.0_fp,   2.0_fp,   2.0_fp/)

      ! Other monoterpenes, methanol, acetone, MBO are each 100% of their
      ! respective categories. The VOCs within the stress category each
      ! account for a specific fraction of emissions across all PFTs
      ! (ethene 58%, toluene 3%, HCN 1.5%). The VOCs within the
      ! other category also account for a given fraction of emissions
      ! across all PFTs (propene 48%, butene 24%, other alkenes 0.2%). But
      ! VOCs in the bidirectional category account for a different amount of
      ! the total flux for the different PFTs. So in this case we define a
      ! vector containing these fractions.

      ! Acetaldehyde: 40% of bidirectional category flux, except 25%
      ! for grasses and crops
      !                 EF1      EF2       EF3       EF4    EF5
      EM_FRAC_ALD2 = (/0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, &
      !                 EF6      EF7       EF8       EF9    EF10
         0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, &
      !                 EF11     EF12      EF13     EF14    EF15
         0.40_fp, 0.25_fp, 0.25_fp, 0.25_fp, 0.25_fp/)

      ! Ethanol: 40% of bidirectional category flux, except 25%
      ! for grasses and crops
      !                 EF1      EF2       EF3       EF4    EF5
      EM_FRAC_EOH  = (/0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, &
      !                 EF6      EF7       EF8       EF9    EF10
         0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, &
      !                 EF11     EF12      EF13     EF14    EF15
         0.40_fp, 0.25_fp, 0.25_fp, 0.25_fp, 0.25_fp/)

      ! Formic acid: 6% of bidirectional category flux, except 15%
      ! for grasses and crops
      !                 EF1      EF2       EF3       EF4    EF5
      EM_FRAC_FAXX = (/0.06_fp, 0.06_fp, 0.06_fp, 0.06_fp, 0.06_fp, &
      !                 EF6      EF7       EF8       EF9    EF10
         0.06_fp, 0.06_fp, 0.06_fp, 0.06_fp, 0.06_fp, &
      !                 EF11     EF12      EF13     EF14    EF15
         0.06_fp, 0.15_fp, 0.15_fp, 0.15_fp, 0.15_fp/)

      ! Acetic acid: 6% of bidirectional category flux, except 15%
      ! for grasses and crops
      !                 EF1      EF2       EF3       EF4    EF5
      EM_FRAC_AAXX = (/0.06_fp, 0.06_fp, 0.06_fp, 0.06_fp, 0.06_fp, &
      !                 EF6      EF7       EF8       EF9    EF10
         0.06_fp, 0.06_fp, 0.06_fp, 0.06_fp, 0.06_fp, &
      !                 EF11     EF12      EF13     EF14    EF15
         0.06_fp, 0.15_fp, 0.15_fp, 0.15_fp, 0.15_fp/)

      ! Formaldehyde: 8% of bidirectional category flux, except 20%
      ! for grasses and crops
      !                 EF1      EF2       EF3       EF4    EF5
      EM_FRAC_CH2O = (/0.08_fp, 0.08_fp, 0.08_fp, 0.08_fp, 0.08_fp, &
      !                 EF6      EF7       EF8       EF9    EF10
         0.08_fp, 0.08_fp, 0.08_fp, 0.08_fp, 0.08_fp, &
      !                 EF11     EF12      EF13     EF14    EF15
         0.08_fp, 0.20_fp, 0.20_fp, 0.20_fp, 0.20_fp/)

      !--------------------------------------------
      ! GET_GAMMAT_T_LI begins here!
      !--------------------------------------------
      AE = 0.0_fp
      ! Loop through plant types
      do P = 1, 15
         ! Add 1 to PFT_16 index to skip bare ground
         ARR_IND = P + 1
         ! Don't need to divide PFT_16 by 100 since it is already fraction
         select case ( TRIM(CMPD) )
            !try to calculate the seven read-in species online now
            ! ISOP: 100% of category
          case ('ISOP')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_ISOP(P)
            ! MBOX: 100% of category
          case ('MBOX')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_MBOX(P)
            ! BPIN: 100% of category
          case ('BPIN')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_BPIN(P)
            ! CARE: 100% of category
          case ('CARE')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_CARE(P)
            ! LIMO: 100% of category
          case ('LIMO')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_LIMO(P)
            ! OCIM: 100% of category
          case ('OCIM')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_OCIM(P)
            ! SABI: 100% of category
          case ('SABI')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_SABI(P)
            ! ---> Now compute EFs for a-pinene and myrcene as well
            ! a-pinene: 100% of category
          case ('APIN')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_APIN(P)
            ! Myrcene: 100% of category
          case ('MYRC')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_MYRC(P)
            ! Other monoterpenes: 100% of category
          case ('OMON')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_OMON(P)
            ! a-Farnesene: 100% of category
          case ('FARN')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_FARN(P)
            ! b-Caryophyllene: 100% of category
          case ('BCAR')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_BCAR(P)
            ! Other sesquiterpenes: 100% of category
          case ('OSQT')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_OSQT(P)
            ! Methanol: 100% of category
          case ('MOH')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_MOH(P)
            ! Acetone: 100% of category
          case ('ACET')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_ACET(P)
            ! Ethanol: variable fraction of category
          case ('EOH')
            AE = AE + PFT_16(ARR_IND) * EM_FRAC_EOH(P) * PFT_EF_BIDR(P)
            ! Formaldehyde: variable fraction of category
          case ('CH2O')
            AE = AE + PFT_16(ARR_IND) * EM_FRAC_CH2O(P) * PFT_EF_BIDR(P)
            ! Acetaldehyde: variable fraction of category
          case ('ALD2')
            AE = AE + PFT_16(ARR_IND) * EM_FRAC_ALD2(P) * PFT_EF_BIDR(P)
            ! Formic acid: variable fraction of category
          case ('FAXX')
            AE = AE + PFT_16(ARR_IND) * EM_FRAC_FAXX(P) * PFT_EF_BIDR(P)
            ! Acetic acid: variable fraction of category
          case ('AAXX')
            AE = AE + PFT_16(ARR_IND) * EM_FRAC_AAXX(P) * PFT_EF_BIDR(P)
            ! Ethene: 58% of "stress" category for all PFTs
          case ('C2H4')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_STRS(P) * 0.58_fp
            ! Toluene: 3% of "stress" category for all PFTs
          case ('TOLU')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_STRS(P) * 0.03_fp
            ! HCN: 1.5% of "stress" category for all PFTs
          case ('HCNX')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_STRS(P) * 0.015_fp
            ! Propene: 48% of "other" category for all PFTs
            ! Butene:  24% of "other" category for all PFTs
            ! Larger alkenes: 0.2% of "other" category for all PFTs
            ! Total: 72.2%
          case ('PRPE')
            AE = AE + PFT_16(ARR_IND) * PFT_EF_OTHR(P) * 0.722_fp
          case default
            RC = CC_FAILURE
            MSG = 'Invalid compound name'
            thisLoc = ' -> at CCPr_bvoc_Common (in process/bvoc/ccpr_bvoc_common_mod.F90)'
            call CC_Error( MSG, RC , thisLoc)
            return
         end select
      enddo

      ! Convert AEF arrays from [ug/m2/hr] to [kg/m2/s]
      FACTOR = 1.0e-9_fp / 3600.0_fp
      AE = AE * Factor
      ! Return w/ success
      RC = CC_SUCCESS
      return
   end subroutine CALC_AEF


end module CCPr_BVOC_Common_Mod
