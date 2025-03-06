!>
!! \file
!! \brief Contains MEGAN2.1 biogenic VOC emission Scheme based on HEMCO and Sam Silva's canopy edits
!!
!! Reference:
!! (1) HEMCO's HCOX_MEGAN_MOD module (https://github.com/geoschem/HEMCO), which is
!!     based on Guenther et al, (GMD 2012) and associated MEGANv2.1 source code
!!     https://doi.org/10.5194/gmd-5-1471-2012
!! (2) Sam Silva's simplified canopy edits (GMD 2020) are also added
!!     https://doi.org/10.5194/gmd-13-2569-2020
!!
!! \author Wei Li
!! \date 07/2024
!! \ingroup catchem_bvoc_process
!!!>

module CCPr_Scheme_MeganV21_Mod

   implicit none

   private

   public :: CCPr_Scheme_MeganV21

contains

   !> \brief Brief description of the subroutine
   !!
   !! \param MeganSpecName     Name of megan species
   !! \param LAI               leaf area index
   !! \param PMISOLAI          LAI of previous month
   !! \param PFT_16(:)         Plant functional type fraction
   !! \param Q_DIR_2           Surface downwelling par diffuse flux
   !! \param Q_DIFF_2          Surface downwelling par beam flux
   !! \param PARDR_LASTXDAYS   Avg. PARDF of last NUM_DAYS
   !! \param PARDF_LASTXDAYS   Avg. PARDR of last NUM_DAYS
   !! \param TS                Surface temperature
   !! \param T_LASTXDAYS       Avg. temperature of last NUM_DAYS
   !! \param T_LAST24H         Avg. temperature of last 24 hours
   !! \param GWETROOT          Root zone soil moisture
   !! \param CO2Inhib          Turn on CO2 inhibition?
   !! \param CO2conc           CO2 concentrations
   !! \param SUNCOS            Cosine of solar zenith angle
   !! \param LAT               Latitude
   !! \param DOY               Day of year
   !! \param LocalHour         Local hour
   !! \param D_BTW_M           Days between mid-months
   !!!! \param AEF_ISOP          Emission factor of ISOP read from file
   !!!! \param AEF_MBOX          Emission factor of MBOX read from file
   !!!! \param AEF_BPIN          Emission factor of BPIN read from file
   !!!! \param AEF_CARE          Emission factor of CARE read from file
   !!!! \param AEF_LIMO          Emission factor of LIMO read from file
   !!!! \param AEF_OCIM          Emission factor of OCIM read from file
   !!!! \param AEF_SABI          Emission factor of SABI read from file
   !! \param EmisPerSpec       Emission per Species
   !! \param RC                Success or Failure
   !!
   !! Note that other state types may be required, e.g. one specific to the process group.
   !!!>
   subroutine CCPr_Scheme_MeganV21( &
      MeganSpecName,       &
      EmisPerSpec,         &
      LAI,                 &
      PFT_16,              &
      PMISOLAI,            &
      Q_DIR_2,             &
      Q_DIFF_2,            &
      PARDR_LASTXDAYS,     &
      PARDF_LASTXDAYS,     &
      TS,                  &
      T_LASTXDAYS,         &
      T_LAST24H,           &
      GWETROOT,            &
      CO2Inhib,            &
      CO2conc,             &
      SUNCOS,              &
      LAT,                 &
      DOY,                 &
      LocalHour,           &
      D_BTW_M,             &
      RC)

      ! Uses
      USE Constants,     Only : PI_180 ! Example to pull in a constant from the CONSTANTS MODULE < Modify as needed >
      use precision_mod, only : fp     ! Example to pull in a precision from the PRECISION MODULE < Modify as needed >
      Use Error_Mod,     Only : CC_SUCCESS    ! Error Check Success
      USE CCPr_BVOC_Common_Mod

      IMPLICIT NONE

      ! Arguments
      !integer,           intent(in)     :: nMeganSpec            !< Number of Megan Species
      !character(len=10), intent(in)     :: MeganSpecName(:)      !< name of megan species
      character(len=10), intent(in)     :: MeganSpecName      !< name of megan species
      real(fp),          intent(in)     :: LAI                   !< leaf area index
      real(fp),          intent(in)     :: PMISOLAI              !< LAI of previous month
      real(fp),          intent(in)     :: PFT_16(:)             !< plant functional type fraction
      real(fp),          intent(in)     :: Q_DIR_2               !< surface downwelling par diffuse flux
      real(fp),          intent(in)     :: Q_DIFF_2              !< surface downwelling par beam flux
      real(fp),          intent(in)     :: PARDR_LASTXDAYS       !< Avg. PARDF of last NUM_DAYS
      real(fp),          intent(in)     :: PARDF_LASTXDAYS       !< Avg. PARDR of last NUM_DAYS
      real(fp),          intent(in)     :: TS                    !< surface temperature
      real(fp),          intent(in)     :: T_LASTXDAYS           !< Avg. temperature of last NUM_DAYS
      real(fp),          intent(in)     :: T_LAST24H             !< Avg. temperature of last 24 hours
      real(fp),          intent(in)     :: GWETROOT              !< Root zone soil moisture
      logical,           intent(in)     :: CO2Inhib              !< turn on CO2 inhibition?
      real(fp),          intent(in)     :: CO2conc               !< CO2 concentrations
      real(fp),          intent(in)     :: SUNCOS                !< Cosine of solar zenith angle
      real(fp),          intent(in)     :: LAT                   !< Latitude
      integer,           intent(in)     :: DOY                   !< Day of year
      real(fp),          intent(in)     :: LocalHour             !< Local hour
      real(fp),          intent(in)     :: D_BTW_M               !< Days between mid-months
      !These seven EFs are changed to be read online now instead of read from file
      !real(fp),          intent(in)     :: AEF_ISOP              !< Emission factor of ISOP read from file
      !real(fp),          intent(in)     :: AEF_MBOX              !< Emission factor of MBOX read from file
      !real(fp),          intent(in)     :: AEF_BPIN              !< Emission factor of BPIN read from file
      !real(fp),          intent(in)     :: AEF_CARE              !< Emission factor of CARE read from file
      !real(fp),          intent(in)     :: AEF_LIMO              !< Emission factor of LIMO read from file
      !real(fp),          intent(in)     :: AEF_OCIM              !< Emission factor of OCIM read from file
      !real(fp),          intent(in)     :: AEF_SABI              !< Emission factor of SABI read from file
      !real(fp),          intent(inout)  :: EmisPerSpec(:)       !< Emission per Species
      real(fp),          intent(inout)  :: EmisPerSpec           !< Emission per Species
      !real(fp),          intent(inout)  :: TotalEmis            !< Total Emission (TODO: not used by now)
      integer,           intent(out)    :: RC                    ! Success or Failure

      ! Local Variables
      character(len=256) :: errMsg
      character(len=256) :: thisLoc

      real(fp)            :: MEGAN_EMIS ! emission for each species [kg/m2/s]
      character(len=256)  :: CMPD      ! Compound name
      REAL(fp)            :: GAMMA_LAI
      REAL(fp)            :: GAMMA_AGE
      REAL(fp)            :: GAMMA_TP  !canopy add
      REAL(fp)            :: CDEA(5)   !canopy add
      REAL(fp)            :: VPGWT(5)  !canopy add
      REAL(fp)            :: GAMMA_PAR_Sun, GAMMA_PAR_Shade !canopy add
      REAL(fp)            :: GAMMA_T_LD_Sun, GAMMA_T_LD_Shade !canopy add
      REAL(fp)            :: GAMMA_T_LI_Sun, GAMMA_T_LI_Shade !canopy add
      !REAL(fp)            :: WINDSP !canopy add
      REAL(fp)            :: GAMMA_PAR
      REAL(fp)            :: GAMMA_T_LD
      REAL(fp)            :: GAMMA_T_LI
      REAL(fp)            :: GAMMA_SM
      REAL(fp)            :: GAMMA_CO2
      REAL(fp)            :: AEF
      !REAL(fp)            :: D_BTW_M
      !REAL(fp)            :: TS, SUNCOS
      !REAL(fp)            :: Q_DIR_2, Q_DIFF_2
      REAL(fp)            :: BETA, LDF, CT1, CEO
      REAL(fp)            :: ANEW, AGRO, AMAT, AOLD
      REAL(fp)            :: MISOLAI, PMISOLAI_NORM
      REAL(fp)            :: PFTSUM
      REAL(fp), parameter :: LAI_MAX = 6.0_fp !Maximum LAI value [cm2/cm2]
      REAL(fp), parameter :: D2RAD  = PI_180  !Degrees to radians
      REAL(fp), parameter :: RAD2D  = 1.0_fp / PI_180  !Radians to degrees
      !REAL(fp)            :: LAT, LocalHour !canopy add
      REAL(fp)            :: PSTD
      REAL(fp)            :: Ea1L, Ea2L, SINbeta, SunF !canopy add
      LOGICAL             :: BIDIR
      INTEGER             :: K !, S, DOY  !canopy add and below
      REAL(fp)            :: T_Leaf_Int_Sun(5)
      REAL(fp)            :: T_Leaf_Int_Shade(5)
      REAL(fp)            :: T_Leaf_Temp_Sun(5)
      REAL(fp)            :: T_Leaf_Temp_Shade(5)
      !REAL(fp)            :: T_Leaf_Wind_Sun(5)
      !REAL(fp)            :: T_Leaf_Wind_Shade(5)
      REAL(fp)            :: P_Leaf_Int_Sun(5)
      REAL(fp)            :: P_Leaf_Int_Shade(5)
      REAL(fp)            :: P_Leaf_LAI_Sun(5)
      REAL(fp)            :: P_Leaf_LAI_Shade(5)
      REAL(fp)            :: Distgauss(5)

      ! Initialize parameters, gamma values, and return value
      errMsg = ''
      thisLoc = ' -> at CCPr_Scheme_MeganV21 (in CCPr_Scheme_MeganV21_mod.F90)'
      RC = CC_SUCCESS

      CDEA       = 0.0_fp  !canopy add
      GAMMA_TP   = 0.0_fp  !canopy add
      MEGAN_EMIS = 0.0_fp
      GAMMA_LAI  = 0.0_fp
      GAMMA_AGE  = 0.0_fp
      GAMMA_T_LD = 0.0_fp
      GAMMA_T_LI = 0.0_fp
      GAMMA_PAR  = 0.0_fp
      GAMMA_SM   = 0.0_fp
      GAMMA_CO2  = 0.0_fp
      BETA       = 0.0_fp
      AEF        = 0.0_fp
      LDF        = 0.0_fp
      CT1        = 0.0_fp
      CEO        = 0.0_fp
      ANEW       = 0.0_fp
      AGRO       = 0.0_fp
      AMAT       = 0.0_fp
      AOLD       = 0.0_fp
      BIDIR      = .FALSE.

      !----------------------------------
      ! Begin SchemeCCPr_Scheme_MeganV21
      !----------------------------------

      EmisPerSpec = 0.0_fp

      !-----------------------------------------------------
      ! Only interested in terrestrial biosphere
      ! If ( local LAI > 0 ) replace the zeros assigned above
      !-----------------------------------------------------
      if ( LAI > 0.0_fp ) then

         !-----------------normalize LAI by total PFT fractions
         PFTSUM = SUM( PFT_16(2:16) )
         MISOLAI  = min(LAI/PFTSUM, LAI_MAX)
         PMISOLAI_NORM= min(PMISOLAI/PFTSUM, LAI_MAX)

         !----------------- %%gamma values not related to compound%% ------------------

         ! --------------------------------------------------
         ! GAMMA_par (light activity factor)
         ! --------------------------------------------------

         ! Calculate GAMMA PAR only during day
         IF ( SUNCOS > 0.0_fp ) THEN

            call GET_GAMMA_PAR_PCEEA( Q_DIR_2,             &
               Q_DIFF_2,              &
               PARDR_LASTXDAYS,       &
               PARDF_LASTXDAYS,       &
               LAT, DOY,              &
               LocalHour,             &
               D2RAD, RAD2D,          &
               GAMMA_PAR)
         ELSE

            ! If night
            GAMMA_PAR = 0.0_fp
         ENDIF

         ! --------------------------------------------------
         ! CO2 inhibition of isoprene (Tai, Jan 2013)
         ! --------------------------------------------------
         IF ( CO2Inhib ) THEN
            call GET_GAMMA_CO2( CO2conc, GAMMA_CO2 )
         ELSE
            GAMMA_CO2 = 1.0_fp
         ENDIF

         !Sam Silva's canopy related coefficients
         T_Leaf_Int_Sun  = (/-13.891_fp, -12.322_fp, -1.032_fp, -5.172_fp, -5.589_fp/)
         T_Leaf_Int_Shade = (/-12.846_fp, -11.343_fp, -1.068_fp,-5.551_fp, -5.955_fp/)
         T_Leaf_Temp_Sun = (/1.064_fp, 1.057_fp, 1.031_fp,  1.050_fp, 1.051_fp/)
         T_Leaf_Temp_Shade = (/1.060_fp, 1.053_fp, 1.031_fp,1.051_fp, 1.052_fp/)
         P_Leaf_Int_Sun  = (/1.0831_fp, 1.0964_fp, 1.1036_fp, 1.0985_fp, 1.0901_fp/)
         P_Leaf_Int_Shade = (/0.8706_fp, 0.8895_fp, 0.9160_fp,0.9407_fp, 0.9564_fp/)
         P_Leaf_LAI_Sun = (/0.0018_fp, -0.1281_fp, -0.2977_fp, -0.4448_fp, -0.5352_fp/)
         P_Leaf_LAI_Shade = (/0.0148_fp, -0.1414_fp, -0.3681_fp,-0.5918_fp, -0.7425_fp/)
         VPGWT = (/0.1184635, 0.2393144, 0.284444444, 0.2393144, 0.1184635/)
         Distgauss = (/0.0469101, 0.2307534, 0.5, 0.7692465, 0.9530899/)

         call  SOLAR_ANGLE(DOY, LocalHour, LAT, D2RAD, SINbeta)
         call GET_CDEA(MISOLAI, CDEA  )

         !--------------------- %%gamma values related to compound%% ----------------------

         !DO S=1, nMeganSpec

         !CMPD = MeganSpecName(S)
         CMPD = MeganSpecName

         ! --------------------------------------------
         ! Get MEGAN parameters for this compound
         ! --------------------------------------------
         CALL GET_MEGAN_PARAMS ( CMPD, BETA, LDF,  CT1,  CEO,      &
            ANEW, AGRO, AMAT, AOLD, BIDIR, RC )

         ! --------------------------------------------------
         ! GAMMA_LAI (leaf area index activity factor)
         ! --------------------------------------------------
         call GET_GAMMA_LAI( MISOLAI, BIDIR, GAMMA_LAI )

         ! --------------------------------------------------
         ! GAMMA_AGE (leaf age activity factor)
         ! --------------------------------------------------
         call GET_GAMMA_AGE(  MISOLAI,           &
            PMISOLAI_NORM,          &
            D_BTW_M,           &
            T_LASTXDAYS,       &
            ANEW, AGRO, AMAT, AOLD,     &
            GAMMA_AGE)

         ! --------------------------------------------------
         ! GAMMA_T_LI (temperature activity factor for
         ! light-independent fraction)
         ! --------------------------------------------------
         !GAMMA_T_LI = GET_GAMMA_T_LI( TS, BETA )

         ! --------------------------------------------------
         ! GAMMA_T_LD (temperature activity factor for
         ! light-dependent fraction)
         ! --------------------------------------------------
         !GAMMA_T_LD = GET_GAMMA_T_LD( TS, Inst%T_LASTXDAYS(I,J), &
         !                             Inst%T_LAST24H(I,J), CT1, CEO )

         ! --------------------------------------------------
         ! Sam Silva's edits to replace GAMMA_T_LD
         !  and GAMMA_T_LI above
         ! --------------------------------------------------
         GAMMA_TP = 0.0_fp

         DO K = 1, 5

            call Calc_Sun_Frac(MISOLAI,SINbeta,Distgauss(K), SunF)

            PSTD = 200_fp
            call GET_GAMMA_PAR_C(Q_DIR_2,             &
               Q_DIFF_2,            &
               PARDR_LASTXDAYS,     &
               PARDF_LASTXDAYS,     &
               P_Leaf_LAI_Sun(K),   &
               P_Leaf_Int_Sun(K),   &
               MISOLAI, PSTD,       &
               GAMMA_PAR_Sun)

            PSTD = 50_fp
            call GET_GAMMA_PAR_C(Q_DIR_2,             &
               Q_DIFF_2,            &
               PARDR_LASTXDAYS,     &
               PARDF_LASTXDAYS,     &
               P_Leaf_LAI_Shade(K), &
               P_Leaf_Int_Shade(K), &
               MISOLAI, PSTD,       &
               GAMMA_PAR_Shade)

            call GET_GAMMA_T_LD_C( TS,                &
            !   T_LASTXDAYS,         &
               T_LAST24H,           &
               CT1, CEO,            &
               T_Leaf_Int_Sun(K),   &
               T_Leaf_Temp_Sun(K),  &
               GAMMA_T_LD_Sun )

            call GET_GAMMA_T_LD_C(TS,                 &
            !   T_LASTXDAYS,         &
               T_LAST24H,           &
               CT1, CEO,            &
               T_Leaf_Int_Shade(K), &
               T_Leaf_Temp_Shade(K),&
               GAMMA_T_LD_Shade )

            call GET_GAMMA_T_LI( TS, BETA,            &
               T_Leaf_Int_Sun(K),   &
               T_Leaf_Temp_Sun(K),  &
               GAMMA_T_LI_Sun )

            call GET_GAMMA_T_LI( TS, BETA,            &
               T_Leaf_Int_Shade(K), &
               T_Leaf_Temp_Shade(K),&
               GAMMA_T_LI_Shade )

            Ea1L  =  CDEA(K) * GAMMA_PAR_Sun * GAMMA_T_LD_Sun * SunF +   &
               GAMMA_PAR_Shade * GAMMA_T_LD_Shade * (1-SunF)

            Ea2L =  GAMMA_T_LI_Sun * SunF +                              &
               GAMMA_T_LI_Shade * (1-SunF)

            GAMMA_TP  = GAMMA_TP +                                       &
               (Ea1L*LDF + Ea2L*(1-LDF))* VPGWT(K)

            !!For test only
            !if ( K==1 ) then
            !write(*,*)'My test1',CMPD, CDEA(1), VPGWT(1),Distgauss(1),GAMMA_PAR_Sun,GAMMA_PAR_Shade,&
            !            GAMMA_T_LD_Shade,GAMMA_T_LD_Sun,GAMMA_T_LI_Shade,GAMMA_T_LI_Sun
            !endif

         ENDDO


         ! --------------------------------------------------
         ! GAMMA_SM (soil moisture activity factor)
         ! --------------------------------------------------
         call GET_GAMMA_SM( GWETROOT, CMPD, GAMMA_SM )

         ! --------------------------------------------------
         ! emission factor (Note: EFs of these seven species are now calculated online instead of reading from file)
         ! --------------------------------------------------
         !select case ( TRIM(CMPD) )
         ! case ('ISOP')
         !   AEF = AEF_ISOP
         ! case ('MBOX')
         !   AEF = AEF_MBOX
         ! case ('BPIN')
         !   AEF = AEF_BPIN
         ! case ('CARE')
         !   AEF = AEF_CARE
         ! case ('LIMO')
         !   AEF = AEF_LIMO
         ! case ('OCIM')
         !   AEF = AEF_OCIM
         ! case ('SABI')
         !   AEF = AEF_SABI
         ! case default !others are calculated inline
         call CALC_AEF(PFT_16, CMPD, AEF, RC)
         !end select

         ! --------------------------------------------------
         ! calculate emission
         ! --------------------------------------------------
         ! Emission is the product of all of these in kg/m2/s.
         ! Normalization factor ensures product of GAMMA values is 1.0 under
         !  standard conditions. Norm_FAC = 0.21. canopy add
         IF ( TRIM(CMPD) == 'ISOP' ) THEN
            ! Only apply CO2 inhibition to isoprene
            ! MEGAN_EMIS = Inst%NORM_FAC(1) * AEF * GAMMA_AGE * GAMMA_SM * &
            !              GAMMA_LAI * ((1.0_fp - LDF) * GAMMA_T_LI +      &
            !              (LDF * GAMMA_PAR * GAMMA_T_LD)) * GAMMA_CO2
            MEGAN_EMIS = MISOLAI * AEF * GAMMA_AGE * GAMMA_SM *  &
               GAMMA_TP*GAMMA_CO2*GAMMA_LAI*0.21_fp
         ELSE
            ! MEGAN_EMIS = Inst%NORM_FAC(1) * AEF * GAMMA_AGE * GAMMA_SM * &
            !              GAMMA_LAI * ((1.0_fp - LDF) * GAMMA_T_LI +      &
            !              (LDF * GAMMA_PAR * GAMMA_T_LD))
            MEGAN_EMIS = MISOLAI * AEF * GAMMA_AGE * GAMMA_SM *  &
               GAMMA_TP * GAMMA_LAI * 0.21_fp

         ENDIF

         !!!For test only
         !write(*,*) 'my test2: ', CMPD, MISOLAI, AEF, GAMMA_AGE, GAMMA_SM, GAMMA_TP, GAMMA_CO2, GAMMA_LAI

         !EmisPerSpec(S) = MEGAN_EMIS
         EmisPerSpec = MEGAN_EMIS

         !ENDDO !each species

      endif

      return

   end subroutine CCPr_Scheme_MeganV21

end module CCPr_Scheme_MeganV21_Mod
