!>
!! \file ccpr_wetdep_common_mod.F90
!! \brief Contains module ccpr_wetdep_common_mod
!!
!! \ingroup catchem_wetdep_process
!!
!! \author Wei Li
!! \date 04/2025
!!!>
module CCPr_Wetdep_Common_Mod
   use precision_mod, only: fp, ZERO, rae, TINY_
   use Error_Mod
   use constants
   implicit none
   private

   public  :: WetDepStateType
   public  :: rainout, rainout_loss
   public  :: washout, washout_loss
   public  :: complete_reevap

   !> \brief Type for CATCHem WetDep Process
   !!
   !! \details Contains all the information needed to run the CATChem WetDep Process
   !!
   !! This type contains the following variables:
   !! - Activate : Activate Process (True/False)
   !! - SchemeOpt : Scheme Option
   !! - WetDep_Flux : WetDep fluxes
   !!
   !! \ingroup catchem_wetdep_process
   !!!>
   TYPE :: WetDepStateType
      ! Generic Variables for Every Process
      Logical                         :: Activate               !< Activate Process (True/False)
      INTEGER                         :: SchemeOpt              !< Scheme Option
      ! Process Specific Parameters
      real, allocatable               :: WetDep_Flux(:,:)         !< WetDep fluxes

   END TYPE WetDepStateType

contains
   !>
   !! \brief Computes RAINFRAC, the fraction of soluble species lost to rainout events in precipitation.
   !!
   !! \param is_aero    aerosol rainout flag
   !! \param efficiency  temperature-dependent scale factor for rainout fraction
   !! \param wd_LidAndGas  ice-to-gas ratio is computed by co-condensation?
   !! \param k0        Henry's solubility constant [M/atm]
   !! \param cr        Henry's volatility constant [K]
   !! \param pKa       Henry's pH correction factor [1]
   !! \param cnvI2G    Conversion factor from ice to gas ratio if wd_LidAndGas is true
   !! \param retfac    Retention factor of species
   !! \param f         Fraction of grid box that is precipiting [unitless]
   !! \param k         Rainout rate constant [1/s]
   !! \param dt        time step [s]
   !! \param tk        temperature [K]
   !! \param c_h2o     Mix ratio of H2O [cm3 H2O/cm3 air]
   !! \param cldice    Precipitable cloud ice mixing ratio [cm3 ice/cm3 air]
   !! \param cldliq    Precipitable cloud liquid mixing ratio [cm3 H2O/cm3 air]
   !! \param spc       Species name
   !! \param lossfrac   Fraction of species lost to rainout [unitless]
   !! \param SO2       SO2 concentration [kg/kg]
   !! \param H2O2      H2O2 concentration [kg/kg]
   !!
   !! \ingroup catchem_wetdep_process
   !!!>
   subroutine rainout( is_aero, efficiency, wd_LidAndGas, k0, cr, pKa, cnvI2G, retfac, f, k, dt, tk, c_h2o, cldice, cldliq, spc, lossfrac, SO2, H2O2 )
      IMPLICIT NONE
      ! Parameters
      !-----------
      logical,    intent(in)  :: is_aero              !< aerosol rainout flag
      real(fp),   intent(in)  :: efficiency(3)        !< temperature-dependent scale factor for rainout fraction
      logical,    intent(in)  :: wd_LidAndGas         !< ice-to-gas ratio is computed by co-condensation?
      real(fp),   intent(in)  :: k0                   !< Henry's solubility constant [M/atm]
      real(fp),   intent(in)  :: cr                   !< Henry's volatility constant [K]
      real(fp),   intent(in)  :: pKa                  !< Henry's pH correction factor [1]
      real(fp),   intent(in)  :: cnvI2G               !< Conversion factor from ice to gas ratio if wd_LidAndGas is true
      real(fp),   intent(in)  :: retfac               !< Retention factor of species
      real(fp),   intent(in)  :: f                    !< Fraction of grid box that is precipiting [unitless]
      real(fp),   intent(in)  :: k                    !< Rainout rate constant [1/s]
      real(fp),   intent(in)  :: dt                   !< time step [s]
      real(fp),   intent(in)  :: tk                   !< temperature [K]
      real(fp),   intent(in)  :: c_h2o                !< Mix ratio of H2O [cm3 H2O/cm3 air]
      real(fp),   intent(in)  :: cldice               !< Precipitable cloud ice mixing ratio [cm3 ice/cm3 air]
      real(fp),   intent(in)  :: cldliq               !< Precipitable cloud liquid mixing ratio [cm3 H2O/cm3 air]
      character(len = 20),  intent(in)  :: spc        !< Species name
      real(fp),   intent(out) :: lossfrac             !< Fraction of species lost to rainout [unitless]
      real(fp),   intent(in) :: SO2                   !< SO2 concentration [kg/kg]
      real(fp),   intent(in) :: H2O2                  !< H2O2 concentration [kg/kg]

      ! Local Variables
      !----------------
      real(fp) :: i2g, l2g, c_tot
      real(fp) :: f_i, f_l, ki
      real(fp) :: SO2s, H2O2s   !after uint conversion to mol/mol
      real(fp) :: SO2LOSS
      real(fp), parameter :: kc = 5.e-3    ! conversion rate from cloud condensate to precip (s-1)

      !--------------------------------------------
      ! main function
      !--------------------------------------------
      !=================================================================
      ! %%% SPECIAL CASE %%%
      ! SO2, HNO3 and H2SO4 scavenges like an aerosol although they are
      ! considered to be a gas-phase species elsewhere
      !=================================================================
      if (is_aero .or. spc == 'SO2' .or. spc == 'HNO3' .or. spc == 'H2SO4') then
         lossfrac = rainfrac( f, k, dt )

         ! -- apply rainout efficiency (simplify from APPLY_RAINOUT_EFF)
         if (tk < 237.0_fp) then
            ! ice
            lossfrac = efficiency(1) * lossfrac
         else if (tk < 258.0_fp) then
            ! snow
            lossfrac = efficiency(2) * lossfrac
         else
            ! liquid rain
            lossfrac = efficiency(3) * lossfrac
         end if

         if (spc == 'SO2') then
            !unit conversion to mol/mol using airMW and SO2MW or H2OMW
            SO2s = SO2 * (AIRMW / 64.04_fp)
            H2O2s = H2O2 * (AIRMW / 34.02_fp)
            ! Update SO2 and H2O2
            if ( SO2s > TINY_ ) then
               ! Limit lossfrac
               SO2LOSS      = MIN( SO2s, H2O2s )
               lossfrac     = SO2LOSS * lossfrac / SO2s
               lossfrac     = MAX( lossfrac, 0e+0_fp )

               ! Update saved H2O2 concentration
               !TODO:: comment out the depletion since we may not have afterchem H2O2 like in GOES-Chem
               !       Otherwise, the depletion will do twice with the later one in washout_loss subroutine
               !H2O2s = H2O2s - ( SO2s * lossfrac )
               !H2O2s = MAX( H2O2s, TINY_ )

            else

               ! If SO2 is not defined then set lossfrac to 0
               lossfrac     = 0.0_fp

            endif

            ! Update saved SO2 concentration
            !TODO:: comment out the depletion since we may not have afterchem SO2 like in GOES-Chem
            !       Otherwise, the depletion will do twice with the later one in washout_loss subroutine
            !SO2s     = SO2s * ( 1.0_fp - lossfrac )
            !SO2s     = MAX( SO2s, TINY_ )
         end if

      else !for gases exceptfor SO2, HNO3 and H2SO4

         ! -- compute ice to gas ratio assuming scavenging by co-condensation
         !    simplified from COMPUTE_Ki
         ! wd_LidAndGas is true only for SO2, NH3 and H2O2 for now
         i2g = zero
         if (wd_LidAndGas) then
            if ( c_h2o > zero ) i2g = cnvI2G * cldice / c_h2o
         end if
         ! -- compute l2g (adopted from COMPUTE_L2G)
         !TODO: seems an error in GOCART version here
         l2g = liq_to_gas_ratio( k0, cr, pKa, tk, cldliq)

         ! -- fraction of species in liquid and ice phases
         c_tot = one + l2g + i2g
         f_l   = l2g / c_tot
         f_i   = i2g / c_tot

         ! -- compute Ki for loss due to scavenging from convective updraft
         if ( tk >= 268.0_fp ) then
            ki = kc * ( f_l + f_i )
         else if ( tk > 248.0_fp ) then
            ki = kc * ( retfac * f_l + f_i )
         else
            ki = kc * f_i
         end if

         ! -- compute rained-out fraction
         lossfrac = rainfrac( f, ki, dt )

      end if

   end subroutine rainout

   !>
   !! \brief Computes WASHFRAC, the fraction of soluble species lost to washout events in precipitation.
   !!
   !! \param radius         Particle radius (um)
   !! \param f              Fraction of grid box that is precipitating [unitless]
   !! \param tk             Temperature in grid box (K)
   !! \param qdwn           Instant precip rate in grid box (cm3 (H2O) / cm2 (air) / s)
   !! \param dz             Height of grid box [cm]
   !! \param dt             Timestep (s)
   !! \param spc            Species name
   !! \param is_aero        aerosol washout flag
   !! \param k0             Henry's solubility constant [M/atm]
   !! \param cr             Henry's volatility constant [K]
   !! \param pKa            Henry's pH correction factor [1]
   !! \param wtune          Washout tuning factor; newly added in GOCART version
   !! \param radius_thr     Fine/coarse particle radius threshold (um); using 1.0 um for now
   !! \param washfrac       Fraction of species lost to washout [unitless]
   !! \param kin            Kinetic process flag [kinetic or equilibrium]
   !! \param SO2            SO2 concentration [kg/kg]
   !! \param H2O2           H2O2 concentration [kg/kg]
   !!
   !! \ingroup catchem_wetdep_process
   !!!>

   subroutine washout( radius, f, tk, qdwn, dz, dt, spc, is_aero, k0, cr, pKa, wtune, radius_thr, washfrac, kin, SO2, H2O2)

      implicit none

      real(fp),  intent(in)  :: radius         !< Particle radius (um)
      real(fp),  intent(in)  :: f              !< Fraction of grid box that is precipitating [unitless]
      real(fp),  intent(in)  :: tk             !< Temperature in grid box (K)
      real(fp),  intent(in)  :: qdwn           !< Instant precip rate in grid box (cm3 (H2O) / cm2 (air) / s)
      real(fp),  intent(in)  :: dz             !< Height of grid box [cm]
      real(fp),  intent(in)  :: dt             !< Timestep (s)
      character(len = 20),  intent(in) :: spc  !< Species name
      logical,   intent(in)  :: is_aero        !< aerosol washout flag
      real(fp),  intent(in) :: k0              !< Henry's solubility constant [M/atm]
      real(fp),  intent(in) :: cr              !< Henry's volatility constant [K]
      real(fp),  intent(in) :: pKa             !< Henry's pH correction factor [1]
      real(fp),  intent(in)  :: wtune          !< Washout tuning factor; newly added in GOCART version
      real(fp),  intent(in)  :: radius_thr     !< Fine/coarse particle radius threshold (um); using 1.0 um for now
      real(fp),  intent(out) :: washfrac       !< Fraction of species lost to washout [unitless]
      logical,   intent(out)  :: kin           !< Kinetic process flag [kinetic or equilibrium]
      real(fp),  intent(in) :: H2O2            !< H2O2 conc [kg/kg]  TODO: conc's after aqueous rxns are applied. These are computed
      real(fp),  intent(in) :: SO2              !< SO2 conc [kg/kg]   in the sulfate chemistry module and passed here (not considered by now)
      ! -- local variables
      real(fp)               :: SO2LOSS
      real(fp)               :: SO2s, H2O2s !unit conversion to mol/mol

      ! -- begin
      washfrac = zero

      !=================================================================
      ! HNO3 scavenges like an aerosol although it is considered
      ! to be a gas-phase species elsewhere (e.g. dry deposition)
      !=================================================================

      IF ( Spc == 'HNO3' ) THEN  !TODO: better way to check species name?

         ! Washout is a kinetic process
         KIN      = .TRUE.

         ! Get washout fraction
         WASHFRAC = WASHFRAC_HNO3(  F, TK, QDWN,  DT )

         !=================================================================
         ! SO2 scavenges like an aerosol although it is considered
         ! to be a gas-phase species elsewhere (e.g. dry deposition)
         !=================================================================
      ELSE IF ( Spc == 'SO2' .or. spc == 'H2SO4' ) THEN  !TODO: better way to check species name?

         ! NOTE: Even though SO2 is not an aerosol we treat it as SO4 in
         ! wet scavenging.  When evaporation occurs, it returns to SO4.
         KIN      = .TRUE.
         !NOTE: Here we use a dummy radius = 0.5 um for SO2/H2SO4 since
         !      it is considered as fine aerosol here
         WASHFRAC = WASHFRAC_AEROSOL( 0.5_fp, F, TK, QDWN, DT, wtune, 1.0_fp )

         IF (Spc == 'SO2') THEN
            ! Use the wet-scavenging following [Chin et al, 1996] such
            ! that a soluble fraction of SO2 is limited by the availability
            ! of H2O2 in the precipitating grid box.  Then scavenge the
            ! soluble SO2 at the same rate as sulfate.

            !unit conversion to mol/mol using airMW and SO2MW or H2OMW
            SO2s = SO2 * (AIRMW / 64.04_fp)
            H2O2s = H2O2 * (AIRMW / 34.02_fp)
            IF ( TK >= 268e+0_fp .AND. SO2s > TINY_ ) THEN

               ! Adjust WASHFRAC
               SO2LOSS  = MIN( SO2s, H2O2s )
               WASHFRAC = SO2LOSS * WASHFRAC / SO2s
               WASHFRAC = MAX( WASHFRAC, 0e+0_fp )

               ! Deplete H2O2s the same as SO2s
               !TODO:: comment out the depletion since we may not have afterchem H2O2 like in GOES-Chem
               !       Otherwise, the depletion will do twice with the later one in washout_loss subroutine
               !H2O2s = H2O2s - ( SO2s * WASHFRAC )
               !H2O2s = MAX( H2O2s, TINY_ )

            ELSE
               WASHFRAC = 0e+0_fp

            ENDIF

            ! Update saved SO2 concentration
            !TODO:: comment out the depletion since we may not have afterchem SO2 like in GOES-Chem
            !       Otherwise, the depletion will do twice with the later one in washout_loss subroutine
            !SO2s = SO2s * ( 1e+0_fp - WASHFRAC )
            !SO2s = MAX( SO2s, TINY_ )
         END IF

         !-----------------------------------------------------------------
         ! Washout for aerosol species
         !-----------------------------------------------------------------
      ELSE IF ( is_aero ) THEN

         ! Washout is a kinetic process
         KIN      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL( radius, F, TK, QDWN, DT, wtune, radius_thr )

         !-----------------------------------------------------------------
         ! Washout for gas-phase species
         ! (except H2SO4, NO3, and SO2, which scavenge like aerosols;
         !-----------------------------------------------------------------
      ELSE
         ! Washout is an equilibrium process
         KIN = .FALSE.
         CALL WASHFRAC_LIQ_GAS( F, TK, QDWN, DZ, DT, K0, CR, PKA, WASHFRAC, KIN )

      END IF

   end subroutine washout

   !>
   !! \brief Computes the fraction of species lost to rainout.
   !!
   !! \param f          Fraction of grid box that is precipiting [unitless]
   !! \param k          Rainout rate constant [1/s]
   !! \param dt         time step [s]
   !!
   !! \ingroup catchem_wetdep_process
   !!!>
   real(fp) function rainfrac( f, k, dt )

      implicit none
      ! INPUT Parameters
      real(fp), intent(in) :: f          !< Fraction of grid box that is precipiting [unitless]
      real(fp), intent(in) :: k          !< Rainout rate constant [1/s]
      real(fp), intent(in) :: dt         !< time step [s]

      rainfrac = f * ( one - exp( -k * dt ) )

   end function rainfrac

   !>
   !! \brief Computes the ratio L2G = Cliq / Cgas, which is the mixing ratio
   !!  of species in the liquid phase, divided by the mixing ratio of species in the gas phase.
   !!
   !! \param k0         Henry's solubility constant [M/atm]
   !! \param cr         Henry's volatility constant [K]
   !! \param pKa        Henry's pH correction factor [1]
   !! \param tk         Temperature [K]
   !! \param qliq       Liquid water content [cm3 H2O/cm3 air]
   !!
   !! \ingroup catchem_wetdep_process
   !!!>
   real(fp) function liq_to_gas_ratio( k0, cr, pKa, tk, qliq )

      real(fp), intent(in) :: k0      !< Henry's solubility constant [M/atm]
      real(fp), intent(in) :: cr      !< Henry's volatility constant [K]
      real(fp), intent(in) :: pKa     !< Henry's pH correction factor [1]
      real(fp), intent(in) :: tk      !< Temperature [K]
      real(fp), intent(in) :: qliq    !< Liquid water content [cm3 H2O/cm3 air]

      ! -- local variables
      real(fp) :: h !, t
      ! -- local parameters
      real(fp), parameter :: cloud_pH = 4.5_fp
      real(fp), parameter  :: Tref = 298.15_fp     ! K
      real(fp), parameter  :: R    = 8.3144598_fp  ! J K-1 mol-1
      real(fp), parameter  :: Pref = 101.325_fp    ! mPa

      ! -- compute Henry's law constant for a given temperature
      !t = real(tk, kind=f8) TODO: not sure why this conversion is needed; do not use it for now
      h = k0 * exp( cr * (1._fp/tk - 1._fp/Tref) ) * R * tk / Pref

      ! -- adjust Henry's law constant for chemical equilibriums in liquid phase
      if ( pKa > -100._fp ) h = h * ( one + 10._fp ** ( cloud_pH - pKa ) )

      liq_to_gas_ratio = h * qliq

   end function liq_to_gas_ratio

   !>
   !! \brief Computes the fraction of soluble aerosol species lost to washout.
   !!
   !! \param radius       Particle radius (um)
   !! \param f            Washout fraction [unitless]
   !! \param tk           Temperature in grid box (K)
   !! \param pdwn         Instant precip rate in grid box (cm3 (H2O) / cm2 (air) / s)
   !! \param dt           Timestep (s)
   !! \param tuning       Washout tuning factor; newly added in GOCART version
   !! \param radius_fine  Fine particle radius threshold (um); using 1.0 um for now
   !!
   !! \ingroup catchem_wetdep_process
   !!!>
   real(fp) function washfrac_aerosol( radius, f, tk, pdwn, dt, tuning, radius_fine )

      implicit none

      real(fp), intent(in) :: radius         !< particle radius (um)
      real(fp), intent(in) :: f              !< washout fraction [unitless]
      real(fp), intent(in) :: tk             !< Temperature in grid box (K)
      real(fp), intent(in) :: pdwn           !< Instant precip rate in grid box (cm3 (H2O) / cm2 (air) / s)
      real(fp), intent(in) :: dt             !< Timestep (s)
      real(fp), intent(in) :: tuning         !< Washout tuning factor; newly added in GOCART version
      real(fp), intent(in) :: radius_fine    !< fine particle radius threshold (um); using 1.0 um for now

      ! -- local variables
      real(fp)          :: dth, pph
      ! -- local parameters
      real(fp), parameter :: k_wash = 1.06e-03_fp
      real(fp), parameter :: h2s = 3600.0_fp ! s-1

      ! -- begin
      washfrac_aerosol = zero

      if ( f > zero ) then
         ! -- convert instant rates (s-1) to hourly rates
         pph = 10. * pdwn * h2s
         dth = dt / h2s

         if ( radius < radius_fine ) then  !for fine aerosol (simplified from WASHFRAC_FINE_AEROSOL)
            if ( tk >= 268e+0_fp  ) then
               washfrac_aerosol = F * ( one  - EXP(-k_wash * tuning * (pph / f ) ** 0.61e+0_fp * dth))
            else
               washfrac_aerosol = F * ( one  - EXP(-2.6e+1_fp * k_wash * tuning  * (pph / f ) ** 0.96e+0_fp * dth))
            endif
         else  !for coarse aerosol (simplified from WASHFRAC_COARSE_AEROSOL)
            if ( tk >= 268e+0_fp  ) then
               washfrac_aerosol = F * ( one  - EXP(-0.92e+0_fp * tuning * (pph / f ) ** 0.79e+0_fp * dth))
            else
               !TODO: GOCART applied a factor of 0.5 to the tuning factor for coarse aerosol????
               !washfrac_aerosol = F * ( one  - EXP(-1.57e+0_fp / 0.5e+0_fp * tuning * (pph / f ) ** 0.96e+0_fp * dth))
               washfrac_aerosol = F * ( one  - EXP(-1.57e+0_fp * tuning * (pph / f ) ** 0.96e+0_fp * dth))
            endif
         endif
      endif

   end function washfrac_aerosol

   !>
   !! \brief Computes the fraction of HNO3 species lost to washout.
   !!
   !! \param f          Fraction of grid box that is precipitating [unitless]
   !! \param tk         Temperature in grid box (K)
   !! \param pdwn       Instant precip rate in grid box (cm3 (H2O) / cm2 (air) / s)
   !! \param dt         Timestep (s)
   !!
   !! \ingroup catchem_wetdep_process
   !!!>
   real(fp) function washfrac_hno3( f, tk, pdwn, dt )

      implicit none

      real(fp), intent(in) :: f      !< Fraction of grid box that is precipitating [unitless]
      real(fp), intent(in) :: tk     !< Temperature in grid box [K]
      real(fp), intent(in) :: pdwn   !< Precip rate thru bottom of grid box (cm3 (H2O) / cm2 (air) / s)
      real(fp), intent(in) :: dt     !< Timestep of washout event (s)

      ! -- local parameters
      real(fp), parameter :: k_wash = 1.0_fp    ! First order washout rate (cm-1)

      ! -- begin
      washfrac_hno3 = zero
      ! -- compute washout fraction only if T >= 268K
      if ( tk >= 268.0_fp .and. f > zero ) then
         washfrac_hno3 = f * ( one - exp( -k_wash * pdwn * dt / f ) )
      end if

   end function washfrac_hno3

   !>
   !! \brief Computes the fraction of soluble liquid/gas phase species lost to washout.
   !!
   !! \param f          Fraction of grid box that is precipitating [unitless]
   !! \param tk         Temperature in grid box (K)
   !! \param pdwn       Instant precip rate in grid box (cm3 (H2O) / cm2 (air) / s)
   !! \param dz         Height of grid box [cm]
   !! \param dt         Timestep of washout event (s)
   !! \param k0         Henry's solubility constant [M/atm]
   !! \param cr         Henry's volatility constant [K]
   !! \param pKa        Henry's pH correction factor [1]
   !! \param washfrac   Fraction of species lost to washout [unitless]
   !! \param kin        Kinetic process flag [kinetic or equilibrium]
   !!
   !! \ingroup catchem_wetdep_process
   !!!>
   subroutine washfrac_liq_gas( f, tk, pdwn, dz, dt, k0, cr, pKa, washfrac, kin)

      implicit none

      real(fp),  intent(in) :: f            !< Fraction of grid box that is precipitating [unitless]
      real(fp),  intent(in) :: tk           !< Temperature in grid box [K]
      real(fp),  intent(in) :: pdwn         !< Precip rate thru bottom of grid box (cm3 (H2O) / cm2 (air) / s)
      real(fp),  intent(in) :: dz           !< Height of grid box [cm]
      real(fp),  intent(in) :: dt           !< Timestep of washout event (s)
      real(fp),  intent(in) :: k0           !< Henry's solubility constant [M/atm]
      real(fp),  intent(in) :: cr           !< Henry's volatility constant [K]
      real(fp),  intent(in) :: pKa          !< Henry's pH correction factor [1]
      real(fp),  intent(out) :: washfrac    !< Fraction of species lost to washout [unitless]
      logical,   intent(out) :: kin         !< Kinetic process flag [kinetic or equilibrium]

      ! -- local variables
      real(fp) :: qliq, l2g, washfrac_kin

      ! -- begin

      ! Start with the assumption that washout will be an
      ! equilibrium process
      kin = .false.

      if ( tk < 268.0_fp ) then
         ! -- no washout
         washfrac = zero

      else

         ! -- compute L2G
         qliq = pdwn * dt / ( f * dz )

         ! Compute liquid to gas ratio
         l2g = liq_to_gas_ratio( k0, cr, pKa, tk, qliq )

         ! -- washout fraction from Henry's Law
         washfrac = l2g / ( one + l2g )

         ! -- washout fraction from kinetic processes (HNO3)
         ! set f = one and call washfrac_hno3 function above
         washfrac_kin = washfrac_hno3( one, tk, pdwn, dt )

         ! -- equilibrium washout must not exceed kinetic washout
         !TODO: GOCART is missing  'washfrac_kin * f' here ??????
         if ( washfrac > washfrac_kin ) then
            washfrac = washfrac_kin * f
            kin = .true. ! washout is a kinetic process
         end if

      end if

   end subroutine washfrac_liq_gas


   !>
   !! \brief Computes the concentrations of species lost to rainout.
   !!
   !! \param k          Layer index
   !! \param lossfrac   Fraction of species lost to washout [unitless]
   !! \param conc       Concentration [kg/m2]
   !! \param dconc      Concentration loss [kg/m2]
   !!
   !! \ingroup catchem_wetdep_process
   !!!>
   subroutine rainout_loss( k, lossfrac, conc, dconc )

      implicit none

      ! -- input/output parameters
      integer,   intent(in)    :: k                     !< layer index
      real(fp),  intent(in)    :: lossfrac              !< fraction of species lost to rainout [unitless]
      real(fp),  dimension(:), intent(inout) :: conc    !< concentration [kg/m2]
      real(fp),  dimension(:), intent(inout) :: dconc   !< concentration loss kg/m2

      ! -- local variables
      real(fp)   :: wetloss  !< wet loss concentration

      ! -- apply loss (simplified from the DO_WASHOUT_ONLY subroutine)
      wetloss = lossfrac * conc(k)
      conc(k) = conc(k) - wetloss
      if ( k == size(conc) ) then !If it is the top layer; we assume the layer index is not reversed
         ! Dconc is an accumulator array for rained-out species.
         ! The species in Dconc are in the liquid phase and will
         ! precipitate to the levels below until a washout occurs.
         dconc(k) = wetloss
      else
         ! Add to Dconc the species lost to rainout in grid box
         ! (I,J,L) plus the species lost to rainout from grid box
         ! (I,J,L+1), which has by now precipitated down into
         ! grid box (I,J,L).  Dconc will continue to accumulate
         ! rained out species in this manner until a washout
         ! event occurs.
         dconc(k) = dconc(k+1) + wetloss
      end if

   end subroutine rainout_loss

   !>
   !! \brief Computes the concentrations of species lost to washout.
   !!
   !! \param k            layer index
   !! \param lossfrac     fraction of species lost to washout [unitless]
   !! \param kin          kinetic process flag [kinetic or equilibrium]
   !! \param f_washout    washout fraction [unitless]
   !! \param f_rainout    rainout fraction [unitless]
   !! \param pdwn         downward flux of precipitation
   !! \param reevap       Precip forming or evaporating [cm3 (h2o)/cm3 (air)]
   !! \param delz_cm      vertical grid spacing [cm]
   !! \param conc         concentration [kg/m2]
   !! \param dconc        concentration loss kg/m2
   !! \param SO4          SO4 concentration [kg/m2]
   !! \param spc          Species name
   !!
   !! \ingroup catchem_wetdep_process
   !!!>
   subroutine washout_loss( k, lossfrac, kin, f_washout, f_rainout, pdwn, reevap, delz_cm, conc, dconc, spc, SO4 )

      implicit none

      ! -- input/output parameters
      integer,   intent(in)    :: k                     !< layer index
      real(fp),  intent(inout) :: lossfrac              !< fraction of species lost to rainout [unitless]
      logical,   intent(in)    :: kin                   !< kinetic process flag [kinetic or equilibrium]
      real(fp),  intent(in)    :: f_washout             !< washout fraction [unitless]
      real(fp),  intent(in)    :: f_rainout             !< rainout fraction [unitless]
      real(fp),  dimension(:), intent(in)    :: pdwn    !< downward flux of precipitation
      !real(fp), dimension(:), intent(in)    :: qq      !< precipatitng water rate [cm3 (h2o) / cm2 (air) / s]
      real(fp)                 :: reevap                !< Precip forming or evaporating [cm3 (h2o)/cm3 (air)]
      real(fp),  dimension(:), intent(in)    :: delz_cm !< vertical grid spacing [cm]
      real(fp),  dimension(:), intent(inout) :: conc    !< concentration [kg/m2]
      real(fp),  dimension(:), intent(inout) :: dconc   !< concentration loss kg/m2
      real(fp),  dimension(:), intent(inout) :: SO4     !< SO4 concentration [kg/m2]
      character(len = 20),  intent(in) :: spc           !< Species name

      ! -- local variables
      integer    :: km1      !< upper one layer index
      real(fp)   :: f        !< washout + rainout fraction [unitless]
      real(fp)   :: wetloss  !< wet loss concentration
      real(fp)   :: alpha    !< re-evaporate fraction [1]
      real(fp)   :: gain     !< washout gain [kg/m2]
      real(fp)   :: washed   !< washout concentration [kg/m2]

      !begins here (simplified from the DO_WASHOUT_ONLY and DO_WASHOUT_AT_SFC subroutines)
      f = f_washout + f_rainout
      km1 = k + 1
      if ( k == 1 ) then !for the surface layer; note this assumes the layer index starts from surface not the top

         ! -- f is included in lossfrac for aerosols and HNO3
         if ( kin ) then
            wetloss = lossfrac * conc(k)
         else
            wetloss = f * lossfrac * conc(k)
         end if
         conc (k) = conc (k) - wetloss
         dconc(k) = dconc(km1) + wetloss

      else  !for middle layers
         if ( kin ) then
            ! -- adjust loss fraction for aerosols
            lossfrac = lossfrac * f_washout / f

            ! Define ALPHA, the fraction of the raindrops that
            ! re-evaporate when falling from (I,J,L+1) to (I,J,L)
            if ( pdwn(km1) - ZERO > ZERO  ) then !avoid divide by zero
               !TODO: is qq(k) right in here?
               !alpha = abs( qq(k) ) * delz_cm(k) / pdwn(km1)
               alpha = abs( reevap ) * delz_cm(k) / pdwn(km1)
            else
               alpha = one
            end if
            ! Restrict ALPHA to be less than 1
            alpha = min( one, alpha )
            ! Assume 50% of the re-evaporated water rains out to aerosols
            ! GAINED is the rained out aerosol coming down from
            ! grid box (I,J,L+1) that will evaporate and re-enter
            ! the atmosphere in the gas phase in grid box (I,J,L).
            gain  = 0.5 * alpha * dconc(km1)
            wetloss  = conc(k) * lossfrac - gain
            ! SO2 in sulfate chemistry is wet-scavenged on the
            ! raindrop and converted to SO4 by aqeuous chem.
            ! If evaporation occurs then SO2 comes back as SO4
            if (spc == 'SO2') then
               SO4(k) = SO4(k) + gain * 96e+0_fp / 64e+0_fp
               conc(k) = conc(k) * (1e+0_fp - lossfrac)
            else
               conc(k) = conc(k) - wetloss
            end if

         else !not kinetic process (mainly for gases)

            !Not sure why dconc is not counted here in wetloss calculation
            washed  = f_washout * conc(k) + dconc(km1)
            wetloss = lossfrac * ( washed - dconc(km1) )
            conc(k) = conc(k) - wetloss

         end if ! kinetic or equilibrium process

         ! Add washout losses in grid box (I,J,L) to dconc
         if ( f_rainout > zero ) then
            dconc(k) = dconc(k) + wetloss
         else
            dconc(k) = dconc(km1) + wetloss
         end if

      end if !surface or middle layer

   end subroutine washout_loss

   !>
   !! \brief Computes the re-evaporation all of the soluble species back into the atmosphere.
   !!
   !! \param k       layer index
   !! \param conc    concentration [kg/m2]
   !! \param dconc   concentration loss [kg/m2]
   !! \param spc     Species name
   !! \param SO4    SO4 concentration [kg/m2]
   !!
   !! \ingroup catchem_wetdep_process
   !!!>
   subroutine complete_reevap( k, conc, dconc, spc, SO4 )

      implicit none

      ! -- input/output parameters
      integer,   intent(in)    :: k                     !< layer index
      real(fp),  dimension(:), intent(inout) :: conc    !< concentration [kg/m2]
      real(fp),  dimension(:), intent(inout) :: dconc   !< concentration loss kg/m2
      character(len = 20),  intent(in) :: spc           !< Species name
      real(fp),  dimension(:), intent(inout)    :: SO4  !< SO4 concentration [kg/m2]

      ! -- local variables
      integer    :: km1      !< upper one layer index
      real(fp)   :: wetloss !< wet loss

      ! -- begin
      km1 = k + 1
      wetloss = -dconc(km1)
      ! All of the rained-out species coming from grid box
      ! (I,J,L+1) goes back into the gas phase at (I,J,L)
      ! In evap, SO2 comes back as SO4
      if (spc == 'SO2') then
         SO4(k) = SO4(k) - (wetloss * 96e+0_fp / 64e+0_fp )
      else
         conc(k) = conc(k) - wetloss
      end if

      ! There is nothing rained out/washed out in grid box
      ! (I,J,L), so set dconc at grid box (I,J,L) to zero.
      dconc(k) = 0e+0_fp

   end subroutine complete_reevap


end module CCPr_WetDep_Common_Mod
