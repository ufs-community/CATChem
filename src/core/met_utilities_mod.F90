!> \file met_utilities_mod.F90
!! \brief Meteorological utility functions for CATChem
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 1.0
!!
!! This module provides meteorological and atmospheric utility functions
!! commonly used in atmospheric chemistry and physics, including calculations
!! for potential temperature, virtual temperature, dew point, relative humidity,
!! saturation vapor pressure, and more.
!!
!! \details
!! The met_utilities module includes:
!! - Potential temperature calculation
!! - Virtual temperature calculation
!! - Dew point calculation
!! - Relative humidity calculation
!! - Saturation vapor pressure (Clausius-Clapeyron)
!! - Mixing ratio and specific humidity conversions
!! - Lapse rate calculations
!!
!! \section met_utilities_usage Usage Example
!! \code{.f90}
!! use met_utilities_mod
!! real(fp) :: T, p, theta, Tv, rh, Td, es
!! theta = potential_temperature(T, p, p0)
!! Tv = virtual_temperature(T, qv)
!! Td = dew_point(T, rh)
!! es = saturation_vapor_pressure(T)
!! \endcode
!!
module Met_Utilities_Mod
   use Precision_Mod
   use Constants
   implicit none
   private

   public :: potential_temperature
   public :: virtual_temperature
   public :: dew_point
   public :: relative_humidity
   public :: saturation_vapor_pressure
   public :: mixing_ratio
   public :: specific_humidity
   public :: dry_adiabatic_lapse_rate
   public :: bulk_richardson_number
   public :: monin_obukhov_length
   public :: friction_velocity
   public :: stability_classification
   public :: saturation_mixing_ratio
   public :: latent_heat_vaporization
   public :: psychrometric_constant
   public :: wind_profile_loglaw
   public :: brunt_vaisala_frequency
   public :: psi_m_businger
   public :: psi_h_businger
   public :: arrhenius_rate
   public :: henrys_law_constant
   public :: photolysis_rate_scaling
   public :: ppm_to_ugm3
   public :: ugm3_to_ppm
   public :: stokes_settling_velocity
   public :: cunningham_correction_factor
   public :: nuclear_decay
   public :: stokes_number
   public :: mean_free_path_air

contains

   !> \brief Calculate potential temperature (theta)
   !! \param[in] T Temperature [K]
   !! \param[in] p Pressure [Pa]
   !! \param[in] p0 Surface pressure [Pa]
   !! \return Potential temperature [K]
   !! \cite WallaceHobbs2006
   function potential_temperature(T, p, p0) result(theta)
      real(fp), intent(in) :: T, p, p0
      real(fp) :: theta
      theta = T * (p0 / p) ** (Rd / Cp)
   end function potential_temperature

   !> \brief Calculate virtual temperature
   !! \param[in] T Temperature [K]
   !! \param[in] qv Water vapor mixing ratio [kg/kg]
   !! \return Virtual temperature [K]
   !! \cite WallaceHobbs2006
   function virtual_temperature(T, qv) result(Tv)
      real(fp), intent(in) :: T, qv
      real(fp) :: Tv
      Tv = T * (1.0_fp + 0.61_fp * qv)
   end function virtual_temperature

   !> \brief Calculate dew point temperature
   !! \param[in] T Temperature [K]
   !! \param[in] rh Relative humidity [0-1]
   !! \return Dew point temperature [K]
   !! \cite Bolton1980
   function dew_point(T, rh) result(Td)
      real(fp), intent(in) :: T, rh
      real(fp) :: Td
      real(fp) :: es, ed
      es = saturation_vapor_pressure(T)
      ed = rh * es
      Td = 243.5_fp / (17.67_fp / log(ed / 611.2_fp) - 1.0_fp) + 273.15_fp
   end function dew_point

   !> \brief Calculate relative humidity
   !! \param[in] T Temperature [K]
   !! \param[in] qv Water vapor mixing ratio [kg/kg]
   !! \param[in] p Pressure [Pa]
   !! \return Relative humidity [0-1]
   !! \cite WallaceHobbs2006
   function relative_humidity(T, qv, p) result(rh)
      real(fp), intent(in) :: T, qv, p
      real(fp) :: rh
      real(fp) :: e, es
      e = qv * p / (0.622_fp + 0.378_fp * qv)
      es = saturation_vapor_pressure(T)
      rh = e / es
      ! Clip to physical limits
      rh = max(0.0_fp, min(1.0_fp, rh))
   end function relative_humidity

   !> \brief Calculate saturation vapor pressure (Clausius-Clapeyron)
   !! \param[in] T Temperature [K]
   !! \return Saturation vapor pressure [Pa]
   !! \cite Bolton1980
   function saturation_vapor_pressure(T) result(es)
      real(fp), intent(in) :: T
      real(fp) :: es
      es = 611.2_fp * exp(17.67_fp * (T - 273.15_fp) / (T - 29.65_fp))
   end function saturation_vapor_pressure

   !> \brief Calculate mixing ratio from specific humidity
   !! \param[in] q Specific humidity [kg/kg]
   !! \return Mixing ratio [kg/kg]
   function mixing_ratio(q) result(r)
      real(fp), intent(in) :: q
      real(fp) :: r
      r = q / (1.0_fp - q)
   end function mixing_ratio

   !> \brief Calculate specific humidity from mixing ratio
   !! \param[in] r Mixing ratio [kg/kg]
   !! \return Specific humidity [kg/kg]
   function specific_humidity(r) result(q)
      real(fp), intent(in) :: r
      real(fp) :: q
      q = r / (1.0_fp + r)
   end function specific_humidity

   !> \brief Calculate dry adiabatic lapse rate
   !! \return Dry adiabatic lapse rate [K/m]
   function dry_adiabatic_lapse_rate() result(gamma_d)
      real(fp) :: gamma_d
      gamma_d = g0 / Cp
   end function dry_adiabatic_lapse_rate

   !> \brief Calculate the bulk Richardson number
   !! \param[in] T0 Surface temperature [K]
   !! \param[in] Tz Temperature at height z [K]
   !! \param[in] u Wind speed at height z [m/s]
   !! \param[in] z Height above ground [m]
   !! \return Bulk Richardson number (dimensionless)
   function bulk_richardson_number(T0, Tz, u, z) result(Ri)
      real(fp), intent(in) :: T0, Tz, u, z
      real(fp) :: Ri
      if (u > 0.0_fp .and. z > 0.0_fp) then
         Ri = (g0 / T0) * (Tz - T0) * z / (u**2)
      else
         Ri = 0.0_fp
      endif
   end function bulk_richardson_number

   !> \brief Calculate the Monin-Obukhov length
   !! \param[in] ustar Friction velocity [m/s]
   !! \param[in] T0 Surface temperature [K]
   !! \param[in] H Sensible heat flux [W/m^2]
   !! \param[in] rho Air density [kg/m^3]
   !! \return Monin-Obukhov length [m]
   function monin_obukhov_length(ustar, T0, H, rho) result(L)
      real(fp), intent(in) :: ustar, T0, H, rho
      real(fp) :: L
      if (ustar > 0.0_fp .and. abs(H) > 0.0_fp) then
         L = - (ustar**3 * rho * Cp * T0) / (VON_KARMAN * g0 * H)
      else
         L = 1.0e5_fp  ! Neutral/very stable default
      endif
   end function monin_obukhov_length

   !> \brief Calculate friction velocity (u*)
   !! \param[in] tau Surface shear stress [N/m^2]
   !! \param[in] rho Air density [kg/m^3]
   !! \return Friction velocity [m/s]
   function friction_velocity(tau, rho) result(ustar)
      real(fp), intent(in) :: tau, rho
      real(fp) :: ustar
      if (rho > 0.0_fp) then
         ustar = sqrt(abs(tau) / rho)
      else
         ustar = 0.0_fp
      endif
   end function friction_velocity

   !> \brief Classify atmospheric stability based on Monin-Obukhov length
   !! \param[in] L Monin-Obukhov length [m]
   !! \return Stability class: -1 (unstable), 0 (neutral), 1 (stable)
   function stability_classification(L) result(class)
      real(fp), intent(in) :: L
      integer :: class
      if (L < -200.0_fp) then
         class = -1  ! Unstable
      else if (L > 200.0_fp) then
         class = 1   ! Stable
      else
         class = 0   ! Neutral
      endif
   end function stability_classification

   !> \brief Calculate saturation mixing ratio
   !! \param[in] p Pressure [Pa]
   !! \param[in] T Temperature [K]
   !! \return Saturation mixing ratio [kg/kg]
   function saturation_mixing_ratio(p, T) result(ws)
      real(fp), intent(in) :: p, T
      real(fp) :: ws
      real(fp) :: es
      es = saturation_vapor_pressure(T)
      ws = 0.622_fp * es / (p - es)
   end function saturation_mixing_ratio

   !> \brief Calculate latent heat of vaporization (temperature dependent)
   !! \param[in] T Temperature [K]
   !! \return Latent heat of vaporization [J/kg]
   function latent_heat_vaporization(T) result(Lv)
      real(fp), intent(in) :: T
      real(fp) :: Lv
      Lv = 2.501e6_fp - 2.361e3_fp * (T - 273.15_fp)
   end function latent_heat_vaporization

   !> \brief Calculate the psychrometric constant
   !! \param[in] p Pressure [Pa]
   !! \param[in] Lv Latent heat of vaporization [J/kg]
   !! \return Psychrometric constant [Pa/K]
   function psychrometric_constant(p, Lv) result(gamma)
      real(fp), intent(in) :: p, Lv
      real(fp) :: gamma
      gamma = Cp * p / (0.622_fp * Lv)
   end function psychrometric_constant

   !> \brief Calculate wind speed at height z using the log-law
   !! \param[in] ustar Friction velocity [m/s]
   !! \param[in] z Height above ground [m]
   !! \param[in] z0 Surface roughness length [m]
   !! \return Wind speed at height z [m/s]
   function wind_profile_loglaw(ustar, z, z0) result(u)
      real(fp), intent(in) :: ustar, z, z0
      real(fp) :: u
      if (z > z0 .and. z0 > 0.0_fp) then
         u = ustar / VON_KARMAN * log(z / z0)
      else
         u = 0.0_fp
      endif
   end function wind_profile_loglaw

   !> \brief Calculate Brunt–Väisälä frequency squared (N^2)
   !! \param[in] T0 Reference temperature [K]
   !! \param[in] dTdz Vertical temperature gradient [K/m]
   !! \return Brunt–Väisälä frequency squared [1/s^2]
   !! \cite WallaceHobbs2006
   function brunt_vaisala_frequency(T0, dTdz) result(N2)
      real(fp), intent(in) :: T0, dTdz
      real(fp) :: N2
      N2 = (g0 / T0) * (dTdz + g0 / Cp)
   end function brunt_vaisala_frequency

   !> \brief Businger-Dyer stability correction for momentum
   !! \param[in] zeta z/L (dimensionless stability parameter)
   !! \return Psi_m (stability correction for momentum)
   function psi_m_businger(zeta) result(psi_m)
      real(fp), intent(in) :: zeta
      real(fp) :: psi_m
      if (zeta < 0.0_fp) then
         psi_m = 2.0_fp * log((1.0_fp + sqrt(1.0_fp - 16.0_fp*zeta)) / 2.0_fp)
      else
         psi_m = -5.0_fp * zeta
      endif
   end function psi_m_businger

   !> \brief Businger-Dyer stability correction for heat
   !! \param[in] zeta z/L (dimensionless stability parameter)
   !! \return Psi_h (stability correction for heat)
   function psi_h_businger(zeta) result(psi_h)
      real(fp), intent(in) :: zeta
      real(fp) :: psi_h
      if (zeta < 0.0_fp) then
         psi_h = 2.0_fp * log((1.0_fp + sqrt(1.0_fp - 16.0_fp*zeta)) / 2.0_fp)
      else
         psi_h = -5.0_fp * zeta
      endif
   end function psi_h_businger

   !> \brief Calculate Arrhenius rate constant
   !! \param[in] A Pre-exponential factor [units vary]
   !! \param[in] Ea Activation energy [J/mol]
   !! \param[in] T Temperature [K]
   !! \return Rate constant [units of A]
   !! \cite SeinfeldPandis2016
   function arrhenius_rate(A, Ea, T) result(k)
      real(fp), intent(in) :: A, Ea, T
      real(fp) :: k
      real(fp), parameter :: R = 8.314462618_fp  ! Gas constant [J/mol/K]
      k = A * exp(-Ea / (R * T))
   end function arrhenius_rate

   !> \brief Calculate Henry's Law constant (temperature dependent)
   !! \param[in] H0 Reference Henry's constant [mol/(m^3*Pa)]
   !! \param[in] dH Enthalpy of solution [J/mol]
   !! \param[in] T Temperature [K]
   !! \param[in] T0 Reference temperature [K]
   !! \return Henry's Law constant at T [mol/(m^3*Pa)]
   !! \cite Sander2015
   function henrys_law_constant(H0, dH, T, T0) result(H)
      real(fp), intent(in) :: H0, dH, T, T0
      real(fp) :: H
      real(fp), parameter :: R = 8.314462618_fp
      H = H0 * exp(-dH/R * (1.0_fp/T - 1.0_fp/T0))
   end function henrys_law_constant

   !> \brief Scale photolysis rate for solar zenith angle
   !! \param[in] J0 Base photolysis rate [1/s]
   !! \param[in] sza Solar zenith angle [degrees]
   !! \return Scaled photolysis rate [1/s]
   function photolysis_rate_scaling(J0, sza) result(J)
      real(fp), intent(in) :: J0, sza
      real(fp) :: J
      J = J0 * max(0.0_fp, cos(sza * 3.141592653589793_fp / 180.0_fp))
   end function photolysis_rate_scaling

   !> \brief Convert ppm to ug/m3
   !! \param[in] ppm Concentration [ppm]
   !! \param[in] M Molar mass [g/mol]
   !! \param[in] T Temperature [K]
   !! \param[in] p Pressure [Pa]
   !! \return Concentration [ug/m3]
   function ppm_to_ugm3(ppm, M, T, p) result(ugm3)
      real(fp), intent(in) :: ppm, M, T, p
      real(fp) :: ugm3
      ugm3 = ppm * 1.0e-6_fp * p * M / (RSTARG * T) * 1.0e3_fp
   end function ppm_to_ugm3

   !> \brief Convert ug/m3 to ppm
   !! \param[in] ugm3 Concentration [ug/m3]
   !! \param[in] M Molar mass [g/mol]
   !! \param[in] T Temperature [K]
   !! \param[in] p Pressure [Pa]
   !! \return Concentration [ppm]
   function ugm3_to_ppm(ugm3, M, T, p) result(ppm)
      real(fp), intent(in) :: ugm3, M, T, p
      real(fp) :: ppm
      ppm = ugm3 * (RSTARG * T) / (p * M * 1.0e3_fp) * 1.0e6_fp
   end function ugm3_to_ppm

   !> \brief Calculate Stokes settling velocity for a particle
   !! \param[in] dp Particle diameter [m]
   !! \param[in] rho_p Particle density [kg/m3]
   !! \param[in] rho_a Air density [kg/m3]
   !! \param[in] mu Air dynamic viscosity [kg/m/s]
   !! \param[in] Cc Cunningham correction factor
   !! \return Settling velocity [m/s]
   function stokes_settling_velocity(dp, rho_p, rho_a, mu, Cc) result(vs)
      real(fp), intent(in) :: dp, rho_p, rho_a, mu, Cc
      real(fp) :: vs
      vs = (dp**2) * (rho_p - rho_a) * g0 * Cc / (18.0_fp * mu)
   end function stokes_settling_velocity

   !> \brief Calculate Cunningham correction factor
   !! \param[in] dp Particle diameter [m]
   !! \param[in] lambda Mean free path of air [m]
   !! \return Cunningham correction factor (dimensionless)
   function cunningham_correction_factor(dp, lambda) result(Cc)
      real(fp), intent(in) :: dp, lambda
      real(fp) :: Cc
      if (dp > 0.0_fp .and. lambda > 0.0_fp) then
         Cc = 1.0_fp + 2.0_fp * lambda / dp * (1.257_fp + 0.4_fp * exp(-1.1_fp * dp / lambda))
      else
         Cc = 1.0_fp
      endif
   end function cunningham_correction_factor

   !> \brief Calculate nuclear decay (first-order)
   !! \param[in] N0 Initial quantity
   !! \param[in] lambda Decay constant [1/s]
   !! \param[in] t Time [s]
   !! \return Remaining quantity after time t
   function nuclear_decay(N0, lambda, t) result(N)
      real(fp), intent(in) :: N0, lambda, t
      real(fp) :: N
      N = N0 * exp(-lambda * t)
   end function nuclear_decay


   !> \brief Calculate Stokes number from base state variables
   !! \param[in] rho_p Particle density [kg/m^3]
   !! \param[in] d_p Particle diameter [m]
   !! \param[in] U Characteristic velocity [m/s]
   !! \param[in] mu Dynamic viscosity [kg/m/s]
   !! \param[in] L Characteristic length scale [m]
   !! \return Stokes number (dimensionless)
   function stokes_number(rho_p, d_p, U, mu, L) result(Stk)
      real(fp), intent(in) :: rho_p, d_p, U, mu, L
      real(fp) :: Stk
      if (mu > 0.0_fp .and. L > 0.0_fp) then
         Stk = (rho_p * d_p**2 * U) / (18.0_fp * mu * L)
      else
         Stk = 0.0_fp
      endif
   end function stokes_number

   !> \brief Calculate the mean free path of air molecules
   !! \param[in] T Temperature [K]
   !! \param[in] p Pressure [Pa]
   !! \return Mean free path [m]
   !! \cite SeinfeldPandis2016
   function mean_free_path_air(T, p) result(lambda)
      real(fp), intent(in) :: T, p
      real(fp) :: lambda
      real(fp), parameter :: d_air = 3.7e-10_fp  ! Effective air molecule diameter [m]
      lambda = BOLTZ * T / (sqrt(2.0_fp) * 3.141592653589793_fp * d_air**2 * p)
   end function mean_free_path_air

end module met_utilities_mod
