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
   public :: hybrid_grid_supported
   public :: get_hybrid_ab
   public :: get_pedge
   public :: get_pmid

   !=========================================================================
   ! Hybrid-sigma vertical grid coefficients
   !
   ! Wet-air pressure at the bottom edge of level L is reconstructed from the
   ! surface pressure following GEOS-Chem (GeosUtil/pressure_mod.F90):
   !     Pedge(L) = Ap(L) + Bp(L) * Psurface        [L = 1 (surface) .. NZ+1 (top)]
   ! Bp is unitless. The stored Ap coefficients keep their native units (GEOS-Chem
   ! 72-level in hPa; GFS 127-level in Pa as read from the model output); the
   ! per-grid unit conversion to Pa is centralized in get_hybrid_ab, which always
   ! returns Ap in Pa. get_pedge therefore does plain Pa math (no *100). Add new
   ! resolutions by defining AP_<n>L/BP_<n>L below and extending the select cases
   ! in get_hybrid_ab and hybrid_grid_supported.
   !=========================================================================
   integer, parameter :: N_HYBRID_72L = 73
   real(fp), parameter :: AP_72L(N_HYBRID_72L) = [ &
      0.000000e+00_fp, 4.804826e-02_fp, 6.593752e+00_fp, 1.313480e+01_fp, &
      1.961311e+01_fp, 2.609201e+01_fp, 3.257081e+01_fp, 3.898201e+01_fp, &
      4.533901e+01_fp, 5.169611e+01_fp, 5.805321e+01_fp, 6.436264e+01_fp, &
      7.062198e+01_fp, 7.883422e+01_fp, 8.909992e+01_fp, 9.936521e+01_fp, &
      1.091817e+02_fp, 1.189586e+02_fp, 1.286959e+02_fp, 1.429100e+02_fp, &
      1.562600e+02_fp, 1.696090e+02_fp, 1.816190e+02_fp, 1.930970e+02_fp, &
      2.032590e+02_fp, 2.121500e+02_fp, 2.187760e+02_fp, 2.238980e+02_fp, &
      2.243630e+02_fp, 2.168650e+02_fp, 2.011920e+02_fp, 1.769300e+02_fp, &
      1.503930e+02_fp, 1.278370e+02_fp, 1.086630e+02_fp, 9.236572e+01_fp, &
      7.851231e+01_fp, 6.660341e+01_fp, 5.638791e+01_fp, 4.764391e+01_fp, &
      4.017541e+01_fp, 3.381001e+01_fp, 2.836781e+01_fp, 2.373041e+01_fp, &
      1.979160e+01_fp, 1.645710e+01_fp, 1.364340e+01_fp, 1.127690e+01_fp, &
      9.292942e+00_fp, 7.619842e+00_fp, 6.216801e+00_fp, 5.046801e+00_fp, &
      4.076571e+00_fp, 3.276431e+00_fp, 2.620211e+00_fp, 2.084970e+00_fp, &
      1.650790e+00_fp, 1.300510e+00_fp, 1.019440e+00_fp, 7.951341e-01_fp, &
      6.167791e-01_fp, 4.758061e-01_fp, 3.650411e-01_fp, 2.785261e-01_fp, &
      2.113490e-01_fp, 1.594950e-01_fp, 1.197030e-01_fp, 8.934502e-02_fp, &
      6.600001e-02_fp, 4.758501e-02_fp, 3.270000e-02_fp, 2.000000e-02_fp, &
      1.000000e-02_fp ]
   real(fp), parameter :: BP_72L(N_HYBRID_72L) = [ &
      1.000000e+00_fp, 9.849520e-01_fp, 9.634060e-01_fp, 9.418650e-01_fp, &
      9.203870e-01_fp, 8.989080e-01_fp, 8.774290e-01_fp, 8.560180e-01_fp, &
      8.346609e-01_fp, 8.133039e-01_fp, 7.919469e-01_fp, 7.706375e-01_fp, &
      7.493782e-01_fp, 7.211660e-01_fp, 6.858999e-01_fp, 6.506349e-01_fp, &
      6.158184e-01_fp, 5.810415e-01_fp, 5.463042e-01_fp, 4.945902e-01_fp, &
      4.437402e-01_fp, 3.928911e-01_fp, 3.433811e-01_fp, 2.944031e-01_fp, &
      2.467411e-01_fp, 2.003501e-01_fp, 1.562241e-01_fp, 1.136021e-01_fp, &
      6.372006e-02_fp, 2.801004e-02_fp, 6.960025e-03_fp, 8.175413e-09_fp, &
      0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, &
      0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, &
      0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, &
      0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, &
      0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, &
      0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, &
      0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, &
      0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, &
      0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, &
      0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, 0.000000e+00_fp, &
      0.000000e+00_fp ]

   !=========================================================================
   ! GFS 127-level hybrid-sigma coefficients (ak/bk global attributes from the
   ! GFS output NetCDF). These are stored VERBATIM in the file's native ordering,
   ! which runs top-of-atmosphere (index 1: bk=0) down to the surface (index 128:
   ! bk=1) -- the reverse of the GEOS-Chem surface-first convention. get_hybrid_ab
   ! reverses BOTH arrays together (they are paired per edge) so the returned
   ! coefficients follow the L=1 (surface) .. NZ+1 (top) convention. ak is already
   ! in Pa in the GFS output, so no unit scaling is applied for this grid.
   !=========================================================================
   integer, parameter :: N_HYBRID_127L = 128
   real(fp), parameter :: AK_127L(N_HYBRID_127L) = [ &
      0.999_fp, 1.605_fp, 2.532_fp, 3.924_fp, &
      5.976_fp, 8.947_fp, 13.177_fp, 19.096_fp, &
      27.243_fp, 38.276_fp, 52.984_fp, 72.293_fp, &
      97.269_fp, 129.11_fp, 169.135_fp, 218.767_fp, &
      279.506_fp, 352.894_fp, 440.481_fp, 543.782_fp, &
      664.236_fp, 803.164_fp, 961.734_fp, 1140.931_fp, &
      1341.538_fp, 1564.119_fp, 1809.028_fp, 2076.415_fp, &
      2366.252_fp, 2678.372_fp, 3012.51_fp, 3368.363_fp, &
      3745.646_fp, 4144.164_fp, 4563.881_fp, 5004.995_fp, &
      5468.017_fp, 5953.848_fp, 6463.864_fp, 7000.0_fp, &
      7563.494_fp, 8150.661_fp, 8756.529_fp, 9376.141_fp, &
      10004.55_fp, 10636.85_fp, 11268.16_fp, 11893.64_fp, &
      12508.52_fp, 13108.09_fp, 13687.73_fp, 14242.89_fp, &
      14769.15_fp, 15262.2_fp, 15717.86_fp, 16132.09_fp, &
      16501.02_fp, 16820.94_fp, 17088.32_fp, 17299.85_fp, &
      17453.08_fp, 17548.35_fp, 17586.77_fp, 17569.7_fp, &
      17498.7_fp, 17375.56_fp, 17202.3_fp, 16981.14_fp, &
      16714.5_fp, 16405.02_fp, 16055.49_fp, 15668.86_fp, &
      15248.25_fp, 14796.87_fp, 14318.04_fp, 13815.15_fp, &
      13291.63_fp, 12750.92_fp, 12196.47_fp, 11631.66_fp, &
      11059.83_fp, 10484.21_fp, 9907.927_fp, 9333.967_fp, &
      8765.155_fp, 8204.142_fp, 7653.387_fp, 7115.147_fp, &
      6591.468_fp, 6084.176_fp, 5594.876_fp, 5124.949_fp, &
      4675.554_fp, 4247.633_fp, 3841.918_fp, 3458.933_fp, &
      3099.01_fp, 2762.297_fp, 2448.768_fp, 2158.238_fp, &
      1890.375_fp, 1644.712_fp, 1420.661_fp, 1217.528_fp, &
      1034.524_fp, 870.778_fp, 725.348_fp, 597.235_fp, &
      485.392_fp, 388.734_fp, 306.149_fp, 236.502_fp, &
      178.651_fp, 131.447_fp, 93.74_fp, 64.392_fp, &
      42.274_fp, 26.274_fp, 15.302_fp, 8.287_fp, &
      4.19_fp, 1.994_fp, 0.81_fp, 0.232_fp, &
      0.029_fp, 0.0_fp, 0.0_fp, 0.0_fp ]
   real(fp), parameter :: BK_127L(N_HYBRID_127L) = [ &
      0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
      0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
      0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
      0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
      0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
      0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
      0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
      0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
      0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
      0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
      1.018e-05_fp, 8.141e-05_fp, 0.00027469_fp, 0.00065078_fp, &
      0.00127009_fp, 0.00219248_fp, 0.00347713_fp, 0.00518228_fp, &
      0.00736504_fp, 0.0100812_fp, 0.01338492_fp, 0.01732857_fp, &
      0.02196239_fp, 0.02733428_fp, 0.03348954_fp, 0.04047056_fp, &
      0.04831661_fp, 0.05706358_fp, 0.06674372_fp, 0.07738548_fp, &
      0.08900629_fp, 0.101594_fp, 0.1151262_fp, 0.1295762_fp, &
      0.1449129_fp, 0.1611008_fp, 0.1780999_fp, 0.195866_fp, &
      0.2143511_fp, 0.2335031_fp, 0.2532663_fp, 0.2735822_fp, &
      0.294389_fp, 0.3156229_fp, 0.337218_fp, 0.3591072_fp, &
      0.3812224_fp, 0.4034951_fp, 0.4258572_fp, 0.4482413_fp, &
      0.4705813_fp, 0.492813_fp, 0.5148743_fp, 0.5367062_fp, &
      0.5582525_fp, 0.5794605_fp, 0.6002815_fp, 0.6206707_fp, &
      0.6405875_fp, 0.6599957_fp, 0.6788633_fp, 0.6971631_fp, &
      0.714872_fp, 0.7319713_fp, 0.7484465_fp, 0.7642871_fp, &
      0.7794867_fp, 0.7940422_fp, 0.8079541_fp, 0.8212263_fp, &
      0.8338652_fp, 0.8458801_fp, 0.8572826_fp, 0.8680866_fp, &
      0.8783077_fp, 0.8879632_fp, 0.8970718_fp, 0.9056532_fp, &
      0.9137284_fp, 0.9213187_fp, 0.9284464_fp, 0.9351338_fp, &
      0.9414037_fp, 0.9472789_fp, 0.9527821_fp, 0.957936_fp, &
      0.962763_fp, 0.9672851_fp, 0.971524_fp, 0.9755009_fp, &
      0.9792364_fp, 0.9827508_fp, 0.9860625_fp, 0.9891851_fp, &
      0.9921299_fp, 0.9949077_fp, 0.9975282_fp, 1.0_fp ]

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

   !> \brief Report whether hybrid-sigma coefficients exist for a level count
   !! \param[in] nlev Number of vertical layers (edges = nlev + 1)
   !! \return .true. if Ap/Bp coefficients are defined for nlev
   !!
   !! Extend get_hybrid_ab (and this test) when adding new resolutions such as
   !! the GFS 127-level grid.
   pure function hybrid_grid_supported(nlev) result(supported)
      integer, intent(in) :: nlev
      logical :: supported
      select case (nlev)
       case (72, 127)
         supported = .true.
       case default
         supported = .false.
      end select
   end function hybrid_grid_supported

   !> \brief Return the hybrid-sigma Ap/Bp coefficients for a given level count
   !! \param[in]  nlev Number of vertical layers
   !! \param[out] ap   Ap coefficients [Pa], size nlev+1 (unallocated if unsupported)
   !! \param[out] bp   Bp coefficients [unitless], size nlev+1 (unallocated if unsupported)
   !! \param[out] ok   .true. when coefficients were returned for nlev
   !!
   !! Coefficients are returned in the surface-first convention: Ap(1)/Bp(1) is the
   !! surface edge and Ap(nlev+1)/Bp(nlev+1) the model top. Ap is always returned
   !! in Pa (the 72-level table is stored in hPa and scaled here; the GFS 127-level
   !! table is already in Pa). The GFS coefficients are stored top-to-surface in
   !! the file, so both arrays are reversed together here to match the convention.
   !! Add a new case here (and in hybrid_grid_supported) to support another grid.
   pure subroutine get_hybrid_ab(nlev, ap, bp, ok)
      integer, intent(in) :: nlev
      real(fp), allocatable, intent(out) :: ap(:)
      real(fp), allocatable, intent(out) :: bp(:)
      logical, intent(out) :: ok

      select case (nlev)
       case (72)
         ! Stored in hPa -> convert to Pa; already surface-first.
         ap = AP_72L * 100.0_fp
         bp = BP_72L
         ok = .true.
       case (127)
         ! GFS ak already in Pa; stored top-to-surface -> reverse both to
         ! surface-first (paired per edge, so reverse ak and bk identically).
         ap = AK_127L(N_HYBRID_127L:1:-1)
         bp = BK_127L(N_HYBRID_127L:1:-1)
         ok = .true.
       case default
         ok = .false.
      end select
   end subroutine get_hybrid_ab

   !> \brief Reconstruct wet-air edge pressures from surface pressure
   !! \param[in] ps   Surface pressure [Pa], shape (nx, ny)
   !! \param[in] nlev Number of vertical layers
   !! \return Edge pressure [Pa], shape (nx, ny, nlev+1); size-zero if nlev
   !!         is unsupported (guard callers with hybrid_grid_supported)
   !! \cite GEOS-Chem GeosUtil/pressure_mod.F90 (GET_PEDGE)
   !!
   !!   Pedge(i,j,L) = Ap(L) + Bp(L) * Psurface(i,j)
   !!
   !! Ap (from get_hybrid_ab) and ps are both in Pa, so the result is in Pa.
   pure function get_pedge(ps, nlev) result(pedge)
      real(fp), intent(in) :: ps(:,:)
      integer,  intent(in) :: nlev
      real(fp), allocatable :: pedge(:,:,:)

      real(fp), allocatable :: ap(:), bp(:)
      logical :: ok
      integer :: i, j, l, nx, ny

      call get_hybrid_ab(nlev, ap, bp, ok)
      if (.not. ok) then
         allocate(pedge(0, 0, 0))
         return
      end if

      nx = size(ps, 1)
      ny = size(ps, 2)
      allocate(pedge(nx, ny, nlev + 1))

      do l = 1, nlev + 1
         do j = 1, ny
            do i = 1, nx
               pedge(i, j, l) = ap(l) + bp(l) * ps(i, j)
            end do
         end do
      end do
   end function get_pedge

   !> \brief Layer mid-point pressure = mean of the bounding edge pressures
   !! \param[in] pedge Edge pressure [Pa], shape (nx, ny, nz+1)
   !! \return Mid-layer pressure [Pa], shape (nx, ny, nz)
   !! \cite GEOS-Chem GeosUtil/pressure_mod.F90 (GET_PCENTER)
   pure function get_pmid(pedge) result(pmid)
      real(fp), intent(in) :: pedge(:,:,:)
      real(fp), allocatable :: pmid(:,:,:)

      integer :: i, j, k, nx, ny, nz

      nx = size(pedge, 1)
      ny = size(pedge, 2)
      nz = size(pedge, 3) - 1
      allocate(pmid(nx, ny, nz))

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               pmid(i, j, k) = 0.5_fp * (pedge(i, j, k) + pedge(i, j, k + 1))
            end do
         end do
      end do
   end function get_pmid

end module met_utilities_mod
