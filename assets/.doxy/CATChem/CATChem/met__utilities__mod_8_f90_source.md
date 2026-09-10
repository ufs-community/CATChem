

# File met\_utilities\_mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**met\_utilities\_mod.F90**](met__utilities__mod_8_f90.md)

[Go to the documentation of this file](met__utilities__mod_8_f90.md)


```Fortran

module met_utilities_mod
   use precision_mod
   use constants
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
   public :: vertical_interp_pressure

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

   function potential_temperature(T, p, p0) result(theta)
      real(fp), intent(in) :: T, p, p0
      real(fp) :: theta
      theta = t * (p0 / p) ** (rd / cp)
   end function potential_temperature

   function virtual_temperature(T, qv) result(Tv)
      real(fp), intent(in) :: T, qv
      real(fp) :: Tv
      tv = t * (1.0_fp + 0.61_fp * qv)
   end function virtual_temperature

   function dew_point(T, rh) result(Td)
      real(fp), intent(in) :: T, rh
      real(fp) :: Td
      real(fp) :: es, ed
      es = saturation_vapor_pressure(t)
      ed = rh * es
      td = 243.5_fp / (17.67_fp / log(ed / 611.2_fp) - 1.0_fp) + 273.15_fp
   end function dew_point

   function relative_humidity(T, qv, p) result(rh)
      real(fp), intent(in) :: T, qv, p
      real(fp) :: rh
      real(fp) :: e, es
      e = qv * p / (0.622_fp + 0.378_fp * qv)
      es = saturation_vapor_pressure(t)
      rh = e / es
      ! Clip to physical limits
      rh = max(0.0_fp, min(1.0_fp, rh))
   end function relative_humidity

   function saturation_vapor_pressure(T) result(es)
      real(fp), intent(in) :: T
      real(fp) :: es
      es = 611.2_fp * exp(17.67_fp * (t - 273.15_fp) / (t - 29.65_fp))
   end function saturation_vapor_pressure

   function mixing_ratio(q) result(r)
      real(fp), intent(in) :: q
      real(fp) :: r
      r = q / (1.0_fp - q)
   end function mixing_ratio

   function specific_humidity(r) result(q)
      real(fp), intent(in) :: r
      real(fp) :: q
      q = r / (1.0_fp + r)
   end function specific_humidity

   function dry_adiabatic_lapse_rate() result(gamma_d)
      real(fp) :: gamma_d
      gamma_d = g0 / cp
   end function dry_adiabatic_lapse_rate

   function bulk_richardson_number(T0, Tz, u, z) result(Ri)
      real(fp), intent(in) :: T0, Tz, u, z
      real(fp) :: Ri
      if (u > 0.0_fp .and. z > 0.0_fp) then
         ri = (g0 / t0) * (tz - t0) * z / (u**2)
      else
         ri = 0.0_fp
      endif
   end function bulk_richardson_number

   function monin_obukhov_length(ustar, T0, H, rho) result(L)
      real(fp), intent(in) :: ustar, T0, H, rho
      real(fp) :: L
      if (ustar > 0.0_fp .and. abs(h) > 0.0_fp) then
         l = - (ustar**3 * rho * cp * t0) / (von_karman * g0 * h)
      else
         l = 1.0e5_fp  ! Neutral/very stable default
      endif
   end function monin_obukhov_length

   function friction_velocity(tau, rho) result(ustar)
      real(fp), intent(in) :: tau, rho
      real(fp) :: ustar
      if (rho > 0.0_fp) then
         ustar = sqrt(abs(tau) / rho)
      else
         ustar = 0.0_fp
      endif
   end function friction_velocity

   function stability_classification(L) result(class)
      real(fp), intent(in) :: L
      integer :: class
      if (l < -200.0_fp) then
         class = -1  ! Unstable
      else if (l > 200.0_fp) then
         class = 1   ! Stable
      else
         class = 0   ! Neutral
      endif
   end function stability_classification

   function saturation_mixing_ratio(p, T) result(ws)
      real(fp), intent(in) :: p, T
      real(fp) :: ws
      real(fp) :: es
      es = saturation_vapor_pressure(t)
      ws = 0.622_fp * es / (p - es)
   end function saturation_mixing_ratio

   function latent_heat_vaporization(T) result(Lv)
      real(fp), intent(in) :: T
      real(fp) :: Lv
      lv = 2.501e6_fp - 2.361e3_fp * (t - 273.15_fp)
   end function latent_heat_vaporization

   function psychrometric_constant(p, Lv) result(gamma)
      real(fp), intent(in) :: p, Lv
      real(fp) :: gamma
      gamma = cp * p / (0.622_fp * lv)
   end function psychrometric_constant

   function wind_profile_loglaw(ustar, z, z0) result(u)
      real(fp), intent(in) :: ustar, z, z0
      real(fp) :: u
      if (z > z0 .and. z0 > 0.0_fp) then
         u = ustar / von_karman * log(z / z0)
      else
         u = 0.0_fp
      endif
   end function wind_profile_loglaw

   function brunt_vaisala_frequency(T0, dTdz) result(N2)
      real(fp), intent(in) :: T0, dTdz
      real(fp) :: N2
      n2 = (g0 / t0) * (dtdz + g0 / cp)
   end function brunt_vaisala_frequency

   function psi_m_businger(zeta) result(psi_m)
      real(fp), intent(in) :: zeta
      real(fp) :: psi_m
      if (zeta < 0.0_fp) then
         psi_m = 2.0_fp * log((1.0_fp + sqrt(1.0_fp - 16.0_fp*zeta)) / 2.0_fp)
      else
         psi_m = -5.0_fp * zeta
      endif
   end function psi_m_businger

   function psi_h_businger(zeta) result(psi_h)
      real(fp), intent(in) :: zeta
      real(fp) :: psi_h
      if (zeta < 0.0_fp) then
         psi_h = 2.0_fp * log((1.0_fp + sqrt(1.0_fp - 16.0_fp*zeta)) / 2.0_fp)
      else
         psi_h = -5.0_fp * zeta
      endif
   end function psi_h_businger

   function arrhenius_rate(A, Ea, T) result(k)
      real(fp), intent(in) :: A, Ea, T
      real(fp) :: k
      real(fp), parameter :: R = 8.314462618_fp  ! Gas constant [J/mol/K]
      k = a * exp(-ea / (r * t))
   end function arrhenius_rate

   function henrys_law_constant(H0, dH, T, T0) result(H)
      real(fp), intent(in) :: H0, dH, T, T0
      real(fp) :: H
      real(fp), parameter :: R = 8.314462618_fp
      h = h0 * exp(-dh/r * (1.0_fp/t - 1.0_fp/t0))
   end function henrys_law_constant

   function photolysis_rate_scaling(J0, sza) result(J)
      real(fp), intent(in) :: J0, sza
      real(fp) :: J
      j = j0 * max(0.0_fp, cos(sza * 3.141592653589793_fp / 180.0_fp))
   end function photolysis_rate_scaling

   function ppm_to_ugm3(ppm, M, T, p) result(ugm3)
      real(fp), intent(in) :: ppm, M, T, p
      real(fp) :: ugm3
      ugm3 = ppm * 1.0e-6_fp * p * m / (rstarg * t) * 1.0e3_fp
   end function ppm_to_ugm3

   function ugm3_to_ppm(ugm3, M, T, p) result(ppm)
      real(fp), intent(in) :: ugm3, M, T, p
      real(fp) :: ppm
      ppm = ugm3 * (rstarg * t) / (p * m * 1.0e3_fp) * 1.0e6_fp
   end function ugm3_to_ppm

   function stokes_settling_velocity(dp, rho_p, rho_a, mu, Cc) result(vs)
      real(fp), intent(in) :: dp, rho_p, rho_a, mu, Cc
      real(fp) :: vs
      vs = (dp**2) * (rho_p - rho_a) * g0 * cc / (18.0_fp * mu)
   end function stokes_settling_velocity

   function cunningham_correction_factor(dp, lambda) result(Cc)
      real(fp), intent(in) :: dp, lambda
      real(fp) :: Cc
      if (dp > 0.0_fp .and. lambda > 0.0_fp) then
         cc = 1.0_fp + 2.0_fp * lambda / dp * (1.257_fp + 0.4_fp * exp(-1.1_fp * dp / lambda))
      else
         cc = 1.0_fp
      endif
   end function cunningham_correction_factor

   function nuclear_decay(N0, lambda, t) result(N)
      real(fp), intent(in) :: N0, lambda, t
      real(fp) :: N
      n = n0 * exp(-lambda * t)
   end function nuclear_decay


   function stokes_number(rho_p, d_p, U, mu, L) result(Stk)
      real(fp), intent(in) :: rho_p, d_p, U, mu, L
      real(fp) :: Stk
      if (mu > 0.0_fp .and. l > 0.0_fp) then
         stk = (rho_p * d_p**2 * u) / (18.0_fp * mu * l)
      else
         stk = 0.0_fp
      endif
   end function stokes_number

   function mean_free_path_air(T, p) result(lambda)
      real(fp), intent(in) :: T, p
      real(fp) :: lambda
      real(fp), parameter :: d_air = 3.7e-10_fp  ! Effective air molecule diameter [m]
      lambda = boltz * t / (sqrt(2.0_fp) * 3.141592653589793_fp * d_air**2 * p)
   end function mean_free_path_air

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

   pure subroutine get_hybrid_ab(nlev, ap, bp, ok)
      integer, intent(in) :: nlev
      real(fp), allocatable, intent(out) :: ap(:)
      real(fp), allocatable, intent(out) :: bp(:)
      logical, intent(out) :: ok

      select case (nlev)
       case (72)
         ! Stored in hPa -> convert to Pa; already surface-first.
         ap = ap_72l * 100.0_fp
         bp = bp_72l
         ok = .true.
       case (127)
         ! GFS ak already in Pa; stored top-to-surface -> reverse both to
         ! surface-first (paired per edge, so reverse ak and bk identically).
         ap = ak_127l(n_hybrid_127l:1:-1)
         bp = bk_127l(n_hybrid_127l:1:-1)
         ok = .true.
       case default
         ok = .false.
      end select
   end subroutine get_hybrid_ab

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

   subroutine vertical_interp_pressure(src_p, src_data, dst_p, dst_data)
      real(fp), intent(in)  :: src_p(:,:,:)
      real(fp), intent(in)  :: src_data(:,:,:)
      real(fp), intent(in)  :: dst_p(:,:,:)
      real(fp), intent(out) :: dst_data(:,:,:)

      integer  :: i, j, k, l, ll, nx, ny, nsrc, ndst
      real(fp) :: pt, p1, p2, w

      nx   = size(dst_data, 1)
      ny   = size(dst_data, 2)
      ndst = size(dst_data, 3)
      nsrc = size(src_p, 3)

      do j = 1, ny
         do i = 1, nx
            do k = 1, ndst
               pt = dst_p(i, j, k)

               if (nsrc == 1) then
                  dst_data(i, j, k) = src_data(i, j, 1)
                  cycle
               end if

               ! Locate the source interval [l, l+1] whose pressures bracket pt.
               l = 0
               bracket_search: do ll = 1, nsrc - 1
                  p1 = src_p(i, j, ll)
                  p2 = src_p(i, j, ll + 1)
                  if ((pt - p1) * (pt - p2) <= 0.0_fp) then
                     l = ll
                     exit bracket_search
                  end if
               end do bracket_search

               if (l == 0) then
                  ! Outside the source column -> nearest-layer (constant) value.
                  if (abs(pt - src_p(i, j, 1)) <= abs(pt - src_p(i, j, nsrc))) then
                     dst_data(i, j, k) = src_data(i, j, 1)
                  else
                     dst_data(i, j, k) = src_data(i, j, nsrc)
                  end if
               else
                  p1 = src_p(i, j, l)
                  p2 = src_p(i, j, l + 1)
                  if (p2 == p1) then
                     w = 0.0_fp
                  else
                     w = (pt - p1) / (p2 - p1)
                  end if
                  dst_data(i, j, k) = src_data(i, j, l) &
                     + w * (src_data(i, j, l + 1) - src_data(i, j, l))
               end if
            end do
         end do
      end do
   end subroutine vertical_interp_pressure

end module met_utilities_mod
```


