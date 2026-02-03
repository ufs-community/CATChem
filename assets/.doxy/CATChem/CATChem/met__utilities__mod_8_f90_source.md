

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

end module met_utilities_mod
```


