!> \file met_utilities_mod.F90
!! \brief Compatibility-preserving Fortran proxy delegating meteorological equations directly to modern C++
!!
module Met_Utilities_Mod
   use Precision_Mod, only: fp
   use Constants, only: Rd, Cp, g0
   use iso_c_binding, only: c_double, c_int, c_char, c_null_char

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
   public :: solar_zenith_angle

   ! C Interoperable Interface Definitions
   interface
      real(c_double) function catchem_met_potential_temperature(temp, press, sfc_press) bind(C, name="catchem_met_potential_temperature")
         import :: c_double
         real(c_double), value :: temp, press, sfc_press
      end function

      real(c_double) function catchem_met_virtual_temperature(temp, qv) bind(C, name="catchem_met_virtual_temperature")
         import :: c_double
         real(c_double), value :: temp, qv
      end function

      real(c_double) function catchem_met_dew_point(temp, rh) bind(C, name="catchem_met_dew_point")
         import :: c_double
         real(c_double), value :: temp, rh
      end function

      real(c_double) function catchem_met_relative_humidity(temp, qv, press) bind(C, name="catchem_met_relative_humidity")
         import :: c_double
         real(c_double), value :: temp, qv, press
      end function

      real(c_double) function catchem_met_saturation_vapor_pressure(temp) bind(C, name="catchem_met_saturation_vapor_pressure")
         import :: c_double
         real(c_double), value :: temp
      end function

      real(c_double) function catchem_met_mixing_ratio(q) bind(C, name="catchem_met_mixing_ratio")
         import :: c_double
         real(c_double), value :: q
      end function

      real(c_double) function catchem_met_specific_humidity(r) bind(C, name="catchem_met_specific_humidity")
         import :: c_double
         real(c_double), value :: r
      end function

      real(c_double) function catchem_met_dry_adiabatic_lapse_rate() bind(C, name="catchem_met_dry_adiabatic_lapse_rate")
         import :: c_double
      end function

      real(c_double) function catchem_met_bulk_richardson_number(t0, tz, u, z) bind(C, name="catchem_met_bulk_richardson_number")
         import :: c_double
         real(c_double), value :: t0, tz, u, z
      end function

      real(c_double) function catchem_met_monin_obukhov_length(ustar, t0, hflux, rho) bind(C, name="catchem_met_monin_obukhov_length")
         import :: c_double
         real(c_double), value :: ustar, t0, hflux, rho
      end function

      real(c_double) function catchem_met_friction_velocity(tau, rho) bind(C, name="catchem_met_friction_velocity")
         import :: c_double
         real(c_double), value :: tau, rho
      end function

      integer(c_int) function catchem_met_stability_classification(l) bind(C, name="catchem_met_stability_classification")
         import :: c_double, c_int
         real(c_double), value :: l
      end function

      real(c_double) function catchem_met_saturation_mixing_ratio(p, t) bind(C, name="catchem_met_saturation_mixing_ratio")
         import :: c_double
         real(c_double), value :: p, t
      end function

      real(c_double) function catchem_met_latent_heat_vaporization(t) bind(C, name="catchem_met_latent_heat_vaporization")
         import :: c_double
         real(c_double), value :: t
      end function

      real(c_double) function catchem_met_psychrometric_constant(p, lv) bind(C, name="catchem_met_psychrometric_constant")
         import :: c_double
         real(c_double), value :: p, lv
      end function

      real(c_double) function catchem_met_wind_profile_loglaw(ustar, z, z0) bind(C, name="catchem_met_wind_profile_loglaw")
         import :: c_double
         real(c_double), value :: ustar, z, z0
      end function

      real(c_double) function catchem_met_brunt_vaisala_frequency(t0, dtdz) bind(C, name="catchem_met_brunt_vaisala_frequency")
         import :: c_double
         real(c_double), value :: t0, dtdz
      end function

      real(c_double) function catchem_met_psi_m_businger(zeta) bind(C, name="catchem_met_psi_m_businger")
         import :: c_double
         real(c_double), value :: zeta
      end function

      real(c_double) function catchem_met_psi_h_businger(zeta) bind(C, name="catchem_met_psi_h_businger")
         import :: c_double
         real(c_double), value :: zeta
      end function

      real(c_double) function catchem_met_arrhenius_rate(a, ea, t) bind(C, name="catchem_met_arrhenius_rate")
         import :: c_double
         real(c_double), value :: a, ea, t
      end function

      real(c_double) function catchem_met_henrys_law_constant(h0, dh, t, t0) bind(C, name="catchem_met_henrys_law_constant")
         import :: c_double
         real(c_double), value :: h0, dh, t, t0
      end function

      real(c_double) function catchem_met_photolysis_rate_scaling(j0, sza) bind(C, name="catchem_met_photolysis_rate_scaling")
         import :: c_double
         real(c_double), value :: j0, sza
      end function

      real(c_double) function catchem_met_ppm_to_ugm3(ppm, m, t, p) bind(C, name="catchem_met_ppm_to_ugm3")
         import :: c_double
         real(c_double), value :: ppm, m, t, p
      end function

      real(c_double) function catchem_met_ugm3_to_ppm(ugm3, m, t, p) bind(C, name="catchem_met_ugm3_to_ppm")
         import :: c_double
         real(c_double), value :: ugm3, m, t, p
      end function

      real(c_double) function catchem_met_stokes_settling_velocity(dp, rho_p, rho_a, mu, cc) bind(C, name="catchem_met_stokes_settling_velocity")
         import :: c_double
         real(c_double), value :: dp, rho_p, rho_a, mu, cc
      end function

      real(c_double) function catchem_met_cunningham_correction_factor(dp, lambda) bind(C, name="catchem_met_cunningham_correction_factor")
         import :: c_double
         real(c_double), value :: dp, lambda
      end function

      real(c_double) function catchem_met_stokes_number(rho_p, d_p, u, mu, l) bind(C, name="catchem_met_stokes_number")
         import :: c_double
         real(c_double), value :: rho_p, d_p, u, mu, l
      end function

      real(c_double) function catchem_met_mean_free_path_air(temp, press) bind(C, name="catchem_met_mean_free_path_air")
         import :: c_double
         real(c_double), value :: temp, press
      end function

      real(c_double) function catchem_met_nuclear_decay(n0, lambda, t) bind(C, name="catchem_met_nuclear_decay")
         import :: c_double
         real(c_double), value :: n0, lambda, t
      end function

      subroutine catchem_met_solar_zenith_angle(jday, xhour, lat_rad, lon_rad, sza_deg, cossza) bind(C, name="catchem_met_solar_zenith_angle")
         import :: c_double, c_int
         integer(c_int), value :: jday
         real(c_double), value :: xhour, lat_rad, lon_rad
         real(c_double), intent(out) :: sza_deg, cossza
      end subroutine
   end interface

contains

   function potential_temperature(T, p, p0) result(theta)
      real(fp), intent(in) :: T, p, p0
      real(fp) :: theta
      theta = real(catchem_met_potential_temperature(real(T, c_double), real(p, c_double), real(p0, c_double)), fp)
   end function potential_temperature

   function virtual_temperature(T, qv) result(Tv)
      real(fp), intent(in) :: T, qv
      real(fp) :: Tv
      Tv = real(catchem_met_virtual_temperature(real(T, c_double), real(qv, c_double)), fp)
   end function virtual_temperature

   function dew_point(T, rh) result(Td)
      real(fp), intent(in) :: T, rh
      real(fp) :: Td
      Td = real(catchem_met_dew_point(real(T, c_double), real(rh, c_double)), fp)
   end function dew_point

   function relative_humidity(T, qv, p) result(rh)
      real(fp), intent(in) :: T, qv, p
      real(fp) :: rh
      rh = real(catchem_met_relative_humidity(real(T, c_double), real(qv, c_double), real(p, c_double)), fp)
   end function relative_humidity

   function saturation_vapor_pressure(T) result(es)
      real(fp), intent(in) :: T
      real(fp) :: es
      es = real(catchem_met_saturation_vapor_pressure(real(T, c_double)), fp)
   end function saturation_vapor_pressure

   function mixing_ratio(q) result(r)
      real(fp), intent(in) :: q
      real(fp) :: r
      r = real(catchem_met_mixing_ratio(real(q, c_double)), fp)
   end function mixing_ratio

   function specific_humidity(r) result(q)
      real(fp), intent(in) :: r
      real(fp) :: q
      q = real(catchem_met_specific_humidity(real(r, c_double)), fp)
   end function specific_humidity

   function dry_adiabatic_lapse_rate() result(gamma_d)
      real(fp) :: gamma_d
      gamma_d = real(catchem_met_dry_adiabatic_lapse_rate(), fp)
   end function dry_adiabatic_lapse_rate

   function bulk_richardson_number(T0, Tz, u, z) result(Ri)
      real(fp), intent(in) :: T0, Tz, u, z
      real(fp) :: Ri
      Ri = real(catchem_met_bulk_richardson_number(real(T0, c_double), real(Tz, c_double), real(u, c_double), real(z, c_double)), fp)
   end function bulk_richardson_number

   function monin_obukhov_length(ustar, T0, H, rho) result(L)
      real(fp), intent(in) :: ustar, T0, H, rho
      real(fp) :: L
      L = real(catchem_met_monin_obukhov_length(real(ustar, c_double), real(T0, c_double), real(H, c_double), real(rho, c_double)), fp)
   end function monin_obukhov_length

   function friction_velocity(tau, rho) result(ustar)
      real(fp), intent(in) :: tau, rho
      real(fp) :: ustar
      ustar = real(catchem_met_friction_velocity(real(tau, c_double), real(rho, c_double)), fp)
   end function friction_velocity

   function stability_classification(L) result(class)
      real(fp), intent(in) :: L
      integer :: class
      class = int(catchem_met_stability_classification(real(L, c_double)))
   end function stability_classification

   function saturation_mixing_ratio(p, T) result(ws)
      real(fp), intent(in) :: p, T
      real(fp) :: ws
      ws = real(catchem_met_saturation_mixing_ratio(real(p, c_double), real(T, c_double)), fp)
   end function saturation_mixing_ratio

   function latent_heat_vaporization(T) result(Lv)
      real(fp), intent(in) :: T
      real(fp) :: Lv
      Lv = real(catchem_met_latent_heat_vaporization(real(T, c_double)), fp)
   end function latent_heat_vaporization

   function psychrometric_constant(p, Lv) result(gamma)
      real(fp), intent(in) :: p, Lv
      real(fp) :: gamma
      gamma = real(catchem_met_psychrometric_constant(real(p, c_double), real(Lv, c_double)), fp)
   end function psychrometric_constant

   function wind_profile_loglaw(ustar, z, z0) result(u)
      real(fp), intent(in) :: ustar, z, z0
      real(fp) :: u
      u = real(catchem_met_wind_profile_loglaw(real(ustar, c_double), real(z, c_double), real(z0, c_double)), fp)
   end function wind_profile_loglaw

   function brunt_vaisala_frequency(T0, dTdz) result(N2)
      real(fp), intent(in) :: T0, dTdz
      real(fp) :: N2
      N2 = real(catchem_met_brunt_vaisala_frequency(real(T0, c_double), real(dTdz, c_double)), fp)
   end function brunt_vaisala_frequency

   function psi_m_businger(zeta) result(psi_m)
      real(fp), intent(in) :: zeta
      real(fp) :: psi_m
      psi_m = real(catchem_met_psi_m_businger(real(zeta, c_double)), fp)
   end function psi_m_businger

   function psi_h_businger(zeta) result(psi_h)
      real(fp), intent(in) :: zeta
      real(fp) :: psi_h
      psi_h = real(catchem_met_psi_h_businger(real(zeta, c_double)), fp)
   end function psi_h_businger

   function arrhenius_rate(A, Ea, T) result(k)
      real(fp), intent(in) :: A, Ea, T
      real(fp) :: k
      k = real(catchem_met_arrhenius_rate(real(A, c_double), real(Ea, c_double), real(T, c_double)), fp)
   end function arrhenius_rate

   function henrys_law_constant(H0, dH, T, T0) result(H)
      real(fp), intent(in) :: H0, dH, T, T0
      real(fp) :: H
      H = real(catchem_met_henrys_law_constant(real(H0, c_double), real(dH, c_double), real(T, c_double), real(T0, c_double)), fp)
   end function henrys_law_constant

   function photolysis_rate_scaling(J0, sza) result(J)
      real(fp), intent(in) :: J0, sza
      real(fp) :: J
      J = real(catchem_met_photolysis_rate_scaling(real(J0, c_double), real(sza, c_double)), fp)
   end function photolysis_rate_scaling

   function ppm_to_ugm3(ppm, M, T, p) result(ugm3)
      real(fp), intent(in) :: ppm, M, T, p
      real(fp) :: ugm3
      ugm3 = real(catchem_met_ppm_to_ugm3(real(ppm, c_double), real(M, c_double), real(T, c_double), real(p, c_double)), fp)
   end function ppm_to_ugm3

   function ugm3_to_ppm(ugm3, M, T, p) result(ppm)
      real(fp), intent(in) :: ugm3, M, T, p
      real(fp) :: ppm
      ppm = real(catchem_met_ugm3_to_ppm(real(ugm3, c_double), real(M, c_double), real(T, c_double), real(p, c_double)), fp)
   end function ugm3_to_ppm

   function stokes_settling_velocity(dp, rho_p, rho_a, mu, Cc) result(vs)
      real(fp), intent(in) :: dp, rho_p, rho_a, mu, Cc
      real(fp) :: vs
      vs = real(catchem_met_stokes_settling_velocity(real(dp, c_double), real(rho_p, c_double), real(rho_a, c_double), real(mu, c_double), real(Cc, c_double)), fp)
   end function stokes_settling_velocity

   function cunningham_correction_factor(dp, lambda) result(Cc)
      real(fp), intent(in) :: dp, lambda
      real(fp) :: Cc
      Cc = real(catchem_met_cunningham_correction_factor(real(dp, c_double), real(lambda, c_double)), fp)
   end function cunningham_correction_factor

   function stokes_number(rho_p, d_p, U, mu, L) result(Stk)
      real(fp), intent(in) :: rho_p, d_p, U, mu, L
      real(fp) :: Stk
      Stk = real(catchem_met_stokes_number(real(rho_p, c_double), real(d_p, c_double), real(U, c_double), real(mu, c_double), real(L, c_double)), fp)
   end function stokes_number

   function mean_free_path_air(T, p) result(lambda)
      real(fp), intent(in) :: T, p
      real(fp) :: lambda
      lambda = real(catchem_met_mean_free_path_air(real(T, c_double), real(p, c_double)), fp)
   end function mean_free_path_air

   function nuclear_decay(N0, lambda, t) result(N)
      real(fp), intent(in) :: N0, lambda, t
      real(fp) :: N
      N = real(catchem_met_nuclear_decay(real(N0, c_double), real(lambda, c_double), real(t, c_double)), fp)
   end function nuclear_decay

   subroutine solar_zenith_angle(jday, xhour, lat_rad, lon_rad, sza_deg, cossza)
      integer, intent(in) :: jday
      real(fp), intent(in) :: xhour, lat_rad, lon_rad
      real(fp), intent(out) :: sza_deg, cossza

      real(c_double) :: s_deg, c_sza

      call catchem_met_solar_zenith_angle(int(jday, c_int), real(xhour, c_double), real(lat_rad, c_double), real(lon_rad, c_double), s_deg, c_sza)
      sza_deg = real(s_deg, fp)
      cossza = real(c_sza, fp)
   end subroutine solar_zenith_angle

end module Met_Utilities_Mod
