#pragma once

#ifdef CATCHEM_ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
namespace catchem {
    namespace math {
        using Kokkos::acos;
        using Kokkos::cos;
        using Kokkos::exp;
        using Kokkos::log;
        using Kokkos::max;
        using Kokkos::min;
        using Kokkos::pow;
        using Kokkos::sin;
        using Kokkos::sqrt;
    } // namespace math
} // namespace catchem
#else
#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif
#include <algorithm>
#include <cmath>
namespace catchem {
    namespace math {
        using std::acos;
        using std::cos;
        using std::exp;
        using std::log;
        using std::max;
        using std::min;
        using std::pow;
        using std::sin;
        using std::sqrt;
    } // namespace math
} // namespace catchem
#endif

#include "catchem_constants.hpp"
#include "catchem_precision.hpp"

namespace catchem {
    namespace met_utilities {

        KOKKOS_INLINE_FUNCTION
        fp potential_temperature(fp temp, fp press, fp sfc_press) {
            return temp * math::pow(sfc_press / press, constants::RD / constants::CP);
        }

        KOKKOS_INLINE_FUNCTION
        fp virtual_temperature(fp temp, fp qv) {
            return temp * (1.0 + 0.61 * qv);
        }

        KOKKOS_INLINE_FUNCTION
        fp saturation_vapor_pressure(fp temp) {
            return 611.2 * math::exp(17.67 * (temp - 273.15) / (temp - 29.65));
        }

        KOKKOS_INLINE_FUNCTION
        fp dew_point(fp temp, fp rh) {
            fp es = saturation_vapor_pressure(temp);
            fp ed = rh * es;
            return 243.5 / (17.67 / math::log(ed / 611.2) - 1.0) + 273.15;
        }

        KOKKOS_INLINE_FUNCTION
        fp relative_humidity(fp T, fp qv, fp p) {
            fp e = qv * p / (0.622 + 0.378 * qv);
            fp es = saturation_vapor_pressure(T);
            fp rh = e / es;
            return math::max(static_cast<fp>(0.0), math::min(static_cast<fp>(1.0), rh));
        }

        KOKKOS_INLINE_FUNCTION
        fp mixing_ratio(fp q) {
            return q / (1.0 - q);
        }

        KOKKOS_INLINE_FUNCTION
        fp specific_humidity(fp r) {
            return r / (1.0 + r);
        }

        KOKKOS_INLINE_FUNCTION
        fp dry_adiabatic_lapse_rate() {
            return constants::G0 / constants::CP;
        }

        KOKKOS_INLINE_FUNCTION
        fp bulk_richardson_number(fp T0, fp Tz, fp u, fp z) {
            if (u > 0.0 && z > 0.0) {
                return (constants::G0 / T0) * (Tz - T0) * z / (u * u);
            }
            return 0.0;
        }

        KOKKOS_INLINE_FUNCTION
        fp monin_obukhov_length(fp ustar, fp T0, fp H, fp rho) {
            if (ustar > 0.0 && (H < 0.0 ? -H : H) > 0.0) {
                return -(ustar * ustar * ustar * rho * constants::CP * T0) / (0.41 * constants::G0 * H);
            }
            return 1.0e5; // Neutral/stable default
        }

        KOKKOS_INLINE_FUNCTION
        fp friction_velocity(fp tau, fp rho) {
            if (rho > 0.0) {
                fp abs_tau = tau < 0.0 ? -tau : tau;
                return math::sqrt(abs_tau / rho);
            }
            return 0.0;
        }

        KOKKOS_INLINE_FUNCTION
        int stability_classification(fp L) {
            if (L < -200.0)
                return -1;
            if (L > 200.0)
                return 1;
            return 0;
        }

        KOKKOS_INLINE_FUNCTION
        fp saturation_mixing_ratio(fp p, fp T) {
            fp es = saturation_vapor_pressure(T);
            return 0.622 * es / (p - es);
        }

        KOKKOS_INLINE_FUNCTION
        fp latent_heat_vaporization(fp T) {
            return 2.501e6 - 2.361e3 * (T - 273.15);
        }

        KOKKOS_INLINE_FUNCTION
        fp psychrometric_constant(fp p, fp Lv) {
            return constants::CP * p / (0.622 * Lv);
        }

        KOKKOS_INLINE_FUNCTION
        fp wind_profile_loglaw(fp ustar, fp z, fp z0) {
            if (z > z0 && z0 > 0.0) {
                return ustar / 0.41 * math::log(z / z0);
            }
            return 0.0;
        }

        KOKKOS_INLINE_FUNCTION
        fp brunt_vaisala_frequency(fp T0, fp dTdz) {
            return (constants::G0 / T0) * (dTdz + constants::G0 / constants::CP);
        }

        KOKKOS_INLINE_FUNCTION
        fp psi_m_businger(fp zeta) {
            if (zeta < 0.0) {
                return 2.0 * math::log((1.0 + math::sqrt(1.0 - 16.0 * zeta)) / 2.0);
            } else {
                return -5.0 * zeta;
            }
        }

        KOKKOS_INLINE_FUNCTION
        fp psi_h_businger(fp zeta) {
            if (zeta < 0.0) {
                return 2.0 * math::log((1.0 + math::sqrt(1.0 - 16.0 * zeta)) / 2.0);
            } else {
                return -5.0 * zeta;
            }
        }

        KOKKOS_INLINE_FUNCTION
        fp arrhenius_rate(fp A, fp Ea, fp T) {
            return A * math::exp(-Ea / (constants::RSTARG * T));
        }

        KOKKOS_INLINE_FUNCTION
        fp henrys_law_constant(fp H0, fp dH, fp T, fp T0) {
            return H0 * math::exp(-dH / constants::RSTARG * (1.0 / T - 1.0 / T0));
        }

        KOKKOS_INLINE_FUNCTION
        fp photolysis_rate_scaling(fp J0, fp sza) {
            return J0 * math::max(static_cast<fp>(0.0), math::cos(sza * constants::PI_180));
        }

        KOKKOS_INLINE_FUNCTION
        fp ppm_to_ugm3(fp ppm, fp M, fp T, fp p) {
            return ppm * 1.0e-6 * p * M / (constants::RSTARG * T) * 1.0e3;
        }

        KOKKOS_INLINE_FUNCTION
        fp ugm3_to_ppm(fp ugm3, fp M, fp T, fp p) {
            return ugm3 * (constants::RSTARG * T) / (p * M * 1.0e3) * 1.0e6;
        }

        KOKKOS_INLINE_FUNCTION
        fp stokes_settling_velocity(fp dp, fp rho_p, fp rho_a, fp mu, fp Cc) {
            return (dp * dp) * (rho_p - rho_a) * constants::G0 * Cc / (18.0 * mu);
        }

        KOKKOS_INLINE_FUNCTION
        fp cunningham_correction_factor(fp dp, fp lambda) {
            if (dp > 0.0 && lambda > 0.0) {
                return 1.0 + 2.0 * lambda / dp * (1.257 + 0.4 * math::exp(-1.1 * dp / lambda));
            }
            return 1.0;
        }

        KOKKOS_INLINE_FUNCTION
        fp stokes_number(fp rho_p, fp d_p, fp U, fp mu, fp L) {
            if (mu > 0.0 && L > 0.0) {
                return (rho_p * d_p * d_p * U) / (18.0 * mu * L);
            }
            return 0.0;
        }

        KOKKOS_INLINE_FUNCTION
        fp mean_free_path_air(fp T, fp p) {
            constexpr fp d_air = 3.7e-10; // Effective air molecule diameter [m]
            return constants::BOLTZ * T / (1.414213562373095 * constants::PI * d_air * d_air * p);
        }

        KOKKOS_INLINE_FUNCTION
        fp nuclear_decay(fp N0, fp lambda, fp t) {
            return N0 * math::exp(-lambda * t);
        }

        KOKKOS_INLINE_FUNCTION
        void solar_zenith_angle(int jday, fp xhour, fp lat_rad, fp lon_rad, fp& sza_deg, fp& cossza) {
            // Solar declination Fourier coefficients
            constexpr fp a0 = 0.006918;
            constexpr fp a1 = 0.399912;
            constexpr fp a2 = 0.006758;
            constexpr fp a3 = 0.002697;
            constexpr fp b1 = 0.070257;
            constexpr fp b2 = 0.000907;
            constexpr fp b3 = 0.000148;

            fp rad2deg = 180.0 / constants::PI;
            fp r = 2.0 * constants::PI * static_cast<fp>(jday - 1) / 365.0;

            fp dec = a0 - a1 * math::cos(r) + b1 * math::sin(r) - a2 * math::cos(2.0 * r) + b2 * math::sin(2.0 * r) -
                     a3 * math::cos(3.0 * r) + b3 * math::sin(3.0 * r);

            fp xlon = lon_rad * rad2deg;
            fp timloc = xhour + xlon / 15.0;
            while (timloc < 0.0)
                timloc += 24.0;
            while (timloc > 24.0)
                timloc -= 24.0;

            fp ahr = (timloc - 12.0) * 15.0;
            if (ahr < 0.0)
                ahr = -ahr;
            fp ahr_rad = ahr * constants::PI / 180.0;

            cossza = math::sin(lat_rad) * math::sin(dec) + math::cos(lat_rad) * math::cos(dec) * math::cos(ahr_rad);
            cossza = math::max(static_cast<fp>(-1.0), math::min(static_cast<fp>(1.0), cossza));

            sza_deg = math::acos(cossza) * rad2deg;

            if (cossza < 0.0)
                cossza = 0.0;
        }

    } // namespace met_utilities
} // namespace catchem
