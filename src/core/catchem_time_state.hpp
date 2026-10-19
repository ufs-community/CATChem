#pragma once
#include "catchem_kokkos_compat.hpp"
#include <algorithm>
#include <cmath>

#include "catchem_constants.hpp"
#include "catchem_precision.hpp"

namespace catchem {

    struct TimeState {
        int year = 2000;
        int month = 1;
        int day = 1;
        int hour = 0;
        int minute = 0;
        int second = 0;
        double timestep = 3600.0;
        double julian_date = 0.0;
        int doy = 1;

        KOKKOS_FUNCTION
        static bool is_leap_year(int y) { return (y % 4 == 0 && y % 100 != 0) || (y % 400 == 0); }

        KOKKOS_FUNCTION
        static int get_days_in_month(int m, int y) {
            const int days_per_month[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
            if (m < 1 || m > 12)
                return 0;
            if (m == 2 && is_leap_year(y))
                return 29;
            return days_per_month[m];
        }

        KOKKOS_FUNCTION
        void calculate_derived_fields() {
            // Julian Day arithmetic
            int a = (14 - month) / 12;
            int y = year + 4800 - a;
            int m = month + 12 * a - 3;
            int jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;

            julian_date = static_cast<double>(jdn) - 0.5 + (hour + minute / 60.0 + second / 3600.0) / 24.0;

            // Calculate Day of Year (DOY)
            const int days_per_month[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
            doy = day;
            for (int i = 1; i < month; ++i) {
                doy += days_per_month[i];
                if (i == 2 && is_leap_year(year)) {
                    doy += 1;
                }
            }
        }

        KOKKOS_FUNCTION
        void advance(double dt) {
            timestep = dt;
            int total_sec = second + static_cast<int>(dt);
            second = total_sec % 60;
            int total_min = minute + total_sec / 60;
            minute = total_min % 60;
            int total_hr = hour + total_min / 60;
            hour = total_hr % 24;
            int extra_days = total_hr / 24;

            while (extra_days > 0) {
                int dim = get_days_in_month(month, year);
                if (day + extra_days <= dim) {
                    day += extra_days;
                    extra_days = 0;
                } else {
                    extra_days -= (dim - day + 1);
                    day = 1;
                    if (month == 12) {
                        month = 1;
                        year += 1;
                    } else {
                        month += 1;
                    }
                }
            }
            calculate_derived_fields();
        }

        KOKKOS_FUNCTION
        double get_cos_sza(double lat_deg, double lon_deg, bool mid_timestep = false) const {
            double lat_rad = lat_deg * constants::PI_180;

            double frac_hour = hour + minute / 60.0 + second / 3600.0;
            if (mid_timestep) {
                frac_hour += (timestep / 2.0) / 3600.0;
            }

            // Day angle [radians]
            double gamma = 2.0 * constants::PI * (doy - 1.0) / 365.0;

            // Solar declination (high-precision Fourier calculation matches GOCART2G)
            double dec = 0.006918 - 0.399912 * std::cos(gamma) + 0.070257 * std::sin(gamma) -
                         0.006758 * std::cos(2.0 * gamma) + 0.000907 * std::sin(2.0 * gamma) -
                         0.002697 * std::cos(3.0 * gamma) + 0.001480 * std::sin(3.0 * gamma);

            // Equation of time
            double eqtime = 229.18 * (0.000075 + 0.001868 * std::cos(gamma) - 0.032077 * std::sin(gamma) -
                                      0.014615 * std::cos(2.0 * gamma) - 0.040849 * std::sin(2.0 * gamma));

            double time_offset = eqtime + 4.0 * lon_deg;
            double true_solar_time = frac_hour * 60.0 + time_offset;

            double hour_angle = (true_solar_time / 4.0) - 180.0;
            double ha_rad = hour_angle * constants::PI_180;

            double cos_sza = std::sin(lat_rad) * std::sin(dec) + std::cos(lat_rad) * std::cos(dec) * std::cos(ha_rad);

            // Clamp output safely to [-1.0, 1.0]
            return std::max(-1.0, std::min(1.0, cos_sza));
        }
    };

} // namespace catchem
