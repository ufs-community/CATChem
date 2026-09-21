# Specification: High-Performance C++ Core Orchestration and Legacy Fortran Elimination

* **Status:** Approved
* **Authors:** Gemini CLI Architect
* **Created:** July 9, 2026
* **Target Version:** 2.0.0
* **Pillars:** High-Performance C++20, Kokkos parallelism, ISO_C_BINDING Memory Interoperability

## 1. Executive Summary & Goals

The CATChem core is undergoing modernization to transition all memory management, physical calculations, configuration loading, and process scheduling from a legacy Fortran structure to a modern, high-performance C++20 standard using Kokkos.

To maintain binary and source-level compatibility with FV3 host models and the ESMF/NUOPC Cap drivers, we are implementing a **C++ Delegate Wrapper Pattern**. This pattern strips the remaining Fortran modules of physical calculations, state storage, and custom allocations, converting them into zero-overhead compatibility proxies.

### Success Criteria
* **0% Business Logic in Fortran:** No physical constants, numerical derivations, or calendar calculations remain in Fortran.
* **100% Core C++ Ownership:** Memory, time-tracking, and unit-conversions live entirely in vectorized, template-optimized C++ structures.
* **100% Compatibility:** Legacy drivers and Caps continue to compile and link flawlessly without modifying a single line of their downstream code.
* **100% Thread-Safety & Exception Shields:** Language boundary transitions are robust, thread-safe, and compile warning-free in GCC/Spack.

---

## 2. Component Design & Dynamic Mappings

### 2.1. TimeState Modernization
`TimeState_Mod.F90` becomes a thin wrapper delegating all date, timezone, leap-year, and Julian Date calculations to C++ (`catchem_time_state.hpp`).

#### Modernized C++ Structure (`src/core/catchem_time_state.hpp`)
```cpp
#pragma once
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
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

        void calculate_derived_fields() {
            int a = (14 - month) / 12;
            int y = year + 4800 - a;
            int m = month + 12 * a - 3;
            int jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;

            julian_date = static_cast<double>(jdn) - 0.5 +
                          (hour + minute / 60.0 + second / 3600.0) / 24.0;

            const int days_per_month[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
            doy = day;
            for (int i = 1; i < month; ++i) {
                doy += days_per_month[i];
                if (i == 2 && is_leap_year(year)) {
                    doy += 1;
                }
            }
        }

        static bool is_leap_year(int y) {
            return (y % 4 == 0 && y % 100 != 0) || (y % 400 == 0);
        }

        static int get_days_in_month(int m, int y) {
            const int days_per_month[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
            if (m < 1 || m > 12) return 0;
            if (m == 2 && is_leap_year(y)) return 29;
            return days_per_month[m];
        }

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

        double get_cos_sza(double lat_deg, double lon_deg, bool mid_timestep = false) const {
            double lat_rad = lat_deg * constants::PI_180;
            double frac_hour = hour + minute / 60.0 + second / 3600.0;
            if (mid_timestep) {
                frac_hour += (timestep / 2.0) / 3600.0;
            }

            double gamma = 2.0 * constants::PI * (doy - 1.0) / 365.0;
            double dec = 0.006918 - 0.399912 * std::cos(gamma) + 0.070257 * std::sin(gamma) -
                         0.006758 * std::cos(2.0 * gamma) + 0.000907 * std::sin(2.0 * gamma) -
                         0.002697 * std::cos(3.0 * gamma) + 0.001480 * std::sin(3.0 * gamma);

            double eqtime = 229.18 * (0.000075 + 0.001868 * std::cos(gamma) - 0.032077 * std::sin(gamma) -
                                      0.014615 * std::cos(2.0 * gamma) - 0.040849 * std::sin(2.0 * gamma));

            double time_offset = eqtime + 4.0 * lon_deg;
            double true_solar_time = frac_hour * 60.0 + time_offset;
            double cos_sza = std::sin(lat_rad) * std::sin(dec) +
                             std::cos(lat_rad) * std::cos(dec) * std::cos(((true_solar_time / 4.0) - 180.0) * constants::PI_180);
            return std::max(-1.0, std::min(1.0, cos_sza));
        }
    };

} // namespace catchem
```

### 2.2. UnitConversion Modernization
All molecular weight evaluation, meteorological density calculations, and chemical unit conversions are transferred from `UnitConversion_Mod.F90` into a parallel-ready template header `catchem_unit_conversion.hpp`.

---

## 3. Flat C API & Exception Shield Interfaces

All BIND(C) export points are added to `src/core/catchem_api.hpp` and `catchem_api.cpp` with strict exception shields to ensure no native C++ exception ever escapes into Fortran.

```cpp
extern "C" {
    // =========================================================================
    // TimeState Core Exports
    // =========================================================================
    void* catchem_time_state_create() {
        try {
            return new catchem::TimeState();
        } catch (...) {
            return nullptr;
        }
    }

    void catchem_time_state_destroy(void* ptr) {
        try {
            delete static_cast<catchem::TimeState*>(ptr);
        } catch (...) {}
    }

    int catchem_time_state_init(void* ptr, int year, int month, int day, int hour, int minute, int second, double timestep) {
        try {
            auto ts = static_cast<catchem::TimeState*>(ptr);
            ts->year = year;
            ts->month = month;
            ts->day = day;
            ts->hour = hour;
            ts->minute = minute;
            ts->second = second;
            ts->timestep = timestep;
            ts->calculate_derived_fields();
            return 0;
        } catch (...) {
            return -1;
        }
    }

    int catchem_time_state_advance(void* ptr, double dt) {
        try {
            static_cast<catchem::TimeState*>(ptr)->advance(dt);
            return 0;
        } catch (...) {
            return -1;
        }
    }
}
```

---

## 4. Testing & Validation Strategy

The complete test suite of 9 targets will be compiled and executed inside the Spack-based Docker container. We will explicitly test:
1. **TimeState Unit Test (`test_TimeState`):** Verifies that advancing, leap year handling, and Julian day arithmetic returned from the C++ delegate is numerically identical to the legacy implementation.
2. **UnitConversion Unit Test (`test_UnitConversion`):** Verifies that molecular weights, ppmv-to-ugm3, mass conversions, and array-vectorized conversions from C++ are identical.
3. **MetState and StateManager Unit Tests:** Validates that meteorological attributes are correctly queried on the unified C++ core.
