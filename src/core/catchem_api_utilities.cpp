#include "catchem_api.hpp"
#include "catchem_api_internal.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_state_manager.hpp"
#include "catchem_unit_conversion.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using namespace catchem::api_internal;

extern "C" {
// =========================================================================
// TimeState C-Linkable API Implementation
// =========================================================================

void* catchem_time_state_create() {
    try {
        return static_cast<void*>(new catchem::TimeState());
    } catch (...) {
        return nullptr;
    }
}

void catchem_time_state_destroy(void* ptr) {
    try {
        delete static_cast<catchem::TimeState*>(ptr);
    } catch (...) {
    }
}

int catchem_time_state_init(void* ptr, int year, int month, int day, int hour, int minute, int second,
                            double timestep) {
    try {
        auto* ts = static_cast<catchem::TimeState*>(ptr);
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

int catchem_time_state_reset(void* ptr) {
    try {
        auto* ts = static_cast<catchem::TimeState*>(ptr);
        ts->year = 2000;
        ts->month = 1;
        ts->day = 1;
        ts->hour = 0;
        ts->minute = 0;
        ts->second = 0;
        ts->timestep = 3600.0;
        ts->calculate_derived_fields();
        return 0;
    } catch (...) {
        return -1;
    }
}

int catchem_time_state_get_year(void* ptr) {
    try {
        return static_cast<catchem::TimeState*>(ptr)->year;
    } catch (...) {
        return 0;
    }
}

int catchem_time_state_get_month(void* ptr) {
    try {
        return static_cast<catchem::TimeState*>(ptr)->month;
    } catch (...) {
        return 0;
    }
}

int catchem_time_state_get_day(void* ptr) {
    try {
        return static_cast<catchem::TimeState*>(ptr)->day;
    } catch (...) {
        return 0;
    }
}

int catchem_time_state_get_hour(void* ptr) {
    try {
        return static_cast<catchem::TimeState*>(ptr)->hour;
    } catch (...) {
        return 0;
    }
}

int catchem_time_state_get_minute(void* ptr) {
    try {
        return static_cast<catchem::TimeState*>(ptr)->minute;
    } catch (...) {
        return 0;
    }
}

int catchem_time_state_get_second(void* ptr) {
    try {
        return static_cast<catchem::TimeState*>(ptr)->second;
    } catch (...) {
        return 0;
    }
}

double catchem_time_state_get_timestep(void* ptr) {
    try {
        return static_cast<catchem::TimeState*>(ptr)->timestep;
    } catch (...) {
        return 0.0;
    }
}

double catchem_time_state_get_julian_date(void* ptr) {
    try {
        return static_cast<catchem::TimeState*>(ptr)->julian_date;
    } catch (...) {
        return 0.0;
    }
}

int catchem_time_state_get_doy(void* ptr) {
    try {
        return static_cast<catchem::TimeState*>(ptr)->doy;
    } catch (...) {
        return 0;
    }
}

double catchem_time_state_get_cos_sza(void* ptr, double lat, double lon, bool mid_timestep) {
    try {
        return static_cast<catchem::TimeState*>(ptr)->get_cos_sza(lat, lon, mid_timestep);
    } catch (...) {
        return 0.0;
    }
}

int catchem_time_state_get_timezone_offset([[maybe_unused]] void* ptr, double lon) {
    try {
        int offset = static_cast<int>(lon / 15.0);
        return std::max(-12, std::min(14, offset));
    } catch (...) {
        return 0;
    }
}

bool catchem_time_state_is_leap_year(int year) {
    return catchem::TimeState::is_leap_year(year);
}

int catchem_time_state_get_days_in_month(int month, int year) {
    return catchem::TimeState::get_days_in_month(month, year);
}

bool catchem_time_state_is_global_holiday(int month, int day) {
    return (month == 1 && day == 1) || (month == 12 && day == 25);
}

bool catchem_time_state_is_us_holiday(int month, int day) {
    return (month == 7 && day == 4) || (month == 11 && day >= 22 && day <= 28);
}

// =========================================================================
// UnitConversion C-Linkable API Implementation
// =========================================================================

double catchem_convert_concentration(double val, const char* from_units, const char* to_units, double mw, double temp,
                                     double press, int* rc) {
    try {
        *rc = 0;
        std::string from = catchem::unit_conversion::to_upper(from_units);
        std::string to = catchem::unit_conversion::to_upper(to_units);
        if (from == to) {
            return val;
        }

        // Identify and extract VMR scale factors relative to PPBV
        auto get_vmr_factor = [](const std::string& unit, bool& is_vmr) -> double {
            is_vmr = true;
            if (unit == "PPMV" || unit == "PPM")
                return 1e3;
            if (unit == "PPBV" || unit == "PPB")
                return 1.0;
            if (unit == "PPTV" || unit == "PPT")
                return 1e-3;
            is_vmr = false;
            return 1.0;
        };

        bool from_is_vmr = false;
        double from_factor = get_vmr_factor(from, from_is_vmr);
        bool to_is_vmr = false;
        double to_factor = get_vmr_factor(to, to_is_vmr);

        // Direct VMR-to-VMR conversion
        if (from_is_vmr && to_is_vmr) {
            return val * (from_factor / to_factor);
        }

        // Normalize VMR conversion from/to mass/volume units
        std::string from_normalized = from_is_vmr ? "PPBV" : from;
        std::string to_normalized = to_is_vmr ? "PPBV" : to;
        double input_val = from_is_vmr ? (val * from_factor) : val;

        double result = 0.0;
        std::string key = from_normalized + " -> " + to_normalized;

        if (key == "PPBV -> UG/M3" || key == "PPBV -> UG M-3" || key == "PPBV -> UG/M^3") {
            result = catchem::unit_conversion::ppbv_to_ugm3(input_val, mw, temp, press);
        } else if (key == "UG/M3 -> PPBV" || key == "UG M-3 -> PPBV" || key == "UG/M^3 -> PPBV") {
            result = catchem::unit_conversion::ugm3_to_ppbv(input_val, mw, temp, press);
        } else if (key == "PPBV -> MG/M3" || key == "PPBV -> MG M-3" || key == "PPBV -> MG/M^3") {
            result = catchem::unit_conversion::ppbv_to_ugm3(input_val, mw, temp, press) * 1e-3;
        } else if (key == "MG/M3 -> PPBV" || key == "MG M-3 -> PPBV" || key == "MG/M^3 -> PPBV") {
            result = catchem::unit_conversion::ugm3_to_ppbv(input_val * 1e3, mw, temp, press);
        } else if (key == "MOLEC/CM3 -> PPBV" || key == "MOLEC CM-3 -> PPBV" || key == "MOLEC/CM^3 -> PPBV") {
            result = catchem::unit_conversion::molcm3_to_ppbv(input_val, temp, press);
        } else if (key == "PPBV -> MOLEC/CM3" || key == "PPBV -> MOLEC CM-3" || key == "PPBV -> MOLEC/CM^3") {
            result = catchem::unit_conversion::ppbv_to_molcm3(input_val, temp, press);
        } else if (key == "MOLEC/CM3 -> UG/M3" || key == "MOLEC CM-3 -> UG/M3" || key == "MOLEC/CM^3 -> UG/M3" ||
                   key == "MOLEC/CM3 -> UG M-3" || key == "MOLEC/CM3 -> UG/M^3") {
            double ppbv = catchem::unit_conversion::molcm3_to_ppbv(input_val, temp, press);
            result = catchem::unit_conversion::ppbv_to_ugm3(ppbv, mw, temp, press);
        } else {
            *rc = -1;
            return val;
        }

        // If target is a VMR, scale output from PPBV to target
        if (to_is_vmr) {
            result /= to_factor;
        }

        return result;
    } catch (...) {
        *rc = -1;
        return val;
    }
}

double catchem_convert_pressure(double val, const char* from_units, const char* to_units, int* rc) {
    try {
        *rc = 0;
        if (catchem::unit_conversion::to_upper(from_units) == catchem::unit_conversion::to_upper(to_units)) {
            return val;
        }
        return catchem::unit_conversion::convert_pressure(val, from_units, to_units, *rc);
    } catch (...) {
        *rc = -1;
        return val;
    }
}

double catchem_convert_temperature(double val, const char* from_units, const char* to_units, int* rc) {
    try {
        *rc = 0;
        if (catchem::unit_conversion::to_upper(from_units) == catchem::unit_conversion::to_upper(to_units)) {
            return val;
        }
        return catchem::unit_conversion::convert_temperature(val, from_units, to_units, *rc);
    } catch (...) {
        *rc = -1;
        return val;
    }
}

double catchem_convert_flux(double val, const char* from_units, const char* to_units, double mw, int* rc) {
    try {
        *rc = 0;
        if (catchem::unit_conversion::to_upper(from_units) == catchem::unit_conversion::to_upper(to_units)) {
            return val;
        }
        return catchem::unit_conversion::convert_flux(val, from_units, to_units, mw, *rc);
    } catch (...) {
        *rc = -1;
        return val;
    }
}

double catchem_convert_rate_constant(double val, const char* from_units, const char* to_units, int* rc) {
    try {
        *rc = 0;
        if (catchem::unit_conversion::to_upper(from_units) == catchem::unit_conversion::to_upper(to_units)) {
            return val;
        }
        return catchem::unit_conversion::convert_rate_constant(val, from_units, to_units, *rc);
    } catch (...) {
        *rc = -1;
        return val;
    }
}

double catchem_convert_mass_units(double val, const char* from_units, const char* to_units, int* rc) {
    try {
        *rc = 0;
        if (catchem::unit_conversion::to_upper(from_units) == catchem::unit_conversion::to_upper(to_units)) {
            return val;
        }
        return catchem::unit_conversion::convert_mass_units(val, from_units, to_units, *rc);
    } catch (...) {
        *rc = -1;
        return val;
    }
}

double catchem_calculate_air_density(double temp, double press, double humidity, bool use_humidity) {
    try {
        return catchem::unit_conversion::calculate_air_density(temp, press, humidity, use_humidity);
    } catch (...) {
        return 0.0;
    }
}

double catchem_calculate_molecular_weight(const char* formula) {
    try {
        return catchem::unit_conversion::calculate_molecular_weight(formula);
    } catch (...) {
        return 0.0;
    }
}

double catchem_convert_imperial(double val, const char* from_units, const char* to_units, const char* category,
                                int* rc) {
    try {
        return catchem::unit_conversion::convert_imperial(val, from_units, to_units, category, *rc);
    } catch (...) {
        *rc = -1;
        return val;
    }
}

int catchem_convert_process_concentration_units(catchem::fp* values, int size, const char* from_units,
                                                const char* to_units, catchem::fp mw, catchem::fp temp,
                                                catchem::fp press) {
    try {
        int rc = 0;
        for (int i = 0; i < size; ++i) {
            values[i] = static_cast<catchem::fp>(catchem_convert_concentration(
                static_cast<double>(values[i]), from_units, to_units, static_cast<double>(mw),
                static_cast<double>(temp), static_cast<double>(press), &rc));
            if (rc != 0)
                return rc;
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

int catchem_convert_process_flux_units(catchem::fp* values, int size, const char* from_units, const char* to_units,
                                       catchem::fp mw) {
    try {
        int rc = 0;
        for (int i = 0; i < size; ++i) {
            values[i] = static_cast<catchem::fp>(catchem_convert_flux(static_cast<double>(values[i]), from_units,
                                                                      to_units, static_cast<double>(mw), &rc));
            if (rc != 0)
                return rc;
        }
        return 0;
    } catch (...) {
        return -1;
    }
}
// =========================================================================
// Meteorological Core Calculation C-API Definitions
// =========================================================================
double catchem_met_potential_temperature(double temp, double press, double sfc_press) {
    return catchem::met_utilities::potential_temperature(temp, press, sfc_press);
}

double catchem_met_virtual_temperature(double temp, double qv) {
    return catchem::met_utilities::virtual_temperature(temp, qv);
}

double catchem_met_dew_point(double temp, double rh) {
    return catchem::met_utilities::dew_point(temp, rh);
}

double catchem_met_relative_humidity(double temp, double qv, double press) {
    return catchem::met_utilities::relative_humidity(temp, qv, press);
}

double catchem_met_saturation_vapor_pressure(double temp) {
    return catchem::met_utilities::saturation_vapor_pressure(temp);
}

double catchem_met_monin_obukhov_length(double ustar, double t0, double hflux, double rho) {
    return catchem::met_utilities::monin_obukhov_length(ustar, t0, hflux, rho);
}

double catchem_met_friction_velocity(double tau, double rho) {
    return catchem::met_utilities::friction_velocity(tau, rho);
}

double catchem_met_cunningham_correction_factor(double dp, double lambda) {
    return catchem::met_utilities::cunningham_correction_factor(dp, lambda);
}

double catchem_met_mean_free_path_air(double temp, double press) {
    return catchem::met_utilities::mean_free_path_air(temp, press);
}

void catchem_met_solar_zenith_angle(int doy, double hour, double lat_rad, double lon_rad, double* sza_deg,
                                    double* cossza) {
    try {
        catchem::fp sza_tmp = 0.0;
        catchem::fp cos_tmp = 0.0;
        catchem::met_utilities::solar_zenith_angle(doy, static_cast<catchem::fp>(hour),
                                                   static_cast<catchem::fp>(lat_rad), static_cast<catchem::fp>(lon_rad),
                                                   sza_tmp, cos_tmp);
        *sza_deg = static_cast<double>(sza_tmp);
        *cossza = static_cast<double>(cos_tmp);
    } catch (...) {
        *sza_deg = 0.0;
        *cossza = 0.0;
    }
}

double catchem_met_mixing_ratio(double q) {
    return catchem::met_utilities::mixing_ratio(q);
}

double catchem_met_specific_humidity(double r) {
    return catchem::met_utilities::specific_humidity(r);
}

double catchem_met_dry_adiabatic_lapse_rate() {
    return catchem::met_utilities::dry_adiabatic_lapse_rate();
}

double catchem_met_bulk_richardson_number(double t0, double tz, double u, double z) {
    return catchem::met_utilities::bulk_richardson_number(t0, tz, u, z);
}

int catchem_met_stability_classification(double l) {
    return catchem::met_utilities::stability_classification(l);
}

double catchem_met_saturation_mixing_ratio(double p, double t) {
    return catchem::met_utilities::saturation_mixing_ratio(p, t);
}

double catchem_met_latent_heat_vaporization(double t) {
    return catchem::met_utilities::latent_heat_vaporization(t);
}

double catchem_met_psychrometric_constant(double p, double lv) {
    return catchem::met_utilities::psychrometric_constant(p, lv);
}

double catchem_met_wind_profile_loglaw(double ustar, double z, double z0) {
    return catchem::met_utilities::wind_profile_loglaw(ustar, z, z0);
}

double catchem_met_brunt_vaisala_frequency(double t0, double dtdz) {
    return catchem::met_utilities::brunt_vaisala_frequency(t0, dtdz);
}

double catchem_met_psi_m_businger(double zeta) {
    return catchem::met_utilities::psi_m_businger(zeta);
}

double catchem_met_psi_h_businger(double zeta) {
    return catchem::met_utilities::psi_h_businger(zeta);
}

double catchem_met_arrhenius_rate(double a, double ea, double t) {
    return catchem::met_utilities::arrhenius_rate(a, ea, t);
}

double catchem_met_henrys_law_constant(double h0, double dh, double t, double t0) {
    return catchem::met_utilities::henrys_law_constant(h0, dh, t, t0);
}

double catchem_met_photolysis_rate_scaling(double j0, double sza) {
    return catchem::met_utilities::photolysis_rate_scaling(j0, sza);
}

double catchem_met_ppm_to_ugm3(double ppm, double m, double t, double p) {
    return catchem::met_utilities::ppm_to_ugm3(ppm, m, t, p);
}

double catchem_met_ugm3_to_ppm(double ugm3, double m, double t, double p) {
    return catchem::met_utilities::ugm3_to_ppm(ugm3, m, t, p);
}

double catchem_met_stokes_settling_velocity(double dp, double rho_p, double rho_a, double mu, double cc) {
    return catchem::met_utilities::stokes_settling_velocity(dp, rho_p, rho_a, mu, cc);
}

double catchem_met_stokes_number(double rho_p, double d_p, double u, double mu, double l) {
    return catchem::met_utilities::stokes_number(rho_p, d_p, u, mu, l);
}

double catchem_met_nuclear_decay(double n0, double lambda, double t) {
    return catchem::met_utilities::nuclear_decay(n0, lambda, t);
}
}
