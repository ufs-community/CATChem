#pragma once

#ifdef CATCHEM_ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#else
#define KOKKOS_INLINE_FUNCTION inline
#endif

#include "catchem_constants.hpp"
#include "catchem_precision.hpp"
#include <algorithm>
#include <string>
#include <vector>

namespace catchem {
    namespace unit_conversion {

        // Numerical constants matching legacy constants.F90
        constexpr double BOLTZ = 1.380649e-23;     // Boltzmann constant [J/K]
        constexpr double AVOGADRO = 6.02214076e23; // Avogadro number [molecules/mol]
        constexpr double RSTARG = 8.314462618;     // Universal gas constant [J/(mol*K)]
        constexpr double AIRMW = 28.9644;          // Dry air molecular weight [g/mol]

        KOKKOS_INLINE_FUNCTION
        double calculate_air_density(double temp, double press, double humidity = 0.0, bool use_humidity = false) {
            if (use_humidity) {
                double rv = 461.5;
                double rd = 287.05;
                double eps = rd / rv;
                double vtemp = temp * (1.0 + humidity * (1.0 - eps) / eps);
                return press / (rd * vtemp);
            } else {
                return press * AIRMW / (RSTARG * temp) * 1.0e-3;
            }
        }

        KOKKOS_INLINE_FUNCTION
        double ppbv_to_ugm3(double ppbv, double mw, double temp, double press) {
            return ppbv * mw * press / (RSTARG * temp) * 1.0e-3;
        }

        KOKKOS_INLINE_FUNCTION
        double ugm3_to_ppbv(double ugm3, double mw, double temp, double press) {
            return ugm3 * RSTARG * temp / (mw * press) * 1.0e3;
        }

        KOKKOS_INLINE_FUNCTION
        double molcm3_to_ppbv(double molcm3, double temp, double press) {
            double number_density = press / (BOLTZ * temp) * 1.0e-6;
            return molcm3 / number_density * 1.0e9;
        }

        KOKKOS_INLINE_FUNCTION
        double ppbv_to_molcm3(double ppbv, double temp, double press) {
            double number_density = press / (BOLTZ * temp) * 1.0e-6;
            return ppbv * number_density * 1.0e-9;
        }

        KOKKOS_INLINE_FUNCTION
        double ppmv_to_mgm3(double ppmv, double mw, double temp, double press) {
            return ppmv * mw * press / (RSTARG * temp);
        }

        KOKKOS_INLINE_FUNCTION
        double mgm3_to_ppmv(double mgm3, double mw, double temp, double press) {
            return mgm3 * RSTARG * temp / (mw * press);
        }

        KOKKOS_INLINE_FUNCTION
        double molcm2s_to_kgm2s(double molcm2s, double mw) {
            return (molcm2s / AVOGADRO) * mw * 1.0e-3 * 1.0e4;
        }

        KOKKOS_INLINE_FUNCTION
        double kgm2s_to_molcm2s(double kgm2s, double mw) {
            return (kgm2s / mw) * AVOGADRO * 1.0e3 * 1.0e-4;
        }

        inline double calculate_molecular_weight(const std::string& formula) {
            std::string f = formula;
            std::transform(f.begin(), f.end(), f.begin(), ::toupper);
            if (f == "O3")
                return 48.0;
            if (f == "NO2")
                return 46.0055;
            if (f == "CO")
                return 28.01;
            if (f == "SO2")
                return 64.066;
            if (f == "HNO3")
                return 63.012;
            if (f == "NH3")
                return 17.031;
            if (f == "HCHO")
                return 30.026;
            if (f == "CH4")
                return 16.04;
            if (f == "CO2")
                return 44.01;
            return AIRMW;
        }

        inline std::string to_upper(const std::string& str) {
            std::string s = str;
            std::transform(s.begin(), s.end(), s.begin(), ::toupper);
            return s;
        }

        inline double convert_pressure(double val, const std::string& from_u, const std::string& to_u, int& rc) {
            rc = 0;
            double val_pa = val;
            std::string from = to_upper(from_u);
            std::string to = to_upper(to_u);

            // Convert to Pa
            if (from == "PA")
                val_pa = val;
            else if (from == "HPA" || from == "MB" || from == "MBAR")
                val_pa = val * 100.0;
            else if (from == "ATM")
                val_pa = val * 101325.0;
            else if (from == "TORR" || from == "MMHG")
                val_pa = val * 133.322387415;
            else if (from == "PSI")
                val_pa = val * 6894.75729317;
            else {
                rc = -1;
                return val;
            }

            // Convert from Pa to Target
            if (to == "PA")
                return val_pa;
            if (to == "HPA" || to == "MB" || to == "MBAR")
                return val_pa / 100.0;
            if (to == "ATM")
                return val_pa / 101325.0;
            if (to == "TORR" || to == "MMHG")
                return val_pa / 133.322387415;
            if (to == "PSI")
                return val_pa / 6894.75729317;

            rc = -1;
            return val;
        }

        inline double convert_temperature(double val, const std::string& from_u, const std::string& to_u, int& rc) {
            rc = 0;
            double val_k = val;
            std::string from = to_upper(from_u);
            std::string to = to_upper(to_u);

            // Convert to K
            if (from == "K" || from == "KELVIN")
                val_k = val;
            else if (from == "C" || from == "CELSIUS" || from == "DEGC" || from == "DEGREE_CELSIUS" ||
                     from == "DEGREE_C")
                val_k = val + 271.15; // match legacy constant
            else if (from == "F" || from == "FAHRENHEIT" || from == "DEGF" || from == "DEGREE_FAHRENHEIT" ||
                     from == "DEGREE_F")
                val_k = (val - 32.0) * 5.0 / 9.0 + 273.15;
            else {
                rc = -1;
                return val;
            }

            // Convert from K to Target
            if (to == "K" || to == "KELVIN")
                return val_k;
            if (to == "C" || to == "CELSIUS" || to == "DEGC" || to == "DEGREE_CELSIUS" || to == "DEGREE_C")
                return val_k - 271.15;
            if (to == "F" || to == "FAHRENHEIT" || to == "DEGF" || to == "DEGREE_FAHRENHEIT" || to == "DEGREE_F")
                return (val_k - 273.15) * 9.0 / 5.0 + 32.0;

            rc = -1;
            return val;
        }

        inline double convert_flux(double val, const std::string& from_u, const std::string& to_u, double mw, int& rc) {
            rc = 0;
            double val_kgm2s = val;
            std::string from = to_upper(from_u);
            std::string to = to_upper(to_u);

            // Convert to kg/m2/s
            if (from == "KG/M2/S" || from == "KG M-2 S-1")
                val_kgm2s = val;
            else if (from == "MOLEC/CM2/S" || from == "MOLEC CM-2 S-1" || from == "MOLECULES/CM2/S")
                val_kgm2s = molcm2s_to_kgm2s(val, mw);
            else if (from == "MOLEC/M2/S" || from == "MOLEC M-2 S-1" || from == "MOLECULES/M2/S")
                val_kgm2s = molcm2s_to_kgm2s(val * 1e-4, mw);
            else {
                rc = -1;
                return val;
            }

            // Convert from kg/m2/s to Target
            if (to == "KG/M2/S" || to == "KG M-2 S-1")
                return val_kgm2s;
            if (to == "MOLEC/CM2/S" || to == "MOLEC CM-2 S-1" || to == "MOLECULES/CM2/S")
                return kgm2s_to_molcm2s(val_kgm2s, mw);
            if (to == "MOLEC/M2/S" || to == "MOLEC M-2 S-1" || to == "MOLECULES/M2/S")
                return kgm2s_to_molcm2s(val_kgm2s, mw) * 1e4;

            rc = -1;
            return val;
        }

        inline double convert_rate_constant(double val, const std::string& from_u, const std::string& to_u, int& rc) {
            rc = 0;
            double val_mks = val;
            std::string from = to_upper(from_u);
            std::string to = to_upper(to_u);

            if (from == "CM3/MOLECULE/S" || from == "CM3 MOLEC-1 S-1")
                val_mks = val * 1e-6;
            else if (from == "M3/MOLECULE/S" || from == "M3 MOLEC-1 S-1")
                val_mks = val;
            else if (from == "1/S" || from == "S-1")
                val_mks = val;
            else {
                rc = -1;
                return val;
            }

            if (to == "CM3/MOLECULE/S" || to == "CM3 MOLEC-1 S-1")
                return val_mks * 1e6;
            if (to == "M3/MOLECULE/S" || to == "M3 MOLEC-1 S-1")
                return val_mks;
            if (to == "1/S" || to == "S-1")
                return val_mks;

            rc = -1;
            return val;
        }

        inline double convert_mass_units(double val, const std::string& from_u, const std::string& to_u, int& rc) {
            rc = 0;
            double val_kg = val;
            std::string from = to_upper(from_u);
            std::string to = to_upper(to_u);

            if (from == "KG" || from == "KILOGRAM")
                val_kg = val;
            else if (from == "G" || from == "GRAM")
                val_kg = val * 1e-3;
            else if (from == "MG" || from == "MILLIGRAM")
                val_kg = val * 1e-6;
            else if (from == "UG" || from == "MICROGRAM")
                val_kg = val * 1e-9;
            else if (from == "LB" || from == "POUND" || from == "LBS")
                val_kg = val * 0.45359237;
            else if (from == "OZ" || from == "OUNCE")
                val_kg = val * 0.028349523125;
            else {
                rc = -1;
                return val;
            }

            if (to == "KG" || to == "KILOGRAM")
                return val_kg;
            if (to == "G" || to == "GRAM")
                return val_kg * 1e3;
            if (to == "MG" || to == "MILLIGRAM")
                return val_kg * 1e6;
            if (to == "UG" || to == "MICROGRAM")
                return val_kg * 1e9;
            if (to == "LB" || to == "POUND" || to == "LBS")
                return val_kg / 0.45359237;
            if (to == "OZ" || to == "OUNCE")
                return val_kg / 0.028349523125;

            rc = -1;
            return val;
        }

        inline double convert_imperial(double val, const std::string& from_u, const std::string& to_u,
                                       const std::string& cat, int& rc) {
            rc = 0;
            std::string c = to_upper(cat);
            std::string from = to_upper(from_u);
            std::string to = to_upper(to_u);

            if (c == "LENGTH") {
                double val_m = val;
                if (from == "M" || from == "METER")
                    val_m = val;
                else if (from == "CM" || from == "CENTIMETER")
                    val_m = val * 1e-2;
                else if (from == "MM" || from == "MILLIMETER")
                    val_m = val * 1e-3;
                else if (from == "IN" || from == "INCH")
                    val_m = val * 0.0254;
                else if (from == "FT" || from == "FOOT" || from == "FEET")
                    val_m = val * 0.3048;
                else if (from == "YD" || from == "YARD")
                    val_m = val * 0.9144;
                else if (from == "MI" || from == "MILE")
                    val_m = val * 1609.344;
                else {
                    rc = -1;
                    return val;
                }

                if (to == "M" || to == "METER")
                    return val_m;
                if (to == "CM" || to == "CENTIMETER")
                    return val_m * 1e2;
                if (to == "MM" || to == "MILLIMETER")
                    return val_m * 1e3;
                if (to == "IN" || to == "INCH")
                    return val_m / 0.0254;
                if (to == "FT" || to == "FOOT" || to == "FEET")
                    return val_m / 0.3048;
                if (to == "YD" || to == "YARD")
                    return val_m / 0.9144;
                if (to == "MI" || to == "MILE")
                    return val_m / 1609.344;
            } else if (c == "AREA") {
                double val_m2 = val;
                if (from == "M2")
                    val_m2 = val;
                else if (from == "CM2")
                    val_m2 = val * 1e-4;
                else if (from == "IN2")
                    val_m2 = val * 0.00064516;
                else if (from == "FT2")
                    val_m2 = val * 0.09290304;
                else if (from == "ACRE" || from == "ACRES")
                    val_m2 = val * 4046.8564224;
                else {
                    rc = -1;
                    return val;
                }

                if (to == "M2")
                    return val_m2;
                if (to == "CM2")
                    return val_m2 * 1e4;
                if (to == "IN2")
                    return val_m2 / 0.00064516;
                if (to == "FT2")
                    return val_m2 / 0.09290304;
                if (to == "ACRE" || to == "ACRES")
                    return val_m2 / 4046.8564224;
            } else if (c == "VOLUME") {
                double val_m3 = val;
                if (from == "M3")
                    val_m3 = val;
                else if (from == "L" || from == "LITER" || from == "LITERS")
                    val_m3 = val * 1e-3;
                else if (from == "CM3" || from == "CC")
                    val_m3 = val * 1e-6;
                else if (from == "GAL" || from == "GALLON" || from == "GALLONS")
                    val_m3 = val * 0.003785411784;
                else if (from == "FT3")
                    val_m3 = val * 0.028316846592;
                else {
                    rc = -1;
                    return val;
                }

                if (to == "M3")
                    return val_m3;
                if (to == "L" || to == "LITER" || to == "LITERS")
                    return val_m3 * 1e3;
                if (to == "CM3" || to == "CC")
                    return val_m3 * 1e6;
                if (to == "GAL" || to == "GALLON" || to == "GALLONS")
                    return val_m3 / 0.003785411784;
                if (to == "FT3")
                    return val_m3 / 0.028316846592;
            } else if (c == "SPEED" || c == "VELOCITY") {
                double val_ms = val;
                if (from == "M/S" || from == "M S-1")
                    val_ms = val;
                else if (from == "KM/H" || from == "KMH")
                    val_ms = val / 3.6;
                else if (from == "MPH" || from == "MILE/H")
                    val_ms = val * 0.44704;
                else if (from == "KTS" || from == "KNOT" || from == "KNOTS")
                    val_ms = val * 0.5144444444;
                else {
                    rc = -1;
                    return val;
                }

                if (to == "M/S" || to == "M S-1")
                    return val_ms;
                if (to == "KM/H" || to == "KMH")
                    return val_ms * 3.6;
                if (to == "MPH" || to == "MILE/H")
                    return val_ms / 0.44704;
                if (to == "KTS" || to == "KNOT" || to == "KNOTS")
                    return val_ms / 0.5144444444;
            } else if (c == "FORCE") {
                double val_n = val;
                if (from == "N" || from == "NEWTON")
                    val_n = val;
                else if (from == "DYNE")
                    val_n = val * 1e-5;
                else if (from == "LBF" || from == "POUND-FORCE")
                    val_n = val * 4.4482216152605;
                else {
                    rc = -1;
                    return val;
                }

                if (to == "N" || to == "NEWTON")
                    return val_n;
                if (to == "DYNE")
                    return val_n * 1e5;
                if (to == "LBF" || to == "POUND-FORCE")
                    return val_n / 4.4482216152605;
            } else if (c == "ENERGY") {
                double val_j = val;
                if (from == "J" || from == "JOULE")
                    val_j = val;
                else if (from == "CAL" || from == "CALORIE")
                    val_j = val * 4.184;
                else if (from == "KCAL")
                    val_j = val * 4184.0;
                else if (from == "BTU")
                    val_j = val * 1055.05585262;
                else if (from == "WH")
                    val_j = val * 3600.0;
                else if (from == "KWH")
                    val_j = val * 3600000.0;
                else {
                    rc = -1;
                    return val;
                }

                if (to == "J" || to == "JOULE")
                    return val_j;
                if (to == "CAL" || to == "CALORIE")
                    return val_j / 4.184;
                if (to == "KCAL")
                    return val_j / 4184.0;
                if (to == "BTU")
                    return val_j / 1055.05585262;
                if (to == "WH")
                    return val_j / 3600.0;
                if (to == "KWH")
                    return val_j / 3600000.0;
            }

            rc = -1;
            return val;
        }

    } // namespace unit_conversion
} // namespace catchem
