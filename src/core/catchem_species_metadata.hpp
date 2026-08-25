#pragma once
#include <string>
#include <vector>

namespace catchem {

    struct SpeciesMetadata {
        // Names
        std::string short_name;
        std::string long_name;
        std::string description;

        // Classification switches
        bool is_gas = false;
        bool is_aerosol = false;
        bool is_tracer = false;
        bool is_advected = true;
        bool is_drydep = false;
        bool is_wetdep = false;
        bool is_photolysis = false;
        bool is_gocart_aero = false;
        bool is_dust = false;
        bool is_seasalt = false;

        // Physical / Numerical properties
        double mw_g = 0.0;
        double density = 0.0;
        double radius = 0.0;
        double lower_radius = 0.0;
        double upper_radius = 0.0;
        double viscosity = 0.0;

        // Dry deposition parameters
        double dd_f0 = 0.0;
        double dd_hstar = 0.0;
        double dd_DvzAerSnow = 0.0;
        double dd_DvzMinVal_snow = 0.0;
        double dd_DvzMinVal_land = 0.0;

        // Wet deposition parameters
        double henry_k0 = 0.0;
        double henry_cr = 0.0;
        double henry_pKa = 0.0;
        double wd_retfactor = 0.0;
        bool wd_LiqAndGas = false;
        double wd_convfacI2G = 0.0;
        std::vector<double> wd_rainouteff = {0.0, 0.0, 0.0};
        double wd_reevap_frac = 0.5;

        // Chemical loss rate and background volume-mixing ratio
        double t_chem_loss = -1.0;
        double BackgroundVV = 1.0e-20;
        std::string mie_name;
    };

} // namespace catchem
