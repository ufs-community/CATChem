#pragma once
#include "catchem_interop_field.hpp"
#include "catchem_species_metadata.hpp"
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace catchem {

    struct ChemState {
        // Single unified 3D View (cols, levels, species)
        std::shared_ptr<InteropField<double, 3>> conc;

        // Species metadata database
        std::vector<SpeciesMetadata> species_list;
        std::unordered_map<std::string, int> species_name_to_index; // 0-based indexing

        // Pre-filtered category lists (0-based)
        std::vector<int> gas_indices;
        std::vector<int> aerosol_indices;
        std::vector<int> tracer_indices;
        std::vector<int> advected_indices;
        std::vector<int> drydep_indices;
        std::vector<int> wetdep_indices;
        std::vector<int> photolysis_indices;
        std::vector<int> dust_indices;
        std::vector<int> seasalt_indices;

        // Cached flat C-character array of short names
        std::vector<char> species_names_c_arr;

        void load_species_config(const std::string& filename) {
            YAML::Node config = YAML::LoadFile(filename);
            species_list.clear();
            species_name_to_index.clear();

            gas_indices.clear();
            aerosol_indices.clear();
            tracer_indices.clear();
            advected_indices.clear();
            drydep_indices.clear();
            wetdep_indices.clear();
            photolysis_indices.clear();
            dust_indices.clear();
            seasalt_indices.clear();

            int index = 0;
            for (auto const& item : config) {
                YAML::Node val = item;
                std::string key = val["name"].as<std::string>();

                SpeciesMetadata meta;
                meta.short_name = key;
                meta.long_name = key;
                meta.description = val["__description"] ? val["__description"].as<std::string>() : "";

                meta.is_gas = val["__is_gas"] ? val["__is_gas"].as<bool>() : false;
                meta.is_aerosol = val["__is_aerosol"] ? val["__is_aerosol"].as<bool>() : false;
                meta.is_tracer = val["__is_tracer"] ? val["__is_tracer"].as<bool>() : false;
                meta.is_advected = val["__is_advected"] ? val["__is_advected"].as<bool>() : true;
                meta.is_drydep = val["__is_drydep"] ? val["__is_drydep"].as<bool>() : false;
                meta.is_wetdep = val["__is_wetdep"] ? val["__is_wetdep"].as<bool>() : false;
                meta.is_photolysis = val["__is_photolysis"] ? val["__is_photolysis"].as<bool>() : false;
                meta.is_dust = val["__is_dust"] ? val["__is_dust"].as<bool>() : false;
                meta.is_seasalt = val["__is_seasalt"] ? val["__is_seasalt"].as<bool>() : false;

                meta.mw_g =
                    val["molecular weight [kg mol-1]"] ? val["molecular weight [kg mol-1]"].as<double>() * 1000.0 : 0.0;
                meta.density = val["__density"] ? val["__density"].as<double>() : 0.0;
                meta.radius = val["__radius"] ? val["__radius"].as<double>() : 0.0;
                meta.lower_radius = val["__lower_radius"] ? val["__lower_radius"].as<double>() : 0.0;
                meta.upper_radius = val["__upper_radius"] ? val["__upper_radius"].as<double>() : 0.0;
                meta.viscosity = val["__viscosity"] ? val["__viscosity"].as<double>() : 0.0;

                meta.dd_f0 = val["__dd_f0"] ? val["__dd_f0"].as<double>() : 0.0;
                meta.dd_hstar = val["__dd_hstar"] ? val["__dd_hstar"].as<double>() : 0.0;
                meta.dd_DvzAerSnow = val["__dd_DvzAerSnow"] ? val["__dd_DvzAerSnow"].as<double>() : 0.0;
                meta.dd_DvzMinVal_snow = val["__dd_DvzMinVal_snow"] ? val["__dd_DvzMinVal_snow"].as<double>() : 0.0;
                meta.dd_DvzMinVal_land = val["__dd_DvzMinVal_land"] ? val["__dd_DvzMinVal_land"].as<double>() : 0.0;

                // Wet deposition parameters
                meta.henry_k0 = val["__henry_k0"] ? val["__henry_k0"].as<double>() : 0.0;
                meta.henry_cr = val["__henry_cr"] ? val["__henry_cr"].as<double>() : 0.0;
                meta.henry_pKa = val["__henry_pKa"] ? val["__henry_pKa"].as<double>() : 0.0;
                meta.wd_retfactor = val["__wd_retfactor"] ? val["__wd_retfactor"].as<double>() : 0.0;
                meta.wd_LiqAndGas = val["__wd_LiqAndGas"] ? val["__wd_LiqAndGas"].as<bool>() : false;
                meta.wd_convfacI2G = val["__wd_convfacI2G"] ? val["__wd_convfacI2G"].as<double>() : 0.0;

                if (val["__wd_rainouteff"]) {
                    meta.wd_rainouteff = val["__wd_rainouteff"].as<std::vector<double>>();
                }
                meta.wd_reevap_frac = val["__wd_reevap_frac"] ? val["__wd_reevap_frac"].as<double>() : 0.5;

                // GOCART carbon chemical loss rate [days]
                meta.t_chem_loss = val["__t_chem_loss"] ? val["__t_chem_loss"].as<double>() : -1.0;

                // Background VMR concentration
                meta.BackgroundVV = val["BackgroundVV"] ? val["BackgroundVV"].as<double>() : 1.0e-20;
                meta.mie_name = val["__mie_name"] ? val["__mie_name"].as<std::string>() : "";

                species_list.push_back(meta);
                species_name_to_index[key] = index;

                // Classify species
                if (meta.is_gas)
                    gas_indices.push_back(index);
                if (meta.is_aerosol)
                    aerosol_indices.push_back(index);
                if (meta.is_tracer)
                    tracer_indices.push_back(index);
                if (meta.is_advected)
                    advected_indices.push_back(index);
                if (meta.is_drydep)
                    drydep_indices.push_back(index);
                if (meta.is_wetdep)
                    wetdep_indices.push_back(index);
                if (meta.is_photolysis)
                    photolysis_indices.push_back(index);
                if (meta.is_dust)
                    dust_indices.push_back(index);
                if (meta.is_seasalt)
                    seasalt_indices.push_back(index);

                index++;
            }

            // Pre-compute and cache flat C-linkable species name character array
            species_names_c_arr.assign(species_list.size() * 32, ' ');
            for (size_t i = 0; i < species_list.size(); ++i) {
                std::string name = species_list[i].short_name;
                for (auto& c : name)
                    c = std::toupper(c);
                for (size_t j = 0; j < name.size() && j < 32; ++j) {
                    species_names_c_arr[i * 32 + j] = name[j];
                }
            }
        }
    };

} // namespace catchem
