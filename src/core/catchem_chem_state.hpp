#pragma once
#include "catchem_config_manager.hpp"
#include "catchem_interop_field.hpp"
#include "catchem_species_metadata.hpp"
#include <algorithm>
#include <array>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace catchem {

    struct ChemState {
        // Single unified 3D View (cols, levels, species)
        std::shared_ptr<InteropField<double, 3>> conc;
        std::shared_ptr<const MechanismDefinition> mechanism;

        // Species metadata database
        std::vector<SpeciesMetadata> species_list;
        std::unordered_map<std::string, int> species_name_to_index; // 0-based indexing

        // Pre-filtered category lists (0-based, in species order), one per
        // SpeciesCategory.  Read them through indices(); add a category by
        // extending the enum and the classification table in
        // load_from_config_manager().
        std::array<std::vector<int>, species_category_count> category_indices;

        const std::vector<int>& indices(SpeciesCategory category) const {
            return category_indices[static_cast<std::size_t>(category)];
        }
        std::vector<int>& indices(SpeciesCategory category) {
            return category_indices[static_cast<std::size_t>(category)];
        }

        // Cached flat C-character array of short names
        std::vector<char> species_names_c_arr;

        void load_from_config_manager(const ConfigManager& config_mgr) {
            species_list.clear();
            species_name_to_index.clear();

            for (auto& list : category_indices)
                list.clear();

            int index = 0;
            auto descriptor = std::make_shared<MechanismDefinition>();
            descriptor->identity =
                config_mgr.data.mechanism_identity.empty() ? "configured" : config_mgr.data.mechanism_identity;
            descriptor->source = config_mgr.config_file_path;
            for (const auto& capability : config_mgr.data.mechanism_capabilities)
                descriptor->capabilities.insert(canonical_species_name(capability));
            for (const auto& sp : config_mgr.data.species) {
                SpeciesMetadata meta;
                meta.short_name = sp.name;
                meta.long_name = sp.long_name.empty() ? sp.name : sp.long_name;
                meta.description = sp.description;
                meta.aliases = sp.aliases;
                meta.roles = sp.roles;

                meta.is_gas = sp.is_gas;
                // A species explicitly declared gas is never passed to
                // particle-only processes, even if a legacy species file
                // also carries a stale aerosol flag.  This preserves the
                // configured gas classification without inventing particle
                // radius or density properties.
                meta.is_aerosol = sp.is_aerosol && !meta.is_gas;
                meta.is_tracer = sp.is_tracer;
                meta.is_advected = sp.is_advected;
                meta.is_drydep = sp.is_drydep;
                meta.is_wetdep = sp.is_wetdep;
                meta.do_settling = sp.do_settling && meta.is_aerosol;
                meta.is_photolysis = sp.is_photolysis;
                meta.is_gocart_aero = sp.is_gocart_aero && !meta.is_gas;
                meta.is_dust = sp.is_dust && !meta.is_gas;
                meta.is_seasalt = sp.is_seasalt && !meta.is_gas;
                meta.is_hydrophilic = sp.is_hydrophilic;

                meta.mw_g = sp.mw_g > 0.0 ? sp.mw_g : sp.molecular_weight_kg_mol * 1000.0;
                meta.density = sp.density;
                meta.radius = sp.radius;
                meta.lower_radius = sp.lower_radius;
                meta.upper_radius = sp.upper_radius;
                meta.viscosity = sp.viscosity;

                meta.dd_f0 = sp.dd_f0;
                meta.dd_hstar = sp.dd_hstar;
                meta.dd_DvzAerSnow = sp.dd_DvzAerSnow;
                meta.dd_DvzMinVal_snow = sp.dd_DvzMinVal_snow;
                meta.dd_DvzMinVal_land = sp.dd_DvzMinVal_land;

                meta.henry_k0 = sp.henry_k0;
                meta.henry_cr = sp.henry_cr;
                meta.henry_pKa = sp.henry_pKa;
                meta.wd_retfactor = sp.wd_retfactor;
                meta.wd_LiqAndGas = sp.wd_LiqAndGas;
                meta.wd_convfacI2G = sp.wd_convfacI2G;
                meta.wd_rainouteff = sp.wd_rainouteff;
                meta.wd_reevap_frac = sp.wd_reevap_frac;

                meta.t_chem_loss = sp.t_chem_loss;
                meta.BackgroundVV = sp.BackgroundVV;
                meta.mie_name = sp.mie_name;

                species_list.push_back(meta);
                descriptor->species.push_back(meta);
                std::string canonical_name = meta.short_name;
                std::transform(canonical_name.begin(), canonical_name.end(), canonical_name.begin(),
                               [](unsigned char c) { return std::toupper(c); });
                if (canonical_name.empty())
                    throw std::invalid_argument("Chemical mechanism contains an empty species name");
                if (!species_name_to_index.emplace(canonical_name, index).second)
                    throw std::invalid_argument("Chemical mechanism contains duplicate species name: " +
                                                meta.short_name);

                // Classify species: one row per category, evaluated in species
                // order so every list keeps the same order as the catalog.
                const std::pair<SpeciesCategory, bool> classification[] = {
                    {SpeciesCategory::Gas, meta.is_gas},
                    {SpeciesCategory::Aerosol, meta.is_aerosol},
                    {SpeciesCategory::Settling, meta.is_aerosol && meta.do_settling},
                    {SpeciesCategory::Tracer, meta.is_tracer},
                    {SpeciesCategory::Advected, meta.is_advected},
                    {SpeciesCategory::DryDep, meta.is_drydep},
                    {SpeciesCategory::WetDep, meta.is_wetdep},
                    {SpeciesCategory::Photolysis, meta.is_photolysis},
                    {SpeciesCategory::Dust, meta.is_dust},
                    {SpeciesCategory::SeaSalt, meta.is_seasalt},
                };
                for (const auto& [category, member] : classification)
                    if (member)
                        indices(category).push_back(index);

                index++;
            }
            descriptor->rebuild_index();
            mechanism = descriptor;

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

        void load_species_config(const std::string& filename, ConfigManager* cfg = nullptr) {
            if (cfg) {
                cfg->load_species_file(filename);
                load_from_config_manager(*cfg);
            } else {
                ConfigManager temp_cfg;
                temp_cfg.load_species_file(filename);
                load_from_config_manager(temp_cfg);
            }
        }
    };

} // namespace catchem
