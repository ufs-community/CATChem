#pragma once
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace catchem {

    struct SpeciesMetadata {
        // Names
        std::string short_name;
        std::string long_name;
        std::string description;
        std::vector<std::string> aliases;
        std::vector<std::string> roles;

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
        bool is_hydrophilic = true;

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

    inline std::string canonical_species_name(std::string name) {
        std::transform(name.begin(), name.end(), name.begin(),
                       [](unsigned char value) { return static_cast<char>(std::toupper(value)); });
        return name;
    }

    class MechanismDefinition {
    public:
        std::string identity;
        std::string source;
        std::size_t generation = 1;
        std::vector<SpeciesMetadata> species;
        std::unordered_set<std::string> capabilities;

        void rebuild_index() {
            name_to_index_.clear();
            role_to_index_.clear();
            if (species.empty())
                throw std::invalid_argument("Chemical mechanism must contain at least one species");
            for (std::size_t index = 0; index < species.size(); ++index) {
                const auto canonical = canonical_species_name(species[index].short_name);
                if (canonical.empty())
                    throw std::invalid_argument("Chemical mechanism contains an empty species name");
                if (!name_to_index_.emplace(canonical, index).second)
                    throw std::invalid_argument("Chemical mechanism contains duplicate species name: " +
                                                species[index].short_name);
                for (const auto& alias : species[index].aliases) {
                    const auto canonical_alias = canonical_species_name(alias);
                    if (canonical_alias.empty() || !name_to_index_.emplace(canonical_alias, index).second)
                        throw std::invalid_argument("Chemical mechanism contains duplicate or empty species alias: " +
                                                    alias);
                }
                for (const auto& role : species[index].roles) {
                    const auto canonical_role = canonical_species_name(role);
                    if (canonical_role.empty() || !role_to_index_.emplace(canonical_role, index).second)
                        throw std::invalid_argument("Chemical mechanism contains duplicate or empty species role: " +
                                                    role);
                }
            }
        }

        std::size_t index_of(const std::string& name) const {
            const auto found = name_to_index_.find(canonical_species_name(name));
            if (found == name_to_index_.end())
                throw std::out_of_range("Species is not present in the active mechanism: " + name);
            return found->second;
        }

        bool contains(const std::string& name) const {
            return name_to_index_.find(canonical_species_name(name)) != name_to_index_.end();
        }

        bool has_role(const std::string& role) const {
            return role_to_index_.find(canonical_species_name(role)) != role_to_index_.end();
        }

        std::size_t index_for_role(const std::string& role) const {
            const auto found = role_to_index_.find(canonical_species_name(role));
            if (found == role_to_index_.end())
                throw std::out_of_range("Required species role is not present in the active mechanism: " + role);
            return found->second;
        }

        bool has_capability(const std::string& capability) const {
            return capabilities.find(canonical_species_name(capability)) != capabilities.end();
        }

    private:
        std::unordered_map<std::string, std::size_t> name_to_index_;
        std::unordered_map<std::string, std::size_t> role_to_index_;
    };

} // namespace catchem
