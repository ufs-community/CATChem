#pragma once
#include "catchem_interop_field.hpp"
#include <algorithm>
#include <cctype>
#include <initializer_list>
#include <memory>
#include <string>
#include <unordered_map>

namespace catchem {

    struct MetState {
        // Dynamic map storage for 2D and 3D meteorological fields
        std::unordered_map<std::string, std::shared_ptr<InteropField<double, 3>>> fields_3d;
        std::unordered_map<std::string, std::shared_ptr<InteropField<double, 2>>> fields_2d;

        // Convenient direct accessors for standard core fields
        std::shared_ptr<InteropField<double, 3>> T;          // Temperature [K]
        std::shared_ptr<InteropField<double, 3>> QV;         // Specific humidity [kg/kg]
        std::shared_ptr<InteropField<double, 3>> RH;         // Relative humidity [0-1]
        std::shared_ptr<InteropField<double, 3>> PMID;       // Mid-level pressure [Pa]
        std::shared_ptr<InteropField<double, 3>> PEDGE;      // Edge-level pressure [Pa]
        std::shared_ptr<InteropField<double, 3>> AIRDEN;     // Wet air density [kg/m³]
        std::shared_ptr<InteropField<double, 3>> AIRDEN_DRY; // Dry air density [kg/m³]
        std::shared_ptr<InteropField<double, 3>> BXHEIGHT;   // Layer thickness height [m]

        std::shared_ptr<InteropField<double, 2>> PS;    // Surface pressure [Pa]
        std::shared_ptr<InteropField<double, 2>> TS;    // Surface temperature [K]
        std::shared_ptr<InteropField<double, 2>> PBLH;  // Boundary layer height [m]
        std::shared_ptr<InteropField<double, 2>> USTAR; // Friction velocity [m/s]
        std::shared_ptr<InteropField<double, 2>> HFLUX; // Sensible heat flux [W/m²]
        std::shared_ptr<InteropField<double, 2>> OBK;   // Monin-Obukhov length [m]
        std::shared_ptr<InteropField<double, 2>> LAT;   // Latitude [deg]
        std::shared_ptr<InteropField<double, 2>> LON;   // Longitude [deg]

        // Dynamic 3D field lookup supporting candidate alias lists
        std::shared_ptr<InteropField<double, 3>> get_3d(std::initializer_list<const char*> candidates) const {
            for (const char* name : candidates) {
                if (!name)
                    continue;
                auto it = fields_3d.find(name);
                if (it != fields_3d.end() && it->second) {
                    return it->second;
                }
            }
            return nullptr;
        }

        // Dynamic 2D field lookup supporting candidate alias lists
        std::shared_ptr<InteropField<double, 2>> get_2d(std::initializer_list<const char*> candidates) const {
            for (const char* name : candidates) {
                if (!name)
                    continue;
                auto it = fields_2d.find(name);
                if (it != fields_2d.end() && it->second) {
                    return it->second;
                }
            }
            return nullptr;
        }

        // Check if any candidate 3D field is bound
        bool has_3d(std::initializer_list<const char*> candidates) const { return get_3d(candidates) != nullptr; }

        // Check if any candidate 2D field is bound
        bool has_2d(std::initializer_list<const char*> candidates) const { return get_2d(candidates) != nullptr; }

        // Register/bind a 2D field dynamically into the map and sync convenience pointers
        void bind_2d_field(const std::string& name, std::shared_ptr<InteropField<double, 2>> field) {
            if (!field || name.empty())
                return;
            fields_2d[name] = field;

            // Generate upper and lower case aliases for key lookup
            std::string upper_name = name;
            std::transform(upper_name.begin(), upper_name.end(), upper_name.begin(),
                           [](unsigned char c) { return std::toupper(c); });
            fields_2d[upper_name] = field;

            std::string lower_name = name;
            std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            fields_2d[lower_name] = field;

            // Sync standard convenience pointers dynamically
            if (upper_name == "PS" || upper_name == "SURFACE_PRESSURE")
                PS = field;
            else if (upper_name == "TS" || upper_name == "SST" || upper_name == "SKIN_TEMPERATURE")
                TS = field;
            else if (upper_name == "PBLH" || upper_name == "HPBL")
                PBLH = field;
            else if (upper_name == "USTAR" || upper_name == "FRICTION_VELOCITY")
                USTAR = field;
            else if (upper_name == "HFLUX")
                HFLUX = field;
            else if (upper_name == "OBK" || upper_name == "OL")
                OBK = field;
            else if (upper_name == "LAT" || upper_name == "LATITUDE")
                LAT = field;
            else if (upper_name == "LON" || upper_name == "LONGITUDE")
                LON = field;
        }

        // Register/bind a 3D field dynamically into the map and sync convenience pointers
        void bind_3d_field(const std::string& name, std::shared_ptr<InteropField<double, 3>> field) {
            if (!field || name.empty())
                return;
            fields_3d[name] = field;

            // Generate upper and lower case aliases for key lookup
            std::string upper_name = name;
            std::transform(upper_name.begin(), upper_name.end(), upper_name.begin(),
                           [](unsigned char c) { return std::toupper(c); });
            fields_3d[upper_name] = field;

            std::string lower_name = name;
            std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            fields_3d[lower_name] = field;

            // Sync standard convenience pointers dynamically
            if (upper_name == "T" || upper_name == "TEMPERATURE" || upper_name == "TEMP")
                T = field;
            else if (upper_name == "QV" || upper_name == "SPHUM")
                QV = field;
            else if (upper_name == "RH" || upper_name == "RELATIVE_HUMIDITY")
                RH = field;
            else if (upper_name == "PMID" || upper_name == "PRESSURE_MID")
                PMID = field;
            else if (upper_name == "PEDGE" || upper_name == "PRESSURE_EDGE")
                PEDGE = field;
            else if (upper_name == "AIRDEN" || upper_name == "AIR_DENSITY") {
                AIRDEN = field;
            } else if (upper_name == "AIRDEN_DRY" || upper_name == "AIR_DENSITY_DRY") {
                AIRDEN_DRY = field;
            } else if (upper_name == "BXHEIGHT" || upper_name == "DZ" || upper_name == "DELZ")
                BXHEIGHT = field;
        }

        void register_fields() {
            if (T)
                bind_3d_field("T", T);
            if (QV)
                bind_3d_field("QV", QV);
            if (RH)
                bind_3d_field("RH", RH);
            if (PMID)
                bind_3d_field("PMID", PMID);
            if (PEDGE)
                bind_3d_field("PEDGE", PEDGE);
            if (AIRDEN)
                bind_3d_field("AIRDEN", AIRDEN);
            if (AIRDEN_DRY)
                bind_3d_field("AIRDEN_DRY", AIRDEN_DRY);
            if (BXHEIGHT)
                bind_3d_field("BXHEIGHT", BXHEIGHT);

            if (PS)
                bind_2d_field("PS", PS);
            if (TS)
                bind_2d_field("TS", TS);
            if (PBLH)
                bind_2d_field("PBLH", PBLH);
            if (USTAR)
                bind_2d_field("USTAR", USTAR);
            if (HFLUX)
                bind_2d_field("HFLUX", HFLUX);
            if (OBK)
                bind_2d_field("OBK", OBK);
            if (LAT)
                bind_2d_field("LAT", LAT);
            if (LON)
                bind_2d_field("LON", LON);
        }
    };

} // namespace catchem
