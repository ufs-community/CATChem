#pragma once
#include "catchem_interop_field.hpp"
#include <memory>
#include <string>
#include <unordered_map>

namespace catchem {

    struct MetState {
        // 3D Volumetric fields
        std::shared_ptr<InteropField<double, 3>> T;          // Temperature [K]
        std::shared_ptr<InteropField<double, 3>> QV;         // Specific humidity [kg/kg]
        std::shared_ptr<InteropField<double, 3>> RH;         // Relative humidity [0-1]
        std::shared_ptr<InteropField<double, 3>> PMID;       // Mid-level pressure [Pa]
        std::shared_ptr<InteropField<double, 3>> PEDGE;      // Edge-level pressure [Pa]
        std::shared_ptr<InteropField<double, 3>> AIRDEN;     // Wet air density [kg/m³]
        std::shared_ptr<InteropField<double, 3>> AIRDEN_DRY; // Dry air density [kg/m³]
        std::shared_ptr<InteropField<double, 3>> BXHEIGHT;   // Layer thickness height [m]

        // 2D Surface fields
        std::shared_ptr<InteropField<double, 2>> PS;    // Surface pressure [Pa]
        std::shared_ptr<InteropField<double, 2>> TS;    // Surface temperature [K]
        std::shared_ptr<InteropField<double, 2>> PBLH;  // Boundary layer height [m]
        std::shared_ptr<InteropField<double, 2>> USTAR; // Friction velocity [m/s]
        std::shared_ptr<InteropField<double, 2>> HFLUX; // Sensible heat flux [W/m²]
        std::shared_ptr<InteropField<double, 2>> OBK;   // Monin-Obukhov length [m]
        std::shared_ptr<InteropField<double, 2>> LAT;   // Latitude [deg]
        std::shared_ptr<InteropField<double, 2>> LON;   // Longitude [deg]

        std::unordered_map<std::string, std::shared_ptr<InteropField<double, 3>>> fields_3d;
        std::unordered_map<std::string, std::shared_ptr<InteropField<double, 2>>> fields_2d;

        void register_fields() {
            fields_3d["T"] = T;
            fields_3d["QV"] = QV;
            fields_3d["RH"] = RH;
            fields_3d["PMID"] = PMID;
            fields_3d["PEDGE"] = PEDGE;
            fields_3d["AIRDEN"] = AIRDEN;
            fields_3d["AIRDEN_DRY"] = AIRDEN_DRY;
            fields_3d["BXHEIGHT"] = BXHEIGHT;

            fields_2d["PS"] = PS;
            fields_2d["TS"] = TS;
            fields_2d["PBLH"] = PBLH;
            fields_2d["USTAR"] = USTAR;
            fields_2d["HFLUX"] = HFLUX;
            fields_2d["OBK"] = OBK;
            fields_2d["LAT"] = LAT;
            fields_2d["LON"] = LON;
        }
    };

} // namespace catchem
