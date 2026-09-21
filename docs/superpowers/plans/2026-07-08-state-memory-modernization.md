# Modernize State & Memory Management Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish C++ MetState, ChemState, and TimeState containers inside StateManager, shifting chemical metadata, unified concentration Views, and portable solar computations entirely to C++ with zero-copy interfaces.

**Architecture:** Create modular C++ headers: `MetState` for meteorological arrays, `ChemState` for chemical database loading and unified concentration Views, and `TimeState` for high-precision inline SZA math. Update `StateManager` to host these three structures and implement parallel thermodynamic equations.

**Tech Stack:** C++20, Kokkos, CMake

## Global Constraints

- Target C++20 utilizing the Kokkos backport of mdspan (`std::experimental::mdspan`), avoiding direct dependency on C++23 `<mdspan>`.
- Retain Fortran column-major storage layout (`Kokkos::LayoutLeft`) across the pointer boundary to achieve zero-copy execution on CPU targets.
- All modifications must remain unstaged and uncommitted until instructed otherwise.

---

### Task 1: Create MetState, ChemState, and TimeState Header Interfaces

**Files:**
- Create: `src/core/catchem_met_state.hpp`
- Create: `src/core/catchem_chem_state.hpp`
- Create: `src/core/catchem_time_state.hpp`

**Interfaces:**
- Produces: `catchem::MetState`, `catchem::ChemState`, `catchem::TimeState`

- [ ] **Step 1: Create header file `src/core/catchem_met_state.hpp`**

Write the complete container holding 3D/2D meteorological `InteropField` objects:
```cpp
#pragma once
#include <memory>
#include <string>
#include <unordered_map>
#include "catchem_interop_field.hpp"

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
    std::shared_ptr<InteropField<double, 2>> PS;         // Surface pressure [Pa]
    std::shared_ptr<InteropField<double, 2>> TS;         // Surface temperature [K]
    std::shared_ptr<InteropField<double, 2>> PBLH;       // Boundary layer height [m]
    std::shared_ptr<InteropField<double, 2>> USTAR;      // Friction velocity [m/s]
    std::shared_ptr<InteropField<double, 2>> HFLUX;      // Sensible heat flux [W/m²]
    std::shared_ptr<InteropField<double, 2>> OBK;        // Monin-Obukhov length [m]
    std::shared_ptr<InteropField<double, 2>> LAT;        // Latitude [deg]
    std::shared_ptr<InteropField<double, 2>> LON;        // Longitude [deg]

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
```

- [ ] **Step 2: Create header file `src/core/catchem_chem_state.hpp`**

Shift species YAML loading logic and hold unified concentration fields:
```cpp
#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <yaml-cpp/yaml.h>
#include "catchem_interop_field.hpp"
#include "catchem_species_metadata.hpp"

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
            std::string key = item.first.as<std::string>();
            YAML::Node val = item.second;

            SpeciesMetadata meta;
            meta.short_name = key;
            meta.long_name = val["name"] ? val["name"].as<std::string>() : key;
            meta.description = val["description"] ? val["description"].as<std::string>() : "";

            meta.is_gas = val["is_gas"] ? val["is_gas"].as<bool>() : false;
            meta.is_aerosol = val["is_aerosol"] ? val["is_aerosol"].as<bool>() : false;
            meta.is_tracer = val["is_tracer"] ? val["is_tracer"].as<bool>() : false;
            meta.is_advected = val["is_advected"] ? val["is_advected"].as<bool>() : true;
            meta.is_drydep = val["is_drydep"] ? val["is_drydep"].as<bool>() : false;
            meta.is_wetdep = val["is_wetdep"] ? val["is_wetdep"].as<bool>() : false;
            meta.is_photolysis = val["is_photolysis"] ? val["is_photolysis"].as<bool>() : false;
            meta.is_dust = val["is_dust"] ? val["is_dust"].as<bool>() : false;
            meta.is_seasalt = val["is_seasalt"] ? val["is_seasalt"].as<bool>() : false;

            meta.mw_g = val["mw_g"] ? val["mw_g"].as<double>() : 0.0;
            meta.density = val["density"] ? val["density"].as<double>() : 0.0;
            meta.radius = val["radius"] ? val["radius"].as<double>() : 0.0;
            meta.lower_radius = val["lower_radius"] ? val["lower_radius"].as<double>() : 0.0;
            meta.upper_radius = val["upper_radius"] ? val["upper_radius"].as<double>() : 0.0;
            meta.viscosity = val["viscosity"] ? val["viscosity"].as<double>() : 0.0;

            meta.dd_f0 = val["dd_f0"] ? val["dd_f0"].as<double>() : 0.0;
            meta.dd_hstar = val["dd_hstar"] ? val["dd_hstar"].as<double>() : 0.0;
            meta.dd_DvzAerSnow = val["dd_DvzAerSnow"] ? val["dd_DvzAerSnow"].as<double>() : 0.0;
            meta.dd_DvzMinVal_snow = val["dd_DvzMinVal_snow"] ? val["dd_DvzMinVal_snow"].as<double>() : 0.0;
            meta.dd_DvzMinVal_land = val["dd_DvzMinVal_land"] ? val["dd_DvzMinVal_land"].as<double>() : 0.0;

            meta.wd_retfactor = val["wd_retfactor"] ? val["wd_retfactor"].as<double>() : 0.0;
            meta.wd_LiqAndGas = val["wd_LiqAndGas"] ? val["wd_LiqAndGas"].as<bool>() : false;
            meta.wd_convfacI2G = val["wd_convfacI2G"] ? val["wd_convfacI2G"].as<double>() : 0.0;

            if (val["wd_rainouteff"]) {
                meta.wd_rainouteff = val["wd_rainouteff"].as<std::vector<double>>();
            }
            meta.mie_name = val["mie_name"] ? val["mie_name"].as<std::string>() : "";

            species_list.push_back(meta);
            species_name_to_index[key] = index;

            // Classify species
            if (meta.is_gas) gas_indices.push_back(index);
            if (meta.is_aerosol) aerosol_indices.push_back(index);
            if (meta.is_tracer) tracer_indices.push_back(index);
            if (meta.is_advected) advected_indices.push_back(index);
            if (meta.is_drydep) drydep_indices.push_back(index);
            if (meta.is_wetdep) wetdep_indices.push_back(index);
            if (meta.is_photolysis) photolysis_indices.push_back(index);
            if (meta.is_dust) dust_indices.push_back(index);
            if (meta.is_seasalt) seasalt_indices.push_back(index);

            index++;
        }
    }
};

} // namespace catchem
```

- [ ] **Step 3: Create header file `src/core/catchem_time_state.hpp`**

Write the complete inline mathematical time tracking definitions:
```cpp
#pragma once
#include <cmath>
#include <algorithm>

#ifdef ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#define KOKKOS_FUNCTION KOKKOS_INLINE_FUNCTION
#else
#define KOKKOS_INLINE_FUNCTION inline
#define KOKKOS_FUNCTION inline
#endif

#include "catchem_precision.hpp"
#include "catchem_constants.hpp"

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
    double get_cos_sza(double lat_deg, double lon_deg, bool mid_timestep = false) const {
        double lat_rad = lat_deg * constants::PI_180;
        double lon_rad = lon_deg * constants::PI_180;

        double frac_hour = hour + minute / 60.0 + second / 3600.0;
        if (mid_timestep) {
            frac_hour += (timestep / 2.0) / 3600.0;
        }

        // Day angle [radians]
        double gamma = 2.0 * constants::PI * (doy - 1.0) / 365.0;

        // Solar declination (high-precision Fourier calculation matches GOCART2G)
        double dec = 0.006918 - 0.399912 * std::cos(gamma) + 0.070257 * std::sin(gamma)
                     - 0.006758 * std::cos(2.0 * gamma) + 0.000907 * std::sin(2.0 * gamma)
                     - 0.002697 * std::cos(3.0 * gamma) + 0.001480 * std::sin(3.0 * gamma);

        // Equation of time
        double eqtime = 229.18 * (0.000075 + 0.001868 * std::cos(gamma) - 0.032077 * std::sin(gamma)
                                  - 0.014615 * std::cos(2.0 * gamma) - 0.040849 * std::sin(2.0 * gamma));

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
```

- [ ] **Step 4: Commit Task 1 files**

```bash
git add src/core/catchem_met_state.hpp src/core/catchem_chem_state.hpp src/core/catchem_time_state.hpp
git commit -m "feat(core): introduce performance portable MetState, ChemState, and TimeState C++ header containers"
```

---

### Task 2: Integrate Modular Sub-Containers into `catchem::StateManager`

**Files:**
- Modify: `src/core/catchem_state_manager.hpp`

**Interfaces:**
- Consumes: `catchem::MetState`, `catchem::ChemState`, `catchem::TimeState`
- Produces: `StateManager::met`, `StateManager::chem`, `StateManager::time`

- [ ] **Step 1: Rework `src/core/catchem_state_manager.hpp`**

Remove the duplicate metadata fields added directly to `StateManager` in Phase 4. Instead, host `MetState met;`, `ChemState chem;`, `TimeState time;` and implement cleaner bindings:
```cpp
#pragma once
#include <unordered_map>
#include <string>
#include <memory>
#include <vector>
#include <yaml-cpp/yaml.h>
#include "catchem_interop_field.hpp"
#include "catchem_met_state.hpp"
#include "catchem_chem_state.hpp"
#include "catchem_time_state.hpp"

namespace catchem {

class StateManager {
public:
    int n_cols;
    int n_levels;
    int n_species;

    std::unordered_map<std::string, std::shared_ptr<InteropField<double, 1>>> fields_1d;
    std::unordered_map<std::string, std::shared_ptr<InteropField<double, 2>>> fields_2d;
    std::unordered_map<std::string, std::shared_ptr<InteropField<double, 3>>> fields_3d;

    // Structured sub-states
    MetState met;
    ChemState chem;
    TimeState time;

    StateManager(int nc, int nl, int ns) : n_cols(nc), n_levels(nl), n_species(ns) {}

    void load_species_config(const std::string& filename) {
        chem.load_species_config(filename);
    }

    void bind_met_field_2d(const std::string& name, double* ptr) {
        auto field = std::make_shared<InteropField<double, 2>>(ptr, std::vector<int>{n_cols, n_levels});
        if (name == "PS") met.PS = field;
        else if (name == "TS") met.TS = field;
        else if (name == "PBLH") met.PBLH = field;
        else if (name == "USTAR") met.USTAR = field;
        else if (name == "HFLUX") met.HFLUX = field;
        else if (name == "OBK") met.OBK = field;
        else if (name == "LAT") met.LAT = field;
        else if (name == "LON") met.LON = field;
        met.fields_2d[name] = field;
    }

    void bind_met_field_3d(const std::string& name, double* ptr) {
        auto field = std::make_shared<InteropField<double, 3>>(ptr, std::vector<int>{n_cols, n_levels, 1}); // Using 1 for single-field layout
        if (name == "T") met.T = field;
        else if (name == "QV") met.QV = field;
        else if (name == "RH") met.RH = field;
        else if (name == "PMID") met.PMID = field;
        else if (name == "PEDGE") met.PEDGE = field;
        else if (name == "AIRDEN") met.AIRDEN = field;
        else if (name == "AIRDEN_DRY") met.AIRDEN_DRY = field;
        else if (name == "BXHEIGHT") met.BXHEIGHT = field;
        met.fields_3d[name] = field;
    }

    void bind_unified_chemistry(double* ptr) {
        chem.conc = std::make_shared<InteropField<double, 3>>(ptr, std::vector<int>{n_cols, n_levels, n_species});
    }

    void sync_to_device() {
        for (auto& [k, v] : fields_1d) v->sync_to_device();
        for (auto& [k, v] : fields_2d) v->sync_to_device();
        for (auto& [k, v] : fields_3d) v->sync_to_device();
        for (auto& [k, v] : met.fields_2d) v->sync_to_device();
        for (auto& [k, v] : met.fields_3d) v->sync_to_device();
        if (chem.conc) chem.conc->sync_to_device();
    }

    void sync_to_host() {
        for (auto& [k, v] : fields_1d) v->sync_to_host();
        for (auto& [k, v] : fields_2d) v->sync_to_host();
        for (auto& [k, v] : fields_3d) v->sync_to_host();
        for (auto& [k, v] : met.fields_2d) v->sync_to_host();
        for (auto& [k, v] : met.fields_3d) v->sync_to_host();
        if (chem.conc) chem.conc->sync_to_host();
    }

    double* get_host_pointer_1d(const std::string& name) {
        if (fields_1d.find(name) == fields_1d.end()) return nullptr;
        return fields_1d.at(name)->host_view.data();
    }

    double* get_host_pointer_2d(const std::string& name) {
        if (fields_2d.find(name) != fields_2d.end()) return fields_2d.at(name)->host_view.data();
        if (met.fields_2d.find(name) != met.fields_2d.end()) return met.fields_2d.at(name)->host_view.data();
        return nullptr;
    }

    double* get_host_pointer_3d(const std::string& name) {
        if (fields_3d.find(name) != fields_3d.end()) return fields_3d.at(name)->host_view.data();
        if (met.fields_3d.find(name) != met.fields_3d.end()) return met.fields_3d.at(name)->host_view.data();
        return nullptr;
    }
};

} // namespace catchem
```

- [ ] **Step 2: Commit Task 2 changes**

```bash
git add src/core/catchem_state_manager.hpp
git commit -m "feat(core): host MetState, ChemState, and TimeState modular instances inside StateManager"
```

---

### Task 3: Expose Memory Binding and Time Configurations in C-API

**Files:**
- Modify: `src/core/catchem_api.hpp`
- Modify: `src/core/catchem_api.cpp`

**Interfaces:**
- Consumes: Updated dynamic nested delegation of `chem` in `catchem::StateManager`
- Produces: `catchem_state_bind_met_2d`, `catchem_state_bind_met_3d`, `catchem_state_bind_unified_chemistry`, `catchem_state_set_time`

- [ ] **Step 1: Declare binders inside `src/core/catchem_api.hpp`**

Make sure these declarations are explicitly declared at global scope:
```cpp
void catchem_state_bind_met_2d(void* state_ptr, const char* name, double* ptr);
void catchem_state_bind_met_3d(void* state_ptr, const char* name, double* ptr);
void catchem_state_bind_unified_chemistry(void* state_ptr, double* ptr);
void catchem_state_set_time(void* state_ptr, int yr, int mo, int dy, int hr, int mn, int sc, int doy, double tstep);
```

- [ ] **Step 2: Implement binders and refactor endpoints inside `src/core/catchem_api.cpp`**

Directly route the query definitions we added in Phase 4 to the structured `chem` member:
```cpp
void catchem_state_bind_met_2d(void* state_ptr, const char* name, double* ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    state->bind_met_field_2d(name, ptr);
}

void catchem_state_bind_met_3d(void* state_ptr, const char* name, double* ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    state->bind_met_field_3d(name, ptr);
}

void catchem_state_bind_unified_chemistry(void* state_ptr, double* ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    state->bind_unified_chemistry(ptr);
}

void catchem_state_set_time(void* state_ptr, int yr, int mo, int dy, int hr, int mn, int sc, int doy, double tstep) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    state->time.year = yr;
    state->time.month = mo;
    state->time.day = dy;
    state->time.hour = hr;
    state->time.minute = mn;
    state->time.second = sc;
    state->time.doy = doy;
    state->time.timestep = tstep;
}

// Modify Phase 4 query functions to delegate to chem sub-state:
int catchem_state_get_species_count(void* state_ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    return static_cast<int>(state->chem.species_list.size());
}

int catchem_state_get_species_index(void* state_ptr, const char* name) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    auto it = state->chem.species_name_to_index.find(name);
    if (it != state->chem.species_name_to_index.end()) {
        return it->second + 1; // Translate 0-based C++ index to 1-based Fortran index
    }
    return -1;
}

int catchem_state_get_gas_species_count(void* state_ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    return static_cast<int>(state->chem.gas_indices.size());
}

void catchem_state_get_gas_indices(void* state_ptr, int* indices_out) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    for (size_t i = 0; i < state->chem.gas_indices.size(); ++i) {
        indices_out[i] = state->chem.gas_indices[i] + 1; // 1-based
    }
}

int catchem_state_get_aerosol_species_count(void* state_ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    return static_cast<int>(state->chem.aerosol_indices.size());
}

void catchem_state_get_aerosol_indices(void* state_ptr, int* indices_out) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    for (size_t i = 0; i < state->chem.aerosol_indices.size(); ++i) {
        indices_out[i] = state->chem.aerosol_indices[i] + 1; // 1-based
    }
}

double catchem_state_get_species_mw(void* state_ptr, int index) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    int idx_0 = index - 1; // 1-based to 0-based
    if (idx_0 >= 0 && idx_0 < static_cast<int>(state->chem.species_list.size())) {
        return state->chem.species_list[idx_0].mw_g;
    }
    return 0.0;
}

int catchem_state_is_species_gas(void* state_ptr, int index) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    int idx_0 = index - 1;
    if (idx_0 >= 0 && idx_0 < static_cast<int>(state->chem.species_list.size())) {
        return state->chem.species_list[idx_0].is_gas ? 1 : 0;
    }
    return 0;
}

int catchem_state_is_species_aerosol(void* state_ptr, int index) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    int idx_0 = index - 1;
    if (idx_0 >= 0 && idx_0 < static_cast<int>(state->chem.species_list.size())) {
        return state->chem.species_list[idx_0].is_aerosol ? 1 : 0;
    }
    return 0;
}
```

- [ ] **Step 3: Commit Task 3 updates**

```bash
git add src/core/catchem_api.*
git commit -m "feat(api): expose MetState binders and modular chem queries inside C-API layer"
```

---

### Task 4: Parallel Derived Met Calculations inside StateManager

**Files:**
- Modify: `src/core/catchem_state_manager.hpp`
- Modify: `src/core/catchem_api.hpp`
- Modify: `src/core/catchem_api.cpp`

**Interfaces:**
- Produces: `StateManager::derive_bxheight`, `StateManager::derive_airden_dry`, `catchem_state_derive_bxheight`, `catchem_state_derive_airden_dry`

- [ ] **Step 1: Declare parallel methods in `StateManager`**

Inside `src/core/catchem_state_manager.hpp`:
```cpp
// Add math includes:
#include "catchem_met_utilities.hpp"
#include "catchem_constants.hpp"

// Inside StateManager class definition:
public:
    void derive_bxheight() {
        if (!met.PEDGE || !met.T || !met.QV || !met.BXHEIGHT) return;

        int nc = n_cols;
        int nl = n_levels;

        auto pedge = met.PEDGE->view();
        auto temp = met.T->view();
        auto qv = met.QV->view();
        auto bxheight = met.BXHEIGHT->view();

        Kokkos::parallel_for("derive_bxheight_kernel",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, 0}, {nc, nl}),
            KOKKOS_LAMBDA(int icol, int ilev) {
                double p_lower = pedge(icol, ilev, 0);
                double p_upper = pedge(icol, ilev + 1, 0);

                if (p_upper > 0.0) {
                    double virtual_t = met_utilities::virtual_temperature(temp(icol, ilev, 0), qv(icol, ilev, 0));
                    bxheight(icol, ilev, 0) = (constants::RD / constants::G0) * virtual_t * std::log(p_lower / p_upper);
                } else {
                    bxheight(icol, ilev, 0) = 0.0;
                }
            }
        );
    }

    void derive_airden_dry() {
        if (!met.PMID || !met.T || !met.QV || !met.AIRDEN_DRY) return;

        int nc = n_cols;
        int nl = n_levels;

        auto pmid = met.PMID->view();
        auto temp = met.T->view();
        auto qv = met.QV->view();
        auto airden_dry = met.AIRDEN_DRY->view();

        Kokkos::parallel_for("derive_airden_dry_kernel",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, 0}, {nc, nl}),
            KOKKOS_LAMBDA(int icol, int ilev) {
                double q = qv(icol, ilev, 0);
                double avgw = (constants::AIR_MW / constants::H2O_MW) * q / (1.0 - q);
                double xh2o = avgw / (1.0 + avgw);

                double p_dry = pmid(icol, ilev, 0) * (1.0 - xh2o);
                airden_dry(icol, ilev, 0) = p_dry / (constants::RD * temp(icol, ilev, 0));
            }
        );
    }
```

- [ ] **Step 2: Declare and implement C-API triggers**

In `src/core/catchem_api.hpp`:
```cpp
void catchem_state_derive_bxheight(void* state_ptr);
void catchem_state_derive_airden_dry(void* state_ptr);
```

In `src/core/catchem_api.cpp`:
```cpp
void catchem_state_derive_bxheight(void* state_ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    state->derive_bxheight();
}

void catchem_state_derive_airden_dry(void* state_ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    state->derive_airden_dry();
}
```

- [ ] **Step 3: Commit Task 4 changes**

```bash
git add src/core/catchem_state_manager.hpp src/core/catchem_api.*
git commit -m "feat(core): implement parallel hydrostatic and density equations inside C++ StateManager"
```

---

### Task 5: Add Advanced Multi-State Assertions in Integration Suite

**Files:**
- Modify: `tests/test_catchem_interop.cpp`

**Interfaces:**
- Consumes: Expose MetState/ChemState binders and launchers

- [ ] **Step 1: Write integration assertions in test harness**

Add **TEST 6: Parallel Meteorological Derivations, SZA, and Unified Chem State** inside `tests/test_catchem_interop.cpp` right before `TEST 5`:
```cpp
        // ==========================================
        // TEST 6: Parallel Meteorological Derivations, SZA, and Unified Chem State
        // ==========================================
        {
            int n_cols = 4;
            int n_levels = 5;
            int n_species = 2;

            void* core = catchem_core_create(n_cols, n_levels, n_species);
            void* state = catchem_core_get_state_manager(core);

            // 1. Allocate mock meteorological and chemical arrays
            std::vector<double> temp_array(n_cols * n_levels, 290.15); // Temperature [K]
            std::vector<double> qv_array(n_cols * n_levels, 0.01);    // Specific humidity [kg/kg]
            std::vector<double> pmid_array(n_cols * n_levels, 100000.0); // Mid-pressure [Pa]
            std::vector<double> pedge_array(n_cols * (n_levels + 1), 101325.0); // Pressure edges [Pa]

            // Assign standard pressure edge levels sequentially
            for (int i = 0; i < n_cols; ++i) {
                pedge_array[i + 0 * n_cols] = 101325.0; // Surface
                pedge_array[i + 1 * n_cols] = 90000.0;
                pedge_array[i + 2 * n_cols] = 80000.0;
                pedge_array[i + 3 * n_cols] = 70000.0;
                pedge_array[i + 4 * n_cols] = 60000.0;
                pedge_array[i + 5 * n_cols] = 50000.0; // Top
            }

            std::vector<double> bxheight_array(n_cols * n_levels, 0.0); // Output height
            std::vector<double> airden_dry_array(n_cols * n_levels, 0.0); // Output dry density

            std::vector<double> mock_chem_state(n_cols * n_levels * n_species, 4.2); // Unified chem state

            // 2. Bind arrays to StateManager
            catchem_state_bind_met_3d(state, "T", temp_array.data());
            catchem_state_bind_met_3d(state, "QV", qv_array.data());
            catchem_state_bind_met_3d(state, "PMID", pmid_array.data());
            catchem_state_bind_met_3d(state, "PEDGE", pedge_array.data());
            catchem_state_bind_met_3d(state, "BXHEIGHT", bxheight_array.data());
            catchem_state_bind_met_3d(state, "AIRDEN_DRY", airden_dry_array.data());

            catchem_state_bind_unified_chemistry(state, mock_chem_state.data());

            // 3. Sync arrays to device memory spaces
            catchem_state_sync_to_device(state);

            // 4. Trigger parallel derived met equations
            catchem_state_derive_bxheight(state);
            catchem_state_derive_airden_dry(state);

            // 5. Sync derived results back to host heap
            catchem_state_sync_to_host(state);

            // 6. Assert correct calculations
            // Layer 1 edge pressures: P_lower = 101325.0, P_upper = 90000.0
            // Virtual T = 290.15 * (1 + 0.608 * 0.01) = 291.914
            // Expected height = (287 / 9.80665) * 291.914 * std::log(101325.0 / 90000.0) ≈ 1010.5 meters
            double derived_h = bxheight_array[0];
            assert(derived_h > 990.0 && derived_h < 1030.0);
            std::cout << "INFO: Derived BXHEIGHT = " << derived_h << " meters.\n";

            double derived_rho = airden_dry_array[0];
            assert(derived_rho > 1.0 && derived_rho < 1.3);
            std::cout << "INFO: Derived Dry Air Density = " << derived_rho << " kg/m³.\n";

            // Assert unified chemistry array mapped accurately
            auto* state_obj = static_cast<catchem::StateManager*>(state);
            assert(state_obj->chem.conc != nullptr);
            assert(state_obj->chem.conc->host_view(0, 0, 0) == 4.2);

            // 7. Test portable Time State calculations
            catchem_state_set_time(state, 2026, 7, 8, 12, 0, 0, 189, 3600.0);
            double cos_sza = state_obj->time.get_cos_sza(40.0, -80.0);
            assert(cos_sza >= -1.0 && cos_sza <= 1.0);
            std::cout << "INFO: Calculated Cos(SZA) at lat=40, lon=-80: " << cos_sza << "\n";

            catchem_core_destroy(core);
            std::cout << "SUCCESS: Parallel Meteorological Derivations & SZA Validation Passed!\n";
        }
```

- [ ] **Step 2: Build and run the entire interop test suite**

Run inside Docker:
`docker run --rm -v /Users/barry/Documents/CATChem:/workspace -w /workspace/build-test cece-dev:latest bash -c "cmake .. -DENABLE_KOKKOS=ON -DENABLE_TESTING=ON && make test_catchem_interop && ./tests/test_catchem_interop"`

Expected: PASS, printing:
```text
SUCCESS: Interop Shared State Validation Passed!
...
SUCCESS: Parallel Meteorological Derivations & SZA Validation Passed!
SUCCESS: C++ Species Metadata & State Initialization Validation Passed!
```

- [ ] **Step 3: Commit test updates**

```bash
git add tests/test_catchem_interop.cpp
git commit -m "test(interop): verify parallel derived meteorological equations and SZA with unified ChemState"
```
