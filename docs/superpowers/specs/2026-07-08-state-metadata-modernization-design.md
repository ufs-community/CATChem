# Technical Design Specification: Phase 4 — C++ State Initialization and Species Metadata Management

**Date:** July 8, 2026  
**Status:** Draft  
**Topic:** Elevating chemical species metadata management, runtime config loading, and indexing to C++20/Kokkos, standardizing configuration via YAML parsing using `yaml-cpp`.

---

## 1. Context & Objectives

In legacy CATChem, species metadata (names, flags like `is_gas` / `is_aerosol`, molecular weights, deposition parameters) were parsed inside Fortran (`ConfigManager_Mod.F90`) and stored across several separate allocated arrays in `ChemStateType` (`chemstate_mod.F90`).

To fully modernize CATChem and enable future Kokkos-based physics/chemistry schemes, we must elevate species metadata ownership to the C++ core. This allows:
1. **Unified Configuration Parsing:** Directly parse `CATChem_species.yml` using C++ `yaml-cpp`.
2. **GPU-Friendly Indexing:** Pre-filter active category indices (e.g., gas, aerosols, tracers, drydep, wetdep) in C++ vectors, facilitating direct capture in parallel Kokkos functors.
3. **Single Source of Truth:** Expose a clean, ISO_C_BINDING-compatible C-API so legacy Fortran processes can query species counts, indices, category filters, and numerical properties (like molecular weights) directly from the C++ `StateManager`.

---

## 2. Component Architecture & Class Designs

### 2.1 Species Metadata representation

A new struct `catchem::SpeciesMetadata` represents the complete set of chemical and physical parameters of a single species.

```cpp
// src/core/catchem_species_metadata.hpp
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
    double wd_retfactor = 0.0;
    bool wd_LiqAndGas = false;
    double wd_convfacI2G = 0.0;
    std::vector<double> wd_rainouteff;
    std::string mie_name;
};

} // namespace catchem
```

---

### 2.2 Extending `catchem::StateManager`

We integrate species database parsing and category mappings directly into `catchem::StateManager`:

```cpp
// src/core/catchem_state_manager.hpp additions
#include "catchem_species_metadata.hpp"
#include <yaml-cpp/yaml.h>

namespace catchem {

class StateManager {
public:
    // ... existing fields ...

    // Species database
    std::vector<SpeciesMetadata> species_list;
    std::unordered_map<std::string, int> species_name_to_index; // 0-based indexing in C++

    // Pre-filtered category indices (extremely helpful for Kokkos loop boundaries)
    std::vector<int> gas_indices;
    std::vector<int> aerosol_indices;
    std::vector<int> tracer_indices;
    std::vector<int> advected_indices;
    std::vector<int> drydep_indices;
    std::vector<int> wetdep_indices;
    std::vector<int> photolysis_indices;
    std::vector<int> dust_indices;
    std::vector<int> seasalt_indices;

    // Loads complete species YAML file and populates the metadata lists
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

---

## 3. C-API Boundary Enhancements

We expose new query functions at the C-API boundary (`src/core/catchem_api.hpp` and `catchem_api.cpp`) to let the Fortran side query metadata and keep indices fully synchronized.

```cpp
// src/core/catchem_api.hpp additions
void catchem_state_load_species_config(void* state_ptr, const char* filename);
int catchem_state_get_species_count(void* state_ptr);
int catchem_state_get_species_index(void* state_ptr, const char* name); // returns 1-based index matching Fortran, or -1 if not found

// Categorized counts and list getters
int catchem_state_get_gas_species_count(void* state_ptr);
void catchem_state_get_gas_indices(void* state_ptr, int* indices_out); // populates 1-based indices
int catchem_state_get_aerosol_species_count(void* state_ptr);
void catchem_state_get_aerosol_indices(void* state_ptr, int* indices_out);

// Individual property getters (by 1-based index)
double catchem_state_get_species_mw(void* state_ptr, int index);
int catchem_state_is_species_gas(void* state_ptr, int index);
int catchem_state_is_species_aerosol(void* state_ptr, int index);
```

---

## 4. CMake Target Configuration

To use `yaml-cpp` within the core C++ target, we explicitly link it to `CATChem_core_cpp` in `src/core/CMakeLists.txt`:

```cmake
target_link_libraries(CATChem_core_cpp PUBLIC Kokkos::kokkos yaml-cpp)
```

---

## 5. Integration Verification & Testing Strategy

To verify this implementation, we will append a dedicated test block `TEST 5: C++ Species Metadata & State Initialization` inside `tests/test_catchem_interop.cpp`. This test will:
1. Initialize the C++ core and retrieve the `StateManager`.
2. Load species metadata using `catchem_state_load_species_config(state, "tests/test_species.yml")`.
3. Assert correct parsing of basic properties for key species (e.g., verifying `so2` is a gas with `mw_g == 64.04`, and `so4` is an aerosol).
4. Assert correct parsing of subcategory counts and indices mapping (e.g. `gas_indices`).

---

## 6. Self-Review Check
* **Placeholder check:** No `TBD` or `TODO` left in design.
* **Consistency:** The C++ index retrieval utilizes 1-based translation to seamlessly align with legacy Fortran without requiring rewriting surrounding array indexing.
* **Scope:** Focused purely on metadata mapping and C++ standard configuration, avoiding unrelated science scheme modifications.
