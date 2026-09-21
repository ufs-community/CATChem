# C++ State Initialization and Species Metadata Management Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Elevate chemical species metadata ownership, runtime configuration loading, and indexing to C++20/Kokkos, standardizing configuration via YAML parsing using `yaml-cpp`.

**Architecture:** Create `SpeciesMetadata` struct. Extend `StateManager` with a species list, string-to-index mappings, category filters, and a YAML-based loader using `yaml-cpp`. Expose querying endpoints via the C-API boundary.

**Tech Stack:** C++20, yaml-cpp, CMake, Kokkos

## Global Constraints

- Target C++20 utilizing the Kokkos backport of mdspan (`std::experimental::mdspan`), avoiding direct dependency on C++23 `<mdspan>`.
- Retain Fortran column-major storage layout (`Kokkos::LayoutLeft`) across the pointer boundary to achieve zero-copy execution on CPU targets.
- All modifications must remain unstaged and uncommitted until instructed otherwise.

---

### Task 1: Create `catchem::SpeciesMetadata` Struct

**Files:**
- Create: `src/core/catchem_species_metadata.hpp`

**Interfaces:**
- Produces: `catchem::SpeciesMetadata`

- [ ] **Step 1: Create header file `src/core/catchem_species_metadata.hpp`**

Write the complete struct declaration:
```cpp
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

- [ ] **Step 2: Commit file**

```bash
git add src/core/catchem_species_metadata.hpp
git commit -m "feat(core): add C++ representation of chemical species metadata"
```

---

### Task 2: Configure CMake and Link `yaml-cpp`

**Files:**
- Modify: `src/core/CMakeLists.txt`

**Interfaces:**
- Consumes: `yaml-cpp` CMake target

- [ ] **Step 1: Link yaml-cpp to CATChem_core_cpp**

Modify `src/core/CMakeLists.txt` to link `yaml-cpp` to `CATChem_core_cpp`:
```cmake
# Replace the target_link_libraries(CATChem_core_cpp PUBLIC Kokkos::kokkos) line with:
target_link_libraries(CATChem_core_cpp PUBLIC Kokkos::kokkos yaml-cpp)
```

- [ ] **Step 2: Commit CMake change**

```bash
git add src/core/CMakeLists.txt
git commit -m "build(cmake): link yaml-cpp target to C++ core library"
```

---

### Task 3: Extend `catchem::StateManager` with Species Loader

**Files:**
- Modify: `src/core/catchem_state_manager.hpp`

**Interfaces:**
- Consumes: `catchem::SpeciesMetadata`
- Produces: `StateManager::species_list`, `StateManager::load_species_config`

- [ ] **Step 1: Write temporary compilation test for StateManager loading**

Modify `tests/test_catchem_interop.cpp` to include the headers and try calling the new method:
```cpp
// Add to the top of tests/test_catchem_interop.cpp:
#include "catchem_state_manager.hpp"

// Inside main(), add:
{
    catchem::StateManager state(4, 5, 2);
    // This will initially fail to compile as load_species_config doesn't exist yet
    // state.load_species_config("tests/test_species.yml");
}
```

- [ ] **Step 2: Run build to verify compilation failure or missing symbol**

Expected: Compilation error or linker error regarding missing `load_species_config` or `<yaml-cpp/yaml.h>` when uncommented.

- [ ] **Step 3: Implement load_species_config in StateManager**

Include headers and define fields/loaders in `src/core/catchem_state_manager.hpp`:
```cpp
// Replace the top of src/core/catchem_state_manager.hpp with:
#pragma once
#include <unordered_map>
#include <string>
#include <memory>
#include <vector>
#include <yaml-cpp/yaml.h>
#include "catchem_interop_field.hpp"
#include "catchem_species_metadata.hpp"

namespace catchem {

class StateManager {
public:
    int n_cols;
    int n_levels;
    int n_species;

    std::unordered_map<std::string, std::shared_ptr<InteropField<double, 1>>> fields_1d;
    std::unordered_map<std::string, std::shared_ptr<InteropField<double, 2>>> fields_2d;
    std::unordered_map<std::string, std::shared_ptr<InteropField<double, 3>>> fields_3d;

    // Species metadata structures
    std::vector<SpeciesMetadata> species_list;
    std::unordered_map<std::string, int> species_name_to_index; // 0-based index

    // Category lists (0-based offsets)
    std::vector<int> gas_indices;
    std::vector<int> aerosol_indices;
    std::vector<int> tracer_indices;
    std::vector<int> advected_indices;
    std::vector<int> drydep_indices;
    std::vector<int> wetdep_indices;
    std::vector<int> photolysis_indices;
    std::vector<int> dust_indices;
    std::vector<int> seasalt_indices;

    StateManager(int nc, int nl, int ns) : n_cols(nc), n_levels(nl), n_species(ns) {}

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
    // ... rest of StateManager methods (bind_field_1d, bind_field_2d, bind_field_3d, etc.) ...
```

- [ ] **Step 4: Verify test compiles successfully**

Expected: Successful compilation when uncommenting test loading code.

- [ ] **Step 5: Commit changes**

```bash
git add src/core/catchem_state_manager.hpp
git commit -m "feat(core): extend StateManager with species metadata loader using yaml-cpp"
```

---

### Task 4: Expose Species Queries inside C-API Boundary

**Files:**
- Modify: `src/core/catchem_api.hpp`
- Modify: `src/core/catchem_api.cpp`

**Interfaces:**
- Produces: `catchem_state_load_species_config`, `catchem_state_get_species_count`, `catchem_state_get_species_index`, `catchem_state_get_gas_species_count`, `catchem_state_get_gas_indices`, `catchem_state_get_aerosol_species_count`, `catchem_state_get_aerosol_indices`, `catchem_state_get_species_mw`, `catchem_state_is_species_gas`, `catchem_state_is_species_aerosol`

- [ ] **Step 1: Declare query endpoints in API header**

Modify `src/core/catchem_api.hpp`:
```cpp
// Append these declarations before the final #ifdef __cplusplus block:
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

- [ ] **Step 2: Implement query endpoints in API implementation**

Modify `src/core/catchem_api.cpp`:
```cpp
// Append these C-API implementations to the extern "C" block in src/core/catchem_api.cpp:

void catchem_state_load_species_config(void* state_ptr, const char* filename) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    state->load_species_config(filename);
}

int catchem_state_get_species_count(void* state_ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    return static_cast<int>(state->species_list.size());
}

int catchem_state_get_species_index(void* state_ptr, const char* name) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    auto it = state->species_name_to_index.find(name);
    if (it != state->species_name_to_index.end()) {
        return it->second + 1; // Translate 0-based C++ index to 1-based Fortran index
    }
    return -1;
}

int catchem_state_get_gas_species_count(void* state_ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    return static_cast<int>(state->gas_indices.size());
}

void catchem_state_get_gas_indices(void* state_ptr, int* indices_out) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    for (size_t i = 0; i < state->gas_indices.size(); ++i) {
        indices_out[i] = state->gas_indices[i] + 1; // 1-based
    }
}

int catchem_state_get_aerosol_species_count(void* state_ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    return static_cast<int>(state->aerosol_indices.size());
}

void catchem_state_get_aerosol_indices(void* state_ptr, int* indices_out) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    for (size_t i = 0; i < state->aerosol_indices.size(); ++i) {
        indices_out[i] = state->aerosol_indices[i] + 1; // 1-based
    }
}

double catchem_state_get_species_mw(void* state_ptr, int index) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    int idx_0 = index - 1; // 1-based to 0-based
    if (idx_0 >= 0 && idx_0 < static_cast<int>(state->species_list.size())) {
        return state->species_list[idx_0].mw_g;
    }
    return 0.0;
}

int catchem_state_is_species_gas(void* state_ptr, int index) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    int idx_0 = index - 1;
    if (idx_0 >= 0 && idx_0 < static_cast<int>(state->species_list.size())) {
        return state->species_list[idx_0].is_gas ? 1 : 0;
    }
    return 0;
}

int catchem_state_is_species_aerosol(void* state_ptr, int index) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    int idx_0 = index - 1;
    if (idx_0 >= 0 && idx_0 < static_cast<int>(state->species_list.size())) {
        return state->species_list[idx_0].is_aerosol ? 1 : 0;
    }
    return 0;
}
```

- [ ] **Step 3: Commit C-API updates**

```bash
git add src/core/catchem_api.hpp src/core/catchem_api.cpp
git commit -m "feat(api): expose complete species metadata query endpoints in C-API layer"
```

---

### Task 5: Add Metadata Tests in Integration Suite

**Files:**
- Modify: `tests/test_catchem_interop.cpp`

**Interfaces:**
- Consumes: C-API metadata querying endpoints

- [ ] **Step 1: Write verification integration assertions**

Add **TEST 5: C++ Species Metadata & State Initialization** inside `tests/test_catchem_interop.cpp`:
```cpp
// Append this test block at the end of main() right before Kokkos::finalize() in tests/test_catchem_interop.cpp:

        // ==========================================
        // TEST 5: C++ Species Metadata & State Initialization
        // ==========================================
        {
            int n_cols = 4;
            int n_levels = 5;
            int n_species = 2;

            void* core = catchem_core_create(n_cols, n_levels, n_species);
            void* state = catchem_core_get_state_manager(core);

            // 1. Load species config from tests/test_species.yml
            catchem_state_load_species_config(state, "tests/test_species.yml");

            // 2. Validate species counts and offsets
            int count = catchem_state_get_species_count(state);
            assert(count > 0);
            std::cout << "INFO: Loaded " << count << " species in integration test.\n";

            // 3. Translate species names to 1-based indices
            int idx_so2 = catchem_state_get_species_index(state, "so2");
            int idx_so4 = catchem_state_get_species_index(state, "so4");
            assert(idx_so2 != -1);
            assert(idx_so4 != -1);

            // 4. Validate physical properties of species
            double mw_so2 = catchem_state_get_species_mw(state, idx_so2);
            assert(mw_so2 == 64.04);

            int is_gas_so2 = catchem_state_is_species_gas(state, idx_so2);
            int is_aero_so2 = catchem_state_is_species_aerosol(state, idx_so2);
            assert(is_gas_so2 == 1);
            assert(is_aero_so2 == 0);

            int is_gas_so4 = catchem_state_is_species_gas(state, idx_so4);
            int is_aero_so4 = catchem_state_is_species_aerosol(state, idx_so4);
            assert(is_gas_so4 == 0);
            assert(is_aero_so4 == 1);

            // 5. Validate category lists (gas / aerosol)
            int gas_count = catchem_state_get_gas_species_count(state);
            assert(gas_count > 0);
            std::vector<int> gas_indices(gas_count);
            catchem_state_get_gas_indices(state, gas_indices.data());

            // Ensure so2 index is present in gas_indices
            bool found_so2 = false;
            for (int idx : gas_indices) {
                if (idx == idx_so2) found_so2 = true;
            }
            assert(found_so2);

            catchem_core_destroy(core);
            std::cout << "SUCCESS: C++ Species Metadata & State Initialization Validation Passed!\n";
        }
```

- [ ] **Step 2: Build and run the entire interop test suite**

Run inside Docker:
`docker run --rm -v /Users/barry/Documents/CATChem:/workspace -w /workspace/build-test cece-dev:latest bash -c "apt-get update -y && apt-get install -y python3 && cmake .. -DENABLE_KOKKOS=ON -DENABLE_TESTING=ON && make test_catchem_interop && ./tests/test_catchem_interop"`

Expected: PASS, printing:
```text
SUCCESS: Interop Shared State Validation Passed!
...
SUCCESS: C++ Species Metadata & State Initialization Validation Passed!
```

- [ ] **Step 3: Commit test updates**

```bash
git add tests/test_catchem_interop.cpp
git commit -m "test(interop): add integration test for C++ species metadata loading and validation"
```
