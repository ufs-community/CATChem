# Core Modernization: ConfigManager and GridManager Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port CATChem's `ConfigManager` and `GridManager` to C++20, and expose their states via the C-API to preserve legacy Fortran interoperability while paving the way for native C++ 3D execution.

**Architecture:** We will create `catchem_config_manager` using `yaml-cpp` to load simulation properties. We will create `catchem_grid_manager` and `catchem_grid_geometry` to manage dimensionality and geography. `catchem::Core` will own these managers. C-API bindings will be created to allow Fortran code to query configurations and grid dimensions directly.

**Tech Stack:** C++20, yaml-cpp, CMake, Kokkos

## Global Constraints

- Target C++20 utilizing the Kokkos backport of mdspan (`std::experimental::mdspan`), avoiding direct dependency on C++23 `<mdspan>`.
- Retain Fortran column-major storage layout (`Kokkos::LayoutLeft`) across the pointer boundary to achieve zero-copy execution on CPU targets.

---

### Task 1: Create `catchem::ConfigManager`

**Files:**
- Create: `src/core/catchem_config_manager.hpp`
- Create: `src/core/catchem_config_manager.cpp`
- Modify: `src/core/CMakeLists.txt`

**Interfaces:**
- Produces: `catchem::ConfigData`, `catchem::ConfigManager`
- Produces: `ConfigManager::load_from_file(const std::string&)`

- [ ] **Step 1: Write the failing compilation test**

```cpp
// tests/test_config_manager_compilation.cpp
#include "catchem_config_manager.hpp"

int main() {
    catchem::ConfigManager config_mgr;
    config_mgr.load_from_file("tests/CATChem_new_config.yml");
    return 0;
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++20 -c tests/test_config_manager_compilation.cpp`
Expected: FAIL (missing header/symbol)

- [ ] **Step 3: Write implementation header**

```cpp
// src/core/catchem_config_manager.hpp
#pragma once
#include <string>
#include <yaml-cpp/yaml.h>

namespace catchem {

struct RuntimeConfig {
    int nx = 1;
    int ny = 1;
    int nz = 1;
    double dt = 3600.0;
    int nsteps = 1;
};

struct ConfigData {
    RuntimeConfig runtime;
    // We can add FilePathConfig, etc. here later
};

class ConfigManager {
public:
    ConfigData data;
    bool is_loaded = false;

    ConfigManager() = default;
    void load_from_file(const std::string& filename);
};

} // namespace catchem
```

- [ ] **Step 4: Write implementation cpp**

```cpp
// src/core/catchem_config_manager.cpp
#include "catchem_config_manager.hpp"
#include <stdexcept>

namespace catchem {

void ConfigManager::load_from_file(const std::string& filename) {
    try {
        YAML::Node config = YAML::LoadFile(filename);
        if (config["simulation"]) {
            auto sim = config["simulation"];
            if (sim["nx"]) data.runtime.nx = sim["nx"].as<int>();
            if (sim["ny"]) data.runtime.ny = sim["ny"].as<int>();
            if (sim["nz"]) data.runtime.nz = sim["nz"].as<int>();
            if (sim["timestep"]) data.runtime.dt = sim["timestep"].as<double>();
            if (sim["nsteps"]) data.runtime.nsteps = sim["nsteps"].as<int>();
        }
        is_loaded = true;
    } catch (const YAML::Exception& e) {
        throw std::runtime_error("Failed to parse config file: " + std::string(e.what()));
    }
}

} // namespace catchem
```

- [ ] **Step 5: Register file in CMake**

Modify `src/core/CMakeLists.txt` to add `catchem_config_manager.cpp` to the `_cpp_core_srcs` list.

- [ ] **Step 6: Remove dummy test and commit**

```bash
rm tests/test_config_manager_compilation.cpp
git add src/core/catchem_config_manager.* src/core/CMakeLists.txt
git commit -m "feat(core): implement C++ ConfigManager using yaml-cpp"
```

---

### Task 2: Create `catchem::GridManager`

**Files:**
- Create: `src/core/catchem_grid_manager.hpp`
- Create: `src/core/catchem_grid_manager.cpp`
- Modify: `src/core/CMakeLists.txt`

**Interfaces:**
- Produces: `catchem::GridGeometry`, `catchem::GridManager`

- [ ] **Step 1: Write the failing compilation test**

```cpp
// tests/test_grid_manager_compilation.cpp
#include "catchem_grid_manager.hpp"

int main() {
    catchem::GridManager grid_mgr(10, 10, 50);
    return 0;
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++20 -c tests/test_grid_manager_compilation.cpp`
Expected: FAIL

- [ ] **Step 3: Write implementation header**

```cpp
// src/core/catchem_grid_manager.hpp
#pragma once
#include <memory>
#include "catchem_interop_field.hpp"

namespace catchem {

struct GridGeometry {
    int nx = 1;
    int ny = 1;
    int nz = 1;

    std::shared_ptr<InteropField<double, 2>> lat;
    std::shared_ptr<InteropField<double, 2>> lon;
    std::shared_ptr<InteropField<double, 2>> grid_area;
    // dz and z_levels can be added as 1D fields if needed globally
};

class GridManager {
public:
    GridGeometry geometry;
    bool is_initialized = false;

    GridManager(int nx, int ny, int nz);

    // Bindings to support Fortran Interop arrays if allocated externally
    void bind_lat(double* ptr);
    void bind_lon(double* ptr);
    void bind_area(double* ptr);
};

} // namespace catchem
```

- [ ] **Step 4: Write implementation cpp**

```cpp
// src/core/catchem_grid_manager.cpp
#include "catchem_grid_manager.hpp"
#include <vector>

namespace catchem {

GridManager::GridManager(int nx, int ny, int nz) {
    geometry.nx = nx;
    geometry.ny = ny;
    geometry.nz = nz;
    is_initialized = true;
}

void GridManager::bind_lat(double* ptr) {
    geometry.lat = std::make_shared<InteropField<double, 2>>(ptr, std::vector<int>{geometry.nx, geometry.ny});
}

void GridManager::bind_lon(double* ptr) {
    geometry.lon = std::make_shared<InteropField<double, 2>>(ptr, std::vector<int>{geometry.nx, geometry.ny});
}

void GridManager::bind_area(double* ptr) {
    geometry.grid_area = std::make_shared<InteropField<double, 2>>(ptr, std::vector<int>{geometry.nx, geometry.ny});
}

} // namespace catchem
```

- [ ] **Step 5: Register file in CMake**

Modify `src/core/CMakeLists.txt` to add `catchem_grid_manager.cpp` to the `_cpp_core_srcs` list.

- [ ] **Step 6: Remove dummy test and commit**

```bash
rm tests/test_grid_manager_compilation.cpp
git add src/core/catchem_grid_manager.* src/core/CMakeLists.txt
git commit -m "feat(core): implement C++ GridManager and GridGeometry"
```

---

### Task 3: Integrate Managers into `catchem::Core`

**Files:**
- Modify: `src/core/catchem_core.hpp`
- Modify: `src/core/catchem_core.cpp`

**Interfaces:**
- Consumes: `catchem::ConfigManager`, `catchem::GridManager`
- Modifies: `catchem::Core` constructor signature to accept config file instead of raw dimensions.

- [ ] **Step 1: Write integration compilation test**

```cpp
// tests/test_core_integration.cpp
#include "catchem_core.hpp"

int main() {
    catchem::Core core("tests/CATChem_new_config.yml");
    return 0;
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++20 -c tests/test_core_integration.cpp`
Expected: FAIL (Core constructor signature mismatch)

- [ ] **Step 3: Modify `catchem_core.hpp`**

```cpp
#pragma once
#include <memory>
#include <vector>
#include <string>
#include "catchem_state_manager.hpp"
#include "catchem_process_interface.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_config_manager.hpp"
#include "catchem_grid_manager.hpp"

namespace catchem {

class Core {
private:
    std::shared_ptr<ConfigManager> config_mgr;
    std::shared_ptr<GridManager> grid_mgr;
    std::shared_ptr<StateManager> state_mgr;
    std::shared_ptr<DiagnosticManager> diag_mgr;
    std::vector<std::shared_ptr<ProcessInterface>> processes;
public:
    Core(const std::string& config_file);
    std::shared_ptr<ConfigManager> get_config_manager();
    std::shared_ptr<GridManager> get_grid_manager();
    std::shared_ptr<StateManager> get_state_manager();
    std::shared_ptr<DiagnosticManager> get_diagnostic_manager();
    void add_process(std::shared_ptr<ProcessInterface> process);
    void run_timestep(); // We will pull dt from config
};

} // namespace catchem
```

- [ ] **Step 4: Modify `catchem_core.cpp`**

```cpp
#include "catchem_core.hpp"

namespace catchem {

Core::Core(const std::string& config_file) {
    config_mgr = std::make_shared<ConfigManager>();
    config_mgr->load_from_file(config_file);

    int nx = config_mgr->data.runtime.nx;
    int ny = config_mgr->data.runtime.ny;
    int nz = config_mgr->data.runtime.nz;

    grid_mgr = std::make_shared<GridManager>(nx, ny, nz);

    // For n_species, we either load it here or default to a safe value until ChemState loads it.
    // Assuming 50 default, or read from config if we add it to ConfigManager later.
    state_mgr = std::make_shared<StateManager>(nx, ny, nz, 50);
    diag_mgr = std::make_shared<DiagnosticManager>();
}

std::shared_ptr<ConfigManager> Core::get_config_manager() {
    return config_mgr;
}

std::shared_ptr<GridManager> Core::get_grid_manager() {
    return grid_mgr;
}

std::shared_ptr<StateManager> Core::get_state_manager() {
    return state_mgr;
}

std::shared_ptr<DiagnosticManager> Core::get_diagnostic_manager() {
    return diag_mgr;
}

void Core::add_process(std::shared_ptr<ProcessInterface> process) {
    processes.push_back(process);
}

void Core::run_timestep() {
    double dt = config_mgr->data.runtime.dt;
    state_mgr->sync_to_device();
    for (auto& process : processes) {
        process->run(state_mgr);
    }
    state_mgr->sync_to_host();
    diag_mgr->sync_to_host();
}

} // namespace catchem
```

**Note:** Ensure `StateManager` constructor in `catchem_state_manager.hpp` is updated if necessary to handle the 4-arg signature used above, or keep it 3-arg and set `n_species` separately. *Fix inline if needed.*

- [ ] **Step 5: Modify `catchem_api.cpp` to use new `Core` signature**

Modify `void* catchem_core_create` to accept a string instead of ints. Wait, changing C-API breaks Fortran tests immediately.
Let's add a new C-API function instead of replacing it, or keep the old one but create a dummy config inside it for backward compatibility.
Let's modify `catchem_api.cpp` to map `catchem_core_create` to a mock config or create a new `catchem_core_create_from_config(const char* config_file)`.

```cpp
// In catchem_api.hpp:
void* catchem_core_create_from_config(const char* config_file);

// In catchem_api.cpp:
void* catchem_core_create_from_config(const char* config_file) {
    return static_cast<void*>(new catchem::Core(config_file));
}

// Update the old create method to use a generic config file string if needed, or leave it if you modified Core to have overloaded constructors.
```
*To be safe, add an overloaded constructor `Core(int nc, int nl, int ns)` to `Core` that initializes default `ConfigManager` and `GridManager` so existing Fortran C-API tests don't break.*

- [ ] **Step 6: Remove test and commit**

```bash
rm tests/test_core_integration.cpp
git add src/core/catchem_core.* src/core/catchem_api.*
git commit -m "feat(core): integrate ConfigManager and GridManager into CATChem Core"
```

---

### Task 4: Expose C-API Endpoints for Fortran Interoperability

**Files:**
- Modify: `src/core/catchem_api.hpp`
- Modify: `src/core/catchem_api.cpp`

**Interfaces:**
- Produces: `catchem_get_grid_dimensions`, `catchem_get_config_timestep`

- [ ] **Step 1: Declare C-API endpoints**

In `src/core/catchem_api.hpp`:
```cpp
void catchem_get_grid_dimensions(void* core_ptr, int* nx, int* ny, int* nz);
double catchem_get_config_timestep(void* core_ptr);
```

- [ ] **Step 2: Implement endpoints**

In `src/core/catchem_api.cpp`:
```cpp
void catchem_get_grid_dimensions(void* core_ptr, int* nx, int* ny, int* nz) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    auto grid = core->get_grid_manager();
    *nx = grid->geometry.nx;
    *ny = grid->geometry.ny;
    *nz = grid->geometry.nz;
}

double catchem_get_config_timestep(void* core_ptr) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return core->get_config_manager()->data.runtime.dt;
}
```

- [ ] **Step 3: Commit**

```bash
git add src/core/catchem_api.*
git commit -m "feat(api): expose Grid and Config queries to C-API for Fortran legacy interop"
```

---

### Task 5: Integration Tests

**Files:**
- Modify: `tests/test_catchem_interop.cpp`

**Interfaces:**
- Consumes: `catchem_core_create_from_config`, `catchem_get_grid_dimensions`, `catchem_get_config_timestep`

- [ ] **Step 1: Write integration assertions**

Add **TEST 7: C++ Config and Grid Initialization** inside `tests/test_catchem_interop.cpp`:
```cpp
        // ==========================================
        // TEST 7: C++ Config and Grid Initialization
        // ==========================================
        {
            void* core = catchem_core_create_from_config("tests/CATChem_new_config.yml");

            int nx, ny, nz;
            catchem_get_grid_dimensions(core, &nx, &ny, &nz);
            assert(nx > 0);
            assert(nz > 0);
            std::cout << "INFO: Loaded grid dimensions from config: " << nx << "x" << ny << "x" << nz << "\n";

            double dt = catchem_get_config_timestep(core);
            assert(dt > 0.0);
            std::cout << "INFO: Loaded timestep from config: " << dt << " s\n";

            catchem_core_destroy(core);
            std::cout << "SUCCESS: C++ Config & Grid Validation Passed!\n";
        }
```

- [ ] **Step 2: Build and run the interop test suite**

Run inside Docker:
`docker run --rm -v $(pwd):/workspace -w /workspace/build-test cece-dev:latest bash -c "cmake .. -DENABLE_KOKKOS=ON -DENABLE_TESTING=ON && make test_catchem_interop && ./tests/test_catchem_interop"`

Expected: PASS, printing:
```text
SUCCESS: C++ Config & Grid Validation Passed!
```

- [ ] **Step 3: Commit test updates**

```bash
git add tests/test_catchem_interop.cpp
git commit -m "test(interop): verify ConfigManager and GridManager C-API bindings"
```
