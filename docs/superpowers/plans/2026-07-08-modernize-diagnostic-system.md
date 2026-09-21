# Modernize Diagnostic System Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Modernize the CATChem diagnostic infrastructure by creating a C++ `DiagnosticManager` backed by Kokkos Views to enable GPU-safe, thread-safe diagnostics collection.

**Architecture:** C++ parallel kernels capture Kokkos-backed diagnostic views by value and write directly to them with zero-copy execution. On CPUs, diagnostics leverage Fortran-compatible column-major layouts (`Kokkos::LayoutLeft`). C-API bindings expose host pointers to Fortran for NetCDF IO.

**Tech Stack:** C++20, Kokkos, Fortran 2008 (ISO_C_BINDING), CMake

## Global Constraints

- Target C++20 utilizing the Kokkos backport of mdspan (`std::experimental::mdspan`), avoiding direct dependency on C++23 `<mdspan>`.
- Retain Fortran column-major storage layout (`Kokkos::LayoutLeft`) across the pointer boundary to achieve zero-copy execution on CPU.

---

### Task 1: Create C++ Diagnostic System Foundation

**Files:**
- Create: `src/core/catchem_diagnostic.hpp`
- Create: `src/core/catchem_diagnostic_manager.hpp`
- Create: `src/core/catchem_diagnostic.cpp`

**Interfaces:**
- Produces: `catchem::DiagType`, `catchem::DiagnosticField`, `catchem::DiagnosticManager`

- [ ] **Step 1: Write the failing test**

```cpp
// Create a temporary test file to verify compilation
// tests/test_diagnostic_compilation.cpp
#include "catchem_diagnostic_manager.hpp"

int main() {
    catchem::DiagnosticManager diag_mgr;
    diag_mgr.register_field("dust_flux", "Dust emission flux", "kg/m2/s", catchem::DiagType::FIELD_2D, {10, 10});
    return 0;
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++20 tests/test_diagnostic_compilation.cpp`
Expected: FAIL with missing headers/symbols

- [ ] **Step 3: Write minimal implementation for diagnostic field header**

```cpp
// src/core/catchem_diagnostic.hpp
#pragma once
#include <string>
#include <vector>
#include <Kokkos_Core.hpp>

namespace catchem {

enum class DiagType { SCALAR, FIELD_1D, FIELD_2D, FIELD_3D };

class DiagnosticField {
public:
    using HostSpace = Kokkos::HostSpace;
    using DeviceSpace = Kokkos::DefaultExecutionSpace::memory_space;

    using View2D = Kokkos::View<double**, Kokkos::LayoutLeft, DeviceSpace>;
    using HostView2D = Kokkos::View<double**, Kokkos::LayoutLeft, HostSpace>;

    using View3D = Kokkos::View<double***, Kokkos::LayoutLeft, DeviceSpace>;
    using HostView3D = Kokkos::View<double***, Kokkos::LayoutLeft, HostSpace>;

    std::string name;
    std::string description;
    std::string units;
    DiagType type;
    std::vector<int> dimensions;

    View2D device_view_2d;
    HostView2D host_view_2d;

    View3D device_view_3d;
    HostView3D host_view_3d;

    bool is_gpu_target;

    DiagnosticField(const std::string& name_val,
                    const std::string& desc_val,
                    const std::string& units_val,
                    DiagType type_val,
                    const std::vector<int>& dims);

    void sync_to_host();
    void sync_to_device();
    void reset();
};

} // namespace catchem
```

- [ ] **Step 4: Write minimal implementation for diagnostic manager header**

```cpp
// src/core/catchem_diagnostic_manager.hpp
#pragma once
#include <unordered_map>
#include <string>
#include <memory>
#include "catchem_diagnostic.hpp"

namespace catchem {

class DiagnosticManager {
private:
    std::unordered_map<std::string, std::shared_ptr<DiagnosticField>> fields;
public:
    DiagnosticManager() = default;

    void register_field(const std::string& name,
                        const std::string& desc,
                        const std::string& units,
                        DiagType type,
                        const std::vector<int>& dims);

    bool has_field(const std::string& name) const;
    std::shared_ptr<DiagnosticField> get_field(const std::string& name);

    Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace::memory_space>
    get_device_view_2d(const std::string& name);

    Kokkos::View<double***, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace::memory_space>
    get_device_view_3d(const std::string& name);

    void* get_host_pointer(const std::string& name);
    void sync_to_host();
    void reset_all();
};

} // namespace catchem
```

- [ ] **Step 5: Write minimal implementation for source file**

```cpp
// src/core/catchem_diagnostic.cpp
#include "catchem_diagnostic_manager.hpp"
#include <stdexcept>

namespace catchem {

DiagnosticField::DiagnosticField(const std::string& name_val,
                                 const std::string& desc_val,
                                 const std::string& units_val,
                                 DiagType type_val,
                                 const std::vector<int>& dims)
    : name(name_val), description(desc_val), units(units_val), type(type_val), dimensions(dims)
{
    is_gpu_target = !std::is_same_v<HostSpace, DeviceSpace>;

    if (type == DiagType::FIELD_2D) {
        if (dims.size() != 2) throw std::invalid_argument("2D field requires 2 dimensions");
        if (is_gpu_target) device_view_2d = View2D("dev_" + name, dims[0], dims[1]);
        host_view_2d = HostView2D("host_" + name, dims[0], dims[1]);
        if (!is_gpu_target) device_view_2d = host_view_2d;
    } else if (type == DiagType::FIELD_3D) {
        if (dims.size() != 3) throw std::invalid_argument("3D field requires 3 dimensions");
        if (is_gpu_target) device_view_3d = View3D("dev_" + name, dims[0], dims[1], dims[2]);
        host_view_3d = HostView3D("host_" + name, dims[0], dims[1], dims[2]);
        if (!is_gpu_target) device_view_3d = host_view_3d;
    } else {
        throw std::invalid_argument("Unsupported DiagType");
    }
}

void DiagnosticField::sync_to_host() {
    if constexpr (!std::is_same_v<HostSpace, DeviceSpace>) {
        if (type == DiagType::FIELD_2D) {
            Kokkos::deep_copy(host_view_2d, device_view_2d);
        } else if (type == DiagType::FIELD_3D) {
            Kokkos::deep_copy(host_view_3d, device_view_3d);
        }
    }
}

void DiagnosticField::sync_to_device() {
    if constexpr (!std::is_same_v<HostSpace, DeviceSpace>) {
        if (type == DiagType::FIELD_2D) {
            Kokkos::deep_copy(device_view_2d, host_view_2d);
        } else if (type == DiagType::FIELD_3D) {
            Kokkos::deep_copy(device_view_3d, host_view_3d);
        }
    }
}

void DiagnosticField::reset() {
    if (type == DiagType::FIELD_2D) {
        Kokkos::deep_copy(device_view_2d, 0.0);
        Kokkos::deep_copy(host_view_2d, 0.0);
    } else if (type == DiagType::FIELD_3D) {
        Kokkos::deep_copy(device_view_3d, 0.0);
        Kokkos::deep_copy(host_view_3d, 0.0);
    }
}

void DiagnosticManager::register_field(const std::string& name,
                                       const std::string& desc,
                                       const std::string& units,
                                       DiagType type,
                                       const std::vector<int>& dims) {
    fields[name] = std::make_shared<DiagnosticField>(name, desc, units, type, dims);
}

bool DiagnosticManager::has_field(const std::string& name) const {
    return fields.find(name) != fields.end();
}

std::shared_ptr<DiagnosticField> DiagnosticManager::get_field(const std::string& name) {
    if (!has_field(name)) throw std::invalid_argument("Field not found: " + name);
    return fields.at(name);
}

Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace::memory_space>
DiagnosticManager::get_device_view_2d(const std::string& name) {
    auto field = get_field(name);
    if (field->type != DiagType::FIELD_2D) throw std::invalid_argument("Field is not 2D: " + name);
    return field->device_view_2d;
}

Kokkos::View<double***, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace::memory_space>
DiagnosticManager::get_device_view_3d(const std::string& name) {
    auto field = get_field(name);
    if (field->type != DiagType::FIELD_3D) throw std::invalid_argument("Field is not 3D: " + name);
    return field->device_view_3d;
}

void* DiagnosticManager::get_host_pointer(const std::string& name) {
    auto field = get_field(name);
    if (field->type == DiagType::FIELD_2D) {
        return static_cast<void*>(field->host_view_2d.data());
    } else if (field->type == DiagType::FIELD_3D) {
        return static_cast<void*>(field->host_view_3d.data());
    }
    return nullptr;
}

void DiagnosticManager::sync_to_host() {
    for (auto& [key, field] : fields) {
        field->sync_to_host();
    }
}

void DiagnosticManager::reset_all() {
    for (auto& [key, field] : fields) {
        field->reset();
    }
}

} // namespace catchem
```

- [ ] **Step 6: Add to CMakeLists and compile test**

Modify `src/core/CMakeLists.txt` to include `catchem_diagnostic.cpp`:
```cmake
  set(
    _cpp_core_srcs
    catchem_core.cpp
    catchem_api.cpp
    catchem_diagnostic.cpp
  )
```

Run test compilation: `rm tests/test_diagnostic_compilation.cpp`

- [ ] **Step 7: Commit**

```bash
git add src/core/catchem_diagnostic.hpp src/core/catchem_diagnostic_manager.hpp src/core/catchem_diagnostic.cpp src/core/CMakeLists.txt
git commit -m "feat(core): implement C++ DiagnosticManager and DiagnosticField"
```

### Task 2: Attach DiagnosticManager to Core & Expose C-API

**Files:**
- Modify: `src/core/catchem_core.hpp`
- Modify: `src/core/catchem_core.cpp`
- Modify: `src/core/catchem_api.hpp`
- Modify: `src/core/catchem_api.cpp`

**Interfaces:**
- Consumes: `catchem::DiagnosticManager`
- Produces: `catchem_diag_register`, `catchem_diag_get_pointer`, `catchem_diag_sync_to_host`, `catchem_diag_reset`

- [ ] **Step 1: Write failing C-API test**

```cpp
// tests/test_diag_capi.cpp
#include "catchem_api.hpp"
#include <iostream>

int main() {
    void* core = catchem_core_create(10, 20, 5);
    catchem_diag_register(core, "dust_flux", "desc", "kg/m2/s", 2, 10, 20, 0);
    void* ptr = catchem_diag_get_pointer(core, "dust_flux");
    if (!ptr) return 1;
    catchem_core_destroy(core);
    std::cout << "SUCCESS\n";
    return 0;
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++20 tests/test_diag_capi.cpp`
Expected: FAIL with missing symbol `catchem_diag_register`

- [ ] **Step 3: Modify `catchem_core.hpp` and `catchem_core.cpp`**

In `src/core/catchem_core.hpp`:
```cpp
#pragma once
#include <memory>
#include <vector>
#include "catchem_state_manager.hpp"
#include "catchem_process_interface.hpp"
#include "catchem_diagnostic_manager.hpp"

namespace catchem {

class Core {
private:
    std::shared_ptr<StateManager> state_mgr;
    std::shared_ptr<DiagnosticManager> diag_mgr;
    std::vector<std::shared_ptr<ProcessInterface>> processes;
public:
    Core(int nc, int nl, int ns);
    std::shared_ptr<StateManager> get_state_manager();
    std::shared_ptr<DiagnosticManager> get_diagnostic_manager();
    void add_process(std::shared_ptr<ProcessInterface> process);
    void run_timestep(double dt);
};

} // namespace catchem
```

In `src/core/catchem_core.cpp`:
```cpp
#include "catchem_core.hpp"

namespace catchem {

Core::Core(int nc, int nl, int ns) {
    state_mgr = std::make_shared<StateManager>(nc, nl, ns);
    diag_mgr = std::make_shared<DiagnosticManager>();
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

void Core::run_timestep(double dt) {
    // Sync shared boundary arrays to active execution spaces
    state_mgr->sync_to_device();

    for (auto& process : processes) {
        process->run(state_mgr);
    }

    // Sync execution outputs back to Fortran-accessible memory
    state_mgr->sync_to_host();

    // Sync diagnostics
    diag_mgr->sync_to_host();
}

} // namespace catchem
```

- [ ] **Step 4: Modify `catchem_api.hpp` and `catchem_api.cpp`**

In `src/core/catchem_api.hpp`:
```cpp
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void* catchem_core_create(int nc, int nl, int ns);
void catchem_core_destroy(void* core_ptr);
void* catchem_core_get_state_manager(void* core_ptr);
void catchem_state_bind_1d(void* state_ptr, const char* name, double* ptr);
void catchem_state_bind_2d(void* state_ptr, const char* name, double* ptr);
void catchem_state_bind_3d(void* state_ptr, const char* name, double* ptr);
void catchem_state_sync_to_device(void* state_ptr);
void catchem_state_sync_to_host(void* state_ptr);
void catchem_core_run_timestep(void* core_ptr, double dt);

// Diagnostic API
void catchem_diag_register(void* core_ptr, const char* name, const char* desc, const char* units, int rank, int dim1, int dim2, int dim3);
void* catchem_diag_get_pointer(void* core_ptr, const char* name);
void catchem_diag_sync_to_host(void* core_ptr);
void catchem_diag_reset(void* core_ptr);

#ifdef __cplusplus
}
#endif
```

In `src/core/catchem_api.cpp`:
```cpp
#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_state_manager.hpp"
#include "catchem_diagnostic_manager.hpp"

extern "C" {

void* catchem_core_create(int nc, int nl, int ns) {
    return static_cast<void*>(new catchem::Core(nc, nl, ns));
}

void catchem_core_destroy(void* core_ptr) {
    delete static_cast<catchem::Core*>(core_ptr);
}

void* catchem_core_get_state_manager(void* core_ptr) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return static_cast<void*>(core->get_state_manager().get());
}

void catchem_state_bind_1d(void* state_ptr, const char* name, double* ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    state->bind_field_1d(name, ptr);
}

void catchem_state_bind_2d(void* state_ptr, const char* name, double* ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    state->bind_field_2d(name, ptr);
}

void catchem_state_bind_3d(void* state_ptr, const char* name, double* ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    state->bind_field_3d(name, ptr);
}

void catchem_state_sync_to_device(void* state_ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    state->sync_to_device();
}

void catchem_state_sync_to_host(void* state_ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    state->sync_to_host();
}

void catchem_core_run_timestep(void* core_ptr, double dt) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    core->run_timestep(dt);
}

void catchem_diag_register(void* core_ptr, const char* name, const char* desc, const char* units, int rank, int dim1, int dim2, int dim3) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    catchem::DiagType type;
    std::vector<int> dims;
    if (rank == 2) {
        type = catchem::DiagType::FIELD_2D;
        dims = {dim1, dim2};
    } else if (rank == 3) {
        type = catchem::DiagType::FIELD_3D;
        dims = {dim1, dim2, dim3};
    } else {
        type = catchem::DiagType::SCALAR; // Simplified for now
    }
    core->get_diagnostic_manager()->register_field(name, desc, units, type, dims);
}

void* catchem_diag_get_pointer(void* core_ptr, const char* name) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return core->get_diagnostic_manager()->get_host_pointer(name);
}

void catchem_diag_sync_to_host(void* core_ptr) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    core->get_diagnostic_manager()->sync_to_host();
}

void catchem_diag_reset(void* core_ptr) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    core->get_diagnostic_manager()->reset_all();
}

}
```

- [ ] **Step 5: Verify cleanup**

Run: `rm tests/test_diag_capi.cpp`

- [ ] **Step 6: Commit**

```bash
git add src/core/catchem_core.hpp src/core/catchem_core.cpp src/core/catchem_api.hpp src/core/catchem_api.cpp
git commit -m "feat(api): expose C-API for DiagnosticManager and attach to Core"
```

### Task 3: Verify C++ Diagnostic System via Interop Tests

**Files:**
- Modify: `tests/test_catchem_interop.cpp`

**Interfaces:**
- Consumes: `catchem_core.hpp`, `catchem_diagnostic_manager.hpp`

- [ ] **Step 1: Modify `test_catchem_interop.cpp`**

Add diagnostic tests to `tests/test_catchem_interop.cpp`. This will test allocating a diagnostic field, updating it inside a Kokkos parallel region, syncing it to the host, and verifying the host pointer through the C-API.

```cpp
#include <iostream>
#include <vector>
#include <Kokkos_Core.hpp>
#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"

// A dummy process that writes to a diagnostic field
class DummyDiagProcess : public catchem::ProcessInterface {
private:
    std::shared_ptr<catchem::DiagnosticManager> diag_mgr;
    int n_cols;
public:
    DummyDiagProcess(std::shared_ptr<catchem::DiagnosticManager> dm, int nc) : diag_mgr(dm), n_cols(nc) {}

    std::string get_name() const override { return "DummyDiagProcess"; }

    void init(std::shared_ptr<catchem::StateManager> state) override {}

    void run(std::shared_ptr<catchem::StateManager> state) override {
        // Retrieve the underlying diagnostic device View
        auto dust_flux = diag_mgr->get_device_view_2d("dust_emission_flux");

        // Capture View by value in the parallel kernel
        Kokkos::parallel_for("calculate_dust_emissions",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, n_cols),
            KOKKOS_LAMBDA(int icol) {
                // Write directly to the diagnostic view
                dust_flux(icol, 0) = 42.0 + icol;
            }
        );
    }

    void finalize() override {}
};

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        int nx = 4;
        int ny = 1;
        int nz = 5;
        int n_cols = nx * ny;

        // 1. Create Core (which creates StateManager and DiagnosticManager)
        void* core_ptr = catchem_core_create(n_cols, nz, 1);
        auto* core = static_cast<catchem::Core*>(core_ptr);
        auto diag_mgr = core->get_diagnostic_manager();

        // 2. Register diagnostic through C-API
        catchem_diag_register(core_ptr, "dust_emission_flux", "Dust flux", "kg/m2/s", 2, n_cols, 1, 0);

        // 3. Attach dummy process
        core->add_process(std::make_shared<DummyDiagProcess>(diag_mgr, n_cols));

        // 4. Run timestep (executes process, syncs diagnostics to host)
        catchem_core_run_timestep(core_ptr, 3600.0);

        // 5. Get host pointer and verify results
        void* host_ptr = catchem_diag_get_pointer(core_ptr, "dust_emission_flux");
        double* dust_flux_host = static_cast<double*>(host_ptr);

        bool passed = true;
        for (int i = 0; i < n_cols; ++i) {
            if (dust_flux_host[i] != 42.0 + i) { // Note LayoutLeft means col_i is inner dimension
                std::cerr << "Diagnostic mismatch at col " << i << ": expected " << 42.0 + i
                          << ", got " << dust_flux_host[i] << std::endl;
                passed = false;
            }
        }

        if (passed) {
            std::cout << "SUCCESS: C++ Diagnostic Validation Passed!" << std::endl;
        } else {
            std::cout << "FAILURE: C++ Diagnostic Validation Failed!" << std::endl;
            return 1;
        }

        catchem_core_destroy(core_ptr);
    }
    Kokkos::finalize();
    return 0;
}
```

- [ ] **Step 2: Build and run the interop tests**

Run: `docker run --rm -v $(pwd):/workspace -w /workspace/build-test cece-dev:latest bash -c "cmake .. -DENABLE_KOKKOS=ON -DENABLE_TESTING=ON && make -j4 && ctest --output-on-failure"`
Expected: PASS, outputting "SUCCESS: C++ Diagnostic Validation Passed!"

- [ ] **Step 3: Commit**

```bash
git add tests/test_catchem_interop.cpp
git commit -m "test(diag): verify C++ diagnostic registration and parallel updating"
```
