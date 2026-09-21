# Technical Design Specification: Modernizing CATChem Diagnostic Infrastructure to C++/Kokkos

## 1. Context & Motivation

In the legacy CATChem framework, diagnostic variables (such as emission fluxes, deposition velocities, and physical intermediate calculations) are managed via `DiagnosticInterface_Mod.F90` and `DiagnosticManager_Mod.F90` in Fortran. This system relies on dynamic allocations of nested union-like types (`DiagnosticDataType`) and requires manual, serial population during the physical processes' run loop.

As CATChem's physical kernels are ported to high-performance C++20/Kokkos (Phase 2), writing to diagnostics within parallel GPU regions (such as `Kokkos::parallel_for`) under the legacy system would introduce massive performance bottlenecks. Specifically, thread-safe updates would either require locking/atomic operations, or force costly host-device-host synchronizations.

By modernizing `DiagnosticInterface` and `DiagnosticManager` into a unified, performance-portable C++ library, we achieve:
*   **Direct GPU-Safe Updates:** C++ parallel kernels can capture Kokkos-backed diagnostic views by value and write directly to them with zero-copy execution.
*   **Zero-Overhead on CPU Spaces:** On CPU execution targets, C++ diagnostics leverage Fortran-compatible column-major layouts (`Kokkos::LayoutLeft`), aligning memory directly with host NetCDF outputters.
*   **Decoupled Metadata & Storage:** Diagnostic variables are managed dynamically via a string-keyed map, avoiding macro include code generation entirely.

---

## 2. Core Architecture

The C++ diagnostic infrastructure resides in `src/core/catchem_diagnostic.hpp` and `src/core/catchem_diagnostic.cpp`, replacing the legacy Fortran files entirely in the core computation block, while preserving a thin Fortran interop shell for standard NetCDF file IO.

```
                +-------------------------------------------------+
                |             catchem::DiagnosticField            |
                | - name, description, units, output_frequency    |
                | - Host/Device Kokkos Views (LayoutLeft)         |
                +-------------------------------------------------+
                                         |
                                         v
                +-------------------------------------------------+
                |            catchem::DiagnosticManager           |
                | - Registry map: std::string -> DiagnosticField  |
                | - Methods: register_field(), get_device_view()  |
                | - sync_to_host() boundaries                      |
                +-------------------------------------------------+
                                         |
                                         v
                +-------------------------------------------------+
                |                  Extern "C" API                 |
                | - catchem_diag_get_pointer_2d()                 |
                | - catchem_diag_get_pointer_3d()                 |
                +-------------------------------------------------+
                                         |
                                         v (c_f_pointer association)
                +-------------------------------------------------+
                |           Thin Preserved Fortran Shell          |
                | - Drives NetCDF outputters using C++ pointers   |
                +-------------------------------------------------+
```

---

## 3. Detailed Class & Struct Interfaces

### 3.1 `catchem::DiagType` Enum
Defines the dimensions of diagnostic variables:
```cpp
namespace catchem {
enum class DiagType {
    SCALAR,
    FIELD_1D,
    FIELD_2D,
    FIELD_3D
};
} // namespace catchem
```

### 3.2 `catchem::DiagnosticField`
Manages the metadata, host/device View structures, and execution-space coherence checks for an individual diagnostic parameter.
*   Uses `Kokkos::LayoutLeft` (Column-Major) layout mapping.
*   Maintains separate device views and host mirror views when targeting discrete accelerators (GPUs), syncing them only at thread-safe boundaries.

```cpp
// src/core/catchem_diagnostic.hpp
#pragma once
#include <string>
#include <vector>
#include <memory>
#include <Kokkos_Core.hpp>

namespace catchem {

class DiagnosticField {
public:
    using HostSpace = Kokkos::HostSpace;
    using DeviceSpace = Kokkos::DefaultExecutionSpace::memory_space;

    // Kokkos Views matching Column-Major layouts (2D and 3D)
    using View2D = Kokkos::View<double**, Kokkos::LayoutLeft, DeviceSpace>;
    using HostView2D = Kokkos::View<double**, Kokkos::LayoutLeft, HostSpace>;

    using View3D = Kokkos::View<double***, Kokkos::LayoutLeft, DeviceSpace>;
    using HostView3D = Kokkos::View<double***, Kokkos::LayoutLeft, HostSpace>;

    std::string name;
    std::string description;
    std::string units;
    DiagType type;
    std::vector<int> dimensions;

    // Allocated View instances
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

### 3.3 `catchem::DiagnosticManager`
Acts as the global, dynamic registry of all registered simulation diagnostics.
*   Integrates with `catchem::Core` lifecycle.
*   Provides dynamic, type-safe getters for processes to capture underlying Views upfront before launching parallel loops.

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

    // Upfront View Getters for capture in Kokkos parallel kernels
    Kokkos::View<double**, Kokkos::LayoutLeft, DefaultExecutionSpace::memory_space>
    get_device_view_2d(const std::string& name);

    Kokkos::View<double***, Kokkos::LayoutLeft, DefaultExecutionSpace::memory_space>
    get_device_view_3d(const std::string& name);

    void* get_host_pointer(const std::string& name);
    void sync_to_host();
    void reset_all();
};

} // namespace catchem
```

---

## 4. Execution Flow & Kernel Captures

### 4.1 Upfront Capture in C++ Parallel Regions (Option 1)
When a concrete C++ process executes, it extracts the target `Kokkos::View` from the manager and captures it by value in its GPU parallel functor.

```cpp
void DustProcess::run(std::shared_ptr<StateManager> state, std::shared_ptr<DiagnosticManager> diag_mgr) {
    int n_cols = state->n_cols;

    // 1. Retrieve the underlying diagnostic device View
    auto dust_flux = diag_mgr->get_device_view_2d("dust_emission_flux");

    // 2. Capture View by value in the parallel kernel
    Kokkos::parallel_for("calculate_dust_emissions",
        Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, n_cols),
        KOKKOS_LAMBDA(int icol) {
            double emission_rate = run_dust_scheme_on_column(icol);

            // 3. Write directly to the diagnostic view with zero-copy execution on GPU space
            dust_flux(icol, 0) = emission_rate;
        }
    );
}
```

### 4.2 Timestep Sync and File Extraction
At the end of a timestep:
1.  **`diag_mgr->sync_to_host()`** is triggered.
    *   If running on CPUs (where `DeviceSpace == HostSpace`), this is compiled away as a **no-op** with zero computational overhead.
    *   If running on GPU accelerators, `Kokkos::deep_copy` is invoked once per diagnostic field to pull device memory mirrors back into unified host memory.
2.  The thin outer Fortran model driver accesses the dynamic bind pointer via the extern "C" API and maps it into its Fortran pointers using standard `iso_c_binding` mechanics.

---

## 5. Extern "C" Binding Boundary

To keep the existing NetCDF output infrastructure completely intact on the host Fortran side, we define thin, extern "C" binders returning raw double pointers from the diagnostic fields' unified host views:

```cpp
// src/core/catchem_api.hpp
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void catchem_diag_register(void* core_ptr, const char* name, const char* desc, const char* units, int rank, int dim1, int dim2, int dim3);
void* catchem_diag_get_pointer(void* core_ptr, const char* name);
void catchem_diag_sync_to_host(void* core_ptr);
void catchem_diag_reset(void* core_ptr);

#ifdef __cplusplus
}
#endif
```

```cpp
// src/core/catchem_api.cpp
#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"

extern "C" {

void catchem_diag_register(void* core_ptr, const char* name, const char* desc, const char* units, int rank, int dim1, int dim2, int dim3) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    catchem::DiagType type;
    std::vector<int> dims;
    if (rank == 2) {
        type = catchem::DiagType::FIELD_2D;
        dims = {dim1, dim2};
    } else {
        type = catchem::DiagType::FIELD_3D;
        dims = {dim1, dim2, dim3};
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

---

## 6. Porting and Integration Roadmap

We integrate this diagnostic modernization plan smoothly into the Phase 2 roadmap:

1.  **Task 1:** Create `catchem_diagnostic.hpp` and `catchem_diagnostic_manager.hpp` with full Kokkos View allocations and host mirror synchronizers.
2.  **Task 2:** Attach `catchem::DiagnosticManager` as a shared member owned by `catchem::Core`.
3.  **Task 3:** Expose the dynamic C-API diagnostic getters (`catchem_diag_register`, `catchem_diag_get_pointer`).
4.  **Task 4:** Refactor `tests/test_catchem_interop.cpp` to verify that capturing diagnostic Views inside a GPU parallel region updates host buffers sequentially and perfectly coordinates with Fortran `c_f_pointer` wrappers.
