# Technical Design: Modernizing CATChem Core with C++ & Kokkos

**Date:** July 8, 2026  
**Status:** Approved  
**Topic:** Transitioning CATChem Core to a C++ modern standard using Kokkos::View and the Kokkos backport of mdspan (`std::experimental::mdspan`) targeting C++20 for high-performance physics/chemistry interop in CCPP and host frameworks.

---

## 1. Executive Summary & Objectives

This design details the replacement of the Fortran operational core (`src/core/**`) with a native, GPU-enabled C++ implementation powered by **Kokkos**.

### Primary Goals:
* **GPU & CPU Interoperability:** Provide compile-time toggles and runtime optimization for CPU (zero-copy) and GPU execution.
* **Modern C++ Semantics:** Transition from Fortran pointer-virtualization to modern multidimensional indexing using the Kokkos backport of `mdspan` (`std::experimental::mdspan`) targeting C++20.
* **CCPP Compatibility:** Retain full compatibility with the CCPP framework and host models (e.g., UFS, FV3) without altering external `.meta` files or CCPP argument tables.
* **Code Consolidation:** Deprecate complex custom build-time Python code generation, virtual columns, and extensive unpacking utility routines.

---

## 2. Memory Management Layer (`InteropField`)

To handle the Fortran-C++ boundary across varying HPC compute node architectures, a template-based `InteropField` class manages both zero-copy CPU views and GPU-mirrored memory.

### Concept:
* **Host Space:** Wraps Fortran array memory pointers without copy using unmanaged `Kokkos::View` with `Kokkos::LayoutLeft` (column-major to preserve Fortran layout).
* **Device Space:** Allocates native managed memory on the active device.
* **Synchronizations:** Synclines (`sync_to_device()`, `sync_to_host()`) resolve to compile-time no-ops when executing on CPU-only spaces.

```cpp
// src/core/catchem_interop_field.hpp
#pragma once
#include <Kokkos_Core.hpp>
#include <vector>
#include <memory>

namespace catchem {

template <typename DataType, int Rank>
class InteropField {
public:
    using HostSpace = Kokkos::HostSpace;
    using DeviceSpace = Kokkos::DefaultExecutionSpace::memory_space;

    // Helper to resolve unmanaged View type based on Rank
    template <typename T, int R, typename Space, bool Unmanaged>
    struct ViewType;

    // Specializations for Rank 1, 2, 3
    template <typename T, typename Space, bool Unmanaged>
    struct ViewType<T, 1, Space, Unmanaged> {
        using type = typename std::conditional_t<Unmanaged,
            Kokkos::View<T*, Kokkos::LayoutLeft, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
            Kokkos::View<T*, Kokkos::LayoutLeft, Space>>;
    };

    template <typename T, typename Space, bool Unmanaged>
    struct ViewType<T, 2, Space, Unmanaged> {
        using type = typename std::conditional_t<Unmanaged,
            Kokkos::View<T**, Kokkos::LayoutLeft, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
            Kokkos::View<T**, Kokkos::LayoutLeft, Space>>;
    };

    template <typename T, typename Space, bool Unmanaged>
    struct ViewType<T, 3, Space, Unmanaged> {
        using type = std::conditional_t<Unmanaged,
            Kokkos::View<T***, Kokkos::LayoutLeft, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
            Kokkos::View<T***, Kokkos::LayoutLeft, Space>>;
    };

    using HostViewType = typename ViewType<DataType, Rank, HostSpace, true>::type;
    using DeviceViewType = typename ViewType<DataType, Rank, DeviceSpace, false>::type;

    HostViewType host_view;
    DeviceViewType device_view;
    bool is_gpu_target;

    InteropField(DataType* ptr, const std::vector<int>& dims) {
        is_gpu_target = !std::is_same_v<HostSpace, DeviceSpace>;

        if constexpr (Rank == 1) {
            host_view = HostViewType(ptr, dims[0]);
            if (is_gpu_target) device_view = DeviceViewType("dev_field_1d", dims[0]);
        } else if constexpr (Rank == 2) {
            host_view = HostViewType(ptr, dims[0], dims[1]);
            if (is_gpu_target) device_view = DeviceViewType("dev_field_2d", dims[0], dims[1]);
        } else if constexpr (Rank == 3) {
            host_view = HostViewType(ptr, dims[0], dims[1], dims[2]);
            if (is_gpu_target) device_view = DeviceViewType("dev_field_3d", dims[0], dims[1], dims[2]);
        }
    }

    void sync_to_device() {
        if (is_gpu_target) {
            Kokkos::deep_copy(device_view, host_view);
        }
    }

    void sync_to_host() {
        if (is_gpu_target) {
            Kokkos::deep_copy(host_view, device_view);
        }
    }

    // Accesses active execution view
    auto view() const {
        if constexpr (!std::is_same_v<HostSpace, DeviceSpace>) {
            return device_view;
        } else {
            return host_view;
        }
    }
};

} // namespace catchem
```

---

## 3. Core Architecture & Interoperability Layer

The C++ core mimics the organizational boundary of CATChem, exposing standard C-linkage entrypoints (`extern "C"`) to the Fortran wrappers.

```
       +------------------------------------+
       |          CCPP Host Model           |
       +-----------------+------------------+
                         |
                         | (Fortran Arrays)
                         v
       +-----------------+------------------+
       |   ccpp_catchem_interface.F90       |  <- Exposes CCPP schemes
       +-----------------+------------------+
                         |
                         | (ISO_C_BINDING raw pointers)
                         v
       +-----------------+------------------+
       |      CATChem_API C-Bindings        |  <- Thin Fortran/C boundary
       +-----------------+------------------+
                         |
                         | (extern "C" C-API)
                         v
       +-----------------+------------------+
       |          C++ Core Engine           |
       |  - catchem::Core                   |
       |  - catchem::StateManager           |
       |  - Maps pointers to InteropFields  |
       +-----------------+------------------+
                         |
                         | (parallel loops)
                         v
       +-----------------+------------------+
       |       Kokkos Kernels (CPU/GPU)     |
       +------------------------------------+
```

### C++ State Manager
Manages the string-to-field map of active met and chemical state arrays:

```cpp
// src/core/catchem_state_manager.hpp
#pragma once
#include <unordered_map>
#include <string>
#include <memory>
#include "catchem_interop_field.hpp"

namespace catchem {

class StateManager {
public:
    int n_cols, n_levels, n_species;

    std::unordered_map<std::string, std::shared_ptr<InteropField<double, 1>>> fields_1d;
    std::unordered_map<std::string, std::shared_ptr<InteropField<double, 2>>> fields_2d;
    std::unordered_map<std::string, std::shared_ptr<InteropField<double, 3>>> fields_3d;

    StateManager(int nc, int nl, int ns) : n_cols(nc), n_levels(nl), n_species(ns) {}

    void bind_field_2d(const std::string& name, double* ptr) {
        fields_2d[name] = std::make_shared<InteropField<double, 2>>(ptr, std::vector<int>{n_cols, n_levels});
    }

    void bind_field_3d(const std::string& name, double* ptr) {
        fields_3d[name] = std::make_shared<InteropField<double, 3>>(ptr, std::vector<int>{n_cols, n_levels, n_species});
    }

    void sync_to_device() {
        for (auto& [k, v] : fields_2d) v->sync_to_device();
        for (auto& [k, v] : fields_3d) v->sync_to_device();
    }

    void sync_to_host() {
        for (auto& [k, v] : fields_2d) v->sync_to_host();
        for (auto& [k, v] : fields_3d) v->sync_to_host();
    }
};

} // namespace catchem
```

---

## 4. Porting/Replacement Plan for `src/core/**`

To achieve complete architectural coverage, this section maps every single one of the 24 files currently in the `src/core/` directory to its destiny in the modernized C++20 / Kokkos code structure:

### Replaced (Rewritten in C++):
* **`CATChemCore_Mod.F90`** $\rightarrow$ `catchem::Core` (Central entry point; orchestrates timestepping, physical scheme dispatch, and diagnostic collection).
* **`StateManager_Mod.F90`** $\rightarrow$ `catchem::StateManager` (Dynamic key-value registering mapping field names to raw multi-dimensional `InteropField` objects).
* **`chemstate_mod.F90`** $\rightarrow$ Consolidated into C++ `catchem::StateManager` and associated `catchem::Species` structures for chemistry concentrations.
* **`metstate_mod.F90`** $\rightarrow$ Dynamic key-value registration maps in C++ `StateManager` holding individual meteorological fields dynamically mapped to string names, removing macro generation entirely.
* **`TimeState_Mod.F90`** $\rightarrow$ `catchem_time.hpp` (defines `catchem::TimeState` holding calendar data and providing device-compilable astronomical solar zenith angle functions via `KOKKOS_INLINE_FUNCTION`).
* **`species_mod.F90`** $\rightarrow$ `catchem_species.hpp` (defines `catchem::Species` metadata structs for molecular weights, physical radii, and classification flags).
* **`GridGeometry_Mod.F90` & `GridManager_Mod.F90`** $\rightarrow$ `catchem_grid.hpp` (defines `catchem::Grid` holding spatial extents, spacing, and computing local cell areas/volumes; utilizes Kokkos subviews for zero-copy vertical profiling).
* **`constants.F90`** $\rightarrow$ `catchem_constants.hpp` (defines compile-time standard values for standard atmospheric pressure, Avogadro's number, thermodynamic heats, and universal physical constants).
* **`error_mod.F90`** $\rightarrow$ `catchem_error.hpp` (defines modern exception handling, thread-safe warning loggers, and context tracing).
* **`UnitConversion_Mod.F90`** $\rightarrow$ `catchem_unit_conversion.hpp` (defines header-only, device-compilable inline routines for converting concentration units [ppbv, ug/m3, molec/cm3] and deposition/emission fluxes).
* **`utilities_mod.F90`** $\rightarrow$ `catchem_utilities.hpp` (defines standard host-side mathematical helpers, safe divisions, and geopotential calculations).
* **`met_utilities_mod.F90`** $\rightarrow$ `catchem_met_utilities.hpp` (defines device-compilable inline thermodynamics including potential/virtual temperature, Monin-Obukhov scale heights, and Cunningham slippage correction factors).
* **`ProcessInterface_Mod.F90`** $\rightarrow$ Abstract virtual base class `catchem::ProcessInterface` defining the execution lifecycle interfaces for physical schemes.
* **`ProcessManager_Mod.F90`** $\rightarrow$ Ported to C++ timeline dispatchers orchestrating sequential or concurrent (Kokkos timeline parallelized) execution steps.
* **`ProcessRegistry_Mod.F90`** $\rightarrow$ Ported to C++ process class registries.
* **`ProcessFactory_Mod.F90`** $\rightarrow$ Reimplemented using C++ creator lambdas and registry factories.
* **`ExtEmisData_Mod.F90`** $\rightarrow$ Ported to C++ emissions handling to cleanly support device-level external emission inputs.
* **`EmissionConfigValidator_Mod.F90`** $\rightarrow$ Reimplemented as static C++ verification routines during emissions registry setup.

### Dropped Entirely (Obsoleted):
* **`VirtualColumn_Mod.F90`** (Obsoleted by zero-copy column indexing via `Kokkos::subview` in C++ / Kokkos kernels).
* **`ChemSpeciesUtils_Mod.F90`** (Obsoleted by modern C++ STL lists, mapping tables, and dynamic lookups via string hash maps).
* **`generate_metstate_macros.py`** & all auto-generated include targets (`virtualmet_type.inc`, `virtualmet_populate.inc`, `virtualmet_cleanup.inc`, `metstate_*.inc`) (Obsoleted; metadata and accessor functions are resolved dynamically via string maps and InteropField wrappers in C++).

### Adapted (Thin Fortran wrappers preserved for interop):
* **`ConfigManager_Mod.F90`** (Retained in Fortran to utilize YAML parsing library, extracting basic dimensions and physical flags, then passing scalar parameters to the C++ core constructor).
* **`DiagnosticManager_Mod.F90` & `DiagnosticInterface_Mod.F90`** (Retained in Fortran to expose diagnosed variables to host-models' NetCDF outputters, syncing modified diagnostic memory buffers back from C++ to Fortran at the end of each timestep).

---

## 5. Phased Incremental Migration Plan

To maintain standard operational readiness and continuous verification during this massive update, we utilize a multi-phase **Strangler Fig Pattern**:

```
+--------------------------------------------------------------+
| Phase 1: Dual Core / Boundary Setup                          |
| - Implement C++ StateManager, InteropField and C-API bindings|
| - Retain legacy Fortran Core                                 |
+-------------------------------------------------------------+
                               |
                               v
+-------------------------------------------------------------+
| Phase 2: Process-by-Process Migration                        |
| - Port Settling & DryDep kernels to Kokkos C++              |
| - Toggle Fortran wrapper to invoke new C-API dispatch       |
| - Automatically sync memory in InteropField                 |
+-------------------------------------------------------------+
                               |
                               v
+-------------------------------------------------------------+
| Phase 3: Final Consolidation                                 |
| - Delete legacy Fortran core files and drop VirtualColumn    |
| - Leave only thin Fortran CCPP/API layer                    |
+-------------------------------------------------------------+
```

### Verification & Testing Strategy:
1. **Bitwise Validation:** Retain the original Fortran processes in a legacy namespace/directory during development. Run differential regression runs to verify C++ and Fortran calculations match down to numerical precision tolerances.
2. **Step-by-Step Porting:** Work on small modules (e.g. *Settling*) first, ensuring C++ executions compile and run on both CPU and GPU before proceeding to more complex chemistry modules.
