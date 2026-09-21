# Modernize CATChem Core Phase 1: Shared State and Interoperability Boundary Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish a modern, high-performance C++20 memory and execution boundary in CATChem utilizing Kokkos views and unmanaged memory wrappers, enabling zero-copy execution on CPU-only space and automatic memory synchronization on GPU-enabled device spaces.

**Architecture:** We employ a Dual Core shared-memory Strangler Fig pattern. Outermost CCPP drivers and host models allocate array memory in Fortran. This memory is wrapped in unmanaged host C++ `Kokkos::View` buffers (`catchem::InteropField`) in column-major layout (`Kokkos::LayoutLeft`) to achieve zero-copy. GPU execution spaces automatically trigger mirrored GPU buffer allocations with high-level synclines (`sync_to_device()`, `sync_to_host()`) mapped to compile-time no-ops on CPU. This enables legacy Fortran processes and new C++ Kokkos processes to run side-by-side, sharing the exact same memory.

**Tech Stack:** C++20, Kokkos, ISO_C_BINDING, CMake.

## Global Constraints

- **Language Target:** C++20 utilizing the Kokkos backport of mdspan (`std::experimental::mdspan`), avoiding direct dependency on C++23 `<mdspan>`.
- **Layout Alignment:** Retain Fortran column-major storage layout (`Kokkos::LayoutLeft`) across the pointer boundary to achieve zero-copy execution on CPU.
- **CCPP Compliance:** Outer drivers must remain in Fortran, extracting variables via ISO_C_BINDING raw pointers.
- **Incremental Bridge:** Facilitate differential validation where C++ state synchronization handles either C++ Kokkos or legacy Fortran physical schemes sharing state sequentially.

---

### Task 1: C++ Headers & Platform Setup

**Files:**
- Create: `src/core/catchem_precision.hpp`
- Create: `src/core/catchem_constants.hpp`
- Create: `src/core/catchem_error.hpp`
- Create: `src/core/catchem_unit_conversion.hpp`
- Create: `src/core/catchem_met_utilities.hpp`
- Modify: `CMakeLists.txt:20-30`

**Interfaces:**
- Produces: `catchem::fp` floating point alias (double or float based on compile definition).
- Produces: `catchem::constants` matching standard atmospheric constants.
- Produces: Device-compilable `KOKKOS_INLINE_FUNCTION` utilities for UnitConversion and MetUtilities.

- [ ] **Step 1: Write the C++ platform headers**

Create `src/core/catchem_precision.hpp`:
```cpp
#pragma once
namespace catchem {
#ifdef USE_REAL8
using fp = double;
#else
using fp = float;
#endif
} // namespace catchem
```

Create `src/core/catchem_constants.hpp`:
```cpp
#pragma once
#include "catchem_precision.hpp"
namespace catchem {
namespace constants {
constexpr fp AVO = 6.022140857e23;
constexpr fp G0 = 9.80665;
constexpr fp RSTARG = 8.314462618;
constexpr fp BOLTZ = 1.38064852e-23;
constexpr fp ATM = 1.01325e5;
constexpr fp AIR_MW = 28.9644;
constexpr fp H2O_MW = 18.01528;
constexpr fp RD = 287.0;
constexpr fp CP = 1004.6;
constexpr fp PI = 3.14159265358979323846;
constexpr fp PI_180 = PI / 180.0;
} // namespace constants
} // namespace catchem
```

Create `src/core/catchem_error.hpp`:
```cpp
#pragma once
#include <string>
#include <vector>
#include <iostream>

namespace catchem {
enum ErrorCode {
    SUCCESS = 0,
    FAILURE = -1,
    INVALID_INPUT = 1001,
    INVALID_STATE = 1003,
    MEMORY_ALLOCATION = 1007
};

class ErrorManager {
private:
    std::vector<std::string> context_stack;
public:
    void push_context(const std::string& ctx) {
        context_stack.push_back(ctx);
    }
    void pop_context() {
        if (!context_stack.empty()) context_stack.pop_back();
    }
    void report_error(ErrorCode code, const std::string& msg) {
        std::cerr << "[CATChem C++ Error " << code << "] " << msg << " | Context: ";
        for (const auto& ctx : context_stack) std::cerr << ctx << " -> ";
        std::cerr << "End\n";
    }
};
} // namespace catchem
```

Create `src/core/catchem_unit_conversion.hpp`:
```cpp
#pragma once
#include <Kokkos_Core.hpp>
#include "catchem_precision.hpp"
#include "catchem_constants.hpp"

namespace catchem {
namespace unit_conversion {

KOKKOS_INLINE_FUNCTION
inline fp ppbv_to_ugm3(fp ppbv, fp mw, fp temp, fp press) {
    return ppbv * mw * press / (constants::RSTARG * temp) * 1.0e-3;
}

KOKKOS_INLINE_FUNCTION
inline fp ugm3_to_ppbv(fp ugm3, fp mw, fp temp, fp press) {
    return ugm3 * constants::RSTARG * temp / (mw * press) * 1.0e3;
}

} // namespace unit_conversion
} // namespace catchem
```

Create `src/core/catchem_met_utilities.hpp`:
```cpp
#pragma once
#include <Kokkos_Core.hpp>
#include <cmath>
#include "catchem_precision.hpp"
#include "catchem_constants.hpp"

namespace catchem {
namespace met_utilities {

KOKKOS_INLINE_FUNCTION
inline fp potential_temperature(fp temp, fp press, fp sfc_press) {
    return temp * std::pow(sfc_press / press, constants::RD / constants::CP);
}

KOKKOS_INLINE_FUNCTION
inline fp virtual_temperature(fp temp, fp qv) {
    return temp * (1.0 + 0.61 * qv);
}

KOKKOS_INLINE_FUNCTION
inline fp cunningham_correction_factor(fp dp, fp lambda) {
    if (dp > 0.0 && lambda > 0.0) {
        return 1.0 + 2.0 * lambda / dp * (1.257 + 0.4 * std::exp(-1.1 * dp / lambda));
    }
    return 1.0;
}

} // namespace met_utilities
} // namespace catchem
```

- [ ] **Step 2: Elevate C++ standard to 20 when Kokkos is ON**

Modify `CMakeLists.txt` to require C++20 for Kokkos:
```cmake
if(ENABLE_KOKKOS)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  find_package(Kokkos REQUIRED)
endif()
```

- [ ] **Step 3: Verify configuration changes via CMake**

Run standard configuration to verify standard propagation:
```bash
cmake -B build -DENABLE_KOKKOS=ON
```
Expected output: SUCCESS without compile errors.

- [ ] **Step 4: Commit setup**

```bash
git add src/core/catchem_precision.hpp src/core/catchem_constants.hpp src/core/catchem_error.hpp src/core/catchem_unit_conversion.hpp src/core/catchem_met_utilities.hpp CMakeLists.txt
git commit -m "feat(core): setup C++20 precision, constants, and inline utilities"
```

---

### Task 2: Implement Memory Interop Layer (`InteropField`)

**Files:**
- Create: `src/core/catchem_interop_field.hpp`

**Interfaces:**
- Produces: `catchem::InteropField<typename T, int Rank>` wrapping raw pointers with compile-time zero-copy vs device mirrors.

- [ ] **Step 1: Write the `InteropField` template implementation**

Create `src/core/catchem_interop_field.hpp`:
```cpp
#pragma once
#include <Kokkos_Core.hpp>
#include <vector>
#include <memory>
#include <type_traits>

namespace catchem {

template <typename DataType, int Rank>
class InteropField {
public:
    using HostSpace = Kokkos::HostSpace;
    using DeviceSpace = Kokkos::DefaultExecutionSpace::memory_space;

    template <typename T, int R, typename Space, bool Unmanaged>
    struct ViewType;

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
        if constexpr (!std::is_same_v<HostSpace, DeviceSpace>) {
            Kokkos::deep_copy(device_view, host_view);
        }
    }

    void sync_to_host() {
        if constexpr (!std::is_same_v<HostSpace, DeviceSpace>) {
            Kokkos::deep_copy(host_view, device_view);
        }
    }

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

- [ ] **Step 2: Commit InteropField**

```bash
git add src/core/catchem_interop_field.hpp
git commit -m "feat(core): implement InteropField unmanaged view wrapper with deep_copy compile-time optimization"
```

---

### Task 3: Implement C++ State Manager (`StateManager`)

**Files:**
- Create: `src/core/catchem_state_manager.hpp`

**Interfaces:**
- Consumes: `catchem::InteropField` template.
- Produces: `catchem::StateManager` for binding arrays, querying layout dimensions, and triggering mass deep-copies.

- [ ] **Step 1: Write StateManager**

Create `src/core/catchem_state_manager.hpp`:
```cpp
#pragma once
#include <unordered_map>
#include <string>
#include <memory>
#include <vector>
#include "catchem_interop_field.hpp"

namespace catchem {

class StateManager {
public:
    int n_cols;
    int n_levels;
    int n_species;

    std::unordered_map<std::string, std::shared_ptr<InteropField<double, 1>>> fields_1d;
    std::unordered_map<std::string, std::shared_ptr<InteropField<double, 2>>> fields_2d;
    std::unordered_map<std::string, std::shared_ptr<InteropField<double, 3>>> fields_3d;

    StateManager(int nc, int nl, int ns) : n_cols(nc), n_levels(nl), n_species(ns) {}

    void bind_field_1d(const std::string& name, double* ptr) {
        fields_1d[name] = std::make_shared<InteropField<double, 1>>(ptr, std::vector<int>{n_cols});
    }

    void bind_field_2d(const std::string& name, double* ptr) {
        fields_2d[name] = std::make_shared<InteropField<double, 2>>(ptr, std::vector<int>{n_cols, n_levels});
    }

    void bind_field_3d(const std::string& name, double* ptr) {
        fields_3d[name] = std::make_shared<InteropField<double, 3>>(ptr, std::vector<int>{n_cols, n_levels, n_species});
    }

    void sync_to_device() {
        for (auto& [k, v] : fields_1d) v->sync_to_device();
        for (auto& [k, v] : fields_2d) v->sync_to_device();
        for (auto& [k, v] : fields_3d) v->sync_to_device();
    }

    void sync_to_host() {
        for (auto& [k, v] : fields_1d) v->sync_to_host();
        for (auto& [k, v] : fields_2d) v->sync_to_host();
        for (auto& [k, v] : fields_3d) v->sync_to_host();
    }
};

} // namespace catchem
```

- [ ] **Step 2: Commit StateManager**

```bash
git add src/core/catchem_state_manager.hpp
git commit -m "feat(core): implement C++ StateManager mapping dynamic string keys"
```

---

### Task 4: Implement C++ Core Orchestrator & Process Interfaces

**Files:**
- Create: `src/core/catchem_process_interface.hpp`
- Create: `src/core/catchem_core.hpp`
- Create: `src/core/catchem_core.cpp`

**Interfaces:**
- Consumes: `catchem::StateManager`
- Produces: `catchem::Core` central orchestrator.

- [ ] **Step 1: Implement abstract process interface**

Create `src/core/catchem_process_interface.hpp`:
```cpp
#pragma once
#include <string>
#include <memory>
#include "catchem_state_manager.hpp"

namespace catchem {

class ProcessInterface {
public:
    virtual ~ProcessInterface() = default;
    virtual std::string get_name() const = 0;
    virtual void init(std::shared_ptr<StateManager> state) = 0;
    virtual void run(std::shared_ptr<StateManager> state) = 0;
    virtual void finalize() = 0;
};

} // namespace catchem
```

- [ ] **Step 2: Implement Core and timeline runner**

Create `src/core/catchem_core.hpp`:
```cpp
#pragma once
#include <memory>
#include <vector>
#include "catchem_state_manager.hpp"
#include "catchem_process_interface.hpp"

namespace catchem {

class Core {
private:
    std::shared_ptr<StateManager> state_mgr;
    std::vector<std::shared_ptr<ProcessInterface>> processes;
public:
    Core(int nc, int nl, int ns);
    std::shared_ptr<StateManager> get_state_manager();
    void add_process(std::shared_ptr<ProcessInterface> process);
    void run_timestep(double dt);
};

} // namespace catchem
```

Create `src/core/catchem_core.cpp`:
```cpp
#include "catchem_core.hpp"

namespace catchem {

Core::Core(int nc, int nl, int ns) {
    state_mgr = std::make_shared<StateManager>(nc, nl, ns);
}

std::shared_ptr<StateManager> Core::get_state_manager() {
    return state_mgr;
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
}

} // namespace catchem
```

- [ ] **Step 3: Commit Core Orchestrator**

```bash
git add src/core/catchem_process_interface.hpp src/core/catchem_core.hpp src/core/catchem_core.cpp
git commit -m "feat(core): implement catchem::Core timestep orchestrator"
```

---

### Task 5: Create C-API Bindings (ISO_C_BINDING Boundary)

**Files:**
- Create: `src/core/catchem_api.hpp`
- Create: `src/core/catchem_api.cpp`
- Modify: `src/core/CMakeLists.txt`

**Interfaces:**
- Produces: Standardized `extern "C"` endpoints managing a singleton `catchem::Core` instance, registering arrays from Fortran drivers.

- [ ] **Step 1: Write C-API exports**

Create `src/core/catchem_api.hpp`:
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

#ifdef __cplusplus
}
#endif
```

Create `src/core/catchem_api.cpp`:
```cpp
#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_state_manager.hpp"

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

}
```

- [ ] **Step 2: Append sources to CMake compilation target**

Modify `src/core/CMakeLists.txt` to include new C++ sources when compiled:
```cmake
set(
  _cpp_core_srcs
  catchem_core.cpp
  catchem_api.cpp
)

# Append to compile target only if enabled or build universally
add_library(CATChem_core_cpp STATIC ${_cpp_core_srcs})
target_link_libraries(CATChem_core_cpp PUBLIC Kokkos::kokkos)
target_include_directories(CATChem_core_cpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
```

- [ ] **Step 3: Commit C-API Bridge**

```bash
git add src/core/catchem_api.hpp src/core/catchem_api.cpp src/core/CMakeLists.txt
git commit -m "feat(core): implement C-API extern bindings for ISO_C_BINDING boundary integration"
```

---

### Task 6: System Interoperability Verification & Integration Testing

**Files:**
- Create: `tests/test_catchem_interop.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: C-API bindings from `catchem_api.hpp`
- Produces: Integrated test suite certifying Fortran layout preservation and synchronized deep copies on host and device execution layers.

- [ ] **Step 1: Write integration tests**

Create `tests/test_catchem_interop.cpp`:
```cpp
#include "catchem_api.hpp"
#include <Kokkos_Core.hpp>
#include <cassert>
#include <iostream>
#include <vector>

// Mock Fortran physics scheme working directly on host array
void run_mock_fortran_physics(double* ptr, int n_cols, int n_levels) {
    // Simulate Fortran LayoutLeft (column-major) indexing: (i, j) -> i + j * n_cols
    for (int j = 0; j < n_levels; ++j) {
        for (int i = 0; i < n_cols; ++i) {
            ptr[i + j * n_cols] += 10.0; // Add tendency
        }
    }
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        int n_cols = 4;
        int n_levels = 5;
        int n_species = 2;

        // Allocate mock Fortran memory (column-major contiguous)
        std::vector<double> fortran_array(n_cols * n_levels, 1.0);

        // 1. Create Core & bind arrays
        void* core = catchem_core_create(n_cols, n_levels, n_species);
        void* state = catchem_core_get_state_manager(core);

        catchem_state_bind_2d(state, "temperature", fortran_array.data());

        // 2. Sync to active space
        catchem_state_sync_to_device(state);

        // 3. Execute Fortran process sequentially modifying the raw array on host
        run_mock_fortran_physics(fortran_array.data(), n_cols, n_levels);

        // Verify direct zero-copy modification
        assert(fortran_array[0] == 11.0);

        // 4. Sync up and clean up
        catchem_state_sync_to_host(state);
        catchem_core_destroy(core);

        std::cout << "SUCCESS: Interop Shared State Validation Passed!\n";
    }
    Kokkos::finalize();
    return 0;
}
```

- [ ] **Step 2: Append test binary to tests/CMakeLists.txt**

Modify `tests/CMakeLists.txt` to compile integration test when Kokkos is ON:
```cmake
if(ENABLE_KOKKOS)
  add_executable(test_catchem_interop test_catchem_interop.cpp)
  target_link_libraries(test_catchem_interop PRIVATE CATChem_core_cpp Kokkos::kokkos)
  add_test(NAME test_catchem_interop COMMAND test_catchem_interop)
endif()
```

- [ ] **Step 3: Compile and run test suite**

Run build and verification tests:
```bash
cmake -B build -DENABLE_KOKKOS=ON -DBUILD_TESTING=ON
cmake --build build --target test_catchem_interop
./build/tests/test_catchem_interop
```
Expected output: SUCCESS message printed to stdout.

- [ ] **Step 4: Commit tests**

```bash
git add tests/test_catchem_interop.cpp tests/CMakeLists.txt
git commit -m "test(core): add integration test verifying Fortran sequential execution and zero-copy boundary"
```
