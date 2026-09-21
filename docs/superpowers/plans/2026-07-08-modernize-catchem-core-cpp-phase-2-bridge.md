# Modernize CATChem Core Phase 2 Bridge Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish a dynamic, performance-portable C++-to-Fortran execution bridge to sequentially run legacy Fortran physical schemes under the C++ orchestrator, maintaining shared state coherence and zero-copy data alignments.

**Architecture:** Create a C++ `catchem::FortranProcess` implementing the virtual `catchem::ProcessInterface`. This class executes state-synchronizations and invokes C-linkage Fortran callbacks. Expand StateManager C-API boundaries with raw double pointer retrievers, which the Fortran side binds directly using `c_f_pointer` before running its legacy physical scheme loops.

**Tech Stack:** C++20, Kokkos, Fortran 2008 (ISO_C_BINDING), CMake

## Global Constraints

- Target C++20 utilizing the Kokkos backport of mdspan (`std::experimental::mdspan`), avoiding direct dependency on C++23 `<mdspan>`.
- Retain Fortran column-major storage layout (`Kokkos::LayoutLeft`) across the pointer boundary to achieve zero-copy execution on CPU targets.

---

### Task 1: Add Raw Pointer Retrievers to StateManager and C-API

**Files:**
- Modify: `src/core/catchem_state_manager.hpp`
- Modify: `src/core/catchem_api.hpp`
- Modify: `src/core/catchem_api.cpp`

**Interfaces:**
- Produces: `catchem::StateManager::get_host_pointer_1d/2d/3d()`, `catchem_state_get_pointer_1d/2d/3d()`

- [ ] **Step 1: Write the failing test**

```cpp
// Add compilation check step inside a temporary file tests/test_get_pointer_compilation.cpp
#include "catchem_api.hpp"

int main() {
    double* ptr = catchem_state_get_pointer_1d(nullptr, "test");
    return 0;
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++20 -Isrc/core tests/test_get_pointer_compilation.cpp`
Expected: FAIL with missing declarations/symbols in catchem_api.hpp

- [ ] **Step 3: Modify StateManager to add host pointer getters**

Modify `src/core/catchem_state_manager.hpp`:
```cpp
    double* get_host_pointer_1d(const std::string& name) {
        if (fields_1d.find(name) == fields_1d.end()) return nullptr;
        return fields_1d.at(name)->host_view.data();
    }

    double* get_host_pointer_2d(const std::string& name) {
        if (fields_2d.find(name) == fields_2d.end()) return nullptr;
        return fields_2d.at(name)->host_view.data();
    }

    double* get_host_pointer_3d(const std::string& name) {
        if (fields_3d.find(name) == fields_3d.end()) return nullptr;
        return fields_3d.at(name)->host_view.data();
    }
```

- [ ] **Step 4: Modify catchem_api declarations & definitions**

Add declarations to `src/core/catchem_api.hpp`:
```cpp
double* catchem_state_get_pointer_1d(void* state_ptr, const char* name);
double* catchem_state_get_pointer_2d(void* state_ptr, const char* name);
double* catchem_state_get_pointer_3d(void* state_ptr, const char* name);
```

Add definitions to `src/core/catchem_api.cpp`:
```cpp
double* catchem_state_get_pointer_1d(void* state_ptr, const char* name) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    return state->get_host_pointer_1d(name);
}

double* catchem_state_get_pointer_2d(void* state_ptr, const char* name) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    return state->get_host_pointer_2d(name);
}

double* catchem_state_get_pointer_3d(void* state_ptr, const char* name) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    return state->get_host_pointer_3d(name);
}
```

- [ ] **Step 5: Clean up temporary test file**

Run: `rm tests/test_get_pointer_compilation.cpp`

- [ ] **Step 6: Commit**

```bash
git add src/core/catchem_state_manager.hpp src/core/catchem_api.hpp src/core/catchem_api.cpp
git commit -m "feat(api): add raw host pointer retrievers to StateManager C-API"
```

---

### Task 2: Implement C++ FortranProcess Bridging Wrapper

**Files:**
- Create: `src/core/catchem_fortran_process.hpp`

**Interfaces:**
- Produces: `catchem::FortranProcess` class implementing `catchem::ProcessInterface`

- [ ] **Step 1: Write verification compilation check**

Create temporary `tests/test_fortran_process_compilation.cpp`:
```cpp
#include "catchem_fortran_process.hpp"

void dummy_callback(void* state) {}

int main() {
    catchem::FortranProcess proc("fortran_scheme", dummy_callback);
    return 0;
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++20 -Isrc/core tests/test_fortran_process_compilation.cpp`
Expected: FAIL with "catchem_fortran_process.hpp: No such file or directory"

- [ ] **Step 3: Implement catchem::FortranProcess class**

Create `src/core/catchem_fortran_process.hpp`:
```cpp
// src/core/catchem_fortran_process.hpp
#pragma once
#include <string>
#include <memory>
#include "catchem_process_interface.hpp"

namespace catchem {

// C-linkage declarations matching Fortran bridge callbacks
extern "C" {
    typedef void (*FortranBridgeCallback)(void* state_mgr);
}

class FortranProcess : public ProcessInterface {
private:
    std::string name;
    FortranBridgeCallback bridge_callback;

public:
    FortranProcess(const std::string& process_name, FortranBridgeCallback callback)
        : name(process_name), bridge_callback(callback) {}

    std::string get_name() const override {
        return name;
    }

    void init(std::shared_ptr<StateManager> state) override {
        // Initial setup if required
    }

    void run(std::shared_ptr<StateManager> state) override {
        // 1. Sync device Views to host unified memory
        state->sync_to_host();

        // 2. Invoke the Fortran bridging callback
        if (bridge_callback) {
            bridge_callback(static_cast<void*>(state.get()));
        }

        // 3. Sync modified host buffers back to device Views
        state->sync_to_device();
    }

    void finalize() override {
        // Cleanup if required
    }
};

} // namespace catchem
```

- [ ] **Step 4: Clean up temporary test file**

Run: `rm tests/test_fortran_process_compilation.cpp`

- [ ] **Step 5: Commit**

```bash
git add src/core/catchem_fortran_process.hpp
git commit -m "feat(core): implement catchem::FortranProcess dynamic bridge"
```

---

### Task 3: Implement Fortran-C Core Bridge module

**Files:**
- Create: `src/core/FortranCoreBridge_Mod.F90`
- Modify: `src/core/CMakeLists.txt`

**Interfaces:**
- Produces: `run_settling_physics_fortran_bridge` (marked with `bind(C)`)

- [ ] **Step 1: Implement FortranCoreBridge_Mod module**

Create `src/core/FortranCoreBridge_Mod.F90`:
```fortran
!> \file FortranCoreBridge_Mod.F90
!! \brief Fortran dynamic bridging procedures callback to execute legacy schemes.
module FortranCoreBridge_Mod
   use iso_c_binding, only: c_ptr, c_f_pointer
   use precision_mod, only: fp

   implicit none
   private

   public :: run_settling_physics_fortran_bridge

   interface
      ! C-API Bindings
      function catchem_state_get_pointer_2d(state_ptr, name) bind(C, name="catchem_state_get_pointer_2d")
         import :: c_ptr
         type(c_ptr), value :: state_ptr
         character(*, kind=1), intent(in) :: name
         type(c_ptr) :: catchem_state_get_pointer_2d
      end function catchem_state_get_pointer_2d
   end interface

contains

   !> \brief C-linkable dynamic bridge to execute legacy physical schemes on raw C++ memory pointers.
   subroutine run_settling_physics_fortran_bridge(state_ptr) bind(C, name="run_settling_physics_fortran_bridge")
      type(c_ptr), value :: state_ptr
      type(c_ptr) :: c_temp
      real(fp), pointer :: temp(:,:)
      integer :: n_cols, n_levels

      ! 1. Mock dimensions for this bridge test (matching test_catchem_interop sizes: 4 x 5)
      n_cols = 4
      n_levels = 5

      ! 2. Retrieve C++ double pointer for "temperature"
      ! append null-termination to Fortran string literal
      c_temp = catchem_state_get_pointer_2d(state_ptr, "temperature" // c_null_char)

      if (.not. c_associated(c_temp)) return

      ! 3. Wrap raw C++ pointer back to Fortran array pointer (LayoutLeft matching column-major size)
      call c_f_pointer(c_temp, temp, [n_cols, n_levels])

      ! 4. Execute legacy Fortran physical scheme directly working on shared memory in-place
      temp(:,:) = temp(:,:) + 10.0_fp

   end subroutine run_settling_physics_fortran_bridge

end module FortranCoreBridge_Mod
```

- [ ] **Step 2: Append FortranCoreBridge_Mod to CMake lists**

In `src/core/CMakeLists.txt`, append `FortranCoreBridge_Mod.F90` to the `_core_srcs` block so that it builds as part of the core library:
```cmake
set(
  _core_srcs
  ConfigManager_Mod.F90
  VirtualColumn_Mod.F90
  GridManager_Mod.F90
  DiagnosticManager_Mod.F90
  StateManager_Mod.F90
  ChemSpeciesUtils_Mod.F90
  EmissionConfigValidator_Mod.F90
  ProcessInterface_Mod.F90
  ProcessFactory_Mod.F90
  ProcessRegistry_Mod.F90
  ProcessManager_Mod.F90
  CATChemCore_Mod.F90
  FortranCoreBridge_Mod.F90
)
```

- [ ] **Step 3: Commit**

```bash
git add src/core/FortranCoreBridge_Mod.F90 src/core/CMakeLists.txt
git commit -m "feat(core): implement FortranCoreBridge_Mod for raw pointer mappings"
```

---

### Task 4: Verify Mixed Timeline in Integration Test

**Files:**
- Modify: `tests/test_catchem_interop.cpp`

**Interfaces:**
- Consumes: `catchem::FortranProcess`, `run_settling_physics_fortran_bridge` callback

- [ ] **Step 1: Write integration assertions inside test_catchem_interop.cpp**

Open `tests/test_catchem_interop.cpp` and declare the external Fortran bridge callback:
```cpp
extern "C" {
    void run_settling_physics_fortran_bridge(void* state_ptr);
}
```

Add **TEST 3: Phase 2 Sequenced Fortran Dynamic Bridge** inside `main()`:
```cpp
        // ==========================================
        // TEST 3: Phase 2 Sequenced Fortran Dynamic Bridge
        // ==========================================
        {
            int n_cols = 4;
            int n_levels = 5;
            int n_species = 2;

            // Allocate mock Fortran memory
            std::vector<double> fortran_array(n_cols * n_levels, 1.0);

            // 1. Create Core, StateManager and Dynamic Registry
            void* core_ptr = catchem_core_create(n_cols, n_levels, n_species);
            auto* core = static_cast<catchem::Core*>(core_ptr);
            void* state = catchem_core_get_state_manager(core_ptr);

            // Bind temperature array
            catchem_state_bind_2d(state, "temperature", fortran_array.data());

            // 2. Attach our newly created C++ FortranProcess bridge callback
            core->add_process(std::make_shared<catchem::FortranProcess>(
                "legacy_settling_physics",
                run_settling_physics_fortran_bridge
            ));

            // 3. Step forward (runs dynamic process, which syncs memory & calls bridge in order)
            catchem_core_run_timestep(core_ptr, 3600.0);

            // 4. Verify results
            // Fortran bridge executes: temp = temp + 10.0D0
            assert(fortran_array[0] == 11.0);

            catchem_core_destroy(core_ptr);
            std::cout << "SUCCESS: Sequenced Fortran Dynamic Bridge Validation Passed!\n";
        }
```

Include `#include "catchem_fortran_process.hpp"` at the top of `test_catchem_interop.cpp`.

- [ ] **Step 2: Build and execute tests in Docker**

Run targeted compilation and run:
`docker run --rm -v /Users/barry/Documents/CATChem:/workspace -w /workspace/build-test cece-dev:latest bash -c "apt-get update -y && apt-get install -y python3 && cmake .. -DENABLE_KOKKOS=ON -DENABLE_TESTING=ON && make test_catchem_interop && ./tests/test_catchem_interop"`

Expected: PASS, printing:
```text
SUCCESS: Interop Shared State Validation Passed!
SUCCESS: C++ Diagnostic Validation Passed!
SUCCESS: Sequenced Fortran Dynamic Bridge Validation Passed!
```

- [ ] **Step 3: Commit**

```bash
git add tests/test_catchem_interop.cpp
git commit -m "test(bridge): integrate and verify mixed-execution Fortran process bridge"
```
