# Technical Design: Phase 2 — Interfacing C++ with Fortran Processes

**Date:** July 8, 2026  
**Status:** Approved  
**Topic:** Establishing a dynamic, dual-execution bridge allowing the modernized C++ core orchestrator to seamlessly invoke legacy Fortran physical schemes on CPU, maintaining shared state coherence and CCPP boundaries.

---

## 1. Executive Summary & Architectural Goals

In the Strangler Fig migration pattern, CATChem physical schemes (e.g., Settling, Dry Deposition, Wet Deposition) will be migrated to native Kokkos kernels step-by-step. To maintain a functional build at each step, legacy Fortran physical schemes must run seamlessly under the centralized C++ execution loop (`catchem::Core`).

### Primary Goals:
* **Option A: Centralized C++ Control:** The C++ core retains absolute ownership of the execution sequence, timelines, and memory states.
* **No Legacy Code Intrusion:** Original Fortran subroutines (with their `.meta` files and argument signatures) must remain completely unmodified. No changes to CCPP meta tables.
* **Shared Memory Coherence:** Synchronize memory automatically across compute devices. CPU executions automatically use zero-copy wrappers, while GPU device mirrors sync host pointers dynamically before Fortran callbacks.

---

## 2. C++ Dynamic Bridging Interface (`catchem::FortranProcess`)

To wrap legacy Fortran subroutines in the C++ process pipeline, we introduce a generic `FortranProcess` class implementing `catchem::ProcessInterface`.

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
        // Initial setup for the Fortran process if required
    }

    void run(std::shared_ptr<StateManager> state) override {
        // 1. Synchronize device views to host to ensure Fortran CPU process gets latest data
        state->sync_to_host();

        // 2. Call the C-linked Fortran bridge callback passing the StateManager pointer
        if (bridge_callback) {
            bridge_callback(static_cast<void*>(state.get()));
        }

        // 3. Post-execution: Sync host views to device so subsequent Kokkos C++ kernels see the updates
        state->sync_to_device();
    }

    void finalize() override {
        // Cleanup resources
    }
};

} // namespace catchem
```

---

## 3. Extension of C-API State Query Interfaces

The Fortran bridging function needs to retrieve raw memory pointers from the `catchem::StateManager` dynamically using string IDs. We expand the C-API boundary (`src/core/catchem_api.hpp`) with pointer extractor endpoints:

```cpp
// src/core/catchem_api.hpp (Additions)
#ifdef __cplusplus
extern "C" {
#endif

// Returns a direct raw double pointer to the host memory bound under 'name'
double* catchem_state_get_pointer_1d(void* state_ptr, const char* name);
double* catchem_state_get_pointer_2d(void* state_ptr, const char* name);
double* catchem_state_get_pointer_3d(void* state_ptr, const char* name);

#ifdef __cplusplus
}
#endif
```

```cpp
// src/core/catchem_api.cpp (Additions)
#include "catchem_api.hpp"
#include "catchem_state_manager.hpp"

extern "C" {

double* catchem_state_get_pointer_1d(void* state_ptr, const char* name) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    auto it = state->fields_1d.find(name);
    if (it != state->fields_1d.end()) {
        return it->second->host_view.data();
    }
    return nullptr;
}

double* catchem_state_get_pointer_2d(void* state_ptr, const char* name) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    auto it = state->fields_2d.find(name);
    if (it != state->fields_2d.end()) {
        return it->second->host_view.data();
    }
    return nullptr;
}

double* catchem_state_get_pointer_3d(void* state_ptr, const char* name) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    auto it = state->fields_3d.find(name);
    if (it != state->fields_3d.end()) {
        return it->second->host_view.data();
    }
    return nullptr;
}

}
```

---

## 4. Fortran Dynamic Slicing & Association Layer

To execute legacy physics subroutines that expect standard multidimensional pointers (e.g. `double precision, intent(inout) :: temp(:,:)`), the bridging routine uses `c_f_pointer` to safely map the raw pointers into standard Fortran shapes:

```fortran
! src/kokkos/KokkosDispatch_Mod.F90 (Additions or thin wrapper module)
module FortranCoreBridge_Mod
   use iso_c_binding
   implicit none
   private

   public :: run_settling_physics_fortran_bridge

   ! Interfaces to C-API getters
   interface
      function catchem_state_get_pointer_2d(state_ptr, name) bind(C, name="catchem_state_get_pointer_2d") result(ptr)
         import :: c_ptr, c_char
         type(c_ptr), value :: state_ptr
         character(kind=c_char), intent(in) :: name(*)
         type(c_ptr) :: ptr
      end function catchem_state_get_pointer_2d
   end interface

contains

   ! C-linked bridging callback
   subroutine run_settling_physics_fortran_bridge(state_ptr) bind(C, name="run_settling_physics_fortran_bridge")
      type(c_ptr), value :: state_ptr

      type(c_ptr) :: c_temp_ptr
      real(c_double), pointer :: fortran_temp(:,:)
      integer :: n_cols, n_levels

      ! 1. Query size bounds (can be queried via separate getters or stored in a common module)
      n_cols = 4
      n_levels = 5

      ! 2. Retrieve the raw C pointer for "temperature"
      c_temp_ptr = catchem_state_get_pointer_2d(state_ptr, "temperature" // c_null_char)

      if (c_associated(c_temp_ptr)) then
         ! 3. Associate raw pointer to Fortran multidimensional pointer (column-major)
         call c_f_pointer(c_temp_ptr, fortran_temp, [n_cols, n_levels])

         ! 4. Call the unchanged legacy physics subroutine
         call legacy_settling_physics_run(fortran_temp, n_cols, n_levels)
      end if
   end subroutine run_settling_physics_fortran_bridge

   ! Example legacy physics representation matching the exact original signatures
   subroutine legacy_settling_physics_run(temp, n_cols, n_levels)
      integer, intent(in) :: n_cols, n_levels
      real(c_double), intent(inout) :: temp(n_cols, n_levels)

      ! Sequentially modify values
      temp(:,:) = temp(:,:) + 10.0D0
   end subroutine legacy_settling_physics_run

end module FortranCoreBridge_Mod
```

---

## 5. Sequential Execution Timeline & Shared Coherence

```
 [ C++ Orchestration Loop ]                [ InteropField ]            [ Fortran Bridge ]
           |                                       |                            |
           v                                       |                            |
 1. fortran_proc->run()                            |                            |
    state->sync_to_host() ------------------------>|                            |
                                                   | (deep_copy if GPU)         |
                                                   |<----------------           |
 2. bridge_callback(state) ---------------------------------------------------->|
                                                                                | catchem_state_get_pointer()
                                                                                | c_f_pointer association
                                                                                | Run Fortran physics
                                                                                | (modifies host memory in place)
                                                                                |<---------
 3. state->sync_to_device() ---------------------->|                            |
                                                   | (deep_copy to GPU)         |
                                                   |<----------------           |
           v                                       |                            |
 [ C++ Kokkos Kernels Continue ]                   |                            |
```

---

## 6. Implementation & Porting Roadmap

To establish this interop framework, we define a highly targeted, 4-task execution plan:
1. **Task 1:** Expand C-API pointer retrievers (`catchem_state_get_pointer_*d`) inside C++ StateManager.
2. **Task 2:** Coder `catchem::FortranProcess` bridging wrapper in C++.
3. **Task 3:** Create `FortranCoreBridge_Mod.F90` to handle raw pointer mappings and `c_f_pointer` shapes.
4. **Task 4:** Extend `tests/test_catchem_interop.cpp` to register and run a mixed timeline (C++ and Fortran processes sequentially), verifying numerical coherence.
