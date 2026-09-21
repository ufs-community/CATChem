# Technical Design Specification: Phase 2 — Modernizing CATChem Processes & Creators to C++/Kokkos

## 1. Context & Objectives

In standard atmospheric physical modeling, physical processes (such as SeaSalt, Settling, Dry Deposition, and Wet Deposition) act as distinct operational nodes on the simulation timeline. Under the legacy CATChem framework, these nodes were governed by polymorphic Fortran objects (`ProcessInterface_Mod.F90`) and allocated sequentially via a Fortran factory registry.

Having successfully established a centralized C++ timestep orchestrator (`catchem::Core`), an unmanaged shared state manager (`catchem::StateManager`), and a dual-sync execution boundary (`catchem::FortranProcess`), we can now port the Process interfaces, factories, and operational nodes themselves directly to C++20/Kokkos.

### Core Objectives:
*   **C++ Execution Control:** Shift structural process interfaces, managers, and factories completely to modern C++20, keeping only the underlying physical computational schemes in Fortran if not yet ported.
*   **Dynamic Process Registry:** Establish a standard, extensible C++ `ProcessRegistry` mapping string identifiers to process instantiators.
*   **Unified Scheme Dispatch:** Ensure C++ process nodes can query runtime configurations (e.g. scheme name) and cleanly choose between highly parallel Kokkos C++ kernels or sequential legacy Fortran callbacks with zero-copy shared memory.

---

## 2. Architecture & Class Mappings

All process operations, from discovery to execution, are shifted to C++ in `src/core/`.

```
                    +------------------------------------+
                    |        catchem::Core Orchestrator  |
                    +------------------------------------+
                                       |
                                       v (sequential timeline list)
                    +------------------------------------+
                    |      catchem::ProcessInterface     | (Virtual Base)
                    +------------------------------------+
                        /              |               \
                       v               v                v
         +-------------------+  +--------------+  +-------------------+
         |  SettlingProcess  |  |SeaSaltProcess|  |   DryDepProcess   | (Concrete Classes)
         +-------------------+  +--------------+  +-------------------+
             /          \
            v (Kokkos)   v (Fortran Bridge)
      [C++ Kernel]     [Legacy Bridge Subroutine]
```

### 2.1 Abstract Base Process Interface
The abstract base interface governs lifecycle hooks, metadata, and parallel execution:

```cpp
// src/core/catchem_process_interface.hpp
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

---

## 3. Dynamic Creator & Process Registry

To avoid hardcoded creation switches and maintain high decouple constraints, we implement a dynamic creator registry inside `src/core/catchem_process_registry.hpp`:

```cpp
// src/core/catchem_process_registry.hpp
#pragma once
#include <string>
#include <memory>
#include <unordered_map>
#include <functional>
#include <stdexcept>
#include "catchem_process_interface.hpp"

namespace catchem {

using ProcessCreator = std::function<std::shared_ptr<ProcessInterface>()>;

class ProcessRegistry {
private:
    std::unordered_map<std::string, ProcessCreator> creators;

    // Singleton private constructor
    ProcessRegistry() = default;

public:
    static ProcessRegistry& get_instance() {
        static ProcessRegistry instance;
        return instance;
    }

    void register_process(const std::string& name, ProcessCreator creator) {
        creators[name] = creator;
    }

    bool has_process(const std::string& name) const {
        return creators.find(name) != creators.end();
    }

    std::shared_ptr<ProcessInterface> create(const std::string& name) {
        if (!has_process(name)) {
            throw std::invalid_argument("Process not registered in C++: " + name);
        }
        return creators.at(name)();
    }

    void clear() {
        creators.clear();
    }
};

} // namespace catchem
```

---

## 4. Concrete Physical Processes Implementation

Each physical process class reads configuration keys at initialization time and implements a scheme dispatch switch inside its `run()` loop.

### 4.1 `catchem::SettlingProcess` Class
Demonstrates how the C++ process handles both native Kokkos GPU/CPU targets and Fortran dynamic bridges side-by-side:

```cpp
// src/core/catchem_process_settling.hpp
#pragma once
#include "catchem_process_interface.hpp"
#include <functional>

namespace catchem {

class SettlingProcess : public ProcessInterface {
private:
    std::string scheme_name;
    std::function<void(void*)> fortran_callback;

public:
    SettlingProcess() : scheme_name("GOCART"), fortran_callback(nullptr) {}

    std::string get_name() const override {
        return "settling";
    }

    void init(std::shared_ptr<StateManager> state) override {
        // Read active scheme from configuration (Mocked as GOCART or Fortran bridge callback registration)
    }

    void set_fortran_bridge_callback(void (*callback)(void*)) {
        fortran_callback_ = callback;
    }

    void run(std::shared_ptr<StateManager> state) override {
        if (scheme_name == "kokkos_gocart") {
            // Native C++ Kokkos path: Run direct device loops
            execute_kokkos_settling_kernel(state);
        } else if (fortran_callback) {
            // Legacy Fortran bridge callback path: Sync memory sequentially and call Fortran
            state->sync_to_host();
            fortran_callback(static_cast<void*>(state.get()));
            state->sync_to_device();
        }
    }

    void finalize() override {}

private:
    void (*fortran_callback_)(void*) = nullptr;

    void execute_kokkos_settling_kernel(std::shared_ptr<StateManager> state) {
        // Parallel Kokkos kernel capture view by value
    }
};

} // namespace catchem
```

### 4.2 Application to All Processes
The same structural C++ class porting applies to remaining categories:
1.  **`catchem::SeaSaltProcess`:** Instantiates `GONG97`, `GONG03`, or `GEOS12` schemes. Calls C++ Kokkos parallel dispatchers, or falls back to legacy Fortran modules.
2.  **`catchem::DryDepProcess`:** Instantiates `WESELY`, `GOCART`, or `ZHANG` deposition schemes.
3.  **`catchem::WetDepProcess`:** Instantiates `JACOB` scavenging schemes.
4.  **`catchem::SO4chemProcess`:** Instantiates GOCART sulfate chemistry.

---

## 5. CMake List and Library Re-Organization

We cleanly update `src/core/CMakeLists.txt` to include the modernized C++ processes and drop legacy macro inclusions.

```cmake
# src/core/CMakeLists.txt
if(ENABLE_KOKKOS)
  set(
    _cpp_core_srcs
    catchem_core.cpp
    catchem_api.cpp
    catchem_diagnostic.cpp
    # Ported dynamic processes
    catchem_process_settling.cpp
    catchem_process_seasalt.cpp
    catchem_process_drydep.cpp
    catchem_process_wetdep.cpp
    catchem_process_so4chem.cpp
  )
  add_library(CATChem_core_cpp STATIC ${_cpp_core_srcs})
...
```

The old Fortran registries (`ProcessRegistry_Mod.F90`, `ProcessFactory_Mod.F90`) are removed, keeping files fully clean and maintainable.

---

## 6. Implementation and Verification Plan

We execute this modernization step-by-step:

*   **Task 1: Dynamic Creator and Registry**  
    Create `src/core/catchem_process_registry.hpp` and verify compilation.
*   **Task 2: Port Concrete Process Classes to C++**  
    Create and compile the class files:
    *   `catchem_process_settling.hpp` / `cpp`
    *   `catchem_process_seasalt.hpp` / `cpp`
    *   `catchem_process_drydep.hpp` / `cpp`
    *   `catchem_process_wetdep.hpp` / `cpp`
    *   `catchem_process_so4chem.hpp` / `cpp`
*   **Task 3: Integrate with central C++ orchestrator**  
    Modify `catchem::Core` to dynamically load registered processes based on keys.
*   **Task 4: Interop Validation Verification**  
    Configure and execute tests inside `cece-dev:latest` Docker verifying that registering a mixed list of native Kokkos and Fortran bridged processes executes sequentially and returns mathematically identical states.
