# Modernize Processes & Creators Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the Process interfaces and Process Creators for all physical processes to C++20, and establish a dynamic `ProcessRegistry` in the C++ orchestrator, keeping only the execution schemes themselves in Fortran.

**Architecture:** Create an abstract virtual `catchem::ProcessInterface`. Implement concrete C++ processes (`SettlingProcess`, `SeaSaltProcess`, `DryDepProcess`, `WetDepProcess`, `SO4chemProcess`) inheriting from it. Establish a singleton `catchem::ProcessRegistry` mapping process string names to creator functions.

**Tech Stack:** C++20, Kokkos, ISO_C_BINDING, CMake

## Global Constraints

- Target C++20 utilizing the Kokkos backport of mdspan (`std::experimental::mdspan`), avoiding direct dependency on C++23 `<mdspan>`.
- Retain Fortran column-major storage layout (`Kokkos::LayoutLeft`) across the pointer boundary to achieve zero-copy execution on CPU targets.

---

### Task 1: Create Dynamic C++ Creator & Registry

**Files:**
- Create: `src/core/catchem_process_registry.hpp`

**Interfaces:**
- Produces: `catchem::ProcessRegistry`, `catchem::ProcessCreator`

- [ ] **Step 1: Write a temporary compilation test**

```cpp
// Create tests/test_registry_compilation.cpp
#include "catchem_process_registry.hpp"
#include <iostream>

class TestProc : public catchem::ProcessInterface {
public:
    std::string get_name() const override { return "test_proc"; }
    void init(std::shared_ptr<catchem::StateManager> state) override {}
    void run(std::shared_ptr<catchem::StateManager> state) override {}
    void finalize() override {}
};

int main() {
    auto& reg = catchem::ProcessRegistry::get_instance();
    reg.register_process("test", []() { return std::make_shared<TestProc>(); });
    if (reg.has_process("test")) {
        std::cout << "SUCCESS\n";
    }
    return 0;
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++20 -Isrc/core tests/test_registry_compilation.cpp`
Expected: FAIL with missing `catchem_process_registry.hpp`

- [ ] **Step 3: Implement catchem::ProcessRegistry**

Create `src/core/catchem_process_registry.hpp`:
```cpp
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

- [ ] **Step 4: Verify test compilation passes**

Run: `g++ -std=c++20 -Isrc/core tests/test_registry_compilation.cpp && ./a.out`
Expected: SUCCESS

- [ ] **Step 5: Clean up temporary test file**

Run: `rm -f tests/test_registry_compilation.cpp a.out`

- [ ] **Step 6: Commit**

```bash
git add src/core/catchem_process_registry.hpp
git commit -m "feat(core): implement C++ ProcessRegistry for dynamic creators"
```

---

### Task 2: Port Concrete Process Classes to C++

**Files:**
- Create: `src/core/catchem_process_settling.hpp` / `cpp`
- Create: `src/core/catchem_process_seasalt.hpp` / `cpp`
- Create: `src/core/catchem_process_drydep.hpp` / `cpp`
- Create: `src/core/catchem_process_wetdep.hpp` / `cpp`
- Create: `src/core/catchem_process_so4chem.hpp` / `cpp`
- Modify: `src/core/CMakeLists.txt`

**Interfaces:**
- Produces: `catchem::SettlingProcess`, `catchem::SeaSaltProcess`, `catchem::DryDepProcess`, `catchem::WetDepProcess`, `catchem::SO4chemProcess` classes

- [ ] **Step 1: Write header class files for Settling, SeaSalt, DryDep, WetDep, and SO4chem**

Create `src/core/catchem_process_settling.hpp`:
```cpp
#pragma once
#include "catchem_process_interface.hpp"
#include <functional>

namespace catchem {

class SettlingProcess : public ProcessInterface {
private:
    std::string active_scheme;
    std::function<void(void*)> fortran_callback;

public:
    SettlingProcess();
    std::string get_name() const override { return "settling"; }
    void init(std::shared_ptr<StateManager> state) override;
    void set_fortran_bridge_callback(std::function<void(void*)> cb);
    void run(std::shared_ptr<StateManager> state) override;
    void finalize() override {}
};

} // namespace catchem
```

Create `src/core/catchem_process_seasalt.hpp`:
```cpp
#pragma once
#include "catchem_process_interface.hpp"
#include <functional>

namespace catchem {

class SeaSaltProcess : public ProcessInterface {
private:
    std::string active_scheme;
    std::function<void(void*)> fortran_callback;

public:
    SeaSaltProcess();
    std::string get_name() const override { return "seasalt"; }
    void init(std::shared_ptr<StateManager> state) override;
    void set_fortran_bridge_callback(std::function<void(void*)> cb);
    void run(std::shared_ptr<StateManager> state) override;
    void finalize() override {}
};

} // namespace catchem
```

Create `src/core/catchem_process_drydep.hpp`:
```cpp
#pragma once
#include "catchem_process_interface.hpp"
#include <functional>

namespace catchem {

class DryDepProcess : public ProcessInterface {
private:
    std::string active_scheme;
    std::function<void(void*)> fortran_callback;

public:
    DryDepProcess();
    std::string get_name() const override { return "drydep"; }
    void init(std::shared_ptr<StateManager> state) override;
    void set_fortran_bridge_callback(std::function<void(void*)> cb);
    void run(std::shared_ptr<StateManager> state) override;
    void finalize() override {}
};

} // namespace catchem
```

Create `src/core/catchem_process_wetdep.hpp`:
```cpp
#pragma once
#include "catchem_process_interface.hpp"
#include <functional>

namespace catchem {

class WetDepProcess : public ProcessInterface {
private:
    std::string active_scheme;
    std::function<void(void*)> fortran_callback;

public:
    WetDepProcess();
    std::string get_name() const override { return "wetdep"; }
    void init(std::shared_ptr<StateManager> state) override;
    void set_fortran_bridge_callback(std::function<void(void*)> cb);
    void run(std::shared_ptr<StateManager> state) override;
    void finalize() override {}
};

} // namespace catchem
```

Create `src/core/catchem_process_so4chem.hpp`:
```cpp
#pragma once
#include "catchem_process_interface.hpp"
#include <functional>

namespace catchem {

class SO4chemProcess : public ProcessInterface {
private:
    std::string active_scheme;
    std::function<void(void*)> fortran_callback;

public:
    SO4chemProcess();
    std::string get_name() const override { return "so4chem"; }
    void init(std::shared_ptr<StateManager> state) override;
    void set_fortran_bridge_callback(std::function<void(void*)> cb);
    void run(std::shared_ptr<StateManager> state) override;
    void finalize() override {}
};

} // namespace catchem
```

- [ ] **Step 2: Write concrete C++ process implementations**

Create `src/core/catchem_process_settling.cpp`:
```cpp
#include "catchem_process_settling.hpp"

namespace catchem {

SettlingProcess::SettlingProcess() : active_scheme("legacy_fortran"), fortran_callback(nullptr) {}

void SettlingProcess::init(std::shared_ptr<StateManager> state) {
    // Register self in global registry on initialization
}

void SettlingProcess::set_fortran_bridge_callback(std::function<void(void*)> cb) {
    fortran_callback = cb;
}

void SettlingProcess::run(std::shared_ptr<StateManager> state) {
    if (fortran_callback) {
        state->sync_to_host();
        fortran_callback(static_cast<void*>(state.get()));
        state->sync_to_device();
    }
}

} // namespace catchem
```

Create `src/core/catchem_process_seasalt.cpp`:
```cpp
#include "catchem_process_seasalt.hpp"

namespace catchem {

SeaSaltProcess::SeaSaltProcess() : active_scheme("legacy_fortran"), fortran_callback(nullptr) {}

void SeaSaltProcess::init(std::shared_ptr<StateManager> state) {}

void SeaSaltProcess::set_fortran_bridge_callback(std::function<void(void*)> cb) {
    fortran_callback = cb;
}

void SeaSaltProcess::run(std::shared_ptr<StateManager> state) {
    if (fortran_callback) {
        state->sync_to_host();
        fortran_callback(static_cast<void*>(state.get()));
        state->sync_to_device();
    }
}

} // namespace catchem
```

Create `src/core/catchem_process_drydep.cpp`:
```cpp
#include "catchem_process_drydep.hpp"

namespace catchem {

DryDepProcess::DryDepProcess() : active_scheme("legacy_fortran"), fortran_callback(nullptr) {}

void DryDepProcess::init(std::shared_ptr<StateManager> state) {}

void DryDepProcess::set_fortran_bridge_callback(std::function<void(void*)> cb) {
    fortran_callback = cb;
}

void DryDepProcess::run(std::shared_ptr<StateManager> state) {
    if (fortran_callback) {
        state->sync_to_host();
        fortran_callback(static_cast<void*>(state.get()));
        state->sync_to_device();
    }
}

} // namespace catchem
```

Create `src/core/catchem_process_wetdep.cpp`:
```cpp
#include "catchem_process_wetdep.hpp"

namespace catchem {

WetDepProcess::WetDepProcess() : active_scheme("legacy_fortran"), fortran_callback(nullptr) {}

void WetDepProcess::init(std::shared_ptr<StateManager> state) {}

void WetDepProcess::set_fortran_bridge_callback(std::function<void(void*)> cb) {
    fortran_callback = cb;
}

void WetDepProcess::run(std::shared_ptr<StateManager> state) {
    if (fortran_callback) {
        state->sync_to_host();
        fortran_callback(static_cast<void*>(state.get()));
        state->sync_to_device();
    }
}

} // namespace catchem
```

Create `src/core/catchem_process_so4chem.cpp`:
```cpp
#include "catchem_process_so4chem.hpp"

namespace catchem {

SO4chemProcess::SO4chemProcess() : active_scheme("legacy_fortran"), fortran_callback(nullptr) {}

void SO4chemProcess::init(std::shared_ptr<StateManager> state) {}

void SO4chemProcess::set_fortran_bridge_callback(std::function<void(void*)> cb) {
    fortran_callback = cb;
}

void SO4chemProcess::run(std::shared_ptr<StateManager> state) {
    if (fortran_callback) {
        state->sync_to_host();
        fortran_callback(static_cast<void*>(state.get()));
        state->sync_to_device();
    }
}

} // namespace catchem
```

- [ ] **Step 3: Modify `src/core/CMakeLists.txt` to include newly ported processes**

Add source files to the static library C++ sources inside `src/core/CMakeLists.txt`:
```cmake
  set(
    _cpp_core_srcs
    catchem_core.cpp
    catchem_api.cpp
    catchem_diagnostic.cpp
    catchem_process_settling.cpp
    catchem_process_seasalt.cpp
    catchem_process_drydep.cpp
    catchem_process_wetdep.cpp
    catchem_process_so4chem.cpp
  )
```

- [ ] **Step 4: Commit**

```bash
git add src/core/catchem_process_settling.* src/core/catchem_process_seasalt.* src/core/catchem_process_drydep.* src/core/catchem_process_wetdep.* src/core/catchem_process_so4chem.* src/core/CMakeLists.txt
git commit -m "feat(core): implement C++ process wrappers for physical schemes"
```

---

### Task 3: Integrate ProcessRegistry with central Core Orchestrator

**Files:**
- Modify: `src/core/catchem_core.hpp`
- Modify: `src/core/catchem_core.cpp`
- Modify: `src/core/catchem_api.hpp`
- Modify: `src/core/catchem_api.cpp`

**Interfaces:**
- Produces: `catchem_core_register_process`, `catchem_core_register_fortran_callback` in C-API layer

- [ ] **Step 1: Modify `catchem_core.hpp` to instantiate processes via ProcessRegistry**

Modify `src/core/catchem_core.hpp`:
```cpp
// Add public methods
    void register_process_cpp(const std::string& name);
    void register_fortran_callback(const std::string& name, void (*callback)(void*));
```

- [ ] **Step 2: Modify `catchem_core.cpp` to implement creation and callback bridges**

Modify `src/core/catchem_core.cpp`:
```cpp
#include "catchem_core.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_process_settling.hpp"
#include "catchem_process_seasalt.hpp"
#include "catchem_process_drydep.hpp"
#include "catchem_process_wetdep.hpp"
#include "catchem_process_so4chem.hpp"

namespace catchem {

Core::Core(int nc, int nl, int ns) {
    state_mgr = std::make_shared<StateManager>(nc, nl, ns);
    diag_mgr = std::make_shared<DiagnosticManager>();

    // Register physical processes into the registry
    auto& reg = ProcessRegistry::get_instance();
    reg.register_process("settling", []() { return std::make_shared<SettlingProcess>(); });
    reg.register_process("seasalt", []() { return std::make_shared<SeaSaltProcess>(); });
    reg.register_process("drydep", []() { return std::make_shared<DryDepProcess>(); });
    reg.register_process("wetdep", []() { return std::make_shared<WetDepProcess>(); });
    reg.register_process("so4chem", []() { return std::make_shared<SO4chemProcess>(); });
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

void Core::register_process_cpp(const std::string& name) {
    add_process(ProcessRegistry::get_instance().create(name));
}

void Core::register_fortran_callback(const std::string& name, void (*callback)(void*)) {
    for (auto& process : processes) {
        if (process->get_name() == name) {
            if (name == "settling") {
                static_cast<SettlingProcess*>(process.get())->set_fortran_bridge_callback(callback);
            } else if (name == "seasalt") {
                static_cast<SeaSaltProcess*>(process.get())->set_fortran_bridge_callback(callback);
            } else if (name == "drydep") {
                static_cast<DryDepProcess*>(process.get())->set_fortran_bridge_callback(callback);
            } else if (name == "wetdep") {
                static_cast<WetDepProcess*>(process.get())->set_fortran_bridge_callback(callback);
            } else if (name == "so4chem") {
                static_cast<SO4chemProcess*>(process.get())->set_fortran_bridge_callback(callback);
            }
            break;
        }
    }
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

- [ ] **Step 3: Expose registration endpoints in the C-API boundary**

Modify `src/core/catchem_api.hpp`:
```cpp
void catchem_core_register_process(void* core_ptr, const char* name);
void catchem_core_register_fortran_callback(void* core_ptr, const char* name, void (*callback)(void*));
```

Modify `src/core/catchem_api.cpp`:
```cpp
void catchem_core_register_process(void* core_ptr, const char* name) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    core->register_process_cpp(name);
}

void catchem_core_register_fortran_callback(void* core_ptr, const char* name, void (*callback)(void*)) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    core->register_fortran_callback(name, callback);
}
```

- [ ] **Step 4: Commit**

```bash
git add src/core/catchem_core.* src/core/catchem_api.*
git commit -m "feat(api): expose C-API endpoints for registering C++ processes and Fortran callbacks"
```

---

### Task 4: Verify Full mixed C++ process lifecycle in interop tests

**Files:**
- Modify: `tests/test_catchem_interop.cpp`

**Interfaces:**
- Consumes: C-API process registries and callbacks

- [ ] **Step 1: Write verification integration assertions**

Add **TEST 4: Phase 3 Unified C++ Process and Creator Registry** sequentially inside `tests/test_catchem_interop.cpp`:
```cpp
        // ==========================================
        // TEST 4: Phase 3 Unified C++ Process and Creator Registry
        // ==========================================
        {
            int n_cols = 4;
            int n_levels = 5;
            int n_species = 2;

            // Allocate mock Fortran memory
            std::vector<double> fortran_array(n_cols * n_levels, 1.0);

            // 1. Create Core Orchestrator
            void* core_ptr = catchem_core_create(n_cols, n_levels, n_species);
            void* state = catchem_core_get_state_manager(core_ptr);

            // Bind temperature array
            catchem_state_bind_2d(state, "temperature", fortran_array.data());

            // 2. Register Process dynamically through Creator registry and C-API
            catchem_core_register_process(core_ptr, "settling");

            // 3. Register Fortran bridge callback dynamically through C-API
            catchem_core_register_fortran_callback(core_ptr, "settling", run_settling_physics_fortran_bridge);

            // 4. Run central timestep
            catchem_core_run_timestep(core_ptr, 3600.0);

            // 5. Assert sequential update in shared memory
            assert(fortran_array[0] == 11.0);

            catchem_core_destroy(core_ptr);
            std::cout << "SUCCESS: Unified C++ Process and Creator Registry Validation Passed!\n";
        }
```

- [ ] **Step 2: Build and run the entire interop test suite**

Run inside Docker:
`docker run --rm -v /Users/barry/Documents/CATChem:/workspace -w /workspace/build-test cece-dev:latest bash -c "apt-get update -y && apt-get install -y python3 && cmake .. -DENABLE_KOKKOS=ON -DENABLE_TESTING=ON && make test_catchem_interop && ./tests/test_catchem_interop"`

Expected: PASS, printing:
```text
SUCCESS: Interop Shared State Validation Passed!
SUCCESS: C++ Diagnostic Validation Passed!
SUCCESS: Sequenced Fortran Dynamic Bridge Validation Passed!
SUCCESS: Unified C++ Process and Creator Registry Validation Passed!
```

- [ ] **Step 3: Commit**

```bash
git add tests/test_catchem_interop.cpp
git commit -m "test(interop): verify unified C++ process registry and Dynamic Creator lifecycle"
```
