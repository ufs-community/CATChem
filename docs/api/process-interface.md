# Process Interface API

This section covers the modernized C++ process interface APIs that define how atmospheric physics and chemistry processes integrate with the centralized CATChem framework.

## Overview

Under the modernized C++20 architecture, the Process Interface system provides:

- **catchem::ProcessInterface**: Virtual base class defining execution lifecycles for all processes.
- **catchem::ProcessRegistry**: Singleton registry utilizing creator lambdas for dynamic process instantiation.
- **catchem::FortranProcess**: Bridging adapter allowing legacy Fortran schemes to run under the C++ orchestrator.
- **Linker-Safe Registration**: BIND(C) callbacks preventing dead-code stripping of static process libraries.

---

## Core Interfaces

### catchem::ProcessInterface

The abstract C++ virtual base class and unified process name identifiers that all native processes implement:

```cpp
#pragma once
#include "catchem_state_manager.hpp"
#include <memory>
#include <string>
#include <string_view>

namespace catchem {

    // Centralized process identifier constants
    struct ProcessNames {
        static constexpr std::string_view GasChem = "gaschem";
        static constexpr std::string_view Photolysis = "photolysis";
        static constexpr std::string_view Settling = "settling";
        static constexpr std::string_view Dust = "dust";
        static constexpr std::string_view SeaSalt = "seasalt";
        static constexpr std::string_view CarbChem = "carbchem";
        static constexpr std::string_view SO4Chem = "so4chem";
    };

    class ProcessInterface {
    public:
        virtual ~ProcessInterface() = default;

        // Returns the lowercase string identifier of the process (using catchem::ProcessNames)
        virtual std::string get_name() const = 0;

        // Allocates local resources and binds config parameters
        virtual void init(std::shared_ptr<StateManager> state) = 0;

        // Executes physical/chemical calculations
        virtual void run(std::shared_ptr<StateManager> state) = 0;

        // Cleans up heap memory and releases resources
        virtual void finalize() = 0;
    };

} // namespace catchem
```

### Column Processing & Kokkos Parallelization

Unlike the legacy architecture which required a separate `ColumnProcessInterface` class, column virtualization and 1D column loops are now handled natively inside `run()` using zero-copy **Kokkos subviews** and parallel execution kernels:

```cpp
void MyProcess::run(std::shared_ptr<StateManager> state) {
    auto n_cols = state->n_cols;
    auto n_levels = state->n_levels;

    // Parallel-for over 1D columns using Kokkos
    Kokkos::parallel_for("MyProcessLoop", Kokkos::RangePolicy<Kokkos::HostSpace>(0, n_cols),
        [=](const int icol) {
            // Slice 3D fields to 1D columns with zero copy
            auto col_temp = Kokkos::subview(state->met.temp, icol, Kokkos::ALL(), 0);

            for (int k = 0; k < n_levels; ++k) {
                // Perform levels calculations
                process_level(col_temp(k));
            }
        });
}
```

---

## Process Registry & Factory Pattern

### ProcessRegistry

The dynamic creator factory maps lowercase process name strings to dynamic creator lambda functions:

```cpp
#pragma once
#include "catchem_process_interface.hpp"
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

namespace catchem {

    using ProcessCreator = std::function<std::shared_ptr<ProcessInterface>()>;

    class ProcessRegistry {
    private:
        std::unordered_map<std::string, ProcessCreator> creators;
        ProcessRegistry() = default;
    public:
        static ProcessRegistry& get_instance() {
            static ProcessRegistry instance;
            return instance;
        }

        void register_process(const std::string& name, ProcessCreator creator) {
            creators[name] = creator;
        }

        std::shared_ptr<ProcessInterface> create(const std::string& name) {
            auto it = creators.find(name);
            if (it == creators.end()) {
                throw std::runtime_error("Process not found in registry: " + name);
            }
            return it->second();
        }
    };

} // namespace catchem
```

---

## Linker-Safe Callback Registration

When physical processes are compiled into distinct static libraries (e.g. `catchem_process_gaschem`), linker optimizations often strip "unused" C++ symbols if they aren't explicitly referenced during main application linking.

To guarantee that creator lambdas are preserved and registered at startup, CATChem implements an explicit, linker-safe dynamic callback mechanism across the mixed-language boundary:

### 1. C++ Registration Hook
The process defines and exports a flat `extern "C"` registration hook:

```cpp
#include "catchem_process_gaschem.hpp"
#include "catchem_process_registry.hpp"

extern "C" void catchem_register_gaschem_cpp() {
    catchem::ProcessRegistry::get_instance().register_process("gaschem", []() {
        return std::make_shared<catchem::GasChemProcess>();
    });
}
```

### 2. Fortran Trigger
The process Fortran initialization wrapper explicitly triggers this dynamic hook, forcing the compiler and linker to preserve the C++ object files:

```fortran
module ProcessGasChemInterface_Mod
    use iso_c_binding
    implicit none

    interface
        subroutine catchem_register_gaschem_cpp() bind(C, name="catchem_register_gaschem_cpp")
        end subroutine
    end interface

contains

    subroutine init_gaschem(rc)
        integer, intent(out) :: rc
        rc = 0
        ! Trigger C++ registry hook explicitly to avoid dead-stripping
        call catchem_register_gaschem_cpp()
    end subroutine
end module
```

---

## C++ Bridge for Legacy Schemes

For physical schemes that have not yet been migrated to native C++ Kokkos, the generic `FortranProcess` bridging adapter wraps Fortran subroutines so they execute seamlessly within the centralized C++ schedule:

```cpp
namespace catchem {

    class FortranProcess : public ProcessInterface {
    private:
        std::string name;
        void (*fortran_run_callback)(void*); // Pointer to raw Fortran subroutine
    public:
        FortranProcess(const std::string& n, void (*cb)(void*))
            : name(n), fortran_run_callback(cb) {}

        std::string get_name() const override { return name; }

        void init(std::shared_ptr<StateManager> state) override {}

        void run(std::shared_ptr<StateManager> state) override {
            // 1. Sync GPU views to host CPU space
            state->sync_to_host();

            // 2. Pass StateManager handle and invoke legacy Fortran calculation
            fortran_run_callback(state.get());

            // 3. Flush CPU changes back to GPU device space
            state->sync_to_device();
        }

        void finalize() override {}
    };

} // namespace catchem
```

---

## Configuration & Diagnostics Integration

### YAML Configuration Loading
Processes load configuration settings in their `init` method by querying the thread-safe `StateManager` and associated YAML files:

```cpp
void GasChemProcess::init(std::shared_ptr<StateManager> state) {
    // Obtain configuration path
    std::string config_dir = state->config_dir;

    // Parse using yaml-cpp
    YAML::Node config = YAML::LoadFile(config_dir + "/micm_config.yaml");
    double absolute_tolerance = config["absolute_tolerance"].as<double>(1e-12);
}
```

### Diagnostics Binding
Processes write diagnostics directly into pre-allocated memory buffers in the `DiagnosticManager`:

```cpp
void PhotolysisProcess::run(std::shared_ptr<StateManager> state) {
    auto diag_mgr = state->get_diagnostic_manager();

    // Retrieve direct pointer address of registered diagnostic field
    double* jrate_ptr = diag_mgr->get_field_pointer("photolysis_rate_jfoo");

    // Write directly into diagnostic buffer
    jrate_ptr[cell_idx] = calculated_jrate;
}
```

---

## Best Practices

### Performance
1. **Zero-Copy Slicing**: Use `Kokkos::subview` to slice 3D arrays to 1D column vectors.
2. **Synchronize Efficiently**: Minimize the frequency of calling `sync_to_host()` and `sync_to_device()`. Keep computations inside Kokkos parallel device loops wherever possible.
3. **Register Fields on Startup**: Avoid dynamic metadata lookup or map queries in the `run()` loop; cache offsets and raw pointer locations inside `init()`.

### Code Quality
1. **Linker Safety**: Always declare and invoke the `extern "C"` registration hook for all C++ process classes.
2. **Exceptional Protection**: Ensure no native C++ exception escapes the dynamic registration hooks or BIND(C) boundaries; wrap code in complete `try-catch` blocks returning integer failure codes.

---

## See Also

- [State Management API](state-management.md) - Host/Device Kokkos View synchronization
- [Column Interface API](column-interface.md) - C++ GridManager and column views
- [Configuration API](configuration.md) - C++ YAML Configuration Manager
- [GasChem Process Documentation](../processes/gaschem/index.md) - Details on C++ MICM solver integration
- [Photolysis Process Documentation](../processes/photolysis/index.md) - Details on C++ TUV-x engine integration
# Process data contracts

Processes request canonical scientific fields or declared aliases. Missing lookup never substitutes an unrelated allocation. In particular, relative humidity cannot resolve to air density, and sulfate chemistry obtains interface geometric height from `Z` (or explicitly reconstructs it from `BXHEIGHT`) rather than pressure edges.

Fortran science bridges mark unified chemistry as host-current after modification. Device writers mark their views device-current. This writer declaration controls the next synchronization direction and prevents stale-side overwrites.
