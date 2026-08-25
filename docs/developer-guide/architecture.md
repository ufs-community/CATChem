# Developer Architecture Guide

This guide provides a comprehensive overview of CATChem's modernized C++ software architecture, design principles, and implementation patterns for developers working on the codebase.

---

## High-Level Architecture Overview

CATChem has transitioned to a high-performance **C++20** and **Kokkos** central orchestration framework. This architecture decouples scientific calculation kernels from legacy memory systems, enabling execution on both CPU multi-core and GPU (CUDA/HIP) high-performance computing (HPC) nodes.

### System Components

```mermaid
graph TB
    A["Host Model / Driver Layer"] -->|BIND(C) flat APIs| B["C++ Delegate Boundary (catchem_api.hpp)"]
    B --> C["C++ Core Orchestrator (catchem::Core)"]

    C --> D["Config Manager (catchem::ConfigManager)"]
    C --> E["State Manager (catchem::StateManager)"]
    C --> F["Grid Manager (catchem::GridManager)"]
    C --> G["Diagnostic Manager (catchem::DiagnosticManager)"]

    D --> H["yaml-cpp Parser"]

    E --> I["Unified Chem Concentrations (Kokkos Views)"]
    E --> J["Dynamic Met Registers (InteropField)"]

    C --> K["C++ Process Registry (catchem::ProcessRegistry)"]
    K --> L["Native C++ Processes"]
    K --> M["Fortran Bridge Processes (catchem::FortranProcess)"]

    L --> N["GasChem (MICM/musica)"]
    L --> O["Photolysis (TUV-x)"]
    L --> P["Settling (Kokkos)"]
```

---

## Layer Responsibilities

### 1. Host Model / Driver Boundary
*   Exposes a BIND(C) flat endpoint interface declared in `catchem_api.hpp`.
*   Passes raw host memory pointers of meteorological and concentration buffers from FV3 / ESMF / NUOPC driver grids.
*   Enforces zero-copy operations on startup by mapping physical pointers directly to pre-allocated C++ memory structures.

### 2. C++ Core Engine (`catchem::Core`)
*   Maintains the single source of truth for the entire physical simulation.
*   Instantiates Config, Grid, State, and Diagnostic managers.
*   Coordinates scheduled atmospheric physical/chemical transformation timelines.

### 3. State Management (`catchem::StateManager`)
*   Allocates and manages multidimensional Kokkos Views (Dual Space layout supporting CPU host and GPU device).
*   Enforces explicit device-to-host and host-to-device memory synchronization loops (`sync_to_host()` and `sync_to_device()`).
*   Uses `InteropField` objects to warp external raw host pointers with zero-copy unmanaged Views.

### 4. Process Layer (`catchem::ProcessInterface`)
*   Implements concrete atmospheric physical/chemical schemes.
*   Leverages Kokkos subviews to slice 1D columns instantly for localized transformation calculations (e.g. settling, chemistry).
*   Utilizes a linker-safe Dynamic Registration callback pattern to prevent optimization dead-code stripping.

---

## Core Design Principles

### 1. Zero-Copy Pointer Mapping
To avoid massive data copy loops and duplicate allocations across the C++/Fortran boundary, physical variables are mapped directly using standard C pointers and `c_f_pointer`:
*   Fortran registers outer model array addresses inside C++ registries on startup.
*   C++ wraps these contiguous host memory addresses as unmanaged host Kokkos Views.
*   GPU backends use device-mirrored spaces, synchronized explicitly using Kokkos memory copy tools.

### 2. Separation of Concerns
*   **Orchestration**: Managed in C++ (`catchem::Core` and managers).
*   **Calculation**: Scientific code resides inside process modules (`src/process/<name>/**`).
*   **Metadata**: Configuration and species chemical properties are read from unified YAML files and queried dynamically.

---

## Directory Organization & Code Structure

The source code is organized into three major functional directories:

*   **`src/core/`**: Central C++ framework files (Core, StateManager, InteropField, ConfigManager, DiagnosticManager) and the BIND(C) export layer (`catchem_api.cpp`).
*   **`src/process/`**: Standalone physical/chemical processes. Each process compiles to a static library (e.g., `catchem_process_gaschem`) and dynamically registers itself to the C++ core registry on startup.
*   **`src/external/`**: External dependencies (such as `musica`, which contains MICM and TUV-x engines).

---

## Timestep Execution Lifecycles

During a simulation timestep, the execution pipeline flows sequentially through the centralized C++ core loop:

```
[Driver / Host Model] ──(run_timestep)──> [catchem::Core]
                                               │
                                               ▼
                                  [1. sync_to_device()]
                                  (Flush Host updates to GPU)
                                               │
                                               ▼
                                  [2. Scheduled Processes]
                                  Loop over active processes:
                                    • Process1->run(state)
                                    • Process2->run(state)
                                               │
                                               ▼
                                  [3. sync_to_host()]
                                  (Pull Device outputs to Host)
                                               │
                                               ▼
                      [Driver / Host retrieves calculated tendencies]
```

---

## Best Practices for Developers

### Thread and GPU Portability
*   Always write process loops using standard **Kokkos Parallel policies** (e.g. `Kokkos::parallel_for`) to ensure performance compatibility across CUDA, HIP, OpenMP, and threads.
*   Avoid standard standard C++ heap allocations (like `new` or `malloc`) or STL vector operations inside the `run()` loops; perform all resource allocations in the `init()` phase.

### Mixed-Language Bridging
*   When integrating or wrapping unported legacy Fortran schemes, declare and execute them under the generic `catchem::FortranProcess` bridge.
*   Avoid exposing raw C++ standard exceptions across the ISO_C_BINDING boundary; wrap all exported flat APIs in complete `try-catch` blocks and return integer error statuses.

---

## See Also

- [State Management API](../api/state-management.md) - Kokkos View management and pointer mappings
- [Process Interface API](../api/process-interface.md) - C++ Processes, Registry, and linker safety
- [Column Interface API](../api/column-interface.md) - Grid layout and column subviews slicing
- [Configuration API](../api/configuration.md) - YAML Config Manager

---
