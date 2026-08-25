# Documentation Update Plan: Modernized C++ Core & Native Processes

This plan outlines the specific updates and additions required to ensure that CATChem's user-facing documentation completely and accurately reflects the modernized C++ core, the native C++ APIs, and the native C++ processes (GasChem and Photolysis).

---

## 1. Objectives

1.  **Eliminate Outdated Information**: Transition core interface, API, state management, and configuration guides from 100% legacy Fortran definitions to modern C++20 and Kokkos structures.
2.  **Document Modernized Architecture**: Fully document the centralized C++ Core orchestrator (`catchem::Core`), the thread-safe `StateManager` (using Kokkos Views), and the `DiagnosticManager`.
3.  **Add Native Process Reference**: Create comprehensive user-facing documentation for the newly-integrated `GasChemProcess` (wrapping Model Independent Chemistry Module - MICM via `musica`) and `PhotolysisProcess` (wrapping TUV-x photolysis solver).
4.  **Preserve Fortran Bridges**: Document how legacy Fortran processes run seamlessly under the centralized C++ execution loop using the `catchem::FortranProcess` wrapper and generic BIND(C) flat APIs.

---

## 2. Update Matrix & Scopes

### A. General Guides
*   **`docs/index.md`**: Update Key Features, Quick Start (C++ and Fortran building), and Integration examples to show both C++ `catchem::Core` and Fortran delegate wrappers.
*   **`docs/developer-guide/architecture.md`**: Update High-Level System Components diagram to include C++ Core, Kokkos device spaces, and BIND(C) boundaries. Rewrite code blocks to C++ namespaces.

### B. Core APIs (`docs/api/`)
*   **`docs/api/index.md`**: Incorporate C++ types (`catchem::Core`, `catchem::StateManager`, `catchem::ProcessInterface`, etc.) and registration/creation common patterns.
*   **`docs/api/process-interface.md`**: Document C++ `ProcessInterface`, dynamic `ProcessRegistry` creator lambdas, linker-safe dynamic registration callbacks, and the `FortranProcess` bridge.
*   **`docs/api/state-management.md`**: Detail `catchem::StateManager`, `InteropField` pointer binding, unified chemical concentration views, and Kokkos host-to-device synchronization.
*   **`docs/api/configuration.md`**: Document the YAML configuration loading using `yaml-cpp` inside the C++ `ConfigManager`.
*   **`docs/api/column-interface.md`**: Explain the `GridManager` and column virtualization via zero-copy Kokkos subviews.

### C. Processes (`docs/processes/`)
*   **`docs/processes/index.md`**: Add `GasChem` and `Photolysis` to the available processes list.
*   **`docs/processes/MODERNIZED_PROCESSES.md`**: Add GasChem C++ MICM solver and Photolysis C++ TUV-x solver as modernized processes, detailing their high-performance features.
*   **`docs/processes/gaschem/index.md` (NEW)**: Write detailed documentation for the `GasChemProcess` (MICM, unit conversions, automatic photolysis coupling, boundary clamping).
*   **`docs/processes/photolysis/index.md` (NEW)**: Write detailed documentation for the `PhotolysisProcess` (TUV-x engine, meteorological mapping, layer-midpoint interpolation).

---

## 3. Verification Strategy

We will verify our documentation changes using the following gate function:
1.  Verify that all modified and newly-created documentation files compile cleanly with Markdown/Doxygen tools.
2.  Assert that no outdated Fortran-only references remain for components that are now fully owned by C++ (or clearly indicate them as legacy wrapper interfaces).
3.  Confirm index correctness by executing the `doc-tools.sh check-freshness` command.
