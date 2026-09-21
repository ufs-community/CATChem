# Core Modernization Design: ConfigManager and GridManager

## Overview
This design outlines the plan for porting the CATChem `ConfigManager` and `GridManager` from Fortran to C++20, resolving the final architectural questions around column virtualization (`VirtualColumn`). This is a crucial step in the broader modernization effort to make the CATChem core a high-performance C++ orchestrator.

## 1. `catchem::ConfigManager`
- **Role:** Replaces the legacy Fortran `ConfigManager_Mod.F90`. Responsible for parsing YAML configuration files directly into structured C++ data.
- **Components:**
  - `src/core/catchem_config_manager.hpp`
  - `src/core/catchem_config_manager.cpp`
  - `catchem::ConfigData` struct to hold parsed data (e.g., `RuntimeConfig`, `FilePathConfig`).
- **Implementation Details:**
  - Will utilize `yaml-cpp` to load `CATChem_new_config.yml`.
  - It extracts core dimensions (`nx`, `ny`, `nz`), simulation runtime parameters (`dt`, `nsteps`), and active process flags.
  - `catchem::Core` will own a `std::shared_ptr<ConfigManager>`. Upon initialization, the core will load the configuration and use the resulting dimensions to size the `GridManager` and `StateManager`.

## 2. `catchem::GridManager` & `catchem::GridGeometry`
- **Role:** Replaces `GridManager_Mod.F90` and `GridGeometry_Mod.F90`. Manages the 3D execution space, dimensions, and geographic coordinate mapping.
- **Components:**
  - `src/core/catchem_grid_manager.hpp`
  - `src/core/catchem_grid_manager.cpp`
- **Implementation Details:**
  - `catchem::GridGeometry` stores the integer grid dimensions (`nx`, `ny`, `nz`).
  - To support coordinate-aware calculations (e.g., solar zenith angle), it will hold `InteropField` bindings for geographic arrays such as `lat`, `lon`, `grid_area`, and `dz`.
  - `catchem::Core` will instantiate `GridManager` based on the dimensions extracted by the `ConfigManager`.

## 3. Physics Execution Architecture & `VirtualColumn`
- **Role:** Formalizes how C++ physical schemes access memory and execute in parallel, explicitly resolving the future of the `VirtualColumn` abstraction.
- **Architectural Decision:**
  - We will **bypass** the creation of a C++ `VirtualColumn` struct.
  - C++ physics schemes will use the idiomatic Kokkos approach: receiving the full 3D `Kokkos::View`s (or the full `StateManager`) and launching `Kokkos::parallel_for` or `Kokkos::MDRangePolicy` kernels directly over `n_cols` (and optionally `n_levels`).
  - This provides the maximum performance on GPUs by minimizing register pressure and overhead that would otherwise be incurred by passing 1D `Kokkos::subview` objects.
  - Legacy Fortran schemes will continue to use the Fortran `VirtualColumn_Mod.F90` through the dynamic C++-to-Fortran bridge.

## 4. Fortran Interoperability and Legacy Support
- **Role:** Ensure that unported legacy Fortran processes continue to function seamlessly alongside the new C++ core components.
- **Implementation Details:**
  - **Configuration & Grid Access:** We will expose C-API endpoints (e.g., `catchem_get_grid_dimensions`, `catchem_get_config_value`) so the legacy Fortran code can query the C++ `ConfigManager` and `GridManager` directly, replacing redundant Fortran parsers.
  - **VirtualColumn Preservation:** Due to the strided memory layout of a 1D column in a `Kokkos::LayoutLeft` 3D array (`nx, ny, nz`), creating a `Kokkos::subview` in C++ yields a non-contiguous slice. Since Fortran's `c_f_pointer` only accepts contiguous blocks, we cannot easily pass C++ subviews across the language boundary.
  - **Execution Flow:** Therefore, when `catchem::Core` executes a legacy Fortran process, it will pass the **raw pointer to the entire 3D array** to the Fortran bridge callback. The Fortran bridging layer will use `c_f_pointer` to construct a 3D Fortran array pointer, and then utilize the legacy `VirtualColumn_Mod.F90` to create the strided 1D slices (`col_ptr => conc(i, j, :)`) that the legacy Fortran schemes expect. This keeps all strided pointer math safely inside the language that owns the object.
