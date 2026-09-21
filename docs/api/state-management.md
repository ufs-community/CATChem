# State Management API

This section covers the modernized C++ state management APIs in CATChem, providing unified, high-performance multidimensional data handling on both CPUs and GPUs.

## Overview

The modernized state management system is centered around C++20 and **Kokkos**. It provides:

- **catchem::StateManager**: Central data repository holding Kokkos Views for chemical species and meteorological attributes.
- **catchem::InteropField**: Dynamic field wrapper supporting raw pointer bindings with zero-copy mapping.
- **Dual-Space Coherence**: Memory tracking and explicit `sync_to_host()` and `sync_to_device()` synchronization boundaries.
- **Zero-Copy Fortran Interoperability**: Direct dynamic binding of C++ allocated buffers to legacy Fortran pointer slices, eliminating double-buffering.

---

## Core Components

### catchem::StateManager

The C++ StateManager coordinates physical grid dimensions and acts as the single source of truth for all simulation states:

```cpp
#pragma once
#include <Kokkos_Core.hpp>
#include <memory>
#include <string>
#include <unordered_map>

namespace catchem {

    class StateManager {
    public:
        const int n_cols;                       ///< Number of horizontal columns.
        const int n_levels;                     ///< Number of vertical levels.
        const int n_species;                    ///< Number of chemical species.

        // Unified 3D chemical concentrations view (Layout: Species, Columns, Levels)
        Kokkos::View<double***, Kokkos::LayoutLeft> chem_conc_device;
        Kokkos::View<double***, Kokkos::LayoutLeft>::HostMirror chem_conc_host;

        // Dynamic Meteorological and Diagnostic state registries (Key-Value map to InteropField)
        std::unordered_map<std::string, InteropField> met_fields;
        std::unordered_map<std::string, InteropField> diag_fields;

        StateManager(int nc, int nl, int ns);

        // Explicitly copies modified host data to the GPU device space
        void sync_to_device();

        // Explicitly copies modified device calculations back to CPU memory buffers
        void sync_to_host();
    };

} // namespace catchem
```

---

## Data Access Patterns

### 1. C++ Native Access (Kokkos parallel kernels)
C++ kernels execute directly on GPU device space or CPU host space using standard Kokkos accessors:

```cpp
// Kernel executing on GPU Device space
auto conc_device = state->chem_conc_device;
Kokkos::parallel_for("ScaleChemistry", Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0}, {n_species, n_cols, n_levels}),
    KOKKOS_LAMBDA(const int s, const int icol, const int k) {
        conc_device(s, icol, k) *= 1.05; // Apply scaling factor in-place
    });
```

### 2. Zero-Copy Fortran Access (C-API boundaries)
Legacy Fortran codes query the StateManager dynamically using BIND(C) flat APIs and construct standard pointer slices without any allocations or copies:

```fortran
use iso_c_binding
use StateManager_Mod
implicit none

type(StateContainerType) :: container
real(dp), pointer :: temp_ptr(:,:,:) => null()
type(c_ptr) :: raw_c_addr
integer :: nx, ny, nz, rc

! 1. Query dimensions
call container%get_dimensions(nx, ny, nz)

! 2. Retrieve C-pointer address of the meteorological Temperature field
raw_c_addr = catchem_state_get_pointer_3d(container%cpp_handle, "temp")

! 3. Standardize raw pointer as a native 3D Fortran slice pointer
call c_f_pointer(raw_c_addr, temp_ptr, [nx, ny, nz])

! 4. Modify memory in-place directly on C++ allocated buffers
temp_ptr(i, j, k) = 298.15
```

---

## Dual-Space Synchronization

To prevent race conditions across heterogeneous systems (CPU host and GPU devices), developers must explicitly schedule memory synchronizations at computational boundaries:

```cpp
// 1. Host-model updates variables on the CPU host (via Fortran or C-API)
catchem_state_bind_3d(state_ptr, "temp", host_temp_buffer);

// 2. Synchronize modified host memory to GPU device before C++ Kokkos execution
state->sync_to_device();

// 3. Execute Kokkos GPU kernel
core->run_timestep(300.0);

// 4. Synchronize GPU results back to CPU host buffers for I/O and diagnostics
state->sync_to_host();
```

---

## Thread Safety

*   **Kokkos Space**: Kokkos Views utilize standard memory traits and layouts ensuring safe concurrent execution inside multi-threaded OpenMP, CUDA, or HIP dispatch zones.
*   **Fortran OMP Boundaries**: Fortran drivers can safely process columns concurrently using standard OpenMP directives because all C++ pointer-bindings are thread-safe and contiguous.

---

## Best Practices

### Performance
1.  **Enforce Zero-Copy**: Never allocate intermediate local arrays inside wrappers or process modules. Retrieve the raw pointer via the dynamic interop utility and write directly to C++ memory buffers.
2.  **Combine Synchronizations**: Group multiple field updates together and call `sync_to_device()` exactly once before running a timestep.
3.  **Use Contiguous Layouts**: Leverage `Kokkos::LayoutLeft` or standard C-row-major contiguous buffers to maximize hardware pre-fetcher cache alignment.

### Safety
1.  **Check Bind Status**: Ensure that physical fields (like `temp`, `pres`, `AIRDEN_DRY`) are successfully bound before invoking calculations; non-bound fields will return a null pointer, raising a segmentation fault.
2.  **Explicit Memory Sync**: Always call `sync_to_host()` before invoking NetCDF output writers or legacy diagnostics.

---

## See Also

- [Process Interface API](process-interface.md) - Physics processes lifecycles
- [Column Interface API](column-interface.md) - GridManager column slicing
- [Configuration API](configuration.md) - YAML configuration loadings
- [GasChem Process Documentation](../processes/gaschem/index.md) - C++ MICM state handling

---

**Auto-Generated Documentation:** [Complete State Management Reference](../CATChem/namespacestatemanager__mod.md)
# Dataflow ownership and freshness

Shared meteorology and unified chemistry use non-owning `InteropField` views. Each view has immutable extents, a generation, availability, and a latest-writer state. A NUOPC import generation invalidates the preceding generation before current fields are rebound. Missing optional fields therefore cannot silently retain a current pointer.

Concentration slabs are mechanism-relative: callers resolve a species name against the active species configuration, validate the returned 1-based index, and only then form the checked `[column, level]` slab. No process or emission mapping may assume a fixed chemical-species slot.

`HostCurrent` permits only host-to-execution synchronization; `DeviceCurrent` permits only execution-to-host synchronization; `Synchronized` performs no copy. Rebinding marks host storage current.
