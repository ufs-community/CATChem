# Column Interface API

This section covers the modernized column virtualization and grid layout APIs that enable high-performance 1D atmospheric column processing in CATChem.

## Overview

Under the C++20 and Kokkos architecture, the Column Interface system has been fully redesigned to achieve zero-overhead, copy-free execution on modern hardware:

- **catchem::GridManager**: Manages grid layout dimensions ($N_{\text{cols}}$, $N_{\text{levels}}$).
- **Kokkos Subviews**: Slices 1D columns from 3D multi-dimensional arrays instantly with zero allocation or data copying.
- **Kokkos Parallel For**: Executes column calculations in parallel across heterogeneous computing nodes (CPUs and GPUs).

---

## Core Concepts

### Column Virtualization via Subviews

Rather than copying full 3D fields to local 1D data structures, CATChem slices 1D column vectors using standard **Kokkos subviews**. Slices point directly to the underlying contiguous multidimensional memory:

```cpp
#include <Kokkos_Core.hpp>

void MyProcess::run(std::shared_ptr<StateManager> state) {
    auto n_cols = state->n_cols;

    // Parallelize column executions
    Kokkos::parallel_for("ProcessColumns", Kokkos::RangePolicy<Kokkos::HostSpace>(0, n_cols),
        [=](const int icol) {
            // Slice the 3D temperature view for column 'icol' with zero-copy
            // Coordinates: (column_index, level_index, species_or_attribute_index)
            auto col_temp = Kokkos::subview(state->met.temp, icol, Kokkos::ALL(), 0);

            // Access and modify elements directly in-place
            for (int k = 0; k < state->n_levels; ++k) {
                double t_k = col_temp(k);
                col_temp(k) = calculate_new_temp(t_k);
            }
        });
}
```

**Benefits**:
*   **Zero Duplication**: No heap allocations or data copies are made during slicing.
*   **Locality Optimization**: Leverages row-major or column-major Kokkos views matching underlying system architectures (e.g. `LayoutLeft` for CUDA/HIP GPUs, `LayoutRight` for OpenMP CPUs).

---

## Grid Layout & GridManager

### GridManager

The C++ class managing physical grid dimensions:

```cpp
#pragma once

namespace catchem {

    class GridManager {
    public:
        const int n_cols;                       ///< Total number of contiguous horizontal columns.
        const int n_levels;                     ///< Total number of vertical levels.

        GridManager(int nc, int nl) : n_cols(nc), n_levels(nl) {}
    };

} // namespace catchem
```

---

## Processing Patterns

### 1. Parallel Column Processing (CPU/GPU)
Kokkos automatically maps parallel column loops to available hardware threads (e.g., OpenMP threads on multi-core CPUs, or threads/blocks on GPUs):

```cpp
Kokkos::parallel_for("ColumnAerosols", Kokkos::RangePolicy<Kokkos::HostSpace>(0, state->n_cols),
    [=](const int icol) {
        // Run aerosol settling scheme on independent 1D column
        run_aerosol_column_kernel(state, icol);
    });
```

### 2. Contiguous 3D Processing (MDRangePolicy)
For operations that are completely independent of vertical column contexts, use Kokkos multi-dimensional range policies to loop over all grid cells concurrently:

```cpp
Kokkos::parallel_for("ScaleMetPres", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {state->n_cols, state->n_levels}),
    KOKKOS_LAMBDA(const int icol, const int k) {
        state->met.pressure(icol, k, 0) *= 1.01;
    });
```

---

## Best Practices

### Performance
1.  **Prefer Kokkos Subviews**: Never use manual copying or temporary `std::vector` allocations to represent columns. Subviews are instant and run on both CPU and GPU memory spaces in-place.
2.  **Align Memory Layouts**: Match Kokkos View layouts with execution backends. Use `LayoutLeft` for GPU acceleration to achieve coalesced global memory access, or `LayoutRight` for sequential CPU cache locality.

### Safety
1.  **Prevent Out-of-Bounds**: Always use grid boundaries from `state->n_cols` and `state->n_levels` to guard level-iteration loops.
2.  **Coordinate Thread Space**: Avoid calling host-only or BIND(C) routines inside parallel Kokkos device lambda loops. Keep parallel lambda bodies purely mathematical and scientific.

---

## See Also

- [State Management API](state-management.md) - Multidimensional Kokkos Views
- [Process Interface API](process-interface.md) - Scheduling and Process Registry
- [GasChem Process Documentation](../processes/gaschem/index.md) - 3D Grid flattening to 1D MICM arrays

---
