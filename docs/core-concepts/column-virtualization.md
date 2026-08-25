# Column Virtualization

This section covers the column virtualization concepts that enable efficient 1D atmospheric processing in CATChem.

## Overview

Column virtualization is a key performance optimization in CATChem. Instead of processing the entire 3D atmospheric grid at once, CATChem treats the grid as a collection of independent 1D vertical columns. This approach has several advantages:

- **Performance**: By processing data in 1D columns, CATChem can take advantage of modern CPU architectures, which are highly optimized for linear data access patterns. This results in better cache utilization and fewer cache misses.
- **Scalability**: Column-based processing is highly scalable. Since each column can be processed independently, the workload can be easily distributed across multiple processors or nodes.
- **Simplicity**: By abstracting the grid into a collection of 1D columns, column virtualization simplifies the development of atmospheric processes. Developers can focus on the physics and chemistry of a single column, without having to worry about the complexity of the 3D grid.

## Core Concepts

### Zero-Copy Subview Column Slicing

Under CATChem's C++ and Kokkos architecture, column virtualization is achieved via zero-copy **Kokkos subviews** (`Kokkos::subview`). Slicing a 3D memory view down to a 1D column vector creates an unmanaged view that points directly to contiguous backing memory without any heap allocation or array copying:

```cpp
// Slice 1D column from 3D temperature view with zero copies
auto col_temp = Kokkos::subview(state->met.T->view(), icol, Kokkos::ALL(), 0);
```

### Fortran ScienceBridge Column Slicing

For physical processes written in Fortran, the `ScienceBridge` module receives raw C pointers (`c_ptr`) from C++, converts them to Fortran array pointers (`c_f_pointer`), and iterates over columns using standard Fortran array section slicing:

```fortran
do icol = 1, n_cols
   call compute_scheme(n_levels, n_species, dt, f_conc(icol, :, :), f_tendency(icol, :, :))
end do
```

## Data Access Patterns

The column virtualization system supports high-performance parallel execution:

- **Kokkos Parallel For**: Executes column kernels concurrently using `Kokkos::parallel_for` across OpenMP, CUDA, or HIP execution spaces.
- **Fortran OpenMP Parallelization**: Process ScienceBridges can process independent columns concurrently across multi-core CPUs.

## Process Integration

Atmospheric processes integrate with the column system by implementing C++ `catchem::ProcessInterface`. Processes retrieve meteorological and chemical view pointers from `StateManager`, perform column calculations (natively in C++ or via a Fortran `ScienceBridge`), and update device views in-place.

## Performance Considerations

Column virtualization is a key performance optimization in CATChem. By processing data in 1D columns, CATChem can take advantage of modern CPU architectures, which are highly optimized for linear data access patterns. This results in better cache utilization and fewer cache misses.

In addition, column virtualization enables natural parallelization. Since each column can be processed independently, the workload can be easily distributed across multiple processors or nodes.
