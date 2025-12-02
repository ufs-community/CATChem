# Column Interface API

This section covers the column virtualization APIs that enable efficient 1D atmospheric processing in CATChem.

## Overview

The Column Interface system provides:

- **ColumnType**: 1D atmospheric column data structure
- **ColumnInterface_Mod**: Column access and manipulation
- **Grid virtualization**: Efficient processing across grid columns
- **Memory optimization**: Cache-friendly data access patterns

## Core Concepts

### Column Virtualization

CATChem processes atmospheric data as independent 1D columns rather than full 3D grids:

```fortran
use ColumnInterface_Mod
type(ColumnType) :: column

! Get column from 3D grid
call container%get_column(i, j, column, rc)

! Process column data efficiently
do k = 1, column%nz
    ! Atmospheric level k processing
    call process_atmospheric_level(column, k)
end do

! Update grid with processed column
call container%update_column(i, j, column, rc)
```

**Benefits:**

- **Scalability**: Linear scaling with grid resolution
- **Cache efficiency**: Better data locality

**Auto-Generated Documentation:** [Column Interface Reference](../CATChem/namespacecolumninterface__mod.md)

## Column Data Structure

### ColumnType

The basic column data structure:

```fortran
type :: ColumnType
    integer :: nz                           ! Number of vertical levels
    real(fp), allocatable :: pressure(:)    ! Pressure levels [Pa]
    real(fp), allocatable :: temperature(:) ! Temperature profile [K]
    real(fp), allocatable :: height(:)      ! Height levels [m]
    real(fp), allocatable :: chem_data(:,:) ! Chemical concentrations
    ! ... additional fields
end type
```

### Column Access

```fortran
! Direct column access
call container%get_column(i, j, column, rc)

! Safe column access (thread-safe)
call container%get_column_safe(i, j, column, rc)

! Batch column access
call container%get_columns_batch(i_start, i_end, j_start, j_end, &
                                columns, rc)
```

## Processing Patterns

### Sequential Processing

```fortran
! Process all columns sequentially
do j = 1, ny
    do i = 1, nx
        call container%get_column(i, j, column, rc)
        call my_process_column(column, rc)
        call container%update_column(i, j, column, rc)
    end do
end do
```

### Parallel Processing

```fortran
! OpenMP parallel column processing
!$OMP PARALLEL DO PRIVATE(column, rc)
do j = 1, ny
    do i = 1, nx
        call container%get_column_safe(i, j, column, rc)
        call my_process_column(column, rc)
        call container%update_column_safe(i, j, column, rc)
    end do
end do
!$OMP END PARALLEL DO
```

### Process Integration

```fortran
! Column-capable process
type, extends(ColumnProcessInterface_t) :: MyProcessType
contains
    procedure :: run_column => my_column_process
end type

subroutine my_column_process(this, column, rc)
    class(MyProcessType), intent(inout) :: this
    type(ColumnType), intent(inout) :: column
    integer, intent(out) :: rc

    ! Column-based atmospheric calculations
    call this%calculate_column_physics(column, rc)
end subroutine
```

## Memory Management

### Efficient Memory Usage

```fortran
! Column data uses minimal memory allocation
type(ColumnType) :: column

! Automatic memory management
call container%get_column(i, j, column, rc)  ! Allocates column data
call container%update_column(i, j, column, rc)  ! Automatically cleans up
```

### Memory Optimization

```fortran
! Pre-allocate column workspace for performance
call container%allocate_column_workspace(nz, n_species, rc)

! Use workspace for repeated operations
do j = 1, ny
    do i = 1, nx
        call container%get_column_workspace(i, j, column, rc)
        call process_column(column, rc)
        call container%return_column_workspace(i, j, column, rc)
    end do
end do
```

## Performance Considerations

### Cache Optimization

Column processing optimizes CPU cache usage:

- **Spatial locality**: Processing contiguous vertical levels
- **Temporal locality**: Reusing column data for multiple operations
- **Reduced cache misses**: 1D access patterns vs. 3D strided access

### Vectorization

Modern compilers can vectorize column operations:

```fortran
! Vectorizable column operation
do k = 1, column%nz
    column%temperature(k) = column%temperature(k) + heating_rate(k) * dt
end do
```

### Parallel Scaling

Column processing enables natural parallelization:

- **Grid-level parallelism**: Process multiple columns simultaneously
- **Process-level parallelism**: Multiple processes per column
- **Thread safety**: Independent column operations

## Advanced Features

### Column Interpolation

```fortran
! Interpolate between pressure levels
call column%interpolate_to_pressure(target_pressure, interpolated_data, rc)

! Interpolate between height levels
call column%interpolate_to_height(target_height, interpolated_data, rc)
```

### Column Diagnostics

```fortran
! Calculate column-integrated quantities
call column%integrate_column('CO', total_co, rc)
call column%calculate_column_burden('aerosol', burden, rc)
```

**Auto-Generated Documentation:** [Diagnostic Manager Reference](../CATChem/namespacediagnosticmanager__mod.md)

**Auto-Generated Documentation:** [Diagnostic Interface Reference](../CATChem/namespacediagnosticinterface__mod.md)

### Boundary Conditions

```fortran
! Apply surface boundary conditions
call column%set_surface_value(species_idx, surface_concentration, rc)

! Apply top boundary conditions
call column%set_top_value(species_idx, top_concentration, rc)
```

## Error Handling

```fortran
use Error_Mod

! Column operation with error handling
call container%get_column(i, j, column, rc)
if (rc /= CC_SUCCESS) then
    call error_mgr%report_error(ERROR_COLUMN_ACCESS, &
                               'Failed to get column data', rc, &
                               additional_info='Grid coordinates: ' // &
                               trim(str(i)) // ', ' // trim(str(j)))
    return
endif
```

**Auto-Generated Documentation:** [Error Handling Reference](../CATChem/namespaceerror__mod.md)

## Best Practices

### Performance

1. **Use column processing** for all atmospheric calculations
2. **Minimize column copies** - work with pointers when possible
3. **Batch operations** when processing multiple columns
4. **Preallocate workspaces** for repeated operations

### Code Quality

1. **Check return codes** for all column operations
2. **Handle edge cases** (surface, top boundary)
3. **Validate column data** before processing
4. **Use appropriate precision** for calculations

### Threading

1. **Use thread-safe methods** in parallel regions
2. **Avoid shared column data** between threads
3. **Use private column variables** in OpenMP regions
4. **Consider NUMA effects** for large grids

## See Also

- [State Management API](state-management.md) - Data container interfaces
- [Process Interface API](process-interface.md) - Column-capable process development
- [Performance Guide](../user-guide/advanced_topics/performance.md) - Optimization strategies
- [Column Virtualization Guide](../user-guide/advanced_topics/column-virtualization.md) - Architecture details

---

**Auto-Generated Documentation:** [Complete Column Interface Reference](../CATChem/namespacecolumninterface__mod.md)
