# Implementation Plan: Zero-Copy ESMF/NUOPC Cap and Memory Interoperability

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fully streamline coupling transformations inside `catchem_nuopc_interface.F90` to bind raw pointers directly to modern C++ unmanaged Views and delete legacy double-buffering.

## Global Constraints
* Internally representation inside C++ uses `Kokkos::View` LayoutLeft host views and C++20 `mdspan` backports mapping contiguous boundary memory.
* High-Performance: 0 copies and 0 allocations in Fortran inside the coupling transformation routines.

---

## Task List & Execution Schedule

### Task 1: Refactor 2D Surface Pointer Bindings
*   **Goal:** Map raw double precision pointers directly via model%bind_met_2d.
*   **Target File:** `drivers/nuopc/catchem_nuopc_interface.F90`
*   **Steps:**
    1. Locate Case(2) inside `transform_field_to_catchem` in `catchem_nuopc_interface.F90`.
    2. Replace local copies for all standard double precision arrays (PS, TS, PBLH, USTAR, HFLUX, OBK, LAT, LON) with direct `call cc_wrap%catchem_model%bind_met_2d(...)` operations using `c_loc(fptr2d(1,1))`.

### Task 2: Refactor 3D Meteorological Pointer Bindings
*   **Goal:** Deallocate legacy buffers and bind volumetric 3D arrays to Kokkos views directly.
*   **Target File:** `drivers/nuopc/catchem_nuopc_interface.F90`
*   **Steps:**
    1. Locate Case(3) inside `transform_field_to_catchem`.
    2. Strip out local allocations for `fptr3d_rev(ni, nj, nk1)` and its deep nested element-copy loops.
    3. Replace with a direct pointer binding: `call cc_wrap%catchem_model%bind_met_3d(trim(field_map%catchem_var), c_loc(fptr3d(1,1,1)))`.

### Task 3: Refactor 4D Chemistry Pointer Bindings
*   **Goal:** Map 4D tracer concentrations in-place.
*   **Target File:** `drivers/nuopc/catchem_nuopc_interface.F90`
*   **Steps:**
    1. Locate Case(4) inside `transform_field_to_catchem`.
    2. Strip out `fptr4d_rev` allocations, zero-out loops, and nested element-copy loops.
    3. Replace with a single, direct, zero-copy binder: `call cc_wrap%catchem_model%bind_unified_chemistry(c_loc(fptr4d(1,1,1,1)))`.

### Task 4: Compile & Run Verification Suite
*   **Goal:** Compile the full codebase in parallel and run tests.
*   **Steps:**
    1. Run compilation: `docker exec dazzling_lehmann cmake --build /workspace/build-cece --parallel`.
    2. Run tests: `docker exec dazzling_lehmann ctest --test-dir /workspace/build-cece --output-on-failure`.
