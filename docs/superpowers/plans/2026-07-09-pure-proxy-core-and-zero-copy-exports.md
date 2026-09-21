# Implementation Plan: Pure Proxy Core and Zero-Copy Coupled Exports

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the Pure Proxy Pattern for `metstate_mod.F90` and `chemstate_mod.F90` to eliminate duplicate Fortran allocations (plugging a massive memory leak), and refactor ESMF exports to execute completely zero-copy.

---

## Task List & Execution Schedule

### Task 1: Convert `MetStateType` to Pure Proxy (Plugs Memory Leak)
*   **Goal:** Deallocate internal Fortran allocations during metstate initialization, forcing pointers to stay null until bound to C++ views.
*   **Target File:** `src/core/metstate_mod.F90`
*   **Steps:**
    1. Strip out the `allocate_metstate_arrays` helper and any case statements inside `allocate_arrays`.
    2. Update `metstate_init` to purely assign geometry dimensions and zero out status variables, avoiding heap array allocations.

### Task 2: Convert `ChemStateType` to Pure Proxy (Plugs Memory Leak)
*   **Goal:** Eradicate duplicate allocations of `unified_conc` and `ChemSpecies%conc` pointers in Fortran.
*   **Target File:** `src/core/chemstate_mod.F90`
*   **Steps:**
    1. Prune local `allocate` statements for `unified_conc` and individual concentrations inside `chemstate_init`.
    2. Retain allocation only for the lightweight descriptors `ChemSpecies(max_species)` structures.
    3. Register pointer sync inside state-retrieval to bind `%conc` slices directly to C++ unified chemistry addresses dynamically.

### Task 3: Streamline Zero-Copy ESMF Exports in Cap
*   **Goal:** Avoid duplicate buffers and loop copying for diagnostic and tracer exports.
*   **Target File:** `drivers/nuopc/catchem_nuopc_interface.F90`
*   **Steps:**
    1. In `transform_catchem_to_field` Case(4), bypass element-by-element copy loops for tracers that are already bound in-place (`fptr4d` already maps `chem.conc`).
    2. Access diagnostics (PM2.5 / PM10) via standard direct references without duplicating arrays.

### Task 4: Compile & Run Verification Suite
*   **Goal:** Confirm successful compilation and 100% green unit tests.
*   **Steps:**
    1. Run build: `docker exec dazzling_lehmann cmake --build /workspace/build-cece --parallel`.
    2. Run tests: `docker exec dazzling_lehmann ctest --test-dir /workspace/build-cece --output-on-failure`.
