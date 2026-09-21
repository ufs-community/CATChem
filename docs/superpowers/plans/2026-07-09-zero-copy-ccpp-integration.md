# Implementation Plan: Zero-Copy CCPP Integration

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Simplify the CCPP Fortran interface (`ccpp_catchem_interface.F90` and `ccpp_catchem_interface.meta`) to completely remove standalone Fortran process states (Dust, Seasalt, Drydep) and bind standard CCPP physical variables and chemistry tracers directly to C++ unmanaged LayoutLeft Views.

---

## Task List & Execution Schedule

### Task 1: Streamline Interface Variables and Meta Table (`ccpp_catchem_interface.meta`)
*   **Goal:** Restructure standard names and lists inside CCPP's metadata description file.
*   **Target File:** `drivers/ccpp/ccpp_catchem_interface.meta`
*   **Steps:**
    1. Delete the CCPP arguments for legacy `DustState`, `SeaSaltState`, and `DryDepState`.
    2. Keep only standard grid parameters, meteorology arrays, and the unified chemical tracer concentration array.

### Task 2: Refactor CCPP Entry Points (`ccpp_catchem_interface.F90`)
*   **Goal:** Strip duplicate state storage and call C++ StateManager binders on-the-fly.
*   **Target File:** `drivers/ccpp/ccpp_catchem_interface.F90`
*   **Steps:**
    1. Delete static variables: `DustState`, `SeaSaltState`, `DryDepState`.
    2. Inside `ccpp_catchem_interface_init`, call `cc_wrap%catchem_model%initialize(...)` directly to load config and spawn the C++ central Core.
    3. Inside `ccpp_catchem_interface_run`, call `cc_wrap%catchem_model%bind_met_3d` and `bind_unified_chemistry` to bind standard CCPP double arrays to C++ views, then execute `catchem_core_run_timestep(core_ptr, dt)`.

### Task 3: Compile & Verification Check
*   **Goal:** Re-build the full codebase in parallel and run tests.
*   **Steps:**
    1. Run: `docker exec dazzling_lehmann cmake --build /workspace/build-cece --parallel`.
    2. Run: `docker exec dazzling_lehmann ctest --test-dir /workspace/build-cece --output-on-failure`.
