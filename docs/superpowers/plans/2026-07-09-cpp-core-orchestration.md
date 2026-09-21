# Implementation Plan: High-Performance C++ Core Orchestration and Legacy Fortran Elimination

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate all meteorological physical derivations, Gregorian calendar arithmetic, and unit conversion business logic from legacy Fortran modules, replacing them with standard C++20 and Kokkos implementations accessed via thin compatibility delegates.

**Architecture:** We use the C++ Delegate Wrapper Pattern. High-performance C++20 structures implement all state and calculation fields, and Fortran modules are empty wrappers forwarding arguments to standard `extern "C"` interfaces in `catchem_api.cpp`.

**Tech Stack:** C++20, Kokkos parallelism, standard Fortran ISO_C_BINDING memory interop.

## Global Constraints
* Language Target: Target C++20 utilizing Kokkos host views and C++20 standard-conforming mdspan backports, avoiding C++23 direct `<mdspan>` dependencies.
* Layout Alignment: Retain Fortran column-major storage layout (Kokkos::LayoutLeft / Kokkos::layout_left) for zero-copy CPU executions.
* Single Source of Truth: Core C++ handles 100% of orchestration, memory management, configuration loading, and diagnostics, completely bypassing legacy Fortran orchestration.
* Flat-Science Adapter Pattern: Cleanly map C++ unmanaged double views to Fortran array slices inside flat BIND(C) bridges (*ScienceBridge.F90) via standard c_f_pointer, avoiding duplicate physics code.
* Interoperable Kind Declarations: All integer boundaries in BIND(C) signatures must utilize explicit integer kinds, specifically `integer(c_int)`, to match standard C/C++ dimensions.
* Language Boundary Exception Checks: Wrap all BIND(C) export endpoints (catchem_api.cpp) in robust C++ try-catch blocks to prevent escaping C++ exceptions.

---

## Task List & Execution Schedule

### Task 1: Rewrite TimeState Arithmetic in C++
*   **Goal:** Re-implement year-month-day advances, leap checks, Julian dates, and holidays in C++.
*   **Target File:** `src/core/catchem_time_state.hpp`
*   **Steps:**
    1. Update the `TimeState` struct to include standard Gregorian arithmetic calculations and Julian date conversions.
    2. Add standard unit test assertions for the newly migrated calendar calculations inside `tests/test_TimeState.f90` (once wrapped).

### Task 2: Flat C-API Export Declarations for TimeState
* **Goal:** Register `TimeState` creation, destruction, advancement, and attribute retrieval endpoints under `extern "C"`.
* **Target File:** `src/core/catchem_api.cpp` / `src/core/catchem_api.hpp`
*   **Signatures to implement:**
    ```cpp
    extern "C" {
        void* catchem_time_state_create(int year, int month, int day, int hour, int minute, int second, double timestep);
        void catchem_time_state_destroy(void* ptr);
        void catchem_time_state_advance(void* ptr, double dt);
        int catchem_time_state_get_year(void* ptr);
        int catchem_time_state_get_month(void* ptr);
        int catchem_time_state_get_day(void* ptr);
        int catchem_time_state_get_hour(void* ptr);
        int catchem_time_state_get_minute(void* ptr);
        int catchem_time_state_get_second(void* ptr);
        double catchem_time_state_get_timestep(void* ptr);
        double catchem_time_state_get_julian_date(void* ptr);
        int catchem_time_state_get_doy(void* ptr);
        bool catchem_time_state_is_leap_year(int year);
        int catchem_time_state_get_days_in_month(int month, int year);
    }
    ```

### Task 3: Strip and Rewrite Fortran TimeState Proxy
*   **Goal:** Replace legacy Fortran calculations in `TimeState_Mod.F90` with direct forwarding calls to the new flat C-API.
*   **Target File:** `src/core/TimeState_Mod.F90`
*   **Steps:**
    1. Modify `type :: TimeStateType` to only contain `type(c_ptr) :: cpp_ptr`.
    2. Implement thin wrappers invoking `catchem_time_state_...` functions.
    3. Verify and compile using `docker exec dazzling_lehmann cmake --build /workspace/build-cece` and run `test_TimeState`.

### Task 4: Rewrite Concentration & Pressure Unit Conversions in C++
*   **Goal:** Port all unit conversions (including mass-to-volume, molecules/cm3-to-ppbv, and temperature-pressure calculations) to template header.
*   **Target File:** `src/core/catchem_unit_conversion.hpp`
*   **Steps:**
    1. Add standard functions to compute densities, molecular weights, and vectorized in-place array scaling in C++.

### Task 5: Strip and Rewrite Fortran UnitConversion Proxy
*   **Goal:** Delegate all Fortran unit conversions to C++.
*   **Target File:** `src/core/UnitConversion_Mod.F90`
*   **Steps:**
    1. Declare corresponding BIND(C) flat APIs inside `UnitConversion_Mod.F90`.
    2. Re-route the standalone conversion interfaces to delegate directly to the C++ core.
    3. Verify compilation and execute `test_UnitConversion` inside the development container.

---

## 5. Verification Checklist

*   [ ] Run `docker exec dazzling_lehmann cmake --build /workspace/build-cece --parallel`
*   [ ] Run `docker exec dazzling_lehmann ctest --test-dir /workspace/build-cece --output-on-failure`
*   [ ] Confirm 9/9 targets are completely green and passing.
