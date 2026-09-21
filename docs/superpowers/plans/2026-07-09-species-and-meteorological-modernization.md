# Implementation Plan: High-Performance C++ Core Modernization (Species and Meteorological Utilities)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fully migrate all chemical species properties parsing and meteorological equations from legacy Fortran files ('species_mod.F90', 'met_utilities_mod.F90') to the modern C++20 core, making Fortran zero-overhead delegation proxies.

**Architecture:** We use the C++ Delegate Wrapper Pattern. C++ `ChemState` and `met_utilities` namespace handle all metadata and math, and Fortran wraps them using standard interoperable C bindings.

**Tech Stack:** C++20, Kokkos parallelism, ISO_C_BINDING interop.

## Global Constraints
* Language Target: Target C++20 utilizing Kokkos host views and C++20 standard-conforming mdspan backports, avoiding C++23 direct `<mdspan>` dependencies.
* Single Source of Truth: Core C++ handles 100% of orchestration, memory management, configuration loading, and diagnostics, completely bypassing legacy Fortran orchestration.
* Language Boundary Exception Checks: Wrap all BIND(C) export endpoints (catchem_api.cpp) in robust C++ try-catch blocks to prevent escaping C++ exceptions.

---

## Task List & Execution Schedule

### Task 1: Extend SpeciesMetadata & YAML Parser in C++
*   **Goal:** Add wet/dry deposition, carbon chemical loss, and background properties to standard C++ loading.
*   **Target File:** `src/core/catchem_species_metadata.hpp`, `src/core/catchem_chem_state.hpp`
*   **Steps:**
    1. Add `henry_k0`, `henry_cr`, `henry_pKa`, `wd_reevap_frac`, `t_chem_loss`, and `BackgroundVV` fields to `SpeciesMetadata` inside `catchem_species_metadata.hpp`.
    2. Update `load_species_config` in `catchem_chem_state.hpp` to parse these fields from `CATChem_species.yml` with double underscore prefixes matching the schema.

### Task 2: Implement Missing Meteorological Utilities in C++
*   **Goal:** Re-implement specialized equations like Cunningham correction, Stokes settling velocity, dew point, relative humidity, Monin-Obukhov lengths, Businger stability, and solar zenith angle declination in C++.
*   **Target File:** `src/core/catchem_met_utilities.hpp`
*   **Steps:**
    1. Move high-precision, GPU-friendly physical calculation functions to standard template inline definitions under the `catchem::met_utilities` namespace.

### Task 3: Declare and Implement C-API Boundary Exports
*   **Goal:** Register all new species-database lookup methods and met-utilities physics routines under `extern "C"`.
*   **Target File:** `src/core/catchem_api.hpp` / `src/core/catchem_api.cpp`
*   **Steps:**
    1. Implement the wide range of species property lookups (MW, classifications, wet/dry dep variables) and physical calculation entry points using robust `try-catch` exception shields.

### Task 4: Strip and Re-route Fortran species_mod to C++
*   **Goal:** Convert `species_mod.F90` into a thin delegator wrapper.
*   **Target File:** `src/core/species_mod.F90`
*   **Steps:**
    1. Purge all local array allocations, loaders, and validation logic.
    2. Implement `populate_species_from_cpp(species, state_ptr, index)` to fetch all metadata once from C++ at startup.
    3. Update `SpeciesManagerType%load_from_file` to delegate directly to C++ loading.

### Task 5: Strip and Re-route Fortran met_utilities_mod to C++
*   **Goal:** Relocate all math and physical equations out of `met_utilities_mod.F90`, converting them into flat BIND(C) proxies.
*   **Target File:** `src/core/met_utilities_mod.F90`
*   **Steps:**
    1. Empty all mathematical formulas, constants, and loops.
    2. Interface and invoke their corresponding `catchem_met_...` C-API functions directly.

---

## 5. Verification Checklist

*   [ ] Run `docker exec dazzling_lehmann cmake --build /workspace/build-cece --parallel`
*   [ ] Run `docker exec dazzling_lehmann ctest --test-dir /workspace/build-cece --output-on-failure`
*   [ ] Verify 100% green status on all 9 compilation targets.
