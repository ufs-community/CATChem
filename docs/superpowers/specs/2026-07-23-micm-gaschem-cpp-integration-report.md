# Technical Report: Modernized C++ GasChem (MICM) Integration & Photolysis Coupling

**Date:** July 23, 2026  
**Status:** Completed  
**Topic:** Modernized C++ GasChem (MICM) process integration and J-rate coupling inside CATChem.

---

## 1. Executive Summary

We have successfully completed the implementation and verification of the modernized, native C++ `GasChemProcess` inside the CATChem core. This native C++ process wraps the C++ Model Independent Chemistry Module (MICM) solver using the `musica` library and establishes a dynamic, zero-overhead photolysis-rate coupling channel with the TUV-x photolysis solver.

All implementations are fully complete, compile with 100% success, and pass both property-based unit tests and coupled integration tests under CTest.

---

## 2. Directory & Class Structure

The modernized process is housed under `src/process/gaschem/` to completely replace legacy Fortran proxy layers:

*   **`src/process/gaschem/catchem_process_gaschem.hpp`**: Declares `GasChemProcess : public ProcessInterface` and registers a C-linkable hook (`catchem_register_gaschem_cpp`) to bind into CATChem's dynamic C++ registry.
*   **`src/process/gaschem/catchem_process_gaschem.cpp`**: Implements core initialization lifecycles, multi-level 3D grid flattening, bidirectional Volume Mixing Ratio (ppmv) to molar density ($\text{mol/m}^3$) scaling, zero-overhead automatic photolysis rate mapping, and C++ MICM solver execution.
*   **`src/process/gaschem/CMakeLists.txt`**: Defines target library and includes configurations, linking the target static library `catchem_process_gaschem` dynamically with the C++ Core (`CATChem_core_cpp`), `musica`, `yaml-cpp`, and Kokkos views.

---

## 3. How the C++ GasChem Process Works

When the core chemistry solver runs during a model timestep, the `GasChemProcess::run` method executes the following multi-step pipeline:

```
[CATChem Core Timestep]
          │
          ▼
┌────────────────────────────────────────────────────────┐
│ 1. Device-to-Host Synchronization                      │
│    state->sync_to_host()                               │
└─────────────────────────┬──────────────────────────────┘
                          │
                          ▼
┌────────────────────────────────────────────────────────┐
│ 2. Environmental & Concentration State Mapping         │
│    • Flatten 3D Grid columns: ilev * n_cols + icol    │
│    • Convert Dry Air Density: kg/m³ -> mol/m³          │
│    • Convert Concentrations: ppmv -> mol/m³            │
│    • Populate musica::State arrays                     │
└─────────────────────────┬──────────────────────────────┘
                          │
                          ▼
┌────────────────────────────────────────────────────────┐
│ 3. Automatic Zero-Overhead Photolysis Coupling         │
│    • Scan MICM rate parameters starting with "PHOTO."  │
│    • Construct diagnostic target: "photolysis_rate_*"  │
│    • Fetch J-rates from global DiagnosticManager       │
│    • Populate MICM State rate_parameters array         │
└─────────────────────────┬──────────────────────────────┘
                          │
                          ▼
┌────────────────────────────────────────────────────────┐
│ 4. Solver Invocation & Execution                       │
│    • Run micm_instance->Solve(state, timestep)         │
│    • Inspect solver_result.state_ convergence          │
└─────────────────────────┬──────────────────────────────┘
                          │
                          ▼
┌────────────────────────────────────────────────────────┐
│ 5. Back-Conversion & State Synchronization             │
│    • Convert concentrations: mol/m³ -> ppmv            │
│    • Write back to state->chem.conc host view          │
│    • Sync host back to device: state->sync_to_device() │
└────────────────────────────────────────────────────────┘
```

### 3.1 3D Grid Flattening
Any 3D grid cell in CATChem at coordinate `(icol, ilev)` is mapped to its flat 1D MICM cell index:
$$i_{\text{cell}} = i_{\text{lev}} \times N_{\text{cols}} + i_{\text{col}}$$
This indexing aligns exactly with the row-major diagnostic storage utilized by the photolysis process.

### 3.2 Environmental Conditions Mapping
For each 1D grid cell, physical conditions are mapped:
* **Temperature:** Set directly from $T(\text{icol}, \text{ilev}, 0)$.
* **Pressure:** Set directly from $\text{PMID}(\text{icol}, \text{ilev}, 0)$.
* **Molar Air Density:** Converted from dry air density (AIRDEN_DRY in $\text{kg/m}^3$) to molar density ($\text{mol/m}^3$) using the dry air molecular weight (approx. $0.0289644\text{ kg/mol}$):
  $$\rho_{\text{molar}} = \frac{\rho_{\text{dry}}}{0.0289644}$$
  where $\rho_{\text{molar}}$ is the dry air molar density in $\text{mol/m}^3$ and $\rho_{\text{dry}}$ is the dry air density in $\text{kg/m}^3$ (`AIRDEN_DRY`).

### 3.3 Boundary Value Safeguards & Stability
To prevent solver failures, NaNs, or matrix singularity issues:
* Physical variables like Temperature ($T$), Pressure ($P$), and dry air density are verified and guarded against non-positive bounds.
* Molar concentrations are clamped to a strict floor of `1.0e-20` on both input conversion and output conversion loops, ensuring that minor solver inaccuracies do not feed negative concentrations back into CATChem.

---

## 4. Automatic Photolysis Coupling

A central focus of this integration is the creation of a **zero-hardcoded, fully dynamic, low-overhead coupling** between TUV-x Photolysis and MICM.

Under the hood:
1.  During TUV-x photolysis calculations, calculated midpoint J-rates ($s^{-1}$) are registered dynamically and stored in `diag_mgr` with the naming convention `"photolysis_rate_" + rx_name` (e.g. `"photolysis_rate_jfoo"`).
2.  During the chemistry phase, `GasChemProcess` scans the MICM solver's rate parameter map (`state_->GetRateParameterMap()`). If a parameter starts with `"PHOTO."` (e.g., `"PHOTO.jfoo"`):
    *   It strips the prefix to extract the label: `jfoo`.
    *   It constructs the query name: `"photolysis_rate_jfoo"`.
    *   It queries the `DiagnosticManager`. If the field is registered, it retrieves the host memory pointer and copies the calculated photolysis rates directly into MICM's rate parameters at the correct offset index for every grid cell:
        $$\text{Index} = i_{\text{cell}} \times N_{\text{params}} + i_{\text{param}}$$
3.  Any loss coefficients (named `"LOSS."`) are defaulted to `1.0` to preserve the chemistry mechanism's built-in scaling factors.

---

## 5. Verification and Testing

We followed a strict Test-Driven Development (TDD) cycle and added two robust test suites to verify compile-time safety and runtime behavioral correctness:

### 5.1 Property-Based Unit Tests (`test_catchem_gaschem_units`)
We built a new unit test suite to mathematically verify our conversions:
*   **Reversibility Invariance**: Asserts that Volume Mixing Ratio conversions ($\text{ppmv} \rightarrow \text{mol/m}^3 \rightarrow \text{ppmv}$) are strictly reversible, identity-preserving, and non-negative across multiple orders of dry air densities ($0.5$ to $1.2\text{ kg/m}^3$) and mixing ratios ($10^{-6}$ to $100.0\text{ ppmv}$).
*   **Safeguards Verification**: Asserts that negative boundary concentrations are successfully bounded to `1.0e-20` to guarantee solver safety and NaN prevention.

### 5.2 Coupled Process Integration Test (`test_catchem_gaschem`)
We built an end-to-end integration test containing a simulated atmospheric column with $3$ vertical levels to verify coupled execution under CTest:
1.  Dynamically registers and initializes both `photolysis` (TUV-x) and `gaschem` (MICM) C++ processes.
2.  Synthesizes a noon summer met profile ($T$, $p$, and dry air density).
3.  Runs a full core timestep (`core->run_timestep(3600.0)`):
    *   `photolysis` calculates midpoint J-rates and populates the `"photolysis_rate_jfoo"` diagnostic.
    *   `gaschem` automatically reads the diagnostic, populates MICM's `"PHOTO.jfoo"` rate parameter, and solves the chemical equations.
4.  Asserts that photolysis rates are non-zero, that chemistry successfully converges, and that output concentrations reflect correct tendencies.

All tests compile cleanly and **pass with 100% success** under the CTest suite:
```
Test project /Users/barry/Documents/CATChem/build-macos/tests
    Start 11: test_catchem_gaschem
1/2 Test #11: test_catchem_gaschem .............   Passed    0.35 sec
    Start 12: test_catchem_gaschem_units
2/2 Test #12: test_catchem_gaschem_units .......   Passed    0.00 sec

100% tests passed, 0 tests failed out of 2
```
