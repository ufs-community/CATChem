# Modernized C++ Gas-Phase Chemistry (GasChem) Process Documentation

## Overview

The modernized `GasChemProcess` is a native C++ physics/chemistry process that acts as a core atmospheric transformation node in CATChem. It wraps the C++ Model Independent Chemistry Module (**MICM**) solver using NCAR's **musica** library, completely replacing the legacy Fortran proxy layers and double-buffering.

*   **Namespace / Class**: `catchem::GasChemProcess : public catchem::ProcessInterface`
*   **Header Location**: `src/process/gaschem/catchem_process_gaschem.hpp`
*   **Source Location**: `src/process/gaschem/catchem_process_gaschem.cpp`

---

## Architecture & How It Works

During a chemistry simulation timestep, `GasChemProcess::run` executes a highly-vectorized, multi-step pipeline on the CPU host:

1.  **Device-to-Host Synchronization**: Retrieves active meteorological and concentration fields from Kokkos device memories:
    ```cpp
    state->sync_to_host();
    ```
2.  **3D Grid Column Flattening**: Converts CATChem's row-major 3D grid layouts into a contiguous 1D array structure suited for MICM:
    $$i_{\text{cell}} = i_{\text{lev}} \times N_{\text{cols}} + i_{\text{col}}$$
3.  **Physical Boundary Safeguards**: Sanitizes input temperature, pressure, and dry air density fields to enforce positive physical bounds.
4.  **Environmental State Mapping**:
    *   Converts dry air density (`AIRDEN_DRY` in $\text{kg/m}^3$) to dry air molar density ($\text{mol/m}^3$) using the molecular weight of dry air ($0.0289644\text{ kg/mol}$):
        $$\rho_{\text{molar}} = \frac{\rho_{\text{dry}}}{0.0289644}$$
    *   Converts gas concentrations from Volume Mixing Ratios (VMR, $\text{ppmv}$) to molar densities ($\text{mol/m}^3$) using dry air molar density:
        $$C_{\text{molar}} = C_{\text{ppmv}} \times 10^{-6} \times \rho_{\text{molar}}$$
5.  **Clamping Stability Bounds**: To prevent solver failure or NaN propagation under extreme conditions, all species concentrations are clamped to a hard floor of `1.0e-20` on both input and output loops.
6.  **Solver Execution**: Runs the C++ MICM Rosenbrock standard-ordered solver:
    ```cpp
    micm_instance->Solve(state, timestep);
    ```
7.  **Back-Conversion & Synchronization**: Converts output concentrations back to VMR ($\text{ppmv}$), updates the host state views, and flushes data to device memories:
    ```cpp
    state->sync_to_device();
    ```

---

## Automatic Zero-Overhead Photolysis Coupling

A standout feature of the C++ core is the dynamic, zero-hardcoded coupling channel with the C++ Photolysis process:

1.  **Scan Rate Parameters**: During initialization and run phases, `GasChemProcess` iterates through the solver's rate parameter map.
2.  **Prefix Extraction**: If a rate parameter matches the pattern `"PHOTO.<label>"` (e.g., `"PHOTO.jfoo"`):
    *   It strips the prefix to extract the raw photolysis reaction identifier: `jfoo`.
    *   It constructs the diagnostic search key: `"photolysis_rate_jfoo"`.
3.  **Dynamic Memory Binding**: It queries the global `DiagnosticManager` for `"photolysis_rate_jfoo"`. If registered, it obtains the raw pointer address and copies the calculated photolysis J-rates directly into MICM's rate parameter array at zero-copy speeds.
4.  **Loss Coefficient Defaults**: Any mechanism loss coefficients matching `"LOSS.<label>"` are initialized to `1.0` to preserve the chemistry mechanism's built-in scaling factors.

---

## Configuration

The process is dynamically added and configured via YAML layout files.

```yaml
processes:
  - name: "gaschem"
    enabled: true
    parameters:
      config_dir: "src/external/musica/configs/tuvx/from_host/"
    diagnostics:
      - "photolysis_rate_jfoo"
```

---

## Verification & Testing Suite

We employ a strict Test-Driven Development (TDD) harness containing two robust verification targets:

### 1. Mathematical Invariant Unit Tests (`test_catchem_gaschem_units`)
Verifies mathematical conversion reversibility and clamping safety across a wide physical envelope:
*   Asserts that $(\text{ppmv} \rightarrow \text{mol/m}^3 \rightarrow \text{ppmv})$ is strictly identity-preserving ($< 1.0\times 10^{-15}$ relative tolerance).
*   Validates conversion safety across a broad range of dry air densities ($0.5$ to $1.2\text{ kg/m}^3$) and mixing ratios ($10^{-6}$ to $100.0\text{ ppmv}$).
*   Enforces concentration flooring of `1.0e-20` on negative out-of-bound inputs.

### 2. End-to-End Simulation Tests (`test_catchem_gaschem`)
Runs a coupled multi-level atmospheric column simulation containing both the C++ `photolysis` and `gaschem` processes:
1.  Registers both processes dynamically inside a C++ `Core` instance.
2.  Synthesizes physical conditions modeling a noon summer atmospheric column.
3.  Executes `core->run_timestep(3600.0)`.
4.  Asserts that photolysis rate diagnostics are successfully populated and mapped to MICM, and that solver concentrations successfully converge with stable mass balance.
