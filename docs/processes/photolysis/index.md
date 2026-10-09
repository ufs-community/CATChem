# Modernized C++ Photolysis Process Documentation

## Overview

The modernized `PhotolysisProcess` is a native C++ physics process that implements radiative transfer calculations to determine atmospheric photolysis rates ($J$-values). It encapsulates NCAR's **TUV-x** (Tropospheric Ultraviolet and Visible) engine, which is part of NCAR's **musica** library submodule, completely replacing legacy Fortran proxy layers.

*   **Namespace / Class**: `catchem::PhotolysisProcess : public catchem::ProcessInterface`
*   **Header Location**: `src/process/photolysis/catchem_process_photolysis.hpp`
*   **Source Location**: `src/process/photolysis/catchem_process_photolysis.cpp`

---

## Architecture & How It Works

During a simulation timestep, `PhotolysisProcess::run` executes a column-wise photolysis calculation pipeline on the CPU host:

1.  **Device-to-Host Synchronization**: Pulls physical and meteorological variables from Kokkos GPU memory space back to CPU buffers:
    ```cpp
    state->sync_to_host();
    ```
2.  **Column-Wise Input Extraction**: Loops over 1D columns, extracting key meteorological fields for every vertical level:
    *   Midpoint Altitudes and Grid-cell Thicknesses
    *   Air Density (`AIRDEN_DRY`)
    *   Solar Zenith Angle (SZA)
    *   Ozone ($O_3$) profile density
3.  **Solver Execution**: Delegates the 1D columns to the scientific TUV-x engine:
    ```cpp
    tuvx_instance->Calculate(photo_mappings, columns);
    ```
4.  **Edge-to-Midpoint Interpolation**: The scientific TUV-x engine calculates photolysis rates at grid-cell edges. `PhotolysisProcess` performs linear/exponential interpolation to map these edge rates to grid-cell midpoints (layer centers) to match CATChem's chemistry layout.
5.  **Diagnostic Publishing**: Midpoint photolysis rates are registered dynamically and written into the global `DiagnosticManager` under the naming convention:
    $$\text{diagnostic\_name} = \text{"photolysis\_rate\_"} + \text{rx\_name}$$
    (e.g., `"photolysis_rate_jfoo"` for reaction `jfoo`).
6.  **Coupling Hand-Off**: These diagnosed photolysis rates are then automatically consumed by downstream chemistry processes (like `GasChemProcess`), requiring zero copies or hardcoding.

---

## Configuration

The process is configured through CATChem's YAML system:

```yaml
processes:
  - name: "photolysis"
    enabled: true
    parameters:
      config_file: "src/external/musica/configs/tuvx/tuv_5_4.yml"
    diagnostics:
      - "photolysis_rate_jfoo"
```

---

## Integration and Verification

Like all modernized processes under the C++ core, the photolysis process is tested and verified under the standard CTest suite. Its dynamic coupling to chemistry is end-to-end verified by `test_catchem_gaschem` to ensure that midpoint-interpolated $J$-values are stable, non-negative, and properly scaled.
