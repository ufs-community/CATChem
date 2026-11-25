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
    *   Number and layer column densities for tagged gas species
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

In the CATChem species file, give each gas species supplied by the host a
`photolysis.<profile>` role. The suffix must match the TUV-x profile name,
including case. For example, use `roles: [photolysis.NO2]` to feed the host
nitrogen dioxide mixing ratio into the TUV-x profile named `NO2`. CATChem
creates the profile if it is absent from the TUV-x configuration. An explicitly
configured profile for that role must be host-updatable and use
`molecule cm-3`.

CATChem converts dry-air ppmv to molecules per cubic centimeter using local dry-air
density and supplies molecules per square centimeter for each layer using its
thickness. Negative mixing ratios contribute zero absorption, and the original
concentration remains available to gas chemistry. The legacy `photolysis.ozone`
role maps to `O3` and preserves any fixed ozone profile in an existing TUV-x
configuration. For untagged `O2` and `O3`, CATChem uses configured profiles or
background mixing ratios of 20.95% and 0.3 ppmv, respectively.

TUV-x reaction names are separate from profile names. For MICM to consume a
photolysis rate, its `PHOTO.<label>` must match the TUV-x reaction label.

Temperature-dependent quantum yields use temperatures at layer interfaces. A
host can supply `TEDGE` in kelvin through the axis-aware meteorology binding
with the interface axis and one more value per column than `T`. When `TEDGE`
is absent from the current input, CATChem linearly interpolates `T` in height
at interior interfaces and uses the nearest layer-center temperature at the
bottom and top. A fixed temperature profile from the TUV-x configuration
keeps its configured values.

---

## Integration and Verification

Like all modernized processes under the C++ core, the photolysis process is tested and verified under the standard CTest suite. Its dynamic coupling to chemistry is end-to-end verified by `test_catchem_gaschem` to ensure that midpoint-interpolated $J$-values are stable, non-negative, and properly scaled.
