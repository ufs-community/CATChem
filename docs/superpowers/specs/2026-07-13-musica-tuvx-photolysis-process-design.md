# Technical Design Specification: MUSICA TUV-x Photolysis Process Integration

**Date:** July 13, 2026  
**Status:** Approved  
**Topic:** Integrating MUSICA TUV-x photolysis calculations into the modernized CATChem C++ framework as a native physical process.

---

## 1. Context & Objectives

In physical atmospheric models, photolysis is a critical physical transformation process that drives chemical kinetics by converting solar radiation into chemical reaction rates ($J$-values). In CATChem, photolysis represents one of the foundational transformation processes.

The **MUlti-Scale Infrastructure for Chemistry and Aerosols (MUSICA)** library, developed by NCAR and included under `src/external/musica/`, contains the standard **TUV-x** photolysis calculator engine.

The objective of this design is to integrate the MUSICA TUV-x engine directly into the modernized C++ CATChem core as a native C++ process (`PhotolysisProcess`) that extends the `catchem::ProcessInterface`.

---

## 2. Architecture & Integration Layout

Under the modernized CATChem core, all physical processes run within the sequential orchestrator loop `catchem::Core` on CPU host memory or GPU device views. To minimize bridging overhead, the photolysis process is implemented entirely in C++20 and registered with the dynamic `catchem::ProcessRegistry`.

```
                    +------------------------------------+
                    |        catchem::Core Orchestrator  |
                    +------------------------------------+
                                       |
                                       v (sequential timeline list)
                    +------------------------------------+
                    |      catchem::ProcessInterface     | (Virtual Base)
                    +------------------------------------+
                                       |
                                       v
                    +------------------------------------+
                    |     catchem::PhotolysisProcess     | (Concrete C++ Process)
                    +------------------------------------+
                                       |
                                       v (Column loop)
                    +------------------------------------+
                    |         musica::TUVX solver        | (External C++ Library)
                    +------------------------------------+
```

### 2.1 File Changes & Additions

A new directory `src/process/photolysis/` will be introduced to house the files:

*   `src/process/photolysis/CMakeLists.txt` - Build details linking the process and its source files to `musica_tuvx` and the CATChem core.
*   `src/process/photolysis/catchem_process_photolysis.hpp` - Header declaring the `PhotolysisProcess` class.
*   `src/process/photolysis/catchem_process_photolysis.cpp` - Detailed implementation of the process lifecycle hooks (`init`, `run`, `finalize`) and C++ registration symbols.

---

## 3. Detailed Specification

### 3.1 YAML Configuration Passing

To dynamically pass the TUV-x configuration file via the main configuration runtime (e.g. `catchem_config.yml`), the following updates will be applied:

1.  **Extend `catchem::StateManager`** with a public string member:
    ```cpp
    std::string config_file_path;
    ```
2.  **Populate `config_file_path`** in `catchem::Core` constructors:
    ```cpp
    state_mgr->config_file_path = config_file;
    ```
3.  **Define the photolysis block** in the main runtime config file:
    ```yaml
    process:
      photolysis:
        activate: true
        config_file: "src/external/musica/configs/tuvx/tuv_5_4.yml"
    ```

### 3.2 Class Declaration

```cpp
// src/process/photolysis/catchem_process_photolysis.hpp
#pragma once
#include "catchem_process_interface.hpp"
#include <musica/tuvx/tuvx.hpp>
#include <memory>
#include <string>

namespace catchem {

    class PhotolysisProcess : public ProcessInterface {
    private:
        std::string config_path;
        std::unique_ptr<musica::TUVX> tuvx_instance;
        musica::Mappings photo_mappings;

    public:
        PhotolysisProcess();
        ~PhotolysisProcess() override;

        std::string get_name() const override { return "photolysis"; }

        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override;
    };

} // namespace catchem
```

### 3.3 Initialization (`init`) Hook

In `init(state)`:
1.  Parse the main configuration file from `state->config_file_path` to extract the photolysis configuration block.
2.  Create instances of `GridMap`, `ProfileMap`, and `RadiatorMap` using the MUSICA API.
3.  Instantiate `musica::TUVX` and call `Create(config_path, ...)` to initialize the scientific TUV-x engine.
4.  Retrieve photolysis reaction listings using `GetPhotolysisRateConstantsOrdering` and map reaction names to internal solver indices.
5.  Dynamically register 2D midpoint-level diagnostic rate fields (e.g. `photolysis_rate_O3_O1D`) inside `DiagnosticManager`.

### 3.4 Execution (`run`) Hook

In `run(state)`:
1.  Synchronize active device views to CPU host memory using `state->sync_to_host()`.
2.  Obtain the `ProfileMap` from the `tuvx_instance` and fetch radiator profiles for `"air"`, `"O2"`, and `"O3"`.
3.  Iterate column-by-column across the horizontal grid. For each column:
    *   Calculate the Cosine of Solar Zenith Angle using `state->time.get_cos_sza(...)`.
    *   Accumulate vertical grid height edges in kilometers from `state->met.BXHEIGHT`.
    *   Populate vertical columns of air density, $O_2$ (standard fraction), and $O_3$ concentrations, converting to `molecules/cm3` units.
    *   Apply midpoint profile updates to TUV-x profiles using `SetMidpointValues(...)`.
    *   Execute the solver by calling `tuvx_instance->Run(...)`, receiving photolysis rates at grid edges.
    *   Interpolate edge values to layer midpoints:
        $$\text{rate\_midpoint}_i = \frac{\text{edge\_rate}_i + \text{edge\_rate}_{i+1}}{2}$$
    *   Populate the interpolated values directly into the registered diagnostic 2D views.
4.  Synchronize the execution outputs and diagnostics back to device memory via `state->sync_to_device()`.

---

## 4. Testing & Verification Plan

Verification of the new process will be handled by a dedicated C++ test file:
*   **Location:** `tests/test_catchem_photolysis.cpp`
*   **Action:** Add tests checking:
    1.  Successful process creation and registration in `catchem::ProcessRegistry`.
    2.  Successful parsing of configuration files and initialization of TUV-x.
    3.  Execution verification across simulated meteorological columns, validating non-zero interpolated photolysis rate diagnostics.
*   **Automation:** Integrate the test file into the `tests/CMakeLists.txt` suite.
