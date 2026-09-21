# Technical Design Specification: Modernized C++ GasChem (MICM) Process Integration

**Date:** July 23, 2026  
**Status:** Under Review  
**Topic:** Integrating the legacy MICM solver (GasChem process) into the modernized C++ core with automatic photolysis rate coupling.

---

## 1. Executive Summary

This design specification details the porting and modernization of the GasChem (Model Independent Chemistry Module - MICM) process from its legacy Fortran wrapper framework (referencing PR #160) into a native C++ implementation within CATChem's modern core.

A primary objective of this integration is the creation of a direct, dynamic coupling channel between the newly-modernized **TUV-x Photolysis Process** and the **GasChem MICM solver**. Under this architecture, the C++ GasChem process dynamically maps MICM's `"PHOTO.<label>"` rate parameters to the photolysis midpoint $J$-rate diagnostics computed and stored in the global `DiagnosticManager` as `"photolysis_rate_<label>"`. This achieves a zero-hardcoding, fully dynamic, low-overhead coupling that supports any photolysis reactions defined in the chemistry mechanism.

---

## 2. Architecture & Component Design

The native C++ GasChem process will be implemented as `catchem::GasChemProcess`, deriving from the base `catchem::ProcessInterface`.

```cpp
namespace catchem {

    class GasChemProcess : public ProcessInterface {
    private:
        std::string config_dir;
        std::unique_ptr<musica::MICM> micm_instance;
        std::unique_ptr<musica::State> micm_state;
        bool initialized = false;

    public:
        GasChemProcess();
        ~GasChemProcess() override;

        std::string get_name() const override { return "gaschem"; }

        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override;
    };

} // namespace catchem
```

### 2.1 Initialization Lifecycle
During `init(state)`:
1. **Config Directory Discovery:** The process inspects the main configuration path (`state->config_file_path`). It extracts the parent directory of this file (e.g. `tests/Configs/Default/`) to use as the base directory for MICM.
2. **MICM Solver Instantiation:** The process instantiates a `musica::MICM` object with this directory path and the standard-ordered Rosenbrock solver type (`musica::RosenbrockStandardOrder`). This dynamically parses `micm_config.yaml`, `phases.yaml`, and `reactions.yaml` located in that directory.
3. **State Allocation:** It allocates a `musica::State` instance for exactly $N_{\text{cells}} = \text{n\_cols} \times \text{n\_levels}$ grid cells.

---

## 3. Data Flow & Unit Conversions

The native solver executes on the host. To support this execution, CATChem's Kokkos views are synchronized to the host, and grid attributes are flattened to a 1D sequence of size `n_cols * n_levels`.

### 3.1 3D Grid Flattening
Any 3D grid cell in CATChem at coordinate `(icol, ilev)` is mapped to its flat 1D MICM cell index:
$$\text{i\_cell} = \text{ilev} \times \text{n\_cols} + \text{icol}$$
This indexing aligns exactly with the row-major diagnostic storage utilized by the photolysis process.

### 3.2 Environmental Conditions Mapping
For each 1D grid cell, physical conditions are mapped:
* **Temperature:** Set directly from $T(\text{icol}, \text{ilev}, 0)$.
* **Pressure:** Set directly from $\text{PMID}(\text{icol}, \text{ilev}, 0)$.
* **Molar Air Density:** Converted from dry air density (AIRDEN_DRY in $\text{kg/m}^3$) to molar density ($\text{mol/m}^3$) using the dry air molecular weight (approx. $0.0289644\text{ kg/mol}$):
  $$\text{air\_density\_mol\_m}^3 = \frac{\text{AIRDEN\_DRY}(\text{icol}, \text{ilev}, 0)}{0.0289644}$$

### 3.3 Species Concentrations Mapping
Species concentrations are mapped bidirectionally:
* **Input Conversion:** Concentrations in CATChem are stored in volume mixing ratio (ppmv). Before calling the solver, they are converted to molecular number density ($\text{mol/m}^3$):
  $$\text{conc\_mol\_m}^3 = \text{conc\_ppmv} \times 10^{-6} \times \text{air\_density\_mol\_m}^3$$
  The resulting value is stored in MICM's flat ordered concentrations vector at index:
  $$\text{idx} = \text{i\_cell} \times N_{\text{species}} + \text{i\_micm\_spec}$$
* **Output Conversion:** After solver execution, the final concentrations are converted back to mixing ratio:
  $$\text{conc\_ppmv} = \frac{\text{conc\_mol\_m}^3}{\text{air\_density\_mol\_m}^3} \times 10^6$$
  and written back into `state->chem.conc`.

---

## 4. Automatic Photolysis Coupling

Dynamic coupling between the photolysis and gaschem processes is handled automatically.

During the run phase, the `GasChemProcess` iterates through the MICM rate parameter map (`state_->GetRateParameterMap()`). If a rate parameter matches the pattern `"PHOTO.<label>"`:
1. It extracts the reaction `label` (e.g. `PHOTO.jfoo` $\rightarrow$ `jfoo`).
2. It constructs the corresponding diagnostic field name: `"photolysis_rate_" + label`.
3. It queries the `DiagnosticManager` for this field name. If the field is registered, it retrieves its midpoint J-rate diagnostic pointer on the host.
4. For each grid cell, the photolysis rate ($s^{-1}$) is read from index `ilev * n_cols + icol` of the diagnostic array and copied into MICM's rate parameters at:
   $$\text{idx} = \text{i\_cell} \times N_{\text{rate\_params}} + \text{i\_param}$$
5. If the parameter is a loss coefficient starting with `"LOSS."`, it is defaulted to $1.0$ to preserve the mechanism's built-in scaling factors.

---

## 5. Testing Strategy

We follow a rigorous development lifecycle utilizing both unit, property-based, and coupled integration tests.

### 5.1 Unit & Property-Based Tests
A new test suite under `tests/test_catchem_gaschem_units.cpp` will verify core mathematical operations and boundary conditions:
* **Molar Density Calculations:** Verify that the dry air density to molar density conversion matches physical laws for standard atmospheric profiles.
* **VMR-to-Molar Concentration Conversions:** Property-based checks ensuring that concentration conversions (ppmv $\leftrightarrow \text{mol/m}^3$) are strictly reversible, identity-preserving, and non-negative.
* **Grid Mapping Invariance:** Confirm that grid cell flattening to 1D index is invariant under varying level-column sizes, and that boundaries are strictly enforced.
* **Negative Value Safeguards:** Property test checking that if negative values exist in input concentrations, they are either safely bounded or handled gracefully without solver NaN propagation.

### 5.2 Coupled Integration Test (`test_catchem_gaschem.cpp`)
An end-to-end integration test will be registered under CTest to verify coupled execution:
1. **Registration:** Dynamically register both `photolysis` and `gaschem` C++ processes.
2. **Setup:** Set up a 1D column with 3 levels, and populate profiles ($T, p, \text{airden\_dry}$).
3. **Execution:** Add both processes to `Core` and run a timestep.
   * `photolysis` calculates midpoint J-rates and populates the `photolysis_rate_jfoo` diagnostic.
   * `gaschem` automatically reads the diagnostic, populates MICM's `PHOTO.jfoo` rate parameter, and solves the chemical equations.
4. **Validation:** Check that:
   * Photolysis rates calculated are non-zero and finite.
   * MICM solver accepts the rates, executes, and advances chemistry successfully.
   * Final species concentrations reflect the correct chemical tendencies.

---

## 6. Build System & Directory Changes

To compile and link the native C++ GasChem process:

* **Add `src/process/gaschem/` Directory:**
  * `CMakeLists.txt` will define the static library `catchem_process_gaschem`, target headers, and link dependencies on `CATChem_core_cpp`, `musica`, `yaml-cpp`, and `Kokkos::kokkos`.
* **Update `src/process/CMakeLists.txt`:**
  * Add `add_subdirectory(gaschem)`.
* **Update `tests/CMakeLists.txt`:**
  * Declare `test_catchem_gaschem` and link to `catchem_process_gaschem` and `catchem_process_photolysis`.
  * Add the test to the CTest suite.
