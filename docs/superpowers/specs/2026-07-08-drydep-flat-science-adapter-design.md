# Spec: C++ to Flat-Fortran Science Adapter for DryDep

## 1. Overview
The goal of this design is to fully modernize the **`DryDep` (Dry Deposition)** process. In the original legacy codebase, `DryDep` relies on complex Fortran abstractions (`ProcessDryDepInterface_Mod.F90` and `DryDepProcessCreator_Mod.F90`), which require the legacy Fortran StateManager, VirtualColumns, and auto-generated meteorological macros to run.

This design implements a **C++ to Flat-Fortran Science Adapter (Approach 1)**. It completely bypasses the legacy Fortran Core. Instead, the centralized C++ Core (`catchem::Core` and `catchem::StateManager`) acts as the single source of truth. When running dry deposition, the C++ wrapper retrieves raw host pointers of its Kokkos Views (and dynamic diagnostic fields), and passes them straight to a thin, C-linkable Fortran bridge (**`DryDepScienceBridge.F90`**). The bridge standardizes the raw pointers as native Fortran arrays, loops over columns, and slices them to dispatch unmodified physical calculations (`compute_wesely`, `compute_gocart`, or `compute_zhang`) in-place.

This pattern preserves **100% of the untouched science modules under `schemes/`**, eliminates duplicate core libraries, and achieves high-performance zero-copy execution on CPU targets.

---

## 2. Component Layout & Files Destiny

The modernization refactors the `src/process/drydep/` directory to have the following file destiny:

| File | Status | Action / Role |
| :--- | :--- | :--- |
| **`schemes/DryDepScheme_WESELY_Mod.F90`** | Retained | Unmodified scientific solver (`compute_wesely`) |
| **`schemes/DryDepScheme_GOCART_Mod.F90`** | Retained | Unmodified scientific solver (`compute_gocart`) |
| **`schemes/DryDepScheme_ZHANG_Mod.F90`** | Retained | Unmodified scientific solver (`compute_zhang`) |
| **`DryDepCommon_Mod.F90`** | Retained | Retains standard types like `DryDepProcessConfig` used in YAML parsing |
| **`DryDepScienceBridge.F90`** | **NEW** | C-bound entrypoint mapping C++ pointers to Fortran slices and executing the column loops |
| **`catchem_process_drydep.hpp`** | Modified | Redesigned standard C++ `ProcessInterface` declaration |
| **`catchem_process_drydep.cpp`** | Modified | Unified C++ process caller pulling Views, registering diagnostics, and invoking the C-Bridge |
| **`ProcessDryDepInterface_Mod.F90`** | **DELETED** | Removed to decouple from legacy `VirtualColumn` and `StateManager` |
| **`DryDepProcessCreator_Mod.F90`** | **DELETED** | Removed as registration is now handled on the C++ side |
| **`CMakeLists.txt`** | Modified | Simplifies targets, links with C++ code, and removes deleted Fortran dependencies |

---

## 3. Data Flow and Interface Definitions

### A. The C-Linkable Bridge: `DryDepScienceBridge.F90`
The bridge defines a standard, flat `BIND(C)` routine that receives the raw C++ memory pointers, constructs Fortran slice pointers, and loops over grid columns:

```fortran
subroutine run_drydep_science_bridge( &
   n_cols, n_levels, n_species, dt, &
   gas_scheme, aero_scheme, diagnostics, &
   ! 3D Met Pointers
   c_bxheight, c_airden, c_t_air, c_z_edges, c_rh, &
   ! 2D/1D Met Pointers
   c_cldfrc, c_frlai, c_frlanduse, c_iland, c_is_ice, c_is_land, c_is_snow, &
   c_lat, c_lon, c_obk, c_ps, c_salinity, c_suncosmid, c_swgdn, c_ts, c_tskin, &
   c_ustar, c_z0, c_frlake, c_gwettop, c_hflux, c_lwi, c_pblh, c_u10m, c_v10m, c_z0h, &
   ! Species Metadata Arrays
   species_names, species_mw_g, species_dd_f0, species_dd_hstar, &
   species_dd_DvzAerSnow, species_dd_DvzMinVal_snow, species_dd_DvzMinVal_land, &
   species_density, species_radius, species_is_seasalt, species_is_dust, &
   species_lower_radius, species_upper_radius, is_gas, &
   ! Concentrations, Tendencies & Diagnostics
   c_conc, c_tendency, c_diag_con, c_diag_vel, &
   diagnostic_species_id, n_diag_species &
) bind(C, name="run_drydep_science_bridge")
```

#### Mapping & Slicing inside the loop:
1. Reconstruct multi-dimensional Fortran pointer slices using `c_f_pointer`:
   ```fortran
   call c_f_pointer(c_conc, conc, [n_cols, n_levels, n_species])
   call c_f_pointer(c_diag_con, diag_con, [n_cols, n_species])
   ```
2. Loop over columns (`do icol = 1, n_cols`):
   * Feed slices straight into unmodified science routines:
     ```fortran
     call compute_wesely(..., conc(icol, 1, :), col_tendencies, ...)
     ```
   * Sliced diagnostics outputs `diag_con(icol, :)` are written directly into the C++ `DiagnosticManager` host buffer in-place.
   * Apply computed tendencies in-place:
     ```fortran
     conc(icol, 1, :) = conc(icol, 1, :) + dt * col_tendencies(1, :)
     ```

---

### B. Modernized C++ Process: `catchem::DryDepProcess`

#### Class Declaration:
```cpp
namespace catchem {

class DryDepProcess : public ProcessInterface {
private:
    std::string gas_scheme;
    std::string aero_scheme;
    bool diagnostics_enabled;
    std::vector<int> diagnostic_species_id;

public:
    DryDepProcess();
    std::string get_name() const override { return "drydep"; }
    void init(std::shared_ptr<StateManager> state) override;
    void run(std::shared_ptr<StateManager> state) override;
    void finalize() override;
};

} // namespace catchem
```

#### Zero-Copy Execution pipeline in `run()`:
1. Performs **`state->sync_to_host()`** to ensure the latest Concentrations and MetState inputs are synchronized to the host CPU array heap.
2. Extracts raw host pointer addresses from central `InteropField` View buffers:
   ```cpp
   double* t_ptr = state->met.T->host_view.data();
   double* conc_ptr = state->chem.conc->host_view.data();
   ```
3. Dynamically fetches raw pointers of C++ allocated diagnostics:
   ```cpp
   double* diag_con = (double*)state->diag_mgr->get_host_pointer("drydep_con_per_species");
   ```
4. Invokes **`run_drydep_science_bridge`** with standard C/C++ primitives.
5. Performs **`state->sync_to_device()`** to flush updated concentrations and populated diagnostics to the active GPU execution View mirrors for subsequent Kokkos solvers.

---

## 4. Error Handling and Verification

### A. Input Verification and Safe Checks
* **Dimension Checks:** The Fortran bridge validates dimensions (`n_cols`, `n_levels`, `n_species`) against expected YAML configuration to prevent out-of-bounds segfaults.
* **Target Pointer Integrity:** Pointers received via `c_ptr` are verified via standard `c_associated()` before calling `c_f_pointer`.

### B. Regression & Equivalence Testing
* **Numeric Integrity:** The integration verification runs all unported schemes (`Wesely`, `GOCART`, `Zhang`) side-by-side inside the compiled test suites and asserts bit-for-bit identity against reference values.

---

## 5. Timeline & Plan

1. **Create `DryDepScienceBridge.F90`** containing the flat BIND(C) dynamic router and array slicing.
2. **Rewrite `catchem_process_drydep.hpp` & `.cpp`** to support configuration checking, diagnostics allocation, and direct pointer extraction.
3. **Delete legacy wrappers** `ProcessDryDepInterface_Mod.F90` and `DryDepProcessCreator_Mod.F90`.
4. **Modify `src/process/drydep/CMakeLists.txt`** to wire up the new source files and targets.
5. **Compile and execute verification** in Docker to confirm 100% correct, zero-copy interop.
