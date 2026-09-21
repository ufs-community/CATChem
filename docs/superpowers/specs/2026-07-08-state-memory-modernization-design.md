# Technical Design Specification: Phase 5 — State and Memory Management Modernization

**Date:** July 8, 2026  
**Status:** Draft  
**Topic:** Transitioning CATChem State, Meteorology, Chemistry, and Time/Solar tracking systems completely to C++20 and Kokkos, establishing zero-copy interfaces for both standalone and coupled execution.

---

## 1. Context & Objectives

In standard coupled modeling configurations, CATChem acts as an incorporated chemical sub-model within a host climate or weather forecast system (e.g., FV3/UFS via CCPP or NUOPC). Consequently, the majority of meteorological state fields are allocated externally by the host model and passed to CATChem as raw pointers.

In legacy CATChem:
*   Meteorological arrays were maintained as direct allocations across 50+ polymorphic Fortran state variables, requiring custom macro generation during compile time.
*   The chemical species were fragmented across separate allocated 3D heaps (`ChemSpecies(s)%conc(:,:,:)`), preventing uniform GPU coalescence.
*   Thermodynamic derivations and SZA solar time tracking calculations were bounded to CPU executions in Fortran.

### Core Objectives:
*   **Dual-Mode Meteorological State (`MetState`):** Establish a structured C++ `MetState` container holding optional, lazy-allocated `InteropField` objects. It must seamlessly run in **Bind mode** (wrapping externally allocated raw Fortran pointers without memory copies) or **Allocate mode** (allocating managed device memory on demand for standalone runs).
*   **Unified Chemistry State View:** Consolidate the entire chemistry state into a single contiguous 3D Kokkos View `(n_cols, n_levels, n_species)` in C++ mapping directly to a standard 4D Fortran array `(nx, ny, nz, n_species)` for coalesced execution.
*   **Performance Portable Time & Solar Engine:** Transition Julian date tracking and high-precision Solar Zenith Angle computations to portable C++ structures decorated with `KOKKOS_INLINE_FUNCTION`.
*   **Thread-Safe Parallel Field Derivations:** Leverage Kokkos `parallel_for` kernels to compute derived met fields (e.g. dry air density, virtual temperature, box heights) on active acceleration targets.

---

## 2. Structured Meteorological State Architecture

We consolidate meteorological tracking under a unified `catchem::MetState` container, located inside `src/core/catchem_state_manager.hpp` or as a standalone module `src/core/catchem_met_state.hpp`.

### 2.1 C++ MetState Class Definition

Each meteorological field is represented as a `std::shared_ptr` to a strongly typed `InteropField`:

```cpp
// src/core/catchem_met_state.hpp
#pragma once
#include <memory>
#include <string>
#include <unordered_map>
#include "catchem_interop_field.hpp"

namespace catchem {

struct MetState {
    // Standard Grid flag fields (2D)
    std::shared_ptr<InteropField<double, 2>> IsLand;
    std::shared_ptr<InteropField<double, 2>> IsWater;
    std::shared_ptr<InteropField<double, 2>> IsIce;
    std::shared_ptr<InteropField<double, 2>> IsSnow;

    // 3D Volumetric fields
    std::shared_ptr<InteropField<double, 3>> T;          // Temperature [K]
    std::shared_ptr<InteropField<double, 3>> QV;         // Specific humidity [kg/kg]
    std::shared_ptr<InteropField<double, 3>> RH;         // Relative humidity [0-1]
    std::shared_ptr<InteropField<double, 3>> PMID;       // Mid-level pressure [Pa]
    std::shared_ptr<InteropField<double, 3>> PEDGE;      // Edge-level pressure [Pa]
    std::shared_ptr<InteropField<double, 3>> AIRDEN;     // Wet air density [kg/m³]
    std::shared_ptr<InteropField<double, 3>> AIRDEN_DRY; // Dry air density [kg/m³]
    std::shared_ptr<InteropField<double, 3>> BXHEIGHT;   // Layer thickness height [m]

    // 2D Surface fields
    std::shared_ptr<InteropField<double, 2>> PS;         // Surface pressure [Pa]
    std::shared_ptr<InteropField<double, 2>> TS;         // Surface temperature [K]
    std::shared_ptr<InteropField<double, 2>> PBLH;       // Boundary layer height [m]
    std::shared_ptr<InteropField<double, 2>> USTAR;      // Friction velocity [m/s]
    std::shared_ptr<InteropField<double, 2>> HFLUX;      // Sensible heat flux [W/m²]
    std::shared_ptr<InteropField<double, 2>> OBK;        // Monin-Obukhov length [m]
    std::shared_ptr<InteropField<double, 2>> LAT;        // Latitude [deg]
    std::shared_ptr<InteropField<double, 2>> LON;        // Longitude [deg]

    // Helper map for dynamic string-based binding
    std::unordered_map<std::string, std::shared_ptr<InteropField<double, 3>>> fields_3d;
    std::unordered_map<std::string, std::shared_ptr<InteropField<double, 2>>> fields_2d;

    void register_fields() {
        fields_3d["T"] = T;
        fields_3d["QV"] = QV;
        fields_3d["RH"] = RH;
        fields_3d["PMID"] = PMID;
        fields_3d["PEDGE"] = PEDGE;
        fields_3d["AIRDEN"] = AIRDEN;
        fields_3d["AIRDEN_DRY"] = AIRDEN_DRY;
        fields_3d["BXHEIGHT"] = BXHEIGHT;

        fields_2d["PS"] = PS;
        fields_2d["TS"] = TS;
        fields_2d["PBLH"] = PBLH;
        fields_2d["USTAR"] = USTAR;
        fields_2d["HFLUX"] = HFLUX;
        fields_2d["OBK"] = OBK;
        fields_2d["LAT"] = LAT;
        fields_2d["LON"] = LON;
    }
};

} // namespace catchem
```

---

### 2.2 Unmanaged Binding Mechanics (Zero-Copy)

For coupled simulations, raw pointer binders at the C-API boundary construct unmanaged View maps instantly:

```cpp
// Within src/core/catchem_api.cpp
void catchem_state_bind_met_3d(void* state_ptr, const char* name, double* ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    std::string key(name);

    // Bind directly to unmanaged InteropField mapping the Fortran heap address
    auto field = std::make_shared<catchem::InteropField<double, 3>>(ptr, std::vector<int>{state->n_cols, state->n_levels, 1});
    state->met.fields_3d[key] = field;

    // Maintain strong pointer bounds
    if (key == "T") state->met.T = field;
    else if (key == "QV") state->met.QV = field;
    else if (key == "RH") state->met.RH = field;
    else if (key == "PMID") state->met.PMID = field;
    else if (key == "PEDGE") state->met.PEDGE = field;
    else if (key == "AIRDEN") state->met.AIRDEN = field;
    else if (key == "AIRDEN_DRY") state->met.AIRDEN_DRY = field;
    else if (key == "BXHEIGHT") state->met.BXHEIGHT = field;
}
```

---

## 3. Unified Contiguous Chemistry State

To achieve uniform, GPU-coalesced access over memory loads during multi-phase chemical updates, we consolidate all species concentration buffers into a single, unified 3D View managed by `catchem::StateManager`.

```cpp
// src/core/catchem_state_manager.hpp additions
namespace catchem {

class StateManager {
public:
    // ... metadata ...
    int n_cols;
    int n_levels;
    int n_species;

    // Single unified 3D View: layout is (cols, levels, species)
    // Matches Fortran 4D LayoutLeft representation: (nx, ny, nz, n_species)
    std::shared_ptr<InteropField<double, 3>> unified_chem_state;

    // Standard allocation interface for standalone executions
    void allocate_unified_chemistry() {
        std::vector<int> dims = {n_cols, n_levels, n_species};
        unified_chem_state = std::make_shared<InteropField<double, 3>>(nullptr, dims);
    }

    // Direct C-API binder for zero-copy coupling
    void bind_unified_chemistry(double* fortran_contiguous_ptr) {
        std::vector<int> dims = {n_cols, n_levels, n_species};
        unified_chem_state = std::make_shared<InteropField<double, 3>>(fortran_contiguous_ptr, dims);
    }
};

} // namespace catchem
```

---

## 4. Performance Portable Time & Solar Zenith Engine

We replace `TimeState_Mod.F90` with a pure C++ struct containing mathematical formulas decorated with Kokkos directives. This enables high-performance solar tracking computations inside parallel compute grids.

```cpp
// src/core/catchem_time_state.hpp
#pragma once
#include <cmath>

#ifdef ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#define KOKKOS_FUNCTION KOKKOS_INLINE_FUNCTION
#else
#define KOKKOS_INLINE_FUNCTION inline
#define KOKKOS_FUNCTION inline
#endif

#include "catchem_precision.hpp"
#include "catchem_constants.hpp"

namespace catchem {

struct TimeState {
    int year = 2000;
    int month = 1;
    int day = 1;
    int hour = 0;
    int minute = 0;
    int second = 0;
    double timestep = 3600.0;
    double julian_date = 0.0;
    int doy = 1;

    KOKKOS_FUNCTION
    double get_cos_sza(double lat_deg, double lon_deg, bool mid_timestep = false) const {
        double lat_rad = lat_deg * constants::PI_180;
        double lon_rad = lon_deg * constants::PI_180;

        double frac_hour = hour + minute / 60.0 + second / 3600.0;
        if (mid_timestep) {
            frac_hour += (timestep / 2.0) / 3600.0;
        }

        // Day angle [radians]
        double gamma = 2.0 * constants::PI * (doy - 1.0) / 365.0;

        // Solar declination (high-precision Fourier calculation matches GOCART2G)
        double dec = 0.006918 - 0.399912 * std::cos(gamma) + 0.070257 * std::sin(gamma)
                     - 0.006758 * std::cos(2.0 * gamma) + 0.000907 * std::sin(2.0 * gamma)
                     - 0.002697 * std::cos(3.0 * gamma) + 0.001480 * std::sin(3.0 * gamma);

        // Equation of time
        double eqtime = 229.18 * (0.000075 + 0.001868 * std::cos(gamma) - 0.032077 * std::sin(gamma)
                                  - 0.014615 * std::cos(2.0 * gamma) - 0.040849 * std::sin(2.0 * gamma));

        double time_offset = eqtime + 4.0 * lon_deg;
        double true_solar_time = frac_hour * 60.0 + time_offset;

        double hour_angle = (true_solar_time / 4.0) - 180.0;
        double ha_rad = hour_angle * constants::PI_180;

        double cos_sza = std::sin(lat_rad) * std::sin(dec) + std::cos(lat_rad) * std::cos(dec) * std::cos(ha_rad);

        // Clamp output safely to [-1.0, 1.0]
        return std::max(-1.0, std::min(1.0, cos_sza));
    }
};

} // namespace catchem
```

---

## 5. Acceleration-Targeted Meteorological Derivations

Deriving meteorological variables (such as air density `AIRDEN_DRY` and layer box height `BXHEIGHT` based on edge pressures) is ported into dynamic C++ Kokkos loops.

### Ported Example: Box Height Calculation (Hydrostatic Equation)
Instead of CPU-bound double loops inside `metstate_mod.F90`, we execute the derivation inside Kokkos parallel region sweeps:

```cpp
// Inside src/core/catchem_state_manager.cpp
void StateManager::derive_bxheight() {
    int nc = n_cols;
    int nl = n_levels;

    // 1. Fetch Views
    auto t_view = met.T->view();
    auto qv_view = met.QV->view();
    auto pedge_view = met.PEDGE->view();
    auto bxheight_view = met.BXHEIGHT->view();

    double Rdg0 = constants::RD / constants::G0;

    // 2. Parallel Hydrostatic Integration
    Kokkos::parallel_for("derive_bxheight",
        Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, nc),
        KOKKOS_LAMBDA(int icol) {
            for (int k = 0; k < nl; ++k) {
                // Calculate virtual temperature: Tv = T * (1 + 0.608 * qv)
                double tv = t_view(icol, k, 0) * (1.0 + 0.608 * qv_view(icol, k, 0));

                // Box height from edge pressures: H_k = (R_d / g) * Tv * log(P_lower_edge / P_upper_edge)
                double p_lower = pedge_view(icol, k, 0);
                double p_upper = pedge_view(icol, k + 1, 0);

                bxheight_view(icol, k, 0) = Rdg0 * tv * std::log(p_lower / p_upper);
            }
        }
    );
}
```

---

## 6. Integration & Testing Blueprint

We will append verification assertions inside `tests/test_catchem_interop.cpp` to validate state transformations:
1. **Mock Fortran Arrays:** Allocate 3D contiguous heaps representing temperatures, moisture levels, pressures, and concentrations.
2. **C-API Bind Assertions:** Invoke C-API binders (`catchem_state_bind_met_3d` and `catchem_state_bind_unified_chemistry`).
3. **Execution space Verification:** Trigger `sync_to_device()` and execute a mock derived kernel. Update assertions after fetching data back to host space using `sync_to_host()`.
4. **Julian Calendar Assertions:** Validate timezone calculations and high-precision SZA values computed in the C++ time engine.
