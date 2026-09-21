# Specification: Species Database and Meteorological Physics Migration to C++

* **Status:** Approved
* **Authors:** Gemini CLI Architect
* **Created:** July 9, 2026
* **Target Version:** 2.1.0
* **Pillars:** Complete C++ Core Ownership, Zero-Overhead Metadata Sync, Thread-Safe Physics Calculations

## 1. Executive Summary & Architecture

To fully eliminate legacy calculations and custom file-I/O from Fortran, we are modernizing both **`species_mod.F90`** and **`met_utilities_mod.F90`**.

1. **`species_mod.F90`:** All species configurations (`CATChem_species.yml`) will be loaded, parsed, and validated solely in C++ using `yaml-cpp` inside `catchem_chem_state.hpp`. At initialization, the Fortran `SpeciesManagerType` queries the C++ core once via Flat BIND(C) APIs to sync all metadata into local Fortran array structs. This ensures **zero-overhead, local memory reads** during nested physical calculations inside standard schemes, while completely eliminating Fortran config parsing and memory allocators.
2. **`met_utilities_mod.F90`:** All high-frequency meteorological equations (e.g. Stokes settling velocity, Businger stability, Monin-Obukhov lengths) are relocated to template-inlined, thread-safe functions inside `catchem_met_utilities.hpp`, enabling full GPU/CPU parallelization in Kokkos. The Fortran module is reduced to a flat forwarding wrapper.

---

## 2. Dynamic Metadata Expansion

We are expanding the C++ `SpeciesMetadata` database to store the remaining dry/wet deposition, carbon loss, and background properties:

```cpp
namespace catchem {
    struct SpeciesMetadata {
        std::string short_name;
        std::string long_name;
        std::string description;

        bool is_gas = false;
        bool is_aerosol = false;
        bool is_tracer = false;
        bool is_advected = true;
        bool is_drydep = false;
        bool is_wetdep = false;
        bool is_photolysis = false;
        bool is_gocart_aero = false;
        bool is_dust = false;
        bool is_seasalt = false;

        double mw_g = 0.0;
        double density = 0.0;
        double radius = 0.0;
        double lower_radius = 0.0;
        double upper_radius = 0.0;
        double viscosity = 0.0;

        // Dry deposition parameters
        double dd_f0 = 0.0;
        double dd_hstar = 0.0;
        double dd_DvzAerSnow = 0.0;
        double dd_DvzMinVal_snow = 0.0;
        double dd_DvzMinVal_land = 0.0;

        // Wet deposition parameters
        double henry_k0 = 0.0;
        double henry_cr = 0.0;
        double henry_pKa = 0.0;
        double wd_retfactor = 0.0;
        bool wd_LiqAndGas = false;
        double wd_convfacI2G = 0.0;
        std::vector<double> wd_rainouteff = {0.0, 0.0, 0.0};
        double wd_reevap_frac = 0.5;

        // Chemical loss rate and background volume-mixing ratio
        double t_chem_loss = -1.0;
        double BackgroundVV = 1.0e-20;
        std::string mie_name;
    };
}
```

Exposing these fields to Fortran via `extern "C"` exports allows `species_mod.F90` to load and validate 100% of species properties dynamically.

---

## 3. Flat C Boundary APIs

The corresponding C-linkable interfaces are declared in `catchem_api.hpp` and protected inside `catchem_api.cpp` with exception shields:

```cpp
extern "C" {
    // =========================================================================
    // Species Database Query Exports
    // =========================================================================
    int catchem_state_get_species_count(void* state_ptr);
    void catchem_state_get_species_name_at(void* state_ptr, int index, char* name_out);
    void catchem_state_get_species_long_name_at(void* state_ptr, int index, char* name_out);
    void catchem_state_get_species_desc_at(void* state_ptr, int index, char* desc_out);

    double catchem_state_get_species_mw(void* state_ptr, int index);
    double catchem_state_get_species_density(void* state_ptr, int index);
    double catchem_state_get_species_radius(void* state_ptr, int index);
    double catchem_state_get_species_lower_radius(void* state_ptr, int index);
    double catchem_state_get_species_upper_radius(void* state_ptr, int index);
    double catchem_state_get_species_viscosity(void* state_ptr, int index);

    int catchem_state_is_species_gas(void* state_ptr, int index);
    int catchem_state_is_species_aerosol(void* state_ptr, int index);
    int catchem_state_get_species_is_tracer(void* state_ptr, int index);
    int catchem_state_get_species_is_advected(void* state_ptr, int index);
    int catchem_state_get_species_is_drydep(void* state_ptr, int index);
    int catchem_state_get_species_is_wetdep(void* state_ptr, int index);
    int catchem_state_get_species_is_photolysis(void* state_ptr, int index);
    int catchem_state_get_species_is_dust(void* state_ptr, int index);
    int catchem_state_get_species_is_seasalt(void* state_ptr, int index);

    double catchem_state_get_species_dd_f0(void* state_ptr, int index);
    double catchem_state_get_species_dd_hstar(void* state_ptr, int index);
    double catchem_state_get_species_dd_DvzAerSnow(void* state_ptr, int index);
    double catchem_state_get_species_dd_DvzMinVal_snow(void* state_ptr, int index);
    double catchem_state_get_species_dd_DvzMinVal_land(void* state_ptr, int index);

    double catchem_state_get_species_henry_k0(void* state_ptr, int index);
    double catchem_state_get_species_henry_cr(void* state_ptr, int index);
    double catchem_state_get_species_henry_pKa(void* state_ptr, int index);
    double catchem_state_get_species_wd_retfactor(void* state_ptr, int index);
    int catchem_state_get_species_wd_LiqAndGas(void* state_ptr, int index);
    double catchem_state_get_species_wd_convfacI2G(void* state_ptr, int index);
    void catchem_state_get_species_wd_rainouteff(void* state_ptr, int index, double* eff_out);
    double catchem_state_get_species_wd_reevap_frac(void* state_ptr, int index);
    double catchem_state_get_species_t_chem_loss(void* state_ptr, int index);
    double catchem_state_get_species_BackgroundVV(void* state_ptr, int index);
    void catchem_state_get_species_mie_name(void* state_ptr, int index, char* name_out);

    // =========================================================================
    // MetUtilities Core Calculation Exports
    // =========================================================================
    double catchem_met_potential_temperature(double temp, double press, double sfc_press);
    double catchem_met_virtual_temperature(double temp, double qv);
    double catchem_met_dew_point(double temp, double rh);
    double catchem_met_relative_humidity(double temp, double qv, double press);
    double catchem_met_saturation_vapor_pressure(double temp);
    double catchem_met_monin_obukhov_length(double ustar, double t0, double hflux, double rho);
    double catchem_met_friction_velocity(double tau, double rho);
    double catchem_met_cunningham_correction_factor(double dp, double lambda);
    double catchem_met_mean_free_path_air(double temp, double press);
    void catchem_met_solar_zenith_angle(int doy, double hour, double lat_rad, double lon_rad, double* sza_deg, double* cossza);
}
```

---

## 4. Verification and Downstream Compatibility

All compiled targets (static library `libCATChem_core.a`, unit tests `test_Precision`, `test_MetState`, `test_TimeState`, and `test_UnitConversion`) will be built and tested within the GCC Spack Container to confirm exact numerical consistency, case insensitivity, and backward compatibility.
