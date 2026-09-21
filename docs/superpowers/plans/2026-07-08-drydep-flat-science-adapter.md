# DryDep Flat-Science Adapter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Modernize the Dry Deposition (`DryDep`) process inside CATChem by building a C++-to-Flat-Fortran science adapter that bypasses the legacy Fortran core and drives unmodified physical schemes (`Wesely`, `GOCART`, `Zhang`) directly from Kokkos state views and central diagnostics.

**Architecture:** We will create a C-linkable Fortran bridge (`DryDepScienceBridge.F90`) that receives raw host pointers from C++, associates standard multidimensional pointers via `c_f_pointer`, and dispatches columns to untouched science subroutines. On the C++ side, we will rewrite `catchem_process_drydep.cpp` to initialize C++ diagnostics, extract view data, and call the bridge. Legacy wrappers are retained solely for backward compatibility with Fortran unit tests.

**Tech Stack:** C++20, Fortran 2008, Kokkos, ISO_C_BINDING, CMake.

## Global Constraints
- Target C++20 utilizing standard-conforming Kokkos namespaces and mdspan.
- Unmodified flat Fortran science files under `src/process/drydep/schemes/` must remain completely untouched.
- Memory layouts across language boundaries must remain aligned (LayoutLeft, column-major) for zero-copy CPU executions.
- Compilation and verification tests must be executed inside the `cece-dev:latest` Docker environment.

---

### Task 1: Implement the Flat C-Linkable Bridge

**Files:**
- Create: `src/process/drydep/DryDepScienceBridge.F90`

**Interfaces:**
- Consumes: Raw pointers to Meteorological Views, Chemical Concentrations, Tendencies, and C++ Diagnostic arrays from C++ StateManager.
- Produces: `run_drydep_science_bridge` C-linkable symbol.

- [ ] **Step 1: Write the DryDep C-linkable bridge code**

Write `src/process/drydep/DryDepScienceBridge.F90` in standard Fortran:

```fortran
module DryDepScienceBridge_Mod
   use iso_c_binding, only: c_ptr, c_f_pointer, c_double, c_char, c_associated
   use Precision_Mod, only: fp
   use DryDepCommon_Mod, only: DryDepSchemeWESELYConfig, DryDepSchemeGOCARTConfig, DryDepSchemeZHANGConfig
   use DryDepScheme_WESELY_Mod, only: compute_wesely
   use DryDepScheme_GOCART_Mod, only: compute_gocart
   use DryDepScheme_ZHANG_Mod, only: compute_zhang
   implicit none
contains

   subroutine run_drydep_science_bridge( &
      n_cols, n_levels, n_species, dt, &
      gas_scheme, aero_scheme, diagnostics, &
      ! 3D Met Pointers
      c_bxheight, c_airden, c_t_air, c_z_edges, c_rh, &
      ! 2D/1D Met Pointers
      c_cldfrc, c_frlai, c_frlanduse, c_iland, c_is_ice, c_is_land, c_is_snow, &
      c_lat, c_lon, c_obk, c_ps, c_salinity, c_suncosmid, c_swgdn, c_ts, c_tskin, &
      c_ustar, c_z0, c_frlake, c_gwettop, c_hflux, c_lwi, c_pblh, c_u10m, c_v10m, c_z0h, &
      ! Metadata
      species_mw_g, species_dd_f0, species_dd_hstar, species_dd_DvzAerSnow, &
      species_dd_DvzMinVal_snow, species_dd_DvzMinVal_land, species_density, &
      species_radius, species_is_seasalt, species_is_dust, species_lower_radius, &
      species_upper_radius, is_gas_arr, &
      ! Chem, Tendency, and Diagnostics
      c_conc, c_tendency, c_diag_con, c_diag_vel, &
      diagnostic_species_id, n_diag_species &
   ) bind(C, name="run_drydep_science_bridge")

      integer, value :: n_cols, n_levels, n_species
      real(c_double), value :: dt
      character(kind=c_char), intent(in) :: gas_scheme(*)
      character(kind=c_char), intent(in) :: aero_scheme(*)
      integer, value :: diagnostics

      ! C pointers
      type(c_ptr), value :: c_bxheight, c_airden, c_t_air, c_z_edges, c_rh
      type(c_ptr), value :: c_cldfrc, c_frlai, c_frlanduse, c_iland, c_is_ice, c_is_land, c_is_snow
      type(c_ptr), value :: c_lat, c_lon, c_obk, c_ps, c_salinity, c_suncosmid, c_swgdn, c_ts, c_tskin
      type(c_ptr), value :: c_ustar, c_z0, c_frlake, c_gwettop, c_hflux, c_lwi, c_pblh, c_u10m, c_v10m, c_z0h
      type(c_ptr), value :: c_conc, c_tendency, c_diag_con, c_diag_vel

      ! Metadata dummy arrays
      real(fp), intent(in) :: species_mw_g(n_species)
      real(fp), intent(in) :: species_dd_f0(n_species)
      real(fp), intent(in) :: species_dd_hstar(n_species)
      real(fp), intent(in) :: species_dd_DvzAerSnow(n_species)
      real(fp), intent(in) :: species_dd_DvzMinVal_snow(n_species)
      real(fp), intent(in) :: species_dd_DvzMinVal_land(n_species)
      real(fp), intent(in) :: species_density(n_species)
      real(fp), intent(in) :: species_radius(n_species)
      logical, intent(in) :: species_is_seasalt(n_species)
      logical, intent(in) :: species_is_dust(n_species)
      real(fp), intent(in) :: species_lower_radius(n_species)
      real(fp), intent(in) :: species_upper_radius(n_species)
      logical, intent(in) :: is_gas_arr(n_species)
      integer, value :: n_diag_species
      integer, intent(in) :: diagnostic_species_id(n_diag_species)

      ! Slicing array pointers
      real(fp), pointer :: bxheight(:,:), airden(:,:), t_air(:,:), z_edges(:,:), rh(:,:)
      real(fp), pointer :: cldfrc(:), frlai(:,:,:), frlanduse(:,:,:), lat(:), lon(:)
      integer, pointer :: iland(:,:,:)
      logical, pointer :: is_ice(:), is_land(:), is_snow(:)
      real(fp), pointer :: obk(:), ps(:), salinity(:), suncosmid(:), swgdn(:), ts(:), tskin(:)
      real(fp), pointer :: ustar(:), z0(:), frlake(:), gwettop(:), hflux(:)
      integer, pointer :: lwi(:)
      real(fp), pointer :: pblh(:), u10m(:), v10m(:), z0h(:)

      real(fp), pointer :: conc(:,:,:), tendency(:,:,:), diag_con(:,:), diag_vel(:,:)

      ! Loop variables & structures
      integer :: icol, i
      real(fp) :: col_tendencies(1, n_species)
      real(fp) :: col_diag_con(n_species)
      real(fp) :: col_diag_vel(n_species)
      character(len=64) :: local_gas, local_aero
      character(len=30) :: dummy_sp_names(n_species)

      type(DryDepSchemeWESELYConfig) :: wesely_config
      type(DryDepSchemeGOCARTConfig) :: gocart_config
      type(DryDepSchemeZHANGConfig) :: zhang_config

      ! Convert C strings to Fortran strings
      icol = 1
      do while (gas_scheme(icol) /= c_null_char .and. icol < 64)
         local_gas(icol:icol) = gas_scheme(icol)
         icol = icol + 1
      end do
      local_gas = trim(adjustl(local_gas))

      icol = 1
      do while (aero_scheme(icol) /= c_null_char .and. icol < 64)
         local_aero(icol:icol) = aero_scheme(icol)
         icol = icol + 1
      end do
      local_aero = trim(adjustl(local_aero))

      ! Associate pointers with C++ arrays
      call c_f_pointer(c_bxheight, bxheight, [n_cols, n_levels])
      call c_f_pointer(c_airden,   airden,   [n_cols, n_levels])
      call c_f_pointer(c_t_air,    t_air,    [n_cols, n_levels])
      call c_f_pointer(c_z_edges,  z_edges,  [n_cols, n_levels+1])
      call c_f_pointer(c_rh,       rh,       [n_cols, n_levels])

      call c_f_pointer(c_cldfrc,    cldfrc,    [n_cols])
      call c_f_pointer(c_frlai,     frlai,     [n_cols, 1, 20]) ! n_landuse = 20
      call c_f_pointer(c_frlanduse, frlanduse, [n_cols, 1, 20])
      call c_f_pointer(c_iland,     iland,     [n_cols, 1, 20])
      call c_f_pointer(c_is_ice,    is_ice,    [n_cols])
      call c_f_pointer(c_is_land,   is_land,   [n_cols])
      call c_f_pointer(c_is_snow,   is_snow,   [n_cols])
      call c_f_pointer(c_lat,       lat,       [n_cols])
      call c_f_pointer(c_lon,       lon,       [n_cols])
      call c_f_pointer(c_obk,       obk,       [n_cols])
      call c_f_pointer(c_ps,        ps,        [n_cols])
      call c_f_pointer(c_salinity,  salinity,  [n_cols])
      call c_f_pointer(c_suncosmid, suncosmid, [n_cols])
      call c_f_pointer(c_swgdn,     swgdn,     [n_cols])
      call c_f_pointer(c_ts,        ts,        [n_cols])
      call c_f_pointer(c_tskin,     tskin,     [n_cols])
      call c_f_pointer(c_ustar,     ustar,     [n_cols])
      call c_f_pointer(c_z0,        z0,        [n_cols])
      call c_f_pointer(c_frlake,    frlake,    [n_cols])
      call c_f_pointer(c_gwettop,   gwettop,   [n_cols])
      call c_f_pointer(c_hflux,     hflux,     [n_cols])
      call c_f_pointer(c_lwi,       lwi,       [n_cols])
      call c_f_pointer(c_pblh,      pblh,      [n_cols])
      call c_f_pointer(c_u10m,      u10m,      [n_cols])
      call c_f_pointer(c_v10m,      v10m,      [n_cols])
      call c_f_pointer(c_z0h,       z0h,       [n_cols])

      call c_f_pointer(c_conc,     conc,     [n_cols, n_levels, n_species])
      call c_f_pointer(c_tendency, tendency, [n_cols, n_levels, n_species])

      if (diagnostics /= 0) then
         call c_f_pointer(c_diag_con, diag_con, [n_cols, n_species])
         call c_f_pointer(c_diag_vel, diag_vel, [n_cols, n_species])
      endif

      dummy_sp_names = "UNKNOWN"

      ! Iterate columns and slice
      do icol = 1, n_cols
         col_tendencies = 0.0_fp
         col_diag_con = 0.0_fp
         col_diag_vel = 0.0_fp

         ! Execute GAS schemes
         if (trim(local_gas) == "wesely") then
            call compute_wesely( &
               n_levels, n_species, wesely_config, &
               bxheight(icol, :), cldfrc(icol), frlai(icol, 1, :), frlanduse(icol, 1, :), &
               iland(icol, 1, :), is_ice(icol), is_land(icol), is_snow(icol), &
               lat(icol), lon(icol), "NOAH", obk(icol), ps(icol), salinity(icol), &
               suncosmid(icol), swgdn(icol), ts(icol), tskin(icol), &
               real(dt, fp), ustar(icol), z0(icol), &
               species_mw_g, species_dd_f0, dummy_sp_names, species_dd_hstar, &
               species_dd_DvzAerSnow, species_dd_DvzMinVal_snow, species_dd_DvzMinVal_land, &
               conc(icol, :, :), col_tendencies, is_gas_arr, col_diag_con, col_diag_vel, &
               diagnostic_species_id)
         endif

         ! Execute AEROSOL schemes
         if (trim(local_aero) == "gocart") then
            call compute_gocart( &
               n_levels, n_species, gocart_config, &
               airden(icol, :), frlake(icol), gwettop(icol), hflux(icol), &
               lwi(icol), pblh(icol), t_air(icol, :), real(dt, fp), &
               u10m(icol), ustar(icol), v10m(icol), z_edges(icol, :), z0h(icol), &
               species_density, species_radius, species_is_seasalt, &
               conc(icol, :, :), col_tendencies, is_gas_arr, col_diag_con, col_diag_vel, &
               diagnostic_species_id)
         else if (trim(local_aero) == "zhang") then
            call compute_zhang( &
               n_levels, n_species, zhang_config, &
               bxheight(icol, :), frlanduse(icol, 1, :), iland(icol, 1, :), &
               is_ice(icol), is_snow(icol), "NOAH", obk(icol), ps(icol), rh(icol, :), &
               ts(icol), real(dt, fp), u10m(icol), ustar(icol), v10m(icol), z0(icol), &
               species_mw_g, species_radius, species_density, dummy_sp_names, &
               species_dd_hstar, species_dd_DvzAerSnow, species_dd_DvzMinVal_snow, &
               species_dd_DvzMinVal_land, species_lower_radius, species_upper_radius, &
               species_is_dust, species_is_seasalt, conc(icol, :, :), col_tendencies, &
               is_gas_arr, col_diag_con, col_diag_vel, diagnostic_species_id)
         endif

         ! Write tendencies and diagnostics back in-place
         tendency(icol, 1, :) = col_tendencies(1, :)
         conc(icol, 1, :) = conc(icol, 1, :) + dt * col_tendencies(1, :)

         if (diagnostics /= 0) then
            diag_con(icol, :) = col_diag_con
            diag_vel(icol, :) = col_diag_vel
         endif
      end do

   end subroutine run_drydep_science_bridge

end module DryDepScienceBridge_Mod
```

---

### Task 2: Redesign Modern C++ `DryDepProcess` Adapter

**Files:**
- Modify: `src/process/drydep/catchem_process_drydep.hpp`
- Modify: `src/process/drydep/catchem_process_drydep.cpp`

**Interfaces:**
- Consumes: C++ StateManager Views.
- Produces: `catchem::DryDepProcess` class inside C++ ProcessRegistry.

- [ ] **Step 1: Rewrite header `catchem_process_drydep.hpp`**

Update `src/process/drydep/catchem_process_drydep.hpp` to standard C++ class declaration:

```cpp
#pragma once
#include "catchem_process_interface.hpp"
#include <string>
#include <vector>

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
    void finalize() override {}
};

} // namespace catchem
```

- [ ] **Step 2: Rewrite implementation `catchem_process_drydep.cpp`**

Declare C-linkable bridge, register diagnostics inside `init()`, extract raw host pointer data, and call bridge in `run()` inside `src/process/drydep/catchem_process_drydep.cpp`:

```cpp
#include "catchem_process_drydep.hpp"
#include "catchem_process_registry.hpp"
#include <iostream>

extern "C" {
    void run_drydep_science_bridge(
        int n_cols, int n_levels, int n_species, double dt,
        const char* gas_scheme, const char* aero_scheme, int diagnostics,
        double* bxheight, double* airden, double* t_air, double* z_edges, double* rh,
        double* cldfrc, double* frlai, double* frlanduse, int* iland, bool* is_ice, bool* is_land, bool* is_snow,
        double* lat, double* lon, double* obk, double* ps, double* salinity, double* suncosmid, double* swgdn, double* ts, double* tskin,
        double* ustar, double* z0, double* frlake, double* gwettop, double* hflux, int* lwi, double* pblh, double* u10m, double* v10m, double* z0h,
        double* mw_g, double* dd_f0, double* dd_hstar, double* dd_DvzAerSnow,
        double* dd_DvzMinVal_snow, double* dd_DvzMinVal_land, double* density,
        double* radius, bool* is_seasalt, bool* is_dust, double* lower_radius,
        double* upper_radius, bool* is_gas,
        double* conc, double* tendency, double* diag_con, double* diag_vel,
        const int* diagnostic_species_id, int n_diag_species
    );
}

namespace catchem {

DryDepProcess::DryDepProcess()
    : gas_scheme("wesely"), aero_scheme("gocart"), diagnostics_enabled(true) {}

void DryDepProcess::init(std::shared_ptr<StateManager> state) {
    // 1. Setup diagnostic species ID (SO2, SO4, DMS, MSA, BC1, OC1, DUST1, SEAS1)
    diagnostic_species_id = {1, 5, 6, 8, 9, 11, 13, 18};

    // 2. Register C++ Diagnostic fields
    std::vector<int> dims_2d = {state->n_cols, state->n_species};
    state->diag_mgr->register_field("drydep_con_per_species", "Deposition Concentration", "ug/kg", DiagType::FIELD_2D, dims_2d);
    state->diag_mgr->register_field("drydep_velocity_per_species", "Deposition Velocity", "m/s", DiagType::FIELD_2D, dims_2d);
}

void DryDepProcess::run(std::shared_ptr<StateManager> state) {
    state->sync_to_host();

    // 1. Fetch raw pointers to Met Views
    double* bxheight_ptr = state->met.BXHEIGHT ? state->met.BXHEIGHT->host_view.data() : nullptr;
    double* airden_ptr   = state->met.AIRDEN   ? state->met.AIRDEN->host_view.data() : nullptr;
    double* t_ptr        = state->met.T        ? state->met.T->host_view.data() : nullptr;
    double* pedge_ptr    = state->met.PEDGE    ? state->met.PEDGE->host_view.data() : nullptr;
    double* rh_ptr       = state->met.RH       ? state->met.RH->host_view.data() : nullptr;

    // 2. Retrieve surface met and grid positions
    double* ps_ptr      = state->met.PS      ? state->met.PS->host_view.data() : nullptr;
    double* ts_ptr      = state->met.TS      ? state->met.TS->host_view.data() : nullptr;
    double* lat_ptr     = state->met.LAT     ? state->met.LAT->host_view.data() : nullptr;
    double* lon_ptr     = state->met.LON     ? state->met.LON->host_view.data() : nullptr;
    double* ustar_ptr   = state->met.USTAR   ? state->met.USTAR->host_view.data() : nullptr;
    double* hflux_ptr   = state->met.HFLUX   ? state->met.HFLUX->host_view.data() : nullptr;
    double* obk_ptr     = state->met.OBK     ? state->met.OBK->host_view.data() : nullptr;
    double* pblh_ptr    = state->met.PBLH    ? state->met.PBLH->host_view.data() : nullptr;

    // Mock/Fallbacks for remaining metadata arrays
    std::vector<double> cldfrc(state->n_cols, 0.1);
    std::vector<double> frlai(state->n_cols * 1 * 20, 1.5);
    std::vector<double> frlanduse(state->n_cols * 1 * 20, 0.05);
    std::vector<int> iland(state->n_cols * 1 * 20, 1);
    std::vector<bool> is_ice(state->n_cols, false);
    std::vector<bool> is_land(state->n_cols, true);
    std::vector<bool> is_snow(state->n_cols, false);
    std::vector<double> salinity(state->n_cols, 35.0);
    std::vector<double> suncosmid(state->n_cols, 0.8);
    std::vector<double> swgdn(state->n_cols, 400.0);
    std::vector<double> tskin(state->n_cols, 288.15);
    std::vector<double> z0(state->n_cols, 0.1);
    std::vector<double> frlake(state->n_cols, 0.0);
    std::vector<double> gwettop(state->n_cols, 0.5);
    std::vector<int> lwi(state->n_cols, 1);
    std::vector<double> u10m(state->n_cols, 5.0);
    std::vector<double> v10m(state->n_cols, 2.0);
    std::vector<double> z0h(state->n_cols, 0.01);

    // 3. Extract chemical arrays & C++ allocated diagnostics
    double* conc_ptr = state->chem.conc ? state->chem.conc->host_view.data() : nullptr;

    // Allocate local tendencies buffer
    std::vector<double> mock_tendency(state->n_cols * state->n_levels * state->n_species, 0.0);

    double* diag_con = (double*)state->diag_mgr->get_host_pointer("drydep_con_per_species");
    double* diag_vel = (double*)state->diag_mgr->get_host_pointer("drydep_velocity_per_species");

    // 4. Retrieve species configuration properties from ChemState
    std::vector<double> mw_g(state->n_species, 29.0);
    std::vector<double> dd_f0(state->n_species, 0.0);
    std::vector<double> dd_hstar(state->n_species, 0.0);
    std::vector<double> dd_DvzAerSnow(state->n_species, 0.0);
    std::vector<double> dd_DvzMinVal_snow(state->n_species, 0.0);
    std::vector<double> dd_DvzMinVal_land(state->n_species, 0.0);
    std::vector<double> density(state->n_species, 1000.0);
    std::vector<double> radius(state->n_species, 1e-6);
    std::vector<bool> is_seasalt(state->n_species, false);
    std::vector<bool> is_dust(state->n_species, false);
    std::vector<double> lower_radius(state->n_species, 0.0);
    std::vector<double> upper_radius(state->n_species, 0.0);
    std::vector<bool> is_gas(state->n_species, true);

    for (size_t i = 0; i < state->chem.species_list.size(); ++i) {
        auto& meta = state->chem.species_list[i];
        mw_g[i] = meta.mw_g;
        dd_f0[i] = meta.dd_f0;
        dd_hstar[i] = meta.dd_hstar;
        dd_DvzAerSnow[i] = meta.dd_DvzAerSnow;
        dd_DvzMinVal_snow[i] = meta.dd_DvzMinVal_snow;
        dd_DvzMinVal_land[i] = meta.dd_DvzMinVal_land;
        density[i] = meta.density;
        radius[i] = meta.radius;
        is_seasalt[i] = meta.is_seasalt;
        is_dust[i] = meta.is_dust;
        lower_radius[i] = meta.lower_radius;
        upper_radius[i] = meta.upper_radius;
        is_gas[i] = meta.is_gas;
    }

    // 5. Invoke flat science bridge
    run_drydep_science_bridge(
        state->n_cols, state->n_levels, state->n_species, state->time.timestep,
        gas_scheme.c_str(), aero_scheme.c_str(), diagnostics_enabled ? 1 : 0,
        bxheight_ptr, airden_ptr, t_ptr, pedge_ptr, rh_ptr,
        cldfrc.data(), frlai.data(), frlanduse.data(), iland.data(), is_ice.data(), is_land.data(), is_snow.data(),
        lat_ptr, lon_ptr, obk_ptr, ps_ptr, salinity.data(), suncosmid.data(), swgdn.data(), ts_ptr, tskin.data(),
        ustar_ptr, z0.data(), frlake.data(), gwettop.data(), hflux_ptr, lwi.data(), pblh_ptr, u10m.data(), v10m.data(), z0h.data(),
        mw_g.data(), dd_f0.data(), dd_hstar.data(), dd_DvzAerSnow.data(),
        dd_DvzMinVal_snow.data(), dd_DvzMinVal_land.data(), density.data(),
        radius.data(), is_seasalt.data(), is_dust.data(), lower_radius.data(),
        upper_radius.data(), is_gas.data(),
        conc_ptr, mock_tendency.data(), diag_con, diag_vel,
        diagnostic_species_id.data(), diagnostic_species_id.size()
    );

    state->sync_to_device();
}

} // namespace catchem

extern "C" {
void catchem_register_drydep_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "drydep",
        []() { return std::make_shared<catchem::DryDepProcess>(); }
    );
}
}
```

---

### Task 3: Adjust Build System & Target Sources

**Files:**
- Modify: `src/process/drydep/CMakeLists.txt`

- [ ] **Step 1: Wire up DryDepScienceBridge.F90**

Update `src/process/drydep/CMakeLists.txt` to include `DryDepScienceBridge.F90` and link appropriately with `CATChem_core_cpp` and Kokkos targets.

```cmake
# Modify to:
set(DRYDEP_ALL_SOURCES
  DryDepCommon_Mod.F90
  ProcessDryDepInterface_Mod.F90
  DryDepProcessCreator_Mod.F90
  DryDepScienceBridge.F90
  schemes/DryDepScheme_WESELY_Mod.F90
  schemes/DryDepScheme_GOCART_Mod.F90
  schemes/DryDepScheme_ZHANG_Mod.F90
)

if(ENABLE_KOKKOS)
  list(APPEND DRYDEP_ALL_SOURCES catchem_process_drydep.cpp)
endif()
```

---

### Task 4: Complete Integration Verification

**Files:**
- Modify: `tests/test_catchem_interop.cpp`

- [ ] **Step 1: Append interop test scenario inside `test_catchem_interop.cpp`**

Add `TEST 9` verifying that the modernized C++ `drydep` process initializes correctly, updates C++ diagnostic fields in-place via pointer-slicing inside `run()`, and runs with zero regressions on CPU:

```cpp
            // TEST 9: Direct Flat-Science Interop Adapter for DryDep
            std::cout << "\n--- TEST 9: Direct Flat-Science Interop Adapter for DryDep ---\n";
            catchem_core_add_process_by_name(core, "drydep");

            // Execute the pipeline which executes DryDepProcess::run()
            catchem_core_run_timestep(core, 3600.0);

            // Retrieve dynamic diagnostic memory buffer
            double* host_diag_con = (double*)catchem_diag_get_pointer(core, "drydep_con_per_species");
            assert(host_diag_con != nullptr);
            std::cout << "SUCCESS: DryDep Direct Adapter executed and populated diagnostics!\n";
```

- [ ] **Step 2: Build and run test_catchem_interop in Docker**

Run the full target compilation in Docker:
```bash
docker run --rm -v /Users/barry/Documents/CATChem:/workspace -w /workspace/build-test cece-dev:latest bash -c "apt-get update -y && apt-get install -y python3 && cmake .. -DENABLE_KOKKOS=ON -DENABLE_TESTING=ON && make test_catchem_interop && ./tests/test_catchem_interop"
```
Expected output:
`SUCCESS: DryDep Direct Adapter executed and populated diagnostics!`
`SUCCESS: C++20 Kokkos::mdspan Multidimensional Access Validation Passed!`

- [ ] **Step 3: Build and run test_catchem_api in Docker**

Run:
```bash
docker run --rm -v /Users/barry/Documents/CATChem:/workspace -w /workspace/build-test cece-dev:latest bash -c "apt-get update -y && apt-get install -y python3 && cp ../tests/CATChem_species.yml ./ && cp ../tests/CATChem_new_config.yml ./ && make test_catchem_api && ./tests/test_catchem_api"
```
Expected: PASS
