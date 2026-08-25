#include "catchem_process_drydep.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_process_registry.hpp"
#include <iostream>

extern "C" {
void run_drydep_science_bridge(int n_cols, int n_levels, int n_species, double dt, const char* gas_scheme,
                               const char* aero_scheme, int diagnostics, double* bxheight, double* airden,
                               double* t_air, double* pedge, double* rh, double* cldfrc, double* frlai,
                               double* frlanduse, int* iland, bool* is_ice, bool* is_land, bool* is_snow, double* lat,
                               double* lon, double* obk, double* ps, double* salinity, double* suncosmid, double* swgdn,
                               double* ts, double* tskin, double* ustar, double* z0, double* frlake, double* gwettop,
                               double* hflux, int* lwi, double* pblh, double* u10m, double* v10m, double* z0h,
                               double* mw_g, double* dd_f0, double* dd_hstar, double* dd_DvzAerSnow,
                               double* dd_DvzMinVal_snow, double* dd_DvzMinVal_land, double* density, double* radius,
                               bool* is_seasalt, bool* is_dust, double* lower_radius, double* upper_radius,
                               bool* is_gas, double* conc, double* tendency, double* diag_con, double* diag_vel,
                               const int* diagnostic_species_id, int n_diag_species);
}

namespace catchem {

    DryDepProcess::DryDepProcess() : gas_scheme("wesely"), aero_scheme("gocart"), diagnostics_enabled(true) {}

    void DryDepProcess::init(std::shared_ptr<StateManager> state) {
        // 1. Setup diagnostic species ID dynamically based on the is_drydep metadata switch
        for (size_t i = 0; i < state->chem.species_list.size(); ++i) {
            if (state->chem.species_list[i].is_drydep) {
                diagnostic_species_id.push_back(i + 1); // 1-based for Fortran bridge
            }
        }

        // 2. Register C++ Diagnostic fields
        std::vector<int> dims_2d = {state->n_cols, state->n_species};
        state->diag_mgr->register_field("drydep_con_per_species", "Deposition Concentration", "ug/kg",
                                        DiagType::FIELD_2D, dims_2d);
        state->diag_mgr->register_field("drydep_velocity_per_species", "Deposition Velocity", "m/s", DiagType::FIELD_2D,
                                        dims_2d);
    }

    void DryDepProcess::run(std::shared_ptr<StateManager> state) {
        state->sync_to_host();

        // 1. Fetch raw pointers to Met Views
        double* bxheight_ptr = state->met.BXHEIGHT ? state->met.BXHEIGHT->host_view.data() : nullptr;
        double* airden_ptr = state->met.AIRDEN ? state->met.AIRDEN->host_view.data() : nullptr;
        double* t_ptr = state->met.T ? state->met.T->host_view.data() : nullptr;
        double* pedge_ptr = state->met.PEDGE ? state->met.PEDGE->host_view.data() : nullptr;
        double* rh_ptr = state->met.RH ? state->met.RH->host_view.data() : nullptr;

        // 2. Retrieve surface met and grid positions
        double* ps_ptr = state->met.PS ? state->met.PS->host_view.data() : nullptr;
        double* ts_ptr = state->met.TS ? state->met.TS->host_view.data() : nullptr;
        double* lat_ptr = state->met.LAT ? state->met.LAT->host_view.data() : nullptr;
        double* lon_ptr = state->met.LON ? state->met.LON->host_view.data() : nullptr;
        double* ustar_ptr = state->met.USTAR ? state->met.USTAR->host_view.data() : nullptr;
        double* hflux_ptr = state->met.HFLUX ? state->met.HFLUX->host_view.data() : nullptr;
        double* obk_ptr = state->met.OBK ? state->met.OBK->host_view.data() : nullptr;
        double* pblh_ptr = state->met.PBLH ? state->met.PBLH->host_view.data() : nullptr;

        // Mock/Fallbacks for remaining metadata arrays - Using char for bool to support standard .data()
        std::vector<double> cldfrc(state->n_cols, 0.1);

        // Multi-dimensional standard C++20 views using standard layout_left to match Fortran column-major
        std::vector<double> frlai_storage(state->n_cols * 1 * 20, 1.5);
        std::vector<double> frlanduse_storage(state->n_cols * 1 * 20, 0.05);
        std::vector<int> iland_storage(state->n_cols * 1 * 20, 1);

        Kokkos::mdspan<double, Kokkos::extents<int, Kokkos::dynamic_extent, 1, 20>, Kokkos::layout_left> frlai(
            frlai_storage.data(), state->n_cols);
        Kokkos::mdspan<double, Kokkos::extents<int, Kokkos::dynamic_extent, 1, 20>, Kokkos::layout_left> frlanduse(
            frlanduse_storage.data(), state->n_cols);
        Kokkos::mdspan<int, Kokkos::extents<int, Kokkos::dynamic_extent, 1, 20>, Kokkos::layout_left> iland(
            iland_storage.data(), state->n_cols);

        std::vector<char> is_ice(state->n_cols, 0);
        std::vector<char> is_land(state->n_cols, 1);
        std::vector<char> is_snow(state->n_cols, 0);
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
        std::vector<char> is_seasalt(state->n_species, 0);
        std::vector<char> is_dust(state->n_species, 0);
        std::vector<double> lower_radius(state->n_species, 0.0);
        std::vector<double> upper_radius(state->n_species, 0.0);
        std::vector<char> is_gas(state->n_species, 1);

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
            is_seasalt[i] = meta.is_seasalt ? 1 : 0;
            is_dust[i] = meta.is_dust ? 1 : 0;
            lower_radius[i] = meta.lower_radius;
            upper_radius[i] = meta.upper_radius;
            is_gas[i] = meta.is_gas ? 1 : 0;
        }

        // 5. Invoke flat science bridge (casting char* vectors to bool* pointers)
        run_drydep_science_bridge(
            state->n_cols, state->n_levels, state->n_species, state->time.timestep, gas_scheme.c_str(),
            aero_scheme.c_str(), diagnostics_enabled ? 1 : 0, bxheight_ptr, airden_ptr, t_ptr, pedge_ptr, rh_ptr,
            cldfrc.data(), frlai.data_handle(), frlanduse.data_handle(), iland.data_handle(), (bool*)is_ice.data(),
            (bool*)is_land.data(), (bool*)is_snow.data(), lat_ptr, lon_ptr, obk_ptr, ps_ptr, salinity.data(),
            suncosmid.data(), swgdn.data(), ts_ptr, tskin.data(), ustar_ptr, z0.data(), frlake.data(), gwettop.data(),
            hflux_ptr, lwi.data(), pblh_ptr, u10m.data(), v10m.data(), z0h.data(), mw_g.data(), dd_f0.data(),
            dd_hstar.data(), dd_DvzAerSnow.data(), dd_DvzMinVal_snow.data(), dd_DvzMinVal_land.data(), density.data(),
            radius.data(), (bool*)is_seasalt.data(), (bool*)is_dust.data(), lower_radius.data(), upper_radius.data(),
            (bool*)is_gas.data(), conc_ptr, mock_tendency.data(), diag_con, diag_vel, diagnostic_species_id.data(),
            diagnostic_species_id.size());

        state->sync_to_device();
    }

} // namespace catchem

extern "C" {
void catchem_register_drydep_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "drydep", []() { return std::make_shared<catchem::DryDepProcess>(); });
}
}
