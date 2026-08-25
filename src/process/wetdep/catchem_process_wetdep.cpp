#include "catchem_process_wetdep.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_process_registry.hpp"
#include <iostream>

extern "C" {
void run_wetdep_science_bridge(int n_cols, int n_levels, int n_species, double dt, int diagnostics, double* airden_dry,
                               double* mairden, double* pedge, double* pfilsan, double* pfllsan, double* reevapls,
                               double* t_air, bool* is_aerosol, double* henry_cr, double* henry_k0, double* henry_pKa,
                               double* wd_retfactor, bool* wd_LiqAndGas, double* wd_convfacI2G, double* wd_rainouteff,
                               double* wd_reevap_frac, double* radius, double* mw_g, const char* species_names,
                               double* conc, double* tendency, double* diag_mass, double* diag_flux,
                               const int* diagnostic_species_id, int n_diag_species);
}

namespace catchem {

    WetDepProcess::WetDepProcess() : active_scheme("jacob"), diagnostics_enabled(true) {}

    void WetDepProcess::init(std::shared_ptr<StateManager> state) {
        if (state->diag_mgr) {
            std::vector<int> dims_2d = {state->n_cols, state->n_levels};
            for (size_t i = 0; i < state->chem.species_list.size(); ++i) {
                auto& meta = state->chem.species_list[i];
                if (meta.is_wetdep) {
                    std::string mass_name = "wetdep_mass_" + meta.short_name;
                    std::string flux_name = "wetdep_flux_" + meta.short_name;
                    state->diag_mgr->register_field(mass_name, "Wet Mass " + meta.short_name, "kg/m2",
                                                    DiagType::FIELD_2D, dims_2d);
                    state->diag_mgr->register_field(flux_name, "Wet Flux " + meta.short_name, "kg/m2/s",
                                                    DiagType::FIELD_2D, dims_2d);

                    // Track diagnostic species index (1-based)
                    diagnostic_species_id.push_back(i + 1);
                }
            }
        }
    }

    void WetDepProcess::run(std::shared_ptr<StateManager> state) {
        state->sync_to_host();

        // 1. Fetch raw pointers to Met Views
        double* airden_dry_ptr = state->met.AIRDEN_DRY ? state->met.AIRDEN_DRY->host_view.data() : nullptr;
        double* mairden_ptr = state->met.AIRDEN ? state->met.AIRDEN->host_view.data() : nullptr;
        double* pedge_ptr = state->met.PEDGE ? state->met.PEDGE->host_view.data() : nullptr;
        double* t_ptr = state->met.T ? state->met.T->host_view.data() : nullptr;

        auto pfilsan_it = state->met.fields_3d.find("PFILSAN");
        double* pfilsan_ptr =
            (pfilsan_it != state->met.fields_3d.end()) ? pfilsan_it->second->host_view.data() : nullptr;

        auto pfllsan_it = state->met.fields_3d.find("PFLLSAN");
        double* pfllsan_ptr =
            (pfllsan_it != state->met.fields_3d.end()) ? pfllsan_it->second->host_view.data() : nullptr;

        auto reevapls_it = state->met.fields_3d.find("REEVAPLS");
        double* reevapls_ptr =
            (reevapls_it != state->met.fields_3d.end()) ? reevapls_it->second->host_view.data() : nullptr;

        // 2. Extract chemical arrays & C++ allocated diagnostics
        double* conc_ptr = state->chem.conc ? state->chem.conc->host_view.data() : nullptr;

        // Allocate local tendencies buffer
        std::vector<double> mock_tendency(state->n_cols * state->n_levels * state->n_species, 0.0);

        // Allocate 3D diagnostic buffers for Fortran bridge
        std::vector<double> diag_mass_bin(state->n_cols * state->n_levels * state->n_species, 0.0);
        std::vector<double> diag_flux_bin(state->n_cols * state->n_levels * state->n_species, 0.0);

        // 3. Extract species configuration properties from ChemState
        std::vector<char> is_aerosol(state->n_species, 0);
        std::vector<double> henry_cr(state->n_species, 0.0);
        std::vector<double> henry_k0(state->n_species, 0.0);
        std::vector<double> henry_pKa(state->n_species, 0.0);
        std::vector<double> wd_retfactor(state->n_species, 0.0);
        std::vector<char> wd_LiqAndGas(state->n_species, 0);
        std::vector<double> wd_convfacI2G(state->n_species, 0.0);
        std::vector<double> wd_reevap_frac(state->n_species, 0.0);
        std::vector<double> wd_rainouteff_storage(state->n_species * 3, 0.0);
        std::vector<double> radius(state->n_species, 1e-6);
        std::vector<double> mw_g(state->n_species, 29.0);

        // Dynamic 2D view for rainouteff with species as dimension 0, and 3-element efficiency as dimension 1
        Kokkos::mdspan<double, Kokkos::extents<int, Kokkos::dynamic_extent, 3>, Kokkos::layout_left> wd_rainouteff(
            wd_rainouteff_storage.data(), state->n_species);

        for (size_t i = 0; i < state->chem.species_list.size(); ++i) {
            auto& meta = state->chem.species_list[i];
            is_aerosol[i] = meta.is_aerosol ? 1 : 0;
            henry_k0[i] = meta.henry_k0;
            henry_cr[i] = meta.henry_cr;
            henry_pKa[i] = meta.henry_pKa;
            wd_retfactor[i] = meta.wd_retfactor;
            wd_LiqAndGas[i] = meta.wd_LiqAndGas ? 1 : 0;
            wd_convfacI2G[i] = meta.wd_convfacI2G;
            wd_reevap_frac[i] = 1.0; // dummy default
            radius[i] = meta.radius;
            mw_g[i] = meta.mw_g;

            // Fill rainouteff safely up to 3 efficiency factors
            for (int k = 0; k < 3; ++k) {
                if (meta.wd_rainouteff.size() > (size_t)k) {
                    wd_rainouteff(i, k) = meta.wd_rainouteff[k];
                } else {
                    wd_rainouteff(i, k) = 0.0;
                }
            }
        }

        // 4. Invoke flat science bridge
        run_wetdep_science_bridge(
            state->n_cols, state->n_levels, state->n_species, state->time.timestep, diagnostics_enabled ? 1 : 0,
            airden_dry_ptr, mairden_ptr, pedge_ptr, pfilsan_ptr, pfllsan_ptr, reevapls_ptr, t_ptr,
            (bool*)is_aerosol.data(), henry_cr.data(), henry_k0.data(), henry_pKa.data(), wd_retfactor.data(),
            (bool*)wd_LiqAndGas.data(), wd_convfacI2G.data(), wd_rainouteff.data_handle(), wd_reevap_frac.data(),
            radius.data(), mw_g.data(), state->chem.species_names_c_arr.data(), conc_ptr, mock_tendency.data(),
            diag_mass_bin.data(), diag_flux_bin.data(), diagnostic_species_id.data(), diagnostic_species_id.size());

        // 5. Map 3D bin diagnostics back to dynamically registered individual C++ diagnostics
        if (state->diag_mgr && diagnostics_enabled) {
            for (size_t i = 0; i < state->chem.species_list.size(); ++i) {
                auto& meta = state->chem.species_list[i];
                if (meta.is_wetdep) {
                    std::string mass_name = "wetdep_mass_" + meta.short_name;
                    std::string flux_name = "wetdep_flux_" + meta.short_name;
                    double* mass_ptr = (double*)state->diag_mgr->get_host_pointer(mass_name);
                    double* num_ptr = (double*)state->diag_mgr->get_host_pointer(flux_name);
                    for (int col = 0; col < state->n_cols; ++col) {
                        for (int lvl = 0; lvl < state->n_levels; ++lvl) {
                            int idx = col + lvl * state->n_cols + i * state->n_cols * state->n_levels;
                            if (mass_ptr)
                                mass_ptr[col + lvl * state->n_cols] = diag_mass_bin[idx];
                            if (num_ptr)
                                num_ptr[col + lvl * state->n_cols] = diag_flux_bin[idx];
                        }
                    }
                }
            }
        }

        state->sync_to_device();
    }

} // namespace catchem

extern "C" {
void catchem_register_wetdep_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "wetdep", []() { return std::make_shared<catchem::WetDepProcess>(); });
}
}
