#include "catchem_process_seasalt.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_process_registry.hpp"
#include <iostream>

extern "C" {
void run_seasalt_science_bridge(int n_cols, int n_levels, int n_species, double dt, const char* active_scheme,
                                int diagnostics, double* frocean, double* frseaice, double* lat, double* lon,
                                double* sst, double* u10m, double* v10m, double* ustar, double* delp, double* density,
                                double* radius, double* lower_radius, double* upper_radius, bool* is_gas, double* mw_g,
                                double* conc, double* tendency, double* diag_mass_total, double* diag_num_total,
                                double* diag_mass_bin, double* diag_num_bin, const int* diagnostic_species_id,
                                int n_diag_species);
}

namespace catchem {

    SeaSaltProcess::SeaSaltProcess() : active_scheme("geos12"), diagnostics_enabled(true) {}

    void SeaSaltProcess::init(std::shared_ptr<StateManager> state) {
        if (state->diag_mgr) {
            std::vector<int> dims_1d = {state->n_cols, 1};
            state->diag_mgr->register_field("seasalt_mass_emission_total", "Total Mass Emission", "kg/m2/s",
                                            DiagType::FIELD_2D, dims_1d);
            state->diag_mgr->register_field("seasalt_number_emission_total", "Total Number Emission", "#/m2/s",
                                            DiagType::FIELD_2D, dims_1d);

            for (size_t i = 0; i < state->chem.species_list.size(); ++i) {
                auto& meta = state->chem.species_list[i];
                if (meta.is_seasalt) {
                    std::string mass_name = "seasalt_mass_emission_" + meta.short_name;
                    std::string num_name = "seasalt_number_emission_" + meta.short_name;
                    state->diag_mgr->register_field(mass_name, "Mass Emission " + meta.short_name, "kg/m2/s",
                                                    DiagType::FIELD_2D, dims_1d);
                    state->diag_mgr->register_field(num_name, "Number Emission " + meta.short_name, "#/m2/s",
                                                    DiagType::FIELD_2D, dims_1d);
                }
            }
        }
    }

    void SeaSaltProcess::run(std::shared_ptr<StateManager> state) {
        state->sync_to_host();

        // 1. Fetch raw pointers to Met fields dynamically from unordered_maps inside state->met
        auto frocean_it = state->met.fields_2d.find("FROCEAN");
        double* frocean_ptr =
            (frocean_it != state->met.fields_2d.end()) ? frocean_it->second->host_view.data() : nullptr;

        auto frseaice_it = state->met.fields_2d.find("FRSEAICE");
        double* frseaice_ptr =
            (frseaice_it != state->met.fields_2d.end()) ? frseaice_it->second->host_view.data() : nullptr;

        auto sst_it = state->met.fields_2d.find("SST");
        double* sst_ptr = (sst_it != state->met.fields_2d.end()) ? sst_it->second->host_view.data() : nullptr;

        auto lat_it = state->met.fields_2d.find("LAT");
        double* lat_ptr = (lat_it != state->met.fields_2d.end()) ? lat_it->second->host_view.data() : nullptr;

        auto lon_it = state->met.fields_2d.find("LON");
        double* lon_ptr = (lon_it != state->met.fields_2d.end()) ? lon_it->second->host_view.data() : nullptr;

        auto delp_it = state->met.fields_3d.find("DELP");
        double* delp_ptr = (delp_it != state->met.fields_3d.end()) ? delp_it->second->host_view.data() : nullptr;

        double* ustar_ptr = state->met.USTAR ? state->met.USTAR->host_view.data() : nullptr;

        // Local fallbacks for winds if not bound
        std::vector<double> u10m(state->n_cols, 5.0);
        std::vector<double> v10m(state->n_cols, 2.0);

        // 2. Identify and slice SeaSalt-only chemical species to comply with science solver limits
        std::vector<int> ss_global_indices;
        std::vector<double> density;
        std::vector<double> radius;
        std::vector<double> lower_radius;
        std::vector<double> upper_radius;
        std::vector<char> is_gas;
        std::vector<double> mw_g;

        for (size_t i = 0; i < state->chem.species_list.size(); ++i) {
            auto& meta = state->chem.species_list[i];
            if (meta.is_seasalt) {
                std::cout << "DEBUG SS: species=" << meta.short_name << ", lower_radius=" << meta.lower_radius
                          << ", upper_radius=" << meta.upper_radius << std::endl;
                ss_global_indices.push_back(i);
                density.push_back(meta.density);
                radius.push_back(meta.radius);
                lower_radius.push_back(meta.lower_radius);
                upper_radius.push_back(meta.upper_radius);
                is_gas.push_back(meta.is_gas ? 1 : 0);
                mw_g.push_back(meta.mw_g);
            }
        }

        int n_seasalt = ss_global_indices.size();
        if (n_seasalt == 0)
            return; // No sea salt species configured

        // Extract chemical concentration raw host pointer
        double* conc_ptr = state->chem.conc ? state->chem.conc->host_view.data() : nullptr;
        if (!conc_ptr)
            return;

        // Allocate contiguous temporary slice for concentrations
        std::vector<double> sliced_conc(state->n_cols * state->n_levels * n_seasalt, 0.0);
        for (int i = 0; i < n_seasalt; ++i) {
            int g_idx = ss_global_indices[i];
            for (int col = 0; col < state->n_cols; ++col) {
                for (int lvl = 0; lvl < state->n_levels; ++lvl) {
                    int src_idx = col + lvl * state->n_cols + g_idx * state->n_cols * state->n_levels;
                    int dest_idx = col + lvl * state->n_cols + i * state->n_cols * state->n_levels;
                    sliced_conc[dest_idx] = conc_ptr[src_idx];
                }
            }
        }

        // Allocate local tendencies buffer for sliced subset
        std::vector<double> mock_tendency(state->n_cols * state->n_levels * n_seasalt, 0.0);

        // 3. Extract diagnostics
        double* diag_mass_total_ptr =
            state->diag_mgr ? (double*)state->diag_mgr->get_host_pointer("seasalt_mass_emission_total") : nullptr;
        double* diag_num_total_ptr =
            state->diag_mgr ? (double*)state->diag_mgr->get_host_pointer("seasalt_number_emission_total") : nullptr;

        std::vector<double> diag_mass_bin(state->n_cols * n_seasalt, 0.0);
        std::vector<double> diag_num_bin(state->n_cols * n_seasalt, 0.0);

        // Dynamic ID array mapping diagnostic species (1-based index in the sliced sea salt subset!)
        std::vector<int> diagnostic_species_id(n_seasalt);
        for (int i = 0; i < n_seasalt; ++i) {
            diagnostic_species_id[i] = i + 1;
        }

        // 5. Invoke flat science bridge
        run_seasalt_science_bridge(state->n_cols, state->n_levels, n_seasalt, state->time.timestep,
                                   active_scheme.c_str(), diagnostics_enabled ? 1 : 0, frocean_ptr, frseaice_ptr,
                                   lat_ptr, lon_ptr, sst_ptr, u10m.data(), v10m.data(), ustar_ptr, delp_ptr,
                                   density.data(), radius.data(), lower_radius.data(), upper_radius.data(),
                                   (bool*)is_gas.data(), mw_g.data(), sliced_conc.data(), mock_tendency.data(),
                                   diag_mass_total_ptr, diag_num_total_ptr, diag_mass_bin.data(), diag_num_bin.data(),
                                   diagnostic_species_id.data(), diagnostic_species_id.size());

        // 6. Copy sliced concentrations back to main unified chemistry state
        for (int i = 0; i < n_seasalt; ++i) {
            int g_idx = ss_global_indices[i];
            for (int col = 0; col < state->n_cols; ++col) {
                for (int lvl = 0; lvl < state->n_levels; ++lvl) {
                    int src_idx = col + lvl * state->n_cols + i * state->n_cols * state->n_levels;
                    int dest_idx = col + lvl * state->n_cols + g_idx * state->n_cols * state->n_levels;
                    conc_ptr[dest_idx] = sliced_conc[src_idx];
                }
            }
        }

        // 7. Map bin diagnostics back to dynamically registered individual C++ diagnostics
        if (state->diag_mgr && diagnostics_enabled) {
            for (int i = 0; i < n_seasalt; ++i) {
                auto& meta = state->chem.species_list[ss_global_indices[i]];
                std::string mass_name = "seasalt_mass_emission_" + meta.short_name;
                std::string num_name = "seasalt_number_emission_" + meta.short_name;
                double* mass_ptr = (double*)state->diag_mgr->get_host_pointer(mass_name);
                double* num_ptr = (double*)state->diag_mgr->get_host_pointer(num_name);
                for (int col = 0; col < state->n_cols; ++col) {
                    if (mass_ptr)
                        mass_ptr[col] = diag_mass_bin[col + i * state->n_cols];
                    if (num_ptr)
                        num_ptr[col] = diag_num_bin[col + i * state->n_cols];
                }
            }
        }

        state->sync_to_device();
    }

} // namespace catchem

extern "C" {
void catchem_register_seasalt_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "seasalt", []() { return std::make_shared<catchem::SeaSaltProcess>(); });
}
}
