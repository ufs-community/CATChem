#include "catchem_process_carbchem.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_process_registry.hpp"
#include <algorithm>
#include <iostream>

extern "C" {
void run_carbchem_science_bridge(int n_cols, int n_levels, int n_species, double dt, const char* active_scheme,
                                 int diagnostics, int year, int month, int day, int hour, int minute, int second,
                                 double* airden, double* delp, double* pmid, double* species_t_chem_loss,
                                 const char* species_names_char, double* conc, double* tendency, double* diag_prod_mass,
                                 double* diag_loss_flux, double* diag_phobic_mass, double* diag_phobic_flux,
                                 const int* diagnostic_species_id, int n_diag_species);
}

namespace catchem {

    CarbChemProcess::CarbChemProcess() : active_scheme("gocart"), diagnostics_enabled(true) {}

    void CarbChemProcess::init(std::shared_ptr<StateManager> state) {
        // 1. Setup diagnostic species ID dynamically (using a dummy is_carbchem flag if we had one, but we map all
        // indices here for simplicity since CarbChem filters internally)
        for (size_t i = 0; i < state->chem.species_list.size(); ++i) {
            diagnostic_species_id.push_back(i + 1); // 1-based for Fortran bridge
        }

        // 2. Register C++ Diagnostic fields
        std::vector<int> dims_3d = {state->n_cols, state->n_levels, state->n_species};
        std::vector<int> dims_2d = {state->n_cols, state->n_species};

        state->diag_mgr->register_field("carbchem_prod_mass", "Carbon Chemistry Production Mass", "kg/kg",
                                        DiagType::FIELD_3D, dims_3d);
        state->diag_mgr->register_field("carbchem_loss_flux", "Carbon Chemistry Loss Flux", "kg/m2/s",
                                        DiagType::FIELD_2D, dims_2d);
        state->diag_mgr->register_field("carbchem_phobic_mass", "Carbon Chemistry Phobic to Philic Mass", "kg/kg",
                                        DiagType::FIELD_3D, dims_3d);
        state->diag_mgr->register_field("carbchem_phobic_flux", "Carbon Chemistry Phobic to Philic Flux", "kg/m2/s",
                                        DiagType::FIELD_2D, dims_2d);
    }

    void CarbChemProcess::run(std::shared_ptr<StateManager> state) {
        state->sync_to_host();

        // 1. Retrieve 3D Meteorological state pointers
        auto airden_ptr_it = state->met.fields_3d.find("AIRDEN_DRY");
        double* airden_ptr =
            (airden_ptr_it != state->met.fields_3d.end()) ? airden_ptr_it->second->host_view.data() : nullptr;
        auto delp_ptr_it = state->met.fields_3d.find("DELP");
        double* delp_ptr =
            (delp_ptr_it != state->met.fields_3d.end()) ? delp_ptr_it->second->host_view.data() : nullptr;
        auto pmid_ptr_it = state->met.fields_3d.find("PMID");
        double* pmid_ptr =
            (pmid_ptr_it != state->met.fields_3d.end()) ? pmid_ptr_it->second->host_view.data() : nullptr;

        if (!airden_ptr || !delp_ptr || !pmid_ptr) {
            std::cerr << "CarbChemProcess: Missing required meteorological fields." << std::endl;
            return;
        }

        // 2. Diagnostic Views
        double* diag_prod_mass = nullptr;
        double* diag_loss_flux = nullptr;
        double* diag_phobic_mass = nullptr;
        double* diag_phobic_flux = nullptr;

        if (state->diag_mgr && diagnostics_enabled) {
            diag_prod_mass = (double*)state->diag_mgr->get_host_pointer("carbchem_prod_mass");
            diag_loss_flux = (double*)state->diag_mgr->get_host_pointer("carbchem_loss_flux");
            diag_phobic_mass = (double*)state->diag_mgr->get_host_pointer("carbchem_phobic_mass");
            diag_phobic_flux = (double*)state->diag_mgr->get_host_pointer("carbchem_phobic_flux");
        }

        double* conc_ptr = state->chem.conc ? state->chem.conc->host_view.data() : nullptr;

        // Allocate local tendencies buffer
        std::vector<double> mock_tendency(state->n_cols * state->n_levels * state->n_species, 0.0);

        // 4. Retrieve species properties from ChemState
        std::vector<double> t_chem_loss(state->n_species, 0.0);
        for (size_t i = 0; i < state->chem.species_list.size(); ++i) {
            // Mock t_chem_loss for carbchem if not directly in metadata struct
            t_chem_loss[i] = 1.0;
        }

        // 5. Invoke flat science bridge
        run_carbchem_science_bridge(
            state->n_cols, state->n_levels, state->n_species, state->time.timestep, active_scheme.c_str(),
            diagnostics_enabled ? 1 : 0, state->time.year, state->time.month, state->time.day, state->time.hour,
            state->time.minute, state->time.second, airden_ptr, delp_ptr, pmid_ptr, t_chem_loss.data(),
            state->chem.species_names_c_arr.data(), conc_ptr, mock_tendency.data(), diag_prod_mass, diag_loss_flux,
            diag_phobic_mass, diag_phobic_flux, diagnostic_species_id.data(), diagnostic_species_id.size());
    }

    void CarbChemProcess::finalize() {}

} // namespace catchem

extern "C" {
void catchem_register_carbchem_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "carbchem", []() { return std::make_shared<catchem::CarbChemProcess>(); });
}
}
