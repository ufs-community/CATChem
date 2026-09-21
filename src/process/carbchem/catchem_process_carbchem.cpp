#include "catchem_process_carbchem.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_error.hpp"
#include "catchem_logger.hpp"
#include "catchem_process_registry.hpp"
#include <algorithm>
#include <cctype>
#include <iostream>

extern "C" {
void run_carbchem_science_bridge(int n_cols, int n_levels, int n_species, double dt, const char* active_scheme,
                                 int diagnostics, double gocart_time_days_hydrophobic_to_hydrophilic, int year,
                                 int month, int day, int hour, int minute, int second, double* airden, double* delp,
                                 double* pmid, double* species_t_chem_loss, const char* species_names_char,
                                 double* conc, double* tendency, double* diag_prod_mass, double* diag_loss_flux,
                                 double* diag_phobic_mass, double* diag_phobic_flux, const int* diagnostic_species_id,
                                 int n_diag_species);
}

namespace catchem {

    ProcessContract CarbChemProcess::get_contract() const {
        return make_contract(get_name(), {host_field_3d("PMID", "Pa"), host_field_3d("DELP", "Pa"),
                                          host_field_3d("AIRDEN", "kg/m3"), host_concentration()});
    }

    CarbChemProcess::CarbChemProcess() : active_scheme("gocart"), diagnostics_enabled(true) {}

    void CarbChemProcess::prepare_inputs(std::shared_ptr<StateManager> state) {
        state->derive_delp();
        state->derive_airden();
    }

    void CarbChemProcess::init(std::shared_ptr<StateManager> state) {
        const auto config = state->config_manager();
        if (!config)
            throw std::invalid_argument("CarbChem requires a runtime YAML configuration");
        const auto configured = config->data.processes.find("carbchem");
        if (configured == config->data.processes.end() || configured->second.scheme != "gocart")
            throw std::invalid_argument("CarbChem requires processes.carbchem.scheme: gocart");
        diagnostics_enabled = configured->second.diagnostics;

        // Read scheme tuning options from the runtime YAML.  The lookup falls
        // back to the compiled default declared in CarbChemCommon_Mod.F90, so
        // a configuration that omits the option keeps current behavior.
        gocart_time_days =
            configured->second.get_double("gocart/time_days_hydrophobic_to_hydrophilic", gocart_time_days);
        if (!(gocart_time_days > 0.0))
            throw std::invalid_argument("CarbChem gocart time_days_hydrophobic_to_hydrophilic must be positive");

        // Surface the effective scheme options so the run log confirms what
        // was parsed from the runtime YAML and will be passed to the bridge.
        Logger::debug(state.get(), "CarbChem scheme options",
                      {{"scheme", active_scheme},
                       {"gocart/time_days_hydrophobic_to_hydrophilic", std::to_string(gocart_time_days)}});

        // 1. Resolve the diagnostic species set (parity with legacy).  The
        // GOCART scheme matches diagnostic_species_id(diag_idx) == species_idx
        // where species_idx is the GLOBAL catalog position (it loops over the
        // full species list it is handed), so ids live in global 1-based
        // space.  Default set = the explicit carbon species the scheme
        // converts (oc1/oc2, bc1/bc2, plus br1/br2 where configured); an
        // explicit diag_species overrides it and must name real mechanism
        // species (fail-loud).  Names resolve case-insensitively against the
        // mechanism's canonical (upper-cased) name map.
        diagnostic_species_id.clear();
        const auto& settings = configured->second;
        const std::vector<std::string> default_carbon = {"oc1", "oc2", "bc1", "bc2", "br1", "br2"};
        const auto& requested = settings.diag_species.empty() ? default_carbon : settings.diag_species;
        for (const auto& name : requested) {
            std::string canonical = name;
            std::transform(canonical.begin(), canonical.end(), canonical.begin(),
                           [](unsigned char c) { return std::toupper(c); });
            const auto found = state->chemistry().species_name_to_index.find(canonical);
            if (found == state->chemistry().species_name_to_index.end()) {
                if (settings.diag_species.empty())
                    continue; // default set: tolerate species absent from this mechanism
                throw std::invalid_argument("CarbChem diag_species names an unknown species: " + name);
            }
            diagnostic_species_id.push_back(found->second + 1); // global, 1-based for the bridge
        }

        if (!diagnostics_enabled)
            return;
        // 2. Register C++ Diagnostic fields.  The species dimension is the
        // (possibly defaulted) diagnostic count, not the full catalog: the
        // scheme scatters each species into its diag_idx slot, so packed
        // [ncol, nz, ndiag] / [ncol, ndiag] layouts match the bridge shape.
        const int n_diag = static_cast<int>(diagnostic_species_id.size());
        // A mechanism may configure none of the default carbon species; with an
        // empty packed axis there is nothing to register (matches the seasalt/
        // settling guards).
        if (n_diag > 0) {
            std::vector<int> dims_3d = {state->column_count(), state->level_count(), n_diag};
            std::vector<int> dims_2d = {state->column_count(), n_diag};

            // Packed-axis labels: diagnostic_species_id holds GLOBAL 1-based
            // catalog positions (the space the GOCART scheme matches), so slot
            // i's label is that species' short_name (FR-006).  The NUOPC driver
            // unpacks each field into one named variable per species (feature 013).
            std::vector<std::string> carbon_labels;
            carbon_labels.reserve(diagnostic_species_id.size());
            for (const int gid : diagnostic_species_id)
                carbon_labels.push_back(state->chemistry().species_list[static_cast<size_t>(gid) - 1].short_name);
            const std::vector<SemanticAxis> axes_3d = {SemanticAxis::Column, SemanticAxis::Level,
                                                       SemanticAxis::Species};
            const std::vector<SemanticAxis> axes_2d = {SemanticAxis::Column, SemanticAxis::Species};

            state->diagnostic_manager()->register_field_contract(
                "carbchem_prod_mass", "Carbon Chemistry Production Mass", "kg/kg", DiagType::FIELD_3D, dims_3d,
                DiagnosticPolicy::Instantaneous, 0.0, axes_3d, carbon_labels);
            state->diagnostic_manager()->register_field_contract(
                "carbchem_loss_flux", "Carbon Chemistry Loss Flux", "kg/m2/s", DiagType::FIELD_2D, dims_2d,
                DiagnosticPolicy::Instantaneous, 0.0, axes_2d, carbon_labels);
            state->diagnostic_manager()->register_field_contract(
                "carbchem_phobic_mass", "Carbon Chemistry Phobic to Philic Mass", "kg/kg", DiagType::FIELD_3D, dims_3d,
                DiagnosticPolicy::Instantaneous, 0.0, axes_3d, carbon_labels);
            state->diagnostic_manager()->register_field_contract(
                "carbchem_phobic_flux", "Carbon Chemistry Phobic to Philic Flux", "kg/m2/s", DiagType::FIELD_2D,
                dims_2d, DiagnosticPolicy::Instantaneous, 0.0, axes_2d, carbon_labels);
        }
    }

    void CarbChemProcess::run(std::shared_ptr<StateManager> state) {

        // The execution plan invokes prepare_inputs before run(), but direct
        // API users and focused process tests may call run() themselves.
        // These derivations are generation-aware and therefore preserve
        // host-provided fields while supplying only absent prerequisites
        // (matches SettlingProcess::run).
        prepare_inputs(state);

        // 1. Retrieve 3D Meteorological state pointers
        double* airden_ptr = state->write_field<3>("AIRDEN");

        double* delp_ptr = state->write_field<3>("DELP");

        double* pmid_ptr = state->write_field<3>("PMID");

        require_field_pointer("CarbChem", "AIRDEN", airden_ptr);
        require_field_pointer("CarbChem", "DELP", delp_ptr);
        require_field_pointer("CarbChem", "PMID", pmid_ptr);

        // 2. Diagnostic Views
        double* diag_prod_mass = nullptr;
        double* diag_loss_flux = nullptr;
        double* diag_phobic_mass = nullptr;
        double* diag_phobic_flux = nullptr;

        if (state->diagnostic_manager() && diagnostics_enabled) {
            diag_prod_mass = (double*)state->diagnostic_manager()->get_host_pointer("carbchem_prod_mass");
            diag_loss_flux = (double*)state->diagnostic_manager()->get_host_pointer("carbchem_loss_flux");
            diag_phobic_mass = (double*)state->diagnostic_manager()->get_host_pointer("carbchem_phobic_mass");
            diag_phobic_flux = (double*)state->diagnostic_manager()->get_host_pointer("carbchem_phobic_flux");
        }

        double* conc_ptr = state->chemistry().conc ? state->chemistry().conc->host_write() : nullptr;
        require_field_pointer("CarbChem", "CHEM_CONC", conc_ptr);

        // Allocate local tendencies buffer
        std::vector<double> mock_tendency(state->column_count() * state->level_count() * state->species_count(), 0.0);

        // 4. Retrieve species properties from ChemState
        // Preserve negative values: GOCART's carbonChemLoss uses tChemLoss<0 as
        // its own "loss disabled" sentinel and early-returns on it; clamping to
        // 0 here instead makes it compute exp(-cdt/0), annihilating the species.
        std::vector<double> t_chem_loss(state->species_count(), -1.0);
        for (size_t i = 0; i < state->chemistry().species_list.size(); ++i) {
            t_chem_loss[i] = state->chemistry().species_list[i].t_chem_loss;
        }

        // 5. Invoke flat science bridge.  diagnostic_species_id (built in
        // init) holds GLOBAL 1-based catalog indices — the space the GOCART
        // scheme matches on.  When diagnostics are off the bridge still
        // forwards the array to the scheme, so pass a valid size-1 dummy
        // holding 0 — it matches no 1-based species_idx and the scheme writes
        // nothing (mirrors the settling/dust/seasalt no_diag_species guard).
        const bool forward_diag = diagnostics_enabled && !diagnostic_species_id.empty();
        static const int no_diag_species = 0;
        const int* diag_ids = forward_diag ? diagnostic_species_id.data() : &no_diag_species;
        const int n_diag_species = forward_diag ? static_cast<int>(diagnostic_species_id.size()) : 0;
        run_carbchem_science_bridge(
            state->column_count(), state->level_count(), state->species_count(), state->clock().timestep,
            active_scheme.c_str(), diagnostics_enabled ? 1 : 0, gocart_time_days, state->clock().year,
            state->clock().month, state->clock().day, state->clock().hour, state->clock().minute, state->clock().second,
            airden_ptr, delp_ptr, pmid_ptr, t_chem_loss.data(), state->chemistry().species_names_c_arr.data(), conc_ptr,
            mock_tendency.data(), diag_prod_mass, diag_loss_flux, diag_phobic_mass, diag_phobic_flux, diag_ids,
            n_diag_species);

        if (state->chemistry().conc)
            state->chemistry().conc->mark_host_modified();
    }

    void CarbChemProcess::finalize() {}

} // namespace catchem

extern "C" {
void catchem_register_carbchem_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "carbchem", []() { return std::make_shared<catchem::CarbChemProcess>(); }, {},
        catchem::make_settings_validator("carbchem", {"gocart/time_days_hydrophobic_to_hydrophilic"}));
}
}
