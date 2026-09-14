// src/process/gaschem/catchem_process_gaschem.cpp
#include "catchem_process_gaschem.hpp"
#include "catchem_constants.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_logger.hpp"
#include "catchem_process_registry.hpp"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>

namespace catchem {

    ProcessContract GasChemProcess::get_contract() const {
        return make_contract(get_name(),
                             {host_field_3d("T", "K"), host_field_3d("PMID", "Pa"),
                              host_field_3d("AIRDEN_DRY", "kg/m3", FieldRequirement::Optional), host_concentration()});
    }

    GasChemProcess::GasChemProcess() = default;
    GasChemProcess::~GasChemProcess() = default;

    void GasChemProcess::init(std::shared_ptr<StateManager> state) {
        Logger::debug(state.get(), "GasChemProcess::init started");

        // 1. Resolve configuration path dynamically via ConfigManager. An explicit
        // "config_file" selects a single-file MUSICA mechanism (v1 format) and takes
        // precedence over the "config_dir" directory (v0 CAMP layout) fallback.
        if (state->config_manager()) {
            std::string file = state->config_manager()->get_string("processes/gaschem/config_file", "");
            if (file.empty()) {
                file = state->config_manager()->get_string("process/gaschem/config_file", "");
            }
            if (!file.empty()) {
                this->config_file = file;
            }

            std::string dir = state->config_manager()->get_string("processes/gaschem/config_dir", "");
            if (dir.empty()) {
                dir = state->config_manager()->get_string("process/gaschem/config_dir", "");
            }
            if (!dir.empty()) {
                this->config_dir = dir;
            }
        }

        if (this->config_dir.empty()) {
            if (!state->configuration_path().empty()) {
                std::string path = state->configuration_path();
                size_t last_slash = path.find_last_of("/\\");
                if (last_slash != std::string::npos) {
                    this->config_dir = path.substr(0, last_slash + 1);
                } else {
                    this->config_dir = "./";
                }
            } else {
                this->config_dir = "tests/Configs/Default/";
            }
        }

        // A directory may also hold a single-file (v1) mechanism; MUSICA's universal
        // parser only routes directories to the v0 CAMP parser, so auto-detect the
        // v1 YAML entry point (config.yaml / config.yml) before falling back to the
        // directory itself. config.json is intentionally skipped: in a v0 CAMP layout
        // it is the version-less camp-files manifest and must stay a directory parse.
        if (this->config_file.empty() && !this->config_dir.empty()) {
            for (const char* candidate : {"config.yaml", "config.yml"}) {
                std::filesystem::path probe = std::filesystem::path(this->config_dir) / candidate;
                if (std::filesystem::is_regular_file(probe)) {
                    this->config_file = probe.string();
                    break;
                }
            }
        }

        const std::string micm_config_path = this->config_file.empty() ? this->config_dir : this->config_file;
        Logger::info(
            state.get(), "GasChemProcess: resolved MICM configuration",
            {{"path", micm_config_path}, {"type", this->config_file.empty() ? "directory (v0 CAMP)" : "file"}});

        // 2. Initialize MICM and State using musica library
        try {
            micm_instance = std::make_unique<musica::MICM>(micm_config_path, musica::RosenbrockStandardOrder);
            micm_state = std::make_unique<musica::State>(*micm_instance, state->column_count() * state->level_count());
            initialized = true;
            std::clog << "[INFO] GasChemProcess: initialized MICM successfully!" << std::endl;

            // Validate that all active CATChem species are mapped inside the MICM solver
            auto variable_map = micm_state->GetVariableMap();
            for (int ispec = 0; ispec < state->species_count(); ++ispec) {
                std::string name = state->chemistry().species_list[ispec].short_name;
                for (auto& c : name)
                    c = std::toupper(c);
                if (variable_map.find(name) == variable_map.end()) {
                    Logger::warn(state.get(), "Active species not found in MICM solver variable map",
                                 {{"species", name}});
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "GasChemProcess: Error: failed to initialize MICM: " << e.what() << std::endl;
            initialized = false;
            throw std::runtime_error(std::string("GasChemProcess: failed to initialize MICM: ") + e.what());
        }
    }

    void GasChemProcess::run(std::shared_ptr<StateManager> state) {
        if (!initialized) {
            std::cerr << "GasChemProcess: Warning: skipped run because solver is not initialized." << std::endl;
            return;
        }

        // 1. Sync device to host

        if (!state->meteorology().AIRDEN_DRY && state->meteorology().PMID && state->meteorology().T) {
            state->derive_airden_dry();
        }

        if (!state->meteorology().T || !state->meteorology().PMID || !state->meteorology().AIRDEN_DRY ||
            !state->chemistry().conc) {
            throw std::runtime_error(
                "GasChem requires current T, PMID, AIRDEN_DRY, and chemistry concentration fields");
        }

        auto temp = state->meteorology().T->host_view;
        auto pmid = state->meteorology().PMID->host_view;
        auto airden_dry = state->meteorology().AIRDEN_DRY->host_view;
        state->chemistry().conc->host_write();
        auto conc = state->chemistry().conc->host_view;

        auto& micm_conditions = micm_state->GetConditions();
        auto& micm_concs = micm_state->GetOrderedConcentrations();
        auto& micm_rate_params = micm_state->GetOrderedRateParameters();

        int nc = state->column_count();
        int nl = state->level_count();
        int ns = state->species_count();

        size_t vector_size_ = micm_instance->GetVectorSize();
        auto variable_map = micm_state->GetVariableMap();
        size_t n_micm_species = variable_map.size();

        auto rate_param_map = micm_state->GetRateParameterMap();
        size_t n_rate_params = rate_param_map.size();

        // Dry air molecular weight in kg/mol (sourced from catchem::constants)
        constexpr double air_mw_kg = catchem::constants::AIR_MW * 1.0e-3;

        // 2. Map environmental variables and input concentrations to state
        for (int ilev = 0; ilev < nl; ++ilev) {
            for (int icol = 0; icol < nc; ++icol) {
                int i_cell = ilev * nc + icol;

                double t_val = temp(icol, ilev, 0);
                double p_val = pmid(icol, ilev, 0);
                double density_dry_kg = airden_dry(icol, ilev, 0);

                if (!std::isfinite(t_val) || t_val <= 0.0 || !std::isfinite(p_val) || p_val <= 0.0 ||
                    !std::isfinite(density_dry_kg) || density_dry_kg <= 0.0)
                    throw std::domain_error("GasChem received non-physical T, PMID, or AIRDEN_DRY");

                // Convert dry air density: kg/m3 to mol/m3
                double air_density_mol = density_dry_kg / air_mw_kg;

                micm_conditions[i_cell].temperature_ = t_val;
                micm_conditions[i_cell].pressure_ = p_val;
                micm_conditions[i_cell].air_density_ = air_density_mol;

                // Copy concentrations: ppmv -> mol/m3
                for (int ispec = 0; ispec < ns; ++ispec) {
                    std::string name = state->chemistry().species_list[ispec].short_name;
                    for (auto& c : name)
                        c = std::toupper(c);

                    auto it = variable_map.find(name);
                    if (it != variable_map.end()) {
                        size_t i_micm_spec = it->second;
                        double ppmv_val = conc(icol, ilev, ispec);
                        if (!std::isfinite(ppmv_val) || ppmv_val < 0.0)
                            throw std::domain_error("GasChem received a negative or non-finite concentration for " +
                                                    name);

                        double conc_molar = ppmv_val * 1.0e-6 * air_density_mol;

                        size_t group_index = i_cell / vector_size_;
                        size_t row_in_group = i_cell % vector_size_;
                        size_t idx = (group_index * n_micm_species + i_micm_spec) * vector_size_ + row_in_group;
                        micm_concs[idx] = conc_molar;
                    }
                }

                // 3. Dynamic Photolysis Mapping (PHOTO.<label> to photolysis_rate_<label>)
                for (const auto& [param_name, i_param] : rate_param_map) {
                    size_t group_index = i_cell / vector_size_;
                    size_t row_in_group = i_cell % vector_size_;
                    size_t idx = (group_index * n_rate_params + i_param) * vector_size_ + row_in_group;

                    if (param_name.rfind("PHOTO.", 0) == 0) {
                        std::string label = param_name.substr(6);
                        std::string diag_name = "photolysis_rate_" + label;

                        double rate_val = 0.0;
                        if (state->diagnostic_manager() && state->diagnostic_manager()->has_field(diag_name)) {
                            double* diag_ptr =
                                static_cast<double*>(state->diagnostic_manager()->get_host_pointer(diag_name));
                            if (diag_ptr) {
                                int diag_idx = ilev * nc + icol;
                                rate_val = diag_ptr[diag_idx];
                            }
                        }
                        micm_rate_params[idx] = rate_val;
                    } else if (param_name.rfind("LOSS.", 0) == 0) {
                        micm_rate_params[idx] = 1.0;
                    }
                }
            }
        }

        // 4. Run standard CPU solver
        double tstep = state->clock().timestep;
        if (tstep <= 0.0) {
            Logger::error(state.get(), "Invalid timestep encountered", {{"timestep", std::to_string(tstep)}});
            throw std::runtime_error("GasChemProcess: timestep must be greater than zero.");
        }
        auto solver_result = micm_instance->Solve(micm_state.get(), tstep);
        if (solver_result.state_ != micm::SolverState::Converged &&
            solver_result.state_ != micm::SolverState::AcceptingUnconvergedIntegration) {
            Logger::error(state.get(), "MICM Solver did not reach convergence!",
                          {{"final_state", micm::SolverStateToString(solver_result.state_)}});
            throw std::runtime_error("GasChemProcess: MICM Solver failed to reach convergence or acceptable state: " +
                                     micm::SolverStateToString(solver_result.state_));
        }

        // 5. Convert output concentrations back: mol/m3 -> ppmv
        for (int ilev = 0; ilev < nl; ++ilev) {
            for (int icol = 0; icol < nc; ++icol) {
                int i_cell = ilev * nc + icol;
                double air_density_mol = micm_conditions[i_cell].air_density_;

                for (int ispec = 0; ispec < ns; ++ispec) {
                    std::string name = state->chemistry().species_list[ispec].short_name;
                    for (auto& c : name)
                        c = std::toupper(c);

                    auto it = variable_map.find(name);
                    if (it != variable_map.end()) {
                        size_t i_micm_spec = it->second;

                        size_t group_index = i_cell / vector_size_;
                        size_t row_in_group = i_cell % vector_size_;
                        size_t idx = (group_index * n_micm_species + i_micm_spec) * vector_size_ + row_in_group;
                        double conc_molar = micm_concs[idx];

                        double ppmv_val = (conc_molar / air_density_mol) * 1.0e6;
                        if (!std::isfinite(ppmv_val) || ppmv_val < 0.0)
                            throw std::domain_error(
                                "GasChem solver returned a negative or non-finite concentration for " + name);
                        conc(icol, ilev, ispec) = ppmv_val;
                    }
                }
            }
        }

        // 6. Sync back to device
        if (state->chemistry().conc)
            state->chemistry().conc->mark_host_modified();
    }

    void GasChemProcess::finalize() {}

} // namespace catchem

extern "C" {
void catchem_register_gaschem_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        std::string(catchem::ProcessNames::GasChem), []() { return std::make_shared<catchem::GasChemProcess>(); });
}
}
