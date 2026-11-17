// src/process/gaschem/catchem_process_gaschem.cpp
#include "catchem_process_gaschem.hpp"
#include "catchem_constants.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_logger.hpp"
#include "catchem_process_registry.hpp"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace catchem {

    namespace {

        // Builds the message for the concentration guards in GasChemProcess::run.
        // The value is written in scientific notation because fixed six-decimal
        // formatting prints every tiny negative as "-0.000000" and hides its
        // magnitude. Column, level, and batch are zero-based; column is the index
        // on this PET, not within the batch.
        std::string concentration_guard_message(const std::string& what, double value, const std::string& species,
                                                int column, int level, int batch) {
            std::ostringstream message;
            message << what << " (" << std::scientific << std::setprecision(6) << value << ") for " << species
                    << " at column=" << column << ", level=" << level << ", batch=" << batch << " (zero-based indices)";
            return message.str();
        }

    } // namespace

    ProcessContract GasChemProcess::get_contract() const {
        return make_contract(get_name(),
                             {host_field_3d("T", "K"), host_field_3d("PMID", "Pa"),
                              host_field_3d("AIRDEN_DRY", "kg/m3", FieldRequirement::Optional), host_concentration()});
    }

    GasChemProcess::GasChemProcess() = default;
    GasChemProcess::~GasChemProcess() = default;

    void GasChemProcess::init(std::shared_ptr<StateManager> state) {
        Logger::debug(*state, "GasChemProcess::init started");

        // 1. Resolve configuration path dynamically via ConfigManager. An explicit
        // "config_file" selects a single-file MUSICA mechanism (v1 format) and takes
        // precedence over the "config_dir" directory (v0 CAMP layout) fallback.
        if (state->config_manager()) {
            std::string file = state->config_manager()->get_string("processes/gaschem/config_file", "");
            if (file.empty()) {
                file = state->config_manager()->get_string("process/gaschem/config_file", "");
            }
            if (!file.empty()) {
                this->config_file_ = file;
            }

            std::string dir = state->config_manager()->get_string("processes/gaschem/config_dir", "");
            if (dir.empty()) {
                dir = state->config_manager()->get_string("process/gaschem/config_dir", "");
            }
            if (!dir.empty()) {
                this->config_dir_ = dir;
            }
        }

        int cols_per_batch = 1;
        if (state->config_manager()) {
            if (!state->config_manager()->get_string("processes/gaschem/columns_per_batch", "").empty()) {
                cols_per_batch = state->config_manager()->get_int("processes/gaschem/columns_per_batch", 1);
            } else if (!state->config_manager()->get_string("process/gaschem/columns_per_batch", "").empty()) {
                cols_per_batch = state->config_manager()->get_int("process/gaschem/columns_per_batch", 1);
            }
        }
        if (cols_per_batch < 1) {
            throw std::invalid_argument("GasChemProcess: columns_per_batch must be >= 1, got " +
                                        std::to_string(cols_per_batch));
        }
        this->columns_per_batch_ = static_cast<std::size_t>(cols_per_batch);

        if (this->config_dir_.empty()) {
            if (!state->configuration_path().empty()) {
                std::string path = state->configuration_path();
                size_t last_slash = path.find_last_of("/\\");
                if (last_slash != std::string::npos) {
                    this->config_dir_ = path.substr(0, last_slash + 1);
                } else {
                    this->config_dir_ = "./";
                }
            } else {
                this->config_dir_ = "tests/Configs/Default/";
            }
        }

        // A directory may also hold a single-file (v1) mechanism; MUSICA's universal
        // parser only routes directories to the v0 CAMP parser, so auto-detect the
        // v1 YAML entry point (config.yaml / config.yml) before falling back to the
        // directory itself. config.json is intentionally skipped: in a v0 CAMP layout
        // it is the version-less camp-files manifest and must stay a directory parse.
        if (this->config_file_.empty() && !this->config_dir_.empty()) {
            for (const char* candidate : {"config.yaml", "config.yml"}) {
                std::filesystem::path probe = std::filesystem::path(this->config_dir_) / candidate;
                if (std::filesystem::is_regular_file(probe)) {
                    this->config_file_ = probe.string();
                    break;
                }
            }
        }

        const std::string micm_config_path = this->config_file_.empty() ? this->config_dir_ : this->config_file_;
        Logger::info(
            *state, "GasChemProcess: resolved MICM configuration",
            {{"path", micm_config_path}, {"type", this->config_file_.empty() ? "directory (v0 CAMP)" : "file"}});

        // 2. Initialize MICM and State using musica library
        try {
            micm_instance_ = std::make_unique<musica::MICM>(micm_config_path, musica::RosenbrockStandardOrder);
            const std::size_t batch_cells = this->columns_per_batch_ * static_cast<std::size_t>(state->level_count());
            micm_state_ = std::make_unique<musica::State>(*micm_instance_, batch_cells);
            initialized_ = true;
            std::clog << "[INFO] GasChemProcess: initialized MICM successfully!" << std::endl;

            // Precompute species and rate-parameter index mappings once so run() performs
            // no per-cell string constructions or hash-map lookups.
            vector_size_ = micm_instance_->GetVectorSize();
            auto variable_map = micm_state_->GetVariableMap();
            n_micm_species_ = variable_map.size();
            species_mappings_.clear();
            species_mappings_.reserve(static_cast<std::size_t>(state->species_count()));

            for (int ispec = 0; ispec < state->species_count(); ++ispec) {
                const std::string& raw_name = state->chemistry().species_list[ispec].short_name;
                std::string upper_name = raw_name;
                std::transform(upper_name.begin(), upper_name.end(), upper_name.begin(),
                               [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
                std::string lower_name = raw_name;
                std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

                auto it = variable_map.find(raw_name);
                if (it == variable_map.end()) {
                    it = variable_map.find(upper_name);
                }
                if (it == variable_map.end()) {
                    it = variable_map.find(lower_name);
                }

                if (it != variable_map.end()) {
                    species_mappings_.push_back(SpeciesMapping{ispec, it->second, it->first});
                } else {
                    Logger::warn(*state, "Active species not found in MICM solver variable map",
                                 {{"species", upper_name}});
                }
            }

            auto rate_param_map = micm_state_->GetRateParameterMap();
            n_rate_params_ = rate_param_map.size();
            rate_param_mappings_.clear();
            rate_param_mappings_.reserve(n_rate_params_);
            for (const auto& [param_name, i_param] : rate_param_map) {
                if (param_name.rfind("PHOTO.", 0) == 0) {
                    rate_param_mappings_.push_back(RateParamMapping{i_param, RateParamKind::Photolysis,
                                                                    "photolysis_rate_" + param_name.substr(6)});
                } else if (param_name.rfind("LOSS.", 0) == 0) {
                    rate_param_mappings_.push_back(RateParamMapping{i_param, RateParamKind::Loss, {}});
                }
            }

            Logger::info(*state, "GasChemProcess: initialized MICM state",
                         {{"column_count", std::to_string(state->column_count())},
                          {"level_count", std::to_string(state->level_count())},
                          {"columns_per_batch", std::to_string(columns_per_batch_)},
                          {"micm_grid_cells", std::to_string(micm_state_->NumberOfGridCells())},
                          {"micm_species_count", std::to_string(micm_state_->NumberOfSpecies())}});
        } catch (const std::exception& e) {
            Logger::error(*state, "GasChemProcess: failed to initialize MICM", {{"error", e.what()}});
            initialized_ = false;
            throw std::runtime_error(std::string("GasChemProcess: failed to initialize MICM: ") + e.what());
        }
    }

    void GasChemProcess::run(std::shared_ptr<StateManager> state) {
        if (!initialized_) {
            Logger::warn(*state, "GasChemProcess: skipped run because solver is not initialized");
            return;
        }

        // 1. Sync device to host

        // Derived fields are retained across timesteps but invalidated when a
        // new NUOPC import generation begins.  Pointer existence therefore
        // does not mean that AIRDEN_DRY is current.
        if (state->meteorology().PMID && state->meteorology().T) {
            state->derive_airden_dry();
        }

        const auto import_generation = state->current_import_generation();
        if (!state->meteorology().T || !state->meteorology().PMID || !state->meteorology().AIRDEN_DRY ||
            !state->meteorology().T->is_current(import_generation) ||
            !state->meteorology().PMID->is_current(import_generation) ||
            !state->meteorology().AIRDEN_DRY->is_current(import_generation) || !state->chemistry().conc) {
            throw std::runtime_error(
                "GasChem requires current T, PMID, AIRDEN_DRY, and chemistry concentration fields");
        }

        auto temp = state->meteorology().T->host_view;
        auto pmid = state->meteorology().PMID->host_view;
        auto airden_dry = state->meteorology().AIRDEN_DRY->host_view;
        state->chemistry().conc->host_write();
        auto conc = state->chemistry().conc->host_view;

        auto& micm_conditions = micm_state_->GetConditions();
        auto& micm_concs = micm_state_->GetOrderedConcentrations();
        auto& micm_rate_params = micm_state_->GetOrderedRateParameters();

        int nc = state->column_count();
        int nl = state->level_count();

        // Dry air molecular weight in kg/mol (sourced from catchem::constants)
        constexpr double air_mw_kg = catchem::constants::AIR_MW * 1.0e-3;

        double tstep = state->clock().timestep;
        if (tstep <= 0.0) {
            Logger::error(*state, "Invalid timestep encountered", {{"timestep", std::to_string(tstep)}});
            throw std::runtime_error("GasChemProcess: timestep must be greater than zero.");
        }

        // Resolve photolysis diagnostic host pointers once per timestep rather than per cell.
        std::vector<const double*> photo_diag_ptrs(rate_param_mappings_.size(), nullptr);
        if (auto* diag_mgr = state->diagnostic_manager().get()) {
            for (std::size_t idx_map = 0; idx_map < rate_param_mappings_.size(); ++idx_map) {
                const auto& mapping = rate_param_mappings_[idx_map];
                if (mapping.kind == RateParamKind::Photolysis && diag_mgr->has_field(mapping.diag_name)) {
                    photo_diag_ptrs[idx_map] =
                        static_cast<const double*>(diag_mgr->get_host_pointer(mapping.diag_name));
                }
            }
        }

        const int batch_cols = static_cast<int>(columns_per_batch_);
        const int n_batches = (nc + batch_cols - 1) / batch_cols;

        for (int batch_idx = 0; batch_idx < n_batches; ++batch_idx) {
            const int col_start = batch_idx * batch_cols;
            const int col_count = std::min(batch_cols, nc - col_start);

            // 2. Map environmental variables and input concentrations to state for this batch.
            // Unused trailing columns in a short final batch (bcol >= col_count) are filled
            // with a copy of the last valid column in the batch so all cells in the fixed-size
            // MICM state hold valid physical state during Solve().
            for (int ilev = 0; ilev < nl; ++ilev) {
                for (int bcol = 0; bcol < batch_cols; ++bcol) {
                    const int icol = col_start + std::min(bcol, col_count - 1);
                    const std::size_t i_cell =
                        static_cast<std::size_t>(ilev) * columns_per_batch_ + static_cast<std::size_t>(bcol);

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

                    const std::size_t group_index = i_cell / vector_size_;
                    const std::size_t row_in_group = i_cell % vector_size_;
                    const std::size_t species_group_base = group_index * n_micm_species_;
                    const std::size_t param_group_base = group_index * n_rate_params_;

                    // Copy concentrations: ppmv -> mol/m3
                    for (const auto& sp : species_mappings_) {
                        const double ppmv_val = conc(icol, ilev, sp.catchem_idx);
                        if (!std::isfinite(ppmv_val) || ppmv_val < 0.0)
                            throw std::domain_error(
                                concentration_guard_message("GasChem received a negative or non-finite concentration",
                                                            ppmv_val, sp.name, icol, ilev, batch_idx));

                        const double conc_molar = ppmv_val * 1.0e-6 * air_density_mol;
                        const std::size_t idx = (species_group_base + sp.micm_idx) * vector_size_ + row_in_group;
                        micm_concs[idx] = conc_molar;
                    }

                    // 3. Dynamic Photolysis Mapping (PHOTO.<label> to photolysis_rate_<label>)
                    const int diag_idx = ilev * nc + icol;
                    for (std::size_t idx_map = 0; idx_map < rate_param_mappings_.size(); ++idx_map) {
                        const auto& param = rate_param_mappings_[idx_map];
                        const std::size_t idx = (param_group_base + param.micm_param_idx) * vector_size_ + row_in_group;
                        if (param.kind == RateParamKind::Photolysis) {
                            const double* diag_ptr = photo_diag_ptrs[idx_map];
                            micm_rate_params[idx] = diag_ptr ? diag_ptr[diag_idx] : 0.0;
                        } else {
                            micm_rate_params[idx] = 1.0;
                        }
                    }
                }
            }

            // 4. Run standard CPU solver for this batch
            auto solver_result = micm_instance_->Solve(micm_state_.get(), tstep);
            if (solver_result.state_ != micm::SolverState::Converged &&
                solver_result.state_ != micm::SolverState::AcceptingUnconvergedIntegration) {
                Logger::error(*state, "MICM Solver did not reach convergence!",
                              {{"batch", std::to_string(batch_idx)},
                               {"col_start", std::to_string(col_start)},
                               {"col_count", std::to_string(col_count)},
                               {"final_state", micm::SolverStateToString(solver_result.state_)}});
                throw std::runtime_error(
                    "GasChemProcess: MICM Solver failed to reach convergence or acceptable state in batch " +
                    std::to_string(batch_idx) + " (col_start=" + std::to_string(col_start) + ", col_count=" +
                    std::to_string(col_count) + "): " + micm::SolverStateToString(solver_result.state_));
            }

            // 5. Convert output concentrations back: mol/m3 -> ppmv (valid columns only)
            for (int ilev = 0; ilev < nl; ++ilev) {
                for (int bcol = 0; bcol < col_count; ++bcol) {
                    const int icol = col_start + bcol;
                    const std::size_t i_cell =
                        static_cast<std::size_t>(ilev) * columns_per_batch_ + static_cast<std::size_t>(bcol);
                    double air_density_mol = micm_conditions[i_cell].air_density_;

                    const std::size_t group_index = i_cell / vector_size_;
                    const std::size_t row_in_group = i_cell % vector_size_;
                    const std::size_t species_group_base = group_index * n_micm_species_;

                    for (const auto& sp : species_mappings_) {
                        const std::size_t idx = (species_group_base + sp.micm_idx) * vector_size_ + row_in_group;
                        const double conc_molar = micm_concs[idx];

                        const double ppmv_val = (conc_molar / air_density_mol) * 1.0e6;
                        if (!std::isfinite(ppmv_val) || ppmv_val < 0.0)
                            throw std::domain_error(concentration_guard_message(
                                "GasChem solver returned a negative or non-finite concentration", ppmv_val, sp.name,
                                icol, ilev, batch_idx));
                        conc(icol, ilev, sp.catchem_idx) = ppmv_val;
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
