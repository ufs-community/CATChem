// src/process/gaschem/catchem_process_gaschem.cpp
#include "catchem_process_gaschem.hpp"
#include "catchem_constants.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_logger.hpp"
#include "catchem_process_registry.hpp"
#include <algorithm>
#include <array>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace catchem {

    namespace {

        // Scientific notation, because fixed six-decimal formatting prints every
        // tiny negative as "-0.000000" and hides its magnitude.
        std::string format_scientific(double value) {
            std::ostringstream text;
            text << std::scientific << std::setprecision(6) << value;
            return text.str();
        }

        // Builds the message for the concentration guards in GasChemProcess::run.
        // Column, level, and batch are zero-based; column is the index on this
        // PET, not within the batch.
        std::string concentration_guard_message(const std::string& what, double value, const std::string& species,
                                                int column, int level, int batch) {
            std::ostringstream message;
            message << what << " (" << format_scientific(value) << ") for " << species << " at column=" << column
                    << ", level=" << level << ", batch=" << batch << " (zero-based indices)";
            return message.str();
        }

        std::string join_names(const std::vector<std::string>& names) {
            std::string joined;
            for (const auto& name : names) {
                if (!joined.empty()) {
                    joined += ',';
                }
                joined += name;
            }
            return joined;
        }

        // The most negative concentration floored during one run() call.
        struct WorstNegative {
            std::size_t count = 0;
            double value = 0.0;
            const std::string* species = nullptr;
            int column = -1;
            int level = -1;

            void record(double candidate, const std::string& name, int icol, int ilev) {
                ++count;
                if (candidate < value) {
                    value = candidate;
                    species = &name;
                    column = icol;
                    level = ilev;
                }
            }
            std::string species_name() const { return species ? *species : std::string("none"); }
        };

    } // namespace

    ProcessContract GasChemProcess::get_contract() const {
        return make_contract(get_name(),
                             {host_field_3d("T", "K"), host_field_3d("PMID", "Pa"),
                              host_field_3d("AIRDEN_DRY", "kg/m3", FieldRequirement::Optional), host_concentration()});
    }

    GasChemProcess::GasChemProcess() = default;
    GasChemProcess::~GasChemProcess() = default;

    GasChemSettings GasChemProcess::parse_settings(const ProcessConfig& config) {
        static constexpr std::array<std::string_view, 4> accepted = {"columns_per_batch", "config_dir", "config_file",
                                                                     "negative_concentration_policy"};

        // Every gaschem option is a flat key, so a nested block is unknown too.
        std::vector<std::string> unknown;
        for (const auto& key : config.flat_option_keys()) {
            if (std::find(accepted.begin(), accepted.end(), key) == accepted.end()) {
                unknown.push_back(key);
            }
        }
        for (const auto& path : config.option_paths()) {
            unknown.push_back(path);
        }
        // gaschem has no scheme to select. The configuration layer treats the
        // scheme selectors as framework keys for every process, so neither
        // listing above reports them and they are checked here.
        if (!config.scheme.empty()) {
            unknown.emplace_back("scheme");
        }
        for (const char* selector : {"gas_scheme", "aero_scheme"}) {
            if (!config.get_string(selector).empty()) {
                unknown.emplace_back(selector);
            }
        }
        if (!unknown.empty()) {
            std::string accepted_names;
            for (const auto& name : accepted) {
                if (!accepted_names.empty()) {
                    accepted_names += ", ";
                }
                accepted_names += name;
            }
            std::sort(unknown.begin(), unknown.end());
            std::string unknown_names;
            for (const auto& name : unknown) {
                if (!unknown_names.empty()) {
                    unknown_names += ", ";
                }
                unknown_names += name;
            }
            throw std::invalid_argument("Unknown option(s) for process 'gaschem': " + unknown_names +
                                        ". Accepted options: " + accepted_names + ".");
        }

        GasChemSettings settings;
        settings.config_file = config.get_string("config_file");
        settings.config_dir = config.get_string("config_dir");

        const int columns_per_batch = config.get_int("columns_per_batch", 1);
        if (columns_per_batch < 1) {
            throw std::invalid_argument("GasChemProcess: columns_per_batch must be >= 1, got " +
                                        std::to_string(columns_per_batch));
        }
        settings.columns_per_batch = static_cast<std::size_t>(columns_per_batch);

        // Negative concentrations are rejected unless the run opts in to flooring
        // them. This key is separate from physical_validation/policy, which governs
        // meteorology, so a run can floor tracer noise from the host and still stop
        // on a bad temperature.
        const std::string negative_policy = config.get_string("negative_concentration_policy", "reject");
        if (negative_policy == "reject") {
            settings.negative_concentration_policy = NegativeConcentrationPolicy::Reject;
        } else if (negative_policy == "clamp") {
            settings.negative_concentration_policy = NegativeConcentrationPolicy::Clamp;
        } else {
            throw std::invalid_argument("GasChemProcess: negative_concentration_policy must be reject or clamp, got " +
                                        negative_policy);
        }
        return settings;
    }

    void GasChemProcess::init(std::shared_ptr<StateManager> state) {
        Logger::debug(*state, "GasChemProcess::init started");

        // 1. Read the options. The core has already run parse_settings as the
        // settings validator when it built this process from the runtime YAML;
        // calling it here keeps the same checks in force when a test or driver
        // initializes the process directly. A legacy top-level "process:" block
        // is parsed into the same ProcessConfig, so it needs no separate lookup.
        // An explicit config_file (single-file v1 mechanism) takes precedence
        // over config_dir.
        GasChemSettings settings;
        if (const auto config = state->config_manager()) {
            const auto configured = config->data.processes.find(get_name());
            if (configured != config->data.processes.end()) {
                settings = parse_settings(configured->second);
            }
        }
        this->config_file_ = settings.config_file;
        this->config_dir_ = settings.config_dir;
        this->columns_per_batch_ = settings.columns_per_batch;
        this->negative_concentration_policy_ = settings.negative_concentration_policy;
        this->negative_inputs_clamped_ = 0;
        this->negative_outputs_clamped_ = 0;
        this->min_negative_input_ = 0.0;
        this->min_negative_output_ = 0.0;
        this->input_clamp_warned_ = false;
        this->output_clamp_warned_ = false;

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
            std::clog << "[INFO] GasChemProcess: initialized MICM successfully!" << std::endl;
        } catch (const std::exception& e) {
            Logger::error(*state, "GasChemProcess: failed to initialize MICM", {{"error", e.what()}});
            initialized_ = false;
            throw std::runtime_error(std::string("GasChemProcess: failed to initialize MICM: ") + e.what());
        }

        // 3. Resolve every MICM variable to a CATChem species through the active
        // mechanism definition, which matches names case-insensitively and honours
        // species aliases. Doing this once here keeps string construction and
        // hash-map lookups out of the cell loops in run().
        vector_size_ = micm_instance_->GetVectorSize();
        const auto variable_map = micm_state_->GetVariableMap();
        n_micm_species_ = variable_map.size();
        species_mappings_.clear();
        species_mappings_.reserve(n_micm_species_);
        std::vector<std::string> unmapped_micm_species;

        const auto& mechanism = state->chemistry().mechanism;
        for (const auto& [micm_name, micm_idx] : variable_map) {
            if (!mechanism || !mechanism->contains(micm_name)) {
                unmapped_micm_species.push_back(micm_name);
                continue;
            }
            species_mappings_.push_back(
                SpeciesMapping{static_cast<int>(mechanism->index_of(micm_name)), micm_idx, micm_name});
        }

        // MICM integrates every species in its mechanism, so a species the host does
        // not supply would start from whatever an earlier solve left in the reused
        // state. Require complete coverage instead of choosing a value for it.
        if (!unmapped_micm_species.empty()) {
            std::sort(unmapped_micm_species.begin(), unmapped_micm_species.end());
            throw std::invalid_argument(
                "GasChemProcess: MICM mechanism species missing from the CATChem species list: " +
                join_names(unmapped_micm_species) +
                ". Add them to the species file; use is_advected: false for species the host should not transport.");
        }

        // The variable map is unordered. Sort by CATChem species index so run()
        // visits species in a fixed order, and so two MICM variables that resolve
        // to the same CATChem species end up next to each other.
        std::sort(species_mappings_.begin(), species_mappings_.end(),
                  [](const SpeciesMapping& a, const SpeciesMapping& b) {
                      return a.catchem_idx != b.catchem_idx ? a.catchem_idx < b.catchem_idx : a.name < b.name;
                  });
        for (std::size_t i = 1; i < species_mappings_.size(); ++i) {
            const auto& previous = species_mappings_[i - 1];
            const auto& current = species_mappings_[i];
            if (previous.catchem_idx == current.catchem_idx) {
                const auto& target = state->chemistry().species_list[static_cast<std::size_t>(current.catchem_idx)];
                throw std::invalid_argument("GasChemProcess: MICM species " + previous.name + " and " + current.name +
                                            " both resolve to CATChem species " + target.short_name);
            }
        }

        // Gas species the host carries that the mechanism does not contain are left
        // untouched by run(). Report them once, in one line.
        std::vector<char> in_mechanism(static_cast<std::size_t>(state->species_count()), 0);
        for (const auto& sp : species_mappings_) {
            in_mechanism[static_cast<std::size_t>(sp.catchem_idx)] = 1;
        }
        std::vector<std::string> gas_outside_mechanism;
        for (const int ispec : state->chemistry().indices(SpeciesCategory::Gas)) {
            if (!in_mechanism[static_cast<std::size_t>(ispec)]) {
                gas_outside_mechanism.push_back(
                    state->chemistry().species_list[static_cast<std::size_t>(ispec)].short_name);
            }
        }
        if (!gas_outside_mechanism.empty()) {
            Logger::warn(*state, "GasChemProcess: CATChem gas species not in the MICM mechanism",
                         {{"count", std::to_string(gas_outside_mechanism.size())},
                          {"species", join_names(gas_outside_mechanism)}});
        }

        const auto rate_param_map = micm_state_->GetRateParameterMap();
        n_rate_params_ = rate_param_map.size();
        rate_param_mappings_.clear();
        rate_param_mappings_.reserve(n_rate_params_);
        for (const auto& [param_name, i_param] : rate_param_map) {
            if (param_name.rfind("PHOTO.", 0) == 0) {
                rate_param_mappings_.push_back(
                    RateParamMapping{i_param, RateParamKind::Photolysis, "photolysis_rate_" + param_name.substr(6)});
            } else if (param_name.rfind("LOSS.", 0) == 0) {
                rate_param_mappings_.push_back(RateParamMapping{i_param, RateParamKind::Loss, {}});
            }
        }

        initialized_ = true;
        Logger::info(*state, "GasChemProcess: initialized MICM state",
                     {{"column_count", std::to_string(state->column_count())},
                      {"level_count", std::to_string(state->level_count())},
                      {"columns_per_batch", std::to_string(columns_per_batch_)},
                      {"negative_concentration_policy",
                       negative_concentration_policy_ == NegativeConcentrationPolicy::Clamp ? "clamp" : "reject"},
                      {"micm_grid_cells", std::to_string(micm_state_->NumberOfGridCells())},
                      {"micm_species_count", std::to_string(micm_state_->NumberOfSpecies())},
                      {"mapped_species_count", std::to_string(species_mappings_.size())}});
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

        const bool reject_negative = negative_concentration_policy_ == NegativeConcentrationPolicy::Reject;
        WorstNegative floored_inputs;
        WorstNegative floored_outputs;

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

                    // Copy concentrations: ppmv -> mol/m3. init() guarantees every MICM
                    // species is mapped, so this overwrites the whole cell and nothing
                    // from the previous batch survives in the reused state.
                    for (const auto& sp : species_mappings_) {
                        double ppmv_val = conc(icol, ilev, sp.catchem_idx);
                        if (!std::isfinite(ppmv_val) || (ppmv_val < 0.0 && reject_negative))
                            throw std::domain_error(
                                concentration_guard_message("GasChem received a negative or non-finite concentration",
                                                            ppmv_val, sp.name, icol, ilev, batch_idx));
                        if (ppmv_val < 0.0) {
                            // Padding cells repeat the last valid column; count a value only in a real cell.
                            if (bcol < col_count) {
                                floored_inputs.record(ppmv_val, sp.name, icol, ilev);
                            }
                            ppmv_val = 0.0;
                        }

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
            adjust_solver_output_for_testing(micm_concs, batch_idx);

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

                        double ppmv_val = (conc_molar / air_density_mol) * 1.0e6;
                        if (!std::isfinite(ppmv_val) || (ppmv_val < 0.0 && reject_negative))
                            throw std::domain_error(concentration_guard_message(
                                "GasChem solver returned a negative or non-finite concentration", ppmv_val, sp.name,
                                icol, ilev, batch_idx));
                        if (ppmv_val < 0.0) {
                            floored_outputs.record(ppmv_val, sp.name, icol, ilev);
                            ppmv_val = 0.0;
                        }
                        conc(icol, ilev, sp.catchem_idx) = ppmv_val;
                    }
                }
            }
        }

        // 6. Report floored values. This PET warns once for inputs and once for
        // solver outputs, the first time each is floored, and logs at debug after
        // that, so a long run does not log one line per PET per timestep. A
        // negative from the solver is a different problem from noise in the host's
        // tracers, so an earlier input warning does not silence it.
        if (floored_inputs.count + floored_outputs.count > 0) {
            negative_inputs_clamped_ += floored_inputs.count;
            negative_outputs_clamped_ += floored_outputs.count;
            min_negative_input_ = std::min(min_negative_input_, floored_inputs.value);
            min_negative_output_ = std::min(min_negative_output_, floored_outputs.value);

            const bool first_input = floored_inputs.count > 0 && !input_clamp_warned_;
            const bool first_output = floored_outputs.count > 0 && !output_clamp_warned_;
            input_clamp_warned_ = input_clamp_warned_ || first_input;
            output_clamp_warned_ = output_clamp_warned_ || first_output;
            const bool warning_due = first_input || first_output;

            // Format the context only when a line will be written.
            if (warning_due || Logger::enabled(Logger::Level::Debug)) {
                const std::string_view first_occurrence = first_input
                                                              ? (first_output ? "inputs,solver_outputs" : "inputs")
                                                              : (first_output ? "solver_outputs" : "none");
                const std::string inputs_floored = std::to_string(floored_inputs.count);
                const std::string min_input = format_scientific(floored_inputs.value);
                const std::string min_input_species = floored_inputs.species_name();
                const std::string min_input_column = std::to_string(floored_inputs.column);
                const std::string min_input_level = std::to_string(floored_inputs.level);
                const std::string outputs_floored = std::to_string(floored_outputs.count);
                const std::string min_output = format_scientific(floored_outputs.value);
                const std::string min_output_species = floored_outputs.species_name();
                const std::string min_output_column = std::to_string(floored_outputs.column);
                const std::string min_output_level = std::to_string(floored_outputs.level);
                const Logger::ContextList context = {{"first_occurrence", first_occurrence},
                                                     {"inputs_floored", inputs_floored},
                                                     {"min_input", min_input},
                                                     {"min_input_species", min_input_species},
                                                     {"min_input_column", min_input_column},
                                                     {"min_input_level", min_input_level},
                                                     {"outputs_floored", outputs_floored},
                                                     {"min_output", min_output},
                                                     {"min_output_species", min_output_species},
                                                     {"min_output_column", min_output_column},
                                                     {"min_output_level", min_output_level}};
                if (warning_due) {
                    Logger::warn(*state, "GasChemProcess: floored negative concentrations to zero", context);
                } else {
                    Logger::debug(*state, "GasChemProcess: floored negative concentrations to zero", context);
                }
            }
        }

        // 7. Sync back to device
        if (state->chemistry().conc)
            state->chemistry().conc->mark_host_modified();
    }

    double GasChemProcess::solver_relative_tolerance() const {
        return micm_instance_ ? micm_instance_->GetRosenbrockSolverParameters().relative_tolerance : 0.0;
    }

    void GasChemProcess::finalize() {
        if (negative_inputs_clamped_ + negative_outputs_clamped_ > 0) {
            Logger::info("GasChemProcess: negative concentrations floored to zero since init",
                         {{"inputs_floored", std::to_string(negative_inputs_clamped_)},
                          {"min_input", format_scientific(min_negative_input_)},
                          {"outputs_floored", std::to_string(negative_outputs_clamped_)},
                          {"min_output", format_scientific(min_negative_output_)}});
        }
    }

} // namespace catchem

extern "C" {
void catchem_register_gaschem_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        std::string(catchem::ProcessNames::GasChem), []() { return std::make_shared<catchem::GasChemProcess>(); }, {},
        [](const catchem::ProcessConfig& config) { (void)catchem::GasChemProcess::parse_settings(config); });
}
}
