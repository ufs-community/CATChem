#include "catchem_core.hpp"
#include "catchem_logger.hpp"
#include "catchem_process_registry.hpp"
#include <algorithm>
#include <cctype>
#include <iostream>
#include <limits>
#include <unordered_set>

namespace {

    // Resolve a (possibly relative) species path against the config file's
    // directory, so YAML-relative paths work regardless of the host's CWD.
    std::string resolve_against_config(const std::string& config_file, const std::string& path) {
        if (path.empty() || path.front() == '/')
            return path;
        auto slash = config_file.find_last_of('/');
        if (slash == std::string::npos)
            return path;
        return config_file.substr(0, slash + 1) + path;
    }

    void load_and_validate_configuration(catchem::ConfigManager& config, const std::string& config_file) {
        if (!config.data.species_filename.empty())
            config.load_species_file(resolve_against_config(config_file, config.data.species_filename));
        if (!config.data.simulation.emission_filename.empty())
            config.load_emission_mapping_file(
                resolve_against_config(config_file, config.data.simulation.emission_filename));
        config.validate_or_throw();
    }

    // Diagnostic-only: reports min/max of a fixed species watch-list after each
    // process, so a runaway species can be bisected to the introducing process
    // without re-running the whole coupled model multiple times.
    void log_watch_species_bounds(catchem::StateManager& state, const std::string& process_name, std::size_t step) {
        // Per-step, per-process, per-species diagnostics: keep them behind the
        // DEBUG threshold and skip the full-array sweep when it is off.
        if (!catchem::Logger::enabled(catchem::Logger::Level::Debug))
            return;
        static const std::vector<std::string> watch = {"so2",   "so4",   "dms",   "msa",   "bc1",   "bc2",
                                                       "oc1",   "oc2",   "dust1", "dust2", "dust3", "dust4",
                                                       "dust5", "seas1", "seas2", "seas3", "seas4", "seas5"};
        if (!state.chemistry().conc)
            return;
        const auto view = state.chemistry().conc->mdspan();
        const int nc = state.column_count();
        const int nl = state.level_count();
        for (std::size_t ispec = 0; ispec < state.chemistry().species_list.size(); ++ispec) {
            std::string lower = state.chemistry().species_list[ispec].short_name;
            std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
            if (std::find(watch.begin(), watch.end(), lower) == watch.end())
                continue;
            double lo = std::numeric_limits<double>::infinity();
            double hi = -std::numeric_limits<double>::infinity();
            double sum = 0.0;
            for (int icol = 0; icol < nc; ++icol)
                for (int ilev = 0; ilev < nl; ++ilev) {
                    const double v = view(icol, ilev, static_cast<int>(ispec));
                    lo = std::min(lo, v);
                    hi = std::max(hi, v);
                    sum += v;
                }
            catchem::Logger::debug(&state, "watch-species bounds after process",
                                   {{"step", std::to_string(step)},
                                    {"process", process_name},
                                    {"species", state.chemistry().species_list[ispec].short_name},
                                    {"min", std::to_string(lo)},
                                    {"max", std::to_string(hi)},
                                    {"mixing_ratio_sum", std::to_string(sum)}});
        }
    }

} // namespace

namespace catchem {

    // Member-wise assignment rather than designated initialisers: the host-only
    // build is C++17 (kept for UFS toolchain compatibility) and designated
    // initialisers are a C++20 feature that -pedantic rejects there.
    CoreCreateOptions CoreCreateOptions::direct_dimensions(int columns, int levels, int species) {
        CoreCreateOptions options;
        options.config_file.clear();
        options.columns = columns;
        options.levels = levels;
        options.species = species;
        options.use_configuration_grid = false;
        return options;
    }

    CoreCreateOptions CoreCreateOptions::configured(std::string config_file) {
        CoreCreateOptions options;
        options.config_file = std::move(config_file);
        options.columns = 1;
        options.levels = 1;
        options.species = 0;
        options.use_configuration_grid = true;
        return options;
    }

    CoreCreateOptions CoreCreateOptions::configured_with_host_grid(std::string config_file, int columns, int levels) {
        CoreCreateOptions options;
        options.config_file = std::move(config_file);
        options.columns = columns;
        options.levels = levels;
        options.species = 0;
        options.use_configuration_grid = false;
        return options;
    }

    Core::~Core() noexcept {
        try {
            shutdown();
        } catch (...) {
        }
    }

    void Core::shutdown() {
        std::lock_guard<std::mutex> lock(lifecycle_mutex_);
        if (shutdown_)
            return;
        std::exception_ptr first_error;
        for (std::size_t i = processes.size(); i-- > 0;) {
            if (i >= initialized_processes_.size() || !initialized_processes_[i])
                continue;
            try {
                processes[i]->finalize();
            } catch (...) {
                if (!first_error)
                    first_error = std::current_exception();
            }
            initialized_processes_[i] = false;
        }
        shutdown_ = true;
        runtime_lease_.release();
        if (first_error)
            std::rethrow_exception(first_error);
    }

    Core::Core(const CoreCreateOptions& options) {
        initialize(options);
    }

    Core::Core(int nc, int nl, int ns) : Core(CoreCreateOptions::direct_dimensions(nc, nl, ns)) {}

    Core::Core(const std::string& config_file) : Core(CoreCreateOptions::configured(config_file)) {}

    Core::Core(const std::string& config_file, int nc, int nl)
        : Core(CoreCreateOptions::configured_with_host_grid(config_file, nc, nl)) {}

    void Core::initialize(const CoreCreateOptions& options) {
        runtime_lease_ = RuntimeLease(RuntimeMode::CATChemOwned);
        try {
            config_mgr = std::make_shared<ConfigManager>();
            const bool configured = !options.config_file.empty();
            int nx = options.columns;
            int ny = 1;
            int nz = options.levels;
            int species_count = options.species;

            if (configured) {
                config_mgr->load_from_file(options.config_file);
                if (options.use_configuration_grid) {
                    nx = config_mgr->data.runtime.nx;
                    ny = config_mgr->data.runtime.ny;
                    nz = config_mgr->data.runtime.nz;
                } else {
                    // A coupled host owns its local grid; YAML still owns all other settings.
                    config_mgr->data.runtime.nx = nx;
                    config_mgr->data.runtime.ny = ny;
                    config_mgr->data.runtime.nz = nz;
                }
                load_and_validate_configuration(*config_mgr, options.config_file);
                species_count = static_cast<int>(config_mgr->data.species.size());
            } else {
                config_mgr->data.runtime.nx = nx;
                config_mgr->data.runtime.ny = ny;
                config_mgr->data.runtime.nz = nz;
            }

            grid_mgr = std::make_shared<GridManager>(nx, ny, nz);
            state_mgr = std::make_shared<StateManager>(nx * ny, nz, species_count);
            if (configured)
                state_mgr->set_validation_policy(config_mgr->data.physical_validation_policy);
            state_mgr->attach_config_manager(config_mgr);
            diag_mgr = std::make_shared<DiagnosticManager>();
            state_mgr->attach_diagnostic_manager(diag_mgr);
            if (configured) {
                state_mgr->set_configuration_path(config_mgr->config_file_path);
                state_mgr->chemistry().load_from_config_manager(*config_mgr);
                add_configured_processes();
            }
        } catch (...) {
            // A constructor that fails after process initialization does not
            // run Core's destructor.  Release initialized processes and the
            // runtime lease explicitly before propagating the configuration
            // or process-registration failure.
            try {
                shutdown();
            } catch (...) {
            }
            throw;
        }
    }

    std::shared_ptr<ConfigManager> Core::get_config_manager() {
        return config_mgr;
    }

    std::shared_ptr<GridManager> Core::get_grid_manager() {
        return grid_mgr;
    }

    std::shared_ptr<StateManager> Core::get_state_manager() {
        return state_mgr;
    }

    std::shared_ptr<DiagnosticManager> Core::get_diagnostic_manager() {
        return diag_mgr;
    }

    std::size_t Core::get_num_processes() const {
        return processes.size();
    }

    std::vector<std::string> Core::get_required_host_fields() const {
        std::vector<std::string> fields;
        std::unordered_set<std::string> seen;
        for (std::size_t index = 0; index < execution_plan_.size(); ++index) {
            for (const auto& access : execution_plan_.contract(index).fields) {
                if (!access.reads() || access.requirement != FieldRequirement::Required ||
                    canonicalize_field_identity(access.canonical_name) == "CONCENTRATION")
                    continue;
                const auto name = StateManager::canonical_field_name(access.canonical_name);
                if (seen.insert(name).second)
                    fields.push_back(name);
            }
        }
        return fields;
    }

    void Core::add_configured_processes() {
        auto& registry = ProcessRegistry::get_instance();
        for (const auto& process_name : config_mgr->data.active_processes) {
            const auto settings = config_mgr->data.processes.find(process_name);
            // run_phases defines schedule order; the process block controls
            // whether that scheduled entry is enabled for this runtime YAML.
            // A missing block defaults to disabled, so mechanisms can retain
            // a common schedule without hardcoding a process selection.
            if (settings == config_mgr->data.processes.end() || !settings->second.activate)
                continue;
            if (settings != config_mgr->data.processes.end())
                registry.validate_settings(process_name, settings->second);
            auto process = registry.create(process_name);
            process->init(state_mgr);
            try {
                add_process(process);
            } catch (...) {
                try {
                    process->finalize();
                } catch (...) {
                }
                throw;
            }
        }
    }

    void Core::add_process(std::shared_ptr<ProcessInterface> process) {
        if (!process)
            throw std::invalid_argument("Cannot add a null process");
        std::lock_guard<std::mutex> lock(lifecycle_mutex_);
        if (shutdown_)
            throw std::logic_error("Cannot add a process after shutdown");
        processes.push_back(process);
        initialized_processes_.push_back(true);
        execution_plan_.compile(processes, state_mgr->chemistry().mechanism.get());
        if (execution_plan_.validation().has_errors()) {
            const auto issue_summary = execution_plan_.validation().format();
            processes.pop_back();
            initialized_processes_.pop_back();
            execution_plan_.compile(processes, state_mgr->chemistry().mechanism.get());
            throw std::invalid_argument("Invalid process schedule: " + issue_summary);
        }
    }

    void Core::run_timestep(double dt) {
        // A process may hold host/device views while it runs.  Serialize the
        // whole step with add_process() and shutdown() so no process can be
        // finalized or schedule-recompiled underneath an active bridge call.
        std::lock_guard<std::mutex> lock(lifecycle_mutex_);
        if (shutdown_)
            throw std::logic_error("Cannot run a timestep after shutdown");
        if (tainted_) {
            if (state_mgr->current_import_generation() > last_outcome_.import_generation) {
                tainted_ = false;
            } else {
                throw std::runtime_error(
                    "Previous timestep partially updated state; a new import generation is required");
            }
        }
        last_outcome_ = {};
        last_outcome_.timestep = ++timestep_counter_;
        last_outcome_.duration = dt;
        last_outcome_.import_generation = state_mgr->current_import_generation();
        if (dt <= 0.0 || dt > 86400.0) {
            last_outcome_.status = TimestepStatus::ValidationFailed;
            last_outcome_.state = StateClassification::RequiresReimport;
            last_outcome_.cause = "invalid timestep duration";
            throw std::out_of_range(
                "Timestep dt must be positive and within a plausible physical daily limit (0 < dt <= 86400).");
        }

        try {
            state_mgr->validate_ready_for_execution();
        } catch (const std::exception& error) {
            last_outcome_.status = TimestepStatus::ValidationFailed;
            last_outcome_.state = StateClassification::RequiresReimport;
            last_outcome_.cause = error.what();
            throw;
        }
        last_outcome_.status = TimestepStatus::Running;
        const bool verbose = config_mgr && config_mgr->data.simulation.verbose_enabled;
        if (verbose) {
            Logger::debug(state_mgr.get(), "Core timestep begin",
                          {{"step", std::to_string(last_outcome_.timestep)},
                           {"dt_s", std::to_string(dt)},
                           {"processes", std::to_string(processes.size())},
                           {"import_generation", std::to_string(last_outcome_.import_generation)}});
        }
        diag_mgr->begin_timestep();

        for (std::size_t index = 0; index < processes.size(); ++index) {
            try {
                if (verbose) {
                    Logger::debug(state_mgr.get(), "Core process prepare",
                                  {{"step", std::to_string(last_outcome_.timestep)},
                                   {"index", std::to_string(index)},
                                   {"process", processes[index]->get_name()}});
                    Logger::debug(state_mgr.get(), "Core process input preparation",
                                  {{"step", std::to_string(last_outcome_.timestep)},
                                   {"index", std::to_string(index)},
                                   {"process", processes[index]->get_name()}});
                }
                processes[index]->prepare_inputs(state_mgr);
                execution_plan_.prepare(index, *state_mgr);
                if (verbose) {
                    Logger::debug(state_mgr.get(), "Core process run",
                                  {{"step", std::to_string(last_outcome_.timestep)},
                                   {"index", std::to_string(index)},
                                   {"process", processes[index]->get_name()}});
                }
                processes[index]->run(state_mgr);
                if (verbose) {
                    log_watch_species_bounds(*state_mgr, processes[index]->get_name(), last_outcome_.timestep);
                    Logger::debug(state_mgr.get(), "Core process bookkeeping",
                                  {{"step", std::to_string(last_outcome_.timestep)},
                                   {"index", std::to_string(index)},
                                   {"process", processes[index]->get_name()}});
                }
                execution_plan_.complete(index, *state_mgr);
                if (verbose) {
                    Logger::debug(state_mgr.get(), "Core process complete",
                                  {{"step", std::to_string(last_outcome_.timestep)},
                                   {"index", std::to_string(index)},
                                   {"process", processes[index]->get_name()}});
                }
            } catch (const std::exception& error) {
                last_outcome_.status = TimestepStatus::PartialUpdate;
                last_outcome_.process_name = processes[index]->get_name();
                last_outcome_.process_index = index;
                last_outcome_.cause = error.what();
                last_outcome_.state = StateClassification::RequiresReimport;
                tainted_ = true;
                diag_mgr->mark_generation_failed();
                Logger::error(state_mgr.get(), "Core process failed",
                              {{"step", std::to_string(last_outcome_.timestep)},
                               {"index", std::to_string(index)},
                               {"process", last_outcome_.process_name},
                               {"cause", last_outcome_.cause}});
                throw;
            } catch (...) {
                last_outcome_.status = TimestepStatus::PartialUpdate;
                last_outcome_.process_name = processes[index]->get_name();
                last_outcome_.process_index = index;
                last_outcome_.cause = "unknown non-standard exception";
                last_outcome_.state = StateClassification::RequiresReinitialize;
                tainted_ = true;
                diag_mgr->mark_generation_failed();
                Logger::error(state_mgr.get(), "Core process failed",
                              {{"step", std::to_string(last_outcome_.timestep)},
                               {"index", std::to_string(index)},
                               {"process", last_outcome_.process_name},
                               {"cause", last_outcome_.cause}});
                throw;
            }
        }

        // Sync diagnostics
        if (verbose)
            Logger::debug(state_mgr.get(), "Core diagnostics sync", {{"step", std::to_string(last_outcome_.timestep)}});
        diag_mgr->sync_to_host();
        last_outcome_.status = TimestepStatus::Succeeded;
        last_outcome_.state = StateClassification::Reusable;
        if (verbose)
            Logger::debug(state_mgr.get(), "Core timestep complete",
                          {{"step", std::to_string(last_outcome_.timestep)}});
    }

    void Core::run_timestep() {
        double dt = config_mgr->data.runtime.dt;
        run_timestep(dt);
    }

} // namespace catchem
