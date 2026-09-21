/**
 * @file catchem_core.hpp
 * @brief Central Orchestration Engine for the CATChem modern C++ framework.
 */

#pragma once
#include "catchem_config_manager.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_execution_plan.hpp"
#include "catchem_grid_manager.hpp"
#include "catchem_process_interface.hpp"
#include "catchem_runtime_lease.hpp"
#include "catchem_state_manager.hpp"
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace catchem {

    /** Inputs for the single staged Core construction pipeline. */
    struct CoreCreateOptions {
        std::string config_file;
        int columns = 1;
        int levels = 1;
        int species = 0;
        bool use_configuration_grid = false;

        static CoreCreateOptions direct_dimensions(int columns, int levels, int species);
        static CoreCreateOptions configured(std::string config_file);
        static CoreCreateOptions configured_with_host_grid(std::string config_file, int columns, int levels);
    };

    /**
     * @class Core
     * @brief High-performance orchestration class managing the atmospheric physics lifecycle.
     *
     * The Core class coordinates grid layout, meteorological states, and schedules
     * chemical processes. It acts as the single source of truth for the physical simulation.
     */
    class Core {
    private:
        RuntimeLease runtime_lease_;                              ///< Process-global execution-runtime lease.
        std::shared_ptr<ConfigManager> config_mgr;                ///< Shared configuration manager.
        std::shared_ptr<GridManager> grid_mgr;                    ///< Shared grid details.
        std::shared_ptr<StateManager> state_mgr;                  ///< Reference-counted memory state.
        std::shared_ptr<DiagnosticManager> diag_mgr;              ///< Global diagnostics registry.
        std::vector<std::shared_ptr<ProcessInterface>> processes; ///< Scheduled physics processes.
        TimestepOutcome last_outcome_;
        std::size_t timestep_counter_ = 0;
        bool tainted_ = false;
        ExecutionPlan execution_plan_;
        mutable std::mutex lifecycle_mutex_;
        std::vector<bool> initialized_processes_;
        bool shutdown_ = false;

        /** @brief Creates and initializes the processes requested by configuration. */
        void add_configured_processes();

        /** @brief Performs the common ordered construction stages. */
        void initialize(const CoreCreateOptions& options);

    public:
        explicit Core(const CoreCreateOptions& options);

        /**
         * @brief Constructs the Core with dimensions.
         * @param nc Number of contiguous columns.
         * @param nl Number of vertical levels.
         * @param ns Number of chemical species.
         */
        Core(int nc, int nl, int ns);

        /**
         * @brief Constructs the Core and loads configuration parameters from a YAML file.
         * @param config_file Path to the YAML configuration file.
         * @throws std::runtime_error If the file cannot be read or parsed.
         */
        Core(const std::string& config_file);

        /**
         * @brief Constructs the Core from a YAML file with host-supplied grid dimensions.
         *
         * Configuration (species, processes, runtime options) comes from the file;
         * the grid is sized by the host — required under domain decomposition
         * (e.g. UFS per-rank tiles), where the YAML grid section does not apply.
         * @param config_file Path to the YAML configuration file.
         * @param nc Number of contiguous columns (host-local).
         * @param nl Number of vertical levels.
         * @throws std::runtime_error If the file cannot be read or parsed.
         */
        Core(const std::string& config_file, int nc, int nl);
        ~Core() noexcept;

        /** Finalizes initialized processes exactly once, in reverse order. */
        void shutdown();

        /** @brief Get the configuration manager. */
        std::shared_ptr<ConfigManager> get_config_manager();

        /** @brief Get the grid manager. */
        std::shared_ptr<GridManager> get_grid_manager();

        /** @brief Get the state manager. */
        std::shared_ptr<StateManager> get_state_manager();

        /** @brief Get the diagnostic manager. */
        std::shared_ptr<DiagnosticManager> get_diagnostic_manager();

        /** @brief Get the number of scheduled physics processes. */
        std::size_t get_num_processes() const;

        /**
         * @brief Return the unique required host-input field names declared by
         *        the active process contracts.
         *
         * Derived state and driver-owned static inputs are intentionally not
         * distinguished here; integrations use their own field contracts to
         * decide which of these requirements are supplied by their imports.
         */
        std::vector<std::string> get_required_host_fields() const;
        const TimestepOutcome& get_timestep_outcome() const noexcept { return last_outcome_; }

        /**
         * @brief Registers a new physics process into the simulation schedule.
         * @param process Shared pointer to the process interface.
         */
        void add_process(std::shared_ptr<ProcessInterface> process);

        /**
         * @brief Executes a single timestep.
         * @param dt Timestep duration in seconds.
         * @throws std::out_of_range If dt is non-positive.
         */
        void run_timestep(double dt);

        /**
         * @brief Executes a single timestep using the default step size from configuration.
         */
        void run_timestep();
    };

} // namespace catchem
