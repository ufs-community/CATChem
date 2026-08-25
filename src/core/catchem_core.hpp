/**
 * @file catchem_core.hpp
 * @brief Central Orchestration Engine for the CATChem modern C++ framework.
 */

#pragma once
#include "catchem_config_manager.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_grid_manager.hpp"
#include "catchem_process_interface.hpp"
#include "catchem_state_manager.hpp"
#include <memory>
#include <string>
#include <vector>

namespace catchem {

    /**
     * @class Core
     * @brief High-performance orchestration class managing the atmospheric physics lifecycle.
     *
     * The Core class coordinates grid layout, meteorological states, and schedules
     * chemical processes. It acts as the single source of truth for the physical simulation.
     */
    class Core {
    private:
        std::shared_ptr<ConfigManager> config_mgr;                ///< Shared configuration manager.
        std::shared_ptr<GridManager> grid_mgr;                    ///< Shared grid details.
        std::shared_ptr<StateManager> state_mgr;                  ///< Reference-counted memory state.
        std::shared_ptr<DiagnosticManager> diag_mgr;              ///< Global diagnostics registry.
        std::vector<std::shared_ptr<ProcessInterface>> processes; ///< Scheduled physics processes.
    public:
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

        /** @brief Get the configuration manager. */
        std::shared_ptr<ConfigManager> get_config_manager();

        /** @brief Get the grid manager. */
        std::shared_ptr<GridManager> get_grid_manager();

        /** @brief Get the state manager. */
        std::shared_ptr<StateManager> get_state_manager();

        /** @brief Get the diagnostic manager. */
        std::shared_ptr<DiagnosticManager> get_diagnostic_manager();

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
