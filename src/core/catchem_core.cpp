#include "catchem_core.hpp"
#include <iostream>

namespace catchem {

    Core::Core(int nc, int nl, int ns) {
        config_mgr = std::make_shared<ConfigManager>();
        config_mgr->data.runtime.nx = nc;
        config_mgr->data.runtime.ny = 1;
        config_mgr->data.runtime.nz = nl;

        grid_mgr = std::make_shared<GridManager>(nc, 1, nl);
        state_mgr = std::make_shared<StateManager>(nc, nl, ns);
        state_mgr->config_mgr = config_mgr;
        diag_mgr = std::make_shared<DiagnosticManager>();
        state_mgr->diag_mgr = diag_mgr;
    }

    Core::Core(const std::string& config_file) {
        config_mgr = std::make_shared<ConfigManager>();
        config_mgr->load_from_file(config_file);

        int nx = config_mgr->data.runtime.nx;
        int ny = config_mgr->data.runtime.ny;
        int nz = config_mgr->data.runtime.nz;

        grid_mgr = std::make_shared<GridManager>(nx, ny, nz);
        state_mgr = std::make_shared<StateManager>(nx * ny, nz, 50);
        state_mgr->config_mgr = config_mgr;
        state_mgr->config_file_path = config_mgr->config_file_path;
        diag_mgr = std::make_shared<DiagnosticManager>();
        state_mgr->diag_mgr = diag_mgr;
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

    void Core::add_process(std::shared_ptr<ProcessInterface> process) {
        processes.push_back(process);
    }

    void Core::run_timestep(double dt) {
        if (dt <= 0.0 || dt > 86400.0) {
            throw std::out_of_range(
                "Timestep dt must be positive and within a plausible physical daily limit (0 < dt <= 86400).");
        }

        // Sync shared boundary arrays to active execution spaces
        state_mgr->sync_to_device();

        for (auto& process : processes) {
            process->run(state_mgr);
        }

        // Sync execution outputs back to Fortran-accessible memory
        state_mgr->sync_to_host();

        // Sync diagnostics
        diag_mgr->sync_to_host();
    }

    void Core::run_timestep() {
        double dt = config_mgr->data.runtime.dt;
        run_timestep(dt);
    }

} // namespace catchem
