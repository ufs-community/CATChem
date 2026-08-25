#include "catchem_config_manager.hpp"
#include <iostream>

namespace catchem {

    void ConfigManager::load_from_file(const std::string& filename) {
        config_file_path = filename;
        try {
            YAML::Node config = YAML::LoadFile(filename);
            root_node = config;
            if (config["simulation"]) {
                auto sim = config["simulation"];
                if (sim["nx"]) {
                    data.runtime.nx = sim["nx"].as<int>();
                }
                if (sim["ny"]) {
                    data.runtime.ny = sim["ny"].as<int>();
                }
                if (sim["nz"]) {
                    data.runtime.nz = sim["nz"].as<int>();
                }
                if (sim["timestep"]) {
                    data.runtime.dt = sim["timestep"].as<double>();
                } else if (sim["dt"]) {
                    data.runtime.dt = sim["dt"].as<double>();
                }
                if (sim["nsteps"]) {
                    data.runtime.nsteps = sim["nsteps"].as<int>();
                }
            }
            is_loaded = true;
        } catch (const std::exception& e) {
            std::cerr << "Error loading configuration file " << filename << ": " << e.what() << std::endl;
            is_loaded = false;
            throw;
        }
    }

} // namespace catchem
