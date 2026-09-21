#include "catchem_config_manager.hpp"
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

int main() {
    catchem::ConfigManager config;
    config.load_from_file("Configs/Default/CATChem_new_config.yml");
    config.load_species_file("Configs/Default/CATChem_species.yml");
    config.validate_or_throw();

    const std::vector<std::string> expected_schedule = {"seasalt", "dust",    "carbchem", "settling",
                                                        "drydep",  "so4chem", "wetdep"};
    assert(config.data.active_processes == expected_schedule);
    assert(config.data.grid.number_of_levels == 64);
    assert(config.data.timesteps.transport_timestep_in_s == 10);
    assert(config.data.timesteps.chemistry_timestep_in_s == 60);

    // This fixture intentionally exercises the complete Default process
    // schedule. Keep its activation state and science selections stable: a
    // legacy-vs-C++ numerical parity runner consumes this exact contract.
    for (const auto& name : expected_schedule) {
        assert(config.is_process_active(name));
    }
    assert(config.data.processes.at("seasalt").scheme == "geos12");
    assert(config.data.processes.at("dust").scheme == "fengsha");
    assert(config.data.processes.at("settling").scheme == "gocart");
    assert(config.data.processes.at("wetdep").scheme == "jacob");
    assert(config.data.processes.at("so4chem").scheme == "gocart");
    assert(config.data.processes.at("carbchem").scheme == "gocart");
    assert(config.data.processes.at("drydep").get_string("gas_scheme") == "wesely");
    assert(config.data.processes.at("drydep").get_string("aero_scheme") == "gocart");

    assert(config.is_process_active("extemis"));
    // Diagnostics output is disabled in the Default configuration (upstream
    // ddf07c50, "Disable diagnostics output in configuration file"); the
    // rest of the output block (compression level etc.) stays as configured.
    assert(!config.data.diagnostics.output.enabled);
    assert(config.data.diagnostics.output.compress_lev == 2);
    assert(config.data.species.size() == 22);
    std::cout << "PASS: Default runtime configuration contract\n";
    return 0;
}
