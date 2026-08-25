#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_property_harness.hpp"
#include "catchem_state_manager.hpp"
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

extern "C" {
void catchem_register_seasalt_cpp();
void catchem_register_drydep_cpp();
void catchem_register_wetdep_cpp();
void catchem_register_settling_cpp();
void catchem_register_so4chem_cpp();
void catchem_register_dust_cpp();
void catchem_register_carbchem_cpp();
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "=== RUNNING ATMOSPHERIC PROPERTY TESTS ===" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        // Verify Trace ID Generation
        {
            auto test_state = std::make_shared<catchem::StateManager>(4, 10, 50);
            assert(test_state->trace_id.length() == 8);
            assert(!test_state->trace_id.empty());
        }

        // Register All C++ Process Handlers
        catchem_register_seasalt_cpp();
        catchem_register_drydep_cpp();
        catchem_register_wetdep_cpp();
        catchem_register_settling_cpp();
        catchem_register_so4chem_cpp();
        catchem_register_dust_cpp();
        catchem_register_carbchem_cpp();

        catchem::test::PropertyHarnessConfig harness_config;
        harness_config.n_cols = 12;
        harness_config.n_levels = 8;
        harness_config.n_species = 22;
        harness_config.iterations_per_scenario = 20;

        // Create Core Orchestration Layer
        void* core_ptr = catchem_core_create(harness_config.n_cols, harness_config.n_levels, harness_config.n_species);
        auto* core = static_cast<catchem::Core*>(core_ptr);
        auto state = core->get_state_manager();

        // Load species metadata configuration
        std::string config_path = "";
        std::vector<std::string> candidates = {"tests/Configs/Default/CATChem_species.yml",
                                               "../tests/Configs/Default/CATChem_species.yml",
                                               "../../tests/Configs/Default/CATChem_species.yml",
                                               "Configs/Default/CATChem_species.yml",
                                               "tests/CATChem_species.yml",
                                               "../tests/CATChem_species.yml"};
        for (const auto& path : candidates) {
            std::ifstream f(path);
            if (f.good()) {
                config_path = path;
                break;
            }
        }
        if (config_path.empty()) {
            std::cerr << "ERROR: Could not find CATChem_species.yml inside test_catchem_properties.cpp\n";
            std::exit(1);
        }
        state->load_species_config(config_path);

        // Schedule All Registered Processes dynamically to test synchronized coupled execution
        auto& registry = catchem::ProcessRegistry::get_instance();
        std::vector<std::string> process_names = {"settling", "drydep", "seasalt", "wetdep",
                                                  "so4chem",  "dust",   "carbchem"};

        for (const auto& name : process_names) {
            auto proc = registry.create(name);
            proc->init(state);
            core->add_process(proc);
        }

        // Execute Property Test Harness across all atmospheric scenarios
        catchem::test::PropertyTestHarness harness(harness_config);
        harness.run_full_suite(core);

        // Finalize lifecycle
        catchem_core_destroy(core_ptr);
    }
    Kokkos::finalize();
    return 0;
}
