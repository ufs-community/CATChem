#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_state_manager.hpp"
#include "catchem_test_config.hpp"
#include <Kokkos_Core.hpp>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

extern "C" {
void catchem_register_photolysis_cpp();
}

// Simple helper to check if a file exists using C++17 std::filesystem
bool file_exists(const std::string& name) {
    return std::filesystem::exists(name);
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "RUNNING TEST: MUSICA TUV-x Photolysis Process Integration" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        // 1. Verify dynamic registration
        catchem_register_photolysis_cpp();
        assert(catchem::ProcessRegistry::get_instance().has_process("photolysis"));
        std::cout << "SUCCESS: Photolysis process successfully registered in registry." << std::endl;

        // 2. Set up core with a single column and vertical grid matching config (3 levels -> 4 edges)
        int n_cols = 1;
        int n_levels = 3;
        int n_species = 5;

        auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
        auto state = core->get_state_manager();

        // Configure the simulation time to be noon (12:00 PM) during summer to ensure positive solar radiation
        state->time.year = 2026;
        state->time.month = 7;
        state->time.day = 13;
        state->time.hour = 12; // Noon
        state->time.minute = 0;
        state->time.second = 0;
        state->time.calculate_derived_fields();

        // 3. Define mock Meteorological profiles
        std::vector<double> lat(n_cols, 40.0);
        std::vector<double> lon(n_cols, -105.0);
        std::vector<double> temperature(n_cols * n_levels, 280.0);
        std::vector<double> airden(n_cols * n_levels, 1.2);
        std::vector<double> pedge(n_cols * (n_levels + 1), 101300.0);
        std::vector<double> bxheight(n_cols * n_levels, 100.0);

        // Populate pedagogical profiles
        for (int i = 0; i < n_levels; ++i) {
            temperature[i] = 280.0 - 0.5 * i;
            airden[i] = 1.2 * std::exp(-i / 10.0);
            bxheight[i] = 1000.0; // 1km thick layers
            pedge[i] = 101300.0 * std::exp(-i / 10.0);
        }
        pedge[n_levels] = 101300.0 * std::exp(-n_levels / 10.0);

        // Bind standard meteorological fields to StateManager
        state->bind_met_field_2d("LAT", lat.data());
        state->bind_met_field_2d("LON", lon.data());
        state->bind_met_field_3d("T", temperature.data());
        state->bind_met_field_3d("AIRDEN", airden.data());
        state->bind_met_field_3d("PEDGE", pedge.data());
        state->bind_met_field_3d("BXHEIGHT", bxheight.data());

        // 4. Resolve the TUV-x configuration file path explicitly using configured header
        std::string config_path =
            std::string(catchem::test::SOURCE_DIR) + "/src/external/musica/configs/tuvx/from_host/config.json";

        std::cout << "DEBUG: Using TUV-x config path: " << config_path << std::endl;
        assert(file_exists(config_path) &&
               "Error: TUV-x configuration file could not be located relative to test runner!");

        // Write a temp main config for propagation test
        std::string temp_main_config = "test_main_config.yml";
        std::ofstream main_conf_writer(temp_main_config);
        main_conf_writer << "process:\n";
        main_conf_writer << "  photolysis:\n";
        main_conf_writer << "    active: true\n";
        main_conf_writer << "    config_file: \"" << config_path << "\"\n";
        main_conf_writer.close();

        // Propagate config file path
        state->config_file_path = temp_main_config;
        if (state->config_mgr) {
            state->config_mgr->load_from_file(temp_main_config);
        }

        // 5. Create and initialize the photolysis process
        auto process = catchem::ProcessRegistry::get_instance().create("photolysis");
        process->init(state);
        core->add_process(process);

        // 6. Run a single timestep (forcing photolysis calculations)
        core->run_timestep(3600.0);

        std::cout << "SUCCESS: Executed run_timestep without error." << std::endl;

        // 7. Verify the output diagnostic photolysis rates are registered and non-zero
        auto diag_mgr = core->get_diagnostic_manager();
        assert(diag_mgr->has_field("photolysis_rate_jfoo"));
        std::cout << "SUCCESS: jfoo photolysis reaction rate diagnostic registered dynamically." << std::endl;

        double* jfoo_rates = static_cast<double*>(diag_mgr->get_host_pointer("photolysis_rate_jfoo"));
        assert(jfoo_rates != nullptr);

        double total_rate = 0.0;
        for (int i = 0; i < n_levels; ++i) {
            double rate = jfoo_rates[i];
            std::cout << "DEBUG: Level " << i << " photolysis rate J(jfoo): " << rate << " s-1" << std::endl;
            assert(std::isfinite(rate) && "Photolysis rate must be finite!");
            assert(rate >= 0.0 && "Photolysis rate must be non-negative!");
            total_rate += rate;
        }

        assert(total_rate > 0.0 && "Photolysis rates should be non-zero for positive solar zenith angles!");
        std::cout << "SUCCESS: Calculated non-zero finite J-rates for all levels." << std::endl;

        // Clean up temp file
        std::remove(temp_main_config.c_str());

        std::cout << "\n==========================================" << std::endl;
        std::cout << "TEST PASSED SUCCESSFULLY!" << std::endl;
        std::cout << "==========================================\n" << std::endl;
    }
    Kokkos::finalize();
    return 0;
}
