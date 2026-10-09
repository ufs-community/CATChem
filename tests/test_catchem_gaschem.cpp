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
void catchem_register_gaschem_cpp();
}

bool file_exists(const std::string& name) {
    return std::filesystem::exists(name);
}

// Runs the full coupled photolysis+gaschem scenario once for a given gaschem
// mechanism configuration (a v0 CAMP directory, a v1 single file, or a v1
// directory relying on config.yaml auto-detection).
static void run_coupled_case(const std::string& label, const std::string& config_key, const std::string& config_value) {
    std::cout << "\n------------------------------------------" << std::endl;
    std::cout << "CASE: gaschem mechanism via " << config_key << " (" << label << ")" << std::endl;
    std::cout << "------------------------------------------\n" << std::endl;

    // 2. Set up core and states
    int n_cols = 1;
    int n_levels = 3;
    int n_species = 4; // matches tests/fixtures/mechanisms/chapman.yml

    auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
    auto state = core->get_state_manager();

    // Time setup (Noon, summer)
    state->clock().year = 2026;
    state->clock().month = 7;
    state->clock().day = 13;
    state->clock().hour = 12;
    state->clock().minute = 0;
    state->clock().second = 0;
    state->clock().calculate_derived_fields();

    // Meteorological arrays
    std::vector<double> lat(n_cols, 40.0);
    std::vector<double> lon(n_cols, -105.0);
    std::vector<double> temperature(n_cols * n_levels, 280.0);
    std::vector<double> airden(n_cols * n_levels, 1.2);
    std::vector<double> pedge(n_cols * (n_levels + 1), 101300.0);
    std::vector<double> bxheight(n_cols * n_levels, 100.0);

    for (int i = 0; i < n_levels; ++i) {
        temperature[i] = 280.0 - 0.5 * i;
        airden[i] = 1.2 * std::exp(-i / 10.0);
        bxheight[i] = 1000.0;
        pedge[i] = 101300.0 * std::exp(-i / 10.0);
    }
    pedge[n_levels] = 101300.0 * std::exp(-n_levels / 10.0);

    state->bind_met_field_2d("LAT", lat.data());
    state->bind_met_field_2d("LON", lon.data());
    state->bind_met_field_3d("T", temperature.data());
    state->bind_met_field_3d("AIRDEN", airden.data());
    state->bind_met_field_3d("AIRDEN_DRY", airden.data());
    state->bind_met_field_3d("PEDGE", pedge.data());
    state->bind_met_field_3d("PMID", pedge.data()); // PMID maps to PMID in tests
    state->bind_met_field_3d("BXHEIGHT", bxheight.data());

    // Load a mechanism whose species mirror chapman v0 (so every CATChem
    // species maps into the MICM variable map) and that carries the
    // photolysis.ozone role the photolysis contract requires.
    std::string species_config = std::string(catchem::test::TEST_DIR) + "/fixtures/mechanisms/chapman.yml";
    assert(file_exists(species_config) && "ERROR: Could not find the chapman mechanism fixture!");
    state->load_species_config(species_config);
    std::vector<double> conc_data(n_cols * n_levels * n_species, 1.0); // 1.0 ppmv initially
    state->bind_unified_chemistry(conc_data.data());

    // 3. Resolve configs explicitly using configured header
    std::string photolysis_config =
        std::string(catchem::test::SOURCE_DIR) + "/src/external/musica/configs/tuvx/from_host/config.json";
    assert(file_exists(photolysis_config) && "ERROR: Could not find photolysis config explicitly!");

    // The gaschem mechanism under test is supplied by the caller: either a v0
    // CAMP directory (config_dir) or a v1 single-file mechanism (config_file).
    std::string temp_main_coupled_config = "test_main_coupled_config.yml";
    std::ofstream main_conf_writer(temp_main_coupled_config);
    main_conf_writer << "process:\n";
    main_conf_writer << "  photolysis:\n";
    main_conf_writer << "    active: true\n";
    main_conf_writer << "    config_file: \"" << photolysis_config << "\"\n";
    main_conf_writer << "  gaschem:\n";
    main_conf_writer << "    active: true\n";
    main_conf_writer << "    " << config_key << ": \"" << config_value << "\"\n";
    main_conf_writer.close();

    state->set_configuration_path(temp_main_coupled_config);
    if (state->config_manager()) {
        state->config_manager()->load_from_file(temp_main_coupled_config);
    }

    // 4. Create and add processes to core pipeline
    auto photolysis_proc = catchem::ProcessRegistry::get_instance().create("photolysis");
    photolysis_proc->init(state);
    core->add_process(photolysis_proc);

    auto gaschem_proc = catchem::ProcessRegistry::get_instance().create("gaschem");
    gaschem_proc->init(state);
    core->add_process(gaschem_proc);

    // 5. Execute timestep
    core->run_timestep(3600.0);
    std::cout << "SUCCESS: Timestep run completed without error." << std::endl;

    // 6. Verify photolysis rate diagnostics calculation
    auto diag_mgr = core->get_diagnostic_manager();
    assert(diag_mgr->has_field("photolysis_rate_jfoo") && "TUV-x jfoo diagnostic must exist!");

    double* jfoo_rates = static_cast<double*>(diag_mgr->get_host_pointer("photolysis_rate_jfoo"));
    assert(jfoo_rates != nullptr);
    assert(jfoo_rates[0] > 0.0 && "Photolysis rates must be non-zero at noon!");

    std::cout << "SUCCESS: Dynamic J-rates calculated and coupled successfully!" << std::endl;

    std::remove(temp_main_coupled_config.c_str());
    std::cout << "CASE PASSED: " << label << std::endl;
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "RUNNING INTEGRATION TEST: Photolysis and GasChem Coupling" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        // 1. Dynamic registrations
        catchem_register_photolysis_cpp();
        catchem_register_gaschem_cpp();
        assert(catchem::ProcessRegistry::get_instance().has_process("photolysis"));
        assert(catchem::ProcessRegistry::get_instance().has_process("gaschem"));

        // Case A: legacy v0 CAMP mechanism directory (multi-file layout).
        std::string v0_dir = std::string(catchem::test::SOURCE_DIR) + "/src/external/musica/configs/v0/chapman/";
        assert((file_exists(v0_dir + "config.json") || file_exists(v0_dir + "config.yaml")) &&
               "ERROR: Could not find gaschem v0 config dir explicitly!");
        run_coupled_case("v0 CAMP directory", "config_dir", v0_dir);

        // Case B: MUSICA v1 single-file mechanism passed via the new config_file key.
        std::string v1_file =
            std::string(catchem::test::SOURCE_DIR) + "/src/external/musica/configs/v1/chapman/config.yaml";
        assert(file_exists(v1_file) && "ERROR: Could not find gaschem v1 config file explicitly!");
        run_coupled_case("v1 single file (config_file)", "config_file", v1_file);

        // Case C: v1 mechanism directory exercising config.yaml auto-detection.
        std::string v1_dir = std::string(catchem::test::SOURCE_DIR) + "/src/external/musica/configs/v1/chapman/";
        assert(file_exists(v1_dir + "config.yaml") && "ERROR: Could not find gaschem v1 config dir explicitly!");
        run_coupled_case("v1 directory (config.yaml auto-detect)", "config_dir", v1_dir);

        std::cout << "\n==========================================" << std::endl;
        std::cout << "COUPLED INTEGRATION TESTS PASSED (v0 dir, v1 file, v1 auto-detect)!" << std::endl;
        std::cout << "==========================================\n" << std::endl;
    }
    Kokkos::finalize();
    return 0;
}
