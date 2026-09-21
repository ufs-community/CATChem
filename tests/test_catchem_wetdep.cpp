#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_kokkos_compat.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_state_manager.hpp"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

extern "C" {
void catchem_register_wetdep_cpp();
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "RUNNING TEST: WetDep Process Unit Test" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        catchem_register_wetdep_cpp();
        assert(catchem::ProcessRegistry::get_instance().has_process("wetdep"));

        int n_cols = 4;
        int n_levels = 5;
        int n_species = 22;

        auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
        auto state = core->get_state_manager();
        auto runtime_config = std::make_shared<catchem::ConfigManager>();
        runtime_config->load_from_file("CATChem_new_config.yml");
        state->attach_config_manager(runtime_config);

        std::string species_path = "CATChem_species.yml";
        std::vector<std::string> candidates = {species_path, "tests/" + species_path, "../tests/" + species_path,
                                               "../../tests/" + species_path};
        for (const auto& candidate : candidates) {
            std::ifstream f(candidate);
            if (f.good()) {
                species_path = candidate;
                break;
            }
        }
        state->load_species_config(species_path);

        std::vector<double> temperature(n_cols * n_levels, 280.0);
        std::vector<double> airden(n_cols * n_levels, 1.2);
        std::vector<double> pedge(n_cols * (n_levels + 1), 101300.0);
        std::vector<double> bxheight(n_cols * n_levels, 100.0);
        std::vector<double> cldf(n_cols * n_levels, 0.3);
        std::vector<double> pfilsan(n_cols * n_levels, 1.0e-5);
        std::vector<double> pfllsan(n_cols * n_levels, 1.0e-5);
        std::vector<double> reevapls(n_cols * n_levels, 0.0);
        std::vector<double> delp(n_cols * n_levels, 1000.0);
        std::vector<double> chem_conc(n_cols * n_levels * n_species, 1.0e-8);

        state->bind_met_field_3d("T", temperature.data());
        state->bind_met_field_3d("AIRDEN", airden.data());
        state->bind_met_field_3d("AIRDEN_DRY", airden.data());
        state->bind_met_field_3d("PEDGE", pedge.data());
        state->bind_met_field_3d("BXHEIGHT", bxheight.data());
        state->bind_met_field_3d("CLDF", cldf.data());
        state->bind_met_field_3d("PFILSAN", pfilsan.data());
        state->bind_met_field_3d("PFLLSAN", pfllsan.data());
        state->bind_met_field_3d("REEVAPLS", reevapls.data());
        state->bind_met_field_3d("DELP", delp.data());
        state->bind_unified_chemistry(chem_conc.data());

        auto wetdep = catchem::ProcessRegistry::get_instance().create("wetdep");
        assert(wetdep != nullptr);
        wetdep->init(state);
        wetdep->run(state);
        state->sync_to_host();

        std::cout << "SUCCESS: WetDep process executed successfully." << std::endl;
    }
    Kokkos::finalize();
    return 0;
}
