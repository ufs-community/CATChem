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
void catchem_register_carbchem_cpp();
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "RUNNING TEST: CarbChem Process Unit Test" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        catchem_register_carbchem_cpp();
        assert(catchem::ProcessRegistry::get_instance().has_process("carbchem"));

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
        std::vector<double> airden_dry(n_cols * n_levels, 1.2);
        std::vector<double> delp(n_cols * n_levels, 1000.0);
        std::vector<double> pmid(n_cols * n_levels, 90000.0);
        std::vector<double> chem_conc(n_cols * n_levels * n_species, 1.0e-8);

        state->bind_met_field_3d("T", temperature.data());
        state->bind_met_field_3d("AIRDEN_DRY", airden_dry.data());
        state->bind_met_field_3d("DELP", delp.data());
        state->bind_met_field_3d("PMID", pmid.data());
        state->bind_unified_chemistry(chem_conc.data());

        auto carbchem = catchem::ProcessRegistry::get_instance().create("carbchem");
        assert(carbchem != nullptr);
        carbchem->init(state);

        // The default diagnostic set is the explicit carbon species (oc1, oc2,
        // bc1, bc2) resolved against the mechanism, NOT the whole catalog.
        // Packed fields therefore carry a species dimension of n_diag, not
        // n_species.  The Default mechanism carries all four carbon species.
        {
            const auto manager = core->get_diagnostic_manager();
            assert(manager->has_field("carbchem_prod_mass"));
            assert(manager->get_field("carbchem_prod_mass")->dimensions ==
                   std::vector<int>({n_cols, n_levels, 4}));
            assert(manager->get_field("carbchem_loss_flux")->dimensions == std::vector<int>({n_cols, 4}));
            assert(manager->get_field("carbchem_phobic_mass")->dimensions ==
                   std::vector<int>({n_cols, n_levels, 4}));
            assert(manager->get_field("carbchem_phobic_flux")->dimensions == std::vector<int>({n_cols, 4}));
            std::cout << "  PASS default carbon set: fields register with species dim = 4" << std::endl;
        }

        // A diag_species subset must shrink the species dimension to 1 and map
        // to the GLOBAL 1-based species index (the space the scheme matches on).
        // register_field throws on different-dims re-registration, so this runs
        // on a fresh Core.  init() reads only config + mechanism.
        {
            auto core2 = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
            auto state2 = core2->get_state_manager();
            auto cfg2 = std::make_shared<catchem::ConfigManager>();
            cfg2->load_from_file("CATChem_new_config.yml");
            cfg2->data.processes["carbchem"].diagnostics = true;
            cfg2->data.processes["carbchem"].diag_species = {"bc2"};
            state2->attach_config_manager(cfg2);
            state2->load_species_config(species_path);
            auto carbchem2 = catchem::ProcessRegistry::get_instance().create("carbchem");
            carbchem2->init(state2);
            const auto mgr2 = core2->get_diagnostic_manager();
            assert(mgr2->get_field("carbchem_prod_mass")->dimensions == std::vector<int>({n_cols, n_levels, 1}));
            assert(mgr2->get_field("carbchem_loss_flux")->dimensions == std::vector<int>({n_cols, 1}));
            std::cout << "  PASS diag_species subset: fields register with species dim = 1" << std::endl;
        }

        carbchem->run(state);
        state->sync_to_host();

        std::cout << "SUCCESS: CarbChem process executed successfully." << std::endl;
    }
    Kokkos::finalize();
    return 0;
}
