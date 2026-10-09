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
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
void catchem_register_seasalt_cpp();
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "RUNNING TEST: SeaSalt Process Unit Test" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        catchem_register_seasalt_cpp();
        assert(catchem::ProcessRegistry::get_instance().has_process("seasalt"));

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

        std::vector<double> u10m(n_cols, 12.0); // 12 m/s wind to excite seasalt emissions
        std::vector<double> v10m(n_cols, 0.0);
        std::vector<double> sst(n_cols, 290.0);
        std::vector<double> frocean(n_cols, 1.0);
        std::vector<double> frseaice(n_cols, 0.0);
        std::vector<double> lat(n_cols, 10.0);
        std::vector<double> lon(n_cols, 100.0);
        std::vector<double> ustar(n_cols, 0.5);
        std::vector<double> airden(n_cols * n_levels, 1.2);
        std::vector<double> bxheight(n_cols * n_levels, 100.0);
        std::vector<double> delp(n_cols * n_levels, 1000.0);
        std::vector<double> chem_conc(n_cols * n_levels * n_species, 0.0);

        state->bind_met_field_2d("U10M", u10m.data());
        state->bind_met_field_2d("V10M", v10m.data());
        state->bind_met_field_2d("SST", sst.data());
        state->bind_met_field_2d("FROCEAN", frocean.data());
        state->bind_met_field_2d("FRSEAICE", frseaice.data());
        state->bind_met_field_2d("LAT", lat.data());
        state->bind_met_field_2d("LON", lon.data());
        state->bind_met_field_2d("USTAR", ustar.data());
        state->bind_met_field_3d("AIRDEN", airden.data());
        state->bind_met_field_3d("BXHEIGHT", bxheight.data());
        state->bind_met_field_3d("DELP", delp.data());
        state->bind_unified_chemistry(chem_conc.data());

        auto seasalt = catchem::ProcessRegistry::get_instance().create("seasalt");
        assert(seasalt != nullptr);
        seasalt->init(state);

        // Per-bin emissions must register as compact [ncols, n_seasalt] fields
        // (one mass + one number array), replacing the former per-species 2D
        // fields, so the NUOPC driver writes a single 3D (nx, ny, nbin)
        // variable.  Totals stay [ncols, 1].
        {
            const auto manager = core->get_diagnostic_manager();
            int n_seasalt = 0;
            for (const auto& meta : state->chemistry().species_list)
                if (meta.is_seasalt)
                    ++n_seasalt;
            assert(n_seasalt > 0);
            assert(manager->has_field("seasalt_mass_emission_bins"));
            assert(manager->get_field("seasalt_mass_emission_bins")->dimensions ==
                   std::vector<int>({n_cols, n_seasalt}));
            assert(manager->has_field("seasalt_number_emission_bins"));
            assert(manager->get_field("seasalt_number_emission_bins")->dimensions ==
                   std::vector<int>({n_cols, n_seasalt}));
            assert(manager->get_field("seasalt_mass_emission_total")->dimensions == std::vector<int>({n_cols, 1}));
        }

        // A diag_species subset must shrink the per-bin fields to [ncols, 1]
        // and map the id to the LOCAL bin position (SEAS3 -> 3), matching the
        // space the schemes search.  register_field throws on different-dims
        // re-registration, so this runs on a fresh Core.  init() reads only
        // the config and the chemistry mechanism, so no met fields are needed.
        {
            auto core2 = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
            auto state2 = core2->get_state_manager();
            auto cfg2 = std::make_shared<catchem::ConfigManager>();
            cfg2->load_from_file("CATChem_new_config.yml");
            cfg2->data.processes["seasalt"].diagnostics = true;
            cfg2->data.processes["seasalt"].diag_species = {"seas3"};
            state2->attach_config_manager(cfg2);
            state2->load_species_config(species_path);
            auto seasalt2 = catchem::ProcessRegistry::get_instance().create("seasalt");
            seasalt2->init(state2);
            const auto mgr2 = core2->get_diagnostic_manager();
            assert(mgr2->get_field("seasalt_mass_emission_bins")->dimensions == std::vector<int>({n_cols, 1}));
            assert(mgr2->get_field("seasalt_number_emission_bins")->dimensions == std::vector<int>({n_cols, 1}));
            assert(mgr2->get_field("seasalt_mass_emission_total")->dimensions == std::vector<int>({n_cols, 1}));
            std::cout << "  PASS diag_species subset: per-bin fields register as [ncols, 1]" << std::endl;
        }

        seasalt->run(state);
        state->sync_to_host();

        std::cout << "SUCCESS: SeaSalt process executed successfully." << std::endl;

        auto missing_field_core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
        auto missing_field_state = missing_field_core->get_state_manager();
        auto missing_runtime_config = std::make_shared<catchem::ConfigManager>();
        missing_runtime_config->load_from_file("CATChem_new_config.yml");
        missing_field_state->attach_config_manager(missing_runtime_config);
        missing_field_state->load_species_config(species_path);

        std::vector<double> missing_sst(n_cols, 290.0);
        std::vector<double> missing_frocean(n_cols, 1.0);
        std::vector<double> missing_frseaice(n_cols, 0.0);
        std::vector<double> missing_lat(n_cols, 10.0);
        std::vector<double> missing_lon(n_cols, 100.0);
        std::vector<double> missing_delp(n_cols * n_levels, 1000.0);
        std::vector<double> missing_chem_conc(n_cols * n_levels * n_species, 0.0);

        missing_field_state->bind_met_field_2d("SST", missing_sst.data());
        missing_field_state->bind_met_field_2d("FROCEAN", missing_frocean.data());
        missing_field_state->bind_met_field_2d("FRSEAICE", missing_frseaice.data());
        missing_field_state->bind_met_field_2d("LAT", missing_lat.data());
        missing_field_state->bind_met_field_2d("LON", missing_lon.data());
        missing_field_state->bind_met_field_3d("DELP", missing_delp.data());
        missing_field_state->bind_unified_chemistry(missing_chem_conc.data());

        auto missing_field_seasalt = catchem::ProcessRegistry::get_instance().create("seasalt");
        assert(missing_field_seasalt != nullptr);
        missing_field_seasalt->init(missing_field_state);

        bool saw_missing_ustar = false;
        try {
            missing_field_seasalt->run(missing_field_state);
        } catch (const std::runtime_error& error) {
            saw_missing_ustar = std::string(error.what()).find("USTAR") != std::string::npos;
        }
        assert(saw_missing_ustar);

        std::cout << "SUCCESS: SeaSalt process rejected missing USTAR." << std::endl;
    }
    Kokkos::finalize();
    return 0;
}
