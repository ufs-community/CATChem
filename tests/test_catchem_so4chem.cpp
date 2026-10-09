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
#include <yaml-cpp/yaml.h>

extern "C" {
void catchem_register_so4chem_cpp();
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "RUNNING TEST: SO4chem Process Unit Test" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        catchem_register_so4chem_cpp();
        assert(catchem::ProcessRegistry::get_instance().has_process("so4chem"));

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

        std::vector<double> lat(n_cols, 40.0);
        std::vector<double> lon(n_cols, -100.0);
        std::vector<double> hflux(n_cols, 10.0);
        std::vector<double> pblh(n_cols, 1000.0);
        std::vector<double> ustar(n_cols, 0.3);
        std::vector<double> u10m(n_cols, 3.0), v10m(n_cols, 1.0), lwi(n_cols, 1.0), z0(n_cols, 0.01);
        std::vector<double> temperature(n_cols * n_levels, 280.0);
        std::vector<double> airden(n_cols * n_levels, 1.2);
        std::vector<double> pmid(n_cols * n_levels, 100000.0);
        std::vector<double> pedge(n_cols * (n_levels + 1), 101300.0);
        std::vector<double> bxheight(n_cols * n_levels, 100.0);
        std::vector<double> cldf(n_cols * n_levels, 0.2);
        std::vector<double> z(n_cols * (n_levels + 1), 0.0);
        std::vector<double> mairden(n_cols * n_levels, 1.2);
        std::vector<double> delp(n_cols * n_levels, 1000.0);
        std::vector<double> chem_conc(n_cols * n_levels * n_species, 1.0e-8);

        state->bind_met_field_2d("LAT", lat.data());
        state->bind_met_field_2d("LON", lon.data());
        state->bind_met_field_2d("HFLUX", hflux.data());
        state->bind_met_field_2d("PBLH", pblh.data());
        state->bind_met_field_2d("USTAR", ustar.data());
        state->bind_met_field_2d("U10M", u10m.data());
        state->bind_met_field_2d("V10M", v10m.data());
        state->bind_met_field_2d("LWI", lwi.data());
        // NUOPC supplies Z0 but not Z0H. SO4chem must derive Z0H = 0.1 * Z0.
        state->bind_met_field_2d("Z0", z0.data());
        state->bind_met_field_3d("T", temperature.data());
        state->bind_met_field_3d("AIRDEN", airden.data());
        state->bind_met_field_3d("AIRDEN_DRY", airden.data());
        state->bind_met_field_3d("PMID", pmid.data());
        state->bind_met_field_3d("PEDGE", pedge.data());
        state->bind_met_field_3d("Z", z.data());
        state->bind_met_field_3d("BXHEIGHT", bxheight.data());
        state->bind_met_field_3d("CLDF", cldf.data());
        state->bind_met_field_3d("MAIRDEN", mairden.data());
        state->bind_met_field_3d("DELP", delp.data());
        state->bind_unified_chemistry(chem_conc.data());

        auto so4chem = catchem::ProcessRegistry::get_instance().create("so4chem");
        assert(so4chem != nullptr);
        so4chem->init(state);
        so4chem->run(state);
        state->sync_to_host();

        // The shared config selects processes.so4chem.diag_species = [so2, so4],
        // so each gets its own Production_rate_<sp> field.  Regression: the
        // bridge used to copy the FIRST diagnostic slot into every field, so
        // Production_rate_so4 duplicated Production_rate_so2.  Each field must
        // carry its OWN species' production rate, so the two buffers must not
        // be identical, and the SO4 field must be non-zero (the scheme fills
        // the SO4 slot from SO2 oxidation).
        {
            const auto manager = core->get_diagnostic_manager();
            assert(manager->has_field("Production_rate_so2"));
            assert(manager->has_field("Production_rate_so4"));
            const double* a = (const double*)manager->get_host_pointer("Production_rate_so2");
            const double* b = (const double*)manager->get_host_pointer("Production_rate_so4");
            assert(a != nullptr && b != nullptr);
            bool differ = false;
            bool so4_nonzero = false;
            for (int i = 0; i < n_cols * n_levels; ++i) {
                if (a[i] != b[i])
                    differ = true;
                if (b[i] != 0.0)
                    so4_nonzero = true;
            }
            assert(differ && "Production_rate fields must not duplicate the first slot");
            assert(so4_nonzero && "Production_rate_so4 must carry the SO4 slot, not zeros");
            std::cout << "  PASS per-slot production-rate fields differ" << std::endl;
        }

        std::cout << "SUCCESS: SO4chem process executed successfully." << std::endl;
    }

    // --- US2 do_drydep routing gate (spec US2.3) -----------------------------
    // Differential runs over an identical fixture: do_drydep=false keeps the
    // HGHTE sentinel that suppresses dry deposition inside SulfateChemDriver
    // (routing stays with the drydep process), while do_drydep=true publishes
    // the physical boundary height so the driver deposits SO2/SO4/MSA itself,
    // GOCART-faithful.  Non-sulfur species must be bit-identical between the
    // two runs, proving the flag only moves sulfate routing.
    {
        int n_cols = 4;
        int n_levels = 5;
        int n_species = 22;

        auto run_case = [&](bool do_drydep, std::vector<double>& conc_out) {
            auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
            auto state = core->get_state_manager();
            auto runtime_config = std::make_shared<catchem::ConfigManager>();
            runtime_config->load_from_file("CATChem_new_config.yml");
            {
                YAML::Node settings;
                settings["gocart"]["update_so2"] = true;
                settings["gocart"]["do_drydep"] = do_drydep;
                auto& proc = runtime_config->data.processes["so4chem"];
                proc.activate = true;
                proc.set_settings_node(settings);
            }
            state->attach_config_manager(runtime_config);

            std::string species_path = "CATChem_species.yml";
            for (const std::string& candidate :
                 {species_path, "tests/" + species_path, "../tests/" + species_path, "../../tests/" + species_path}) {
                std::ifstream f(candidate);
                if (f.good()) {
                    species_path = candidate;
                    break;
                }
            }
            state->load_species_config(species_path);

            std::vector<double> lat(n_cols, 40.0);
            std::vector<double> lon(n_cols, -100.0);
            std::vector<double> hflux(n_cols, 10.0);
            std::vector<double> pblh(n_cols, 1000.0);
            std::vector<double> ustar(n_cols, 0.3);
            std::vector<double> u10m(n_cols, 3.0), v10m(n_cols, 1.0), lwi(n_cols, 1.0), z0(n_cols, 0.01);
            std::vector<double> temperature(n_cols * n_levels, 280.0);
            std::vector<double> airden(n_cols * n_levels, 1.2);
            std::vector<double> pmid(n_cols * n_levels, 100000.0);
            std::vector<double> pedge(n_cols * (n_levels + 1), 101300.0);
            std::vector<double> bxheight(n_cols * n_levels, 100.0);
            std::vector<double> cldf(n_cols * n_levels, 0.2);
            // Physically ascending interface heights: DryDeposition derives the
            // surface-layer dz from HGHTE, and an all-zero column would make
            // the deposition frequency identically zero in both runs, hiding
            // the do_drydep difference.  The sentinel trick (do_drydep=false)
            // then inflates only the bottom dz, suppressing sulfate drydep.
            std::vector<double> z(n_cols * (n_levels + 1));
            for (int lev = 0; lev <= n_levels; ++lev)
                for (int c = 0; c < n_cols; ++c)
                    z[static_cast<size_t>(c) + static_cast<size_t>(lev) * n_cols] = 100.0 * lev;
            std::vector<double> mairden(n_cols * n_levels, 1.2);
            std::vector<double> delp(n_cols * n_levels, 1000.0);
            std::vector<double> chem_conc(n_cols * n_levels * n_species, 1.0e-8);

            state->bind_met_field_2d("LAT", lat.data());
            state->bind_met_field_2d("LON", lon.data());
            state->bind_met_field_2d("HFLUX", hflux.data());
            state->bind_met_field_2d("PBLH", pblh.data());
            state->bind_met_field_2d("USTAR", ustar.data());
            state->bind_met_field_2d("U10M", u10m.data());
            state->bind_met_field_2d("V10M", v10m.data());
            state->bind_met_field_2d("LWI", lwi.data());
            state->bind_met_field_2d("Z0", z0.data());
            state->bind_met_field_3d("T", temperature.data());
            state->bind_met_field_3d("AIRDEN", airden.data());
            state->bind_met_field_3d("AIRDEN_DRY", airden.data());
            state->bind_met_field_3d("PMID", pmid.data());
            state->bind_met_field_3d("PEDGE", pedge.data());
            state->bind_met_field_3d("Z", z.data());
            state->bind_met_field_3d("BXHEIGHT", bxheight.data());
            state->bind_met_field_3d("CLDF", cldf.data());
            state->bind_met_field_3d("MAIRDEN", mairden.data());
            state->bind_met_field_3d("DELP", delp.data());
            state->bind_unified_chemistry(chem_conc.data());

            auto so4chem = catchem::ProcessRegistry::get_instance().create("so4chem");
            so4chem->init(state);
            so4chem->run(state);
            state->sync_to_host();
            conc_out = chem_conc;
        };

        std::vector<double> off_conc, on_conc;
        run_case(false, off_conc);
        run_case(true, on_conc);

        // Default catalog slots: so2=0, h2o2=1, oh=2, no3=3, so4=4, msa=7.
        // Chemistry layout is species-outermost column-major: element
        // (col, level, species) = species*n_cols*n_levels + level*n_cols + col,
        // so each species occupies a contiguous block of n_cols*n_levels.
        const int slab = n_cols * n_levels;
        auto species_total = [&](const std::vector<double>& c, int species) {
            double total = 0.0;
            for (int i = 0; i < slab; ++i)
                total += c[static_cast<size_t>(species) * slab + static_cast<size_t>(i)];
            return total;
        };
        for (int species : {0, 4, 7}) {
            const bool deposited = species_total(on_conc, species) < species_total(off_conc, species);
            std::cout << (deposited ? "  PASS: " : "  FAIL: ") << "do_drydep=true deposits species slot " << species
                      << std::endl;
            assert(deposited);
        }
        // OH (2) and NO3 (3) are the oxidants the scheme echoes unchanged every
        // step.  H2O2 (1) is deliberately written back each step — it is
        // depleted by aqueous SO2 oxidation inside the driver regardless of
        // do_drydep — so it is excluded from the bit-identical set.
        for (int species : {2, 3}) {
            // Compare the full species block bit-exactly, not just the sum.
            bool slot_exact = true;
            for (int i = 0; i < slab; ++i)
                slot_exact = slot_exact && (on_conc[static_cast<size_t>(species) * slab + i] ==
                                            off_conc[static_cast<size_t>(species) * slab + i]);
            std::cout << (slot_exact ? "  PASS: " : "  FAIL: ") << "do_drydep leaves non-sulfur slot " << species
                      << " bit-identical" << std::endl;
            assert(slot_exact);
        }
    }
    Kokkos::finalize();
    return 0;
}
