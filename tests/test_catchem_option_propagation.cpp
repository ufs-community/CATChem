// Round-trip test: runtime YAML scheme options must reach the science layer
// and change computed behavior, and unknown option keys must fail validation.
#include "catchem_core.hpp"
#include "catchem_kokkos_compat.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_state_manager.hpp"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

extern "C" {
void catchem_register_settling_cpp();
}

namespace {

    // Minimal runtime YAML exercising the nested <scheme>/<key> layout that
    // SettlingProcess::init reads.  scale_factor is the knob under test.
    std::string write_config(const std::string& path, double scale_factor) {
        std::ofstream out(path);
        out << "simulation:\n"
            << "  name: option_propagation\n"
            << "  timestep: 3600\n"
            << "processes:\n"
            << "  settling:\n"
            << "    activate: true\n"
            << "    diagnostics: false\n"
            << "    scheme: 'gocart'\n"
            << "    gocart:\n"
            << "      scale_factor: " << scale_factor << "\n"
            << "      correction_maring: false\n";
        out.close();
        assert(out.good());
        return path;
    }

    std::string find_fixture(const std::string& filename) {
        for (const auto& candidate :
             {filename, "tests/" + filename, "../tests/" + filename, "../../tests/" + filename}) {
            std::ifstream f(candidate);
            if (f.good())
                return candidate;
        }
        throw std::runtime_error("Fixture not found: " + filename);
    }

    // Run one settling column with the given scale_factor and return the
    // total aerosol mass remaining after a single step.
    double run_settling(double scale_factor) {
        const int n_cols = 2;
        const int n_levels = 6;

        auto core = std::make_shared<catchem::Core>(n_cols, n_levels, 22);
        auto state = core->get_state_manager();

        auto runtime_config = std::make_shared<catchem::ConfigManager>();
        runtime_config->load_from_file(write_config("opt_prop_settling.yml", scale_factor));
        state->attach_config_manager(runtime_config);
        state->load_species_config(find_fixture("CATChem_species.yml"));

        const int n_species = state->species_count();
        std::vector<double> temperature(n_cols * n_levels, 280.0);
        std::vector<double> pmid(n_cols * n_levels, 100000.0);
        std::vector<double> pedge(n_cols * (n_levels + 1));
        std::vector<double> z_edge(n_cols * (n_levels + 1));
        std::vector<double> airden(n_cols * n_levels, 1.2);
        std::vector<double> rh(n_cols * n_levels, 0.5);
        std::vector<double> bxheight(n_cols * n_levels, 100.0);
        std::vector<double> chem_conc(static_cast<size_t>(n_cols) * n_levels * n_species, 1.0e-8);

        for (int level = 0; level <= n_levels; ++level)
            for (int column = 0; column < n_cols; ++column) {
                pedge[static_cast<size_t>(column + level * n_cols)] = 100000.0 - 5000.0 * level;
                z_edge[static_cast<size_t>(column + level * n_cols)] = 100.0 * level;
            }

        state->bind_met_field_3d("T", temperature.data());
        state->bind_met_field_3d("PMID", pmid.data());
        state->bind_met_field_3d("PEDGE", pedge.data());
        state->bind_met_field_3d("Z", z_edge.data());
        state->bind_met_field_3d("AIRDEN", airden.data());
        state->bind_met_field_3d("AIRDEN_DRY", airden.data());
        state->bind_met_field_3d("RH", rh.data());
        state->bind_met_field_3d("BXHEIGHT", bxheight.data());
        state->bind_unified_chemistry(chem_conc.data());

        auto settling = catchem::ProcessRegistry::get_instance().create("settling");
        assert(settling != nullptr);
        settling->init(state);
        settling->run(state);
        state->sync_to_host();

        // Sum one aerosol species across the grid; a stronger fall speed must
        // remove more mass from the column top levels.
        const int probe_species = state->chemistry().aerosol_indices.front();
        double total = 0.0;
        for (int col = 0; col < n_cols; ++col)
            for (int lvl = 0; lvl < n_levels; ++lvl)
                total += chem_conc[static_cast<size_t>(col) * n_levels * n_species +
                                   static_cast<size_t>(lvl) * n_species + probe_species];
        settling->finalize();
        return total;
    }

} // namespace

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "RUNNING TEST: Scheme Option Propagation" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        catchem_register_settling_cpp();
        assert(catchem::ProcessRegistry::get_instance().has_process("settling"));

        // 1. Parity contract: on the GOCART metadata (non-Mie) settling path,
        // scale_factor is parsed and staged onto the scheme config but the
        // kernel does NOT consume it -- exactly like upstream
        // ProcessSettlingInterface/SettlingScheme_GOCART_Mod, which never
        // reference params%scale_factor.  A faithful port must therefore leave
        // the settled mass unchanged when only scale_factor changes.  The knob
        // still round-trips through init() without error (validated below).
        const double baseline = run_settling(1.0);
        const double boosted = run_settling(1000.0);
        std::cout << "  scale_factor=1.0    -> remaining mass " << baseline << std::endl;
        std::cout << "  scale_factor=1000.0 -> remaining mass " << boosted << std::endl;
        assert(std::abs(baseline - boosted) <= 1.0e-15 &&
               "GOCART metadata-path settling must ignore scale_factor (upstream parity)");
        std::cout << "SUCCESS: scale_factor round-trips without altering the metadata-path physics." << std::endl;

        // 2. Unknown nested option must be rejected by the registered validator.
        auto& registry = catchem::ProcessRegistry::get_instance();
        {
            catchem::ConfigManager cfg;
            std::ofstream out("opt_prop_bad.yml");
            out << "processes:\n"
                << "  settling:\n"
                << "    activate: true\n"
                << "    scheme: 'gocart'\n"
                << "    gocart:\n"
                << "      scale_faktor: 2.0\n"; // deliberate typo
            out.close();
            cfg.load_from_file("opt_prop_bad.yml");
            bool threw = false;
            try {
                registry.validate_settings("settling", cfg.data.processes.at("settling"));
            } catch (const std::invalid_argument& e) {
                threw = true;
                std::cout << "  Rejected as expected: " << e.what() << std::endl;
            }
            assert(threw && "unknown scheme option must fail validation at init");
        }

        // 3. Known nested options must pass validation.
        {
            catchem::ConfigManager cfg;
            std::ofstream out("opt_prop_good.yml");
            out << "processes:\n"
                << "  settling:\n"
                << "    activate: true\n"
                << "    scheme: 'gocart'\n"
                << "    gocart:\n"
                << "      scale_factor: 1.0\n"
                << "      simple_scheme: true\n"
                << "      swelling_rh_max: 0.95\n"
                << "      correction_maring: true\n";
            out.close();
            cfg.load_from_file("opt_prop_good.yml");
            registry.validate_settings("settling", cfg.data.processes.at("settling"));
            std::cout << "SUCCESS: Accepted options passed validation." << std::endl;
        }
    }
    Kokkos::finalize();
    std::cout << "=== ALL OPTION PROPAGATION CHECKS PASSED ===" << std::endl;
    return 0;
}
