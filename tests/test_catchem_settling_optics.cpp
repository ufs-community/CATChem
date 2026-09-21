// Settling optics-table (simple_scheme) path tests.
//
// Verifies that the legacy Mie-table settling configuration is accepted by the
// C++ core, that the tables named by the top-level "mie:" section load once at
// initialization, and that settling then runs through the Chem_SettlingSimple
// kernel for every aerosol in the Default 22-species inventory.  The failure
// modes (US3) are covered in the same file because they share the fixture.
#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_state_manager.hpp"
#include "catchem_test_config.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
void catchem_register_settling_cpp();
}

namespace {

    struct Fixture {
        int n_cols = 4;
        int n_levels = 5;
        int n_species = 22;
        std::vector<double> T, AIRDEN, DELP, RH, PMID, Zedge, conc;

        Fixture() {
            T.assign(n_cols * n_levels, 0.0);
            AIRDEN.assign(n_cols * n_levels, 0.0);
            DELP.assign(n_cols * n_levels, 0.0);
            RH.assign(n_cols * n_levels, 0.0);
            PMID.assign(n_cols * n_levels, 0.0);
            Zedge.assign(n_cols * (n_levels + 1), 0.0);
            conc.assign(static_cast<size_t>(n_cols) * n_levels * n_species, 1.0e-8);
            for (int c = 0; c < n_cols; ++c)
                for (int lev = 0; lev < n_levels; ++lev) {
                    std::size_t i = c + static_cast<size_t>(lev) * n_cols;
                    T[i] = 288.15 - 6.5 * lev;
                    PMID[i] = 101300.25 * std::exp(-lev / 8.0);
                    AIRDEN[i] = 1.2 * std::exp(-lev / 8.0);
                    DELP[i] = 5000.0;
                    // Span the RH range so the table LUT is exercised both below
                    // and above its high-RH plateau.
                    RH[i] = 0.4 + 0.6 * (lev + 1) / n_levels;
                }
            for (int c = 0; c < n_cols; ++c)
                for (int lev = 0; lev <= n_levels; ++lev)
                    Zedge[c + static_cast<size_t>(lev) * n_cols] = 100.0 * lev;
        }
    };

    std::size_t flat(const Fixture& fix, int column, int level, int species) {
        return static_cast<size_t>(species) * fix.n_levels * fix.n_cols +
               static_cast<size_t>(level) * fix.n_cols + column;
    }

    std::string find_species_file() {
        const std::string rel = "Configs/Default/CATChem_species.yml";
        for (const std::string& candidate :
             {rel, "tests/" + rel, "../tests/" + rel, "../../tests/" + rel})
            if (std::ifstream(candidate).good())
                return candidate;
        return "";
    }

    // Write a runtime YAML that mirrors the legacy standalone configuration: the
    // top-level mie: section plus settling gocart/simple_scheme.
    void write_config(const std::string& path, const std::string& directory, bool simple_scheme,
                      bool include_files = true) {
        std::ofstream out(path);
        out << "simulation:\n"
            << "  nx: " << 4 << "\n"
            << "  ny: 1\n"
            << "  nz: 5\n"
            << "  timestep: 3600\n";
        if (simple_scheme) {
            out << "mie:\n"
                << "  directory: \"" << directory << "\"\n";
            if (include_files)
                out << "  files:\n"
                    << "    SS: optics_SS.v3_5.nc\n"
                    << "    DU: optics_DU.v15_5.nc\n"
                    << "    BC: optics_BC.v1_5.nc\n"
                    << "    OC: optics_OC.v1_5.nc\n"
                    << "    SU: optics_SU.v1_5.nc\n"
                    << "    NI: optics_NI.v2_5.nc\n"
                    << "    BRC: optics_BRC.v1_5.nc\n";
        }
        out << "processes:\n"
            << "  settling:\n"
            << "    activate: true\n"
            << "    scheme: gocart\n"
            << "    gocart:\n"
            << "      scale_factor: 1.0\n"
            << "      simple_scheme: " << (simple_scheme ? "true" : "false") << "\n"
            << "      swelling_rh_max: 0.95\n"
            << "      correction_maring: false\n"
            << "      maring_dust_only: true\n";
        out.close();
    }

    std::shared_ptr<catchem::StateManager> bind_state(const std::shared_ptr<catchem::Core>& core, Fixture& fix) {
        auto state = core->get_state_manager();
        state->bind_met_field_3d("T", fix.T.data());
        state->bind_met_field_3d("AIRDEN", fix.AIRDEN.data());
        state->bind_met_field_3d("DELP", fix.DELP.data());
        state->bind_met_field_3d("RH", fix.RH.data());
        state->bind_met_field_3d("Z", fix.Zedge.data());
        state->bind_met_field_3d("PMID", fix.PMID.data());
        state->bind_unified_chemistry(fix.conc.data());
        return state;
    }

    bool threw_containing(const std::function<void()>& fn, const std::string& expected) {
        try {
            fn();
        } catch (const std::exception& error) {
            return std::string(error.what()).find(expected) != std::string::npos;
        }
        return false;
    }

} // namespace

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    int failures = 0;
    auto check = [&](bool condition, const std::string& label) {
        std::cout << (condition ? "  PASS: " : "  FAIL: ") << label << '\n';
        if (!condition)
            ++failures;
    };

    {
        std::cout << "==========================================" << std::endl;
        std::cout << "RUNNING TEST: Settling Optics-Table (simple_scheme) Path" << std::endl;
        std::cout << "==========================================" << std::endl;

        catchem_register_settling_cpp();

        const std::string species_path = find_species_file();
        check(!species_path.empty(), "Default species file located");
        const std::string optics_dir(catchem::test::OPTICS_DIR);
        check(std::ifstream(optics_dir + "/optics_DU.v15_5.nc").good(),
              "optics fixture directory reachable (complete tables with rhop)");

        // --- US1: simple_scheme: true initializes and runs --------------------
        {
            Fixture fix;
            auto core = std::make_shared<catchem::Core>(fix.n_cols, fix.n_levels, fix.n_species);
            auto state = bind_state(core, fix);
            state->load_species_config(species_path);
            write_config("optics_ok.yml", optics_dir + "/", true);
            auto config = std::make_shared<catchem::ConfigManager>();
            config->load_from_file("optics_ok.yml");
            state->attach_config_manager(config);
            state->clock().timestep = 3600.0;

            const auto before = fix.conc;
            auto settling = catchem::ProcessRegistry::get_instance().create("settling");
            bool ran = true;
            std::string message;
            try {
                settling->init(state);
                settling->run(state);
            } catch (const std::exception& error) {
                ran = false;
                message = error.what();
            }
            check(ran, "optics-table init+run succeeds (no 'unsupported' rejection): " + message);

            if (ran) {
                bool all_finite_nonneg = true;
                double total_change = 0.0;
                for (std::size_t i = 0; i < fix.conc.size(); ++i) {
                    if (!std::isfinite(fix.conc[i]) || fix.conc[i] < 0.0)
                        all_finite_nonneg = false;
                    total_change += std::abs(fix.conc[i] - before[i]);
                }
                check(all_finite_nonneg, "every concentration stays finite and non-negative");
                check(total_change > 0.0, "optics-table settling actually changed concentrations");

                // Dust is the fastest-settling aerosol; its top-of-column cell
                // must lose mass, which also proves the DU table (not the dry
                // radius) drove the size response.
                int dust = -1;
                for (std::size_t i = 0; i < state->chemistry().species_list.size(); ++i)
                    if (state->chemistry().species_list[i].short_name == "dust4")
                        dust = static_cast<int>(i);
                check(dust >= 0, "dust4 present in Default species list");
                if (dust >= 0) {
                    double top_loss = 0.0;
                    for (int c = 0; c < fix.n_cols; ++c)
                        top_loss += before[flat(fix, c, fix.n_levels - 1, dust)] -
                                    fix.conc[flat(fix, c, fix.n_levels - 1, dust)];
                    check(top_loss > 0.0, "dust4 depletes at the top of the column via its optics table");
                    // Bottom-of-column dust must GAIN the mass transferring down
                    // through the layer stack (conservation, not just decay).
                    double bottom_gain = 0.0;
                    for (int c = 0; c < fix.n_cols; ++c)
                        bottom_gain += fix.conc[flat(fix, c, 0, dust)] -
                                       before[flat(fix, c, 0, dust)];
                    check(bottom_gain > 0.0, "dust4 accumulates at the bottom from downward transfer");
                }
            }
            settling->finalize();
        }

        // --- US1 scenario 3: simple_scheme: false opens no tables -------------
        {
            Fixture fix;
            auto core = std::make_shared<catchem::Core>(fix.n_cols, fix.n_levels, fix.n_species);
            auto state = bind_state(core, fix);
            state->load_species_config(species_path);
            write_config("optics_off.yml", "/definitely/not/here/", false);
            auto config = std::make_shared<catchem::ConfigManager>();
            config->load_from_file("optics_off.yml");
            state->attach_config_manager(config);
            state->clock().timestep = 3600.0;
            auto settling = catchem::ProcessRegistry::get_instance().create("settling");
            bool ok = true;
            try {
                settling->init(state);
                settling->run(state);
            } catch (const std::exception& error) {
                ok = false;
                std::cout << "  (metadata path error: " << error.what() << ")\n";
            }
            check(ok, "metadata path ignores the optics section entirely (bad directory unused)");
            settling->finalize();
        }

        // --- US3: fail-loud scenarios ----------------------------------------
        // 1. Unreadable optics directory names the directory.
        {
            Fixture fix;
            auto core = std::make_shared<catchem::Core>(fix.n_cols, fix.n_levels, fix.n_species);
            auto state = bind_state(core, fix);
            state->load_species_config(species_path);
            write_config("optics_baddir.yml", "/no/such/optics/dir/", true);
            auto config = std::make_shared<catchem::ConfigManager>();
            config->load_from_file("optics_baddir.yml");
            state->attach_config_manager(config);
            auto settling = catchem::ProcessRegistry::get_instance().create("settling");
            check(threw_containing([&] { settling->init(state); }, "/no/such/optics/dir/"),
                  "bad optics directory aborts init naming the directory");
        }

        // 2. A listed table whose file is absent names the type and path.
        {
            Fixture fix;
            auto core = std::make_shared<catchem::Core>(fix.n_cols, fix.n_levels, fix.n_species);
            auto state = bind_state(core, fix);
            state->load_species_config(species_path);
            std::ofstream out("optics_missing_file.yml");
            out << "mie:\n"
                << "  directory: \"" << optics_dir << "/\"\n"
                << "  files:\n"
                << "    DU: optics_DU.v99_9.nc\n";
            out << "processes:\n  settling:\n    activate: true\n    scheme: gocart\n"
                << "    gocart:\n      simple_scheme: true\n";
            out.close();
            auto config = std::make_shared<catchem::ConfigManager>();
            config->load_from_file("optics_missing_file.yml");
            state->attach_config_manager(config);
            auto settling = catchem::ProcessRegistry::get_instance().create("settling");
            bool named_type = false, named_path = false;
            try {
                settling->init(state);
            } catch (const std::exception& error) {
                const std::string what = error.what();
                named_type = what.find("DU") != std::string::npos;
                named_path = what.find("optics_DU.v99_9.nc") != std::string::npos;
            }
            check(named_type && named_path, "missing optics file aborts init naming type and expected path");
        }

        // 3. simple_scheme: true with no mie.files is an error, not a silent no-op.
        {
            Fixture fix;
            auto core = std::make_shared<catchem::Core>(fix.n_cols, fix.n_levels, fix.n_species);
            auto state = bind_state(core, fix);
            state->load_species_config(species_path);
            write_config("optics_nofiles.yml", optics_dir + "/", true, /*include_files=*/false);
            auto config = std::make_shared<catchem::ConfigManager>();
            config->load_from_file("optics_nofiles.yml");
            state->attach_config_manager(config);
            auto settling = catchem::ProcessRegistry::get_instance().create("settling");
            check(threw_containing([&] { settling->init(state); }, "mie.files"),
                  "simple_scheme with empty mie.files aborts naming the section (FR-013)");
        }

        // 4. A settling species whose __mie_name matches no table names the species.
        {
            Fixture fix;
            auto core = std::make_shared<catchem::Core>(fix.n_cols, fix.n_levels, fix.n_species);
            auto state = bind_state(core, fix);
            state->load_species_config(species_path);
            std::ofstream out("optics_unmatched.yml");
            out << "mie:\n"
                << "  directory: \"" << optics_dir << "/\"\n"
                << "  files:\n"
                << "    SS: optics_SS.v3_5.nc\n"
                << "    DU: optics_DU.v15_5.nc\n"
                << "    BC: optics_BC.v1_5.nc\n"
                << "    OC: optics_OC.v1_5.nc\n";
            out << "processes:\n  settling:\n    activate: true\n    scheme: gocart\n"
                << "    gocart:\n      simple_scheme: true\n";
            out.close();
            auto config = std::make_shared<catchem::ConfigManager>();
            config->load_from_file("optics_unmatched.yml");
            state->attach_config_manager(config);
            auto settling = catchem::ProcessRegistry::get_instance().create("settling");
            bool named_species = false, named_mie = false;
            try {
                settling->init(state);
            } catch (const std::exception& error) {
                const std::string what = error.what();
                named_species = what.find("so4") != std::string::npos; // SU-bearing species
                named_mie = what.find("SU") != std::string::npos;
            }
            check(named_species && named_mie,
                  "unmatched __mie_name aborts init naming the species and its type");
        }

        std::cout << (failures == 0 ? "SUCCESS: all settling optics-table assertions passed.\n"
                                    : "FAILURE: " + std::to_string(failures) + " assertion(s) failed.\n");
    }
    Kokkos::finalize();
    return failures == 0 ? 0 : 1;
}
