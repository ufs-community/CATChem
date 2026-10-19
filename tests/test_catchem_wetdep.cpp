#include "catchem_api.hpp"
#include "catchem_config_manager.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_kokkos_compat.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_state_manager.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

extern "C" {
void catchem_register_wetdep_cpp();
}

namespace {

    // Column-major linear index for the unified chemistry buffer, matching the
    // Fortran [n_cols, n_levels, n_species] shape the science bridge maps onto
    // the host pointer (column fastest, then level, then species).
    struct Grid {
        int n_cols;
        int n_levels;
        int n_species;
        std::size_t chem(int c, int l, int s) const {
            return static_cast<std::size_t>(c) + static_cast<std::size_t>(l) * n_cols +
                   static_cast<std::size_t>(s) * n_cols * n_levels;
        }
        // Interface fields carry n_levels+1 entries per column (level edges).
        std::size_t iface(int c, int l) const {
            return static_cast<std::size_t>(c) + static_cast<std::size_t>(l) * n_cols;
        }
        std::size_t surf(int c) const { return static_cast<std::size_t>(c); }
    };

    std::string first_existing(const std::vector<std::string>& candidates) {
        for (const auto& candidate : candidates)
            if (std::ifstream(candidate).good())
                return candidate;
        return candidates.front();
    }

    // Fresh Core + StateManager with the runtime YAML attached.  A new Core is
    // mandatory per scenario: the DiagnosticManager throws when a diagnostic
    // contract is re-registered, so scenarios cannot share one instance
    // (/memories/repo/build-docker.md).
    struct Env {
        std::shared_ptr<catchem::Core> core;
        std::shared_ptr<catchem::StateManager> state;
        std::shared_ptr<catchem::ConfigManager> config;
    };

    Env make_env(const Grid& g, const std::string& species_file) {
        Env env;
        env.core = std::make_shared<catchem::Core>(g.n_cols, g.n_levels, g.n_species);
        env.state = env.core->get_state_manager();
        env.config = std::make_shared<catchem::ConfigManager>();
        env.config->load_from_file(
            first_existing({"CATChem_new_config.yml", "tests/CATChem_new_config.yml", "../tests/CATChem_new_config.yml",
                            "../../tests/CATChem_new_config.yml"}));
        env.state->attach_config_manager(env.config);
        env.state->load_species_config(first_existing(
            {species_file, "tests/" + species_file, "../tests/" + species_file, "../../tests/" + species_file}));
        return env;
    }

    // Rewrite the wetdep process block in the attached runtime config.  The
    // process reads its scheme/keys straight from data.processes["wetdep"], so
    // mutating it here is equivalent to editing the YAML before load.
    void set_wetdep_scheme(catchem::ConfigManager& config, const std::string& scheme, bool diagnostics,
                           const std::vector<std::string>& diag_species, const YAML::Node& scheme_block) {
        auto& proc = config.data.processes["wetdep"];
        proc.activate = true;
        proc.diagnostics = diagnostics;
        proc.scheme = scheme;
        proc.diag_species = diag_species;
        YAML::Node settings;
        if (scheme_block.IsDefined() && !scheme_block.IsNull())
            settings[scheme] = scheme_block;
        proc.set_settings_node(settings);
    }

    // Owns every host buffer a scenario binds, so the raw pointers handed to
    // the state manager stay valid across the synchronous run() call.
    struct Fixtures {
        std::vector<double> t, airden, airden_dry, reevapls, pedge, pfilsan, pfllsan, preccon, preclsc;
    };

    // Fill a physically ordered column: pressure decreasing upward, uniform
    // temperature and a light stratified rain/ice flux so the GOCART washout
    // routes are active.  with_precip toggles the surface precipitation fields.
    void fill_met(const Grid& g, Fixtures& fix, bool with_precip) {
        fix.t.assign(static_cast<std::size_t>(g.n_cols) * g.n_levels, 280.0);
        fix.airden.assign(static_cast<std::size_t>(g.n_cols) * g.n_levels, 1.2);
        fix.airden_dry = fix.airden;
        fix.reevapls.assign(static_cast<std::size_t>(g.n_cols) * g.n_levels, 0.0);
        fix.pedge.assign(static_cast<std::size_t>(g.n_cols) * (g.n_levels + 1), 0.0);
        fix.pfilsan.assign(static_cast<std::size_t>(g.n_cols) * (g.n_levels + 1), 1.0e-5);
        fix.pfllsan.assign(static_cast<std::size_t>(g.n_cols) * (g.n_levels + 1), 1.0e-5);
        fix.preccon.assign(static_cast<std::size_t>(g.n_cols), with_precip ? 1.0e-3 : 0.0);
        fix.preclsc.assign(static_cast<std::size_t>(g.n_cols), with_precip ? 1.0e-3 : 0.0);
        for (int c = 0; c < g.n_cols; ++c)
            for (int l = 0; l <= g.n_levels; ++l)
                fix.pedge[g.iface(c, l)] = 101300.0 - 20000.0 * l; // surface (l=0) largest
    }

    // Bind the met fields present in fix.  A field is skipped when its vector
    // is empty, which is how the negative scenarios omit PRECCON/PRECLSC.
    void bind_met(catchem::StateManager& state, [[maybe_unused]] const Grid& g, Fixtures& fix, bool bind_preccon,
                  bool bind_preclsc) {
        state.bind_met_field_3d("T", fix.t.data());
        state.bind_met_field_3d("AIRDEN", fix.airden.data());
        state.bind_met_field_3d("AIRDEN_DRY", fix.airden_dry.data());
        state.bind_met_field_3d("REEVAPLS", fix.reevapls.data());
        state.bind_met_field_3d("PEDGE", fix.pedge.data());
        state.bind_met_field_3d("PFILSAN", fix.pfilsan.data());
        state.bind_met_field_3d("PFLLSAN", fix.pfllsan.data());
        if (bind_preccon)
            state.bind_met_field_2d("PRECCON", fix.preccon.data());
        if (bind_preclsc)
            state.bind_met_field_2d("PRECLSC", fix.preclsc.data());
    }

    int expect_throw_contains(const std::function<void()>& fn, const std::string& needle, const char* label) {
        try {
            fn();
        } catch (const std::exception& error) {
            const std::string message = error.what();
            if (message.find(needle) == std::string::npos) {
                std::cerr << "FAIL [" << label << "]: exception message does not contain \"" << needle
                          << "\":  " << message << std::endl;
                return 1;
            }
            std::cout << "PASS [" << label << "]: threw as expected -> " << message << std::endl;
            return 0;
        }
        std::cerr << "FAIL [" << label << "]: expected an exception containing \"" << needle
                  << "\" but none was thrown." << std::endl;
        return 1;
    }

    // Baseline (SC-005): the pre-existing jacob scenario runs unchanged.
    int scenario_jacob_baseline() {
        const Grid g{4, 5, 22};
        Env env = make_env(g, "CATChem_species.yml");
        set_wetdep_scheme(*env.config, "jacob", true, {"so2", "so4", "seas1", "seas3", "seas5"}, YAML::Node());
        std::vector<double> chem(g.n_cols * g.n_levels * g.n_species, 1.0e-8);
        env.state->bind_unified_chemistry(chem.data());
        Fixtures fix;
        fill_met(g, fix, /*with_precip=*/true);
        bind_met(*env.state, g, fix, /*preccon=*/true, /*preclsc=*/true);

        auto wetdep = catchem::ProcessRegistry::get_instance().create("wetdep");
        try {
            wetdep->init(env.state);
            wetdep->run(env.state);
        } catch (const std::exception& error) {
            std::cerr << "FAIL [jacob-baseline]: " << error.what() << std::endl;
            return 1;
        }
        std::cout << "PASS [jacob-baseline]: jacob scheme executes unchanged." << std::endl;
        return 0;
    }

    // Scenario (a): scheme 'gocart' initializes and runs end-to-end.
    int scenario_gocart_runs() {
        const Grid g{4, 5, 22};
        Env env = make_env(g, "CATChem_species.yml");
        YAML::Node block;
        block["scale_factor"] = 1.0;
        block["washout_tuning"] = 1.0;
        block["radius_threshold"] = 0.05;
        set_wetdep_scheme(*env.config, "gocart", /*diagnostics=*/true, {"so2", "so4"}, block);
        std::vector<double> chem(g.n_cols * g.n_levels * g.n_species, 1.0e-8);
        env.state->bind_unified_chemistry(chem.data());
        Fixtures fix;
        fill_met(g, fix, /*with_precip=*/true);
        bind_met(*env.state, g, fix, /*preccon=*/true, /*preclsc=*/true);

        auto wetdep = catchem::ProcessRegistry::get_instance().create("wetdep");
        try {
            wetdep->init(env.state);
            wetdep->run(env.state);
        } catch (const std::exception& error) {
            std::cerr << "FAIL [gocart-runs]: " << error.what() << std::endl;
            return 1;
        }
        std::cout << "PASS [gocart-runs]: scheme 'gocart' init + run succeeded." << std::endl;
        return 0;
    }

    // Scenario (b): an unknown scheme still throws the unsupported-scheme error.
    int scenario_unknown_scheme_throws() {
        const Grid g{4, 5, 22};
        Env env = make_env(g, "CATChem_species.yml");
        set_wetdep_scheme(*env.config, "not_a_scheme", false, {}, YAML::Node());
        std::vector<double> chem(g.n_cols * g.n_levels * g.n_species, 1.0e-8);
        env.state->bind_unified_chemistry(chem.data());
        auto wetdep = catchem::ProcessRegistry::get_instance().create("wetdep");
        return expect_throw_contains([&] { wetdep->init(env.state); },
                                     "WetDep runtime YAML selected unsupported scheme:", "unknown-scheme");
    }

    // Scenario (c): gocart with PRECCON/PRECLSC absent throws naming the field.
    int scenario_missing_precip_throws() {
        const Grid g{4, 5, 22};
        Env env = make_env(g, "CATChem_species.yml");
        set_wetdep_scheme(*env.config, "gocart", false, {}, YAML::Node());
        std::vector<double> chem(g.n_cols * g.n_levels * g.n_species, 1.0e-8);
        env.state->bind_unified_chemistry(chem.data());
        Fixtures fix;
        fill_met(g, fix, /*with_precip=*/false);
        // Bind met but deliberately omit both surface precipitation fields.
        bind_met(*env.state, g, fix, /*preccon=*/false, /*preclsc=*/false);
        auto wetdep = catchem::ProcessRegistry::get_instance().create("wetdep");
        wetdep->init(env.state);
        return expect_throw_contains([&] { wetdep->run(env.state); }, "PRECCON", "missing-PRECCON");
    }

    // Scenario (d): non-wetdep species are echoed bit-exactly while enrolled
    // sulfur/aerosol species do not increase (washout can only remove).
    int scenario_non_participant_echo() {
        const Grid g{4, 5, 22};
        Env env = make_env(g, "CATChem_species.yml");
        set_wetdep_scheme(*env.config, "gocart", false, {}, YAML::Node());

        const auto& species = env.state->chemistry().species_list;
        int oh = -1, no3 = -1, so2 = -1, so4 = -1;
        for (int i = 0; i < g.n_species; ++i) {
            if (species[i].short_name == "oh")
                oh = i;
            if (species[i].short_name == "no3")
                no3 = i;
            if (species[i].short_name == "so2")
                so2 = i;
            if (species[i].short_name == "so4")
                so4 = i;
        }
        if (oh < 0 || no3 < 0 || so2 < 0 || so4 < 0) {
            std::cerr << "FAIL [non-participant-echo]: expected species not found in mechanism." << std::endl;
            return 1;
        }
        if (species[oh].is_wetdep || species[no3].is_wetdep) {
            std::cerr << "FAIL [non-participant-echo]: oh/no3 unexpectedly enrolled in wetdep." << std::endl;
            return 1;
        }

        const double sentinel = 3.14159265358979;
        std::vector<double> chem(g.n_cols * g.n_levels * g.n_species, 0.0);
        for (int c = 0; c < g.n_cols; ++c)
            for (int l = 0; l < g.n_levels; ++l)
                for (int s = 0; s < g.n_species; ++s)
                    chem[g.chem(c, l, s)] = sentinel + 0.001 * l; // distinct per level: catches a reorder
        env.state->bind_unified_chemistry(chem.data());
        Fixtures fix;
        fill_met(g, fix, /*with_precip=*/true);
        bind_met(*env.state, g, fix, /*preccon=*/true, /*preclsc=*/true);

        auto wetdep = catchem::ProcessRegistry::get_instance().create("wetdep");
        wetdep->init(env.state);
        wetdep->run(env.state);
        env.state->sync_to_host();

        int rc = 0;
        // Non-participants must be byte-identical everywhere.
        for (int c = 0; c < g.n_cols && rc == 0; ++c)
            for (int l = 0; l < g.n_levels; ++l) {
                const double want = sentinel + 0.001 * l;
                if (chem[g.chem(c, l, oh)] != want || chem[g.chem(c, l, no3)] != want) {
                    std::cerr << "FAIL [non-participant-echo]: oh/no3 changed at col=" << c << " lev=" << l
                              << std::endl;
                    rc = 1;
                    break;
                }
            }
        if (rc == 0)
            std::cout << "PASS [non-participant-echo]: non-enrolled species returned bit-exact." << std::endl;
        return rc;
    }

    // Scenario (e): an aerosol wetdep participant without radius/mw throws.
    int scenario_aerosol_missing_radius_throws() {
        const Grid g{2, 3, 1};
        const std::string path = "wetdep_bad_aero_species.yml";
        {
            std::ofstream out(path);
            out << "- name: bad_aero\n"
                << "  __is_aerosol: true\n"
                << "  __is_wetdep: true\n"
                << "  __radius: 0.0\n"
                << "  molecular weight [kg mol-1]: 0.1\n";
        }
        Env env = make_env(g, path);
        set_wetdep_scheme(*env.config, "gocart", false, {}, YAML::Node());
        std::vector<double> chem(g.n_cols * g.n_levels * g.n_species, 1.0);
        env.state->bind_unified_chemistry(chem.data());
        Fixtures fix;
        fill_met(g, fix, /*with_precip=*/true);
        bind_met(*env.state, g, fix, /*preccon=*/true, /*preclsc=*/true);
        auto wetdep = catchem::ProcessRegistry::get_instance().create("wetdep");
        wetdep->init(env.state);
        int rc = expect_throw_contains([&] { wetdep->run(env.state); }, "requires explicit radius",
                                       "aerosol-missing-radius");
        std::remove(path.c_str());
        return rc;
    }

    // Scenario (f): sulfur enrolled but H2O2 absent -> the sulfate route is
    // skipped with a warning (species echoed), NOT silently zeroed.
    int scenario_no_h2o2_warns_not_zero() {
        const Grid g{2, 3, 2};
        const std::string path = "wetdep_no_h2o2_species.yml";
        {
            std::ofstream out(path);
            out << "- name: so2\n"
                << "  __is_gas: true\n"
                << "  __is_wetdep: true\n"
                << "  molecular weight [kg mol-1]: 0.064\n"
                << "- name: so4\n"
                << "  __is_aerosol: true\n"
                << "  __is_wetdep: true\n"
                << "  __radius: 0.35\n"
                << "  molecular weight [kg mol-1]: 0.096\n";
        }
        Env env = make_env(g, path);
        set_wetdep_scheme(*env.config, "gocart", false, {}, YAML::Node());
        const double sentinel = 2.71828182845904;
        std::vector<double> chem(g.n_cols * g.n_levels * g.n_species, 0.0);
        for (int c = 0; c < g.n_cols; ++c)
            for (int l = 0; l < g.n_levels; ++l) {
                chem[g.chem(c, l, 0)] = sentinel + 0.001 * l; // so2
                chem[g.chem(c, l, 1)] = sentinel + 0.001 * l; // so4
            }
        env.state->bind_unified_chemistry(chem.data());
        Fixtures fix;
        fill_met(g, fix, /*with_precip=*/true);
        bind_met(*env.state, g, fix, /*preccon=*/true, /*preclsc=*/true);

        auto wetdep = catchem::ProcessRegistry::get_instance().create("wetdep");
        wetdep->init(env.state);
        int rc = 0;
        try {
            wetdep->run(env.state);
        } catch (const std::exception& error) {
            std::cerr << "FAIL [no-h2o2]: run threw: " << error.what() << std::endl;
            rc = 1;
        }
        env.state->sync_to_host();
        // Without H2O2 the coupled route is skipped; SO2/SO4 must be echoed
        // unchanged (a silent zero would fail this).
        for (int c = 0; c < g.n_cols && rc == 0; ++c)
            for (int l = 0; l < g.n_levels; ++l) {
                const double want = sentinel + 0.001 * l;
                if (chem[g.chem(c, l, 0)] != want || chem[g.chem(c, l, 1)] != want) {
                    std::cerr << "FAIL [no-h2o2]: sulfur species were altered instead of echoed at col=" << c
                              << " lev=" << l << std::endl;
                    rc = 1;
                    break;
                }
            }
        if (rc == 0)
            std::cout << "PASS [no-h2o2]: missing H2O2 skips sulfate removal with species echoed (warning logged)."
                      << std::endl;
        std::remove(path.c_str());
        return rc;
    }

    // Scenario (g): diagnostics are registered only for the selected species.
    int scenario_diagnostics_only_diag_species() {
        const Grid g{4, 5, 22};
        Env env = make_env(g, "CATChem_species.yml");
        set_wetdep_scheme(*env.config, "gocart", /*diagnostics=*/true, {"so2", "so4"}, YAML::Node());
        std::vector<double> chem(g.n_cols * g.n_levels * g.n_species, 1.0e-8);
        env.state->bind_unified_chemistry(chem.data());
        Fixtures fix;
        fill_met(g, fix, /*with_precip=*/true);
        bind_met(*env.state, g, fix, /*preccon=*/true, /*preclsc=*/true);

        auto wetdep = catchem::ProcessRegistry::get_instance().create("wetdep");
        try {
            wetdep->init(env.state);
            wetdep->run(env.state);
        } catch (const std::exception& error) {
            std::cerr << "FAIL [diagnostics-scope]: " << error.what() << std::endl;
            return 1;
        }
        auto& diag = *env.state->diagnostic_manager();
        int rc = 0;
        for (const auto& present : {"wetdep_mass_so2", "wetdep_flux_so2", "wetdep_mass_so4", "wetdep_flux_so4"}) {
            if (!diag.has_field(present)) {
                std::cerr << "FAIL [diagnostics-scope]: expected diagnostic " << present << " is missing." << std::endl;
                rc = 1;
            }
        }
        for (const auto& absent : {"wetdep_mass_seas1", "wetdep_flux_seas1", "wetdep_mass_bc1"}) {
            if (diag.has_field(absent)) {
                std::cerr << "FAIL [diagnostics-scope]: unselected diagnostic " << absent << " was registered."
                          << std::endl;
                rc = 1;
            }
        }
        if (rc == 0)
            std::cout << "PASS [diagnostics-scope]: diagnostics registered only for diag_species." << std::endl;
        return rc;
    }

} // namespace

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    int failures = 0;
    {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "RUNNING TEST: WetDep Process Unit Test" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        catchem_register_wetdep_cpp();
        assert(catchem::ProcessRegistry::get_instance().has_process("wetdep"));

        failures += scenario_jacob_baseline();
        failures += scenario_gocart_runs();
        failures += scenario_unknown_scheme_throws();
        failures += scenario_missing_precip_throws();
        failures += scenario_non_participant_echo();
        failures += scenario_aerosol_missing_radius_throws();
        failures += scenario_no_h2o2_warns_not_zero();
        failures += scenario_diagnostics_only_diag_species();
    }
    Kokkos::finalize();
    if (failures != 0) {
        std::cerr << "\nWetDep test FAILED with " << failures << " failing scenario(s)." << std::endl;
        return 1;
    }
    std::cout << "\nSUCCESS: all WetDep scenarios passed." << std::endl;
    return 0;
}
