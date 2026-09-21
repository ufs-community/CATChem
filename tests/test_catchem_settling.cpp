#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_kokkos_compat.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_state_manager.hpp"
#include <cassert>
#include <cmath>
#include <cstring>
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

    // A small vertical profile that keeps every settling input finite and
    // physically ordered (temperature lapse, decreasing pressure, ascending
    // geometric height) so the GOCART2G kernel produces a genuine tendency.
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
                    RH[i] = 0.9 * std::exp(-lev / 5.0);
                }
            for (int c = 0; c < n_cols; ++c)
                for (int lev = 0; lev <= n_levels; ++lev)
                    Zedge[c + static_cast<size_t>(lev) * n_cols] = 100.0 * lev;
        }
    };

    // Bit flags selecting which host fields to bind.  AIRDEN/DELP/RH are bound
    // directly (their derivation prerequisites PEDGE/QV are intentionally left
    // unbound) so that prepare_inputs() short-circuits and each field's own
    // guard becomes observable.
    enum FieldFlag : unsigned {
        F_T = 1u << 0,
        F_AIRDEN = 1u << 1,
        F_DELP = 1u << 2,
        F_RH = 1u << 3,
        F_Z = 1u << 4,
        F_PMID = 1u << 5,
        F_CONC = 1u << 6,
        F_ALL = 0x7Fu,
    };

    std::shared_ptr<catchem::StateManager> bind_state(const std::shared_ptr<catchem::Core>& core, Fixture& fix,
                                                      unsigned flags) {
        auto state = core->get_state_manager();
        auto config = std::make_shared<catchem::ConfigManager>();
        config->load_from_file("CATChem_new_config.yml");
        state->attach_config_manager(config);

        // The parity oracle uses the Default 22-species inventory; the local
        // 49-species test file would change the aerosol set and slab count.
        const std::string species_name = "Configs/Default/CATChem_species.yml";
        std::string species_path = "CATChem_species.yml";
        for (const std::string& candidate : {species_name, "tests/" + species_name, "../tests/" + species_name,
                                             "../../tests/" + species_name, std::string("CATChem_species.yml")}) {
            if (std::ifstream(candidate).good()) {
                species_path = candidate;
                break;
            }
        }
        state->load_species_config(species_path);
        if (static_cast<int>(state->chemistry().species_list.size()) != fix.n_species)
            throw std::runtime_error("settling fixture species count mismatch: " + species_path);

        if (flags & F_T)
            state->bind_met_field_3d("T", fix.T.data());
        if (flags & F_AIRDEN)
            state->bind_met_field_3d("AIRDEN", fix.AIRDEN.data());
        if (flags & F_DELP)
            state->bind_met_field_3d("DELP", fix.DELP.data());
        if (flags & F_RH)
            state->bind_met_field_3d("RH", fix.RH.data());
        if (flags & F_Z)
            state->bind_met_field_3d("Z", fix.Zedge.data());
        if (flags & F_PMID)
            state->bind_met_field_3d("PMID", fix.PMID.data());
        if (flags & F_CONC) {
            fix.conc.assign(static_cast<size_t>(fix.n_cols) * fix.n_levels * fix.n_species, 1.0e-8);
            state->bind_unified_chemistry(fix.conc.data());
        }
        return state;
    }

    std::shared_ptr<catchem::ConfigManager> settling_config(bool maring_dust_only, bool correction_maring) {
        auto config = std::make_shared<catchem::ConfigManager>();
        config->load_from_file("CATChem_new_config.yml");
        YAML::Node gocart;
        gocart["scale_factor"] = 1.0;
        gocart["simple_scheme"] = false;
        gocart["swelling_rh_max"] = 0.95;
        gocart["correction_maring"] = correction_maring;
        gocart["maring_dust_only"] = maring_dust_only;
        YAML::Node settings;
        settings["gocart"] = gocart;
        auto& proc = config->data.processes["settling"];
        proc.activate = true;
        proc.scheme = "gocart";
        proc.set_settings_node(settings);
        return config;
    }

    int species_index(const std::shared_ptr<catchem::StateManager>& state, const std::string& name) {
        const auto& list = state->chemistry().species_list;
        for (std::size_t i = 0; i < list.size(); ++i)
            if (list[i].short_name == name)
                return static_cast<int>(i);
        return -1;
    }

    // Flat index into the bound (n_cols, n_levels, n_species) concentration
    // buffer: column fastest, then level, then species.
    std::size_t flat(const Fixture& fix, int column, int level, int species) {
        return static_cast<size_t>(species) * fix.n_levels * fix.n_cols + static_cast<size_t>(level) * fix.n_cols +
               column;
    }

    bool threw_message(const std::function<void()>& fn, const std::string& expected) {
        try {
            fn();
        } catch (const std::exception& error) {
            return std::string(error.what()).find(expected) != std::string::npos;
        }
        return false;
    }

    // Diagnostics: enabling them registers velocity/flux fields sized to the
    // aerosol subset, and the filled velocity must be non-negative, finite, and
    // (for a loaded column with mass) non-zero somewhere.
    void test_settling_diagnostics() {
        Fixture fix;
        auto core = std::make_shared<catchem::Core>(fix.n_cols, fix.n_levels, fix.n_species);
        auto state = bind_state(core, fix, FieldFlag::F_ALL);
        auto config = settling_config(true, false);
        // Turn diagnostics on for the settling process.  The shared test YAML
        // carries a settling diag_species subset; clear it so the default
        // (every aerosol) applies and the field dims match the aerosol subset.
        config->data.processes["settling"].diagnostics = true;
        config->data.processes["settling"].diag_species.clear();
        state->attach_config_manager(config);

        auto settling = catchem::ProcessRegistry::get_instance().create("settling");
        settling->init(state);
        settling->run(state);
        state->sync_to_host();

        const auto mgr = core->get_diagnostic_manager();
        const int n_aerosol = static_cast<int>(state->chemistry().aerosol_indices.size());
        assert(n_aerosol > 0);
        assert(mgr->has_field("settling_velocity_per_species_per_level"));
        assert(mgr->has_field("settling_flux_per_species"));
        assert(mgr->get_field("settling_velocity_per_species_per_level")->dimensions ==
               std::vector<int>({fix.n_cols, fix.n_levels, n_aerosol}));
        assert(mgr->get_field("settling_flux_per_species")->dimensions ==
               std::vector<int>({fix.n_cols, n_aerosol}));

        const double* vel =
            static_cast<const double*>(mgr->get_host_pointer("settling_velocity_per_species_per_level"));
        assert(vel != nullptr);
        const int n = fix.n_cols * fix.n_levels * n_aerosol;
        for (int i = 0; i < n; ++i) {
            assert(std::isfinite(vel[i]));
            assert(vel[i] >= 0.0); // GOCART clamps negative vsettle to zero
        }
        bool any_positive = false;
        for (int i = 0; i < n; ++i)
            if (vel[i] > 0.0) { any_positive = true; break; }
        assert(any_positive && "settling velocity diagnostics must be populated by the bridge");
        std::cout << "  PASS settling_diagnostics: fields present, velocity finite & >= 0" << std::endl;
    }

    // Regression: enabling diagnostics must not perturb the science.  Run the
    // identical fixture twice — once with diagnostics off, once on — and assert
    // the resulting concentration buffers are bit-identical.
    void test_settling_diag_off_matches() {
        auto run_once = [](bool diagnostics, std::vector<double>& conc_out) {
            Fixture fix;
            auto core = std::make_shared<catchem::Core>(fix.n_cols, fix.n_levels, fix.n_species);
            auto state = bind_state(core, fix, FieldFlag::F_ALL);
            auto config = settling_config(true, false);
            config->data.processes["settling"].diagnostics = diagnostics;
            state->attach_config_manager(config);
            state->clock().timestep = 3600.0;

            auto settling = catchem::ProcessRegistry::get_instance().create("settling");
            settling->init(state);
            settling->run(state);
            state->sync_to_host();
            conc_out = fix.conc;
        };

        std::vector<double> off, on;
        run_once(false, off);
        run_once(true, on);
        assert(off.size() == on.size() && !off.empty());
        for (size_t i = 0; i < off.size(); ++i) {
            // Bit-identical: diagnostics are a read-only side product of the
            // same compute_gocart call, so any difference is a wiring bug.
            assert(std::memcmp(&off[i], &on[i], sizeof(double)) == 0);
        }
        std::cout << "  PASS settling_diag_off_matches: tendencies bit-identical with diagnostics on" << std::endl;
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
        std::cout << "RUNNING TEST: Settling Process Unit Test" << std::endl;
        std::cout << "==========================================" << std::endl;

        catchem_register_settling_cpp();
        assert(catchem::ProcessRegistry::get_instance().has_process("settling"));

        // --- Contract now requires PMID (US1) --------------------------------
        {
            auto settling = catchem::ProcessRegistry::get_instance().create("settling");
            auto contract = settling->get_contract();
            for (const char* name : {"T", "AIRDEN", "DELP", "RH", "Z", "PMID"}) {
                bool present = false;
                for (const auto& field : contract.fields) {
                    if (field.canonical_name == name) {
                        present = true;
                        check(field.execution_space == catchem::ExecutionSpaceIntent::Host,
                              std::string("Fortran settling bridge declares host access for ") + name);
                    }
                }
                check(present, std::string("contract requires ") + name);
            }
            for (const auto& field : contract.fields)
                if (field.canonical_name == "CONCENTRATION")
                    check(field.execution_space == catchem::ExecutionSpaceIntent::Host,
                          "Fortran settling bridge declares host access for concentration");
        }

        // --- Absent contracted field aborts with a named error ---------------
        // Host-authoritative inputs (T, Z, PMID, CHEM_CONC) trip the strict
        // require_field_pointer guard; the derived inputs (AIRDEN, DELP, RH)
        // trip their derivation guard, which also names the field.
        struct Case {
            unsigned keep;
            std::string expected;
            std::string label;
        };
        const Case cases[] = {
            {FieldFlag::F_ALL & ~FieldFlag::F_T, "Settling process missing required field T", "T"},
            {FieldFlag::F_ALL & ~FieldFlag::F_Z, "Settling process missing required field Z", "Z"},
            {FieldFlag::F_ALL & ~FieldFlag::F_PMID, "Settling process missing required field PMID", "PMID"},
            {FieldFlag::F_ALL & ~FieldFlag::F_CONC, "Settling process missing required field CHEM_CONC", "CHEM_CONC"},
            {FieldFlag::F_ALL & ~(FieldFlag::F_AIRDEN | FieldFlag::F_PMID), "derive AIRDEN", "AIRDEN"},
            {FieldFlag::F_ALL & ~FieldFlag::F_DELP, "derive DELP", "DELP"},
            {FieldFlag::F_ALL & ~FieldFlag::F_RH, "derive RH", "RH"},
        };
        for (const auto& test : cases) {
            Fixture fix;
            auto core = std::make_shared<catchem::Core>(fix.n_cols, fix.n_levels, fix.n_species);
            bool named = threw_message(
                [&] {
                    auto state = bind_state(core, fix, test.keep);
                    auto settling = catchem::ProcessRegistry::get_instance().create("settling");
                    settling->init(state);
                    settling->run(state);
                },
                test.expected);
            check(named, "missing " + test.label + " aborts with a named error");
        }

        // --- Dust-only Maring gating + µm radius scale -----------------------
        // correction_maring=true with maring_dust_only=true leaves sea salt
        // uncorrected; with maring_dust_only=false sea salt is corrected.
        // is corrected in both, so its tendency is identical.  Settling is
        // mass-conserving within a column, so the probe compares the surface
        // cell (which gains flux from above), not the column sum.  The coarse
        // seas5 bin (r≈7.8 µm) settles fast enough for a measurable surface
        // gain; the Maring 0.0033 m/s subtraction changes that gain.
        auto run_variant = [&](Fixture& fix, bool maring_dust_only, bool correction_maring, double dt) {
            auto core = std::make_shared<catchem::Core>(fix.n_cols, fix.n_levels, fix.n_species);
            auto state = bind_state(core, fix, FieldFlag::F_ALL);
            state->attach_config_manager(settling_config(maring_dust_only, correction_maring));
            state->clock().timestep = dt;
            auto settling = catchem::ProcessRegistry::get_instance().create("settling");
            settling->init(state);
            settling->run(state);
            state->sync_to_host();
        };

        Fixture dust_only, all_corrected, probe;
        auto probe_core = std::make_shared<catchem::Core>(probe.n_cols, probe.n_levels, probe.n_species);
        auto probe_state = bind_state(probe_core, probe, FieldFlag::F_ALL);
        int dust = species_index(probe_state, "dust1");
        int seas = species_index(probe_state, "seas5");
        check(dust >= 0 && seas >= 0, "dust1/seas5 present in Default species list");

        run_variant(dust_only, /*maring_dust_only=*/true, /*correction_maring=*/true, /*dt=*/3600.0);
        run_variant(all_corrected, /*maring_dust_only=*/false, /*correction_maring=*/true, /*dt=*/3600.0);

        double top_loss = 0.0, seas_diff = 0.0, dust_diff = 0.0;
        const int top = probe.n_levels - 1;
        for (int c = 0; c < probe.n_cols; ++c) {
            top_loss += 1.0e-8 - all_corrected.conc[flat(probe, c, top, seas)];
            seas_diff +=
                std::abs(all_corrected.conc[flat(probe, c, top, seas)] - dust_only.conc[flat(probe, c, top, seas)]);
            dust_diff +=
                std::abs(all_corrected.conc[flat(probe, c, top, dust)] - dust_only.conc[flat(probe, c, top, dust)]);
        }
        check(std::isfinite(top_loss) && top_loss > 0.0,
              "seas5 depletes at the top of the column (radius reached kernel in µm)");
        check(seas_diff > 0.0, "sea-salt tendency changes with maring_dust_only");
        check(dust_diff == 0.0, "dust tendency is unchanged by maring_dust_only (always corrected)");

        test_settling_diagnostics();
        test_settling_diag_off_matches();

        std::cout << (failures == 0 ? "SUCCESS: all settling assertions passed.\n"
                                    : "FAILURE: " + std::to_string(failures) + " settling assertion(s) failed.\n");
    }
    Kokkos::finalize();
    return failures == 0 ? 0 : 1;
}
