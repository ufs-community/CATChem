#include "catchem_api.hpp"
#include "catchem_config_manager.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_kokkos_compat.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_state_manager.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

extern "C" {
void catchem_register_drydep_cpp();
}

namespace {

    // Bit flags selecting which host fields to bind.  The negative cases drop
    // one contracted field at a time so the require_field_pointer guard for it
    // becomes observable.
    enum FieldFlag : unsigned {
        F_2D = 1u << 0,
        F_3D = 1u << 1,
        F_PEDGES = 1u << 2,
        F_Z = 1u << 3,
        F_CONC = 1u << 4,
        F_ALL = 0x1Fu,
    };

    struct Fixture {
        int n_cols = 4;
        int n_levels = 5;
        int n_species = 22;
        std::vector<double> lat, lon, ps, ustar, ts, pblh, z0h, hflux, obk, dluse, lai, frsno, swgdn, frlake, gwettop,
            lwi, u10m, v10m, z0, salinity, cldfrc, suncosmid;
        std::vector<double> temperature, airden, qv, pmid, cldf, bxheight, rh, delp;
        std::vector<double> pedge, zedge, chem_conc;

        Fixture()
            : lat(n_cols, 40.0), lon(n_cols, -100.0), ps(n_cols, 101325.0), ustar(n_cols, 0.5), ts(n_cols, 290.0),
              pblh(n_cols, 1000.0), z0h(n_cols, 0.01), hflux(n_cols, 50.0), obk(n_cols, 100.0), dluse(n_cols, 1.0),
              lai(n_cols, 1.0), frsno(n_cols, 0.0), swgdn(n_cols, 100.0), frlake(n_cols, 0.0), gwettop(n_cols, 0.2),
              lwi(n_cols, 1.0), u10m(n_cols, 3.0), v10m(n_cols, 1.0), z0(n_cols, 0.01), salinity(n_cols, 35.0),
              cldfrc(n_cols, 0.2), suncosmid(n_cols, 0.5), temperature(n_cols * n_levels, 290.0),
              airden(n_cols * n_levels, 1.2), qv(n_cols * n_levels, 0.01), pmid(n_cols * n_levels, 90000.0),
              cldf(n_cols * n_levels, 0.2), bxheight(n_cols * n_levels, 100.0), rh(n_cols * n_levels, 50.0),
              delp(n_cols * n_levels, 1000.0), pedge(n_cols * (n_levels + 1), 101300.0), zedge(n_cols * (n_levels + 1)),
              chem_conc(n_cols * n_levels * n_species, 1.0e-8) {
            // Geometric height interfaces ascending from the surface: the
            // GOCART aero scheme's hghte slot.  Feeding pressure here produced
            // NaN deposition velocities (spec 011 US2).
            for (int c = 0; c < n_cols; ++c)
                for (int lev = 0; lev <= n_levels; ++lev)
                    zedge[c + static_cast<size_t>(lev) * n_cols] = 100.0 * lev;
        }
    };

    std::shared_ptr<catchem::StateManager> bind_state(const std::shared_ptr<catchem::Core>& core, Fixture& fix,
                                                      unsigned flags) {
        auto state = core->get_state_manager();
        auto runtime_config = std::make_shared<catchem::ConfigManager>();
        runtime_config->load_from_file("CATChem_new_config.yml");
        // The Default inventory exercises the GOCART aerosol scheme, whose
        // hghte slot is the z_edges bridge argument.  The shared test config
        // may select zhang; force gocart so the US2 regression is observable.
        {
            YAML::Node settings;
            settings["gas_scheme"] = "wesely";
            settings["aero_scheme"] = "gocart";
            auto& proc = runtime_config->data.processes["drydep"];
            proc.activate = true;
            proc.diagnostics = true;
            proc.set_settings_node(settings);
        }
        state->attach_config_manager(runtime_config);

        // The parity oracle uses the Default 22-species inventory.
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
            throw std::runtime_error("drydep fixture species count mismatch: " + species_path);

        if (flags & F_2D) {
            state->bind_met_field_2d("LAT", fix.lat.data());
            state->bind_met_field_2d("LON", fix.lon.data());
            state->bind_met_field_2d("PS", fix.ps.data());
            state->bind_met_field_2d("USTAR", fix.ustar.data());
            state->bind_met_field_2d("TS", fix.ts.data());
            state->bind_met_field_2d("PBLH", fix.pblh.data());
            state->bind_met_field_2d("Z0H", fix.z0h.data());
            state->bind_met_field_2d("HFLUX", fix.hflux.data());
            state->bind_met_field_2d("OBK", fix.obk.data());
            state->bind_met_field_2d("DLUSE", fix.dluse.data());
            state->bind_met_field_2d("LAI", fix.lai.data());
            state->bind_met_field_2d("FRSNO", fix.frsno.data());
            state->bind_met_field_2d("SWGDN", fix.swgdn.data());
            state->bind_met_field_2d("FRLAKE", fix.frlake.data());
            state->bind_met_field_2d("GWETTOP", fix.gwettop.data());
            state->bind_met_field_2d("LWI", fix.lwi.data());
            state->bind_met_field_2d("U10M", fix.u10m.data());
            state->bind_met_field_2d("V10M", fix.v10m.data());
            state->bind_met_field_2d("Z0", fix.z0.data());
            state->bind_met_field_2d("SALINITY", fix.salinity.data());
            state->bind_met_field_2d("CLDFRC", fix.cldfrc.data());
            state->bind_met_field_2d("SUNCOSMID", fix.suncosmid.data());
        }
        if (flags & F_3D) {
            state->bind_met_field_3d("T", fix.temperature.data());
            state->bind_met_field_3d("AIRDEN", fix.airden.data());
            state->bind_met_field_3d("AIRDEN_DRY", fix.airden.data());
            state->bind_met_field_3d("QV", fix.qv.data());
            state->bind_met_field_3d("PMID", fix.pmid.data());
            state->bind_met_field_3d("CLDF", fix.cldf.data());
            state->bind_met_field_3d("BXHEIGHT", fix.bxheight.data());
            state->bind_met_field_3d("RH", fix.rh.data());
            state->bind_met_field_3d("DELP", fix.delp.data());
        }
        if (flags & F_PEDGES)
            state->bind_met_field_3d("PEDGE", fix.pedge.data());
        if (flags & F_Z)
            state->bind_met_field_3d("Z", fix.zedge.data());
        if (flags & F_CONC)
            state->bind_unified_chemistry(fix.chem_conc.data());
        return state;
    }

    bool threw_message(const std::function<void()>& fn, const std::string& expected) {
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
        std::cout << "RUNNING TEST: DryDep Process Unit Test" << std::endl;
        std::cout << "==========================================" << std::endl;

        catchem_register_drydep_cpp();
        assert(catchem::ProcessRegistry::get_instance().has_process("drydep"));

        // --- Contract requires geometric height Z (US2) ----------------------
        {
            auto drydep = catchem::ProcessRegistry::get_instance().create("drydep");
            auto contract = drydep->get_contract();
            bool z_required = false, pedge_required = false;
            for (const auto& field : contract.fields) {
                if (field.canonical_name == "Z" && field.requirement == catchem::FieldRequirement::Required)
                    z_required = true;
                if (field.canonical_name == "PEDGE")
                    pedge_required = true;
            }
            check(z_required, "contract requires Z (geometric height interface)");
            check(pedge_required, "contract still binds PEDGE for pressure consumers");
        }

        // --- Missing Z aborts with a named error ------------------------------
        {
            Fixture fix;
            auto core = std::make_shared<catchem::Core>(fix.n_cols, fix.n_levels, fix.n_species);
            bool named = threw_message(
                [&] {
                    auto state = bind_state(core, fix, F_ALL & ~FieldFlag::F_Z);
                    auto drydep = catchem::ProcessRegistry::get_instance().create("drydep");
                    drydep->prepare_inputs(state);
                    drydep->init(state);
                    drydep->run(state);
                },
                "DryDep process missing required field Z");
            check(named, "missing Z aborts with a named error");
        }

        // --- GOCART aerosol deposition velocities are finite ------------------
        // Regression for the PEDGE-in-z-slot bug: pressure (Pa) fed into the
        // scheme's hghte (m) slot produced NaN deposition velocities.
        {
            Fixture fix;
            auto core = std::make_shared<catchem::Core>(fix.n_cols, fix.n_levels, fix.n_species);
            auto state = bind_state(core, fix, F_ALL);
            state->clock().timestep = 3600.0;
            auto drydep = catchem::ProcessRegistry::get_instance().create("drydep");
            drydep->prepare_inputs(state);
            drydep->init(state);
            drydep->run(state);
            state->sync_to_host();

            auto* diag = state->diagnostic_manager().get();
            check(diag != nullptr && diag->has_field("drydep_velocity_per_species"),
                  "drydep_velocity_per_species diagnostic registered");
            bool all_finite = true;
            std::size_t nan_count = 0;
            if (diag && diag->has_field("drydep_velocity_per_species")) {
                const double* vel =
                    static_cast<const double*>(diag->get_host_read_pointer("drydep_velocity_per_species"));
                const std::size_t n = static_cast<std::size_t>(fix.n_cols) * fix.n_species;
                for (std::size_t i = 0; i < n; ++i) {
                    if (!std::isfinite(vel[i])) {
                        all_finite = false;
                        ++nan_count;
                    }
                }
            }
            check(all_finite, "all deposition velocities finite (no NaN from z-slot mixup)");
            if (!all_finite)
                std::cout << "  (non-finite velocity cells: " << nan_count << ")\n";

            // Concentrations must also stay finite after the scheme applies.
            const double* conc = fix.chem_conc.data();
            bool conc_finite = std::all_of(conc, conc + static_cast<size_t>(fix.n_cols) * fix.n_levels * fix.n_species,
                                           [](double v) { return std::isfinite(v); });
            check(conc_finite, "concentrations remain finite after drydep");
        }

        std::cout << (failures == 0 ? "SUCCESS: all drydep assertions passed.\n"
                                    : "FAILURE: " + std::to_string(failures) + " drydep assertion(s) failed.\n");
    }
    Kokkos::finalize();
    return failures == 0 ? 0 : 1;
}
