#include "catchem_process_settling.hpp"
#include "catchem_error.hpp"
#include "catchem_logger.hpp"
#include "catchem_process_registry.hpp"
#include <iostream>
#include <stdexcept>

namespace catchem {

    extern "C" void run_settling_science_bridge(int n_columns, int n_levels, int n_aerosols, int n_total_species,
                                                double dt, double scale_factor, double swelling_rh_max,
                                                int correction_maring, int maring_dust_only, double* airden,
                                                double* delp, const double* pmid, double* rh, double* temperature,
                                                double* z_edge, const char* aerosol_species_names,
                                                const char* species_names, const int* species_is_dust,
                                                const int* species_is_hydrophilic, const double* radius,
                                                const double* density, double* concentration, int* bridge_rc);

    ProcessContract SettlingProcess::get_contract() const {
        // The current GOCART implementation invokes a Fortran science bridge
        // with host pointers.  Advertising these fields as device accesses
        // makes ExecutionPlan::complete() mark the device copy as the latest
        // writer after the bridge has modified host concentration memory.
        // A later host synchronization can then overwrite the settling result
        // (and the concentration presented to a coupled host) with stale data.
        return make_contract(get_name(),
                             {host_field_3d("T", "K"), host_field_3d("AIRDEN", "kg/m3"), host_field_3d("DELP", "Pa"),
                              host_field_3d("PMID", "Pa"), host_field_3d("RH", "1"), host_field_interface("Z", "m"),
                              host_concentration()});
    }

    SettlingProcess::SettlingProcess() : active_scheme("c++_kokkos"), fortran_callback(nullptr) {}

    void SettlingProcess::prepare_inputs(std::shared_ptr<StateManager> state) {
        // DELP, AIRDEN, and RH are optional host products. GOCART settling
        // requires all three, so make them current before the execution plan
        // validates this process contract.
        state->derive_delp();
        state->derive_airden();
        state->derive_relative_humidity();
    }

    void SettlingProcess::init(std::shared_ptr<StateManager> state) {
        const auto config = state->config_manager();
        if (!config)
            throw std::invalid_argument("Settling requires a runtime YAML configuration");
        const auto configured = config->data.processes.find("settling");
        if (configured == config->data.processes.end() || configured->second.scheme != "gocart")
            throw std::invalid_argument("Settling requires processes.settling.scheme: gocart");
        active_scheme = configured->second.scheme;

        // Read scheme tuning options from the runtime YAML.  Each lookup falls
        // back to the compiled default declared in SettlingCommon_Mod.F90.
        gocart_scale_factor = configured->second.get_double("gocart/scale_factor", gocart_scale_factor);
        gocart_simple_scheme = configured->second.get_bool("gocart/simple_scheme", gocart_simple_scheme);
        gocart_swelling_rh_max = configured->second.get_double("gocart/swelling_rh_max", gocart_swelling_rh_max);
        gocart_correction_maring = configured->second.get_bool("gocart/correction_maring", gocart_correction_maring);
        gocart_maring_dust_only = configured->second.get_bool("gocart/maring_dust_only", gocart_maring_dust_only);
        if (!(gocart_scale_factor > 0.0))
            throw std::invalid_argument("Settling gocart scale_factor must be positive");
        // A non-positive cap disables the clamp; otherwise it must be a valid RH fraction.
        if (gocart_swelling_rh_max > 1.0)
            throw std::invalid_argument(
                "Settling gocart swelling_rh_max must be <= 1.0 (RH fraction), or <= 0 to disable");

        // Surface the effective scheme options so the run log confirms what
        // was parsed from the runtime YAML and reaches the settling kernel.
        // The metadata path corresponds to legacy simple_scheme: false.
        // Mie-table settling remains intentionally unsupported by C++.
        if (gocart_simple_scheme)
            throw std::invalid_argument(
                "Settling simple_scheme requires Mie tables and is unsupported by the C++ core");
        Logger::debug(state.get(), "Settling scheme options",
                      {{"scheme", active_scheme},
                       {"gocart/scale_factor", std::to_string(gocart_scale_factor)},
                       {"gocart/correction_maring", gocart_correction_maring ? "true" : "false"},
                       {"gocart/maring_dust_only", gocart_maring_dust_only ? "true" : "false"},
                       {"gocart/simple_scheme", "false (species metadata)"},
                       {"gocart/swelling", "per-species __hydrophilic (Gerber when true)"},
                       {"gocart/swelling_rh_max", std::to_string(gocart_swelling_rh_max)}});

        int num_aerosols = state->chemistry().aerosol_indices.size();
        if (num_aerosols > 0) {
            aerosol_species_names.assign(static_cast<size_t>(num_aerosols) * 32, ' ');
            host_radius_dry.assign(num_aerosols, 0.0);
            host_rhop_dry.assign(num_aerosols, 0.0);
            host_is_dust.assign(num_aerosols, 0);
            host_is_hydrophilic.assign(num_aerosols, 1);

            for (int i = 0; i < num_aerosols; ++i) {
                int ispec = state->chemistry().aerosol_indices[i];
                double r_val = state->chemistry().species_list[ispec].radius;
                double d_val = state->chemistry().species_list[ispec].density;
                if (!(r_val > 0.0 && d_val > 0.0))
                    throw std::runtime_error("Settling aerosol '" + state->chemistry().species_list[ispec].short_name +
                                             "' requires explicit radius and density");
                // Radii are configured in micrometres and cross the bridge in
                // micrometres; the legacy scheme performs the µm -> m conversion.
                host_radius_dry[i] = r_val;
                host_rhop_dry[i] = d_val;
                host_is_dust[i] = state->chemistry().species_list[ispec].is_dust ? 1 : 0;
                // Per-species hygroscopicity drives wet-particle swelling: a
                // hydrophilic aerosol grows with RH (Gerber), a hydrophobic one
                // settles at its dry size.  Replaces the old global
                // swelling_method knob.
                host_is_hydrophilic[i] = state->chemistry().species_list[ispec].is_hydrophilic ? 1 : 0;
                std::copy_n(state->chemistry().species_names_c_arr.data() + static_cast<size_t>(ispec) * 32, 32,
                            aerosol_species_names.data() + static_cast<size_t>(i) * 32);
            }
        }
    }

    void SettlingProcess::set_fortran_bridge_callback(std::function<void(void*)> cb) {
        fortran_callback = cb;
    }

    void SettlingProcess::run(std::shared_ptr<StateManager> state) {
        if (fortran_callback) {
            // Fallback for tests explicitly requesting the Fortran bridge
            fortran_callback(static_cast<void*>(state.get()));
            if (state->chemistry().conc)
                state->chemistry().conc->mark_host_modified();
            return;
        }

        // The execution plan invokes prepare_inputs before run(), but direct
        // API users and focused process tests may call run() themselves.
        // These derivations are generation-aware and therefore preserve
        // host-provided fields while supplying only absent prerequisites.
        prepare_inputs(state);

        int num_aerosols = state->chemistry().aerosol_indices.size();
        if (num_aerosols == 0) {
            Logger::info(state.get(), "Settling skipped: no aerosol species registered", {});
            return;
        }

        require_field_pointer("Settling", "T", state->meteorology().T ? state->meteorology().T->host_data() : nullptr);
        require_field_pointer("Settling", "AIRDEN",
                              state->meteorology().AIRDEN ? state->meteorology().AIRDEN->host_data() : nullptr);
        double* delp = state->write_field<3>("DELP");
        double* z_edge = state->write_field<3>("Z");
        const double* pmid = state->read_field<3>("PMID");
        require_field_pointer("Settling", "DELP", delp);
        require_field_pointer("Settling", "PMID", pmid);
        require_field_pointer("Settling", "RH",
                              state->meteorology().RH ? state->meteorology().RH->host_data() : nullptr);
        require_field_pointer("Settling", "Z", z_edge);
        require_field_pointer("Settling", "CHEM_CONC",
                              state->chemistry().conc ? state->chemistry().conc->host_data() : nullptr);

        int bridge_rc = 0;
        run_settling_science_bridge(state->column_count(), state->level_count(), num_aerosols, state->species_count(),
                                    state->clock().timestep, gocart_scale_factor, gocart_swelling_rh_max,
                                    gocart_correction_maring ? 1 : 0, gocart_maring_dust_only ? 1 : 0,
                                    state->meteorology().AIRDEN->host_data(), delp, pmid,
                                    state->meteorology().RH->host_data(), state->meteorology().T->host_data(), z_edge,
                                    aerosol_species_names.data(), state->chemistry().species_names_c_arr.data(),
                                    host_is_dust.data(), host_is_hydrophilic.data(), host_radius_dry.data(),
                                    host_rhop_dry.data(), state->chemistry().conc->host_write(), &bridge_rc);
        if (bridge_rc != 0)
            throw std::runtime_error("Settling science bridge failed with status " + std::to_string(bridge_rc));
        state->chemistry().conc->mark_host_modified();
    }

    void SettlingProcess::finalize() {
        // Kokkos views will be deallocated automatically when their reference count goes to zero.
    }

} // namespace catchem

extern "C" {
void catchem_register_settling_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "settling", []() { return std::make_shared<catchem::SettlingProcess>(); }, {},
        catchem::make_settings_validator("settling",
                                         {"gocart/scale_factor", "gocart/simple_scheme", "gocart/swelling_rh_max",
                                          "gocart/correction_maring", "gocart/maring_dust_only"}));
}
}
