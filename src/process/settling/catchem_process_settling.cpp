#include "catchem_process_settling.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_error.hpp"
#include "catchem_logger.hpp"
#include "catchem_process_registry.hpp"
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace catchem {

    extern "C" void run_settling_science_bridge(
        int n_columns, int n_levels, int n_aerosols, int n_total_species, double dt, double scale_factor,
        double swelling_rh_max, int correction_maring, int maring_dust_only, double* airden, double* delp,
        const double* pmid, double* rh, double* temperature, double* z_edge, const char* aerosol_species_names,
        const char* species_names, const int* species_is_dust, const int* species_is_hydrophilic, const double* radius,
        const double* density, double* concentration, int simple_scheme, const char* aerosol_mie_names,
        double* diag_velocity, double* diag_flux, const int* diagnostic_species_id, int n_diag_species, int* bridge_rc);

    // Loads the aerosol optics (Mie) tables named by the top-level "mie:" section
    // into the Fortran bridge store.  type_names/file_paths are fixed-width,
    // blank-padded character buffers packed contiguously (32 and 512 bytes per
    // entry).  init_rc receives the 1-based index of the first file that failed
    // to open, or 0 on success.  Idempotent: a second call replaces any
    // previously loaded set.
    extern "C" void run_settling_mie_init(int n_files, const char* type_names, const char* file_paths, int* init_rc);

    namespace {
        // Right-trim the blank padding of a fixed-width Fortran character buffer.
        std::string trim_trailing_spaces(std::string text) {
            text.erase(text.find_last_not_of(' ') + 1);
            return text;
        }
    } // namespace

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

        diagnostics_enabled = configured->second.diagnostics;

        // Surface the effective scheme options so the run log confirms what
        // was parsed from the runtime YAML and reaches the settling kernel.
        // The metadata path corresponds to legacy simple_scheme: false.
        // simple_scheme: true restores the legacy optics-table path: the
        // NetCDF Mie tables named by the top-level "mie:" section are loaded
        // once here and the scheme selects Chem_SettlingSimple per species via
        // each species' __mie_name (specs/012).
        std::vector<std::string> mie_type_labels; // "<TYPE> -> <file>" for the log
        if (gocart_simple_scheme) {
            const auto& mie = config->data.mie;
            if (mie.files.empty())
                throw std::invalid_argument("Settling gocart/simple_scheme requires a non-empty 'mie.files' section "
                                            "mapping aerosol types to optics tables");
            const std::filesystem::path dir(mie.directory);
            if (mie.directory != "./" && !std::filesystem::is_directory(dir))
                throw std::invalid_argument("Settling optics directory is not a readable directory: '" + mie.directory +
                                            "'");
            std::vector<char> type_names, file_paths;
            type_names.reserve(mie.files.size() * 32);
            file_paths.reserve(mie.files.size() * 512);
            for (const auto& [type, file] : mie.files) {
                const std::filesystem::path joined = mie.directory == "./" ? std::filesystem::path(file) : dir / file;
                std::error_code ec;
                if (!std::filesystem::exists(joined, ec))
                    throw std::invalid_argument("Settling optics table for aerosol type '" + type + "' not found at '" +
                                                joined.string() + "'");
                type_names.insert(type_names.end(), 32, ' ');
                std::copy_n(type.begin(), std::min<size_t>(type.size(), 32), type_names.end() - 32);
                file_paths.insert(file_paths.end(), 512, ' ');
                const std::string path_str = joined.string();
                std::copy_n(path_str.begin(), std::min<size_t>(path_str.size(), 511), file_paths.end() - 512);
                mie_type_labels.push_back(type + " -> " + path_str);
            }
            int init_rc = 0;
            run_settling_mie_init(static_cast<int>(mie.files.size()), type_names.data(), file_paths.data(), &init_rc);
            if (init_rc != 0) {
                const auto& failed = mie.files[static_cast<size_t>(init_rc) - 1];
                throw std::runtime_error("Settling optics table load failed for file #" + std::to_string(init_rc) +
                                         " ('" + failed.first + "': " + failed.second + ")");
            }
            mie_initialized = true;
        }
        Logger::debug(
            state.get(), "Settling scheme options",
            {{"scheme", active_scheme},
             {"gocart/scale_factor", std::to_string(gocart_scale_factor)},
             {"gocart/correction_maring", gocart_correction_maring ? "true" : "false"},
             {"gocart/maring_dust_only", gocart_maring_dust_only ? "true" : "false"},
             {"gocart/simple_scheme", gocart_simple_scheme ? "true (optics tables)" : "false (species metadata)"},
             {"gocart/swelling",
              gocart_simple_scheme ? "optics-table rEff(rh) LUT" : "per-species __hydrophilic (Gerber when true)"},
             {"gocart/swelling_rh_max", std::to_string(gocart_swelling_rh_max)},
             {"gocart/optics_directory", gocart_simple_scheme ? config->data.mie.directory : "(unused)"},
             {"gocart/optics_tables", gocart_simple_scheme ? std::to_string(mie_type_labels.size()) : "0"}});
        for (const auto& label : mie_type_labels)
            Logger::debug(state.get(), "Settling optics table", {{"binding", label}});

        int num_aerosols = state->chemistry().aerosol_indices.size();
        if (num_aerosols > 0) {
            aerosol_species_names.assign(static_cast<size_t>(num_aerosols) * 32, ' ');
            aerosol_mie_names.assign(static_cast<size_t>(num_aerosols) * 32, ' ');
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
                // swelling_method knob.  Ignored on the optics-table path where
                // the table itself encodes the size response.
                host_is_hydrophilic[i] = state->chemistry().species_list[ispec].is_hydrophilic ? 1 : 0;
                std::copy_n(state->chemistry().species_names_c_arr.data() + static_cast<size_t>(ispec) * 32, 32,
                            aerosol_species_names.data() + static_cast<size_t>(i) * 32);
                const std::string& mie_name = state->chemistry().species_list[ispec].mie_name;
                std::copy_n(mie_name.begin(), std::min<size_t>(mie_name.size(), 32),
                            aerosol_mie_names.data() + static_cast<size_t>(i) * 32);
            }

            if (gocart_simple_scheme) {
                // Defensive C++-side pre-resolution so misconfiguration names the
                // species at init instead of surfacing as the Fortran bridge_rc==2
                // backstop at the first step (specs/012 FR-009).
                for (int i = 0; i < num_aerosols; ++i) {
                    const int ispec = state->chemistry().aerosol_indices[i];
                    const std::string species_name = state->chemistry().species_list[ispec].short_name;
                    const std::string trimmed =
                        trim_trailing_spaces(std::string(aerosol_mie_names.data() + static_cast<size_t>(i) * 32, 32));
                    if (trimmed.empty())
                        throw std::invalid_argument("Settling simple_scheme requires __mie_name on species '" +
                                                    species_name + "' but it is empty");
                    bool resolved = false;
                    for (const auto& [type, file] : config->data.mie.files)
                        if (type == trimmed) {
                            resolved = true;
                            break;
                        }
                    if (!resolved)
                        throw std::invalid_argument("Settling species '" + species_name + "' has __mie_name '" +
                                                    trimmed + "' which matches no configured mie.files entry");
                }
            }
        } else if (gocart_simple_scheme) {
            throw std::invalid_argument("Settling simple_scheme is enabled but no aerosol species are configured");
        }

        // --- Per-process diagnostics (parity with legacy) ------------------
        // diagnostic_species_id indexes the aerosol subset (1..num_aerosols),
        // which is the species_idx space compute_gocart iterates.
        {
            const auto& aerosol_idx = state->chemistry().aerosol_indices;
            const auto& settings = configured->second;
            std::vector<int> selected_local; // 1-based positions into aerosol subset
            if (!settings.diag_species.empty()) {
                for (const auto& name : settings.diag_species) {
                    int found = -1;
                    for (size_t a = 0; a < aerosol_idx.size(); ++a) {
                        if (state->chemistry().species_list[aerosol_idx[a]].short_name == name) {
                            found = static_cast<int>(a) + 1;
                            break;
                        }
                    }
                    if (found < 0)
                        throw std::invalid_argument("Settling diag_species names a non-settling species: " + name);
                    selected_local.push_back(found);
                }
            } else {
                for (size_t a = 0; a < aerosol_idx.size(); ++a)
                    selected_local.push_back(static_cast<int>(a) + 1);
            }
            diagnostic_species_id = selected_local;
            diagnostic_species_names.clear();
            for (int local : selected_local)
                diagnostic_species_names.push_back(
                    state->chemistry().species_list[aerosol_idx[static_cast<size_t>(local) - 1]].short_name);

            if (diagnostics_enabled && state->diagnostic_manager() && !selected_local.empty()) {
                const int ndiag = static_cast<int>(selected_local.size());
                std::vector<int> dims_vel = {state->column_count(), state->level_count(), ndiag};
                std::vector<int> dims_flux = {state->column_count(), ndiag};
                // diagnostic_species_names is the packed-axis label list, built in
                // the same aerosol-subset order the scheme iterates (FR-006); the
                // NUOPC driver unpacks each field into one named variable per
                // species (feature 013).
                const std::vector<SemanticAxis> axes_vel = {SemanticAxis::Column, SemanticAxis::Level,
                                                            SemanticAxis::Species};
                const std::vector<SemanticAxis> axes_flux = {SemanticAxis::Column, SemanticAxis::Species};
                state->diagnostic_manager()->register_field_contract(
                    "settling_velocity_per_species_per_level", "Settling velocity", "m/s", DiagType::FIELD_3D, dims_vel,
                    DiagnosticPolicy::Instantaneous, 0.0, axes_vel, diagnostic_species_names);
                state->diagnostic_manager()->register_field_contract("settling_flux_per_species", "Settling column flux",
                                                                    "kg/m2/s", DiagType::FIELD_2D, dims_flux,
                                                                    DiagnosticPolicy::Instantaneous, 0.0, axes_flux,
                                                                    diagnostic_species_names);
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
        double* diag_velocity = nullptr;
        double* diag_flux = nullptr;
        if (diagnostics_enabled && state->diagnostic_manager() && !diagnostic_species_id.empty()) {
            diag_velocity =
                static_cast<double*>(state->diagnostic_manager()->get_host_pointer(
                    "settling_velocity_per_species_per_level"));
            diag_flux =
                static_cast<double*>(state->diagnostic_manager()->get_host_pointer("settling_flux_per_species"));
        }
        const int n_diag_species = diagnostics_enabled ? static_cast<int>(diagnostic_species_id.size()) : 0;
        // The Fortran dummy argument is declared diagnostic_species_id(max(n_diag_species,1)),
        // so the bridge always receives a size-1 intent(in) array even when the count is zero
        // (it is never dereferenced in that case).  std::vector::data() of an empty vector may
        // be nullptr, which would form a Fortran pointer to nothing; pass a valid dummy instead.
        static const int no_diag_species = 0;
        const int* diag_ids =
            diagnostic_species_id.empty() ? &no_diag_species : diagnostic_species_id.data();

        run_settling_science_bridge(
            state->column_count(), state->level_count(), num_aerosols, state->species_count(), state->clock().timestep,
            gocart_scale_factor, gocart_swelling_rh_max, gocart_correction_maring ? 1 : 0,
            gocart_maring_dust_only ? 1 : 0, state->meteorology().AIRDEN->host_data(), delp, pmid,
            state->meteorology().RH->host_data(), state->meteorology().T->host_data(), z_edge,
            aerosol_species_names.data(), state->chemistry().species_names_c_arr.data(), host_is_dust.data(),
            host_is_hydrophilic.data(), host_radius_dry.data(), host_rhop_dry.data(),
            state->chemistry().conc->host_write(), gocart_simple_scheme ? 1 : 0, aerosol_mie_names.data(),
            diag_velocity, diag_flux, diag_ids, n_diag_species, &bridge_rc);
        if (bridge_rc == 2)
            throw std::runtime_error(
                "Settling optics-table mapping failed: a settling species did not resolve to a loaded Mie table");
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
