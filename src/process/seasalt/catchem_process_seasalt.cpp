#include "catchem_process_seasalt.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_error.hpp"
#include "catchem_logger.hpp"
#include "catchem_process_registry.hpp"
#include <array>
#include <cctype>
#include <iostream>

extern "C" {
void run_seasalt_science_bridge(int n_cols, int n_levels, int n_species, int n_total_species, double dt,
                                const char* active_scheme, int diagnostics, double gong97_scale_factor,
                                int gong97_weibull_flag, double gong03_scale_factor, int gong03_weibull_flag,
                                double geos12_scale_factor, int geos12_weibull_flag, double* frocean, double* frseaice,
                                double* lat, double* lon, double* sst, double* u10m, double* v10m, double* ustar,
                                double* delp, double* density, double* radius, double* lower_radius,
                                double* upper_radius, bool* is_gas, double* mw_g, const char* bin_species_names,
                                const char* species_names, double* conc, double* tendency, double* diag_mass_total,
                                double* diag_num_total, double* diag_mass_bin, double* diag_num_bin,
                                const int* diagnostic_species_id, int n_diag_species);
}

namespace catchem {

    namespace {

        /// Resolve the canonical sea-salt bins (SEAS1..SEAS5) to global
        /// species indices.  Shared by init() and run() so the diagnostic bin
        /// order and the physics bin order are guaranteed identical: the
        /// schemes' local species_idx space (1..n_seasalt) is defined by THIS
        /// sequence.  Every bin must exist by canonical name and carry valid
        /// physical properties (same throws run() has always enforced).
        std::vector<int> resolved_seasalt_bins(const std::shared_ptr<StateManager>& state) {
            // Route legacy size bins by canonical name rather than declaration
            // position.  Their ordering is part of the GEOS/Gong scheme contract.
            constexpr std::array<const char*, 5> seasalt_bin_names = {"SEAS1", "SEAS2", "SEAS3", "SEAS4", "SEAS5"};
            std::vector<int> ss_global_indices;
            ss_global_indices.reserve(seasalt_bin_names.size());
            for (const char* name : seasalt_bin_names) {
                const auto found = state->chemistry().species_name_to_index.find(name);
                if (found == state->chemistry().species_name_to_index.end())
                    throw std::runtime_error("SeaSalt requires canonical species '" + std::string(name) + "'");
                const int index = found->second;
                const auto& meta = state->chemistry().species_list[index];
                if (!meta.is_seasalt)
                    throw std::runtime_error("SeaSalt species '" + meta.short_name + "' must set is_seasalt: true");
                if (!(meta.density > 0.0 && meta.radius > 0.0 && meta.lower_radius > 0.0 &&
                      meta.upper_radius > meta.lower_radius && meta.mw_g > 0.0))
                    throw std::runtime_error("SeaSalt species '" + meta.short_name +
                                             "' requires explicit density, radius bounds, and molecular weight");
                ss_global_indices.push_back(index);
            }
            return ss_global_indices;
        }

    } // namespace

    ProcessContract SeaSaltProcess::get_contract() const {
        std::vector<FieldAccessContract> fields{host_field_interface("PEDGE", "Pa"),
                                                host_field_3d("DELP", "Pa"),
                                                host_field_2d("FROCEAN", "frac"),
                                                host_field_2d("FRSEAICE", "frac"),
                                                host_field_2d("TS", "K"),
                                                host_field_2d("LAT", "degrees", FieldRequirement::Required,
                                                              AccessIntent::Read, PersistencePolicy::Persistent),
                                                host_field_2d("LON", "degrees", FieldRequirement::Required,
                                                              AccessIntent::Read, PersistencePolicy::Persistent),
                                                host_field_2d("U10M", "m/s"),
                                                host_field_2d("V10M", "m/s"),
                                                host_concentration()};
        if (active_scheme == "geos12")
            fields.push_back(host_field_2d("USTAR", "m/s"));
        return make_contract(get_name(), std::move(fields));
    }

    SeaSaltProcess::SeaSaltProcess() : active_scheme("geos12"), diagnostics_enabled(true) {}

    void SeaSaltProcess::prepare_inputs(std::shared_ptr<StateManager> state) {
        state->derive_delp();
    }

    void SeaSaltProcess::init(std::shared_ptr<StateManager> state) {
        const auto config = state->config_manager();
        if (!config)
            throw std::invalid_argument("SeaSalt requires a runtime YAML configuration");
        const auto configured = config->data.processes.find("seasalt");
        if (configured == config->data.processes.end() || configured->second.scheme.empty())
            throw std::invalid_argument("SeaSalt requires processes.seasalt.scheme in the runtime YAML");
        active_scheme = configured->second.scheme;
        diagnostics_enabled = configured->second.diagnostics;
        if (active_scheme != "gong97" && active_scheme != "gong03" && active_scheme != "geos12")
            throw std::invalid_argument("SeaSalt runtime YAML selected unsupported scheme: " + active_scheme);

        // Read per-scheme tuning options before the diagnostics early-return:
        // they must reach the bridge even when diagnostics are off.
        const auto& settings = configured->second;
        gong97_scale_factor = settings.get_double("gong97/scale_factor", gong97_scale_factor);
        gong97_weibull_flag = settings.get_bool("gong97/weibull_flag", gong97_weibull_flag);
        gong03_scale_factor = settings.get_double("gong03/scale_factor", gong03_scale_factor);
        gong03_weibull_flag = settings.get_bool("gong03/weibull_flag", gong03_weibull_flag);
        geos12_scale_factor = settings.get_double("geos12/scale_factor", geos12_scale_factor);
        geos12_weibull_flag = settings.get_bool("geos12/weibull_flag", geos12_weibull_flag);

        // Surface the effective scheme options so the run log confirms what
        // was parsed from the runtime YAML and will be passed to the bridge.
        Logger::debug(state.get(), "SeaSalt scheme options",
                      {{"scheme", active_scheme},
                       {"gong97/scale_factor", std::to_string(gong97_scale_factor)},
                       {"gong97/weibull_flag", gong97_weibull_flag ? "true" : "false"},
                       {"gong03/scale_factor", std::to_string(gong03_scale_factor)},
                       {"gong03/weibull_flag", gong03_weibull_flag ? "true" : "false"},
                       {"geos12/scale_factor", std::to_string(geos12_scale_factor)},
                       {"geos12/weibull_flag", geos12_weibull_flag ? "true" : "false"}});

        // 3. Per-process diagnostics (parity with legacy): the schemes match
        // diagnostic_species_id(diag_idx) == species_idx, where species_idx is
        // the LOCAL bin position (1..n_seasalt) in the canonical SEAS1..SEAS5
        // order run() feeds the bridge.  Build the ids in THAT space, honoring
        // a diag_species subset (names matched against short_name exactly as
        // dust/settling do).
        const auto ss_global_indices = resolved_seasalt_bins(state);
        {
            std::vector<int> selected_local; // 1-based positions in the bin list
            if (!settings.diag_species.empty()) {
                for (const auto& name : settings.diag_species) {
                    int found = -1;
                    for (size_t bin = 0; bin < ss_global_indices.size(); ++bin) {
                        if (state->chemistry().species_list[ss_global_indices[bin]].short_name == name) {
                            found = static_cast<int>(bin) + 1;
                            break;
                        }
                    }
                    if (found < 0)
                        throw std::invalid_argument("SeaSalt diag_species names a non-sea-salt species: " + name);
                    selected_local.push_back(found);
                }
            } else {
                for (size_t bin = 0; bin < ss_global_indices.size(); ++bin)
                    selected_local.push_back(static_cast<int>(bin) + 1);
            }
            diagnostic_species_id = selected_local;
        }

        if (!diagnostics_enabled)
            return;
        if (state->diagnostic_manager()) {
            std::vector<int> dims_1d = {state->column_count(), 1};
            state->diagnostic_manager()->register_field("seasalt_mass_emission_total", "Total Mass Emission", "kg/m2/s",
                                                        DiagType::FIELD_2D, dims_1d);
            state->diagnostic_manager()->register_field("seasalt_number_emission_total", "Total Number Emission",
                                                        "#/m2/s", DiagType::FIELD_2D, dims_1d);

            // Per-bin emissions register as one compact [ncols, n_diag] field
            // each with a Category axis and per-slot labels; the NUOPC driver
            // unpacks them into one named 2D variable per size bin (feature 013).
            // n_diag is the (possibly subsetted) diagnostic bin count;
            // column-major layout matches the bridge's [n_cols, n_diag_species].
            const int n_diag = static_cast<int>(diagnostic_species_id.size());
            if (n_diag > 0) {
                std::vector<int> dims_bins = {state->column_count(), n_diag};
                // diagnostic_species_id holds LOCAL bin positions (1-based) into
                // the canonical resolved_seasalt_bins order run() feeds the
                // bridge, so slot i's label is that bin's short_name (FR-006).
                std::vector<std::string> ss_bin_labels;
                ss_bin_labels.reserve(ss_global_indices.size());
                for (const int local : diagnostic_species_id)
                    ss_bin_labels.push_back(
                        state->chemistry().species_list[ss_global_indices[static_cast<size_t>(local) - 1]].short_name);
                const std::vector<SemanticAxis> axes_bin = {SemanticAxis::Column, SemanticAxis::Category};
                state->diagnostic_manager()->register_field_contract("seasalt_mass_emission_bins", "Mass Emission Per Bin",
                                                                    "kg/m2/s", DiagType::FIELD_2D, dims_bins,
                                                                    DiagnosticPolicy::Instantaneous, 0.0, axes_bin,
                                                                    ss_bin_labels);
                state->diagnostic_manager()->register_field_contract("seasalt_number_emission_bins", "Number Emission Per Bin",
                                                                    "#/m2/s", DiagType::FIELD_2D, dims_bins,
                                                                    DiagnosticPolicy::Instantaneous, 0.0, axes_bin,
                                                                    ss_bin_labels);
            }
        }
    }

    void SeaSaltProcess::run(std::shared_ptr<StateManager> state) {

        double* frocean_ptr = state->write_field<2>("FROCEAN");
        double* frseaice_ptr = state->write_field<2>("FRSEAICE");
        double* sst_ptr = state->write_field<2>("TS");
        double* lat_ptr = state->write_field<2>("LAT");
        double* lon_ptr = state->write_field<2>("LON");
        double* ustar_ptr = state->write_field<2>("USTAR");
        double* u10m_ptr = state->write_field<2>("U10M");
        double* v10m_ptr = state->write_field<2>("V10M");

        double* delp_ptr = state->write_field<3>("DELP");

        require_field_pointer("SeaSalt", "FROCEAN", frocean_ptr);
        require_field_pointer("SeaSalt", "FRSEAICE", frseaice_ptr);
        require_field_pointer("SeaSalt", "LAT", lat_ptr);
        require_field_pointer("SeaSalt", "LON", lon_ptr);
        require_field_pointer("SeaSalt", "SST", sst_ptr);
        if (active_scheme == "geos12")
            require_field_pointer("SeaSalt", "USTAR", ustar_ptr);
        require_field_pointer("SeaSalt", "U10M", u10m_ptr);
        require_field_pointer("SeaSalt", "V10M", v10m_ptr);
        require_field_pointer("SeaSalt", "DELP", delp_ptr);

        // 2. Identify and slice SeaSalt-only chemical species to comply with science solver limits
        std::vector<double> density;
        std::vector<double> radius;
        std::vector<double> lower_radius;
        std::vector<double> upper_radius;
        std::vector<char> is_gas;
        std::vector<double> mw_g;

        // Canonical SEAS1..SEAS5 resolution is shared with init() so the
        // diagnostic ids (built there) index the same local bin order the
        // physics arrays use here.
        const auto ss_global_indices = resolved_seasalt_bins(state);
        for (const int index : ss_global_indices) {
            const auto& meta = state->chemistry().species_list[index];
            density.push_back(meta.density);
            radius.push_back(meta.radius);
            lower_radius.push_back(meta.lower_radius);
            upper_radius.push_back(meta.upper_radius);
            is_gas.push_back(meta.is_gas ? 1 : 0);
            mw_g.push_back(meta.mw_g);
        }

        int n_seasalt = ss_global_indices.size();
        if (n_seasalt == 0)
            return; // No sea salt species configured

        // Extract chemical concentration raw host pointer
        double* conc_ptr = state->chemistry().conc ? state->chemistry().conc->host_write() : nullptr;
        require_field_pointer("SeaSalt", "CHEM_CONC", conc_ptr);

        // The bridge now operates directly on the full unified concentration
        // array and resolves each sea-salt bin to its slot by name.  No slice
        // buffers and no manual copy-back (mirrors settling/drydep/wetdep).
        const int n_total_species = state->species_count();

        // Pack the sea-salt bin names into the 32-char-per-name flat layout the
        // bridge expects, in the same order as the metadata arrays.
        // Upper-case to match the catalog convention in species_names_c_arr
        // (ChemState stores short names upper-cased); the bridge compares by
        // trimmed name, so both sides must use the same case.
        std::vector<char> bin_species_names(static_cast<size_t>(32) * n_seasalt, ' ');
        for (int i = 0; i < n_seasalt; ++i) {
            const std::string& nm = state->chemistry().species_list[ss_global_indices[i]].short_name;
            for (size_t c = 0; c < nm.size() && c < 32; ++c)
                bin_species_names[static_cast<size_t>(i) * 32 + c] =
                    static_cast<char>(std::toupper(static_cast<unsigned char>(nm[c])));
        }

        // Full-width tendency scratch (the bridge writes only the sea-salt slots).
        std::vector<double> full_tendency(
            static_cast<size_t>(state->column_count()) * state->level_count() * n_total_species, 0.0);

        // 3. Extract diagnostics
        double* diag_mass_total_ptr =
            diagnostics_enabled && state->diagnostic_manager()
                ? (double*)state->diagnostic_manager()->get_host_pointer("seasalt_mass_emission_total")
                : nullptr;
        double* diag_num_total_ptr =
            diagnostics_enabled && state->diagnostic_manager()
                ? (double*)state->diagnostic_manager()->get_host_pointer("seasalt_number_emission_total")
                : nullptr;

        // Per-bin diagnostics write straight into the compact registered
        // fields ([ncols, n_diag], column-major).  The bridge's Fortran shape
        // is [n_cols, n_diag_species], so the memory layout matches with no
        // gather/scatter step.
        double* diag_mass_bin_ptr =
            diagnostics_enabled && state->diagnostic_manager()
                ? (double*)state->diagnostic_manager()->get_host_pointer("seasalt_mass_emission_bins")
                : nullptr;
        double* diag_num_bin_ptr =
            diagnostics_enabled && state->diagnostic_manager()
                ? (double*)state->diagnostic_manager()->get_host_pointer("seasalt_number_emission_bins")
                : nullptr;

        // diagnostic_species_id (built in init) holds 1-based LOCAL bin
        // positions, optionally subset by processes.seasalt.diag_species: the
        // scheme runs over all n_seasalt bins internally but only fills the
        // diagnostic slots named in the subset.  The Fortran dummy is
        // declared diagnostic_species_id(max(n_diag_species,1)); when
        // diagnostics are off the bridge still forwards the array to the
        // scheme, so pass a valid size-1 dummy holding 0 — it matches no
        // 1-based species_idx and the scheme writes nothing (mirrors
        // settling/dust's no_diag_species guard).
        const bool forward_diag = diagnostics_enabled && !diagnostic_species_id.empty();
        static const int no_diag_species = 0;
        const int* diag_ids = forward_diag ? diagnostic_species_id.data() : &no_diag_species;
        const int n_diag_species = forward_diag ? static_cast<int>(diagnostic_species_id.size()) : 0;
        if (forward_diag) {
            const int registered_bins =
                diag_mass_bin_ptr
                    ? static_cast<int>(
                          state->diagnostic_manager()->get_field("seasalt_mass_emission_bins")->dimensions[1])
                    : 0;
            if (registered_bins < n_diag_species)
                throw std::runtime_error("SeaSalt diagnostic bin capacity (" + std::to_string(registered_bins) +
                                         ") is smaller than the " + std::to_string(n_diag_species) +
                                         " selected diagnostic sea-salt bins");
        }

        // 5. Invoke flat science bridge (operates on the full conc array;
        // resolves sea-salt bins by name internally).
        run_seasalt_science_bridge(
            state->column_count(), state->level_count(), n_seasalt, n_total_species, state->clock().timestep,
            active_scheme.c_str(), diagnostics_enabled ? 1 : 0, gong97_scale_factor, gong97_weibull_flag ? 1 : 0,
            gong03_scale_factor, gong03_weibull_flag ? 1 : 0, geos12_scale_factor, geos12_weibull_flag ? 1 : 0,
            frocean_ptr, frseaice_ptr, lat_ptr, lon_ptr, sst_ptr, u10m_ptr, v10m_ptr, ustar_ptr, delp_ptr,
            density.data(), radius.data(), lower_radius.data(), upper_radius.data(), (bool*)is_gas.data(), mw_g.data(),
            bin_species_names.data(), state->chemistry().species_names_c_arr.data(), conc_ptr, full_tendency.data(),
            diag_mass_total_ptr, diag_num_total_ptr, diag_mass_bin_ptr, diag_num_bin_ptr, diag_ids, n_diag_species);

        if (state->chemistry().conc)
            state->chemistry().conc->mark_host_modified();
    }

} // namespace catchem

extern "C" {
void catchem_register_seasalt_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "seasalt", []() { return std::make_shared<catchem::SeaSaltProcess>(); }, {},
        catchem::make_settings_validator("seasalt",
                                         {"gong97/scale_factor", "gong97/weibull_flag", "gong03/scale_factor",
                                          "gong03/weibull_flag", "geos12/scale_factor", "geos12/weibull_flag"}));
}
}
