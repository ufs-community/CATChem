#include "catchem_process_dust.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_error.hpp"
#include "catchem_logger.hpp"
#include "catchem_process_registry.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

extern "C" {
void run_dust_science_bridge(int n_cols, int n_levels, int n_species, int n_total_species, int n_soil, double dt,
                             const char* active_scheme, int diagnostics, double fengsha_alpha, double fengsha_gamma,
                             double fengsha_drylimit_factor, double fengsha_moisture_factor, double fengsha_kvhmax,
                             int fengsha_drag_option, int fengsha_horizflux_option, int fengsha_moist_option,
                             int fengsha_distribution_option, const double* ginoux_ch_du, int n_ginoux_ch_du,
                             const double* airden, const double* delp, const double* clayfrac, const double* frlake,
                             const double* frsno, const double* gvf, const double* lai, int* lwi, const double* rdrag,
                             const double* sandfrac, const double* soilm, const double* gwettop, const double* ssm,
                             const double* tskin, const double* u10m, const double* v10m, const double* ustar,
                             const double* ustar_threshold, const double* z0, const double* species_density,
                             const double* species_radius, const double* species_lower_radius,
                             const double* species_upper_radius, const char* bin_species_names,
                             const char* species_names, double* conc, double* tendency, double* diag_emission_total,
                             double* diag_emission_bin, double* diag_horizontal_flux, double* diag_moisture_correction,
                             double* diag_effective_threshold, double* diag_utar_threshold,
                             const int* diagnostic_species_id, int n_diag_species);
}

namespace catchem {

    namespace {

        /// Resolve the canonical dust bins (DUST1..DUSTn) to global species
        /// indices.  Shared by init() and run() so the diagnostic bin order
        /// and the physics bin order are guaranteed identical: the schemes'
        /// local species_idx space (1..n_dust) is defined by THIS sequence.
        /// The bin count is derived from the is_dust metadata, never
        /// hardcoded; every bin must exist by canonical name and carry valid
        /// physical properties (same throws run() has always enforced).
        std::vector<int> resolved_dust_bins(const std::shared_ptr<StateManager>& state) {
            std::size_t bin_count = 0;
            for (const auto& meta : state->chemistry().species_list)
                if (meta.is_dust)
                    ++bin_count;
            std::vector<int> dust_global_indices;
            dust_global_indices.reserve(bin_count);
            for (std::size_t bin = 1; bin <= bin_count; ++bin) {
                const std::string name = "DUST" + std::to_string(bin);
                const auto found = state->chemistry().species_name_to_index.find(name);
                if (found == state->chemistry().species_name_to_index.end())
                    throw std::runtime_error("Dust requires canonical species '" + name + "'");
                const int index = found->second;
                const auto& meta = state->chemistry().species_list[index];
                if (!meta.is_dust)
                    throw std::runtime_error("Dust species '" + meta.short_name + "' must set is_dust: true");
                if (!(meta.density > 0.0 && meta.radius > 0.0 && meta.lower_radius > 0.0 &&
                      meta.upper_radius > meta.lower_radius))
                    throw std::runtime_error("Dust species '" + meta.short_name +
                                             "' requires explicit positive density, radius, lower_radius, and "
                                             "upper_radius");
                dust_global_indices.push_back(index);
            }
            return dust_global_indices;
        }

    } // namespace

    ProcessContract DustProcess::get_contract() const {
        std::vector<FieldAccessContract> fields{host_field_interface("PEDGE", "Pa"), host_field_3d("DELP", "Pa"),
                                                host_field_3d("AIRDEN", "kg/m3"), host_concentration()};
        if (active_scheme == "fengsha") {
            fields.insert(fields.end(),
                          {host_field_soil_layer("SOILM", "m3/m3"), host_field_2d("CLAYFRAC", "1"),
                           host_field_2d("FRLAKE", "1"), host_field_2d("FRSNO", "frac"), host_field_2d("GVF", "frac"),
                           host_field_2d("LAI", "m2/m2", FieldRequirement::Optional), host_field_2d("LWI", "1"),
                           host_field_2d("RDRAG", "1"), host_field_2d("SNDFRC", "1"), host_field_2d("SSM", "1"),
                           host_field_2d("TS", "K"), host_field_2d("USTAR", "m/s"),
                           host_field_2d("USTAR_THRESHOLD", "m/s"), host_field_2d("Z0", "m")});
        } else if (active_scheme == "ginoux") {
            fields.insert(fields.end(),
                          {host_field_2d("FRLAKE", "1"), host_field_2d("FRSNO", "frac"), host_field_2d("GWETTOP", "1"),
                           host_field_2d("LWI", "1"), host_field_2d("SSM", "1"), host_field_2d("TS", "K"),
                           host_field_2d("U10M", "m/s"), host_field_2d("V10M", "m/s")});
        }
        return make_contract(get_name(), std::move(fields));
    }

    DustProcess::DustProcess() : active_scheme("fengsha"), diagnostics_enabled(true) {}

    void DustProcess::prepare_inputs(std::shared_ptr<StateManager> state) {
        state->derive_delp();
        // derive_airden preserves a current host field and rebuilds its
        // process-owned field after begin_import_generation() invalidates
        // timestep data.  Testing pointer existence here misses that latter
        // case on the second and later coupling steps.
        state->derive_airden();
    }

    void DustProcess::init(std::shared_ptr<StateManager> state) {
        const auto config = state->config_manager();
        if (!config)
            throw std::invalid_argument("Dust requires a runtime YAML configuration");
        const auto configured = config->data.processes.find("dust");
        if (configured == config->data.processes.end() || configured->second.scheme.empty())
            throw std::invalid_argument("Dust requires processes.dust.scheme in the runtime YAML");
        active_scheme = configured->second.scheme;
        diagnostics_enabled = configured->second.diagnostics;
        if (active_scheme != "fengsha" && active_scheme != "ginoux")
            throw std::invalid_argument("Dust runtime YAML selected unsupported scheme: " + active_scheme);

        // 2. Read scheme tuning options from the runtime YAML.  Each lookup
        // falls back to the compiled default declared in DustCommon_Mod.F90,
        // so a configuration that omits the option keeps current behavior.
        // Core validates the option names against the registered schema, so a
        // misspelled key fails at initialization rather than being dropped.
        const auto& settings = configured->second;
        fengsha_alpha = settings.get_double("fengsha/alpha", fengsha_alpha);
        fengsha_gamma = settings.get_double("fengsha/gamma", fengsha_gamma);
        fengsha_drylimit_factor = settings.get_double("fengsha/drylimit_factor", fengsha_drylimit_factor);
        fengsha_moist_correction_factor =
            settings.get_double("fengsha/moist_correction_factor", fengsha_moist_correction_factor);
        fengsha_kvhmax = settings.get_double("fengsha/kvhmax", fengsha_kvhmax);
        fengsha_drag_option = settings.get_int("fengsha/drag_option", fengsha_drag_option);
        fengsha_horizflux_option = settings.get_int("fengsha/horizflux_option", fengsha_horizflux_option);
        fengsha_moist_option = settings.get_int("fengsha/moist_option", fengsha_moist_option);
        fengsha_distribution_option = settings.get_int("fengsha/distribution_option", fengsha_distribution_option);

        // Ch_DU carries one multiplier per dust size bin; the scheme type
        // declares exactly five bins, so a provided sequence must match.
        auto ch_du = settings.get_vector("ginoux/Ch_DU");
        if (!ch_du.empty()) {
            if (ch_du.size() != ginoux_ch_du.size())
                throw std::invalid_argument("Dust ginoux Ch_DU must declare " + std::to_string(ginoux_ch_du.size()) +
                                            " values, one per dust size bin");
            ginoux_ch_du = std::move(ch_du);
        }

        // Surface the effective scheme options so the run log confirms what
        // was parsed from the runtime YAML and will be passed to the bridge.
        {
            std::string ch_du_joined;
            for (size_t i = 0; i < ginoux_ch_du.size(); ++i) {
                if (i)
                    ch_du_joined += ",";
                ch_du_joined += std::to_string(ginoux_ch_du[i]);
            }
            Logger::debug(state.get(), "Dust scheme options",
                          {{"scheme", active_scheme},
                           {"fengsha/alpha", std::to_string(fengsha_alpha)},
                           {"fengsha/gamma", std::to_string(fengsha_gamma)},
                           {"fengsha/drylimit_factor", std::to_string(fengsha_drylimit_factor)},
                           {"fengsha/moist_correction_factor", std::to_string(fengsha_moist_correction_factor)},
                           {"fengsha/kvhmax", std::to_string(fengsha_kvhmax)},
                           {"fengsha/drag_option", std::to_string(fengsha_drag_option)},
                           {"fengsha/horizflux_option", std::to_string(fengsha_horizflux_option)},
                           {"fengsha/moist_option", std::to_string(fengsha_moist_option)},
                           {"fengsha/distribution_option", std::to_string(fengsha_distribution_option)},
                           {"ginoux/Ch_DU", ch_du_joined}});
        }

        // 3. Per-process diagnostics (parity with legacy): the schemes match
        // diagnostic_species_id(diag_idx) == species_idx, where species_idx is
        // the LOCAL bin position (1..n_dust) in the canonical-name bin order
        // run() feeds the bridge.  Build the ids in THAT space, honoring a
        // diag_species subset (names matched against short_name exactly as
        // settling does).
        const auto dust_global_indices = resolved_dust_bins(state);
        {
            std::vector<int> selected_local; // 1-based positions in the bin list
            if (!settings.diag_species.empty()) {
                for (const auto& name : settings.diag_species) {
                    int found = -1;
                    for (size_t bin = 0; bin < dust_global_indices.size(); ++bin) {
                        if (state->chemistry().species_list[dust_global_indices[bin]].short_name == name) {
                            found = static_cast<int>(bin) + 1;
                            break;
                        }
                    }
                    if (found < 0)
                        throw std::invalid_argument("Dust diag_species names a non-dust species: " + name);
                    selected_local.push_back(found);
                }
            } else {
                for (size_t bin = 0; bin < dust_global_indices.size(); ++bin)
                    selected_local.push_back(static_cast<int>(bin) + 1);
            }
            diagnostic_species_id = selected_local;
        }

        if (!diagnostics_enabled)
            return;

        // 4. Register C++ Diagnostic fields (registering 1D fields as 2D with second dimension of 1).
        // Per-bin fields use a compact [ncols, n_diag] layout with a Category
        // axis and per-slot labels; the NUOPC driver unpacks them into one
        // named 2D variable per bin (feature 013).  The science bridge already
        // writes exactly that column-major shape.
        const int n_diag = static_cast<int>(diagnostic_species_id.size());
        std::vector<int> dims_1d_as_2d = {state->column_count(), 1};
        std::vector<int> dims_bins = {state->column_count(), n_diag};

        // Per-bin labels: diagnostic_species_id holds LOCAL bin positions (1-based)
        // into the canonical resolved_dust_bins order run() feeds the bridge, so
        // slot i's label is the short_name of that bin — never a global catalog
        // index (index-space rule, FR-006).
        std::vector<std::string> dust_bin_labels;
        dust_bin_labels.reserve(dust_global_indices.size());
        for (const int local : diagnostic_species_id)
            dust_bin_labels.push_back(
                state->chemistry().species_list[dust_global_indices[static_cast<size_t>(local) - 1]].short_name);
        const std::vector<SemanticAxis> axes_bin = {SemanticAxis::Column, SemanticAxis::Category};

        state->diagnostic_manager()->register_field("dust_emission_total", "Total Dust Emission", "kg/m2/s",
                                                    DiagType::FIELD_2D, dims_1d_as_2d);
        state->diagnostic_manager()->register_field_contract("dust_emission_bin", "Dust Emission Per Bin", "kg/m2/s",
                                                             DiagType::FIELD_2D, dims_bins,
                                                             DiagnosticPolicy::Instantaneous, 0.0, axes_bin,
                                                             dust_bin_labels);
        state->diagnostic_manager()->register_field("dust_horizontal_flux", "Dust Horizontal Flux", "kg/m/s",
                                                    DiagType::FIELD_2D, dims_1d_as_2d);
        state->diagnostic_manager()->register_field("dust_moisture_correction", "Dust Moisture Correction", "unitless",
                                                    DiagType::FIELD_2D, dims_1d_as_2d);
        state->diagnostic_manager()->register_field("dust_effective_threshold", "Dust Effective Threshold", "m/s",
                                                    DiagType::FIELD_2D, dims_1d_as_2d);
        state->diagnostic_manager()->register_field_contract("dust_utar_threshold", "Dust Ustar Threshold Per Bin",
                                                             "m/s", DiagType::FIELD_2D, dims_bins,
                                                             DiagnosticPolicy::Instantaneous, 0.0, axes_bin,
                                                             dust_bin_labels);
    }

    void DustProcess::run(std::shared_ptr<StateManager> state) {

        // The execution plan invokes prepare_inputs before run(), but direct
        // API users and focused process tests may call run() themselves.
        // These derivations are generation-aware and therefore preserve
        // host-provided fields while supplying only absent prerequisites
        // (matches SettlingProcess::run).
        prepare_inputs(state);

        // 1. Retrieve Meteorological state pointers
        const double* airden_ptr = state->read_field<3>("AIRDEN");
        const double* delp_ptr = state->read_field<3>("DELP");
        const double* clayfrac_ptr = state->read_field<2>("CLAYFRAC");
        const double* frlake_ptr = state->read_field<2>("FRLAKE");
        const double* frsno_ptr = state->read_field<2>("FRSNO");
        const double* gvf_ptr = state->read_field<2>("GVF");
        const double* lai_ptr = state->read_field<2>("LAI");
        const double* lwi_double_ptr = state->read_field<2>("LWI");
        std::vector<int> lwi(state->column_count());
        require_field_pointer("Dust", "LWI", lwi_double_ptr);
        for (int col = 0; col < state->column_count(); ++col)
            lwi[col] = static_cast<int>(lwi_double_ptr[col]);
        const double* rdrag_ptr = state->read_field<2>("RDRAG");
        const double* sandfrac_ptr = state->read_field<2>("SNDFRC");
        const double* soilm_ptr = state->read_field<3>("SOILM");
        const auto soilm_field = state->find_field<3>("SOILM");
        const int n_soil = soilm_field ? static_cast<int>(soilm_field->extent(1)) : 0;
        const double* gwettop_ptr = state->read_field<2>("GWETTOP");
        const double* ssm_ptr = state->read_field<2>("SSM");
        const double* tskin_ptr = state->read_field<2>("TS");
        const double* u10m_ptr = state->read_field<2>("U10M");
        const double* v10m_ptr = state->read_field<2>("V10M");
        const double* ustar_ptr = state->read_field<2>("USTAR");
        const double* ustar_th_ptr = state->read_field<2>("USTAR_THRESHOLD");
        const double* z0_ptr = state->read_field<2>("Z0");

        require_field_pointer("Dust", "AIRDEN", airden_ptr);
        require_field_pointer("Dust", "DELP", delp_ptr);
        if (active_scheme == "fengsha") {
            require_field_pointer("Dust", "CLAYFRAC", clayfrac_ptr);
            require_field_pointer("Dust", "FRLAKE", frlake_ptr);
            require_field_pointer("Dust", "FRSNO", frsno_ptr);
            require_field_pointer("Dust", "GVF", gvf_ptr);
            require_field_pointer("Dust", "LAI", lai_ptr);
            require_field_pointer("Dust", "RDRAG", rdrag_ptr);
            require_field_pointer("Dust", "SNDFRC", sandfrac_ptr);
            require_field_pointer("Dust", "SOILM", soilm_ptr);
            if (n_soil <= 0)
                throw std::runtime_error("Dust: SOILM has no soil-layer extent");
            require_field_pointer("Dust", "SSM", ssm_ptr);
            require_field_pointer("Dust", "TS", tskin_ptr);
            require_field_pointer("Dust", "USTAR", ustar_ptr);
            require_field_pointer("Dust", "USTAR_THRESHOLD", ustar_th_ptr);
            require_field_pointer("Dust", "Z0", z0_ptr);
        } else {
            require_field_pointer("Dust", "FRLAKE", frlake_ptr);
            require_field_pointer("Dust", "FRSNO", frsno_ptr);
            require_field_pointer("Dust", "GWETTOP", gwettop_ptr);
            require_field_pointer("Dust", "SSM", ssm_ptr);
            require_field_pointer("Dust", "TS", tskin_ptr);
            require_field_pointer("Dust", "U10M", u10m_ptr);
            require_field_pointer("Dust", "V10M", v10m_ptr);
        }

        const auto config = state->config_manager();
        if (config && config->data.simulation.verbose_enabled && Logger::enabled(Logger::Level::Debug)) {
            // DELP is level-major (index = column + level*n_cols).  Because
            // derive_delp() takes pressure_thickness(PEDGE[L], PEDGE[L+1]) with
            // no abs(), a non-zero DELP field only exists when pressure
            // descends with increasing level index, i.e. level 0 is the surface.
            // Log the surface (level 0) and top (last level) layer thicknesses
            // explicitly so a run can confirm both the vertical ordering and the
            // surface-layer depth without inferring them from the array max.
            const int n_levels = state->level_count();
            const int n_columns = state->column_count();
            double min_delp = std::numeric_limits<double>::infinity();
            double max_delp = 0.0;
            double surf_min = std::numeric_limits<double>::infinity();
            double surf_max = 0.0;
            double top_min = std::numeric_limits<double>::infinity();
            double top_max = 0.0;
            for (int level = 0; level < n_levels; ++level) {
                for (int column = 0; column < n_columns; ++column) {
                    const std::size_t index =
                        static_cast<std::size_t>(column) + static_cast<std::size_t>(level) * n_columns;
                    if (!std::isfinite(delp_ptr[index]) || delp_ptr[index] <= 0.0)
                        continue;
                    min_delp = std::min(min_delp, delp_ptr[index]);
                    max_delp = std::max(max_delp, delp_ptr[index]);
                    if (level == 0) {
                        surf_min = std::min(surf_min, delp_ptr[index]);
                        surf_max = std::max(surf_max, delp_ptr[index]);
                    } else if (level == n_levels - 1) {
                        top_min = std::min(top_min, delp_ptr[index]);
                        top_max = std::max(top_max, delp_ptr[index]);
                    }
                }
            }
            Logger::debug(state.get(), "Dust layer-mass conversion inputs",
                          {{"delp_pa_min", std::to_string(min_delp)},
                           {"delp_pa_max", std::to_string(max_delp)},
                           {"delp_surface_pa_min", std::to_string(surf_min)},
                           {"delp_surface_pa_max", std::to_string(surf_max)},
                           {"delp_top_pa_min", std::to_string(top_min)},
                           {"delp_top_pa_max", std::to_string(top_max)},
                           {"conversion", "flux*g/DELP*1e9 kg/kg-to-ug/kg"}});
        }

        // 2. Diagnostic Views
        double* diag_emission_total = nullptr;
        double* diag_emission_bin = nullptr;
        double* diag_horizontal_flux = nullptr;
        double* diag_moisture_correction = nullptr;
        double* diag_effective_threshold = nullptr;
        double* diag_utar_threshold = nullptr;

        if (state->diagnostic_manager() && diagnostics_enabled) {
            diag_emission_total = (double*)state->diagnostic_manager()->get_host_pointer("dust_emission_total");
            diag_emission_bin = (double*)state->diagnostic_manager()->get_host_pointer("dust_emission_bin");
            diag_horizontal_flux = (double*)state->diagnostic_manager()->get_host_pointer("dust_horizontal_flux");
            diag_moisture_correction =
                (double*)state->diagnostic_manager()->get_host_pointer("dust_moisture_correction");
            diag_effective_threshold =
                (double*)state->diagnostic_manager()->get_host_pointer("dust_effective_threshold");
            diag_utar_threshold = (double*)state->diagnostic_manager()->get_host_pointer("dust_utar_threshold");
        }

        double* conc_ptr = state->chemistry().conc ? state->chemistry().conc->host_write() : nullptr;
        require_field_pointer("Dust", "CHEM_CONC", conc_ptr);

        // 4. Route the legacy dust bins by their canonical names.  The science
        // kernel assigns bin physics positionally, so deriving that position
        // from YAML declaration order would silently map a reordered species
        // list onto the wrong bins.  init() resolves the SAME sequence for
        // the diagnostic bin order (resolved_dust_bins), keeping the local
        // species_idx space identical on both sides of the bridge.
        const auto dust_global_indices = resolved_dust_bins(state);
        std::vector<double> density;
        std::vector<double> radius;
        std::vector<double> lower_radius;
        std::vector<double> upper_radius;
        density.reserve(dust_global_indices.size());
        radius.reserve(dust_global_indices.size());
        lower_radius.reserve(dust_global_indices.size());
        upper_radius.reserve(dust_global_indices.size());
        for (const int index : dust_global_indices) {
            const auto& meta = state->chemistry().species_list[index];
            density.push_back(meta.density);
            radius.push_back(meta.radius);
            lower_radius.push_back(meta.lower_radius);
            upper_radius.push_back(meta.upper_radius);
        }

        const int n_dust = static_cast<int>(dust_global_indices.size());
        if (n_dust == 0) {
            return;
        }

        // The bridge now operates directly on the full unified concentration
        // array and resolves each dust bin to its slot by name.  No slice
        // buffers and no manual copy-back (mirrors settling/drydep/wetdep).
        const int n_total_species = state->species_count();

        // Pack the dust bin names into the 32-char-per-name flat layout the
        // bridge expects, in the same order as the metadata arrays.
        // Upper-case to match the catalog convention in species_names_c_arr
        // (ChemState stores short names upper-cased); the bridge compares by
        // trimmed name, so both sides must use the same case.
        std::vector<char> bin_species_names(static_cast<size_t>(32) * n_dust, ' ');
        for (int local_idx = 0; local_idx < n_dust; ++local_idx) {
            const std::string& nm = state->chemistry().species_list[dust_global_indices[local_idx]].short_name;
            for (size_t c = 0; c < nm.size() && c < 32; ++c)
                bin_species_names[static_cast<size_t>(local_idx) * 32 + c] =
                    static_cast<char>(std::toupper(static_cast<unsigned char>(nm[c])));
        }

        // Full-width tendency scratch (the bridge writes only the dust slots).
        std::vector<double> full_tendency(
            static_cast<size_t>(state->column_count()) * state->level_count() * n_total_species, 0.0);

        // diagnostic_species_id (built in init) holds 1-based LOCAL bin
        // positions, optionally subset by processes.dust.diag_species: the
        // scheme runs over all n_dust bins internally but only fills the
        // diagnostic slots named in the subset.  The Fortran dummy is
        // declared diagnostic_species_id(max(n_diag_species,1)); when
        // diagnostics are off the bridge still forwards the array to the
        // scheme, so pass a valid size-1 dummy holding 0 — it matches no
        // 1-based species_idx and the scheme writes nothing (mirrors
        // settling's no_diag_species guard).
        const bool forward_diag = diagnostics_enabled && !diagnostic_species_id.empty();
        static const int no_diag_species = 0;
        const int* diag_ids = forward_diag ? diagnostic_species_id.data() : &no_diag_species;
        const int n_diag_species = forward_diag ? static_cast<int>(diagnostic_species_id.size()) : 0;

        // 5. Invoke flat science bridge
        run_dust_science_bridge(
            state->column_count(), state->level_count(), n_dust, n_total_species, n_soil, state->clock().timestep,
            active_scheme.c_str(), diagnostics_enabled ? 1 : 0, fengsha_alpha, fengsha_gamma, fengsha_drylimit_factor,
            fengsha_moist_correction_factor, fengsha_kvhmax, fengsha_drag_option, fengsha_horizflux_option,
            fengsha_moist_option, fengsha_distribution_option, ginoux_ch_du.data(),
            static_cast<int>(ginoux_ch_du.size()), airden_ptr, delp_ptr, clayfrac_ptr, frlake_ptr, frsno_ptr, gvf_ptr,
            lai_ptr, lwi.data(), rdrag_ptr, sandfrac_ptr, soilm_ptr, gwettop_ptr, ssm_ptr, tskin_ptr, u10m_ptr,
            v10m_ptr, ustar_ptr, ustar_th_ptr, z0_ptr, density.data(), radius.data(), lower_radius.data(),
            upper_radius.data(), bin_species_names.data(), state->chemistry().species_names_c_arr.data(), conc_ptr,
            full_tendency.data(), diag_emission_total, diag_emission_bin, diag_horizontal_flux,
            diag_moisture_correction, diag_effective_threshold, diag_utar_threshold, diag_ids, n_diag_species);

        if (state->chemistry().conc)
            state->chemistry().conc->mark_host_modified();
    }

    void DustProcess::finalize() {}

} // namespace catchem

extern "C" {
void catchem_register_dust_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "dust", []() { return std::make_shared<catchem::DustProcess>(); }, {},
        catchem::make_settings_validator(
            "dust", {"fengsha/alpha", "fengsha/gamma", "fengsha/drylimit_factor", "fengsha/moist_correction_factor",
                     "fengsha/kvhmax", "fengsha/drag_option", "fengsha/horizflux_option", "fengsha/moist_option",
                     "fengsha/distribution_option", "ginoux/Ch_DU"}));
}
}
