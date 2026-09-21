#include "catchem_process_so4chem.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_error.hpp"
#include "catchem_logger.hpp"
#include "catchem_process_registry.hpp"
#include <array>
#include <iostream>

extern "C" {
void run_so4chem_science_bridge(int n_cols, int n_levels, int n_species, double dt, int diagnostics,
                                int gocart_update_so2, int year, int month, int day, int hour, int minute, int second,
                                double* airden, double* cldf, double* delp, double* pmid, double* t_air,
                                double* z_edges, double* hflux, double* lat, double* lon, int* lwi, double* pblh,
                                double* u10m, double* ustar, double* v10m, double* z0h, double* species_mw_g,
                                const char* species_names, double* conc, double* tendency, bool* c_firsttime,
                                int* c_nymd_last, int* c_nhms_last_recycle, double* c_xh2o2_init,
                                double* c_pso4_g_so2, double* c_pso4_aq_so2, double* c_pso2_dms, double* c_dms_flux,
                                double* c_diag_prod_rate, const int* diagnostic_species_id, int n_diag_species);
}

namespace catchem {

    ProcessContract SO4chemProcess::get_contract() const {
        return make_contract(get_name(), {host_field_3d("T", "K"), host_field_3d("PMID", "Pa"),
                                          host_field_interface("PEDGE", "Pa"), host_field_interface("Z", "m"),
                                          host_field_3d("DELP", "Pa"), host_field_3d("AIRDEN", "kg/m3"),
                                          host_field_3d("CLDF", "1"), host_field_2d("HFLUX", "W/m2"),
                                          host_field_2d("LAT", "degrees", FieldRequirement::Required,
                                                        AccessIntent::Read, PersistencePolicy::Persistent),
                                          host_field_2d("LON", "degrees", FieldRequirement::Required,
                                                        AccessIntent::Read, PersistencePolicy::Persistent),
                                          host_field_2d("PBLH", "m"), host_field_2d("USTAR", "m/s"),
                                          host_field_2d("U10M", "m/s"), host_field_2d("V10M", "m/s"),
                                          host_field_2d("LWI", "1"), host_field_2d("Z0", "m"), host_concentration()});
    }

    SO4chemProcess::SO4chemProcess() : active_scheme("gocart"), diagnostics_enabled(true) {}

    void SO4chemProcess::prepare_inputs(std::shared_ptr<StateManager> state) {
        state->derive_delp();
        state->derive_airden();
    }

    void SO4chemProcess::init(std::shared_ptr<StateManager> state) {
        const auto config = state->config_manager();
        if (!config)
            throw std::invalid_argument("SO4Chem requires a runtime YAML configuration");
        const auto configured = config->data.processes.find("so4chem");
        if (configured == config->data.processes.end() || configured->second.scheme != "gocart")
            throw std::invalid_argument("SO4Chem requires processes.so4chem.scheme: gocart");
        diagnostics_enabled = configured->second.diagnostics;

        // Read scheme tuning options from the runtime YAML.  Each lookup falls
        // back to the compiled default declared in SO4chemCommon_Mod.F90, so a
        // configuration that omits the option keeps current behavior.
        const auto& settings = configured->second;
        gocart_update_so2 = settings.get_bool("gocart/update_so2", gocart_update_so2);

        // Surface the effective scheme options so the run log confirms what
        // was parsed from the runtime YAML and will be passed to the bridge.
        Logger::debug(state.get(), "SO4Chem scheme options",
                      {{"scheme", active_scheme}, {"gocart/update_so2", gocart_update_so2 ? "true" : "false"}});

        // Preserve the unit contract of ProcessSO4chemInterface_Mod and
        // SO4chemScheme_GOCART_Mod: gases are carried in ppmv, while SO4 and
        // MSA are aerosol mass in ug/kg.  The science scheme has fixed
        // conversions for these four species, so accepting a different phase
        // classification would silently corrupt source strengths and lifetimes.
        struct SpeciesUnitContract {
            const char* name;
            bool is_gas;
        };
        constexpr std::array<SpeciesUnitContract, 4> unit_contract = {
            {{"dms", true}, {"so2", true}, {"so4", false}, {"msa", false}}};
        const auto& chemistry = state->chemistry();
        if (!chemistry.mechanism)
            throw std::invalid_argument("SO4Chem requires a loaded species mechanism");
        for (const auto& expected : unit_contract) {
            if (!chemistry.mechanism->contains(expected.name))
                throw std::invalid_argument(std::string("SO4Chem requires species '") + expected.name + "'");
            const auto index = chemistry.mechanism->index_of(expected.name);
            const auto& metadata = chemistry.species_list[index];
            if (metadata.is_gas != expected.is_gas || metadata.is_aerosol == expected.is_gas) {
                throw std::invalid_argument(std::string("SO4Chem requires '") + expected.name + "' to be a " +
                                            (expected.is_gas ? "gas (ppmv)" : "aerosol (ug/kg)"));
            }
        }

        // 1. Allocate persistent states
        firsttime.assign(state->column_count(), 1);
        nymd_last.assign(state->column_count(), 0);
        nhms_last_recycle.assign(state->column_count(), 0);
        xh2o2_init.assign(state->column_count() * state->level_count(), 0.0);
        pso4_g_so2.assign(state->column_count() * state->level_count(), 0.0);
        pso4_aq_so2.assign(state->column_count() * state->level_count(), 0.0);
        pso2_dms.assign(state->column_count() * state->level_count(), 0.0);
        dms_flux.assign(state->column_count(), 0.0);

        // 2. Resolve the diagnostic species set (parity with legacy).  The
        // GOCART scheme matches diagnostic_species_id(diag_idx) == species_idx
        // where species_idx is the GLOBAL catalog position, so ids live in
        // global 1-based space.  Default set = the species the legacy so4chem
        // process is configured with (the sulfur chain plus its oxidants); the
        // scheme only fills MSA/SO2/SO4 slots, so the others register as zeros
        // exactly as legacy does.  An explicit diag_species overrides the
        // default and must name a real mechanism species (fail-loud).  Names
        // resolve case-insensitively against the mechanism.
        // Built before the diagnostics gate so the ids are always available to
        // run() (mirrors dust/seasalt/carbchem).
        diagnostic_species_id.clear();
        {
            const std::vector<std::string> default_sulfur = {"dms", "so2", "so4", "msa", "h2o2", "oh", "no3",
                                                             "dms_in"};
            const auto& requested = settings.diag_species.empty() ? default_sulfur : settings.diag_species;
            for (const auto& species_name : requested) {
                if (!chemistry.mechanism->contains(species_name)) {
                    if (settings.diag_species.empty())
                        continue; // default set: tolerate species absent from this mechanism
                    throw std::invalid_argument("SO4Chem diag_species names an unknown species: " + species_name);
                }
                diagnostic_species_id.push_back(static_cast<int>(chemistry.mechanism->index_of(species_name)) + 1);
            }
        }

        if (!diagnostics_enabled)
            return;

        // 3. Register diagnostics.  Axes are explicit so the writer never
        // shape-guesses: a single-level column would otherwise mis-resolve a
        // _per_level field as a singleton (feature 013).
        if (state->diagnostic_manager()) {
            std::vector<int> dims_2d = {state->column_count(), state->level_count()};
            std::vector<int> dims_1d = {state->column_count(), 1};
            const std::vector<SemanticAxis> axes_level = {SemanticAxis::Column, SemanticAxis::Level};
            const std::vector<SemanticAxis> axes_single = {SemanticAxis::Column, SemanticAxis::Singleton};

            state->diagnostic_manager()->register_field_contract("PSO4_from_gaseous_SO2_per_level", "PSO4 gas source",
                                                                 "kg/kg/s", DiagType::FIELD_2D, dims_2d,
                                                                 DiagnosticPolicy::Instantaneous, 0.0, axes_level);
            state->diagnostic_manager()->register_field_contract("PSO4_from_aqueous_SO2_per_level", "PSO4 aq source",
                                                                 "kg/kg/s", DiagType::FIELD_2D, dims_2d,
                                                                 DiagnosticPolicy::Instantaneous, 0.0, axes_level);
            state->diagnostic_manager()->register_field_contract("DMS_emission_flux", "DMS emission surface flux",
                                                                 "kg/m2/s", DiagType::FIELD_2D, dims_1d,
                                                                 DiagnosticPolicy::Instantaneous, 0.0, axes_single);

            // One Production_rate_<sp> field per selected species (names
            // unchanged from the legacy per-species convention).  run() fills
            // each field from its OWN slot of the bridge's packed buffer.
            for (const int global_index : diagnostic_species_id) {
                const auto& meta = state->chemistry().species_list[static_cast<std::size_t>(global_index - 1)];
                std::string diag_name = "Production_rate_" + meta.short_name;
                state->diagnostic_manager()->register_field_contract(
                    diag_name, "Production rate " + meta.short_name, "kg/kg/s", DiagType::FIELD_2D, dims_2d,
                    DiagnosticPolicy::Instantaneous, 0.0, axes_level);
            }
        }
    }

    void SO4chemProcess::run(std::shared_ptr<StateManager> state) {

        // 1. Retrieve 3D Meteorological variables
        double* airden_ptr = state->write_field<3>("AIRDEN");

        double* pmid_ptr = state->write_field<3>("PMID");
        double* t_ptr = state->write_field<3>("T");
        double* z_ptr = state->write_field<3>("Z");
        double* cldf_ptr = state->write_field<3>("CLDF");

        double* delp_ptr = state->write_field<3>("DELP");

        // 2. Retrieve 2D Surface Met variables
        double* hflux_ptr = state->write_field<2>("HFLUX");

        double* lat_ptr = state->write_field<2>("LAT");

        double* lon_ptr = state->write_field<2>("LON");

        double* pblh_ptr = state->write_field<2>("PBLH");

        double* ustar_ptr = state->write_field<2>("USTAR");

        double* u10m_ptr = state->write_field<2>("U10M");

        double* v10m_ptr = state->write_field<2>("V10M");

        double* lwi_ptr = state->write_field<2>("LWI");
        std::vector<int> lwi(state->column_count());
        require_field_pointer("SO4chem", "LWI", lwi_ptr);
        for (int col = 0; col < state->column_count(); ++col)
            lwi[col] = static_cast<int>(lwi_ptr[col]);

        require_field_pointer("SO4chem", "AIRDEN", airden_ptr);
        require_field_pointer("SO4chem", "PMID", pmid_ptr);
        require_field_pointer("SO4chem", "T", t_ptr);
        require_field_pointer("SO4chem", "Z", z_ptr);
        require_field_pointer("SO4chem", "DELP", delp_ptr);
        require_field_pointer("SO4chem", "CLDF", cldf_ptr);
        require_field_pointer("SO4chem", "HFLUX", hflux_ptr);
        require_field_pointer("SO4chem", "LAT", lat_ptr);
        require_field_pointer("SO4chem", "LON", lon_ptr);
        require_field_pointer("SO4chem", "PBLH", pblh_ptr);
        require_field_pointer("SO4chem", "USTAR", ustar_ptr);
        require_field_pointer("SO4chem", "U10M", u10m_ptr);
        require_field_pointer("SO4chem", "V10M", v10m_ptr);

        const double* z0_ptr = state->read_field<2>("Z0");
        require_field_pointer("SO4chem", "Z0", z0_ptr);

        // 3. Chemical and Tendency Views
        double* conc_ptr = state->chemistry().conc ? state->chemistry().conc->host_write() : nullptr;
        require_field_pointer("SO4chem", "CHEM_CONC", conc_ptr);

        // Allocate local tendencies buffer
        std::vector<double> mock_tendency(state->column_count() * state->level_count() * state->species_count(), 0.0);

        // 4. Retrieve species properties from ChemState
        std::vector<double> mw_g(state->species_count(), 0.0);
        for (size_t i = 0; i < state->chemistry().species_list.size(); ++i) {
            if (!(state->chemistry().species_list[i].mw_g > 0.0))
                throw std::runtime_error("SO4Chem species '" + state->chemistry().species_list[i].short_name +
                                         "' requires an explicit molecular weight");
            mw_g[i] = state->chemistry().species_list[i].mw_g;
        }

        // 5. Invoke flat science bridge.  diagnostic_species_id (built in
        // init) holds GLOBAL 1-based catalog indices — the space the GOCART
        // scheme matches on.  The packed diag_prod_rate scratch receives one
        // [ncol, nz] slab per selected species so each Production_rate_<sp>
        // field gets its OWN slot (the former pso4_so2 buffer collapsed every
        // species onto slot 1).  When diagnostics are off, forward a valid
        // size-1 id dummy holding 0 and a zero count (mirrors the other
        // bridges' no_diag_species guard).
        const bool forward_diag = diagnostics_enabled && !diagnostic_species_id.empty();
        static const int no_diag_species = 0;
        const int* diag_ids = forward_diag ? diagnostic_species_id.data() : &no_diag_species;
        const int n_diag_species = forward_diag ? static_cast<int>(diagnostic_species_id.size()) : 0;
        const int slab = state->column_count() * state->level_count();
        std::vector<double> diag_prod_rate(static_cast<size_t>(slab) * (forward_diag ? n_diag_species : 1), 0.0);
        run_so4chem_science_bridge(
            state->column_count(), state->level_count(), state->species_count(), state->clock().timestep,
            diagnostics_enabled ? 1 : 0, gocart_update_so2 ? 1 : 0, state->clock().year, state->clock().month,
            state->clock().day, state->clock().hour, state->clock().minute, state->clock().second, airden_ptr, cldf_ptr,
            delp_ptr, pmid_ptr, t_ptr, z_ptr, hflux_ptr, lat_ptr, lon_ptr, lwi.data(), pblh_ptr, u10m_ptr, ustar_ptr,
            v10m_ptr, const_cast<double*>(z0_ptr), mw_g.data(), state->chemistry().species_names_c_arr.data(), conc_ptr,
            mock_tendency.data(), (bool*)firsttime.data(), nymd_last.data(), nhms_last_recycle.data(),
            xh2o2_init.data(), pso4_g_so2.data(), pso4_aq_so2.data(), pso2_dms.data(), dms_flux.data(),
            diag_prod_rate.data(), diag_ids, n_diag_species);

        // 6. Map persistent column diagnostics straight to registered C++ Diagnostics Views
        if (state->diagnostic_manager() && diagnostics_enabled) {
            double* diag_pso4_g =
                (double*)state->diagnostic_manager()->get_host_pointer("PSO4_from_gaseous_SO2_per_level");
            double* diag_pso4_aq =
                (double*)state->diagnostic_manager()->get_host_pointer("PSO4_from_aqueous_SO2_per_level");
            double* diag_dms_flux = (double*)state->diagnostic_manager()->get_host_pointer("DMS_emission_flux");

            if (diag_pso4_g)
                std::copy(pso4_g_so2.begin(), pso4_g_so2.end(), diag_pso4_g);
            if (diag_pso4_aq)
                std::copy(pso4_aq_so2.begin(), pso4_aq_so2.end(), diag_pso4_aq);
            if (diag_dms_flux)
                std::copy(dms_flux.begin(), dms_flux.end(), diag_dms_flux);

            // Scatter each packed slot into its own Production_rate_<sp> field
            // (field names unchanged from the legacy per-species convention).
            for (int d = 0; d < n_diag_species; ++d) {
                const auto& meta = state->chemistry().species_list[static_cast<std::size_t>(diagnostic_species_id[d]) - 1];
                std::string diag_name = "Production_rate_" + meta.short_name;
                double* diag_prod = (double*)state->diagnostic_manager()->get_host_pointer(diag_name);
                if (diag_prod)
                    std::copy(diag_prod_rate.begin() + static_cast<std::ptrdiff_t>(d) * slab,
                              diag_prod_rate.begin() + static_cast<std::ptrdiff_t>(d + 1) * slab, diag_prod);
            }
        }

        if (state->chemistry().conc)
            state->chemistry().conc->mark_host_modified();
    }

} // namespace catchem

extern "C" {
void catchem_register_so4chem_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "so4chem", []() { return std::make_shared<catchem::SO4chemProcess>(); }, {},
        catchem::make_settings_validator("so4chem", {"gocart/update_so2"}));
}
}
