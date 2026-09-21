#include "catchem_process_drydep.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_error.hpp"
#include "catchem_logger.hpp"
#include "catchem_process_registry.hpp"
#include <iostream>

extern "C" {
void run_drydep_science_bridge(int n_cols, int n_levels, int n_species, double dt, const char* gas_scheme,
                               const char* aero_scheme, int diagnostics, double wesely_scale_factor,
                               int wesely_co2_effect, double wesely_co2_level, double wesely_co2_reference,
                               double gocart_scale_factor, int gocart_resuspension, int gocart_dust_resusp_only,
                               double zhang_scale_factor, double* bxheight, double* airden, double* t_air,
                               double* z_edges, double* rh, double* cldfrc, double* frlai, double* frlanduse,
                               int* iland, bool* is_ice, bool* is_land, bool* is_snow, double* lat, double* lon,
                               double* obk, double* ps, double* salinity, double* suncosmid, double* swgdn, double* ts,
                               double* tskin, double* ustar, double* z0, double* frlake, double* gwettop, double* hflux,
                               int* lwi, double* pblh, double* u10m, double* v10m, double* z0h, double* mw_g,
                               double* dd_f0, double* dd_hstar, double* dd_DvzAerSnow, double* dd_DvzMinVal_snow,
                               double* dd_DvzMinVal_land, double* density, double* radius, bool* is_seasalt,
                               bool* is_dust, double* lower_radius, double* upper_radius, bool* is_gas, double* conc,
                               double* tendency, const char* species_names, double* diag_con, double* diag_vel,
                               const int* diagnostic_species_id, int n_diag_species);
}

namespace catchem {

    ProcessContract DryDepProcess::get_contract() const {
        return make_contract(get_name(), {host_field_3d("T", "K"),
                                          host_field_3d("QV", "kg/kg"),
                                          host_field_3d("PMID", "Pa"),
                                          host_field_interface("PEDGE", "Pa"),
                                          host_field_interface("Z", "m"),
                                          host_field_3d("BXHEIGHT", "m"),
                                          host_field_3d("AIRDEN", "kg/m3"),
                                          host_field_3d("RH", "1"),
                                          host_field_3d("CLDF", "1"),
                                          host_field_2d("PS", "Pa"),
                                          host_field_2d("TS", "K"),
                                          host_field_2d("LAT", "degrees", FieldRequirement::Required,
                                                        AccessIntent::Read, PersistencePolicy::Persistent),
                                          host_field_2d("LON", "degrees", FieldRequirement::Required,
                                                        AccessIntent::Read, PersistencePolicy::Persistent),
                                          host_field_2d("USTAR", "m/s"),
                                          host_field_2d("HFLUX", "W/m2"),
                                          host_field_2d("OBK", "m"),
                                          host_field_2d("PBLH", "m"),
                                          host_field_2d("DLUSE", "1"),
                                          host_field_2d("LAI", "m2/m2"),
                                          host_field_2d("FRSNO", "frac"),
                                          host_field_2d("SWGDN", "W/m2"),
                                          host_field_2d("Z0", "m"),
                                          host_field_2d("SALINITY", "1", FieldRequirement::Optional),
                                          host_field_2d("FRLAKE", "1"),
                                          host_field_2d("GWETTOP", "1"),
                                          host_field_2d("LWI", "1"),
                                          host_field_2d("U10M", "m/s"),
                                          host_field_2d("V10M", "m/s"),
                                          host_concentration()});
    }

    DryDepProcess::DryDepProcess() : gas_scheme("wesely"), aero_scheme("gocart"), diagnostics_enabled(true) {}

    void DryDepProcess::prepare_inputs(std::shared_ptr<StateManager> state) {
        state->derive_bxheight();
        state->derive_airden();
        state->derive_relative_humidity();
        state->derive_salinity();
        state->derive_landuse_categories();
        state->derive_obk();
        state->derive_column_cloud_fraction();
        state->derive_suncosmid();
    }

    void DryDepProcess::init(std::shared_ptr<StateManager> state) {
        const auto config = state->config_manager();
        if (!config)
            throw std::invalid_argument("DryDep requires a runtime YAML configuration");
        const auto configured = config->data.processes.find("drydep");
        if (configured == config->data.processes.end())
            throw std::invalid_argument("DryDep requires a processes.drydep block in the runtime YAML");
        gas_scheme = configured->second.get_string("gas_scheme");
        aero_scheme = configured->second.get_string("aero_scheme");
        diagnostics_enabled = configured->second.diagnostics;
        if (gas_scheme != "wesely" || (aero_scheme != "gocart" && aero_scheme != "zhang"))
            throw std::invalid_argument("DryDep runtime YAML selected an unsupported gas/aerosol scheme combination");

        // Read per-scheme tuning options from the runtime YAML.  Defaults
        // mirror DryDepCommon_Mod.F90; the registered validator rejects any
        // option name the schemes do not declare.
        const auto& settings = configured->second;
        wesely_scale_factor = settings.get_double("wesely/scale_factor", wesely_scale_factor);
        wesely_co2_effect = settings.get_bool("wesely/co2_effect", wesely_co2_effect);
        wesely_co2_level = settings.get_double("wesely/co2_level", wesely_co2_level);
        wesely_co2_reference = settings.get_double("wesely/co2_reference", wesely_co2_reference);
        gocart_scale_factor = settings.get_double("gocart/scale_factor", gocart_scale_factor);
        gocart_resuspension = settings.get_bool("gocart/resuspension", gocart_resuspension);
        gocart_dust_resuspension_only =
            settings.get_bool("gocart/dust_resuspension_only", gocart_dust_resuspension_only);
        zhang_scale_factor = settings.get_double("zhang/scale_factor", zhang_scale_factor);

        // Surface the effective scheme options so the run log confirms what
        // was parsed from the runtime YAML and will be passed to the bridge.
        Logger::debug(state.get(), "DryDep scheme options",
                      {{"gas_scheme", gas_scheme},
                       {"aero_scheme", aero_scheme},
                       {"wesely/scale_factor", std::to_string(wesely_scale_factor)},
                       {"wesely/co2_effect", wesely_co2_effect ? "true" : "false"},
                       {"wesely/co2_level", std::to_string(wesely_co2_level)},
                       {"wesely/co2_reference", std::to_string(wesely_co2_reference)},
                       {"gocart/scale_factor", std::to_string(gocart_scale_factor)},
                       {"gocart/resuspension", gocart_resuspension ? "true" : "false"},
                       {"gocart/dust_resuspension_only", gocart_dust_resuspension_only ? "true" : "false"},
                       {"zhang/scale_factor", std::to_string(zhang_scale_factor)}});

        // 1. Setup diagnostic species ID dynamically based on the is_drydep metadata switch
        for (size_t i = 0; i < state->chemistry().species_list.size(); ++i) {
            if (state->chemistry().species_list[i].is_drydep) {
                diagnostic_species_id.push_back(i + 1); // 1-based for Fortran bridge
            }
        }

        if (!diagnostics_enabled)
            return;
        // 2. Register C++ Diagnostic fields.  The packed species axis is the
        // diagnostic (is_drydep) subset, NOT the full catalog: the schemes
        // scatter each species into its diag_idx slot in [1, n_diag_species],
        // so a [ncol, n_species] buffer leaves the non-drydep tail unwritten.
        // Registering at [ncol, n_diag] with a Species axis and per-slot labels
        // lets the NUOPC driver unpack one named variable per drydep species
        // (feature 013).  diagnostic_species_id holds GLOBAL 1-based catalog
        // positions, so slot i's label is species_list[id-1].short_name (FR-006).
        const int n_diag = static_cast<int>(diagnostic_species_id.size());
        // A mechanism may configure no is_drydep species; with an empty packed
        // axis there is nothing to register (matches the seasalt/settling guards).
        if (n_diag > 0) {
            std::vector<int> dims_2d = {state->column_count(), n_diag};
            std::vector<std::string> drydep_labels;
            drydep_labels.reserve(diagnostic_species_id.size());
            for (const int gid : diagnostic_species_id)
                drydep_labels.push_back(state->chemistry().species_list[static_cast<size_t>(gid) - 1].short_name);
            const std::vector<SemanticAxis> axes_2d = {SemanticAxis::Column, SemanticAxis::Species};
            state->diagnostic_manager()->register_field_contract(
                "drydep_con_per_species", "Deposition Concentration", "ug/kg", DiagType::FIELD_2D, dims_2d,
                DiagnosticPolicy::Instantaneous, 0.0, axes_2d, drydep_labels);
            state->diagnostic_manager()->register_field_contract(
                "drydep_velocity_per_species", "Deposition Velocity", "m/s", DiagType::FIELD_2D, dims_2d,
                DiagnosticPolicy::Instantaneous, 0.0, axes_2d, drydep_labels);
        }
    }

    void DryDepProcess::run(std::shared_ptr<StateManager> state) {

        // All derived fields are established in prepare_inputs.  Reads must
        // not mark host-owned NUOPC imports modified.
        const double* bxheight_ptr = state->read_field<3>("BXHEIGHT");
        const double* airden_ptr = state->read_field<3>("AIRDEN");
        const double* t_ptr = state->read_field<3>("T");
        const double* pedge_ptr = state->read_field<3>("PEDGE");
        const double* z_ptr = state->read_field<3>("Z");
        const double* rh_ptr = state->read_field<3>("RH");

        // 2. Retrieve surface met and grid positions
        const double* ps_ptr = state->read_field<2>("PS");
        const double* ts_ptr = state->read_field<2>("TS");
        const double* lat_ptr = state->read_field<2>("LAT");
        const double* lon_ptr = state->read_field<2>("LON");
        const double* ustar_ptr = state->read_field<2>("USTAR");
        const double* hflux_ptr = state->read_field<2>("HFLUX");
        const double* obk_ptr = state->read_field<2>("OBK");
        const double* pblh_ptr = state->read_field<2>("PBLH");

        require_field_pointer("DryDep", "BXHEIGHT", bxheight_ptr);
        require_field_pointer("DryDep", "AIRDEN", airden_ptr);
        require_field_pointer("DryDep", "T", t_ptr);
        require_field_pointer("DryDep", "PEDGE", pedge_ptr);
        // The GOCART aero scheme consumes geometric height in its z_edges slot;
        // feeding PEDGE (Pa) there produced NaN deposition velocities.
        require_field_pointer("DryDep", "Z", z_ptr);
        require_field_pointer("DryDep", "PS", ps_ptr);
        require_field_pointer("DryDep", "TS", ts_ptr);
        require_field_pointer("DryDep", "RH", rh_ptr);
        require_field_pointer("DryDep", "LAT", lat_ptr);
        require_field_pointer("DryDep", "LON", lon_ptr);
        require_field_pointer("DryDep", "USTAR", ustar_ptr);
        require_field_pointer("DryDep", "HFLUX", hflux_ptr);
        require_field_pointer("DryDep", "OBK", obk_ptr);
        require_field_pointer("DryDep", "PBLH", pblh_ptr);

        const double* cldf_ptr = state->read_field<3>("CLDF");
        const double* dluse_ptr = state->read_field<2>("DLUSE");
        const double* lai_ptr = state->read_field<2>("LAI");
        const double* frsno_ptr = state->read_field<2>("FRSNO");
        const double* swgdn_ptr = state->read_field<2>("SWGDN");
        const double* z0_ptr = state->read_field<2>("Z0");
        const double* salinity_ptr = state->read_field<2>("SALINITY");
        const double* frlake_ptr = state->read_field<2>("FRLAKE");
        const double* gwettop_ptr = state->read_field<2>("GWETTOP");
        const double* lwi_ptr = state->read_field<2>("LWI");
        const double* u10m_ptr = state->read_field<2>("U10M");
        const double* v10m_ptr = state->read_field<2>("V10M");
        require_field_pointer("DryDep", "CLDF", cldf_ptr);
        require_field_pointer("DryDep", "DLUSE", dluse_ptr);
        require_field_pointer("DryDep", "LAI", lai_ptr);
        require_field_pointer("DryDep", "FRSNO", frsno_ptr);
        require_field_pointer("DryDep", "SWGDN", swgdn_ptr);
        require_field_pointer("DryDep", "Z0", z0_ptr);
        require_field_pointer("DryDep", "SALINITY", salinity_ptr);
        require_field_pointer("DryDep", "FRLAKE", frlake_ptr);
        require_field_pointer("DryDep", "GWETTOP", gwettop_ptr);
        require_field_pointer("DryDep", "LWI", lwi_ptr);
        require_field_pointer("DryDep", "U10M", u10m_ptr);
        require_field_pointer("DryDep", "V10M", v10m_ptr);
        const double* cldfrc = state->read_field<2>("CLDFRC");
        const double* suncosmid = state->read_field<2>("SUNCOSMID");
        require_field_pointer("DryDep", "CLDFRC", cldfrc);
        require_field_pointer("DryDep", "SUNCOSMID", suncosmid);

        // Multi-dimensional standard C++20 views using standard layout_left to match Fortran column-major
        // These categorical fields are imported explicitly.  Reconstructing
        // them from a scalar LAI/DLUSE changes the Wesely calculation and is
        // not equivalent to the legacy core's MetState inputs.
        const double* frlai_ptr = state->read_field<3>("FRLAI");
        const double* frlanduse_ptr = state->read_field<3>("FRLANDUSE");
        const double* iland_ptr = state->read_field<3>("ILAND");
        require_field_pointer("DryDep", "FRLAI", frlai_ptr);
        require_field_pointer("DryDep", "FRLANDUSE", frlanduse_ptr);
        require_field_pointer("DryDep", "ILAND", iland_ptr);
        std::vector<double> frlai_storage(state->column_count() * 20);
        std::vector<double> frlanduse_storage(state->column_count() * 20);
        std::vector<int> iland_storage(state->column_count() * 20, 0);

        Kokkos::mdspan<double, Kokkos::extents<int, Kokkos::dynamic_extent, 1, 20>, Kokkos::layout_left> frlai(
            frlai_storage.data(), state->column_count());
        Kokkos::mdspan<double, Kokkos::extents<int, Kokkos::dynamic_extent, 1, 20>, Kokkos::layout_left> frlanduse(
            frlanduse_storage.data(), state->column_count());
        Kokkos::mdspan<int, Kokkos::extents<int, Kokkos::dynamic_extent, 1, 20>, Kokkos::layout_left> iland(
            iland_storage.data(), state->column_count());

        std::vector<char> is_ice(state->column_count()), is_land(state->column_count()), is_snow(state->column_count());
        std::vector<int> lwi(state->column_count());
        for (int c = 0; c < state->column_count(); ++c) {
            for (int category = 0; category < 20; ++category) {
                const std::size_t slot =
                    static_cast<std::size_t>(c) + static_cast<std::size_t>(state->column_count()) * category;
                frlai_storage[slot] = frlai_ptr[slot];
                frlanduse_storage[slot] = frlanduse_ptr[slot];
                iland_storage[slot] = static_cast<int>(iland_ptr[slot]);
            }
            lwi[c] = static_cast<int>(lwi_ptr[c]);
            is_land[c] = lwi[c] == 1;
            is_ice[c] = lwi[c] == 2;
            is_snow[c] = frsno_ptr[c] >= 0.5;
        }

        // 3. Extract chemical arrays & C++ allocated diagnostics
        double* conc_ptr = state->chemistry().conc ? state->chemistry().conc->host_write() : nullptr;
        require_field_pointer("DryDep", "CHEM_CONC", conc_ptr);

        // Allocate local tendencies buffer
        std::vector<double> mock_tendency(state->column_count() * state->level_count() * state->species_count(), 0.0);

        // Diagnostic fields are registered only when enabled.  The science
        // bridge accepts null diagnostic storage in concentration-only runs.
        double* diag_con = nullptr;
        double* diag_vel = nullptr;
        if (diagnostics_enabled) {
            diag_con = (double*)state->diagnostic_manager()->get_host_pointer("drydep_con_per_species");
            diag_vel = (double*)state->diagnostic_manager()->get_host_pointer("drydep_velocity_per_species");
        }

        // 4. Retrieve species configuration properties from ChemState
        std::vector<double> mw_g(state->species_count(), 0.0);
        std::vector<double> dd_f0(state->species_count(), 0.0);
        std::vector<double> dd_hstar(state->species_count(), 0.0);
        std::vector<double> dd_DvzAerSnow(state->species_count(), 0.0);
        std::vector<double> dd_DvzMinVal_snow(state->species_count(), 0.0);
        std::vector<double> dd_DvzMinVal_land(state->species_count(), 0.0);
        std::vector<double> density(state->species_count(), 0.0);
        std::vector<double> radius(state->species_count(), 0.0);
        std::vector<char> is_seasalt(state->species_count(), 0);
        std::vector<char> is_dust(state->species_count(), 0);
        std::vector<double> lower_radius(state->species_count(), 0.0);
        std::vector<double> upper_radius(state->species_count(), 0.0);
        std::vector<char> is_gas(state->species_count(), 1);

        for (size_t i = 0; i < state->chemistry().species_list.size(); ++i) {
            auto& meta = state->chemistry().species_list[i];
            if (!(meta.mw_g > 0.0))
                throw std::runtime_error("DryDep species '" + meta.short_name +
                                         "' requires an explicit molecular weight");
            mw_g[i] = meta.mw_g;
            dd_f0[i] = meta.dd_f0;
            dd_hstar[i] = meta.dd_hstar;
            dd_DvzAerSnow[i] = meta.dd_DvzAerSnow;
            dd_DvzMinVal_snow[i] = meta.dd_DvzMinVal_snow;
            dd_DvzMinVal_land[i] = meta.dd_DvzMinVal_land;
            if (meta.is_aerosol && (!(meta.density > 0.0) || !(meta.radius > 0.0) || !(meta.lower_radius > 0.0) ||
                                    !(meta.upper_radius > meta.lower_radius)))
                throw std::runtime_error("DryDep aerosol '" + meta.short_name +
                                         "' requires explicit density and radius bounds");
            density[i] = meta.density;
            radius[i] = meta.radius;
            is_seasalt[i] = meta.is_seasalt ? 1 : 0;
            is_dust[i] = meta.is_dust ? 1 : 0;
            lower_radius[i] = meta.lower_radius;
            upper_radius[i] = meta.upper_radius;
            is_gas[i] = meta.is_gas ? 1 : 0;
        }

        // 5. Invoke flat science bridge (casting char* vectors to bool* pointers)
        run_drydep_science_bridge(
            state->column_count(), state->level_count(), state->species_count(), state->clock().timestep,
            gas_scheme.c_str(), aero_scheme.c_str(), diagnostics_enabled ? 1 : 0, wesely_scale_factor,
            wesely_co2_effect ? 1 : 0, wesely_co2_level, wesely_co2_reference, gocart_scale_factor,
            gocart_resuspension ? 1 : 0, gocart_dust_resuspension_only ? 1 : 0, zhang_scale_factor,
            const_cast<double*>(bxheight_ptr), const_cast<double*>(airden_ptr), const_cast<double*>(t_ptr),
            const_cast<double*>(z_ptr), const_cast<double*>(rh_ptr), const_cast<double*>(cldfrc), frlai.data_handle(),
            frlanduse.data_handle(), iland.data_handle(), (bool*)is_ice.data(), (bool*)is_land.data(),
            (bool*)is_snow.data(), const_cast<double*>(lat_ptr), const_cast<double*>(lon_ptr),
            const_cast<double*>(obk_ptr), const_cast<double*>(ps_ptr), const_cast<double*>(salinity_ptr),
            const_cast<double*>(suncosmid), const_cast<double*>(swgdn_ptr), const_cast<double*>(ts_ptr),
            const_cast<double*>(ts_ptr), const_cast<double*>(ustar_ptr), const_cast<double*>(z0_ptr),
            const_cast<double*>(frlake_ptr), const_cast<double*>(gwettop_ptr), const_cast<double*>(hflux_ptr),
            lwi.data(), const_cast<double*>(pblh_ptr), const_cast<double*>(u10m_ptr), const_cast<double*>(v10m_ptr),
            const_cast<double*>(z0_ptr), mw_g.data(), dd_f0.data(), dd_hstar.data(), dd_DvzAerSnow.data(),
            dd_DvzMinVal_snow.data(), dd_DvzMinVal_land.data(), density.data(), radius.data(), (bool*)is_seasalt.data(),
            (bool*)is_dust.data(), lower_radius.data(), upper_radius.data(), (bool*)is_gas.data(), conc_ptr,
            mock_tendency.data(), state->chemistry().species_names_c_arr.data(), diag_con, diag_vel,
            diagnostic_species_id.data(), diagnostic_species_id.size());

        if (state->chemistry().conc)
            state->chemistry().conc->mark_host_modified();
    }

} // namespace catchem

extern "C" {
void catchem_register_drydep_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "drydep", []() { return std::make_shared<catchem::DryDepProcess>(); }, {},
        catchem::make_settings_validator("drydep",
                                         {"wesely/scale_factor", "wesely/co2_effect", "wesely/co2_level",
                                          "wesely/co2_reference", "gocart/scale_factor", "gocart/resuspension",
                                          "gocart/dust_resuspension_only", "zhang/scale_factor"}));
}
}
