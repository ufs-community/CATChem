// src/process/photolysis/catchem_process_photolysis.cpp
#include "catchem_process_photolysis.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_error.hpp"
#include "catchem_logger.hpp"
#include "catchem_process_registry.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <musica/tuvx/grid.hpp>
#include <musica/tuvx/grid_map.hpp>
#include <musica/tuvx/profile.hpp>
#include <musica/tuvx/profile_map.hpp>
#include <musica/tuvx/radiator_map.hpp>
#include <musica/tuvx/tuvx_c_interface.hpp>
#include <sstream>
#include <unordered_set>
#include <yaml-cpp/yaml.h>

namespace catchem {

    ProcessContract PhotolysisProcess::get_contract() const {
        return make_contract(get_name(),
                             {host_field_2d("LAT", "degrees", FieldRequirement::Required, AccessIntent::Read,
                                            PersistencePolicy::Persistent),
                              host_field_2d("LON", "degrees", FieldRequirement::Required, AccessIntent::Read,
                                            PersistencePolicy::Persistent),
                              host_field_3d("T", "K"), host_field_3d("PMID", "Pa", FieldRequirement::Optional),
                              host_field_interface("PEDGE", "Pa", FieldRequirement::Optional),
                              host_field_3d("BXHEIGHT", "m"),
                              host_field_3d("AIRDEN", "kg/m3", FieldRequirement::Optional),
                              host_field_3d("AIRDEN_DRY", "kg/m3", FieldRequirement::Optional), host_concentration()},
                             {{"photolysis.ozone", "", true}});
    }

    PhotolysisProcess::PhotolysisProcess() : config_path("") {}
    PhotolysisProcess::~PhotolysisProcess() = default;

    void PhotolysisProcess::init(std::shared_ptr<StateManager> state) {
        Logger::debug(state.get(), "PhotolysisProcess::init started");
        if (state->config_manager()) {
            std::string cfg = state->config_manager()->get_string("processes/photolysis/config_file", "");
            if (cfg.empty()) {
                cfg = state->config_manager()->get_string("process/photolysis/config_file", "");
            }
            if (!cfg.empty()) {
                this->config_path = cfg;
            }
        }

        if (this->config_path.empty()) {
            this->config_path = "src/external/musica/configs/tuvx/tuv_5_4.yml";
        }

        Logger::debug(state.get(), "Parsing TUV-x config file to check profiles", {{"config", config_path}});
        std::unordered_set<std::string> config_defined_profiles;
        try {
            YAML::Node tuvx_config = YAML::LoadFile(config_path);
            if (tuvx_config["profiles"]) {
                for (const auto& p : tuvx_config["profiles"]) {
                    if (p["name"]) {
                        config_defined_profiles.insert(p["name"].as<std::string>());
                    }
                }
            }
        } catch (const std::exception& e) {
            Logger::debug(state.get(), "YAML parse warning", {{"error", e.what()}});
        }

        Logger::debug(state.get(), "Creating GridMap, ProfileMap, and RadiatorMap");
        musica::Error err;
        grids = musica::CreateGridMap(&err);
        profiles = musica::CreateProfileMap(&err);
        radiators = musica::CreateRadiatorMap(&err);

        Logger::debug(state.get(), "Creating height grid");
        musica::Grid* height_grid = musica::CreateGrid("height", "km", state->level_count(), &err);
        std::vector<double> dummy_edges(state->level_count() + 1, 0.0);
        std::vector<double> dummy_mids(state->level_count(), 0.0);
        for (int i = 0; i <= state->level_count(); ++i) {
            dummy_edges[i] = i * 1.0;
            if (i < state->level_count()) {
                dummy_mids[i] = i * 1.0 + 0.5;
            }
        }
        musica::SetGridEdges(height_grid, dummy_edges.data(), dummy_edges.size(), &err);
        musica::SetGridMidpoints(height_grid, dummy_mids.data(), dummy_mids.size(), &err);

        Logger::debug(state.get(), "Adding height grid to GridMap");
        musica::AddGrid(grids, height_grid, &err);

        Logger::debug(state.get(), "Selecting and configuring wavelength grid");

        static constexpr std::array<double, 6> edges_from_host = {300.0, 400.0, 500.0, 600.0, 700.0, 800.0};

        static constexpr std::array<double, 157> edges_standard = {
            120.0000, 121.4000, 121.9000, 122.3000, 123.1000, 123.8000, 124.6000, 125.4000, 126.2000, 127.0000,
            128.6000, 129.4000, 130.3000, 132.0000, 135.0000, 137.0000, 145.0000, 155.0000, 165.0000, 170.0000,
            175.4000, 177.0000, 178.6000, 180.2000, 181.8000, 183.5000, 185.2000, 186.9000, 188.7000, 190.5000,
            192.3000, 194.2000, 196.1000, 198.0000, 200.0000, 202.0000, 204.1000, 206.2000, 208.3330, 210.5260,
            212.7660, 215.0540, 217.3910, 219.7800, 222.2220, 224.7190, 227.2730, 229.8850, 232.5580, 235.2940,
            238.0950, 240.9640, 243.9020, 246.9140, 250.0000, 253.1650, 256.4100, 259.7400, 263.1580, 266.6670,
            270.2700, 273.9730, 277.7780, 281.6900, 285.7140, 289.8550, 294.1180, 298.5000, 302.5000, 303.5000,
            304.5000, 305.5000, 306.5000, 307.5000, 308.5000, 309.5000, 310.5000, 311.5000, 312.5000, 313.5000,
            314.5000, 317.5000, 322.5000, 327.5000, 332.5000, 337.5000, 342.5000, 347.5000, 352.5000, 357.5000,
            362.5000, 367.5000, 372.5000, 377.5000, 382.5000, 387.5000, 392.5000, 397.5000, 402.5000, 407.5000,
            412.5000, 417.5000, 422.5000, 427.5000, 432.5000, 437.5000, 442.5000, 447.5000, 452.5000, 457.5000,
            462.5000, 467.5000, 472.5000, 477.5000, 482.5000, 487.5000, 492.5000, 497.5000, 502.5000, 507.5000,
            512.5000, 517.5000, 522.5000, 527.5000, 532.5000, 537.5000, 542.5000, 547.5000, 552.5000, 557.5000,
            562.5000, 567.5000, 572.5000, 577.5000, 582.5000, 587.5000, 592.5000, 597.5000, 602.5000, 607.5000,
            612.5000, 617.5000, 622.5000, 627.5000, 632.5000, 637.5000, 642.5000, 647.1000, 655.0000, 665.0000,
            675.0000, 685.0000, 695.0000, 705.0000, 715.0000, 725.0000, 735.0000};

        int wl_sections = 156;
        std::vector<double> wl_edges;

        if (config_path.find("from_host") != std::string::npos ||
            config_path.find("config.json") != std::string::npos) {
            wl_sections = 5;
            wl_edges.assign(edges_from_host.begin(), edges_from_host.end());
        } else {
            wl_sections = 156;
            wl_edges.assign(edges_standard.begin(), edges_standard.end());
        }

        Logger::debug(state.get(), "Creating wavelength grid with specified sections",
                      {{"sections", std::to_string(wl_sections)}});
        musica::Grid* wl_grid = musica::CreateGrid("wavelength", "nm", wl_sections, &err);
        std::vector<double> wl_mids(wl_sections, 0.0);
        for (int i = 0; i < wl_sections; ++i) {
            wl_mids[i] = 0.5 * (wl_edges[i] + wl_edges[i + 1]);
        }
        musica::SetGridEdges(wl_grid, wl_edges.data(), wl_edges.size(), &err);
        musica::SetGridMidpoints(wl_grid, wl_mids.data(), wl_mids.size(), &err);

        Logger::debug(state.get(), "Adding wavelength grid to GridMap");
        musica::AddGrid(grids, wl_grid, &err);

        // 3. Register profiles safely only if missing from the config file definition
        register_profile_if_missing(state.get(), config_defined_profiles, "temperature", "K", height_grid, 280.0,
                                    state->level_count(), &err);
        register_profile_if_missing(state.get(), config_defined_profiles, "air", "molecule cm-3", height_grid, 1e12,
                                    state->level_count(), &err);
        register_profile_if_missing(state.get(), config_defined_profiles, "O2", "molecule cm-3", height_grid, 1e12,
                                    state->level_count(), &err);
        register_profile_if_missing(state.get(), config_defined_profiles, "O3", "molecule cm-3", height_grid, 1e12,
                                    state->level_count(), &err);
        register_profile_if_missing(state.get(), config_defined_profiles, "surface albedo", "none", wl_grid, 0.1,
                                    wl_sections, &err);
        register_profile_if_missing(state.get(), config_defined_profiles, "extraterrestrial flux", "photon cm-2 s-1",
                                    wl_grid, 1.5e14, wl_sections, &err);

        // 4. Safely delete local grids as they are cloned/owned inside the GridMap
        Logger::debug(state.get(), "Deleting local height and wavelength grid pointers");
        musica::DeleteGrid(height_grid, &err);
        musica::DeleteGrid(wl_grid, &err);

        // 5. Create TUVX instance using C API
        Logger::debug(state.get(), "Calling musica::CreateTuvx", {{"config", config_path}});
        tuvx_instance = musica::CreateTuvx(config_path.c_str(), grids, profiles, radiators, &err);

        if (err.code_ != 0) {
            std::cerr << "PhotolysisProcess: Error: Failed to initialize TUV-x! "
                      << (err.message_.value_ ? err.message_.value_ : "Unknown Error") << std::endl;
            return;
        }

        Logger::info(state.get(), "PhotolysisProcess: initialized TUV-x successfully!");

        Logger::debug(state.get(), "Getting Photolysis rate constants ordering");
        musica::GetPhotolysisRateConstantsOrdering(tuvx_instance, &photo_mappings, &err);

        Logger::debug(state.get(), "Dynamic diagnostic field registration");
        if (state->diagnostic_manager()) {
            std::vector<int> dims_2d = {state->column_count(), state->level_count()};
            // J-rate fields are rewritten in full (every column x level) on
            // every step, so the blanket per-step reset that register_field()
            // (DiagnosticPolicy::Instantaneous) performs is pure overhead.  They
            // are registered Persistent instead; run() zeroes a column only on a
            // TUV-x solver error, which preserves the Instantaneous semantics
            // (unwritten cells read back as 0) without the per-step memset.
            const std::vector<SemanticAxis> axes_2d = {SemanticAxis::Column, SemanticAxis::Level};
            for (size_t i = 0; i < photo_mappings.size_; ++i) {
                std::string rx_name =
                    photo_mappings.mappings_[i].name_.value_ ? photo_mappings.mappings_[i].name_.value_ : "";
                std::string diag_name = "photolysis_rate_" + rx_name;
                state->diagnostic_manager()->register_field_contract(diag_name, "Photolysis rate for " + rx_name, "s-1",
                                                                     DiagType::FIELD_2D, dims_2d,
                                                                     DiagnosticPolicy::Persistent, 0.0, axes_2d);
                // Persistent fields never see the blanket reset, so seed the
                // storage explicitly once; the per-step reset this replaces
                // would otherwise have been the only zeroing pass.
                state->diagnostic_manager()->get_field(diag_name)->reset();
            }
        }
        Logger::debug(state.get(), "PhotolysisProcess::init complete");
    }

    void PhotolysisProcess::run(std::shared_ptr<StateManager> state) {
        Logger::debug(state.get(), "PhotolysisProcess::run started");
        if (!tuvx_instance) {
            Logger::debug(state.get(), "tuvx_instance is NULL, returning");
            return;
        }

        Logger::debug(state.get(), "Syncing state to host");

        Logger::debug(state.get(), "Resolving configured ozone role");
        if (!state->chemistry().mechanism || !state->chemistry().mechanism->has_role("photolysis.ozone"))
            throw std::runtime_error("Photolysis requires mechanism role photolysis.ozone");
        const int i_o3 = static_cast<int>(state->chemistry().mechanism->index_for_role("photolysis.ozone"));
        Logger::debug(state.get(), "Ozone index resolved", {{"index", std::to_string(i_o3)}});

        musica::Error err;
        int num_reactions = photo_mappings.size_;
        Logger::debug(state.get(), "Number of photolysis reactions mapped",
                      {{"reactions", std::to_string(num_reactions)}});

        Logger::debug(state.get(), "Fetching ProfileMap from TUVX instance");
        musica::ProfileMap* loaded_profiles = musica::GetProfileMap(tuvx_instance, &err);
        if (err.code_ != 0) {
            std::cerr << "PhotolysisProcess: Error getting ProfileMap: "
                      << (err.message_.value_ ? err.message_.value_ : "Unknown Error") << std::endl;
            return;
        }

        Logger::debug(state.get(), "Retrieving individual Profile pointers");
        musica::Profile* profile_air = musica::GetProfile(loaded_profiles, "air", "molecule cm-3", &err);
        musica::Profile* profile_o2 = musica::GetProfile(loaded_profiles, "O2", "molecule cm-3", &err);
        musica::Profile* profile_o3 = musica::GetProfile(loaded_profiles, "O3", "molecule cm-3", &err);
        musica::Profile* profile_temp = musica::GetProfile(loaded_profiles, "temperature", "K", &err);

        std::ostringstream ptr_ss;
        ptr_ss << "air=" << profile_air << ", o2=" << profile_o2 << ", o3=" << profile_o3 << ", temp=" << profile_temp;
        Logger::debug(state.get(), "Profile pointers retrieved", {{"profiles", ptr_ss.str()}});

        if (profile_air) {
            Logger::debug(state.get(), "profile_air name", {{"name", profile_air->GetName(&err)}});
        }
        if (profile_o2) {
            Logger::debug(state.get(), "profile_o2 name", {{"name", profile_o2->GetName(&err)}});
        }

        Logger::debug(state.get(), "Retrieving GridMap from TUVX instance");
        musica::GridMap* loaded_grids = musica::GetGridMap(tuvx_instance, &err);
        musica::Grid* height_grid = musica::GetGrid(loaded_grids, "height", "km", &err);

        std::ostringstream grid_ss;
        grid_ss << height_grid;
        Logger::debug(state.get(), "height_grid retrieved", {{"ptr", grid_ss.str()}});

        std::vector<double> height_edges(state->level_count() + 1, 0.0);
        std::vector<double> air_profile(state->level_count(), 0.0);
        std::vector<double> o2_profile(state->level_count(), 0.0);
        std::vector<double> o3_profile(state->level_count(), 0.0);
        std::vector<double> temp_profile(state->level_count(), 0.0);

        // Derived fields remain allocated between imports, so require the
        // current generation rather than using pointer existence as the
        // freshness test.
        if (state->meteorology().PEDGE && state->meteorology().T) {
            state->derive_bxheight();
        }
        if (state->meteorology().PMID && state->meteorology().T) {
            state->derive_airden_dry();
        }

        const auto import_generation = state->current_import_generation();

        require_field_pointer("Photolysis", "LAT",
                              state->meteorology().LAT ? state->meteorology().LAT->host_write() : nullptr);
        require_field_pointer("Photolysis", "LON",
                              state->meteorology().LON ? state->meteorology().LON->host_write() : nullptr);
        require_field_pointer("Photolysis", "BXHEIGHT",
                              state->meteorology().BXHEIGHT &&
                                      state->meteorology().BXHEIGHT->is_current(import_generation)
                                  ? state->meteorology().BXHEIGHT->host_read()
                                  : nullptr);
        require_field_pointer("Photolysis", "AIRDEN_DRY",
                              state->meteorology().AIRDEN_DRY &&
                                      state->meteorology().AIRDEN_DRY->is_current(import_generation)
                                  ? state->meteorology().AIRDEN_DRY->host_read()
                                  : nullptr);
        require_field_pointer("Photolysis", "T",
                              state->meteorology().T && state->meteorology().T->is_current(import_generation)
                                  ? state->meteorology().T->host_read()
                                  : nullptr);
        require_field_pointer("Photolysis", "CHEM_CONC",
                              state->chemistry().conc ? state->chemistry().conc->host_write() : nullptr);

        // Column-level tracing is only useful when CATCHEM_LOG_LEVEL=DEBUG; the
        // flag is hoisted out of the loop so the per-column string build and the
        // ten Logger calls below cost nothing in nominal runs.
        const bool dbg_column = Logger::enabled(Logger::Level::Debug);
        Logger::debug(state.get(), "Starting column-wise calculation loop");
        for (int i_col = 0; i_col < state->column_count(); ++i_col) {
            std::string col_str = dbg_column ? std::to_string(i_col) : std::string();
            if (dbg_column)
                Logger::debug(state.get(), "Calculating SZA for column", {{"col", col_str}});
            double lat_deg = state->meteorology().LAT->host_view(i_col, 0);
            double lon_deg = state->meteorology().LON->host_view(i_col, 0);
            double cos_sza = state->clock().get_cos_sza(lat_deg, lon_deg, true);
            double sza_rad = std::acos(std::max(-1.0, std::min(1.0, cos_sza)));

            if (dbg_column)
                Logger::debug(state.get(), "Updating grid height edges for column", {{"col", col_str}});
            height_edges[0] = 0.0;
            for (int i_lvl = 0; i_lvl < state->level_count(); ++i_lvl) {
                double dz_m = state->meteorology().BXHEIGHT->host_view(i_col, i_lvl, 0);
                height_edges[i_lvl + 1] = height_edges[i_lvl] + dz_m / 1000.0;
            }
            if (height_grid) {
                musica::SetGridEdges(height_grid, height_edges.data(), height_edges.size(), &err);
            }

            if (dbg_column)
                Logger::debug(state.get(), "Populating profile midpoint vectors for column", {{"col", col_str}});
            for (int i_lvl = 0; i_lvl < state->level_count(); ++i_lvl) {
                double airden_kg_m3 = state->meteorology().AIRDEN_DRY->host_view(i_col, i_lvl, 0);
                air_profile[i_lvl] = airden_kg_m3 * 2.079153e19;
                o2_profile[i_lvl] = air_profile[i_lvl] * 0.2095;
                temp_profile[i_lvl] = state->meteorology().T->host_view(i_col, i_lvl, 0);

                if (i_o3 >= 0 && state->chemistry().conc) {
                    o3_profile[i_lvl] = state->chemistry().conc->host_view(i_col, i_lvl, i_o3);
                } else {
                    o3_profile[i_lvl] = air_profile[i_lvl] * 3e-7;
                }
            }

            if (dbg_column)
                Logger::debug(state.get(), "Updating profiles in TUVX for column", {{"col", col_str}});
            if (profile_air) {
                if (dbg_column)
                    Logger::debug(state.get(), "SetProfileMidpointValues for air in column", {{"col", col_str}});
                musica::SetProfileMidpointValues(profile_air, air_profile.data(), state->level_count(), &err);
            }
            if (profile_o2) {
                if (dbg_column)
                    Logger::debug(state.get(), "SetProfileMidpointValues for O2 in column", {{"col", col_str}});
                musica::SetProfileMidpointValues(profile_o2, o2_profile.data(), state->level_count(), &err);
            }
            if (profile_o3) {
                if (dbg_column)
                    Logger::debug(state.get(), "SetProfileMidpointValues for O3 in column", {{"col", col_str}});
                musica::SetProfileMidpointValues(profile_o3, o3_profile.data(), state->level_count(), &err);
            }
            if (profile_temp) {
                if (dbg_column)
                    Logger::debug(state.get(), "SetProfileMidpointValues for temperature in column",
                                  {{"col", col_str}});
                musica::SetProfileMidpointValues(profile_temp, temp_profile.data(), state->level_count(), &err);
            }

            std::vector<double> edge_photolysis_rates((state->level_count() + 1) * num_reactions, 0.0);
            std::vector<double> edge_heating_rates((state->level_count() + 1) * tuvx_instance->GetHeatingRateCount(),
                                                   0.0);

            if (dbg_column)
                Logger::debug(state.get(), "Calling musica::RunTuvx for column", {{"col", col_str}});
            musica::RunTuvx(tuvx_instance, sza_rad, 1.0, edge_photolysis_rates.data(), edge_heating_rates.data(),
                            nullptr, nullptr, nullptr, &err);

            if (err.code_ != 0) {
                std::cerr << "PhotolysisProcess: Solver error in column " << i_col << ": "
                          << (err.message_.value_ ? err.message_.value_ : "Unknown Error") << std::endl;
                // The write loop below is skipped for this column, so clear its
                // slice of every J-rate field.  This reproduces the zeroed cell
                // that Instantaneous reset used to guarantee, now that the
                // fields are Persistent and carry the previous step's values.
                if (state->diagnostic_manager()) {
                    for (size_t rx_idx = 0; rx_idx < photo_mappings.size_; ++rx_idx) {
                        std::string rx_name = photo_mappings.mappings_[rx_idx].name_.value_
                                                  ? photo_mappings.mappings_[rx_idx].name_.value_
                                                  : "";
                        double* diag_ptr = static_cast<double*>(
                            state->diagnostic_manager()->get_host_pointer("photolysis_rate_" + rx_name));
                        if (diag_ptr) {
                            for (int i_lvl = 0; i_lvl < state->level_count(); ++i_lvl)
                                diag_ptr[i_lvl * state->column_count() + i_col] = 0.0;
                        }
                    }
                }
                continue;
            }

            if (dbg_column)
                Logger::debug(state.get(), "Copying midpoint-interpolated J-rates to diagnostics for column",
                              {{"col", col_str}});
            if (state->diagnostic_manager()) {
                for (size_t rx_idx = 0; rx_idx < photo_mappings.size_; ++rx_idx) {
                    std::string rx_name = photo_mappings.mappings_[rx_idx].name_.value_
                                              ? photo_mappings.mappings_[rx_idx].name_.value_
                                              : "";
                    std::string diag_name = "photolysis_rate_" + rx_name;
                    double* diag_ptr = static_cast<double*>(state->diagnostic_manager()->get_host_pointer(diag_name));

                    if (diag_ptr) {
                        for (int i_lvl = 0; i_lvl < state->level_count(); ++i_lvl) {
                            int idx_edge1 = rx_idx * (state->level_count() + 1) + i_lvl;
                            int idx_edge2 = rx_idx * (state->level_count() + 1) + (i_lvl + 1);

                            double rate_midpoint =
                                0.5 * (edge_photolysis_rates[idx_edge1] + edge_photolysis_rates[idx_edge2]);

                            int diag_idx = i_lvl * state->column_count() + i_col;
                            diag_ptr[diag_idx] = rate_midpoint;
                        }
                    }
                }
            }
        }

        Logger::debug(state.get(), "Syncing diagnostics and state to device");
        Logger::debug(state.get(), "PhotolysisProcess::run complete");
    }

    void PhotolysisProcess::finalize() {
        musica::Error err;
        if (tuvx_instance) {
            musica::DeleteTuvx(tuvx_instance, &err);
            tuvx_instance = nullptr;
        }
        if (grids) {
            musica::DeleteGridMap(grids, &err);
            grids = nullptr;
        }
        if (profiles) {
            musica::DeleteProfileMap(profiles, &err);
            profiles = nullptr;
        }
        if (radiators) {
            musica::DeleteRadiatorMap(radiators, &err);
            radiators = nullptr;
        }
    }

    void PhotolysisProcess::register_profile_if_missing(const StateManager* state,
                                                        const std::unordered_set<std::string>& config_defined_profiles,
                                                        const char* name, const char* units, musica::Grid* grid,
                                                        double default_val, std::size_t num_vals, musica::Error* err) {
        if (config_defined_profiles.find(name) == config_defined_profiles.end()) {
            Logger::debug(state, "Pre-registering missing profile", {{"name", name}, {"units", units}});
            musica::Profile* new_prof = musica::CreateProfile(name, units, grid, err);
            std::vector<double> dummy(num_vals, default_val);
            musica::SetProfileMidpointValues(new_prof, dummy.data(), num_vals, err);
            musica::AddProfile(profiles, new_prof, err);
            musica::DeleteProfile(new_prof, err);
        }
    }

} // namespace catchem

extern "C" {
void catchem_register_photolysis_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(std::string(catchem::ProcessNames::Photolysis), []() {
        return std::make_shared<catchem::PhotolysisProcess>();
    });
}
}
