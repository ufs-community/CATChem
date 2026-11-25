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
#include <musica/tuvx/radiator.hpp>
#include <musica/tuvx/radiator_map.hpp>
#include <musica/tuvx/tuvx_c_interface.hpp>
#include <sstream>
#include <string_view>
#include <unordered_set>
#include <yaml-cpp/yaml.h>

namespace catchem {

    ProcessContract PhotolysisProcess::get_contract() const {
        auto contract = make_contract(
            get_name(), {host_field_2d("LAT", "degrees", FieldRequirement::Required, AccessIntent::Read,
                                       PersistencePolicy::Persistent),
                         host_field_2d("LON", "degrees", FieldRequirement::Required, AccessIntent::Read,
                                       PersistencePolicy::Persistent),
                         host_field_3d("T", "K"), host_field_3d("PMID", "Pa", FieldRequirement::Optional),
                         host_field_interface("TEDGE", "K", FieldRequirement::Optional),
                         host_field_interface("PEDGE", "Pa", FieldRequirement::Optional),
                         host_field_3d("BXHEIGHT", "m"), host_field_3d("AIRDEN", "kg/m3", FieldRequirement::Optional),
                         host_field_3d("AIRDEN_DRY", "kg/m3", FieldRequirement::Optional), host_concentration()});
        if (cloud_feedback_) {
            contract.fields.push_back(host_field_3d("CLDF", "1"));
            contract.fields.push_back(host_field_3d("CLDLIQ", "kg/kg"));
            contract.fields.push_back(host_field_3d("DELP", "Pa"));
        }
        return contract;
    }

    PhotolysisProcess::PhotolysisProcess() : config_path_("") {}
    PhotolysisProcess::~PhotolysisProcess() = default;

    void PhotolysisProcess::init(std::shared_ptr<StateManager> state) {
        Logger::debug(*state, "PhotolysisProcess::init started");
        if (state->config_manager()) {
            cloud_feedback_ = state->config_manager()->get_bool(
                "processes/photolysis/cloud_feedback",
                state->config_manager()->get_bool("process/photolysis/cloud_feedback", false));
            std::string cfg = state->config_manager()->get_string("processes/photolysis/config_file", "");
            if (cfg.empty()) {
                cfg = state->config_manager()->get_string("process/photolysis/config_file", "");
            }
            if (!cfg.empty()) {
                this->config_path_ = cfg;
            }
        }

        if (this->config_path_.empty()) {
            this->config_path_ = "src/external/musica/configs/tuvx/tuv_5_4.yml";
        }

        const auto& mechanism = state->chemistry().mechanism;
        if (!mechanism)
            throw std::runtime_error("Photolysis requires an active chemical mechanism");
        species_profiles_.clear();
        std::unordered_set<std::string> mapped_profiles;
        constexpr std::string_view role_prefix = "PHOTOLYSIS.";
        for (const auto& species : mechanism->species) {
            for (const auto& role : species.roles) {
                const auto canonical_role = canonical_species_name(role);
                if (!canonical_role.starts_with(role_prefix))
                    continue;
                const bool legacy_ozone = canonical_role == "PHOTOLYSIS.OZONE";
                std::string name = legacy_ozone ? "O3" : role.substr(role_prefix.size());
                const auto canonical_name = canonical_species_name(name);
                if (!species.is_gas || name.empty() || canonical_name == "AIR" || canonical_name == "TEMPERATURE" ||
                    canonical_name == "SURFACE ALBEDO" || canonical_name == "EXTRATERRESTRIAL FLUX" ||
                    !mapped_profiles.insert(name).second)
                    throw std::invalid_argument("Invalid or duplicate photolysis species profile for role: " + role);
                species_profiles_.push_back({role, std::move(name), legacy_ozone});
            }
        }

        Logger::debug(*state, "Parsing TUV-x config file to check profiles", {{"config", config_path_}});
        std::unordered_set<std::string> config_defined_profiles;
        try {
            YAML::Node tuvx_config = YAML::LoadFile(config_path_);
            if (tuvx_config["profiles"]) {
                for (const auto& p : tuvx_config["profiles"]) {
                    if (p["name"]) {
                        config_defined_profiles.insert(p["name"].as<std::string>());
                    }
                }
            }
        } catch (const std::exception& e) {
            Logger::debug(*state, "YAML parse warning", {{"error", e.what()}});
        }

        Logger::debug(*state, "Creating GridMap, ProfileMap, and RadiatorMap");
        musica::Error err;
        grids_ = musica::CreateGridMap(&err);
        profiles_ = musica::CreateProfileMap(&err);
        radiators_ = musica::CreateRadiatorMap(&err);

        Logger::debug(*state, "Creating height grid");
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

        Logger::debug(*state, "Adding height grid to GridMap");
        musica::AddGrid(grids_, height_grid, &err);

        Logger::debug(*state, "Selecting and configuring wavelength grid");

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

        if (config_path_.find("from_host") != std::string::npos ||
            config_path_.find("config.json") != std::string::npos) {
            wl_sections = 5;
            wl_edges.assign(edges_from_host.begin(), edges_from_host.end());
        } else {
            wl_sections = 156;
            wl_edges.assign(edges_standard.begin(), edges_standard.end());
        }

        Logger::debug(*state, "Creating wavelength grid with specified sections",
                      {{"sections", std::to_string(wl_sections)}});
        musica::Grid* wl_grid = musica::CreateGrid("wavelength", "nm", wl_sections, &err);
        std::vector<double> wl_mids(wl_sections, 0.0);
        for (int i = 0; i < wl_sections; ++i) {
            wl_mids[i] = 0.5 * (wl_edges[i] + wl_edges[i + 1]);
        }
        musica::SetGridEdges(wl_grid, wl_edges.data(), wl_edges.size(), &err);
        musica::SetGridMidpoints(wl_grid, wl_mids.data(), wl_mids.size(), &err);

        Logger::debug(*state, "Adding wavelength grid to GridMap");
        musica::AddGrid(grids_, wl_grid, &err);

        // 3. Register profiles safely only if missing from the config file definition
        register_profile_if_missing(*state, config_defined_profiles, "temperature", "K", height_grid, 280.0,
                                    state->level_count(), &err);
        register_profile_if_missing(*state, config_defined_profiles, "air", "molecule cm-3", height_grid, 1e12,
                                    state->level_count(), &err);
        for (const auto& binding : species_profiles_) {
            register_profile_if_missing(*state, config_defined_profiles, binding.name.c_str(), "molecule cm-3",
                                        height_grid, 1e12, state->level_count(), &err);
        }
        for (const auto* background : {"O2", "O3"}) {
            if (!mapped_profiles.contains(background))
                register_profile_if_missing(*state, config_defined_profiles, background, "molecule cm-3", height_grid,
                                            1e12, state->level_count(), &err);
        }
        register_profile_if_missing(*state, config_defined_profiles, "surface albedo", "none", wl_grid, 0.1,
                                    wl_sections, &err);
        register_profile_if_missing(*state, config_defined_profiles, "extraterrestrial flux", "photon cm-2 s-1",
                                    wl_grid, 1.5e14, wl_sections, &err);

        if (cloud_feedback_) {
            std::unique_ptr<musica::Radiator> cloud(musica::CreateRadiator("clouds", height_grid, wl_grid, &err));
            if (err.code_ != 0 || !cloud)
                throw std::runtime_error("Photolysis failed to create cloud radiator");
            musica::AddRadiator(radiators_, cloud.get(), &err);
            if (err.code_ != 0)
                throw std::runtime_error("Photolysis failed to register cloud radiator");
        }

        // 4. Safely delete local grids as they are cloned/owned inside the GridMap
        Logger::debug(*state, "Deleting local height and wavelength grid pointers");
        musica::DeleteGrid(height_grid, &err);
        musica::DeleteGrid(wl_grid, &err);

        // 5. Create TUVX instance using C API
        Logger::debug(*state, "Calling musica::CreateTuvx", {{"config", config_path_}});
        tuvx_instance_ = musica::CreateTuvx(config_path_.c_str(), grids_, profiles_, radiators_, &err);

        if (err.code_ != 0) {
            Logger::error(*state, "PhotolysisProcess: failed to initialize TUV-x",
                          {{"error", err.message_.value_ ? err.message_.value_ : "Unknown Error"}});
            return;
        }

        Logger::info(*state, "PhotolysisProcess: initialized TUV-x successfully!");

        Logger::debug(*state, "Getting Photolysis rate constants ordering");
        musica::GetPhotolysisRateConstantsOrdering(tuvx_instance_, &photo_mappings_, &err);

        Logger::debug(*state, "Dynamic diagnostic field registration");
        if (state->diagnostic_manager()) {
            std::vector<int> dims_2d = {state->column_count(), state->level_count()};
            // J-rate fields are rewritten in full (every column x level) on
            // every step, so the blanket per-step reset that register_field()
            // (DiagnosticPolicy::Instantaneous) performs is pure overhead.  They
            // are registered Persistent instead; run() zeroes a column only on a
            // TUV-x solver error, which preserves the Instantaneous semantics
            // (unwritten cells read back as 0) without the per-step memset.
            const std::vector<SemanticAxis> axes_2d = {SemanticAxis::Column, SemanticAxis::Level};
            for (size_t i = 0; i < photo_mappings_.size_; ++i) {
                std::string rx_name =
                    photo_mappings_.mappings_[i].name_.value_ ? photo_mappings_.mappings_[i].name_.value_ : "";
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
        Logger::debug(*state, "PhotolysisProcess::init complete");
    }

    void PhotolysisProcess::run(std::shared_ptr<StateManager> state) {
        Logger::debug(*state, "PhotolysisProcess::run started");
        if (!tuvx_instance_) {
            Logger::debug(*state, "tuvx_instance is NULL, returning");
            return;
        }

        Logger::debug(*state, "Syncing state to host");

        musica::Error err;
        int num_reactions = photo_mappings_.size_;
        Logger::debug(*state, "Number of photolysis reactions mapped", {{"reactions", std::to_string(num_reactions)}});

        Logger::debug(*state, "Fetching ProfileMap from TUVX instance");
        const std::unique_ptr<musica::ProfileMap> loaded_profiles(musica::GetProfileMap(tuvx_instance_, &err));
        if (err.code_ != 0) {
            Logger::error(*state, "PhotolysisProcess: error getting ProfileMap",
                          {{"error", err.message_.value_ ? err.message_.value_ : "Unknown Error"}});
            return;
        }

        Logger::debug(*state, "Retrieving individual Profile pointers");
        const auto density_profile = [&](const char* name) {
            return std::unique_ptr<musica::Profile>(
                musica::GetProfile(loaded_profiles.get(), name, "molecule cm-3", &err));
        };
        const auto is_mapped = [&](std::string_view name) {
            return std::ranges::any_of(species_profiles_, [&](const auto& binding) { return binding.name == name; });
        };
        const auto profile_air = density_profile("air");
        const auto profile_o2 = is_mapped("O2") ? nullptr : density_profile("O2");
        const auto profile_o3 = is_mapped("O3") ? nullptr : density_profile("O3");
        const std::unique_ptr<musica::Profile> profile_temp(
            musica::GetProfile(loaded_profiles.get(), "temperature", "K", &err));

        std::vector<std::pair<std::size_t, std::unique_ptr<musica::Profile>>> host_species;
        for (const auto& binding : species_profiles_) {
            auto profile = density_profile(binding.name.c_str());
            if (!profile && !binding.legacy_ozone)
                throw std::invalid_argument("Photolysis role " + binding.role +
                                            " requires a host-updatable TUV-x profile in molecule cm-3");
            if (profile) {
                const auto index = state->chemistry().mechanism->index_for_role(binding.role);
                host_species.emplace_back(index, std::move(profile));
            }
        }

        std::ostringstream ptr_ss;
        ptr_ss << "air=" << profile_air.get() << ", o2=" << profile_o2.get() << ", o3=" << profile_o3.get()
               << ", temp=" << profile_temp.get();
        Logger::debug(*state, "Profile pointers retrieved", {{"profiles", ptr_ss.str()}});

        if (profile_air) {
            Logger::debug(*state, "profile_air name", {{"name", profile_air->GetName(&err)}});
        }
        if (profile_o2) {
            Logger::debug(*state, "profile_o2 name", {{"name", profile_o2->GetName(&err)}});
        }

        Logger::debug(*state, "Retrieving GridMap from TUVX instance");
        musica::GridMap* loaded_grids = musica::GetGridMap(tuvx_instance_, &err);
        musica::Grid* height_grid = musica::GetGrid(loaded_grids, "height", "km", &err);

        std::ostringstream grid_ss;
        grid_ss << height_grid;
        Logger::debug(*state, "height_grid retrieved", {{"ptr", grid_ss.str()}});

        std::vector<double> height_edges(state->level_count() + 1, 0.0);
        std::vector<double> air_profile(state->level_count(), 0.0);
        std::vector<double> o2_profile(state->level_count(), 0.0);
        std::vector<double> o3_profile(state->level_count(), 0.0);
        std::vector<double> species_density(state->level_count());
        std::vector<double> temp_profile(state->level_count(), 0.0);
        std::vector<double> temp_edges(state->level_count() + 1);
        std::vector<double> layer_density(state->level_count());

        const auto update_density_profile = [&](musica::Profile* profile, std::vector<double>& midpoint_density) {
            if (!profile)
                return;
            // TUV-x uses the integrated column in each layer for optical depth.
            for (std::size_t level = 0; level < midpoint_density.size(); ++level) {
                const double thickness_cm = (height_edges[level + 1] - height_edges[level]) * 1.0e5;
                layer_density[level] = midpoint_density[level] * thickness_cm;
            }
            musica::SetProfileMidpointValues(profile, midpoint_density.data(), midpoint_density.size(), &err);
            musica::SetProfileLayerDensities(profile, layer_density.data(), layer_density.size(), &err);
        };

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
        std::unique_ptr<musica::RadiatorMap> cloud_map;
        std::unique_ptr<musica::Radiator> cloud_radiator;
        std::shared_ptr<InteropField<double, 3>> cloud_fraction, cloud_liquid, pressure_thickness;
        if (cloud_feedback_) {
            cloud_fraction = state->meteorology().get_3d({"CLDF"});
            cloud_liquid = state->meteorology().get_3d({"CLDLIQ"});
            pressure_thickness = state->meteorology().get_3d({"DELP"});
            for (const auto& field : {cloud_fraction, cloud_liquid, pressure_thickness}) {
                if (!field || !field->is_current(import_generation))
                    throw std::runtime_error("Cloud feedback requires current layer fields CLDF, CLDLIQ and DELP");
                field->sync_to_host();
            }
            cloud_map.reset(musica::GetRadiatorMap(tuvx_instance_, &err));
            if (err.code_ != 0 || !cloud_map)
                throw std::runtime_error("Photolysis failed to get cloud radiator map");
            cloud_radiator.reset(musica::GetRadiator(cloud_map.get(), "clouds", &err));
            if (err.code_ != 0 || !cloud_radiator)
                throw std::runtime_error("Photolysis failed to get cloud radiator");
        }

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

        const auto host_temp_edges = state->meteorology().get_3d({"TEDGE"});
        const bool use_host_temp_edges =
            profile_temp && host_temp_edges && host_temp_edges->is_current(import_generation);
        if (use_host_temp_edges) {
            if (host_temp_edges->extent(1) != temp_edges.size() ||
                host_temp_edges->contract.axes[1] != SemanticAxis::Interface)
                throw std::invalid_argument("Photolysis expects TEDGE on the vertical interface axis");
            require_field_pointer("Photolysis", "TEDGE", host_temp_edges->host_read());
        }

        // Column-level tracing is only useful when CATCHEM_LOG_LEVEL=DEBUG; the
        // flag is hoisted out of the loop so the per-column string build and the
        // ten Logger calls below cost nothing in nominal runs.
        const bool dbg_column = Logger::enabled(Logger::Level::Debug);
        Logger::debug(*state, "Starting column-wise calculation loop");
        for (int i_col = 0; i_col < state->column_count(); ++i_col) {
            std::string col_str = dbg_column ? std::to_string(i_col) : std::string();
            if (dbg_column)
                Logger::debug(*state, "Calculating SZA for column", {{"col", col_str}});
            double lat_deg = state->meteorology().LAT->host_view(i_col, 0);
            double lon_deg = state->meteorology().LON->host_view(i_col, 0);
            double cos_sza = state->clock().get_cos_sza(lat_deg, lon_deg, true);
            double sza_rad = std::acos(std::max(-1.0, std::min(1.0, cos_sza)));

            if (dbg_column)
                Logger::debug(*state, "Updating grid height edges for column", {{"col", col_str}});
            height_edges[0] = 0.0;
            for (int i_lvl = 0; i_lvl < state->level_count(); ++i_lvl) {
                double dz_m = state->meteorology().BXHEIGHT->host_view(i_col, i_lvl, 0);
                height_edges[i_lvl + 1] = height_edges[i_lvl] + dz_m / 1000.0;
            }
            if (height_grid) {
                musica::SetGridEdges(height_grid, height_edges.data(), height_edges.size(), &err);
            }

            if (cloud_feedback_) {
                const int wavelength_count = tuvx_instance_->GetNumberOfWavelengthMidpoints();
                std::vector<double> optical_depth(state->level_count() * wavelength_count);
                for (int level = 0; level < state->level_count(); ++level) {
                    const double fraction = cloud_fraction->host_view(i_col, level, 0);
                    const double liquid = cloud_liquid->host_view(i_col, level, 0);
                    const double pressure = pressure_thickness->host_view(i_col, level, 0);
                    if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0 || !std::isfinite(liquid) ||
                        liquid < 0.0 || !std::isfinite(pressure) || pressure <= 0.0)
                        throw std::domain_error(
                            "Photolysis received invalid cloud fraction, liquid water or pressure thickness");
                    const double tau = 155.0 * liquid * pressure / constants::G0 * std::sqrt(fraction);
                    for (int wavelength = 0; wavelength < wavelength_count; ++wavelength)
                        optical_depth[wavelength * state->level_count() + level] = tau;
                }
                musica::SetRadiatorOpticalDepths(cloud_radiator.get(), optical_depth.data(), state->level_count(),
                                                 wavelength_count, &err);
                if (err.code_ != 0)
                    throw std::runtime_error("Photolysis failed to update cloud optical depth");
            }

            if (dbg_column)
                Logger::debug(*state, "Populating profile midpoint vectors for column", {{"col", col_str}});
            for (int i_lvl = 0; i_lvl < state->level_count(); ++i_lvl) {
                double airden_kg_m3 = state->meteorology().AIRDEN_DRY->host_view(i_col, i_lvl, 0);
                air_profile[i_lvl] = airden_kg_m3 * 2.079153e19;
                o2_profile[i_lvl] = air_profile[i_lvl] * 0.2095;
                temp_profile[i_lvl] = state->meteorology().T->host_view(i_col, i_lvl, 0);

                o3_profile[i_lvl] = air_profile[i_lvl] * 3e-7;
            }

            if (dbg_column)
                Logger::debug(*state, "Updating profiles in TUVX for column", {{"col", col_str}});
            update_density_profile(profile_air.get(), air_profile);
            update_density_profile(profile_o2.get(), o2_profile);
            update_density_profile(profile_o3.get(), o3_profile);
            for (const auto& [index, profile] : host_species) {
                for (int level = 0; level < state->level_count(); ++level) {
                    // A negative mixing ratio cannot contribute negative optical depth.
                    const double ppmv = state->chemistry().conc->host_view(i_col, level, index);
                    if (!std::isfinite(ppmv))
                        throw std::domain_error("Photolysis received a non-finite mixing ratio for species " +
                                                state->chemistry().mechanism->species[index].short_name);
                    species_density[level] = std::max(ppmv, 0.0) * 1.0e-6 * air_profile[level];
                }
                update_density_profile(profile.get(), species_density);
            }
            if (profile_temp) {
                if (dbg_column)
                    Logger::debug(*state, "Updating temperature profile in column", {{"col", col_str}});
                if (use_host_temp_edges) {
                    for (std::size_t edge = 0; edge < temp_edges.size(); ++edge)
                        temp_edges[edge] = host_temp_edges->host_view(i_col, edge, 0);
                } else {
                    // Hold the nearest midpoint value at the boundaries. Interior
                    // edges divide the midpoint interval in proportion to layer thickness.
                    temp_edges.front() = temp_profile.front();
                    temp_edges.back() = temp_profile.back();
                    for (std::size_t edge = 1; edge < temp_profile.size(); ++edge) {
                        const double lower_thickness = height_edges[edge] - height_edges[edge - 1];
                        const double span = height_edges[edge + 1] - height_edges[edge - 1];
                        temp_edges[edge] =
                            std::lerp(temp_profile[edge - 1], temp_profile[edge], lower_thickness / span);
                    }
                }
                if (std::ranges::any_of(temp_edges, [](double temp) { return !std::isfinite(temp) || temp <= 0.0; }))
                    throw std::domain_error("Photolysis requires finite, positive interface temperatures in kelvin");
                musica::SetProfileMidpointValues(profile_temp.get(), temp_profile.data(), state->level_count(), &err);
                musica::SetProfileEdgeValues(profile_temp.get(), temp_edges.data(), temp_edges.size(), &err);
            }

            std::vector<double> edge_photolysis_rates((state->level_count() + 1) * num_reactions, 0.0);
            std::vector<double> edge_heating_rates((state->level_count() + 1) * tuvx_instance_->GetHeatingRateCount(),
                                                   0.0);

            if (dbg_column)
                Logger::debug(*state, "Calling musica::RunTuvx for column", {{"col", col_str}});
            musica::RunTuvx(tuvx_instance_, sza_rad, 1.0, edge_photolysis_rates.data(), edge_heating_rates.data(),
                            nullptr, nullptr, nullptr, &err);

            if (err.code_ != 0) {
                Logger::error(*state, "PhotolysisProcess: solver error in column",
                              {{"column", std::to_string(i_col)},
                               {"error", err.message_.value_ ? err.message_.value_ : "Unknown Error"}});
                // The write loop below is skipped for this column, so clear its
                // slice of every J-rate field.  This reproduces the zeroed cell
                // that Instantaneous reset used to guarantee, now that the
                // fields are Persistent and carry the previous step's values.
                if (state->diagnostic_manager()) {
                    for (size_t rx_idx = 0; rx_idx < photo_mappings_.size_; ++rx_idx) {
                        std::string rx_name = photo_mappings_.mappings_[rx_idx].name_.value_
                                                  ? photo_mappings_.mappings_[rx_idx].name_.value_
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
                Logger::debug(*state, "Copying midpoint-interpolated J-rates to diagnostics for column",
                              {{"col", col_str}});
            if (state->diagnostic_manager()) {
                for (size_t rx_idx = 0; rx_idx < photo_mappings_.size_; ++rx_idx) {
                    std::string rx_name = photo_mappings_.mappings_[rx_idx].name_.value_
                                              ? photo_mappings_.mappings_[rx_idx].name_.value_
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

        Logger::debug(*state, "Syncing diagnostics and state to device");
        Logger::debug(*state, "PhotolysisProcess::run complete");
    }

    void PhotolysisProcess::finalize() {
        musica::Error err;
        if (tuvx_instance_) {
            musica::DeleteTuvx(tuvx_instance_, &err);
            tuvx_instance_ = nullptr;
        }
        if (grids_) {
            musica::DeleteGridMap(grids_, &err);
            grids_ = nullptr;
        }
        if (profiles_) {
            musica::DeleteProfileMap(profiles_, &err);
            profiles_ = nullptr;
        }
        if (radiators_) {
            musica::DeleteRadiatorMap(radiators_, &err);
            radiators_ = nullptr;
        }
    }

    void PhotolysisProcess::register_profile_if_missing(const StateManager& state,
                                                        const std::unordered_set<std::string>& config_defined_profiles,
                                                        const char* name, const char* units, musica::Grid* grid,
                                                        double default_val, std::size_t num_vals, musica::Error* err) {
        if (config_defined_profiles.find(name) == config_defined_profiles.end()) {
            Logger::debug(state, "Pre-registering missing profile", {{"name", name}, {"units", units}});
            musica::Profile* new_prof = musica::CreateProfile(name, units, grid, err);
            std::vector<double> dummy(num_vals, default_val);
            musica::SetProfileMidpointValues(new_prof, dummy.data(), num_vals, err);
            musica::AddProfile(profiles_, new_prof, err);
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
