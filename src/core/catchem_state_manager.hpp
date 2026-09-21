#pragma once
#include "catchem_chem_state.hpp"
#include "catchem_config_manager.hpp"
#include "catchem_constants.hpp"
#include "catchem_interop_field.hpp"
#include "catchem_met_state.hpp"
#include "catchem_met_utilities.hpp"
#include "catchem_physical_validation.hpp"
#include "catchem_time_state.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace catchem {

    class DiagnosticManager;

    class StateManager {
    public:
        StateManager(int nc, int nl, int ns);

        int column_count() const noexcept { return n_cols; }
        int level_count() const noexcept { return n_levels; }
        int species_count() const noexcept { return n_species; }
        std::size_t current_import_generation() const noexcept { return import_generation; }

        MetState& meteorology() noexcept { return met; }
        const MetState& meteorology() const noexcept { return met; }
        ChemState& chemistry() noexcept { return chem; }
        const ChemState& chemistry() const noexcept { return chem; }
        TimeState& clock() noexcept { return time; }
        const TimeState& clock() const noexcept { return time; }

        std::shared_ptr<ConfigManager> config_manager() const noexcept { return config_mgr; }
        std::shared_ptr<DiagnosticManager> diagnostic_manager() const noexcept { return diag_mgr; }
        void attach_config_manager(std::shared_ptr<ConfigManager> manager) { config_mgr = std::move(manager); }
        void attach_diagnostic_manager(std::shared_ptr<DiagnosticManager> manager) { diag_mgr = std::move(manager); }

        const std::string& configuration_path() const noexcept { return config_file_path; }
        void set_configuration_path(std::string path) { config_file_path = std::move(path); }
        const std::string& runtime_trace_id() const noexcept { return trace_id; }
        void set_runtime_trace_id(std::string id) { trace_id = std::move(id); }

        PhysicalValidationPolicy validation_policy() const noexcept { return physical_validation_policy; }
        void set_validation_policy(PhysicalValidationPolicy policy) noexcept { physical_validation_policy = policy; }
        const PhysicalValidationReport& validation_report() const noexcept { return physical_validation_report; }

        static std::string canonical_field_name(std::string name) {
            auto key = canonicalize_field_identity(std::move(name));
            static const std::unordered_map<std::string, std::string> aliases{
                {"TEMPERATURE", "T"},
                {"TEMP", "T"},
                {"SPHUM", "QV"},
                {"RELATIVE_HUMIDITY", "RH"},
                {"PRESSURE_MID", "PMID"},
                {"PRESSURE_EDGE", "PEDGE"},
                {"AIR_DENSITY", "AIRDEN"},
                {"AIR_DENSITY_DRY", "AIRDEN_DRY"},
                {"BOX_HEIGHT", "BXHEIGHT"},
                {"DZ", "BXHEIGHT"},
                {"DELZ", "BXHEIGHT"},
                {"LATITUDE", "LAT"},
                {"LONGITUDE", "LON"},
                {"SURFACE_PRESSURE", "PS"},
                {"SKIN_TEMPERATURE", "TS"},
                {"SST", "TS"},
                {"FRICTION_VELOCITY", "USTAR"},
                {"PRESSURE_THICKNESS_OF_ATMOSPHERIC_LAYER", "DELP"},
                {"CLAY_FRACTION", "CLAYFRAC"},
                {"LAKE_FRACTION", "FRLAKE"},
                {"SNOW_FRACTION", "FRSNO"},
                {"VEGETATION_FRACTION", "GVF"},
                {"LEAF_AREA_INDEX", "LAI"},
                {"LAND_WATER_ICE_MASK", "LWI"},
                {"DRAG_COEFFICIENT", "CMM"},
                {"SAND_FRACTION", "SNDFRC"},
                {"SOIL_MOISTURE", "SOILM"},
                {"SURFACE_SOIL_MOISTURE", "GWETTOP"},
                {"U_10M", "U10M"},
                {"V_10M", "V10M"},
                {"THRESHOLD_FRICTION_VELOCITY", "USTAR_THRESHOLD"},
                {"ROUGHNESS_LENGTH", "Z0"},
                {"HPBL", "PBLH"},
                {"OL", "OBK"},
                {"GEOMETRIC_HEIGHT_EDGE", "Z"}};
            const auto found = aliases.find(key);
            return found == aliases.end() ? key : found->second;
        }

        template <int Rank> std::shared_ptr<InteropField<double, Rank>> find_field(const std::string& name) const {
            const auto key = canonical_field_name(name);
            if constexpr (Rank == 1) {
                auto found = fields_1d.find(key);
                return found == fields_1d.end() ? nullptr : found->second;
            } else if constexpr (Rank == 2) {
                auto found = fields_2d.find(key);
                if (found != fields_2d.end())
                    return found->second;
                return met.get_2d({key.c_str()});
            } else {
                auto found = fields_3d.find(key);
                if (found != fields_3d.end())
                    return found->second;
                return met.get_3d({key.c_str()});
            }
        }

        bool prepare_field_access(const FieldAccessContract& access) {
            const auto contract_matches = [&access](const auto& field) {
                return field && field->contract.units == access.units && field->contract.axes == access.axes &&
                       field->contract.persistence == access.persistence;
            };
            const int rank = static_cast<int>(access.axes.size());
            if (canonical_field_name(access.canonical_name) == "CONCENTRATION") {
                if (!contract_matches(chem.conc))
                    return false;
                if (access.execution_space == ExecutionSpaceIntent::Device)
                    chem.conc->sync_to_device();
                else
                    chem.conc->sync_to_host();
                return true;
            }
            if (rank == 1) {
                auto field = find_field<1>(access.canonical_name);
                if (!contract_matches(field) ||
                    (access.persistence == PersistencePolicy::Timestep && !field->is_current(import_generation)))
                    return false;
                if (access.execution_space == ExecutionSpaceIntent::Device)
                    field->sync_to_device();
                else
                    field->sync_to_host();
                return true;
            }
            if (rank == 2) {
                auto field = find_field<2>(access.canonical_name);
                if (!contract_matches(field) ||
                    (access.persistence == PersistencePolicy::Timestep && !field->is_current(import_generation)))
                    return false;
                if (access.execution_space == ExecutionSpaceIntent::Device)
                    field->sync_to_device();
                else
                    field->sync_to_host();
                return true;
            }
            if (rank == 3) {
                auto field = find_field<3>(access.canonical_name);
                if (!contract_matches(field) ||
                    (access.persistence == PersistencePolicy::Timestep && !field->is_current(import_generation)))
                    return false;
                if (access.execution_space == ExecutionSpaceIntent::Device)
                    field->sync_to_device();
                else
                    field->sync_to_host();
                return true;
            }
            return false;
        }

        void complete_field_access(const FieldAccessContract& access) {
            if (!access.writes())
                return;
            const bool device = access.execution_space == ExecutionSpaceIntent::Device;
            if (canonical_field_name(access.canonical_name) == "CONCENTRATION") {
                if (chem.conc)
                    device ? chem.conc->mark_device_modified() : chem.conc->mark_host_modified();
                return;
            }
            const int rank = static_cast<int>(access.axes.size());
            if (rank == 1) {
                auto field = find_field<1>(access.canonical_name);
                if (field)
                    device ? field->mark_device_modified() : field->mark_host_modified();
            } else if (rank == 2) {
                auto field = find_field<2>(access.canonical_name);
                if (field)
                    device ? field->mark_device_modified() : field->mark_host_modified();
            } else if (rank == 3) {
                auto field = find_field<3>(access.canonical_name);
                if (field)
                    device ? field->mark_device_modified() : field->mark_host_modified();
            }
        }

        template <int Rank>
        const double* read_field(const std::string& name, std::size_t required_generation = 0) const {
            auto field = find_field<Rank>(name);
            if (!field)
                return nullptr;
            if (required_generation != 0 && !field->is_current(required_generation))
                return nullptr;
            return field->host_read();
        }

        template <int Rank> double* write_field(const std::string& name, std::size_t required_generation = 0) {
            auto field = find_field<Rank>(name);
            if (!field)
                return nullptr;
            if (required_generation != 0 && !field->is_current(required_generation))
                return nullptr;
            return field->host_write();
        }

        static PersistencePolicy persistence_for(const std::string& name) {
            const auto key = canonical_field_name(name);
            return (key == "LAT" || key == "LON" || key == "AREA_M2") ? PersistencePolicy::Persistent
                                                                      : PersistencePolicy::Timestep;
        }

        static std::string units_for(const std::string& name) {
            const auto key = canonical_field_name(name);
            // Temperature
            if (key == "T" || key == "TS" || key == "T2M" || key == "SOILT")
                return "K";
            // Water vapor
            if (key == "QV" || key == "QV2M")
                return "kg/kg";
            // Unitless variables
            if (key == "RH" || key == "CLDFRC" || key == "SUNCOSMID")
                return "1";
            // Pressure related fields
            if (key == "PMID" || key == "PEDGE" || key == "PS" || key == "DELP")
                return "Pa";
            // Height related fields
            if (key == "Z" || key == "ZMID" || key == "BXHEIGHT" || key == "PBLH" || key == "ORO")
                return "m";
            // Density
            if (key == "AIRDEN" || key == "AIRDEN_DRY" || key == "MAIRDEN")
                return "kg/m3";
            // Coordinates
            if (key == "LAT" || key == "LON")
                return "degrees";
            // Cell area
            if (key == "AREA_M2")
                return "m2";
            // Wind fields
            if (key == "U" || key == "V" || key == "U10M" || key == "V10M" || key == "USTAR" ||
                key == "USTAR_THRESHOLD")
                return "m/s";
            // Precipitation
            if (key == "PRECLSC" || key == "PRECCON")
                return "m";
            if (key == "WCA")
                return "m3";
            // Rain / ice flux rates
            if (key == "PFILSAN" || key == "PFLLSAN")
                return "kg/m2/s";
            // Re-evaporation tendency
            if (key == "REEVAPLS")
                return "kg/kg/s";
            // Green Vegetative Fraction
            if (key == "GVF")
                return "frac";
            // Fractional ocean, seaice, snow
            if (key == "FROCEAN" || key == "FRSEAICE" || key == "FRSNO")
                return "frac";
            // Leaf Area index
            if (key == "LAI" || key == "FRLAI")
                return "m2/m2";
            // Land Parameters
            if (key == "CLAYFRAC" || key == "FRLAKE" || key == "LWI" || key == "DLUSE" || key == "DSOILTYPE" ||
                key == "SNDFRC" || key == "GWETTOP" || key == "CLDF" || key == "SSM" || key == "RDRAG" ||
                key == "FRLANDUSE" || key == "ILAND")
                return "1";
            // Sea Surface temperature
            if (key == "SST")
                return "K";
            // Canopy fields
            if (key == "CMM" || key == "RCA")
                return "s/m";
            // Soil Moisture
            if (key == "SOILM")
                return "m3/m3";
            // Roughness and Obhokov length
            if (key == "Z0" || key == "Z0H" || key == "OBK")
                return "m";
            // Heat fluxes and downward shortwave
            if (key == "HFLUX" || key == "EFLUX" || key == "HFLUX_UP" || key == "SWGDN")
                return "W/m2";
            return "";
        }

        static std::vector<std::string> aliases_for(const std::string& name) {
            const auto key = canonical_field_name(name);
            if (key == "T")
                return {"TEMPERATURE", "TEMP"};
            if (key == "QV")
                return {"SPHUM"};
            if (key == "RH")
                return {"RELATIVE_HUMIDITY"};
            if (key == "PMID")
                return {"PRESSURE_MID"};
            if (key == "PEDGE")
                return {"PRESSURE_EDGE"};
            if (key == "Z")
                return {"GEOMETRIC_HEIGHT_EDGE"};
            if (key == "AIRDEN")
                return {"AIR_DENSITY"};
            if (key == "AIRDEN_DRY")
                return {"AIR_DENSITY_DRY"};
            if (key == "BXHEIGHT")
                return {"DZ", "DELZ"};
            if (key == "LAT")
                return {"LATITUDE"};
            if (key == "LON")
                return {"LONGITUDE"};
            return {};
        }

        template <int Rank>
        static void populate_contract(const std::string& name, InteropField<double, Rank>& field,
                                      std::vector<SemanticAxis> axes) {
            field.contract.canonical_name = canonical_field_name(name);
            field.contract.units = units_for(name);
            field.contract.axes = std::move(axes);
            field.contract.aliases = aliases_for(name);
            field.contract.persistence = persistence_for(name);
        }

        void begin_import_generation() {
            ++import_generation;
            for (auto& [name, field] : met.fields_2d)
                if (field && field->contract.persistence == PersistencePolicy::Timestep)
                    field->invalidate();
            for (auto& [name, field] : met.fields_3d)
                if (field && field->contract.persistence == PersistencePolicy::Timestep)
                    field->invalidate();
        }

        void load_species_config(const std::string& filename) { chem.load_species_config(filename, config_mgr.get()); }

        void bind_met_field_2d(const std::string& name, double* ptr) {
            const auto canonical_name = canonical_field_name(name);
            auto field = met.get_2d({canonical_name.c_str()});
            if (field) {
                field->update_host_pointer(ptr);
                field->set_generation(import_generation);
            } else {
                auto new_field = std::make_shared<InteropField<double, 2>>(ptr, std::vector<int>{n_cols, 1});
                populate_contract(canonical_name, *new_field, {SemanticAxis::Column, SemanticAxis::Singleton});
                met.bind_2d_field(canonical_name, new_field);
                new_field->set_generation(import_generation);
            }
        }

        void bind_met_field_3d(const std::string& name, double* ptr) {
            std::string upper_name = canonical_field_name(name);
            const bool is_interface = upper_name == "PEDGE" || upper_name == "Z";
            bind_met_field_3d_contract(name, ptr, is_interface ? n_levels + 1 : n_levels,
                                       is_interface ? SemanticAxis::Interface : SemanticAxis::Level);
        }

        void bind_met_field_3d_contract(const std::string& name, double* ptr, int vertical_extent,
                                        SemanticAxis vertical_axis) {
            const auto canonical_name = canonical_field_name(name);
            auto field = met.get_3d({canonical_name.c_str()});
            if (field) {
                if (field->extent(1) != static_cast<std::size_t>(vertical_extent) || field->contract.axes.size() != 3 ||
                    field->contract.axes[1] != vertical_axis)
                    throw std::invalid_argument("Rebinding field with a different vertical contract: " + name);
                field->update_host_pointer(ptr);
                field->set_generation(import_generation);
            } else {
                auto new_field =
                    std::make_shared<InteropField<double, 3>>(ptr, std::vector<int>{n_cols, vertical_extent, 1});
                populate_contract(canonical_name, *new_field,
                                  {SemanticAxis::Column, vertical_axis, SemanticAxis::Singleton});
                met.bind_3d_field(canonical_name, new_field);
                new_field->set_generation(import_generation);
            }
        }

        double* find_3d_ptr(std::initializer_list<const char*> names) const {
            auto f = met.get_3d(names);
            if (f && f->availability == AvailabilityState::Current)
                return f->host_write();
            return nullptr;
        }

        double* find_2d_ptr(std::initializer_list<const char*> names) const {
            auto f = met.get_2d(names);
            if (f && f->availability == AvailabilityState::Current)
                return f->host_write();
            return nullptr;
        }

        void bind_unified_chemistry(double* ptr) {
            if (chem.conc) {
                chem.conc->update_host_pointer(ptr);
            } else {
                chem.conc =
                    std::make_shared<InteropField<double, 3>>(ptr, std::vector<int>{n_cols, n_levels, n_species});
                chem.conc->contract.canonical_name = "CONCENTRATION";
                chem.conc->contract.units = "mol/mol";
                chem.conc->contract.axes = {SemanticAxis::Column, SemanticAxis::Level, SemanticAxis::Species};
                chem.conc->contract.persistence = PersistencePolicy::Persistent;
            }
        }

        /**
         * @brief Derives layer thicknesses (BXHEIGHT / dz) hydrostatically.
         *
         * Computes the vertical distance of each grid cell based on pressure edges, temperature, and moisture content.
         */
        void derive_bxheight() {
            // Use a current host-provided layer thickness when available.
            // Hydrostatic reconstruction is only the fallback for hosts that
            // do not exchange BXHEIGHT.
            if (met.BXHEIGHT && met.BXHEIGHT->is_current(import_generation))
                return;
            if (!met.PEDGE || !met.T || met.PEDGE->availability != AvailabilityState::Current ||
                met.T->availability != AvailabilityState::Current)
                return;

            physical_validation_report.clear();
            met.PEDGE->sync_to_host();
            met.T->sync_to_host();
            if (met.QV)
                met.QV->sync_to_host();
            const double* physical_pedge = met.PEDGE->host_data();
            const double* physical_temperature = met.T->host_data();
            for (int icol = 0; icol < n_cols; ++icol)
                for (int ilev = 0; ilev < n_levels; ++ilev) {
                    const auto location = static_cast<std::size_t>(icol + ilev * n_cols);
                    const double lower = physical_pedge[icol + ilev * n_cols];
                    const double upper = physical_pedge[icol + (ilev + 1) * n_cols];
                    const double temperature = physical_temperature[location];
                    if (!std::isfinite(lower) || !std::isfinite(upper) || lower <= 0.0 || upper <= 0.0 ||
                        lower <= upper)
                        physical_validation_report.observe("PEDGE", "positive-descending", lower, location,
                                                           "invalid layers produce zero thickness");
                    if (!std::isfinite(temperature) || temperature <= 0.0)
                        physical_validation_report.observe("T", "finite-positive", temperature, location,
                                                           "temperature is clamped to 1 K");
                }
            if (physical_validation_policy == PhysicalValidationPolicy::Reject && !physical_validation_report.empty())
                throw std::domain_error("Physical validation failed before BXHEIGHT mutation:\n" +
                                        physical_validation_report.format());

            if (!met.BXHEIGHT) {
                auto buf = std::make_shared<std::vector<double>>(n_cols * n_levels, 0.0);
                owned_buffers.push_back(buf);
                bind_met_field_3d("BXHEIGHT", buf->data());
            }

            met.PEDGE->sync_to_device();
            met.T->sync_to_device();
            if (met.QV)
                met.QV->sync_to_device();
            int nc = n_cols;
            int nl = n_levels;

            auto pedge = met.PEDGE->view();
            auto temp = met.T->view();
            auto bxheight = met.BXHEIGHT->view();
            // Resolve QV outside the kernel: capturing `met` would capture `this`
            // (deprecated under C++20 with [=], and not device-safe). When QV is
            // absent, `qv` aliases `temp` purely to have a view of the right type;
            // it is never read in that case.
            const bool have_qv = met.QV && met.QV->availability == AvailabilityState::Current;
            auto qv = have_qv ? met.QV->view() : temp;

#ifdef CATCHEM_ENABLE_KOKKOS
            Kokkos::parallel_for(
                "derive_bxheight_kernel",
                Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, 0}, {nc, nl}),
                KOKKOS_LAMBDA(int icol, int ilev) {
                    double p_lower = pedge(icol, ilev, 0);
                    double p_upper = pedge(icol, ilev + 1, 0);

                    if (p_upper > 0.0 && p_lower > 0.0 && p_lower > p_upper) {
                        double q_val = have_qv ? qv(icol, ilev, 0) : 0.0;
                        if (!(q_val >= 0.0 && q_val < 1.0))
                            q_val = q_val == q_val ? (q_val < 0.0 ? 0.0 : 0.9999) : 0.0;
                        double t_val = temp(icol, ilev, 0);
                        if (!(t_val > 0.0))
                            t_val = 1.0;
                        bxheight(icol, ilev, 0) =
                            met_utilities::hydrostatic_layer_thickness(p_lower, p_upper, t_val, q_val);
                    } else {
                        bxheight(icol, ilev, 0) = 0.0;
                    }
                });
#else
            for (int icol = 0; icol < nc; ++icol) {
                for (int ilev = 0; ilev < nl; ++ilev) {
                    double p_lower = pedge(icol, ilev, 0);
                    double p_upper = pedge(icol, ilev + 1, 0);

                    if (p_upper > 0.0 && p_lower > 0.0 && p_lower > p_upper) {
                        double q_val = have_qv ? qv(icol, ilev, 0) : 0.0;
                        if (!std::isfinite(q_val))
                            q_val = 0.0;
                        q_val = std::clamp(q_val, 0.0, 0.9999);
                        double t_val = temp(icol, ilev, 0);
                        if (!std::isfinite(t_val) || t_val <= 0.0)
                            t_val = 1.0;
                        bxheight(icol, ilev, 0) =
                            met_utilities::hydrostatic_layer_thickness(p_lower, p_upper, t_val, q_val);
                    } else {
                        bxheight(icol, ilev, 0) = 0.0;
                    }
                }
            }
#endif
            met.BXHEIGHT->mark_device_modified();
            met.BXHEIGHT->set_generation(import_generation);
        }

        /**
         * @brief Derives dry air density (AIRDEN_DRY) using the Ideal Gas Law.
         *
         * Calculates dry air mass densities based on mid-point pressures, temperature, and specific humidity.
         */
        void derive_airden_dry() {
            // A current host field is authoritative.  Derive only when the
            // host did not supply AIRDEN_DRY for this import generation.
            if (met.AIRDEN_DRY && met.AIRDEN_DRY->is_current(import_generation))
                return;
            if (!met.PMID || !met.T || met.PMID->availability != AvailabilityState::Current ||
                met.T->availability != AvailabilityState::Current)
                return;

            physical_validation_report.clear();
            met.PMID->sync_to_host();
            met.T->sync_to_host();
            if (met.QV)
                met.QV->sync_to_host();
            const double* physical_pmid = met.PMID->host_data();
            const double* physical_temperature = met.T->host_data();
            const double* physical_humidity = met.QV ? met.QV->host_data() : nullptr;
            for (int icol = 0; icol < n_cols; ++icol)
                for (int ilev = 0; ilev < n_levels; ++ilev) {
                    const auto location = static_cast<std::size_t>(icol + ilev * n_cols);
                    const double pressure = physical_pmid[location];
                    const double temperature = physical_temperature[location];
                    const double humidity = physical_humidity ? physical_humidity[location] : 0.0;
                    if (!std::isfinite(pressure) || pressure <= 0.0)
                        physical_validation_report.observe("PMID", "finite-positive", pressure, location);
                    if (!std::isfinite(temperature) || temperature <= 0.0)
                        physical_validation_report.observe("T", "finite-positive", temperature, location,
                                                           "temperature is clamped to 1 K");
                    if (!std::isfinite(humidity) || humidity < 0.0 || humidity >= 1.0)
                        physical_validation_report.observe("QV", "finite-[0,1)", humidity, location,
                                                           "humidity is clamped to [0,0.9999]");
                }
            if (physical_validation_policy == PhysicalValidationPolicy::Reject && !physical_validation_report.empty())
                throw std::domain_error("Physical validation failed before AIRDEN_DRY mutation:\n" +
                                        physical_validation_report.format());

            if (!met.AIRDEN_DRY) {
                auto buf = std::make_shared<std::vector<double>>(n_cols * n_levels, 0.0);
                owned_buffers.push_back(buf);
                bind_met_field_3d("AIRDEN_DRY", buf->data());
            }

            met.PMID->sync_to_device();
            met.T->sync_to_device();
            if (met.QV)
                met.QV->sync_to_device();
            int nc = n_cols;
            int nl = n_levels;

            auto pmid = met.PMID->view();
            auto temp = met.T->view();
            auto airden_dry = met.AIRDEN_DRY->view();
            // See derive_bxheight: resolve QV outside the kernel.
            const bool have_qv = met.QV && met.QV->availability == AvailabilityState::Current;
            auto qv = have_qv ? met.QV->view() : temp;

#ifdef CATCHEM_ENABLE_KOKKOS
            Kokkos::parallel_for(
                "derive_airden_dry_kernel",
                Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, 0}, {nc, nl}),
                KOKKOS_LAMBDA(int icol, int ilev) {
                    double q = have_qv ? qv(icol, ilev, 0) : 0.0;
                    if (!(q >= 0.0 && q < 1.0))
                        q = q == q ? (q < 0.0 ? 0.0 : 0.9999) : 0.0;
                    double pressure = pmid(icol, ilev, 0);
                    if (!(pressure > 0.0))
                        pressure = 1.0;
                    double t_val = temp(icol, ilev, 0);
                    if (!(t_val > 0.0))
                        t_val = 1.0;
                    airden_dry(icol, ilev, 0) = met_utilities::dry_air_density(pressure, t_val, q);
                });
#else
            for (int icol = 0; icol < nc; ++icol) {
                for (int ilev = 0; ilev < nl; ++ilev) {
                    double q = have_qv ? qv(icol, ilev, 0) : 0.0;
                    if (!std::isfinite(q))
                        q = 0.0;
                    q = std::clamp(q, 0.0, 0.9999);
                    double pressure = pmid(icol, ilev, 0);
                    if (!std::isfinite(pressure) || pressure <= 0.0)
                        pressure = 1.0;
                    double t_val = temp(icol, ilev, 0);
                    if (!std::isfinite(t_val) || t_val <= 0.0)
                        t_val = 1.0;
                    airden_dry(icol, ilev, 0) = met_utilities::dry_air_density(pressure, t_val, q);
                }
            }
#endif
            met.AIRDEN_DRY->mark_device_modified();
            met.AIRDEN_DRY->set_generation(import_generation);
        }

        // Maintain the legacy metstate definition of AIRDEN: pressure divided
        // by dry-air gas constant and temperature.  AIRDEN_DRY is the
        // humidity-corrected quantity and is derived separately above.
        void derive_airden() {
            if (met.AIRDEN && met.AIRDEN->is_current(import_generation))
                return;
            if (!met.PMID || !met.T || !met.PMID->is_current(import_generation) ||
                !met.T->is_current(import_generation))
                throw std::runtime_error("Cannot derive AIRDEN: requires current PMID and T");
            auto airden = find_field<3>("AIRDEN");
            if (!airden) {
                auto buffer = std::make_shared<std::vector<double>>(static_cast<std::size_t>(n_cols) * n_levels, 0.0);
                owned_buffers.push_back(buffer);
                bind_met_field_3d("AIRDEN", buffer->data());
                airden = find_field<3>("AIRDEN");
            }
            met.PMID->sync_to_host();
            met.T->sync_to_host();
            const double* p = met.PMID->host_data();
            const double* t = met.T->host_data();
            double* output = airden->host_write();
            for (int i = 0; i < n_cols * n_levels; ++i)
                output[i] = p[i] / (constants::RD * t[i]);
            airden->mark_host_modified();
            airden->set_generation(import_generation);
        }

        // Derive DELP from the host pressure interfaces.  The sign is not
        // hidden with abs(): a non-descending host vertical coordinate is an
        // invalid layer and produces zero, making the contract error visible.
        void derive_delp() {
            // An imported DELP is authoritative.  Only derive it from PEDGE
            // when the caller did not provide a current layer thickness.
            if (const auto delp = find_field<3>("DELP"); delp && delp->is_current(import_generation))
                return;
            if (!met.PEDGE || !met.PEDGE->is_current(import_generation))
                throw std::runtime_error("Cannot derive DELP: PEDGE is not current");
            auto delp = find_field<3>("DELP");
            if (!delp) {
                auto buffer = std::make_shared<std::vector<double>>(static_cast<std::size_t>(n_cols) * n_levels, 0.0);
                owned_buffers.push_back(buffer);
                bind_met_field_3d("DELP", buffer->data());
                delp = find_field<3>("DELP");
            }
            met.PEDGE->sync_to_host();
            const double* pedge = met.PEDGE->host_data();
            double* output = delp->host_write();
            for (int level = 0; level < n_levels; ++level)
                for (int column = 0; column < n_cols; ++column) {
                    const std::size_t index = static_cast<std::size_t>(column + level * n_cols);
                    output[index] = met_utilities::pressure_thickness(pedge[index], pedge[index + n_cols]);
                }
            delp->mark_host_modified();
            delp->set_generation(import_generation);
        }

        void derive_obk() {
            if (const auto obk = find_field<2>("OBK"); obk && obk->is_current(import_generation))
                return;
            if (!met.USTAR || !met.TS || !met.HFLUX || !met.PMID || !met.T ||
                !met.USTAR->is_current(import_generation) || !met.TS->is_current(import_generation) ||
                !met.HFLUX->is_current(import_generation) || !met.PMID->is_current(import_generation) ||
                !met.T->is_current(import_generation))
                throw std::runtime_error("Cannot derive OBK: requires current USTAR, TS, HFLUX, PMID, and T");
            auto obk = find_field<2>("OBK");
            if (!obk) {
                auto buffer = std::make_shared<std::vector<double>>(n_cols, 0.0);
                owned_buffers.push_back(buffer);
                bind_met_field_2d("OBK", buffer->data());
                obk = find_field<2>("OBK");
            }
            met.USTAR->sync_to_host();
            met.TS->sync_to_host();
            met.HFLUX->sync_to_host();
            met.PMID->sync_to_host();
            met.T->sync_to_host();
            const double* ustar = met.USTAR->host_data();
            const double* ts = met.TS->host_data();
            const double* hflux = met.HFLUX->host_data();
            const double* p = met.PMID->host_data();
            const double* t = met.T->host_data();
            double* output = obk->host_write();
            for (int c = 0; c < n_cols; ++c)
                output[c] =
                    met_utilities::monin_obukhov_length(ustar[c], ts[c], hflux[c], p[c] / (constants::RD * t[c]));
            obk->mark_host_modified();
            obk->set_generation(import_generation);
        }

        // Legacy metstate_mod defaults missing ocean salinity to zero. This
        // disables the optional salinity-dependent ocean chemistry treatment
        // while retaining a host-provided salinity field when one is available.
        void derive_salinity() {
            if (const auto salinity = find_field<2>("SALINITY"); salinity && salinity->is_current(import_generation))
                return;
            auto buffer = std::make_shared<std::vector<double>>(n_cols, 0.0);
            owned_buffers.push_back(buffer);
            bind_met_field_2d("SALINITY", buffer->data());
            auto salinity = find_field<2>("SALINITY");
            salinity->mark_host_modified();
            salinity->set_generation(import_generation);
        }

        // Reconstruct the legacy 20-category land-use fields from UFS's
        // dominant land-use index and leaf-area index. DLUSE=0 is water,
        // which legacy CATChem places in category 17 (one-based).  The derived
        // FRLAI/FRLANDUSE/ILAND fields are generic surface descriptors; any
        // process that declares them in its contract may consume them.
        void derive_landuse_categories() {
            constexpr int n_landuse = 20;
            auto dluse = find_field<2>("DLUSE");
            if (!dluse || !dluse->is_current(import_generation))
                throw std::runtime_error("Cannot derive dry-deposition land use: requires current DLUSE");

            if (const auto frlai = find_field<3>("FRLAI"); !frlai || !frlai->is_current(import_generation)) {
                auto lai = find_field<2>("LAI");
                if (!lai || !lai->is_current(import_generation))
                    throw std::runtime_error("Cannot derive FRLAI: requires current LAI and DLUSE");
                auto buffer = std::make_shared<std::vector<double>>(static_cast<std::size_t>(n_cols) * n_landuse, 0.0);
                owned_buffers.push_back(buffer);
                dluse->sync_to_host();
                lai->sync_to_host();
                const double* land_type = dluse->host_data();
                const double* leaf_area = lai->host_data();
                for (int column = 0; column < n_cols; ++column) {
                    const int category = static_cast<int>(land_type[column]);
                    if (category >= 1 && category <= n_landuse)
                        buffer->at(static_cast<std::size_t>(column + (category - 1) * n_cols)) = leaf_area[column];
                    for (int slot = 14; slot <= 16; ++slot)
                        buffer->at(static_cast<std::size_t>(column + slot * n_cols)) = 0.0;
                }
                bind_met_field_3d_contract("FRLAI", buffer->data(), n_landuse, SemanticAxis::SoilLayer);
                auto out = find_field<3>("FRLAI");
                out->mark_host_modified();
                out->set_generation(import_generation);
            }

            if (const auto frlanduse = find_field<3>("FRLANDUSE");
                !frlanduse || !frlanduse->is_current(import_generation)) {
                auto buffer = std::make_shared<std::vector<double>>(static_cast<std::size_t>(n_cols) * n_landuse, 0.0);
                owned_buffers.push_back(buffer);
                dluse->sync_to_host();
                const double* land_type = dluse->host_data();
                for (int column = 0; column < n_cols; ++column) {
                    int category = static_cast<int>(land_type[column]);
                    if (category == 0)
                        category = 17;
                    if (category >= 1 && category <= n_landuse)
                        buffer->at(static_cast<std::size_t>(column + (category - 1) * n_cols)) = 1.0;
                }
                bind_met_field_3d_contract("FRLANDUSE", buffer->data(), n_landuse, SemanticAxis::SoilLayer);
                auto out = find_field<3>("FRLANDUSE");
                out->mark_host_modified();
                out->set_generation(import_generation);
            }

            if (const auto iland = find_field<3>("ILAND"); !iland || !iland->is_current(import_generation)) {
                auto buffer = std::make_shared<std::vector<double>>(static_cast<std::size_t>(n_cols) * n_landuse, 0.0);
                owned_buffers.push_back(buffer);
                for (int category = 0; category < n_landuse; ++category)
                    for (int column = 0; column < n_cols; ++column)
                        buffer->at(static_cast<std::size_t>(column + category * n_cols)) = category + 1;
                bind_met_field_3d_contract("ILAND", buffer->data(), n_landuse, SemanticAxis::SoilLayer);
                auto out = find_field<3>("ILAND");
                out->mark_host_modified();
                out->set_generation(import_generation);
            }
        }

        void derive_relative_humidity() {
            if (const auto rh = find_field<3>("RH"); rh && rh->is_current(import_generation))
                return;
            if (!met.T || !met.QV || !met.PMID)
                throw std::runtime_error("Cannot derive RH: requires T, QV, and PMID");
            auto buffer = std::make_shared<std::vector<double>>(static_cast<std::size_t>(n_cols) * n_levels, 0.0);
            owned_buffers.push_back(buffer);
            bind_met_field_3d("RH", buffer->data());
            met.T->sync_to_host();
            met.QV->sync_to_host();
            met.PMID->sync_to_host();
            const double* t = met.T->host_data();
            const double* qv = met.QV->host_data();
            const double* p = met.PMID->host_data();
            for (int i = 0; i < n_cols * n_levels; ++i)
                buffer->at(i) = met_utilities::relative_humidity(t[i], qv[i], p[i]);
            met.RH->mark_host_modified();
            met.RH->set_generation(import_generation);
        }

        void derive_column_cloud_fraction() {
            if (const auto out = find_field<2>("CLDFRC"); out && out->is_current(import_generation))
                return;
            auto cldf = find_field<3>("CLDF");
            if (!cldf)
                throw std::runtime_error("Cannot derive CLDFRC: requires CLDF");
            auto buffer = std::make_shared<std::vector<double>>(n_cols, 0.0);
            owned_buffers.push_back(buffer);
            bind_met_field_2d("CLDFRC", buffer->data());
            cldf->sync_to_host();
            const double* src = cldf->host_data();
            // Legacy parity: the upstream metstate_mod derives CLDFRC as the
            // column-total cloud fraction (`CLDFRC(:,:) = SUM(CLDF, DIM=3)`),
            // not the surface layer.  The sum may exceed one where clouds
            // overlap multiple layers; upstream applies no clamp, so neither
            // do we.
            for (int c = 0; c < n_cols; ++c) {
                double total = 0.0;
                for (int lev = 0; lev < n_levels; ++lev)
                    total += src[static_cast<std::size_t>(c) + static_cast<std::size_t>(lev) * n_cols];
                buffer->at(c) = total;
            }
            auto derived = find_field<2>("CLDFRC");
            derived->mark_host_modified();
            derived->set_generation(import_generation);
        }

        void derive_suncosmid() {
            if (const auto out = find_field<2>("SUNCOSMID"); out && out->is_current(import_generation))
                return;
            if (!met.LAT || !met.LON)
                throw std::runtime_error("Cannot derive SUNCOSMID: requires LAT and LON");
            auto buffer = std::make_shared<std::vector<double>>(n_cols, 0.0);
            owned_buffers.push_back(buffer);
            bind_met_field_2d("SUNCOSMID", buffer->data());
            met.LAT->sync_to_host();
            met.LON->sync_to_host();
            const double* lat = met.LAT->host_data();
            const double* lon = met.LON->host_data();
            for (int c = 0; c < n_cols; ++c)
                buffer->at(c) = time.get_cos_sza(lat[c], lon[c], true);
            auto derived = find_field<2>("SUNCOSMID");
            derived->mark_host_modified();
            derived->set_generation(import_generation);
        }

        // Derive large-scale/anvil precipitation re-evaporation [kg/kg/s].
        // This is a process-owned physical diagnostic, not a missing-input
        // fallback.  Its source fields are validated before any mutation.
        void derive_reevapls() {
            // Preserve a current host-provided re-evaporation tendency.
            // The diagnostic calculation below is the fallback for hosts
            // that do not exchange REEVAPLS.
            if (const auto reevapls = find_field<3>("REEVAPLS"); reevapls && reevapls->is_current(import_generation))
                return;
            auto pfilsan = find_field<3>("PFILSAN");
            auto pfllsan = find_field<3>("PFLLSAN");
            if (!met.T || !met.QV || !met.PMID || !met.PEDGE || !pfilsan || !pfllsan ||
                !met.T->is_current(import_generation) || !met.QV->is_current(import_generation) ||
                !met.PMID->is_current(import_generation) || !met.PEDGE->is_current(import_generation) ||
                !pfilsan->is_current(import_generation) || !pfllsan->is_current(import_generation))
                throw std::runtime_error(
                    "WetDep cannot derive REEVAPLS: requires current T, QV, PMID, PEDGE, PFILSAN, and PFLLSAN");

            auto reevapls = find_field<3>("REEVAPLS");
            if (!reevapls) {
                auto buffer = std::make_shared<std::vector<double>>(static_cast<std::size_t>(n_cols) * n_levels, 0.0);
                owned_buffers.push_back(buffer);
                bind_met_field_3d("REEVAPLS", buffer->data());
                reevapls = find_field<3>("REEVAPLS");
            }
            met.T->sync_to_host();
            met.QV->sync_to_host();
            met.PMID->sync_to_host();
            met.PEDGE->sync_to_host();
            pfilsan->sync_to_host();
            pfllsan->sync_to_host();
            double* output = reevapls->host_write();
            const double* temperature = met.T->host_data();
            const double* humidity = met.QV->host_data();
            const double* pmid = met.PMID->host_data();
            const double* pedge = met.PEDGE->host_data();
            const double* ice_flux = pfilsan->host_data();
            const double* liquid_flux = pfllsan->host_data();
            for (int level = 0; level < n_levels; ++level)
                for (int column = 0; column < n_cols; ++column) {
                    const std::size_t index = static_cast<std::size_t>(column + level * n_cols);
                    output[index] = met_utilities::large_scale_reevaporation(
                        temperature[index], humidity[index], pmid[index], pedge[index], pedge[index + n_cols],
                        // Match legacy metstate_mod: the layer-k re-evaporation
                        // diagnostic uses the flux at interface k.  The
                        // nlev+1 extent is retained for consumers that need
                        // the adjacent interface (for example flux divergence).
                        ice_flux[index], liquid_flux[index]);
                }
            reevapls->mark_host_modified();
            reevapls->set_generation(import_generation);
        }

        void sync_to_device() {
            for (auto& [k, v] : fields_1d)
                v->sync_to_device();
            for (auto& [k, v] : fields_2d)
                v->sync_to_device();
            for (auto& [k, v] : fields_3d)
                v->sync_to_device();
            for (auto& [k, v] : met.fields_2d)
                v->sync_to_device();
            for (auto& [k, v] : met.fields_3d)
                v->sync_to_device();
            if (chem.conc)
                chem.conc->sync_to_device();
        }

        void validate_ready_for_execution() const {
            if (chem.conc && (chem.conc->extent(0) != static_cast<std::size_t>(n_cols) ||
                              chem.conc->extent(1) != static_cast<std::size_t>(n_levels) ||
                              chem.conc->extent(2) != static_cast<std::size_t>(n_species)))
                throw std::runtime_error("Chemistry concentration contract does not match the active grid/mechanism");
        }

        void sync_to_host() {
            for (auto& [k, v] : fields_1d)
                v->sync_to_host();
            for (auto& [k, v] : fields_2d)
                v->sync_to_host();
            for (auto& [k, v] : fields_3d)
                v->sync_to_host();
            for (auto& [k, v] : met.fields_2d)
                v->sync_to_host();
            for (auto& [k, v] : met.fields_3d)
                v->sync_to_host();
            if (chem.conc)
                chem.conc->sync_to_host();
        }

        void bind_field_1d(const std::string& name, double* ptr) {
            auto it = fields_1d.find(name);
            if (it != fields_1d.end() && it->second) {
                it->second->update_host_pointer(ptr);
            } else {
                fields_1d[name] = std::make_shared<InteropField<double, 1>>(ptr, std::vector<int>{n_cols});
            }
        }

        void bind_field_2d(const std::string& name, double* ptr) {
            auto it = fields_2d.find(name);
            if (it != fields_2d.end() && it->second) {
                it->second->update_host_pointer(ptr);
            } else {
                fields_2d[name] = std::make_shared<InteropField<double, 2>>(ptr, std::vector<int>{n_cols, 1});
            }
        }

        void bind_field_3d(const std::string& name, double* ptr) {
            auto it = fields_3d.find(name);
            if (it != fields_3d.end() && it->second) {
                it->second->update_host_pointer(ptr);
            } else {
                fields_3d[name] =
                    std::make_shared<InteropField<double, 3>>(ptr, std::vector<int>{n_cols, n_levels, n_species});
            }
        }

        double* get_host_pointer_1d(const std::string& name) {
            if (fields_1d.find(name) == fields_1d.end())
                return nullptr;
            return fields_1d.at(name)->host_write();
        }

        const double* get_host_read_pointer_1d(const std::string& name) {
            const auto found = fields_1d.find(name);
            return found == fields_1d.end() ? nullptr : found->second->host_read();
        }

        double* get_host_pointer_2d(const std::string& name) {
            if (fields_2d.find(name) != fields_2d.end())
                return fields_2d.at(name)->host_write();
            if (met.fields_2d.find(name) != met.fields_2d.end())
                return met.fields_2d.at(name)->host_write();
            return nullptr;
        }

        const double* get_host_read_pointer_2d(const std::string& name) {
            if (fields_2d.find(name) != fields_2d.end())
                return fields_2d.at(name)->host_read();
            if (met.fields_2d.find(name) != met.fields_2d.end())
                return met.fields_2d.at(name)->host_read();
            return nullptr;
        }

        double* get_host_pointer_3d(const std::string& name) {
            if (fields_3d.find(name) != fields_3d.end())
                return fields_3d.at(name)->host_write();
            if (met.fields_3d.find(name) != met.fields_3d.end())
                return met.fields_3d.at(name)->host_write();
            return nullptr;
        }

        const double* get_host_read_pointer_3d(const std::string& name) {
            if (fields_3d.find(name) != fields_3d.end())
                return fields_3d.at(name)->host_read();
            if (met.fields_3d.find(name) != met.fields_3d.end())
                return met.fields_3d.at(name)->host_read();
            return nullptr;
        }

        std::pair<std::size_t, std::size_t> transfer_statistics() const {
            std::size_t transfers = 0;
            std::size_t bytes = 0;
            const auto collect = [&transfers, &bytes](const auto& fields) {
                for (const auto& [name, field] : fields) {
                    (void)name;
                    const auto count = field->host_to_device_sync_count + field->device_to_host_sync_count;
                    std::size_t elements = 1;
                    for (const auto extent : field->immutable_extents)
                        elements *= extent;
                    transfers += count;
                    bytes += count * elements * sizeof(double);
                }
            };
            collect(fields_1d);
            collect(fields_2d);
            collect(fields_3d);
            collect(met.fields_2d);
            collect(met.fields_3d);
            if (chem.conc) {
                const auto count = chem.conc->host_to_device_sync_count + chem.conc->device_to_host_sync_count;
                std::size_t elements = 1;
                for (const auto extent : chem.conc->immutable_extents)
                    elements *= extent;
                transfers += count;
                bytes += count * elements * sizeof(double);
            }
            return {transfers, bytes};
        }

    private:
        int n_cols;
        int n_levels;
        int n_species;
        std::size_t import_generation = 0;
        PhysicalValidationPolicy physical_validation_policy = PhysicalValidationPolicy::WarnAndClamp;
        PhysicalValidationReport physical_validation_report;

        std::string config_file_path;
        std::string trace_id;
        std::shared_ptr<ConfigManager> config_mgr;
        std::shared_ptr<DiagnosticManager> diag_mgr;

        std::unordered_map<std::string, std::shared_ptr<InteropField<double, 1>>> fields_1d;
        std::unordered_map<std::string, std::shared_ptr<InteropField<double, 2>>> fields_2d;
        std::unordered_map<std::string, std::shared_ptr<InteropField<double, 3>>> fields_3d;
        std::vector<std::shared_ptr<std::vector<double>>> owned_buffers;

        MetState met;
        ChemState chem;
        TimeState time;
    };

} // namespace catchem
