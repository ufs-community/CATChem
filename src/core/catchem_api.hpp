/**
 * @file catchem_api.hpp
 * @brief Flat, BIND(C) linkable API endpoints for CATChem host model integration.
 */

#pragma once

#include "catchem_precision.hpp"

#ifdef __cplusplus
extern "C" {
#endif

enum catchem_dataflow_status {
    CATCHEM_SUCCESS = 0,
    CATCHEM_NULL_ARGUMENT = 1,
    CATCHEM_MISSING_FIELD = 2,
    CATCHEM_RANK_MISMATCH = 3,
    CATCHEM_EXTENT_MISMATCH = 4,
    CATCHEM_INVALID_INDEX = 5,
    CATCHEM_STALE_GENERATION = 6,
    CATCHEM_DUPLICATE_MAPPING = 7,
    CATCHEM_INVALID_STATE = 8,
    CATCHEM_INTERNAL_ERROR = 9,
    CATCHEM_INVALID_HANDLE = 10,
    CATCHEM_WRONG_HANDLE_TYPE = 11,
    CATCHEM_RUNTIME_UNAVAILABLE = 12,
    CATCHEM_CONTRACT_VIOLATION = 13,
    CATCHEM_INVALID_CONFIGURATION = 14,
    CATCHEM_PROCESS_FAILURE = 15,
    CATCHEM_SHUTDOWN_FAILURE = 16,
    CATCHEM_PHYSICAL_VALIDATION_FAILURE = 17
};

int catchem_core_create_checked(int nc, int nl, int ns, void** core_out);
int catchem_core_create_from_config_checked(const char* config_file, void** core_out);
int catchem_core_create_from_config_with_grid_checked(const char* config_file, int ncols, int nlevels, void** core_out);
int catchem_core_destroy_checked(void* core_ptr);
int catchem_core_get_state_manager_checked(void* core_ptr, void** state_out);
int catchem_get_last_error(char* buffer, int max_len);

/**
 * @brief Creates the C++ Core orchestrator instance.
 * @param nc Number of horizontal columns.
 * @param nl Number of vertical levels.
 * @param ns Number of chemical species.
 * @return Opacity-wrapped void* handle pointing to catchem::Core.
 */
void* catchem_core_create(int nc, int nl, int ns);

/**
 * @brief Creates the C++ Core orchestrator from a YAML config file.
 * @param config_file Null-terminated filesystem path.
 * @return Opacity-wrapped void* handle pointing to catchem::Core.
 */
void* catchem_core_create_from_config(const char* config_file);

/**
 * @brief Creates the Core from a YAML config file with host-supplied grid dimensions.
 *
 * Configuration comes from the file; the grid is sized by the host (required
 * under domain decomposition, e.g. UFS per-rank tiles).
 * @param config_file Null-terminated filesystem path.
 * @param ncols Host-local number of columns (nx*ny).
 * @param nlevels Number of vertical levels.
 * @return Opacity-wrapped void* handle pointing to catchem::Core, or NULL on failure.
 */
void* catchem_core_create_from_config_with_grid(const char* config_file, int ncols, int nlevels);

/**
 * @brief Destroys the Core orchestrator instance and releases heap memory.
 * @param core_ptr Pointer to the active catchem::Core instance.
 */
void catchem_core_destroy(void* core_ptr);

/**
 * @brief Extracts the underlying StateManager handle from Core.
 * @param core_ptr Pointer to the active catchem::Core instance.
 * @return Pointer to catchem::StateManager.
 */
void* catchem_core_get_state_manager(void* core_ptr);

/** @brief Bind a 1D double field to the StateManager registry. */
void catchem_state_bind_1d(void* state_ptr, const char* name, double* ptr);
int catchem_state_bind_1d_checked(void* state_ptr, const char* name, double* ptr, int dim1);

/** @brief Bind a 2D double field to the StateManager registry. */
void catchem_state_bind_2d(void* state_ptr, const char* name, double* ptr);
int catchem_state_bind_2d_checked(void* state_ptr, const char* name, double* ptr, int dim1, int dim2);

/** @brief Bind a 3D double field to the StateManager registry. */
void catchem_state_bind_3d(void* state_ptr, const char* name, double* ptr);
int catchem_state_bind_3d_checked(void* state_ptr, const char* name, double* ptr, int dim1, int dim2, int dim3);

/** @brief Bind a 2D meteorological field by name. */
void catchem_state_bind_met_2d(void* state_ptr, const char* name, double* ptr);

/** @brief Bind a 3D meteorological field by name. */
void catchem_state_bind_met_3d(void* state_ptr, const char* name, double* ptr);
int catchem_state_begin_import_generation(void* state_ptr);
int catchem_state_set_physical_validation_policy_checked(void* state_ptr, int policy);
int catchem_state_get_physical_validation_report_checked(void* state_ptr, int* issue_count, char* detail,
                                                         int detail_length);
int catchem_state_bind_met_2d_checked(void* state_ptr, const char* name, double* ptr, int dim1, int dim2);
int catchem_state_bind_met_3d_checked(void* state_ptr, const char* name, double* ptr, int dim1, int dim2, int dim3);
int catchem_state_bind_met_3d_axis_checked(void* state_ptr, const char* name, double* ptr, int dim1, int dim2, int dim3,
                                           int semantic_axis);

/** @brief Binds the contiguous, multi-species unified chemistry concentrations array. */
void catchem_state_bind_unified_chemistry(void* state_ptr, double* ptr);
int catchem_state_bind_unified_chemistry_checked(void* state_ptr, double* ptr, int dim1, int dim2, int dim3);
int catchem_state_mark_chem_host_modified(void* state_ptr);

/** @brief Sets current simulation time within the state. */
void catchem_state_set_time(void* state_ptr, int yr, int mo, int dy, int hr, int mn, int sc, int doy, double tstep);
int catchem_state_set_time_checked(void* state_ptr, int yr, int mo, int dy, int hr, int mn, int sc, int doy,
                                   double tstep);

/** @brief Synchronizes registered host pointers to Kokkos device memory space. */
void catchem_state_sync_to_device(void* state_ptr);
int catchem_state_sync_to_device_checked(void* state_ptr);

/** @brief Synchronizes Kokkos device calculations back to host buffers. */
void catchem_state_sync_to_host(void* state_ptr);
int catchem_state_sync_to_host_checked(void* state_ptr);

/** @brief Retrieves direct host pointers from the 1D, 2D, or 3D fields. */
double* catchem_state_get_pointer_1d(void* state_ptr, const char* name);
double* catchem_state_get_pointer_2d(void* state_ptr, const char* name);
double* catchem_state_get_pointer_3d(void* state_ptr, const char* name);
int catchem_state_get_pointer_3d_checked(void* state_ptr, const char* name, void** ptr_out);
double* catchem_state_get_species_conc_pointer(void* state_ptr, int species_index);
int catchem_state_get_species_conc_pointer_checked(void* state_ptr, int species_index, int dim1, int dim2,
                                                   double** ptr_out);

/**
 * @brief Executes a single timestepped execution over scheduled processes.
 * @param core_ptr Core pointer.
 * @param dt Step size in seconds.
 */
int catchem_core_run_timestep(void* core_ptr, double dt);
int catchem_core_get_timestep_outcome(void* core_ptr, int* status, long long* timestep, double* duration,
                                      long long* import_generation, int* process_index, int* state_classification,
                                      char* process_name, int process_name_len, char* cause, int cause_len);

/** @brief Registers and attaches an active physics process handler. */
void catchem_core_add_process_by_name(void* core_ptr, const char* name);

/** @brief Returns the number of active physics processes scheduled on the Core. */
int catchem_core_get_num_processes(void* core_ptr);
int catchem_core_get_num_processes_checked(void* core_ptr, int* count_out);
int catchem_core_get_required_host_field_count_checked(void* core_ptr, int* count_out);
int catchem_core_get_required_host_field_name_checked(void* core_ptr, int index, char* name_out, int name_out_len);

// Grid and Configuration API
void catchem_get_grid_dimensions(void* core_ptr, int* nx, int* ny, int* nz);
double catchem_get_config_timestep(void* core_ptr);
int catchem_config_get_output_frequency(void* core_ptr);
int catchem_config_get_compress_level(void* core_ptr);
void catchem_config_get_output_directory(void* core_ptr, char* buffer, int max_len);
void catchem_config_get_output_prefix(void* core_ptr, char* buffer, int max_len);
int catchem_config_get_latlon_output(void* core_ptr);
int catchem_config_get_diag_enabled(void* core_ptr);
int catchem_config_get_process_diagnostics_enabled(void* core_ptr);
int catchem_config_get_diag_species_count(void* core_ptr);
void catchem_config_get_diag_species_at(void* core_ptr, int index, char* buffer, int max_len);
/// diagnostics.output.attributes (feature 013, FR-011): run-level NetCDF global
/// attributes.  Iteration order is the map's key order, so it is deterministic.
int catchem_config_get_output_attribute_count(void* core_ptr);
void catchem_config_get_output_attribute_key_at(void* core_ptr, int index, char* buffer, int max_len);
void catchem_config_get_output_attribute_value_at(void* core_ptr, int index, char* buffer, int max_len);
/// Path of the YAML the configuration was loaded from ("" before load).
void catchem_config_get_config_file_path(void* core_ptr, char* buffer, int max_len);
/// Build provenance baked in at configure time (feature 013, FR-011).
void catchem_get_build_version(char* buffer, int max_len);
void catchem_get_build_commit(char* buffer, int max_len);
int catchem_config_get_process_active(void* core_ptr, const char* process_name);
int catchem_config_has_emission_mapping(void* core_ptr);
int catchem_config_get_emission_category_count(void* core_ptr);
void catchem_config_get_emission_category_name_at(void* core_ptr, int index, char* name_out, int max_len);
int catchem_config_is_emission_category_active(void* core_ptr, const char* category_name);
int catchem_config_get_emission_field_count(void* core_ptr, const char* category_name);
void catchem_config_get_emission_field_name_at(void* core_ptr, const char* category_name, int field_idx, char* name_out,
                                               int max_len);
void catchem_config_get_emission_field_units(void* core_ptr, const char* category_name, const char* field_name,
                                             char* units_out, int max_len);
int catchem_config_get_emission_species_map_count(void* core_ptr, const char* category_name, const char* field_name);
void catchem_config_get_emission_species_map_at(void* core_ptr, const char* category_name, const char* field_name,
                                                int map_idx, char* target_species_out, int max_len, double* scale_out,
                                                int* species_idx_out);

/**
 * @brief Safely queries a boolean configuration setting by path.
 * @param core_ptr Pointer to catchem::Core instance.
 * @param yaml_path Slash-separated path in config (e.g. "processes/extemis/activate").
 * @param default_val Fallback value if key is missing or incompatible.
 * @return 1 for true, 0 for false, or default_val.
 */
int catchem_config_get_yaml_bool(void* core_ptr, const char* yaml_path, int default_val);

/**
 * @brief Safely queries a double configuration setting by path.
 * @param core_ptr Pointer to catchem::Core instance.
 * @param yaml_path Slash-separated path in config.
 * @param default_val Fallback value if key is missing.
 * @return Double configuration value or default_val.
 */
double catchem_config_get_yaml_double(void* core_ptr, const char* yaml_path, double default_val);

/**
 * @brief Safely queries an integer configuration setting by path.
 * @param core_ptr Pointer to catchem::Core instance.
 * @param yaml_path Slash-separated path in config.
 * @param default_val Fallback value if key is missing.
 * @return Integer configuration value or default_val.
 */
int catchem_config_get_yaml_int(void* core_ptr, const char* yaml_path, int default_val);

/**
 * @brief Safely queries a string configuration setting by path, resolving relative file paths.
 * @param core_ptr Pointer to catchem::Core instance.
 * @param yaml_path Slash-separated path in config.
 * @param val_out Output buffer for null-terminated string.
 * @param max_len Maximum length of output buffer.
 * @param default_val Fallback string if key is missing.
 */
void catchem_config_get_yaml_string(void* core_ptr, const char* yaml_path, char* val_out, int max_len,
                                    const char* default_val);

/**
 * @brief Locates static FENGSHA/dust input file path via ConfigManager process settings.
 * @param core_ptr Pointer to catchem::Core instance.
 * @param val_out Output buffer for null-terminated string.
 * @param max_len Maximum length of output buffer.
 */
void catchem_config_find_fengsha_static_file(void* core_ptr, char* val_out, int max_len);

/**
 * @brief Returns sequence length for a list configuration path.
 * @param core_ptr Pointer to catchem::Core instance.
 * @param yaml_path Slash-separated path in config.
 * @return Number of elements in sequence, or 0.
 */
int catchem_config_get_yaml_list_count(void* core_ptr, const char* yaml_path);

/**
 * @brief Retrieves a string element at index from a list configuration path.
 * @param core_ptr Pointer to catchem::Core instance.
 * @param yaml_path Slash-separated path in config.
 * @param index 0-based element index.
 * @param val_out Output buffer for string element.
 * @param max_len Maximum length of output buffer.
 */
void catchem_config_get_yaml_list_at(void* core_ptr, const char* yaml_path, int index, char* val_out, int max_len);

/**
 * @brief Enable or disable echoing the parsed YAML to stdout on configuration load.
 * @param enabled Non-zero to echo (the default), zero to suppress.
 *
 * The NUOPC driver calls this with zero on non-root PETs so the effective
 * configuration appears once per coupled run instead of once per PET.
 */
void catchem_set_config_echo_enabled(int enabled);

// Diagnostic API
void catchem_diag_register(void* core_ptr, const char* name, const char* desc, const char* units, int rank, int dim1,
                           int dim2, int dim3);
int catchem_diag_register_checked(void* core_ptr, const char* name, const char* desc, const char* units, int rank,
                                  int dim1, int dim2, int dim3);
int catchem_diag_register_contract_checked(void* core_ptr, const char* name, const char* desc, const char* units,
                                           int rank, const int* dims, const int* axes, int policy, double reset_value);
int catchem_diag_get_contract(void* core_ptr, const char* name, int* generation, int* availability, int* latest_writer,
                              int* policy);
void* catchem_diag_get_pointer(void* core_ptr, const char* name);
int catchem_diag_get_rank(void* core_ptr, const char* name);
int catchem_diag_get_rank_checked(void* core_ptr, const char* name, int* rank_out);
void catchem_diag_get_dims(void* core_ptr, const char* name, int* dims_out);
int catchem_diag_get_dims_checked(void* core_ptr, const char* name, int* dims_out, int dims_length);
int catchem_diag_get_pointer_checked(void* core_ptr, const char* name, int rank, const int* dims, void** ptr_out);
int catchem_diag_mark_host_modified(void* core_ptr, const char* name);
int catchem_diag_mark_device_modified(void* core_ptr, const char* name);
void catchem_diag_sync_to_host(void* core_ptr);
void catchem_diag_reset(void* core_ptr);
int catchem_diag_get_count(void* core_ptr);
int catchem_diag_get_count_checked(void* core_ptr, int* count_out);
void catchem_diag_get_name_at(void* core_ptr, int index, char* name_out);
int catchem_diag_get_name_at_checked(void* core_ptr, int index, char* name_out, int name_length);
int catchem_diag_get_units_checked(void* core_ptr, const char* name, char* units_out, int units_length);
int catchem_diag_get_description_checked(void* core_ptr, const char* name, char* desc_out, int desc_length);
/// SemanticAxis ordinals per dimension (length == rank). axes_length must be >= rank.
int catchem_diag_get_axes_checked(void* core_ptr, const char* name, int* axes_out, int axes_length);
/// Label for slot `slot` (0-based) of the field's packed (Species/Category) dimension.
/// Non-zero when the field has no packed dimension or `slot` is out of range.
int catchem_diag_get_unpack_label_at_checked(void* core_ptr, const char* name, int slot, char* label_out,
                                             int label_length);

// YAML Species Metadata
void catchem_state_load_species_config(void* state_ptr, const char* filename);
int catchem_state_get_species_count(void* state_ptr);
int catchem_state_get_species_count_checked(void* state_ptr, int* count_out);
int catchem_state_get_species_index(void* state_ptr,
                                    const char* name); // returns 1-based index matching Fortran, or -1 if not found
int catchem_state_get_species_index_checked(void* state_ptr, const char* name, int* index_out);

// Categorized counts and list getters
int catchem_state_get_gas_species_count(void* state_ptr);
void catchem_state_get_gas_indices(void* state_ptr, int* indices_out); // populates 1-based indices
int catchem_state_get_aerosol_species_count(void* state_ptr);
void catchem_state_get_aerosol_indices(void* state_ptr, int* indices_out);

// Individual property getters (by 1-based index)
double catchem_state_get_species_mw(void* state_ptr, int index);
int catchem_state_get_species_mw_checked(void* state_ptr, int index, double* molecular_weight_out);
int catchem_state_is_species_gas(void* state_ptr, int index);
int catchem_state_is_species_gas_checked(void* state_ptr, int index, int* value_out);
int catchem_state_is_species_aerosol(void* state_ptr, int index);
int catchem_state_is_species_aerosol_checked(void* state_ptr, int index, int* value_out);
void catchem_state_get_species_name_at(void* state_ptr, int index, char* name_out);
int catchem_state_get_species_name_at_checked(void* state_ptr, int index, char* name_out, int name_length);
int catchem_state_get_species_is_advected_checked(void* state_ptr, int index, int* value_out);
void catchem_state_get_species_long_name_at(void* state_ptr, int index, char* name_out);
void catchem_state_get_species_desc_at(void* state_ptr, int index, char* desc_out);
void catchem_state_get_species_mie_name(void* state_ptr, int index, char* mie_out);

// Physics derivations
void catchem_state_derive_bxheight(void* state_ptr);
int catchem_state_derive_bxheight_checked(void* state_ptr);
void catchem_state_derive_airden_dry(void* state_ptr);
int catchem_state_derive_airden_dry_checked(void* state_ptr);

// TimeState C-Linkable API
void* catchem_time_state_create();
void catchem_time_state_destroy(void* ptr);
int catchem_time_state_init(void* ptr, int year, int month, int day, int hour, int minute, int second, double timestep);
int catchem_time_state_advance(void* ptr, double dt);
int catchem_time_state_reset(void* ptr);
int catchem_time_state_get_year(void* ptr);
int catchem_time_state_get_month(void* ptr);
int catchem_time_state_get_day(void* ptr);
int catchem_time_state_get_hour(void* ptr);
int catchem_time_state_get_minute(void* ptr);
int catchem_time_state_get_second(void* ptr);
double catchem_time_state_get_timestep(void* ptr);
double catchem_time_state_get_julian_date(void* ptr);
int catchem_time_state_get_doy(void* ptr);
double catchem_time_state_get_cos_sza(void* ptr, double lat, double lon, bool mid_timestep);
int catchem_time_state_get_timezone_offset(void* ptr, double lon);
bool catchem_time_state_is_leap_year(int year);
int catchem_time_state_get_days_in_month(int month, int year);
bool catchem_time_state_is_global_holiday(int month, int day);
bool catchem_time_state_is_us_holiday(int month, int day);

// UnitConversion C-Linkable API
double catchem_convert_concentration(double val, const char* from_units, const char* to_units, double mw, double temp,
                                     double press, int* rc);
double catchem_convert_pressure(double val, const char* from_units, const char* to_units, int* rc);
double catchem_convert_temperature(double val, const char* from_units, const char* to_units, int* rc);
double catchem_convert_flux(double val, const char* from_units, const char* to_units, double mw, int* rc);
double catchem_convert_rate_constant(double val, const char* from_units, const char* to_units, int* rc);
double catchem_convert_mass_units(double val, const char* from_units, const char* to_units, int* rc);
double catchem_calculate_air_density(double temp, double press, double humidity, bool use_humidity);
double catchem_calculate_molecular_weight(const char* formula);
double catchem_convert_imperial(double val, const char* from_units, const char* to_units, const char* category,
                                int* rc);
int catchem_convert_process_concentration_units(catchem::fp* values, int size, const char* from_units,
                                                const char* to_units, catchem::fp mw, catchem::fp temp,
                                                catchem::fp press);
int catchem_convert_process_flux_units(catchem::fp* values, int size, const char* from_units, const char* to_units,
                                       catchem::fp mw);

// =========================================================================
// Species Metadata and Property Query C-API
// =========================================================================
void catchem_state_get_species_name_at(void* state_ptr, int index, char* name_out);
void catchem_state_get_species_long_name_at(void* state_ptr, int index, char* name_out);
void catchem_state_get_species_desc_at(void* state_ptr, int index, char* desc_out);
#define CATCHEM_SPECIES_DOUBLE_PROPERTY(api, member, fallback)                                                         \
    double catchem_state_get_species_##api(void* state_ptr, int index);
#define CATCHEM_SPECIES_BOOL_PROPERTY(api, member) int catchem_state_get_species_##api(void* state_ptr, int index);
#define CATCHEM_SPECIES_LEGACY_BOOL_PROPERTY(api, member)                                                              \
    int catchem_state_is_species_##api(void* state_ptr, int index);
#include "catchem_species_properties.def"
#undef CATCHEM_SPECIES_LEGACY_BOOL_PROPERTY
#undef CATCHEM_SPECIES_BOOL_PROPERTY
#undef CATCHEM_SPECIES_DOUBLE_PROPERTY

void catchem_state_get_species_wd_rainouteff(void* state_ptr, int index, double* eff_out);
void catchem_state_get_species_mie_name(void* state_ptr, int index, char* name_out);

// =========================================================================
// Meteorological Core Calculation C-API
// =========================================================================
double catchem_met_potential_temperature(double temp, double press, double sfc_press);
double catchem_met_virtual_temperature(double temp, double qv);
double catchem_met_dew_point(double temp, double rh);
double catchem_met_relative_humidity(double temp, double qv, double press);
double catchem_met_saturation_vapor_pressure(double temp);
double catchem_met_monin_obukhov_length(double ustar, double t0, double hflux, double rho);
double catchem_met_friction_velocity(double tau, double rho);
double catchem_met_cunningham_correction_factor(double dp, double lambda);
double catchem_met_mean_free_path_air(double temp, double press);
void catchem_met_solar_zenith_angle(int doy, double hour, double lat_rad, double lon_rad, double* sza_deg,
                                    double* cossza);
double catchem_met_mixing_ratio(double q);
double catchem_met_specific_humidity(double r);
double catchem_met_dry_adiabatic_lapse_rate();
double catchem_met_bulk_richardson_number(double t0, double tz, double u, double z);
int catchem_met_stability_classification(double l);
double catchem_met_saturation_mixing_ratio(double p, double t);
double catchem_met_latent_heat_vaporization(double t);
double catchem_met_psychrometric_constant(double p, double lv);
double catchem_met_wind_profile_loglaw(double ustar, double z, double z0);
double catchem_met_brunt_vaisala_frequency(double t0, double dtdz);
double catchem_met_psi_m_businger(double zeta);
double catchem_met_psi_h_businger(double zeta);
double catchem_met_arrhenius_rate(double a, double ea, double t);
double catchem_met_henrys_law_constant(double h0, double dh, double t, double t0);
double catchem_met_photolysis_rate_scaling(double j0, double sza);
double catchem_met_ppm_to_ugm3(double ppm, double m, double t, double p);
double catchem_met_ugm3_to_ppm(double ugm3, double m, double t, double p);
double catchem_met_stokes_settling_velocity(double dp, double rho_p, double rho_a, double mu, double cc);
double catchem_met_stokes_number(double rho_p, double d_p, double u, double mu, double l);
double catchem_met_nuclear_decay(double n0, double lambda, double t);

/** @brief Get grid X dimension. */
int catchem_state_get_nx(void* state_ptr);

/** @brief Get grid Y dimension. */
int catchem_state_get_ny(void* state_ptr);

/** @brief Get grid Z dimension. */
int catchem_state_get_nz(void* state_ptr);

#ifdef __cplusplus
}
#endif
