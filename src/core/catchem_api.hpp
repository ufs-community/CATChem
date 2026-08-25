/**
 * @file catchem_api.hpp
 * @brief Flat, BIND(C) linkable API endpoints for CATChem host model integration.
 */

#pragma once

#include "catchem_precision.hpp"

#ifdef __cplusplus
extern "C" {
#endif

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

/** @brief Bind a 2D double field to the StateManager registry. */
void catchem_state_bind_2d(void* state_ptr, const char* name, double* ptr);

/** @brief Bind a 3D double field to the StateManager registry. */
void catchem_state_bind_3d(void* state_ptr, const char* name, double* ptr);

/** @brief Bind a 2D meteorological field by name. */
void catchem_state_bind_met_2d(void* state_ptr, const char* name, double* ptr);

/** @brief Bind a 3D meteorological field by name. */
void catchem_state_bind_met_3d(void* state_ptr, const char* name, double* ptr);

/** @brief Binds the contiguous, multi-species unified chemistry concentrations array. */
void catchem_state_bind_unified_chemistry(void* state_ptr, double* ptr);

/** @brief Sets current simulation time within the state. */
void catchem_state_set_time(void* state_ptr, int yr, int mo, int dy, int hr, int mn, int sc, int doy, double tstep);

/** @brief Synchronizes registered host pointers to Kokkos device memory space. */
void catchem_state_sync_to_device(void* state_ptr);

/** @brief Synchronizes Kokkos device calculations back to host buffers. */
void catchem_state_sync_to_host(void* state_ptr);

/** @brief Retrieves direct host pointers from the 1D, 2D, or 3D fields. */
double* catchem_state_get_pointer_1d(void* state_ptr, const char* name);
double* catchem_state_get_pointer_2d(void* state_ptr, const char* name);
double* catchem_state_get_pointer_3d(void* state_ptr, const char* name);

/**
 * @brief Executes a single timestepped execution over scheduled processes.
 * @param core_ptr Core pointer.
 * @param dt Step size in seconds.
 */
void catchem_core_run_timestep(void* core_ptr, double dt);

/** @brief Registers and attaches an active physics process handler. */
void catchem_core_add_process_by_name(void* core_ptr, const char* name);

// Grid and Configuration API
void catchem_get_grid_dimensions(void* core_ptr, int* nx, int* ny, int* nz);
double catchem_get_config_timestep(void* core_ptr);

// Diagnostic API
void catchem_diag_register(void* core_ptr, const char* name, const char* desc, const char* units, int rank, int dim1,
                           int dim2, int dim3);
void* catchem_diag_get_pointer(void* core_ptr, const char* name);
void catchem_diag_sync_to_host(void* core_ptr);
void catchem_diag_reset(void* core_ptr);
int catchem_diag_get_count(void* core_ptr);
void catchem_diag_get_name_at(void* core_ptr, int index, char* name_out);

// YAML Species Metadata
void catchem_state_load_species_config(void* state_ptr, const char* filename);
int catchem_state_get_species_count(void* state_ptr);
int catchem_state_get_species_index(void* state_ptr,
                                    const char* name); // returns 1-based index matching Fortran, or -1 if not found

// Categorized counts and list getters
int catchem_state_get_gas_species_count(void* state_ptr);
void catchem_state_get_gas_indices(void* state_ptr, int* indices_out); // populates 1-based indices
int catchem_state_get_aerosol_species_count(void* state_ptr);
void catchem_state_get_aerosol_indices(void* state_ptr, int* indices_out);

// Individual property getters (by 1-based index)
double catchem_state_get_species_mw(void* state_ptr, int index);
int catchem_state_is_species_gas(void* state_ptr, int index);
int catchem_state_is_species_aerosol(void* state_ptr, int index);

// Physics derivations
void catchem_state_derive_bxheight(void* state_ptr);
void catchem_state_derive_airden_dry(void* state_ptr);

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
double catchem_state_get_species_density(void* state_ptr, int index);
double catchem_state_get_species_radius(void* state_ptr, int index);
double catchem_state_get_species_lower_radius(void* state_ptr, int index);
double catchem_state_get_species_upper_radius(void* state_ptr, int index);
double catchem_state_get_species_viscosity(void* state_ptr, int index);
int catchem_state_get_species_is_tracer(void* state_ptr, int index);
int catchem_state_get_species_is_advected(void* state_ptr, int index);
int catchem_state_get_species_is_drydep(void* state_ptr, int index);
int catchem_state_get_species_is_wetdep(void* state_ptr, int index);
int catchem_state_get_species_is_photolysis(void* state_ptr, int index);
int catchem_state_get_species_is_dust(void* state_ptr, int index);
int catchem_state_get_species_is_seasalt(void* state_ptr, int index);

double catchem_state_get_species_dd_f0(void* state_ptr, int index);
double catchem_state_get_species_dd_hstar(void* state_ptr, int index);
double catchem_state_get_species_dd_DvzAerSnow(void* state_ptr, int index);
double catchem_state_get_species_dd_DvzMinVal_snow(void* state_ptr, int index);
double catchem_state_get_species_dd_DvzMinVal_land(void* state_ptr, int index);

double catchem_state_get_species_henry_k0(void* state_ptr, int index);
double catchem_state_get_species_henry_cr(void* state_ptr, int index);
double catchem_state_get_species_henry_pKa(void* state_ptr, int index);
double catchem_state_get_species_wd_retfactor(void* state_ptr, int index);
int catchem_state_get_species_wd_LiqAndGas(void* state_ptr, int index);
double catchem_state_get_species_wd_convfacI2G(void* state_ptr, int index);
void catchem_state_get_species_wd_rainouteff(void* state_ptr, int index, double* eff_out);
double catchem_state_get_species_wd_reevap_frac(void* state_ptr, int index);
double catchem_state_get_species_t_chem_loss(void* state_ptr, int index);
double catchem_state_get_species_BackgroundVV(void* state_ptr, int index);
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

#ifdef __cplusplus
}
#endif
