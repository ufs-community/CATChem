#include "catchem_api.hpp"
#include "catchem_api_internal.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_state_manager.hpp"
#include "catchem_unit_conversion.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

using namespace catchem::api_internal;

extern "C" {
void catchem_state_load_species_config(void* state_ptr, const char* filename) {
    try {
        auto* state = static_cast<catchem::StateManager*>(state_ptr);
        state->load_species_config(filename);
    } catch (const std::exception& e) {
        std::cerr << "CATChem API Error: Failed to load species configuration '" << filename
                  << "'. Details: " << e.what() << std::endl;
    }
}

int catchem_state_get_species_count(void* state_ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    return static_cast<int>(state->chemistry().species_list.size());
}

int catchem_state_get_species_count_checked(void* state_ptr, int* count_out) {
    if (!count_out)
        return fail(catchem::BoundaryStatus::NullArgument, "state_get_species_count", "count_out",
                    "output pointer is null");
    *count_out = 0;
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, "state_get_species_count", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    *count_out = static_cast<int>(static_cast<catchem::StateManager*>(state_ptr)->chemistry().species_list.size());
    return CATCHEM_SUCCESS;
}

int catchem_state_get_species_index(void* state_ptr, const char* name) {
    if (!state_ptr || !name)
        return -1;
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    std::string canonical_name(name);
    std::transform(canonical_name.begin(), canonical_name.end(), canonical_name.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    auto it = state->chemistry().species_name_to_index.find(canonical_name);
    if (it != state->chemistry().species_name_to_index.end()) {
        return it->second + 1; // Translate 0-based C++ index to 1-based Fortran index
    }
    return -1;
}

int catchem_state_get_species_index_checked(void* state_ptr, const char* name, int* index_out) {
    if (index_out)
        *index_out = 0;
    if (!name || !index_out)
        return fail(catchem::BoundaryStatus::NullArgument, "state_get_species_index", "argument",
                    "species name and index output are required");
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, "state_get_species_index", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    std::string key(name);
    std::transform(key.begin(), key.end(), key.begin(), [](unsigned char c) { return std::toupper(c); });
    const auto& indices = static_cast<catchem::StateManager*>(state_ptr)->chemistry().species_name_to_index;
    const auto found = indices.find(key);
    if (found == indices.end())
        return fail(catchem::BoundaryStatus::InvalidIndex, "state_get_species_index", name,
                    "species is not present in the runtime mechanism");
    *index_out = found->second + 1;
    return CATCHEM_SUCCESS;
}

int catchem_state_get_gas_species_count(void* state_ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    return static_cast<int>(state->chemistry().gas_indices.size());
}

void catchem_state_get_gas_indices(void* state_ptr, int* indices_out) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    for (size_t i = 0; i < state->chemistry().gas_indices.size(); ++i) {
        indices_out[i] = state->chemistry().gas_indices[i] + 1; // 1-based
    }
}

int catchem_state_get_aerosol_species_count(void* state_ptr) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    return static_cast<int>(state->chemistry().aerosol_indices.size());
}

void catchem_state_get_aerosol_indices(void* state_ptr, int* indices_out) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    for (size_t i = 0; i < state->chemistry().aerosol_indices.size(); ++i) {
        indices_out[i] = state->chemistry().aerosol_indices[i] + 1; // 1-based
    }
}

double catchem_state_get_species_mw(void* state_ptr, int index) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    int idx_0 = index - 1; // 1-based to 0-based
    if (idx_0 >= 0 && idx_0 < static_cast<int>(state->chemistry().species_list.size())) {
        return state->chemistry().species_list[idx_0].mw_g;
    }
    return 0.0;
}

int catchem_state_get_species_mw_checked(void* state_ptr, int index, double* molecular_weight_out) {
    if (!molecular_weight_out)
        return fail(catchem::BoundaryStatus::NullArgument, "state_get_species_mw", "output",
                    "molecular-weight output is null");
    *molecular_weight_out = 0.0;
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, "state_get_species_mw", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    const auto& species = static_cast<catchem::StateManager*>(state_ptr)->chemistry().species_list;
    if (index < 1 || index > static_cast<int>(species.size()))
        return fail(catchem::BoundaryStatus::InvalidIndex, "state_get_species_mw", "index",
                    "species index is outside the runtime mechanism");
    *molecular_weight_out = species[static_cast<std::size_t>(index - 1)].mw_g;
    return CATCHEM_SUCCESS;
}

static int checked_species_classification(void* state_ptr, int index, int* value_out, const char* operation,
                                          bool catchem::SpeciesMetadata::*member) {
    if (!value_out)
        return fail(catchem::BoundaryStatus::NullArgument, operation, "output", "output pointer is null");
    *value_out = 0;
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, operation, admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    const auto& species = static_cast<catchem::StateManager*>(state_ptr)->chemistry().species_list;
    if (index < 1 || index > static_cast<int>(species.size()))
        return fail(catchem::BoundaryStatus::InvalidIndex, operation, "index",
                    "species index is outside the runtime mechanism");
    *value_out = species[static_cast<std::size_t>(index - 1)].*member ? 1 : 0;
    return CATCHEM_SUCCESS;
}

int catchem_state_is_species_gas_checked(void* state_ptr, int index, int* value_out) {
    return checked_species_classification(state_ptr, index, value_out, "state_is_species_gas",
                                          &catchem::SpeciesMetadata::is_gas);
}

int catchem_state_is_species_aerosol_checked(void* state_ptr, int index, int* value_out) {
    return checked_species_classification(state_ptr, index, value_out, "state_is_species_aerosol",
                                          &catchem::SpeciesMetadata::is_aerosol);
}

int catchem_state_get_species_is_advected_checked(void* state_ptr, int index, int* value_out) {
    return checked_species_classification(state_ptr, index, value_out, "state_get_species_is_advected",
                                          &catchem::SpeciesMetadata::is_advected);
}

int catchem_state_is_species_gas(void* state_ptr, int index) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    int idx_0 = index - 1;
    if (idx_0 >= 0 && idx_0 < static_cast<int>(state->chemistry().species_list.size())) {
        return state->chemistry().species_list[idx_0].is_gas ? 1 : 0;
    }
    return 0;
}

int catchem_state_is_species_aerosol(void* state_ptr, int index) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    int idx_0 = index - 1;
    if (idx_0 >= 0 && idx_0 < static_cast<int>(state->chemistry().species_list.size())) {
        return state->chemistry().species_list[idx_0].is_aerosol ? 1 : 0;
    }
    return 0;
}

void catchem_state_get_species_name_at(void* state_ptr, int index, char* name_out) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    int idx_0 = index - 1;
    if (idx_0 >= 0 && idx_0 < static_cast<int>(state->chemistry().species_list.size())) {
        copy_string_to_buffer(state->chemistry().species_list[idx_0].short_name, name_out, 64);
    } else {
        copy_string_to_buffer("", name_out, 64);
    }
}

int catchem_state_get_species_name_at_checked(void* state_ptr, int index, char* name_out, int name_length) {
    if (name_out && name_length > 0)
        name_out[0] = '\0';
    if (!name_out || name_length <= 0)
        return fail(catchem::BoundaryStatus::NullArgument, "state_get_species_name", "name_out",
                    "a positive-length name output is required");
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, "state_get_species_name", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    const auto& species = static_cast<catchem::StateManager*>(state_ptr)->chemistry().species_list;
    if (index < 1 || index > static_cast<int>(species.size()))
        return fail(catchem::BoundaryStatus::InvalidIndex, "state_get_species_name", "index",
                    "species index is outside the runtime mechanism");
    copy_string_to_buffer(species[static_cast<std::size_t>(index - 1)].short_name, name_out, name_length);
    return CATCHEM_SUCCESS;
}

void catchem_state_get_species_long_name_at(void* state_ptr, int index, char* name_out) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    int idx_0 = index - 1;
    if (idx_0 >= 0 && idx_0 < static_cast<int>(state->chemistry().species_list.size())) {
        copy_string_to_buffer(state->chemistry().species_list[idx_0].long_name, name_out, 128);
    } else {
        copy_string_to_buffer("", name_out, 128);
    }
}

void catchem_state_get_species_desc_at(void* state_ptr, int index, char* desc_out) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    int idx_0 = index - 1;
    if (idx_0 >= 0 && idx_0 < static_cast<int>(state->chemistry().species_list.size())) {
        copy_string_to_buffer(state->chemistry().species_list[idx_0].description, desc_out, 256);
    } else {
        copy_string_to_buffer("", desc_out, 256);
    }
}

void catchem_state_get_species_mie_name(void* state_ptr, int index, char* mie_out) {
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    int idx_0 = index - 1;
    if (idx_0 >= 0 && idx_0 < static_cast<int>(state->chemistry().species_list.size())) {
        copy_string_to_buffer(state->chemistry().species_list[idx_0].mie_name, mie_out, 64);
    } else {
        copy_string_to_buffer("", mie_out, 64);
    }
}

void catchem_get_grid_dimensions(void* core_ptr, int* nx, int* ny, int* nz) {
    if (core_ptr == nullptr)
        return;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    auto grid = core->get_grid_manager();
    if (nx)
        *nx = grid->geometry.nx;
    if (ny)
        *ny = grid->geometry.ny;
    if (nz)
        *nz = grid->geometry.nz;
}

double catchem_get_config_timestep(void* core_ptr) {
    if (core_ptr == nullptr)
        return 0.0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return core->get_config_manager()->data.runtime.dt;
}

int catchem_config_get_output_frequency(void* core_ptr) {
    if (core_ptr == nullptr)
        return 0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return core->get_config_manager()->data.diagnostics.output.frequency;
}

int catchem_config_get_compress_level(void* core_ptr) {
    if (core_ptr == nullptr)
        return 0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return core->get_config_manager()->data.diagnostics.output.compress_lev;
}

void catchem_config_get_output_directory(void* core_ptr, char* buffer, int max_len) {
    if (core_ptr == nullptr) {
        copy_string_to_buffer("./", buffer, max_len);
        return;
    }
    auto* core = static_cast<catchem::Core*>(core_ptr);
    copy_string_to_buffer(core->get_config_manager()->data.diagnostics.output.directory, buffer, max_len);
}

void catchem_config_get_output_prefix(void* core_ptr, char* buffer, int max_len) {
    if (core_ptr == nullptr) {
        copy_string_to_buffer("catchem_diag", buffer, max_len);
        return;
    }
    auto* core = static_cast<catchem::Core*>(core_ptr);
    copy_string_to_buffer(core->get_config_manager()->data.diagnostics.output.prefix, buffer, max_len);
}

int catchem_config_get_latlon_output(void* core_ptr) {
    if (core_ptr == nullptr)
        return 0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return core->get_config_manager()->data.diagnostics.output.format == "latlon" ? 1 : 0;
}

void catchem_set_config_echo_enabled(int enabled) {
    catchem::ConfigManager::echo_config_to_stdout = (enabled != 0);
}

int catchem_config_get_diag_enabled(void* core_ptr) {
    if (core_ptr == nullptr)
        return 0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return core->get_config_manager()->data.diagnostics.output.enabled ? 1 : 0;
}

int catchem_config_get_process_diagnostics_enabled(void* core_ptr) {
    if (core_ptr == nullptr)
        return 0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return core->get_config_manager()->data.diagnostics.output.process_diagnostics ? 1 : 0;
}

int catchem_config_get_diag_species_count(void* core_ptr) {
    if (core_ptr == nullptr)
        return 0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return static_cast<int>(core->get_config_manager()->data.diagnostics.output.diag_list.size());
}

void catchem_config_get_diag_species_at(void* core_ptr, int index, char* buffer, int max_len) {
    if (core_ptr == nullptr)
        return;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    const auto& diag_list = core->get_config_manager()->data.diagnostics.output.diag_list;
    if (index >= 0 && index < static_cast<int>(diag_list.size())) {
        copy_string_to_buffer(diag_list[index], buffer, max_len);
    } else {
        copy_string_to_buffer("", buffer, max_len);
    }
}

int catchem_config_get_output_attribute_count(void* core_ptr) {
    if (core_ptr == nullptr)
        return 0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return static_cast<int>(core->get_config_manager()->data.diagnostics.output.attributes.size());
}

void catchem_config_get_output_attribute_key_at(void* core_ptr, int index, char* buffer, int max_len) {
    if (core_ptr == nullptr)
        return;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    const auto& attributes = core->get_config_manager()->data.diagnostics.output.attributes;
    auto it = attributes.begin();
    std::advance(it, index);
    if (index >= 0 && it != attributes.end()) {
        copy_string_to_buffer(it->first, buffer, max_len);
    } else {
        copy_string_to_buffer("", buffer, max_len);
    }
}

void catchem_config_get_output_attribute_value_at(void* core_ptr, int index, char* buffer, int max_len) {
    if (core_ptr == nullptr)
        return;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    const auto& attributes = core->get_config_manager()->data.diagnostics.output.attributes;
    auto it = attributes.begin();
    std::advance(it, index);
    if (index >= 0 && it != attributes.end()) {
        copy_string_to_buffer(it->second, buffer, max_len);
    } else {
        copy_string_to_buffer("", buffer, max_len);
    }
}

void catchem_config_get_config_file_path(void* core_ptr, char* buffer, int max_len) {
    if (core_ptr == nullptr)
        return;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    copy_string_to_buffer(core->get_config_manager()->config_file_path, buffer, max_len);
}

// Build provenance is injected by CMake (see src/core/CMakeLists.txt).  A source
// tree configured without those definitions still links; the attributes degrade
// to "unknown" instead of an empty string.
#ifndef CATCHEM_BUILD_VERSION
#define CATCHEM_BUILD_VERSION "unknown"
#endif
#ifndef CATCHEM_BUILD_COMMIT
#define CATCHEM_BUILD_COMMIT "unknown"
#endif

void catchem_get_build_version(char* buffer, int max_len) {
    copy_string_to_buffer(CATCHEM_BUILD_VERSION, buffer, max_len);
}

void catchem_get_build_commit(char* buffer, int max_len) {
    copy_string_to_buffer(CATCHEM_BUILD_COMMIT, buffer, max_len);
}

int catchem_config_get_process_active(void* core_ptr, const char* process_name) {
    if (core_ptr == nullptr || process_name == nullptr)
        return 0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    const auto& processes = core->get_config_manager()->data.processes;
    auto it = processes.find(std::string(process_name));
    if (it != processes.end()) {
        return it->second.activate ? 1 : 0;
    }
    return 0;
}

int catchem_config_has_emission_mapping(void* core_ptr) {
    if (core_ptr == nullptr)
        return 0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return !core->get_config_manager()->data.emission_mappings.empty() ? 1 : 0;
}

int catchem_config_get_emission_category_count(void* core_ptr) {
    if (core_ptr == nullptr)
        return 0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return static_cast<int>(core->get_config_manager()->data.emission_mappings.size());
}

void catchem_config_get_emission_category_name_at(void* core_ptr, int index, char* name_out, int max_len) {
    if (core_ptr == nullptr) {
        copy_string_to_buffer("", name_out, max_len);
        return;
    }
    auto* core = static_cast<catchem::Core*>(core_ptr);
    const auto& mappings = core->get_config_manager()->data.emission_mappings;
    if (index >= 0 && index < static_cast<int>(mappings.size())) {
        auto it = mappings.begin();
        std::advance(it, index);
        copy_string_to_buffer(it->first, name_out, max_len);
    } else {
        copy_string_to_buffer("", name_out, max_len);
    }
}

int catchem_config_is_emission_category_active(void* core_ptr, const char* category_name) {
    if (core_ptr == nullptr || category_name == nullptr)
        return 0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return core->get_config_manager()->is_category_active(category_name) ? 1 : 0;
}

int catchem_config_get_emission_field_count(void* core_ptr, const char* category_name) {
    if (core_ptr == nullptr || category_name == nullptr)
        return 0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    const auto& mappings = core->get_config_manager()->data.emission_mappings;
    auto it = mappings.find(std::string(category_name));
    if (it != mappings.end()) {
        return static_cast<int>(it->second.fields.size());
    }
    return 0;
}

void catchem_config_get_emission_field_name_at(void* core_ptr, const char* category_name, int field_idx, char* name_out,
                                               int max_len) {
    if (core_ptr == nullptr || category_name == nullptr) {
        copy_string_to_buffer("", name_out, max_len);
        return;
    }
    auto* core = static_cast<catchem::Core*>(core_ptr);
    const auto& mappings = core->get_config_manager()->data.emission_mappings;
    auto it = mappings.find(std::string(category_name));
    if (it != mappings.end() && field_idx >= 0 && field_idx < static_cast<int>(it->second.fields.size())) {
        auto fit = it->second.fields.begin();
        std::advance(fit, field_idx);
        copy_string_to_buffer(fit->first, name_out, max_len);
    } else {
        copy_string_to_buffer("", name_out, max_len);
    }
}

void catchem_config_get_emission_field_units(void* core_ptr, const char* category_name, const char* field_name,
                                             char* units_out, int max_len) {
    if (core_ptr == nullptr || category_name == nullptr || field_name == nullptr) {
        copy_string_to_buffer("", units_out, max_len);
        return;
    }
    auto* core = static_cast<catchem::Core*>(core_ptr);
    const auto& mappings = core->get_config_manager()->data.emission_mappings;
    auto it = mappings.find(std::string(category_name));
    if (it != mappings.end()) {
        auto fit = it->second.fields.find(std::string(field_name));
        if (fit != it->second.fields.end()) {
            copy_string_to_buffer(fit->second.units, units_out, max_len);
            return;
        }
    }
    copy_string_to_buffer("", units_out, max_len);
}

int catchem_config_get_emission_species_map_count(void* core_ptr, const char* category_name, const char* field_name) {
    if (core_ptr == nullptr || category_name == nullptr || field_name == nullptr)
        return 0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    const auto& mappings = core->get_config_manager()->data.emission_mappings;
    auto it = mappings.find(std::string(category_name));
    if (it != mappings.end()) {
        auto fit = it->second.fields.find(std::string(field_name));
        if (fit != it->second.fields.end()) {
            return static_cast<int>(fit->second.map.size());
        }
    }
    return 0;
}

void catchem_config_get_emission_species_map_at(void* core_ptr, const char* category_name, const char* field_name,
                                                int map_idx, char* target_species_out, int max_len, double* scale_out,
                                                int* species_idx_out) {
    if (core_ptr == nullptr || category_name == nullptr || field_name == nullptr) {
        copy_string_to_buffer("", target_species_out, max_len);
        if (scale_out)
            *scale_out = 1.0;
        if (species_idx_out)
            *species_idx_out = -1;
        return;
    }
    auto* core = static_cast<catchem::Core*>(core_ptr);
    const auto& mappings = core->get_config_manager()->data.emission_mappings;
    auto it = mappings.find(std::string(category_name));
    if (it != mappings.end()) {
        auto fit = it->second.fields.find(std::string(field_name));
        if (fit != it->second.fields.end() && map_idx >= 0 && map_idx < static_cast<int>(fit->second.map.size())) {
            copy_string_to_buffer(fit->second.map[map_idx], target_species_out, max_len);
            if (scale_out) {
                *scale_out = (map_idx < static_cast<int>(fit->second.scale.size())) ? fit->second.scale[map_idx] : 1.0;
            }
            if (species_idx_out) {
                *species_idx_out =
                    catchem_state_get_species_index(core->get_state_manager().get(), fit->second.map[map_idx].c_str());
            }
            return;
        }
    }
    copy_string_to_buffer("", target_species_out, max_len);
    if (scale_out)
        *scale_out = 1.0;
    if (species_idx_out)
        *species_idx_out = -1;
}

int catchem_config_get_yaml_bool(void* core_ptr, const char* yaml_path, int default_val) {
    if (core_ptr == nullptr || yaml_path == nullptr)
        return default_val;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return core->get_config_manager()->get_bool(yaml_path, default_val != 0) ? 1 : 0;
}

double catchem_config_get_yaml_double(void* core_ptr, const char* yaml_path, double default_val) {
    if (core_ptr == nullptr || yaml_path == nullptr)
        return default_val;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return core->get_config_manager()->get_double(yaml_path, default_val);
}

int catchem_config_get_yaml_int(void* core_ptr, const char* yaml_path, int default_val) {
    if (core_ptr == nullptr || yaml_path == nullptr)
        return default_val;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return core->get_config_manager()->get_int(yaml_path, default_val);
}

void catchem_config_get_yaml_string(void* core_ptr, const char* yaml_path, char* val_out, int max_len,
                                    const char* default_val) {
    if (core_ptr == nullptr || yaml_path == nullptr) {
        copy_string_to_buffer(default_val ? default_val : "", val_out, max_len);
        return;
    }
    auto* core = static_cast<catchem::Core*>(core_ptr);
    std::string result = core->get_config_manager()->get_string(yaml_path, default_val ? default_val : "");
    copy_string_to_buffer(result, val_out, max_len);
}

void catchem_config_find_fengsha_static_file(void* core_ptr, char* val_out, int max_len) {
    if (core_ptr == nullptr) {
        copy_string_to_buffer("", val_out, max_len);
        return;
    }
    auto* core = static_cast<catchem::Core*>(core_ptr);
    std::string result = core->get_config_manager()->find_process_file_setting("fengsha");
    copy_string_to_buffer(result, val_out, max_len);
}

int catchem_config_get_yaml_list_count(void* core_ptr, const char* yaml_path) {
    if (core_ptr == nullptr || yaml_path == nullptr)
        return 0;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return static_cast<int>(core->get_config_manager()->get_string_list(yaml_path).size());
}

void catchem_config_get_yaml_list_at(void* core_ptr, const char* yaml_path, int index, char* val_out, int max_len) {
    if (core_ptr == nullptr || yaml_path == nullptr || index < 0) {
        copy_string_to_buffer("", val_out, max_len);
        return;
    }
    auto* core = static_cast<catchem::Core*>(core_ptr);
    auto list = core->get_config_manager()->get_string_list(yaml_path);
    if (index < static_cast<int>(list.size())) {
        copy_string_to_buffer(list[index], val_out, max_len);
    } else {
        copy_string_to_buffer("", val_out, max_len);
    }
}
// =========================================================================
// Species Metadata and Property Query C-API
// =========================================================================
#define CATCHEM_SPECIES_DOUBLE_PROPERTY(api, member, fallback)                                                         \
    double catchem_state_get_species_##api(void* state_ptr, int index) {                                               \
        try {                                                                                                          \
            auto* state = static_cast<catchem::StateManager*>(state_ptr);                                              \
            const int idx_0 = index - 1;                                                                               \
            if (idx_0 >= 0 && idx_0 < static_cast<int>(state->chemistry().species_list.size()))                        \
                return state->chemistry().species_list[idx_0].member;                                                  \
        } catch (...) {                                                                                                \
        }                                                                                                              \
        return fallback;                                                                                               \
    }
#define CATCHEM_SPECIES_BOOL_PROPERTY(api, member)                                                                     \
    int catchem_state_get_species_##api(void* state_ptr, int index) {                                                  \
        try {                                                                                                          \
            auto* state = static_cast<catchem::StateManager*>(state_ptr);                                              \
            const int idx_0 = index - 1;                                                                               \
            if (idx_0 >= 0 && idx_0 < static_cast<int>(state->chemistry().species_list.size()))                        \
                return state->chemistry().species_list[idx_0].member ? 1 : 0;                                          \
        } catch (...) {                                                                                                \
        }                                                                                                              \
        return 0;                                                                                                      \
    }
#define CATCHEM_SPECIES_LEGACY_BOOL_PROPERTY(api, member)                                                              \
    int catchem_state_is_species_##api(void* state_ptr, int index) {                                                   \
        try {                                                                                                          \
            auto* state = static_cast<catchem::StateManager*>(state_ptr);                                              \
            const int idx_0 = index - 1;                                                                               \
            if (idx_0 >= 0 && idx_0 < static_cast<int>(state->chemistry().species_list.size()))                        \
                return state->chemistry().species_list[idx_0].member ? 1 : 0;                                          \
        } catch (...) {                                                                                                \
        }                                                                                                              \
        return 0;                                                                                                      \
    }
#include "catchem_species_properties.def"
#undef CATCHEM_SPECIES_LEGACY_BOOL_PROPERTY
#undef CATCHEM_SPECIES_BOOL_PROPERTY
#undef CATCHEM_SPECIES_DOUBLE_PROPERTY

void catchem_state_get_species_wd_rainouteff(void* state_ptr, int index, double* eff_out) {
    try {
        auto* state = static_cast<catchem::StateManager*>(state_ptr);
        int idx_0 = index - 1;
        if (idx_0 >= 0 && idx_0 < static_cast<int>(state->chemistry().species_list.size())) {
            auto& eff = state->chemistry().species_list[idx_0].wd_rainouteff;
            for (size_t i = 0; i < 3; ++i) {
                eff_out[i] = i < eff.size() ? eff[i] : 0.0;
            }
        } else {
            eff_out[0] = eff_out[1] = eff_out[2] = 0.0;
        }
    } catch (...) {
        eff_out[0] = eff_out[1] = eff_out[2] = 0.0;
    }
}
}
