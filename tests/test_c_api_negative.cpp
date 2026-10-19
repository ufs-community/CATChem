#include "catchem_api.hpp"
#include <cassert>
#include <cstdint>
#include <string>
#include <vector>

int main() {
    void* core = nullptr;
    assert(catchem_core_create_checked(2, 3, 4, &core) == CATCHEM_SUCCESS);
    assert(core != nullptr);
    void* state = nullptr;
    assert(catchem_core_get_state_manager_checked(core, &state) == CATCHEM_SUCCESS);
    assert(state != nullptr);
    int required_count = -1;
    assert(catchem_core_get_required_host_field_count_checked(core, &required_count) == CATCHEM_SUCCESS);
    assert(required_count == 0);
    char required_name[8] = {'x'};
    assert(catchem_core_get_required_host_field_name_checked(core, 0, required_name, sizeof(required_name)) ==
           CATCHEM_INVALID_INDEX);
    assert(required_name[0] == '\0');
    assert(catchem_core_get_required_host_field_count_checked(state, &required_count) == CATCHEM_WRONG_HANDLE_TYPE);
    assert(catchem_state_set_physical_validation_policy_checked(state, 0) == CATCHEM_SUCCESS);
    assert(catchem_state_set_physical_validation_policy_checked(state, 99) == CATCHEM_INVALID_STATE);
    int physical_issues = -1;
    char physical_detail[8] = {'x'};
    assert(catchem_state_get_physical_validation_report_checked(state, &physical_issues, physical_detail,
                                                                sizeof(physical_detail)) == CATCHEM_SUCCESS);
    assert(physical_issues == 0 && physical_detail[0] == '\0');

    std::vector<double> one(2, 0.0);
    std::vector<double> two(2, 0.0);
    std::vector<double> three(6, 0.0);
    std::vector<double> chemistry(24, 0.0);
    void* random = reinterpret_cast<void*>(static_cast<std::uintptr_t>(0x12345));
    assert(catchem_core_create_from_config_checked(nullptr, &random) == CATCHEM_NULL_ARGUMENT);
    assert(random == nullptr);
    random = reinterpret_cast<void*>(static_cast<std::uintptr_t>(0x12345));
    assert(catchem_core_create_from_config_with_grid_checked(nullptr, 2, 3, &random) == CATCHEM_NULL_ARGUMENT);
    assert(random == nullptr);
    random = reinterpret_cast<void*>(static_cast<std::uintptr_t>(0x12345));
    assert(catchem_core_create_from_config_checked("does-not-exist.yml", &random) == CATCHEM_INVALID_STATE);
    assert(random == nullptr);
    int output = 77;
    assert(catchem_core_get_num_processes_checked(state, &output) == CATCHEM_WRONG_HANDLE_TYPE);
    assert(output == 0);
    assert(catchem_core_get_num_processes_checked(core, nullptr) == CATCHEM_NULL_ARGUMENT);
    assert(catchem_state_sync_to_device_checked(core) == CATCHEM_WRONG_HANDLE_TYPE);
    assert(catchem_state_sync_to_host_checked(nullptr) == CATCHEM_NULL_ARGUMENT);
    assert(catchem_state_set_time_checked(core, 2026, 1, 1, 0, 0, 0, 1, 60.0) == CATCHEM_WRONG_HANDLE_TYPE);
    output = 77;
    assert(catchem_diag_get_count_checked(state, &output) == CATCHEM_WRONG_HANDLE_TYPE);
    assert(output == 0);
    output = 77;
    assert(catchem_diag_get_rank_checked(core, "missing", &output) == CATCHEM_MISSING_FIELD);
    assert(output == 0);
    int dims[3] = {7, 7, 7};
    assert(catchem_diag_get_dims_checked(core, "missing", dims, 3) == CATCHEM_MISSING_FIELD);
    assert(dims[0] == 0 && dims[1] == 0 && dims[2] == 0);
    char missing_name[8] = {'x'};
    assert(catchem_diag_get_name_at_checked(core, 0, missing_name, sizeof(missing_name)) == CATCHEM_INVALID_INDEX);
    assert(missing_name[0] == '\0');
    assert(catchem_diag_register_checked(core, "bad", "bad", "1", 1, 2, 1, 1) == CATCHEM_RANK_MISMATCH);
    int contract_dims[2] = {2, 1};
    int contract_axes[2] = {0, 4};
    assert(catchem_diag_register_contract_checked(core, "bad", "bad", "1", 2, nullptr, contract_axes, 0, 0.0) ==
           CATCHEM_NULL_ARGUMENT);
    void* diagnostic_pointer = random;
    assert(catchem_diag_get_pointer_checked(core, "missing", 2, contract_dims, &diagnostic_pointer) ==
           CATCHEM_MISSING_FIELD);
    assert(diagnostic_pointer == nullptr);
    output = 77;
    assert(catchem_state_get_species_count_checked(core, &output) == CATCHEM_WRONG_HANDLE_TYPE);
    assert(output == 0);
    double molecular_weight = 77.0;
    assert(catchem_state_get_species_mw_checked(state, 1, &molecular_weight) == CATCHEM_INVALID_INDEX);
    assert(molecular_weight == 0.0);
    output = 77;
    assert(catchem_state_is_species_gas_checked(state, 1, &output) == CATCHEM_INVALID_INDEX);
    assert(output == 0);
    char species_name[8] = {'x'};
    assert(catchem_state_get_species_name_at_checked(state, 1, species_name, sizeof(species_name)) ==
           CATCHEM_INVALID_INDEX);
    assert(species_name[0] == '\0');
    output = 77;
    assert(catchem_state_get_species_index_checked(state, "missing", &output) == CATCHEM_INVALID_INDEX);
    assert(output == 0);
    output = 77;
    assert(catchem_state_get_species_is_advected_checked(state, 1, &output) == CATCHEM_INVALID_INDEX);
    assert(output == 0);
    output = 77;
    assert(catchem_state_is_species_aerosol_checked(state, 1, &output) == CATCHEM_INVALID_INDEX);
    assert(output == 0);
    void* field_pointer = random;
    assert(catchem_state_get_pointer_3d_checked(state, "missing", &field_pointer) == CATCHEM_MISSING_FIELD);
    assert(field_pointer == nullptr);

    std::vector<double> invalid_temperature(6, -1.0);
    std::vector<double> valid_pressure(6, 90000.0);
    std::vector<double> valid_humidity(6, 0.01);
    assert(catchem_state_bind_met_3d_checked(state, "T", invalid_temperature.data(), 2, 3, 1) == CATCHEM_SUCCESS);
    assert(catchem_state_bind_met_3d_checked(state, "PMID", valid_pressure.data(), 2, 3, 1) == CATCHEM_SUCCESS);
    assert(catchem_state_bind_met_3d_checked(state, "QV", valid_humidity.data(), 2, 3, 1) == CATCHEM_SUCCESS);
    assert(catchem_state_derive_airden_dry_checked(state) == CATCHEM_PHYSICAL_VALIDATION_FAILURE);
    char derivation_error[256] = {};
    assert(catchem_get_last_error(derivation_error, sizeof(derivation_error)) == CATCHEM_PHYSICAL_VALIDATION_FAILURE);
    assert(std::string(derivation_error).find("Physical validation") != std::string::npos);
    assert(catchem_state_derive_bxheight_checked(core) == CATCHEM_WRONG_HANDLE_TYPE);
    assert(catchem_state_bind_met_2d_checked(state, "LAT", one.data(), 3, 1) == CATCHEM_EXTENT_MISMATCH);
    assert(catchem_state_bind_met_3d_axis_checked(state, "soil", three.data(), 2, 3, 1, 99) == CATCHEM_INVALID_STATE);
    double* concentration_pointer = reinterpret_cast<double*>(static_cast<std::uintptr_t>(0x12345));
    assert(catchem_state_get_species_conc_pointer_checked(state, 9, 2, 3, &concentration_pointer) ==
           CATCHEM_INVALID_INDEX);
    assert(concentration_pointer == nullptr);

    int misuse_cases = 0;
    for (int repetition = 0; repetition < 10; ++repetition) {
        assert(catchem_core_create_checked(0, 3, 4, &random) == CATCHEM_EXTENT_MISMATCH);
        ++misuse_cases;
        assert(catchem_core_get_state_manager_checked(state, &random) == CATCHEM_WRONG_HANDLE_TYPE);
        ++misuse_cases;
        assert(catchem_state_bind_1d_checked(core, "x", one.data(), 2) == CATCHEM_WRONG_HANDLE_TYPE);
        ++misuse_cases;
        assert(catchem_state_bind_2d_checked(state, "x", two.data(), 3, 1) == CATCHEM_EXTENT_MISMATCH);
        ++misuse_cases;
        assert(catchem_state_bind_3d_checked(state, "x", three.data(), 2, 4, 1) == CATCHEM_EXTENT_MISMATCH);
        ++misuse_cases;
        assert(catchem_state_bind_met_3d_checked(state, nullptr, three.data(), 2, 3, 1) == CATCHEM_NULL_ARGUMENT);
        ++misuse_cases;
        assert(catchem_state_bind_met_3d_checked(state, "T", nullptr, 2, 3, 1) == CATCHEM_NULL_ARGUMENT);
        ++misuse_cases;
        assert(catchem_state_bind_unified_chemistry_checked(state, chemistry.data(), 2, 3, 5) ==
               CATCHEM_EXTENT_MISMATCH);
        ++misuse_cases;
    }
    assert(misuse_cases >= 50);

    void* stale_state = state;
    assert(catchem_core_destroy_checked(core) == CATCHEM_SUCCESS);
    assert(catchem_core_destroy_checked(core) == CATCHEM_INVALID_HANDLE);
    assert(catchem_state_begin_import_generation(stale_state) == CATCHEM_INVALID_HANDLE);
    assert(catchem_state_begin_import_generation(nullptr) == CATCHEM_NULL_ARGUMENT);
    assert(catchem_state_begin_import_generation(reinterpret_cast<void*>(static_cast<std::uintptr_t>(0x12345))) ==
           CATCHEM_INVALID_HANDLE);

    char error[256] = {};
    assert(catchem_get_last_error(error, sizeof(error)) == CATCHEM_INVALID_HANDLE);
    assert(std::string(error).find("handle") != std::string::npos);
    return 0;
}
