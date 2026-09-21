#include "catchem_api.hpp"
#include "catchem_test_config.hpp"

#include "catchem_core.hpp"
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

// Unlike assert(), CHECK always evaluates its expression: several calls below
// populate pointers consumed on later lines, and NDEBUG builds must not skip
// the side effects (or the checks themselves).
#define CHECK(expr)                                                                                                    \
    do {                                                                                                               \
        if (!(expr)) {                                                                                                 \
            std::cerr << "CHECK failed: " << #expr << " at " << __FILE__ << ':' << __LINE__ << std::endl;              \
            std::abort();                                                                                              \
        }                                                                                                              \
    } while (0)

// The standalone, NUOPC, and CCPP adapters all terminate in these checked
// operations. This fixture locks their common failure taxonomy and axis values.
int main() {
    const std::string conformance_config =
        std::string(catchem::test::TEST_DIR) + "/fixtures/host_conformance/CATChem_config.yml";
    catchem::Core configured(conformance_config);
    const auto& ordered_species = configured.get_state_manager()->chemistry().species_list;
    CHECK(ordered_species.size() == 3);
    CHECK(ordered_species[0].short_name == "unfamiliar_alpha");
    CHECK(ordered_species[1].short_name == "unfamiliar_beta");
    CHECK(ordered_species[2].short_name == "unfamiliar_gamma");

    void* core = nullptr;
    CHECK(catchem_core_create_checked(2, 3, 2, &core) == CATCHEM_SUCCESS);
    void* state = nullptr;
    CHECK(catchem_core_get_state_manager_checked(core, &state) == CATCHEM_SUCCESS);

    std::vector<double> layers(6), interfaces(8), soil(4), chemistry(12);
    CHECK(catchem_state_begin_import_generation(state) == CATCHEM_SUCCESS);
    CHECK(catchem_state_bind_met_3d_axis_checked(state, "T", layers.data(), 2, 3, 1, 0) == CATCHEM_SUCCESS);
    CHECK(catchem_state_bind_met_3d_axis_checked(state, "PEDGE", interfaces.data(), 2, 4, 1, 1) == CATCHEM_SUCCESS);
    CHECK(catchem_state_bind_met_3d_axis_checked(state, "soil_moisture", soil.data(), 2, 2, 1, 2) == CATCHEM_SUCCESS);
    CHECK(catchem_state_bind_unified_chemistry_checked(state, chemistry.data(), 2, 3, 2) == CATCHEM_SUCCESS);

    CHECK(catchem_state_bind_met_3d_axis_checked(state, "T", layers.data(), 2, 4, 1, 0) == CATCHEM_EXTENT_MISMATCH);
    CHECK(catchem_state_bind_met_3d_axis_checked(state, "T", layers.data(), 2, 3, 1, 99) == CATCHEM_INVALID_STATE);
    CHECK(catchem_state_bind_unified_chemistry_checked(state, chemistry.data(), 2, 3, 3) == CATCHEM_EXTENT_MISMATCH);

    const std::string mechanism_file = std::string(catchem::test::TEST_DIR) + "/fixtures/mechanisms/unfamiliar.yml";
    catchem_state_load_species_config(state, mechanism_file.c_str());
    CHECK(catchem_state_get_species_count(state) == 3);
    CHECK(catchem_state_get_species_index(state, "aerosol_c") == 3);
    CHECK(catchem_state_get_species_index(state, "XENON_A") == 1);

    // The same checked policy/report boundary is consumed by direct, CCPP, and
    // NUOPC hosts. Exercise its stable status, count, and detail representation.
    CHECK(catchem_state_set_physical_validation_policy_checked(state, 1) == CATCHEM_SUCCESS);
    std::vector<double> bad_temperature(6, -1.0);
    std::vector<double> pressure_mid(6, 90000.0);
    std::vector<double> humidity(6, 0.01);
    CHECK(catchem_state_bind_met_3d_checked(state, "T", bad_temperature.data(), 2, 3, 1) == CATCHEM_SUCCESS);
    CHECK(catchem_state_bind_met_3d_checked(state, "PMID", pressure_mid.data(), 2, 3, 1) == CATCHEM_SUCCESS);
    CHECK(catchem_state_bind_met_3d_checked(state, "QV", humidity.data(), 2, 3, 1) == CATCHEM_SUCCESS);
    catchem_state_derive_airden_dry(state);
    int issue_count = -1;
    std::array<char, 512> report{};
    CHECK(catchem_state_get_physical_validation_report_checked(state, &issue_count, report.data(),
                                                               static_cast<int>(report.size())) == CATCHEM_SUCCESS);
    CHECK(issue_count == 1);
    CHECK(std::strstr(report.data(), "T") != nullptr);
    issue_count = -1;
    report[0] = 'x';
    CHECK(catchem_state_get_physical_validation_report_checked(
              nullptr, &issue_count, report.data(), static_cast<int>(report.size())) == CATCHEM_NULL_ARGUMENT);
    CHECK(issue_count == 0 && report[0] == '\0');

    void* stale = state;
    CHECK(catchem_core_destroy_checked(core) == CATCHEM_SUCCESS);
    CHECK(catchem_state_begin_import_generation(stale) == CATCHEM_INVALID_HANDLE);
    return 0;
}
