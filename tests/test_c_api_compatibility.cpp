#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include <cstdlib>
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

int main() {
    void* legacy = catchem_core_create(2, 3, 1);
    void* checked = nullptr;
    CHECK(catchem_core_create_checked(2, 3, 1, &checked) == CATCHEM_SUCCESS);
    CHECK(legacy && checked);
    void* checked_state = nullptr;
    void* legacy_state = catchem_core_get_state_manager(legacy);
    CHECK(catchem_core_get_state_manager_checked(checked, &checked_state) == CATCHEM_SUCCESS);
    CHECK(legacy_state && checked_state);
    std::vector<double> legacy_temperature(6, 280.0), checked_temperature(6, 280.0);
    catchem_state_bind_3d(legacy_state, "T", legacy_temperature.data());
    CHECK(catchem_state_bind_3d_checked(checked_state, "T", checked_temperature.data(), 2, 3, 1) == CATCHEM_SUCCESS);
    CHECK(catchem_state_get_pointer_3d(legacy_state, "T") == legacy_temperature.data());
    CHECK(catchem_state_get_pointer_3d(checked_state, "T") == checked_temperature.data());
    catchem_core_destroy(legacy);
    CHECK(catchem_core_destroy_checked(checked) == CATCHEM_SUCCESS);

    // --- feature 013: axes + unpack-label C getters ---
    void* diag_core = nullptr;
    CHECK(catchem_core_create_checked(2, 3, 1, &diag_core) == CATCHEM_SUCCESS);
    auto diag_manager = static_cast<catchem::Core*>(diag_core)->get_diagnostic_manager();
    diag_manager->register_field_contract(
        "velocity_per_species", "settling velocity", "m/s", catchem::DiagType::FIELD_3D, std::vector<int>{2, 3, 2},
        catchem::DiagnosticPolicy::Instantaneous, 0.0,
        std::vector<catchem::SemanticAxis>{catchem::SemanticAxis::Column, catchem::SemanticAxis::Level,
                                           catchem::SemanticAxis::Species},
        std::vector<std::string>{"SO4", "BC2"});
    diag_manager->register_field_contract(
        "total", "single value", "kg", catchem::DiagType::FIELD_2D, std::vector<int>{2, 1},
        catchem::DiagnosticPolicy::Instantaneous, 0.0,
        std::vector<catchem::SemanticAxis>{catchem::SemanticAxis::Column, catchem::SemanticAxis::Singleton});

    int axes_out[4] = {-1, -1, -1, -1};
    CHECK(catchem_diag_get_axes_checked(diag_core, "velocity_per_species", axes_out, 4) == CATCHEM_SUCCESS);
    CHECK(axes_out[0] == 0 && axes_out[1] == 1 && axes_out[2] == 4 && axes_out[3] == 0);
    // Output shorter than the rank is rejected rather than silently truncated.
    CHECK(catchem_diag_get_axes_checked(diag_core, "velocity_per_species", axes_out, 2) != CATCHEM_SUCCESS);

    char label[16];
    CHECK(catchem_diag_get_unpack_label_at_checked(diag_core, "velocity_per_species", 0, label, sizeof(label)) ==
          CATCHEM_SUCCESS);
    CHECK(std::string(label) == "SO4");
    CHECK(catchem_diag_get_unpack_label_at_checked(diag_core, "velocity_per_species", 1, label, sizeof(label)) ==
          CATCHEM_SUCCESS);
    CHECK(std::string(label) == "BC2");
    // Out-of-range slot and a field with no packed dimension both fail loudly.
    CHECK(catchem_diag_get_unpack_label_at_checked(diag_core, "velocity_per_species", 2, label, sizeof(label)) !=
          CATCHEM_SUCCESS);
    CHECK(catchem_diag_get_unpack_label_at_checked(diag_core, "total", 0, label, sizeof(label)) != CATCHEM_SUCCESS);
    // Unknown field names are reported, not answered with an empty success.
    CHECK(catchem_diag_get_axes_checked(diag_core, "missing_field", axes_out, 4) != CATCHEM_SUCCESS);

    CHECK(catchem_core_destroy_checked(diag_core) == CATCHEM_SUCCESS);
    return 0;
}
