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
#include <sstream>
#include <vector>

using namespace catchem::api_internal;

extern "C" {
void catchem_state_bind_1d(void* state_ptr, const char* name, double* ptr) {
    catchem::AdmissionLease admission;
    if (admit_handle(state_ptr, catchem::HandleType::State, "state_bind_1d", admission) != CATCHEM_SUCCESS)
        return;
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    (void)catchem_state_bind_1d_checked(state_ptr, name, ptr, state->column_count());
}

int catchem_state_bind_1d_checked(void* state_ptr, const char* name, double* ptr, int dim1) {
    catchem::clear_boundary_error();
    if (!name || !ptr)
        return fail(catchem::BoundaryStatus::NullArgument, "state_bind_1d", name, "null name or data");
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, "state_bind_1d", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    if (dim1 != state->column_count())
        return fail(catchem::BoundaryStatus::ExtentMismatch, "state_bind_1d", name, "extent mismatch");
    try {
        state->bind_field_1d(name, ptr);
        return CATCHEM_SUCCESS;
    } catch (const std::exception& error) {
        return fail(catchem::BoundaryStatus::InternalError, "state_bind_1d", name, error.what());
    }
}

void catchem_state_bind_2d(void* state_ptr, const char* name, double* ptr) {
    catchem::AdmissionLease admission;
    if (admit_handle(state_ptr, catchem::HandleType::State, "state_bind_2d", admission) != CATCHEM_SUCCESS)
        return;
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    (void)catchem_state_bind_2d_checked(state_ptr, name, ptr, state->column_count(), 1);
}

int catchem_state_bind_2d_checked(void* state_ptr, const char* name, double* ptr, int dim1, int dim2) {
    catchem::clear_boundary_error();
    if (!name || !ptr)
        return fail(catchem::BoundaryStatus::NullArgument, "state_bind_2d", name, "null name or data");
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, "state_bind_2d", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    if (dim1 != state->column_count() || dim2 != 1)
        return fail(catchem::BoundaryStatus::ExtentMismatch, "state_bind_2d", name, "extent mismatch");
    try {
        state->bind_field_2d(name, ptr);
        return CATCHEM_SUCCESS;
    } catch (const std::exception& error) {
        return fail(catchem::BoundaryStatus::InternalError, "state_bind_2d", name, error.what());
    }
}

void catchem_state_bind_3d(void* state_ptr, const char* name, double* ptr) {
    catchem::AdmissionLease admission;
    if (admit_handle(state_ptr, catchem::HandleType::State, "state_bind_3d", admission) != CATCHEM_SUCCESS)
        return;
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    (void)catchem_state_bind_3d_checked(state_ptr, name, ptr, state->column_count(), state->level_count(), 1);
}

int catchem_state_bind_3d_checked(void* state_ptr, const char* name, double* ptr, int dim1, int dim2, int dim3) {
    catchem::clear_boundary_error();
    if (!name || !ptr)
        return fail(catchem::BoundaryStatus::NullArgument, "state_bind_3d", name, "null name or data");
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, "state_bind_3d", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    if (dim1 != state->column_count() || dim2 != state->level_count() || dim3 != 1)
        return fail(catchem::BoundaryStatus::ExtentMismatch, "state_bind_3d", name, "extent mismatch");
    try {
        state->bind_field_3d(name, ptr);
        return CATCHEM_SUCCESS;
    } catch (const std::exception& error) {
        return fail(catchem::BoundaryStatus::InternalError, "state_bind_3d", name, error.what());
    }
}

void catchem_state_bind_met_2d(void* state_ptr, const char* name, double* ptr) {
    catchem::AdmissionLease admission;
    if (admit_handle(state_ptr, catchem::HandleType::State, "state_bind_met_2d", admission) != CATCHEM_SUCCESS)
        return;
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    state->bind_met_field_2d(name, ptr);
}

void catchem_state_bind_met_3d(void* state_ptr, const char* name, double* ptr) {
    catchem::AdmissionLease admission;
    if (admit_handle(state_ptr, catchem::HandleType::State, "state_bind_met_3d", admission) != CATCHEM_SUCCESS)
        return;
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    state->bind_met_field_3d(name, ptr);
}

int catchem_state_begin_import_generation(void* state_ptr) {
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, "state_begin_import_generation", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    try {
        static_cast<catchem::StateManager*>(state_ptr)->begin_import_generation();
        return CATCHEM_SUCCESS;
    } catch (...) {
        return CATCHEM_INTERNAL_ERROR;
    }
}

int catchem_state_set_physical_validation_policy_checked(void* state_ptr, int policy) {
    catchem::AdmissionLease admission;
    const int status =
        admit_handle(state_ptr, catchem::HandleType::State, "state_set_physical_validation_policy", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    if (policy < 0 || policy > static_cast<int>(catchem::PhysicalValidationPolicy::CountAndContinue))
        return fail(catchem::BoundaryStatus::InvalidState, "state_set_physical_validation_policy", "policy",
                    "policy is outside the supported enumeration");
    static_cast<catchem::StateManager*>(state_ptr)->set_validation_policy(
        static_cast<catchem::PhysicalValidationPolicy>(policy));
    return CATCHEM_SUCCESS;
}

int catchem_state_get_physical_validation_report_checked(void* state_ptr, int* issue_count, char* detail,
                                                         int detail_length) {
    if (issue_count)
        *issue_count = 0;
    if (detail && detail_length > 0)
        detail[0] = '\0';
    if (!issue_count || !detail || detail_length <= 0)
        return fail(catchem::BoundaryStatus::NullArgument, "state_get_physical_validation_report", "output",
                    "issue count and a positive-length detail buffer are required");
    catchem::AdmissionLease admission;
    const int status =
        admit_handle(state_ptr, catchem::HandleType::State, "state_get_physical_validation_report", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    const auto& report = static_cast<catchem::StateManager*>(state_ptr)->validation_report();
    *issue_count = static_cast<int>(report.issue_count());
    copy_string_to_buffer(report.format(), detail, detail_length);
    return CATCHEM_SUCCESS;
}

int catchem_state_bind_met_2d_checked(void* state_ptr, const char* name, double* ptr, int dim1, int dim2) {
    if (!name || !ptr)
        return CATCHEM_NULL_ARGUMENT;
    catchem::AdmissionLease admission;
    const int handle_status = admit_handle(state_ptr, catchem::HandleType::State, "state_bind_met_2d", admission);
    if (handle_status != CATCHEM_SUCCESS)
        return handle_status;
    try {
        auto* state = static_cast<catchem::StateManager*>(state_ptr);
        if (dim1 != state->column_count() || dim2 != 1)
            return CATCHEM_EXTENT_MISMATCH;
        state->bind_met_field_2d(name, ptr);
        return CATCHEM_SUCCESS;
    } catch (...) {
        return CATCHEM_INTERNAL_ERROR;
    }
}

int catchem_state_bind_met_3d_checked(void* state_ptr, const char* name, double* ptr, int dim1, int dim2, int dim3) {
    if (!name || !ptr)
        return fail(catchem::BoundaryStatus::NullArgument, "state_bind_met_3d", "argument",
                    "field name and data pointer are required");
    catchem::AdmissionLease admission;
    const int handle_status = admit_handle(state_ptr, catchem::HandleType::State, "state_bind_met_3d", admission);
    if (handle_status != CATCHEM_SUCCESS)
        return handle_status;
    try {
        auto* state = static_cast<catchem::StateManager*>(state_ptr);
        std::string key(name);
        std::transform(key.begin(), key.end(), key.begin(), [](unsigned char c) { return std::toupper(c); });
        const int expected_levels = (key == "PEDGE" || key == "Z") ? state->level_count() + 1 : state->level_count();
        if (dim1 != state->column_count() || dim2 != expected_levels || dim3 != 1)
            return fail(catchem::BoundaryStatus::ExtentMismatch, "state_bind_met_3d", name,
                        "field extents do not match the state contract");
        state->bind_met_field_3d(name, ptr);
        return CATCHEM_SUCCESS;
    } catch (const std::exception& error) {
        return fail(catchem::BoundaryStatus::InternalError, "state_bind_met_3d", name, error.what());
    } catch (...) {
        return fail(catchem::BoundaryStatus::InternalError, "state_bind_met_3d", name, "unknown exception");
    }
}

int catchem_state_bind_met_3d_axis_checked(void* state_ptr, const char* name, double* ptr, int dim1, int dim2, int dim3,
                                           int semantic_axis) {
    if (!name || !ptr)
        return CATCHEM_NULL_ARGUMENT;
    catchem::AdmissionLease admission;
    const int handle_status = admit_handle(state_ptr, catchem::HandleType::State, "state_bind_met_3d_axis", admission);
    if (handle_status != CATCHEM_SUCCESS)
        return handle_status;
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    catchem::SemanticAxis axis;
    int expected_vertical = dim2;
    switch (semantic_axis) {
    case 0:
        axis = catchem::SemanticAxis::Level;
        expected_vertical = state->level_count();
        break;
    case 1:
        axis = catchem::SemanticAxis::Interface;
        expected_vertical = state->level_count() + 1;
        break;
    case 2:
        axis = catchem::SemanticAxis::SoilLayer;
        break;
    default:
        return CATCHEM_INVALID_STATE;
    }
    if (dim1 != state->column_count() || dim2 <= 0 || dim2 != expected_vertical || dim3 != 1)
        return CATCHEM_EXTENT_MISMATCH;
    try {
        state->bind_met_field_3d_contract(name, ptr, dim2, axis);
        return CATCHEM_SUCCESS;
    } catch (...) {
        return CATCHEM_INTERNAL_ERROR;
    }
}

void catchem_state_bind_unified_chemistry(void* state_ptr, double* ptr) {
    catchem::AdmissionLease admission;
    if (admit_handle(state_ptr, catchem::HandleType::State, "state_bind_chemistry", admission) != CATCHEM_SUCCESS)
        return;
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    (void)catchem_state_bind_unified_chemistry_checked(state_ptr, ptr, state->column_count(), state->level_count(),
                                                       state->species_count());
}

int catchem_state_bind_unified_chemistry_checked(void* state_ptr, double* ptr, int dim1, int dim2, int dim3) {
    if (!ptr)
        return CATCHEM_NULL_ARGUMENT;
    catchem::AdmissionLease admission;
    const int handle_status = admit_handle(state_ptr, catchem::HandleType::State, "state_bind_chemistry", admission);
    if (handle_status != CATCHEM_SUCCESS)
        return handle_status;
    try {
        auto* state = static_cast<catchem::StateManager*>(state_ptr);
        if (dim1 != state->column_count() || dim2 != state->level_count() || dim3 != state->species_count())
            return CATCHEM_EXTENT_MISMATCH;
        state->bind_unified_chemistry(ptr);
        return CATCHEM_SUCCESS;
    } catch (...) {
        return CATCHEM_INTERNAL_ERROR;
    }
}

int catchem_state_mark_chem_host_modified(void* state_ptr) {
    catchem::AdmissionLease admission;
    const int handle_status =
        admit_handle(state_ptr, catchem::HandleType::State, "state_mark_chem_host_modified", admission);
    if (handle_status != CATCHEM_SUCCESS)
        return handle_status;
    try {
        auto* state = static_cast<catchem::StateManager*>(state_ptr);
        if (!state->chemistry().conc)
            return CATCHEM_INVALID_STATE;
        state->chemistry().conc->mark_host_modified();
        return CATCHEM_SUCCESS;
    } catch (...) {
        return CATCHEM_INTERNAL_ERROR;
    }
}

void catchem_state_set_time(void* state_ptr, int yr, int mo, int dy, int hr, int mn, int sc, int doy, double tstep) {
    (void)catchem_state_set_time_checked(state_ptr, yr, mo, dy, hr, mn, sc, doy, tstep);
}

int catchem_state_set_time_checked(void* state_ptr, int yr, int mo, int dy, int hr, int mn, int sc, int doy,
                                   double tstep) {
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, "state_set_time", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    try {
        auto* state = static_cast<catchem::StateManager*>(state_ptr);
        state->clock().year = yr;
        state->clock().month = mo;
        state->clock().day = dy;
        state->clock().hour = hr;
        state->clock().minute = mn;
        state->clock().second = sc;
        state->clock().doy = doy;
        state->clock().timestep = tstep;
        return CATCHEM_SUCCESS;
    } catch (const std::exception& error) {
        return fail(catchem::BoundaryStatus::InternalError, "state_set_time", "time", error.what());
    } catch (...) {
        return fail(catchem::BoundaryStatus::InternalError, "state_set_time", "time", "unknown exception");
    }
}

void catchem_state_sync_to_device(void* state_ptr) {
    (void)catchem_state_sync_to_device_checked(state_ptr);
}

int catchem_state_sync_to_device_checked(void* state_ptr) {
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, "state_sync_to_device", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    try {
        static_cast<catchem::StateManager*>(state_ptr)->sync_to_device();
        return CATCHEM_SUCCESS;
    } catch (const std::exception& error) {
        return fail(catchem::BoundaryStatus::InternalError, "state_sync_to_device", "state", error.what());
    } catch (...) {
        return fail(catchem::BoundaryStatus::InternalError, "state_sync_to_device", "state", "unknown exception");
    }
}

void catchem_state_sync_to_host(void* state_ptr) {
    (void)catchem_state_sync_to_host_checked(state_ptr);
}

int catchem_state_sync_to_host_checked(void* state_ptr) {
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, "state_sync_to_host", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    try {
        static_cast<catchem::StateManager*>(state_ptr)->sync_to_host();
        return CATCHEM_SUCCESS;
    } catch (const std::exception& error) {
        return fail(catchem::BoundaryStatus::InternalError, "state_sync_to_host", "state", error.what());
    } catch (...) {
        return fail(catchem::BoundaryStatus::InternalError, "state_sync_to_host", "state", "unknown exception");
    }
}

double* catchem_state_get_pointer_1d(void* state_ptr, const char* name) {
    catchem::AdmissionLease admission;
    if (admit_handle(state_ptr, catchem::HandleType::State, "state_get_pointer_1d", admission) != CATCHEM_SUCCESS)
        return nullptr;
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    return state->get_host_pointer_1d(name);
}

double* catchem_state_get_pointer_2d(void* state_ptr, const char* name) {
    catchem::AdmissionLease admission;
    if (admit_handle(state_ptr, catchem::HandleType::State, "state_get_pointer_2d", admission) != CATCHEM_SUCCESS)
        return nullptr;
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    return state->get_host_pointer_2d(name);
}

double* catchem_state_get_pointer_3d(void* state_ptr, const char* name) {
    void* pointer = nullptr;
    (void)catchem_state_get_pointer_3d_checked(state_ptr, name, &pointer);
    return static_cast<double*>(pointer);
}

int catchem_state_get_pointer_3d_checked(void* state_ptr, const char* name, void** ptr_out) {
    if (ptr_out)
        *ptr_out = nullptr;
    if (!name || !ptr_out)
        return fail(catchem::BoundaryStatus::NullArgument, "state_get_pointer_3d", "argument",
                    "field name and pointer output are required");
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, "state_get_pointer_3d", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    try {
        *ptr_out = static_cast<catchem::StateManager*>(state_ptr)->get_host_pointer_3d(name);
        if (!*ptr_out)
            return fail(catchem::BoundaryStatus::MissingField, "state_get_pointer_3d", name,
                        "field is not registered or current");
        return CATCHEM_SUCCESS;
    } catch (const std::exception& error) {
        return fail(catchem::BoundaryStatus::InternalError, "state_get_pointer_3d", name, error.what());
    } catch (...) {
        return fail(catchem::BoundaryStatus::InternalError, "state_get_pointer_3d", name, "unknown exception");
    }
}

double* catchem_state_get_species_conc_pointer(void* state_ptr, int species_index) {
    catchem::AdmissionLease admission;
    if (admit_handle(state_ptr, catchem::HandleType::State, "state_get_species_concentration", admission) !=
        CATCHEM_SUCCESS)
        return nullptr;
    auto* state = static_cast<catchem::StateManager*>(state_ptr);
    int idx_0 = species_index - 1;
    if (state->chemistry().conc && state->chemistry().conc->host_data() && idx_0 >= 0 &&
        idx_0 < state->species_count()) {
        return state->chemistry().conc->host_data() +
               static_cast<size_t>(idx_0) * state->column_count() * state->level_count();
    }
    return nullptr;
}

int catchem_state_get_species_conc_pointer_checked(void* state_ptr, int species_index, int dim1, int dim2,
                                                   double** ptr_out) {
    if (!ptr_out)
        return CATCHEM_NULL_ARGUMENT;
    *ptr_out = nullptr;
    catchem::AdmissionLease admission;
    const int handle_status =
        admit_handle(state_ptr, catchem::HandleType::State, "state_get_species_concentration", admission);
    if (handle_status != CATCHEM_SUCCESS)
        return handle_status;
    try {
        auto* state = static_cast<catchem::StateManager*>(state_ptr);
        if (dim1 != state->column_count() || dim2 != state->level_count())
            return CATCHEM_EXTENT_MISMATCH;
        if (species_index < 1 || species_index > state->species_count())
            return CATCHEM_INVALID_INDEX;
        if (!state->chemistry().conc)
            return CATCHEM_INVALID_STATE;
        // This is a read boundary used by NUOPC export and diagnostics.  Make
        // the host mirror current before exposing its raw pointer; otherwise a
        // device-resident process can appear to have produced no concentration
        // change at the coupling boundary.
        const double* concentration = state->chemistry().conc->host_read();
        if (!concentration)
            return CATCHEM_INVALID_STATE;
        *ptr_out = const_cast<double*>(concentration) + static_cast<std::size_t>(species_index - 1) * dim1 * dim2;
        return CATCHEM_SUCCESS;
    } catch (...) {
        return CATCHEM_INTERNAL_ERROR;
    }
}

void catchem_state_derive_bxheight(void* state_ptr) {
    (void)catchem_state_derive_bxheight_checked(state_ptr);
}

int catchem_state_derive_bxheight_checked(void* state_ptr) {
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, "state_derive_bxheight", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    try {
        static_cast<catchem::StateManager*>(state_ptr)->derive_bxheight();
        return CATCHEM_SUCCESS;
    } catch (const std::domain_error& error) {
        return fail(catchem::BoundaryStatus::PhysicalValidationFailure, "state_derive_bxheight", "inputs",
                    error.what());
    } catch (const std::exception& error) {
        return fail(catchem::BoundaryStatus::InternalError, "state_derive_bxheight", "state", error.what());
    } catch (...) {
        return fail(catchem::BoundaryStatus::InternalError, "state_derive_bxheight", "state", "unknown exception");
    }
}

void catchem_state_derive_airden_dry(void* state_ptr) {
    (void)catchem_state_derive_airden_dry_checked(state_ptr);
}

int catchem_state_derive_airden_dry_checked(void* state_ptr) {
    catchem::AdmissionLease admission;
    const int status = admit_handle(state_ptr, catchem::HandleType::State, "state_derive_airden_dry", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    try {
        static_cast<catchem::StateManager*>(state_ptr)->derive_airden_dry();
        return CATCHEM_SUCCESS;
    } catch (const std::domain_error& error) {
        return fail(catchem::BoundaryStatus::PhysicalValidationFailure, "state_derive_airden_dry", "inputs",
                    error.what());
    } catch (const std::exception& error) {
        return fail(catchem::BoundaryStatus::InternalError, "state_derive_airden_dry", "state", error.what());
    } catch (...) {
        return fail(catchem::BoundaryStatus::InternalError, "state_derive_airden_dry", "state", "unknown exception");
    }
}

int catchem_state_get_nx(void* state_ptr) {
    if (!state_ptr)
        return 0;
    return static_cast<catchem::StateManager*>(state_ptr)->column_count();
}

int catchem_state_get_ny(void* state_ptr) {
    if (!state_ptr)
        return 0;
    return 1; // 1D column arrays natively, but exposing ny=1 for 2D interfaces
}

int catchem_state_get_nz(void* state_ptr) {
    if (!state_ptr)
        return 0;
    return static_cast<catchem::StateManager*>(state_ptr)->level_count();
}
}
