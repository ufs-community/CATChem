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
void catchem_diag_register(void* core_ptr, const char* name, const char* desc, const char* units, int rank, int dim1,
                           int dim2, int dim3) {
    (void)catchem_diag_register_checked(core_ptr, name, desc, units, rank, dim1, dim2, dim3);
}

int catchem_diag_register_checked(void* core_ptr, const char* name, const char* desc, const char* units, int rank,
                                  int dim1, int dim2, int dim3) {
    if (!name || !desc || !units)
        return CATCHEM_NULL_ARGUMENT;
    catchem::AdmissionLease admission;
    const int handle_status = admit_handle(core_ptr, catchem::HandleType::Core, "diagnostic_register", admission);
    if (handle_status != CATCHEM_SUCCESS)
        return handle_status;
    try {
        auto* core = static_cast<catchem::Core*>(core_ptr);
        catchem::DiagType type;
        std::vector<int> dims;
        if (rank == 2) {
            type = catchem::DiagType::FIELD_2D;
            dims = {dim1, dim2};
        } else if (rank == 3) {
            type = catchem::DiagType::FIELD_3D;
            dims = {dim1, dim2, dim3};
        } else
            return CATCHEM_RANK_MISMATCH;
        core->get_diagnostic_manager()->register_field(name, desc, units, type, dims);
        return CATCHEM_SUCCESS;
    } catch (const std::invalid_argument&) {
        return CATCHEM_EXTENT_MISMATCH;
    } catch (...) {
        return CATCHEM_INTERNAL_ERROR;
    }
}

int catchem_diag_register_contract_checked(void* core_ptr, const char* name, const char* desc, const char* units,
                                           int rank, const int* dims, const int* axes, int policy, double reset_value) {
    if (!name || !desc || !units || !dims || !axes)
        return CATCHEM_NULL_ARGUMENT;
    catchem::AdmissionLease admission;
    const int handle_status =
        admit_handle(core_ptr, catchem::HandleType::Core, "diagnostic_register_contract", admission);
    if (handle_status != CATCHEM_SUCCESS)
        return handle_status;
    try {
        if (rank != 2 && rank != 3)
            return CATCHEM_RANK_MISMATCH;
        std::vector<int> dimensions(dims, dims + rank);
        std::vector<catchem::SemanticAxis> semantic_axes;
        for (int index = 0; index < rank; ++index) {
            if (axes[index] < 0 || axes[index] > static_cast<int>(catchem::SemanticAxis::Singleton))
                return CATCHEM_INVALID_STATE;
            semantic_axes.push_back(static_cast<catchem::SemanticAxis>(axes[index]));
        }
        if (policy < 0 || policy > static_cast<int>(catchem::DiagnosticPolicy::Persistent))
            return CATCHEM_INVALID_STATE;
        const auto type = rank == 2 ? catchem::DiagType::FIELD_2D : catchem::DiagType::FIELD_3D;
        // This C entry point carries no per-slot label argument, so it cannot
        // satisfy the INV-3..7 unpacking contract for a packed (Species or
        // Category) axis.  Register it without the strict label checks; the
        // axes-driven writer still fails loud (FR-008) at output time if such
        // a field ever needs to unpack.  Process diagnostics that DO unpack
        // register through the C++ DiagnosticManager::register_field_contract
        // API with labels, which keeps the strict validation.
        static_cast<catchem::Core*>(core_ptr)->get_diagnostic_manager()->register_field_contract(
            name, desc, units, type, dimensions, static_cast<catchem::DiagnosticPolicy>(policy), reset_value,
            semantic_axes, {}, false);
        return CATCHEM_SUCCESS;
    } catch (const std::invalid_argument&) {
        return CATCHEM_EXTENT_MISMATCH;
    } catch (...) {
        return CATCHEM_INTERNAL_ERROR;
    }
}

int catchem_diag_get_contract(void* core_ptr, const char* name, int* generation, int* availability, int* latest_writer,
                              int* policy) {
    if (!name || !generation || !availability || !latest_writer || !policy)
        return CATCHEM_NULL_ARGUMENT;
    catchem::AdmissionLease admission;
    const int handle_status = admit_handle(core_ptr, catchem::HandleType::Core, "diagnostic_get_contract", admission);
    if (handle_status != CATCHEM_SUCCESS)
        return handle_status;
    try {
        auto manager = static_cast<catchem::Core*>(core_ptr)->get_diagnostic_manager();
        if (!manager->has_field(name))
            return CATCHEM_MISSING_FIELD;
        const auto field = manager->get_field(name);
        *generation = static_cast<int>(field->generation);
        *availability = static_cast<int>(field->availability);
        *latest_writer = static_cast<int>(field->latest_writer);
        *policy = static_cast<int>(field->reset_policy);
        return CATCHEM_SUCCESS;
    } catch (...) {
        return CATCHEM_INTERNAL_ERROR;
    }
}

void* catchem_diag_get_pointer(void* core_ptr, const char* name) {
    catchem::AdmissionLease admission;
    if (admit_handle(core_ptr, catchem::HandleType::Core, "diagnostic_get_pointer", admission) != CATCHEM_SUCCESS ||
        !name)
        return nullptr;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    return core->get_diagnostic_manager()->get_host_pointer(name);
}

int catchem_diag_get_rank(void* core_ptr, const char* name) {
    int rank = 0;
    (void)catchem_diag_get_rank_checked(core_ptr, name, &rank);
    return rank;
}

int catchem_diag_get_rank_checked(void* core_ptr, const char* name, int* rank_out) {
    if (rank_out)
        *rank_out = 0;
    if (!name || !rank_out)
        return fail(catchem::BoundaryStatus::NullArgument, "diagnostic_get_rank", "argument",
                    "field name and rank output are required");
    catchem::AdmissionLease admission;
    const int status = admit_handle(core_ptr, catchem::HandleType::Core, "diagnostic_get_rank", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    auto manager = static_cast<catchem::Core*>(core_ptr)->get_diagnostic_manager();
    if (!manager || !manager->has_field(name))
        return fail(catchem::BoundaryStatus::MissingField, "diagnostic_get_rank", name, "field is not registered");
    *rank_out = static_cast<int>(manager->get_field(name)->dimensions.size());
    return CATCHEM_SUCCESS;
}

void catchem_diag_get_dims(void* core_ptr, const char* name, int* dims_out) {
    (void)catchem_diag_get_dims_checked(core_ptr, name, dims_out, 3);
}

int catchem_diag_get_dims_checked(void* core_ptr, const char* name, int* dims_out, int dims_length) {
    if (dims_out && dims_length > 0)
        std::fill(dims_out, dims_out + dims_length, 0);
    if (!name || !dims_out || dims_length <= 0)
        return fail(catchem::BoundaryStatus::NullArgument, "diagnostic_get_dims", "argument",
                    "field name and a positive-length dimensions output are required");
    catchem::AdmissionLease admission;
    const int status = admit_handle(core_ptr, catchem::HandleType::Core, "diagnostic_get_dims", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    auto manager = static_cast<catchem::Core*>(core_ptr)->get_diagnostic_manager();
    if (!manager || !manager->has_field(name))
        return fail(catchem::BoundaryStatus::MissingField, "diagnostic_get_dims", name, "field is not registered");
    const auto field = manager->get_field(name);
    if (dims_length < static_cast<int>(field->dimensions.size()))
        return fail(catchem::BoundaryStatus::ExtentMismatch, "diagnostic_get_dims", name,
                    "dimensions output is shorter than the field rank");
    std::copy(field->dimensions.begin(), field->dimensions.end(), dims_out);
    return CATCHEM_SUCCESS;
}

int catchem_diag_get_pointer_checked(void* core_ptr, const char* name, int rank, const int* dims, void** ptr_out) {
    if (ptr_out)
        *ptr_out = nullptr;
    if (!name || !dims || !ptr_out)
        return CATCHEM_NULL_ARGUMENT;
    catchem::AdmissionLease admission;
    const int handle_status = admit_handle(core_ptr, catchem::HandleType::Core, "diagnostic_get_pointer", admission);
    if (handle_status != CATCHEM_SUCCESS)
        return handle_status;
    try {
        auto mgr = static_cast<catchem::Core*>(core_ptr)->get_diagnostic_manager();
        if (!mgr || !mgr->has_field(name))
            return CATCHEM_MISSING_FIELD;
        auto field = mgr->get_field(name);
        if (rank != static_cast<int>(field->dimensions.size()))
            return CATCHEM_RANK_MISMATCH;
        for (int i = 0; i < rank; ++i)
            if (dims[i] != field->dimensions[static_cast<std::size_t>(i)])
                return CATCHEM_EXTENT_MISMATCH;
        *ptr_out = const_cast<void*>(mgr->get_host_read_pointer(name));
        return *ptr_out ? CATCHEM_SUCCESS : CATCHEM_INVALID_STATE;
    } catch (...) {
        return CATCHEM_INTERNAL_ERROR;
    }
}

int catchem_diag_mark_host_modified(void* core_ptr, const char* name) {
    if (!core_ptr || !name)
        return CATCHEM_NULL_ARGUMENT;
    try {
        static_cast<catchem::Core*>(core_ptr)->get_diagnostic_manager()->mark_host_modified(name);
        return CATCHEM_SUCCESS;
    } catch (const std::invalid_argument&) {
        return CATCHEM_MISSING_FIELD;
    } catch (...) {
        return CATCHEM_INTERNAL_ERROR;
    }
}

int catchem_diag_mark_device_modified(void* core_ptr, const char* name) {
    if (!core_ptr || !name)
        return CATCHEM_NULL_ARGUMENT;
    try {
        static_cast<catchem::Core*>(core_ptr)->get_diagnostic_manager()->mark_device_modified(name);
        return CATCHEM_SUCCESS;
    } catch (const std::invalid_argument&) {
        return CATCHEM_MISSING_FIELD;
    } catch (...) {
        return CATCHEM_INTERNAL_ERROR;
    }
}

void catchem_diag_sync_to_host(void* core_ptr) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    core->get_diagnostic_manager()->sync_to_host();
}

void catchem_diag_reset(void* core_ptr) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    core->get_diagnostic_manager()->reset_all();
}

int catchem_diag_get_count(void* core_ptr) {
    int count = 0;
    (void)catchem_diag_get_count_checked(core_ptr, &count);
    return count;
}

int catchem_diag_get_count_checked(void* core_ptr, int* count_out) {
    if (!count_out)
        return fail(catchem::BoundaryStatus::NullArgument, "diagnostic_get_count", "count_out",
                    "output pointer is null");
    *count_out = 0;
    catchem::AdmissionLease admission;
    const int status = admit_handle(core_ptr, catchem::HandleType::Core, "diagnostic_get_count", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    *count_out = static_cast<int>(
        static_cast<catchem::Core*>(core_ptr)->get_diagnostic_manager()->get_registered_names().size());
    return CATCHEM_SUCCESS;
}

void catchem_diag_get_name_at(void* core_ptr, int index, char* name_out) {
    (void)catchem_diag_get_name_at_checked(core_ptr, index, name_out, 64);
}

int catchem_diag_get_name_at_checked(void* core_ptr, int index, char* name_out, int name_length) {
    if (name_out && name_length > 0)
        name_out[0] = '\0';
    if (!name_out || name_length <= 0)
        return fail(catchem::BoundaryStatus::NullArgument, "diagnostic_get_name", "name_out",
                    "a positive-length name output is required");
    catchem::AdmissionLease admission;
    const int status = admit_handle(core_ptr, catchem::HandleType::Core, "diagnostic_get_name", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    const auto names = static_cast<catchem::Core*>(core_ptr)->get_diagnostic_manager()->get_registered_names();
    if (index < 0 || index >= static_cast<int>(names.size()))
        return fail(catchem::BoundaryStatus::InvalidIndex, "diagnostic_get_name", "index",
                    "diagnostic index is outside the registered range");
    copy_string_to_buffer(names[static_cast<std::size_t>(index)], name_out, name_length);
    return CATCHEM_SUCCESS;
}

int catchem_diag_get_units_checked(void* core_ptr, const char* name, char* units_out, int units_length) {
    if (units_out && units_length > 0)
        units_out[0] = '\0';
    if (!name || !units_out || units_length <= 0)
        return fail(catchem::BoundaryStatus::NullArgument, "diagnostic_get_units", "name",
                    "field name and a positive-length units output are required");
    catchem::AdmissionLease admission;
    const int status = admit_handle(core_ptr, catchem::HandleType::Core, "diagnostic_get_units", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    auto manager = static_cast<catchem::Core*>(core_ptr)->get_diagnostic_manager();
    if (!manager || !manager->has_field(name))
        return fail(catchem::BoundaryStatus::MissingField, "diagnostic_get_units", name, "field is not registered");
    copy_string_to_buffer(manager->get_field(name)->units, units_out, units_length);
    return CATCHEM_SUCCESS;
}

int catchem_diag_get_description_checked(void* core_ptr, const char* name, char* desc_out, int desc_length) {
    if (desc_out && desc_length > 0)
        desc_out[0] = '\0';
    if (!name || !desc_out || desc_length <= 0)
        return fail(catchem::BoundaryStatus::NullArgument, "diagnostic_get_description", "name",
                    "field name and a positive-length description output are required");
    catchem::AdmissionLease admission;
    const int status = admit_handle(core_ptr, catchem::HandleType::Core, "diagnostic_get_description", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    auto manager = static_cast<catchem::Core*>(core_ptr)->get_diagnostic_manager();
    if (!manager || !manager->has_field(name))
        return fail(catchem::BoundaryStatus::MissingField, "diagnostic_get_description", name,
                    "field is not registered");
    copy_string_to_buffer(manager->get_field(name)->description, desc_out, desc_length);
    return CATCHEM_SUCCESS;
}

int catchem_diag_get_axes_checked(void* core_ptr, const char* name, int* axes_out, int axes_length) {
    if (!name || !axes_out || axes_length <= 0)
        return fail(catchem::BoundaryStatus::NullArgument, "diagnostic_get_axes", "name",
                    "field name and a positive-length axes output are required");
    catchem::AdmissionLease admission;
    const int status = admit_handle(core_ptr, catchem::HandleType::Core, "diagnostic_get_axes", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    auto manager = static_cast<catchem::Core*>(core_ptr)->get_diagnostic_manager();
    if (!manager || !manager->has_field(name))
        return fail(catchem::BoundaryStatus::MissingField, "diagnostic_get_axes", name, "field is not registered");
    const auto& axes = manager->get_axes(name);
    if (axes_length < static_cast<int>(axes.size()))
        return fail(catchem::BoundaryStatus::ExtentMismatch, "diagnostic_get_axes", name,
                    "axes output is shorter than the field rank");
    std::fill(axes_out, axes_out + axes_length, 0);
    for (std::size_t index = 0; index < axes.size(); ++index)
        axes_out[index] = static_cast<int>(axes[index]);
    return CATCHEM_SUCCESS;
}

int catchem_diag_get_unpack_label_at_checked(void* core_ptr, const char* name, int slot, char* label_out,
                                             int label_length) {
    if (label_out && label_length > 0)
        label_out[0] = '\0';
    if (!name || !label_out || label_length <= 0)
        return fail(catchem::BoundaryStatus::NullArgument, "diagnostic_get_unpack_label_at", "name",
                    "field name and a positive-length label output are required");
    catchem::AdmissionLease admission;
    const int status = admit_handle(core_ptr, catchem::HandleType::Core, "diagnostic_get_unpack_label_at", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    auto manager = static_cast<catchem::Core*>(core_ptr)->get_diagnostic_manager();
    if (!manager || !manager->has_field(name))
        return fail(catchem::BoundaryStatus::MissingField, "diagnostic_get_unpack_label_at", name,
                    "field is not registered");
    const auto& labels = manager->get_unpack_labels(name);
    if (labels.empty())
        return fail(catchem::BoundaryStatus::InvalidState, "diagnostic_get_unpack_label_at", name,
                    "field has no packed dimension to unpack");
    if (slot < 0 || slot >= static_cast<int>(labels.size()))
        return fail(catchem::BoundaryStatus::InvalidIndex, "diagnostic_get_unpack_label_at", name,
                    "label slot is outside the packed dimension");
    copy_string_to_buffer(labels[static_cast<std::size_t>(slot)], label_out, label_length);
    return CATCHEM_SUCCESS;
}
}
