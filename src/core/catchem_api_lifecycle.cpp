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
void* catchem_core_create(int nc, int nl, int ns) {
    void* core = nullptr;
    (void)catchem_core_create_checked(nc, nl, ns, &core);
    return core;
}

int catchem_core_create_checked(int nc, int nl, int ns, void** core_out) {
    catchem::clear_boundary_error();
    if (!core_out)
        return fail(catchem::BoundaryStatus::NullArgument, "core_create", "core_out", "output pointer is null");
    *core_out = nullptr;
    if (nc <= 0 || nl <= 0 || ns <= 0)
        return fail(catchem::BoundaryStatus::ExtentMismatch, "core_create", "dimensions",
                    "column, level, and species extents must be positive");
    try {
        auto* core = new catchem::Core(nc, nl, ns);
        if (!catchem::HandleRegistry::instance().add(
                core, {catchem::HandleType::Core, 1, 1, nullptr, catchem::HandleOwnership::Owned})) {
            delete core;
            return fail(catchem::BoundaryStatus::InternalError, "core_create", "registry",
                        "could not register core handle");
        }
        *core_out = core;
        return CATCHEM_SUCCESS;
    } catch (const std::exception& error) {
        return fail(catchem::BoundaryStatus::InternalError, "core_create", "core", error.what());
    } catch (...) {
        return fail(catchem::BoundaryStatus::InternalError, "core_create", "core", "unknown exception");
    }
}

void* catchem_core_create_from_config(const char* config_file) {
    void* core = nullptr;
    (void)catchem_core_create_from_config_checked(config_file, &core);
    return core;
}

int catchem_core_create_from_config_checked(const char* config_file, void** core_out) {
    catchem::clear_boundary_error();
    if (core_out)
        *core_out = nullptr;
    if (!config_file || !core_out)
        return fail(catchem::BoundaryStatus::NullArgument, "core_create_from_config", "argument",
                    "configuration path and output pointer are required");
    try {
        register_builtin_processes();
        auto* core = new catchem::Core(config_file);
        catchem::HandleRegistry::instance().add(
            core, {catchem::HandleType::Core, 1, 1, nullptr, catchem::HandleOwnership::Owned});
        *core_out = core;
        return CATCHEM_SUCCESS;
    } catch (const std::exception& error) {
        return fail(catchem::BoundaryStatus::InvalidState, "core_create_from_config", config_file, error.what());
    } catch (...) {
        return fail(catchem::BoundaryStatus::InternalError, "core_create_from_config", config_file,
                    "unknown exception");
    }
}

void* catchem_core_create_from_config_with_grid(const char* config_file, int ncols, int nlevels) {
    void* core = nullptr;
    (void)catchem_core_create_from_config_with_grid_checked(config_file, ncols, nlevels, &core);
    return core;
}

int catchem_core_create_from_config_with_grid_checked(const char* config_file, int ncols, int nlevels,
                                                      void** core_out) {
    catchem::clear_boundary_error();
    if (core_out)
        *core_out = nullptr;
    if (!config_file || !core_out)
        return fail(catchem::BoundaryStatus::NullArgument, "core_create_from_config_with_grid", "argument",
                    "configuration path and output pointer are required");
    if (ncols <= 0 || nlevels <= 0)
        return fail(catchem::BoundaryStatus::ExtentMismatch, "core_create_from_config_with_grid", "grid",
                    "column and level extents must be positive");
    try {
        register_builtin_processes();
        auto* core = new catchem::Core(config_file, ncols, nlevels);
        catchem::HandleRegistry::instance().add(
            core, {catchem::HandleType::Core, 1, 1, nullptr, catchem::HandleOwnership::Owned});
        *core_out = core;
        return CATCHEM_SUCCESS;
    } catch (const std::exception& error) {
        return fail(catchem::BoundaryStatus::InvalidState, "core_create_from_config_with_grid", config_file,
                    error.what());
    } catch (...) {
        return fail(catchem::BoundaryStatus::InternalError, "core_create_from_config_with_grid", config_file,
                    "unknown exception");
    }
}

void catchem_core_destroy(void* core_ptr) {
    (void)catchem_core_destroy_checked(core_ptr);
}

int catchem_core_destroy_checked(void* core_ptr) {
    catchem::clear_boundary_error();
    const auto status = catchem::HandleRegistry::instance().close_and_wait(core_ptr, catchem::HandleType::Core);
    if (status != catchem::BoundaryStatus::Success)
        return fail(status, "core_destroy", "core", "null, stale, destroyed, or wrong-type handle");
    catchem::HandleRegistry::instance().invalidate_children(core_ptr);
    auto* core = static_cast<catchem::Core*>(core_ptr);
    try {
        core->shutdown();
        delete core;
        return CATCHEM_SUCCESS;
    } catch (const std::exception& error) {
        const std::string cause = error.what();
        delete core;
        return fail(catchem::BoundaryStatus::ShutdownFailure, "core_destroy", "process cleanup", cause.c_str());
    } catch (...) {
        delete core;
        return fail(catchem::BoundaryStatus::ShutdownFailure, "core_destroy", "process cleanup",
                    "unknown cleanup failure");
    }
}

void* catchem_core_get_state_manager(void* core_ptr) {
    void* state = nullptr;
    (void)catchem_core_get_state_manager_checked(core_ptr, &state);
    return state;
}

int catchem_core_get_state_manager_checked(void* core_ptr, void** state_out) {
    catchem::clear_boundary_error();
    if (!state_out)
        return fail(catchem::BoundaryStatus::NullArgument, "core_get_state", "state_out", "output pointer is null");
    *state_out = nullptr;
    catchem::AdmissionLease admission;
    const int status = admit_handle(core_ptr, catchem::HandleType::Core, "core_get_state", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    auto* core = static_cast<catchem::Core*>(core_ptr);
    auto* state = core->get_state_manager().get();
    const auto existing = catchem::HandleRegistry::instance().validate(state, catchem::HandleType::State, core_ptr);
    if (existing != catchem::BoundaryStatus::Success) {
        catchem::HandleRegistry::instance().remove(state);
        catchem::HandleRegistry::instance().add(
            state, {catchem::HandleType::State, 1, 1, core_ptr, catchem::HandleOwnership::Borrowed});
    }
    *state_out = state;
    return CATCHEM_SUCCESS;
}

int catchem_get_last_error(char* buffer, int max_len) {
    if (!buffer || max_len <= 0)
        return CATCHEM_NULL_ARGUMENT;
    const auto& error = catchem::last_boundary_error();
    std::string message = error.operation;
    if (!error.object.empty())
        message += " [" + error.object + "]";
    if (!error.cause.empty())
        message += ": " + error.cause;
    copy_string_to_buffer(message, buffer, max_len);
    return status_code(error.status);
}

int catchem_core_run_timestep(void* core_ptr, double dt) {
    catchem::AdmissionLease admission;
    const int handle_status = admit_handle(core_ptr, catchem::HandleType::Core, "core_run_timestep", admission);
    if (handle_status != CATCHEM_SUCCESS)
        return handle_status;
    try {
        auto* core = static_cast<catchem::Core*>(core_ptr);
        core->run_timestep(dt);
        return 0;
    } catch (const std::exception& error) {
        return fail(catchem::BoundaryStatus::InvalidState, "core_run_timestep", "timestep", error.what());
    } catch (...) {
        return fail(catchem::BoundaryStatus::InternalError, "core_run_timestep", "timestep", "unknown exception");
    }
}

int catchem_core_get_timestep_outcome(void* core_ptr, int* status, long long* timestep, double* duration,
                                      long long* import_generation, int* process_index, int* state_classification,
                                      char* process_name, int process_name_len, char* cause, int cause_len) {
    if (!status || !timestep || !duration || !import_generation || !process_index || !state_classification ||
        !process_name || !cause)
        return CATCHEM_NULL_ARGUMENT;
    catchem::AdmissionLease admission;
    const int handle_status = admit_handle(core_ptr, catchem::HandleType::Core, "core_get_timestep_outcome", admission);
    if (handle_status != CATCHEM_SUCCESS)
        return handle_status;
    const auto& outcome = static_cast<catchem::Core*>(core_ptr)->get_timestep_outcome();
    *status = static_cast<int>(outcome.status);
    *timestep = static_cast<long long>(outcome.timestep);
    *duration = outcome.duration;
    *import_generation = static_cast<long long>(outcome.import_generation);
    *process_index = static_cast<int>(outcome.process_index);
    *state_classification = static_cast<int>(outcome.state);
    copy_string_to_buffer(outcome.process_name, process_name, process_name_len);
    copy_string_to_buffer(outcome.cause, cause, cause_len);
    return CATCHEM_SUCCESS;
}

void catchem_core_add_process_by_name(void* core_ptr, const char* name) {
    register_builtin_processes();
    auto* core = static_cast<catchem::Core*>(core_ptr);
    auto process = catchem::ProcessRegistry::get_instance().create(name);
    process->init(core->get_state_manager());
    core->add_process(process);
}

int catchem_core_get_num_processes(void* core_ptr) {
    int count = 0;
    (void)catchem_core_get_num_processes_checked(core_ptr, &count);
    return count;
}

int catchem_core_get_num_processes_checked(void* core_ptr, int* count_out) {
    if (!count_out)
        return fail(catchem::BoundaryStatus::NullArgument, "core_get_num_processes", "count_out",
                    "output pointer is null");
    *count_out = 0;
    catchem::AdmissionLease admission;
    const int status = admit_handle(core_ptr, catchem::HandleType::Core, "core_get_num_processes", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    *count_out = static_cast<int>(static_cast<catchem::Core*>(core_ptr)->get_num_processes());
    return CATCHEM_SUCCESS;
}

int catchem_core_get_required_host_field_count_checked(void* core_ptr, int* count_out) {
    if (!count_out)
        return fail(catchem::BoundaryStatus::NullArgument, "core_get_required_host_field_count", "count_out",
                    "output pointer is null");
    *count_out = 0;
    catchem::AdmissionLease admission;
    const int status =
        admit_handle(core_ptr, catchem::HandleType::Core, "core_get_required_host_field_count", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    *count_out = static_cast<int>(static_cast<catchem::Core*>(core_ptr)->get_required_host_fields().size());
    return CATCHEM_SUCCESS;
}

int catchem_core_get_required_host_field_name_checked(void* core_ptr, int index, char* name_out, int name_out_len) {
    if (!name_out || name_out_len <= 0)
        return fail(catchem::BoundaryStatus::NullArgument, "core_get_required_host_field_name", "name_out",
                    "a non-empty output buffer is required");
    name_out[0] = '\0';
    catchem::AdmissionLease admission;
    const int status =
        admit_handle(core_ptr, catchem::HandleType::Core, "core_get_required_host_field_name", admission);
    if (status != CATCHEM_SUCCESS)
        return status;
    const auto fields = static_cast<catchem::Core*>(core_ptr)->get_required_host_fields();
    if (index < 0 || static_cast<std::size_t>(index) >= fields.size())
        return fail(catchem::BoundaryStatus::InvalidIndex, "core_get_required_host_field_name", "index",
                    "required host field index is out of range");
    copy_string_to_buffer(fields[static_cast<std::size_t>(index)], name_out, name_out_len);
    return CATCHEM_SUCCESS;
}
}
