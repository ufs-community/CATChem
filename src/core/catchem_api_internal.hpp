#pragma once

#include "catchem_boundary_guard.hpp"
#include "catchem_error.hpp"
#include "catchem_handle_registry.hpp"
#include <cstring>
#include <mutex>
#include <string>

#if defined(__GNUC__) || defined(__clang__)
#define CATCHEM_WEAK_SYMBOL __attribute__((weak))
#else
#define CATCHEM_WEAK_SYMBOL
#endif

extern "C" {
void catchem_register_carbchem_cpp() CATCHEM_WEAK_SYMBOL;
void catchem_register_drydep_cpp() CATCHEM_WEAK_SYMBOL;
void catchem_register_dust_cpp() CATCHEM_WEAK_SYMBOL;
void catchem_register_seasalt_cpp() CATCHEM_WEAK_SYMBOL;
void catchem_register_settling_cpp() CATCHEM_WEAK_SYMBOL;
void catchem_register_so4chem_cpp() CATCHEM_WEAK_SYMBOL;
void catchem_register_wetdep_cpp() CATCHEM_WEAK_SYMBOL;
}

namespace catchem::api_internal {

    inline int status_code(BoundaryStatus status) {
        return static_cast<int>(status);
    }

    inline int admit_handle(void* handle, HandleType type, const char* operation, AdmissionLease& admission) {
        return status_code(BoundaryGuard::admit(handle, type, operation, admission));
    }

    inline int fail(BoundaryStatus status, const char* operation, const char* object, const char* cause) {
        return status_code(set_boundary_error(status, operation, object ? object : "", cause));
    }

    inline void copy_string_to_buffer(const std::string& source, char* buffer, int max_length) {
        if (!buffer || max_length <= 0)
            return;
        std::strncpy(buffer, source.c_str(), static_cast<std::size_t>(max_length - 1));
        buffer[max_length - 1] = '\0';
    }

    inline void register_builtin_processes() {
        static std::once_flag registered;
        std::call_once(registered, [] {
            if (catchem_register_carbchem_cpp)
                catchem_register_carbchem_cpp();
            if (catchem_register_drydep_cpp)
                catchem_register_drydep_cpp();
            if (catchem_register_dust_cpp)
                catchem_register_dust_cpp();
            if (catchem_register_seasalt_cpp)
                catchem_register_seasalt_cpp();
            if (catchem_register_settling_cpp)
                catchem_register_settling_cpp();
            if (catchem_register_so4chem_cpp)
                catchem_register_so4chem_cpp();
            if (catchem_register_wetdep_cpp)
                catchem_register_wetdep_cpp();
        });
    }

} // namespace catchem::api_internal
