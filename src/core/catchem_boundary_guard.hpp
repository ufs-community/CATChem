#pragma once
#include "catchem_error.hpp"
#include "catchem_handle_registry.hpp"
#include <string>
#include <utility>

namespace catchem {

    class BoundaryGuard {
    public:
        static BoundaryStatus admit(const void* handle, HandleType type, const char* operation,
                                    AdmissionLease& lease) noexcept {
            auto acquired = HandleRegistry::instance().acquire(handle, type);
            if (acquired.first != BoundaryStatus::Success)
                return set_boundary_error(acquired.first, operation ? operation : "", "handle",
                                          "null, stale, destroyed, closing, or wrong-type handle");
            lease = std::move(acquired.second);
            return BoundaryStatus::Success;
        }

        template <typename Output>
        static bool initialize_output(Output* output, Output safe_value, const char* operation,
                                      const char* name) noexcept {
            if (!output) {
                set_boundary_error(BoundaryStatus::NullArgument, operation ? operation : "", name ? name : "output",
                                   "output pointer is null");
                return false;
            }
            *output = std::move(safe_value);
            return true;
        }

        template <typename Function>
        static BoundaryStatus invoke(const char* operation, const char* object, Function&& function) noexcept {
            clear_boundary_error();
            try {
                std::forward<Function>(function)();
                return BoundaryStatus::Success;
            } catch (const std::invalid_argument& error) {
                return set_boundary_error(BoundaryStatus::ContractViolation, operation ? operation : "",
                                          object ? object : "", error.what());
            } catch (const std::exception& error) {
                return set_boundary_error(BoundaryStatus::InternalError, operation ? operation : "",
                                          object ? object : "", error.what());
            } catch (...) {
                return set_boundary_error(BoundaryStatus::InternalError, operation ? operation : "",
                                          object ? object : "", "unknown non-standard exception");
            }
        }
    };

} // namespace catchem
