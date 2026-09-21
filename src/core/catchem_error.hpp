#pragma once

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace catchem {

    enum class BoundaryStatus : int {
        Success = 0,
        NullArgument = 1,
        MissingField = 2,
        RankMismatch = 3,
        ExtentMismatch = 4,
        InvalidIndex = 5,
        StaleGeneration = 6,
        DuplicateMapping = 7,
        InvalidState = 8,
        InternalError = 9,
        InvalidHandle = 10,
        WrongHandleType = 11,
        RuntimeUnavailable = 12,
        ContractViolation = 13,
        InvalidConfiguration = 14,
        ProcessFailure = 15,
        ShutdownFailure = 16,
        PhysicalValidationFailure = 17
    };

    struct BoundaryError {
        BoundaryStatus status = BoundaryStatus::Success;
        std::string operation;
        std::string object;
        std::string cause;
    };

    inline thread_local BoundaryError boundary_error;

    inline void clear_boundary_error() noexcept {
        boundary_error = {};
    }

    inline BoundaryStatus set_boundary_error(BoundaryStatus status, std::string operation, std::string object,
                                             std::string cause) noexcept {
        boundary_error = {status, std::move(operation), std::move(object), std::move(cause)};
        return status;
    }

    inline const BoundaryError& last_boundary_error() noexcept {
        return boundary_error;
    }

    template <typename Function>
    BoundaryStatus guard_boundary_call(const char* operation, const char* object, Function&& function) noexcept {
        clear_boundary_error();
        try {
            function();
            return BoundaryStatus::Success;
        } catch (const std::invalid_argument& error) {
            return set_boundary_error(BoundaryStatus::InvalidState, operation ? operation : "", object ? object : "",
                                      error.what());
        } catch (const std::exception& error) {
            return set_boundary_error(BoundaryStatus::InternalError, operation ? operation : "", object ? object : "",
                                      error.what());
        } catch (...) {
            return set_boundary_error(BoundaryStatus::InternalError, operation ? operation : "", object ? object : "",
                                      "unknown non-standard exception");
        }
    }

    enum ErrorCode { SUCCESS = 0, FAILURE = -1, INVALID_INPUT = 1001, INVALID_STATE = 1003, MEMORY_ALLOCATION = 1007 };

    class ErrorManager {
    private:
        std::vector<std::string> context_stack;

    public:
        void push_context(const std::string& ctx) { context_stack.push_back(ctx); }
        void pop_context() {
            if (!context_stack.empty()) {
                context_stack.pop_back();
            }
        }
        void report_error(ErrorCode code, const std::string& msg) {
            std::cerr << "[CATChem C++ Error " << code << "] " << msg << " | Context: ";
            for (const auto& ctx : context_stack) {
                std::cerr << ctx << " -> ";
            }
            std::cerr << "End\n";
        }
    };

    inline void require_field_pointer(const char* process_name, const char* field_name, const void* ptr) {
        if (ptr == nullptr) {
            throw std::runtime_error(std::string("FATAL ERROR: ") + process_name + " process missing required field " +
                                     field_name);
        }
    }

} // namespace catchem
