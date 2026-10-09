#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_core_architecture_test_helpers.hpp"
#include <cassert>
#include <string>
#include <thread>

class CleanupFailure final : public catchem::test::RecordingProcess {
public:
    using RecordingProcess::RecordingProcess;
    void finalize() override {
        RecordingProcess::finalize();
        throw std::runtime_error("cleanup sentinel");
    }
};

int main() {
    std::string first_error, second_error;
    std::thread first([&] {
        assert(catchem_state_begin_import_generation(nullptr) == CATCHEM_NULL_ARGUMENT);
        char detail[128] = {};
        assert(catchem_get_last_error(detail, sizeof(detail)) == CATCHEM_NULL_ARGUMENT);
        first_error = detail;
    });
    std::thread second([&] {
        void* output = reinterpret_cast<void*>(1);
        assert(catchem_core_get_state_manager_checked(reinterpret_cast<void*>(0x12345), &output) ==
               CATCHEM_INVALID_HANDLE);
        assert(output == nullptr);
        char detail[128] = {};
        assert(catchem_get_last_error(detail, sizeof(detail)) == CATCHEM_INVALID_HANDLE);
        second_error = detail;
    });
    first.join();
    second.join();
    assert(first_error.find("state_begin_import_generation") != std::string::npos);
    assert(second_error.find("core_get_state") != std::string::npos);

    assert(catchem_state_begin_import_generation(nullptr) == CATCHEM_NULL_ARGUMENT);
    char truncated[5] = {'x', 'x', 'x', 'x', 'x'};
    assert(catchem_get_last_error(truncated, sizeof(truncated)) == CATCHEM_NULL_ARGUMENT);
    assert(truncated[4] == '\0');

    void* handle = catchem_core_create(1, 1, 1);
    auto* core = static_cast<catchem::Core*>(handle);
    std::vector<std::string> events;
    auto process = std::make_shared<CleanupFailure>("cleanup", events);
    process->init(core->get_state_manager());
    core->add_process(process);
    assert(catchem_core_destroy_checked(handle) == CATCHEM_SHUTDOWN_FAILURE);
    char detail[128] = {};
    assert(catchem_get_last_error(detail, sizeof(detail)) == CATCHEM_SHUTDOWN_FAILURE);
    assert(std::string(detail).find("cleanup sentinel") != std::string::npos);
    return 0;
}
