#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>

class CountingProcess final : public catchem::ProcessInterface {
public:
    int runs = 0;
    std::string get_name() const override { return "counting"; }
    void init(std::shared_ptr<catchem::StateManager>) override {}
    void run(std::shared_ptr<catchem::StateManager>) override { ++runs; }
    void finalize() override {}
};

class FailOnceProcess final : public catchem::ProcessInterface {
public:
    bool fail = true;
    std::string get_name() const override { return "injected_failure"; }
    void init(std::shared_ptr<catchem::StateManager>) override {}
    void run(std::shared_ptr<catchem::StateManager>) override {
        if (fail) {
            fail = false;
            throw std::runtime_error("injected process failure");
        }
    }
    void finalize() override {}
};

int main() {
    void* handle = catchem_core_create(2, 3, 1);
    assert(handle);
    auto* core = static_cast<catchem::Core*>(handle);
    auto first = std::make_shared<CountingProcess>();
    auto failing = std::make_shared<FailOnceProcess>();
    auto never = std::make_shared<CountingProcess>();
    core->add_process(first);
    core->add_process(failing);
    core->add_process(never);

    assert(catchem_core_run_timestep(handle, 60.0) == CATCHEM_INVALID_STATE);
    assert(first->runs == 1);
    assert(never->runs == 0);
    int status = -1, process_index = -1, classification = -1;
    long long timestep = -1, generation = -1;
    double duration = 0.0;
    char process[64] = {}, cause[128] = {};
    assert(catchem_core_get_timestep_outcome(handle, &status, &timestep, &duration, &generation, &process_index,
                                             &classification, process, sizeof(process), cause,
                                             sizeof(cause)) == CATCHEM_SUCCESS);
    assert(status == static_cast<int>(catchem::TimestepStatus::PartialUpdate));
    assert(timestep == 1 && duration == 60.0 && process_index == 1);
    assert(classification == static_cast<int>(catchem::StateClassification::RequiresReimport));
    assert(std::string(process) == "injected_failure");
    assert(std::string(cause).find("injected") != std::string::npos);

    assert(catchem_core_run_timestep(handle, 60.0) == CATCHEM_INVALID_STATE);
    assert(first->runs == 1);
    void* state = catchem_core_get_state_manager(handle);
    assert(catchem_state_begin_import_generation(state) == CATCHEM_SUCCESS);
    assert(catchem_core_run_timestep(handle, 60.0) == CATCHEM_SUCCESS);
    assert(first->runs == 2 && never->runs == 1);
    assert(core->get_timestep_outcome().status == catchem::TimestepStatus::Succeeded);
    assert(catchem_core_destroy_checked(handle) == CATCHEM_SUCCESS);
    return 0;
}
