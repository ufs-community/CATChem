#include "catchem_core.hpp"
#include "catchem_core_architecture_test_helpers.hpp"
#include <cassert>
#include <stdexcept>
#include <vector>

class FinalizeProcess final : public catchem::test::RecordingProcess {
public:
    FinalizeProcess(std::string name, std::vector<std::string>& events, bool fail)
        : RecordingProcess(std::move(name), events), fail_(fail) {}
    void finalize() override {
        RecordingProcess::finalize();
        if (fail_)
            throw std::runtime_error("injected finalize failure");
    }

private:
    bool fail_;
};

int main() {
    std::vector<std::string> events;
    catchem::Core core(1, 1, 1);
    auto state = core.get_state_manager();
    auto first = std::make_shared<FinalizeProcess>("first", events, false);
    auto second = std::make_shared<FinalizeProcess>("second", events, true);
    auto third = std::make_shared<FinalizeProcess>("third", events, false);
    for (const auto& process : {first, second, third}) {
        process->init(state);
        core.add_process(process);
    }
    bool failed = false;
    try {
        core.shutdown();
    } catch (const std::runtime_error&) {
        failed = true;
    }
    assert(failed);
    assert((events == std::vector<std::string>{"init:first", "init:second", "init:third", "finalize:third",
                                               "finalize:second", "finalize:first"}));
    core.shutdown();
    assert(events.size() == 6);

    std::vector<std::string> partial_events;
    catchem::Core partial(1, 1, 1);
    auto initialized = std::make_shared<catchem::test::RecordingProcess>("initialized", partial_events);
    initialized->init(partial.get_state_manager());
    partial.add_process(initialized);
    auto rejected = std::make_shared<catchem::test::RecordingProcess>("rejected", partial_events, true);
    try {
        rejected->init(partial.get_state_manager());
    } catch (const std::runtime_error&) {
    }
    partial.shutdown();
    assert((partial_events == std::vector<std::string>{"init:initialized", "init:rejected", "finalize:initialized"}));
    return 0;
}
