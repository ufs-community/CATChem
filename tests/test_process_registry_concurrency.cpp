#include "catchem_core_architecture_test_helpers.hpp"
#include "catchem_process_registry.hpp"
#include <atomic>
#include <cassert>
#include <thread>

int main() {
    auto& registry = catchem::ProcessRegistry::get_instance();
    std::atomic<int> created{0};
    registry.register_process("concurrent-test", [&] {
        ++created;
        static thread_local std::vector<std::string> events;
        return std::make_shared<catchem::test::RecordingProcess>("concurrent-test", events);
    });
    std::vector<std::thread> workers;
    for (int thread = 0; thread < 8; ++thread) {
        workers.emplace_back([&] {
            for (int operation = 0; operation < 1250; ++operation) {
                assert(registry.has_process("concurrent-test"));
                assert(registry.create("concurrent-test"));
            }
        });
    }
    for (auto& worker : workers)
        worker.join();
    assert(created == 10000);
    return 0;
}
