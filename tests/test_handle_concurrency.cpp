#include "catchem_api.hpp"
#include <atomic>
#include <cassert>
#include <thread>
#include <vector>

int main() {
    void* core = catchem_core_create(2, 3, 1);
    assert(core);
    std::atomic<bool> start{false};
    std::atomic<int> ready{0};
    std::atomic<int> successes{0};
    std::atomic<int> rejections{0};
    std::vector<std::thread> workers;
    for (int thread = 0; thread < 8; ++thread) {
        workers.emplace_back([&] {
            ++ready;
            while (!start.load(std::memory_order_acquire))
                std::this_thread::yield();
            for (int operation = 0; operation < 1250; ++operation) {
                void* state = nullptr;
                const int status = catchem_core_get_state_manager_checked(core, &state);
                if (status == CATCHEM_SUCCESS) {
                    assert(state);
                    ++successes;
                } else {
                    assert(status == CATCHEM_INVALID_HANDLE);
                    assert(state == nullptr);
                    ++rejections;
                }
            }
        });
    }
    while (ready.load() != 8)
        std::this_thread::yield();
    start.store(true, std::memory_order_release);
    assert(catchem_core_destroy_checked(core) == CATCHEM_SUCCESS);
    for (auto& worker : workers)
        worker.join();
    assert(successes + rejections == 10000);
    void* state = reinterpret_cast<void*>(1);
    assert(catchem_core_get_state_manager_checked(core, &state) == CATCHEM_INVALID_HANDLE);
    assert(state == nullptr);
    return 0;
}
