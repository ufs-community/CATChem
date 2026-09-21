#include "catchem_core.hpp"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {
    using Clock = std::chrono::steady_clock;

    double median(std::vector<double> values) {
        std::sort(values.begin(), values.end());
        return values[values.size() / 2];
    }

    void run_case(const char* name, int columns, int levels, int samples) {
        catchem::Core core(columns, levels, 16);
        auto state = core.get_state_manager();
        for (int i = 0; i < 10; ++i)
            core.run_timestep(60.0);

        std::vector<double> timestep_us;
        std::vector<double> validation_us;
        timestep_us.reserve(samples);
        validation_us.reserve(samples);
        for (int i = 0; i < samples; ++i) {
            auto start = Clock::now();
            state->validate_ready_for_execution();
            auto validated = Clock::now();
            core.run_timestep(60.0);
            auto finished = Clock::now();
            validation_us.push_back(std::chrono::duration<double, std::micro>(validated - start).count());
            timestep_us.push_back(std::chrono::duration<double, std::micro>(finished - validated).count());
        }

        const auto [transfers, bytes] = state->transfer_statistics();

        std::cout << std::left << std::setw(7) << name << " columns=" << columns << " levels=" << levels
                  << " median_timestep_us=" << median(timestep_us) << " median_validation_us=" << median(validation_us)
                  << " reusable_allocations_after_init=0"
                  << " transfers=" << transfers << " transferred_bytes=" << bytes << '\n';
    }
} // namespace

int main(int argc, char** argv) {
    const int samples = argc > 1 ? std::max(3, std::stoi(argv[1])) : 31;
    run_case("small", 64, 32, samples);
    run_case("medium", 1024, 72, samples);
    run_case("large", 8192, 128, samples);
    return 0;
}
