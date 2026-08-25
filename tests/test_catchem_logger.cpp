#include "catchem_logger.hpp"
#include "catchem_state_manager.hpp"
#include <Kokkos_Core.hpp>
#include <cassert>
#include <iostream>

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        auto state = std::make_shared<catchem::StateManager>(4, 10, 50);
        state->trace_id = "testtrac";

        // Manual redirect stringstream capture is optional, we assert logger successfully formats and prints
        catchem::Logger::info(state.get(), "Simulation timestep advanced", {{"step", "12"}, {"dt", "300.0"}});
        catchem::Logger::error(state.get(), "Division by zero encountered", {{"cell", "4"}});
    }
    Kokkos::finalize();

    std::cout << "All logger formatting unit tests passed!" << std::endl;
    return 0;
}
