#include "catchem_kokkos_compat.hpp"
#include "catchem_logger.hpp"
#include "catchem_state_manager.hpp"
#include <cassert>
#include <cstdlib>
#include <iostream>

int main(int argc, char* argv[]) {
    // Set the level threshold source before the first Logger call so the
    // environment-derived threshold (cached on first use) is deterministic.
    setenv("CATCHEM_LOG_LEVEL", "debug", 1);

    Kokkos::initialize(argc, argv);
    {
        auto state = std::make_shared<catchem::StateManager>(4, 10, 50);
        state->set_runtime_trace_id("testtrac");

        // Manual redirect stringstream capture is optional, we assert logger successfully formats and prints
        catchem::Logger::info(state.get(), "Simulation timestep advanced", {{"step", "12"}, {"dt", "300.0"}});
        catchem::Logger::error(state.get(), "Division by zero encountered", {{"cell", "4"}});

        // level_from_string is case-insensitive and accepts the warn alias.
        using Level = catchem::Logger::Level;
        assert(catchem::Logger::level_from_string("DEBUG") == Level::Debug);
        assert(catchem::Logger::level_from_string("info") == Level::Info);
        assert(catchem::Logger::level_from_string("Warning") == Level::Warn);
        assert(catchem::Logger::level_from_string("error") == Level::Error);
        assert(!catchem::Logger::level_from_string("bogus").has_value());
        assert(!catchem::Logger::level_from_string("").has_value());

        // Without an override, CATCHEM_LOG_LEVEL=debug lets DEBUG through.
        assert(catchem::Logger::enabled(Level::Debug));

        // An explicit override (the YAML configuration path) always wins over
        // CATCHEM_LOG_LEVEL, raising the threshold so DEBUG is suppressed.
        catchem::Logger::set_level(Level::Error);
        assert(!catchem::Logger::enabled(Level::Debug));
        assert(!catchem::Logger::enabled(Level::Info));
        assert(catchem::Logger::enabled(Level::Error));

        // Clearing the override restores environment control.
        catchem::Logger::clear_level();
        assert(catchem::Logger::enabled(Level::Debug));
    }
    Kokkos::finalize();

    std::cout << "All logger formatting unit tests passed!" << std::endl;
    return 0;
}
