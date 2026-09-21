#include "catchem_state_manager.hpp"
#include <cstdlib>

namespace catchem {

    StateManager::StateManager(int nc, int nl, int ns) : n_cols(nc), n_levels(nl), n_species(ns) {
        // Generate standard 8-character trace id (alphanumeric random)
        static const char alphanum[] = "0123456789abcdefghijklmnopqrstuvwxyz";
        trace_id = "";
        for (int i = 0; i < 8; ++i) {
            trace_id += alphanum[rand() % (sizeof(alphanum) - 1)];
        }
    }

} // namespace catchem
