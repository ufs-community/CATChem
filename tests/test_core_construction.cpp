#include "catchem_core.hpp"
#include "catchem_test_config.hpp"
#include <cassert>
#include <string>

namespace {
    void assert_common_wiring(catchem::Core& core, int nx, int ny, int nz, int species) {
        const auto config = core.get_config_manager();
        const auto grid = core.get_grid_manager();
        const auto state = core.get_state_manager();
        assert(config && grid && state && core.get_diagnostic_manager());
        assert(config->data.runtime.nx == nx);
        assert(config->data.runtime.ny == ny);
        assert(config->data.runtime.nz == nz);
        assert(grid->geometry.nx == nx);
        assert(grid->geometry.ny == ny);
        assert(grid->geometry.nz == nz);
        assert(state->column_count() == nx * ny);
        assert(state->level_count() == nz);
        assert(state->species_count() == species);
        assert(state->config_manager() == config);
        assert(state->diagnostic_manager() == core.get_diagnostic_manager());
    }
} // namespace

int main() {
    catchem::Core direct(2, 3, 3);
    assert_common_wiring(direct, 2, 1, 3, 3);

    const std::string config = std::string(catchem::test::SOURCE_DIR) + "/tests/fixtures/platform_integrity_valid.yml";
    catchem::Core configured(config);
    assert_common_wiring(configured, 2, 1, 3, 3);

    catchem::Core host_grid(config, 5, 4);
    assert_common_wiring(host_grid, 5, 1, 4, 3);

    catchem::Core options(catchem::CoreCreateOptions::configured_with_host_grid(config, 5, 4));
    assert_common_wiring(options, 5, 1, 4, 3);

    const std::string disabled =
        std::string(catchem::test::SOURCE_DIR) + "/tests/fixtures/configured_process_disabled.yml";
    catchem::Core disabled_process(disabled);
    assert_common_wiring(disabled_process, 2, 1, 3, 3);
    assert(disabled_process.get_num_processes() == 0);
    return 0;
}
