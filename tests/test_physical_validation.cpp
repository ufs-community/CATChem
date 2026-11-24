#include "catchem_kokkos_compat.hpp"
#include "catchem_state_manager.hpp"
#include <cassert>
#include <limits>
#include <vector>

int main(int argc, char* argv[]) {
    // The core's derive_* kernels run through Kokkos when CATCHEM_ENABLE_KOKKOS
    // is compiled in (MUSICA builds force it on); the compat header no-ops
    // these calls for host-only builds.
    Kokkos::initialize(argc, argv);
    {
        catchem::StateManager state(1, 2, 1);
        std::vector<double> edge{80000.0, 90000.0, 70000.0};
        std::vector<double> temperature{280.0, -5.0};
        state.bind_met_field_3d("PEDGE", edge.data());
        state.bind_met_field_3d("T", temperature.data());
        state.set_validation_policy(catchem::PhysicalValidationPolicy::Reject);
        bool rejected = false;
        try {
            state.derive_bxheight();
        } catch (const std::domain_error&) {
            rejected = true;
        }
        assert(rejected);
        assert(!state.meteorology().BXHEIGHT); // validation precedes output allocation/mutation
        assert(state.validation_report().issue_count() == 2);

        state.set_validation_policy(catchem::PhysicalValidationPolicy::WarnAndClamp);
        state.derive_bxheight();
        assert(state.meteorology().BXHEIGHT);
        assert(!state.validation_report().empty());
        assert(state.validation_report().issues().front().locations.size() <= 16);

        catchem::StateManager density(1, 1, 1);
        std::vector<double> pressure{90000.0};
        std::vector<double> bad_temperature{std::numeric_limits<double>::quiet_NaN()};
        std::vector<double> humidity{1.5};
        density.bind_met_field_3d("PMID", pressure.data());
        density.bind_met_field_3d("T", bad_temperature.data());
        density.bind_met_field_3d("QV", humidity.data());
        density.set_validation_policy(catchem::PhysicalValidationPolicy::Reject);
        rejected = false;
        try {
            density.derive_airden_dry();
        } catch (const std::domain_error&) {
            rejected = true;
        }
        assert(rejected);
        assert(!density.meteorology().AIRDEN_DRY);
        density.set_validation_policy(catchem::PhysicalValidationPolicy::CountAndContinue);
        density.derive_airden_dry();
        assert(density.meteorology().AIRDEN_DRY);
        assert(density.validation_report().issue_count() == 2);
    }
    Kokkos::finalize();
    return 0;
}
