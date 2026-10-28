#include "catchem_state_manager.hpp"
#include <cassert>
#include <vector>

int main() {
    catchem::StateManager state(2, 3, 1);
    std::vector<double> temperature(6, 280.0);
    state.bind_met_field_3d("T", temperature.data());
    auto field = state.meteorology().T;
    catchem::FieldAccessContract access{
        "T",
        "K",
        {catchem::SemanticAxis::Column, catchem::SemanticAxis::Level, catchem::SemanticAxis::Singleton},
        catchem::PersistencePolicy::Timestep,
        catchem::FieldRequirement::Required,
        catchem::AccessIntent::Read,
        catchem::ExecutionSpaceIntent::Device};
    assert(state.prepare_field_access(access));
    assert(field->host_to_device_sync_count == 1);
    assert(state.prepare_field_access(access));
    assert(field->host_to_device_sync_count == 1);
    return 0;
}
