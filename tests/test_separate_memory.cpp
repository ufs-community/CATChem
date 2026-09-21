#include "catchem_interop_field.hpp"
#include "catchem_state_manager.hpp"
#include <cassert>
#include <vector>

int main() {
    std::vector<double> host(6, 1.0);
    catchem::InteropField<double, 2> field(host.data(), {2, 3});
    (void)field.device_read();
    assert(field.latest_writer == catchem::LatestWriter::Synchronized);
    auto device = field.device_write();
    device(0, 0) = 42.0;
    assert(field.latest_writer == catchem::LatestWriter::DeviceCurrent);
    const double* current_host = field.host_read();
    assert(current_host[0] == 42.0);
    assert(field.latest_writer == catchem::LatestWriter::Synchronized);
    current_host = field.host_read();
    assert(current_host[0] == 42.0);
    assert(field.latest_writer == catchem::LatestWriter::Synchronized);

    // Schedule-directed access moves only the declared field and does not
    // repeat a transfer while the requested execution space remains current.
    catchem::StateManager state(2, 3, 1);
    std::vector<double> temperature(6, 280.0);
    std::vector<double> pressure(6, 90000.0);
    state.bind_met_field_3d("T", temperature.data());
    state.bind_met_field_3d("PMID", pressure.data());
    const catchem::FieldAccessContract temperature_on_device{
        "T",
        "K",
        {catchem::SemanticAxis::Column, catchem::SemanticAxis::Level, catchem::SemanticAxis::Singleton},
        catchem::PersistencePolicy::Timestep,
        catchem::FieldRequirement::Required,
        catchem::AccessIntent::Read,
        catchem::ExecutionSpaceIntent::Device};
    assert(state.prepare_field_access(temperature_on_device));
    assert(state.meteorology().T->host_to_device_sync_count == 1);
    assert(state.meteorology().PMID->host_to_device_sync_count == 0);
    assert(state.prepare_field_access(temperature_on_device));
    assert(state.meteorology().T->host_to_device_sync_count == 1);
    return 0;
}
