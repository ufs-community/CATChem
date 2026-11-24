#include "catchem_api.hpp"
#include "catchem_interop_field.hpp"
#include "catchem_runtime_lease.hpp"
#include <cassert>
#include <vector>

int main() {
    for (int cycle = 0; cycle < 100; ++cycle) {
        void* first = catchem_core_create(2, 3, 1);
        void* second = catchem_core_create(2, 3, 1);
        assert(first && second);
        void* child = catchem_core_get_state_manager(first);
        assert(child);
        assert(catchem::RuntimeLease::snapshot().lease_count == 2);
        assert(catchem_core_destroy_checked(first) == CATCHEM_SUCCESS);
        assert(catchem_state_begin_import_generation(child) == CATCHEM_INVALID_HANDLE);
        assert(catchem::RuntimeLease::snapshot().lease_count == 1);
        assert(catchem_core_destroy_checked(second) == CATCHEM_SUCCESS);
        assert(catchem::RuntimeLease::snapshot().lease_count == 0);
    }
#ifdef CATCHEM_ENABLE_KOKKOS
    assert(catchem::RuntimeLease::request_finalize() == catchem::BoundaryStatus::Success);
#endif

    std::vector<double> values(6, 1.0);
    catchem::InteropField<double, 2> field(values.data(), {2, 3});
    assert(field.latest_writer == catchem::LatestWriter::HostCurrent);
    assert(field.host_read() == values.data());
    assert(field.latest_writer == catchem::LatestWriter::HostCurrent);
    assert(field.host_write() == values.data());
    assert(field.latest_writer == catchem::LatestWriter::HostCurrent);
    (void)field.device_read();
    assert(field.latest_writer == catchem::LatestWriter::Synchronized);
    (void)field.device_write();
    assert(field.latest_writer == catchem::LatestWriter::DeviceCurrent);
    (void)field.host_read();
    assert(field.latest_writer == catchem::LatestWriter::Synchronized);
    return 0;
}
