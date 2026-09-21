#include "catchem_runtime_lease.hpp"
#include "catchem_kokkos_compat.hpp"
#include <mutex>
#include <stdexcept>

namespace catchem {
    namespace {
        std::mutex runtime_mutex;
        RuntimeSnapshot runtime_snapshot;

#ifdef CATCHEM_ENABLE_KOKKOS
        // Only the Kokkos branch of the RuntimeLease constructor calls these;
        // the host-only build has no backend to query or start.
        bool backend_initialized() noexcept {
            return Kokkos::is_initialized();
        }
        void initialize_backend() {
            Kokkos::initialize();
        }
#endif
        void finalize_backend() noexcept {
#ifdef CATCHEM_ENABLE_KOKKOS
            if (Kokkos::is_initialized() && !Kokkos::is_finalized())
                Kokkos::finalize();
#endif
        }
    } // namespace

    RuntimeLease::RuntimeLease(RuntimeMode requested) {
        std::lock_guard<std::mutex> lock(runtime_mutex);
#ifndef CATCHEM_ENABLE_KOKKOS
        (void)requested;
        mode_ = RuntimeMode::Disabled;
        runtime_snapshot.mode = RuntimeMode::Disabled;
        runtime_snapshot.state = RuntimeState::Initialized;
#else
        if (runtime_snapshot.state == RuntimeState::Finalized)
            throw std::runtime_error("CATChem execution runtime has already been finalized");
        if (runtime_snapshot.lease_count == 0 && runtime_snapshot.state == RuntimeState::Unavailable) {
            if (requested == RuntimeMode::Disabled)
                mode_ = RuntimeMode::Disabled;
            else if (backend_initialized())
                mode_ = RuntimeMode::HostOwned;
            else if (requested == RuntimeMode::CATChemOwned) {
                initialize_backend();
                mode_ = RuntimeMode::CATChemOwned;
                runtime_snapshot.initialized_by_catchem = true;
            } else
                throw std::runtime_error("HostOwned runtime requested before the host initialized Kokkos");
            runtime_snapshot.mode = mode_;
            runtime_snapshot.state = RuntimeState::Initialized;
        } else {
            mode_ = runtime_snapshot.mode;
            if (requested == RuntimeMode::Disabled && mode_ != RuntimeMode::Disabled)
                throw std::runtime_error("Disabled runtime mode conflicts with active runtime leases");
        }
#endif
        ++runtime_snapshot.lease_count;
        active_ = true;
    }

    RuntimeLease::RuntimeLease(RuntimeLease&& other) noexcept : mode_(other.mode_), active_(other.active_) {
        other.active_ = false;
    }
    RuntimeLease& RuntimeLease::operator=(RuntimeLease&& other) noexcept {
        if (this != &other) {
            release();
            mode_ = other.mode_;
            active_ = other.active_;
            other.active_ = false;
        }
        return *this;
    }
    RuntimeLease::~RuntimeLease() {
        release();
    }

    void RuntimeLease::release() noexcept {
        if (!active_)
            return;
        std::lock_guard<std::mutex> lock(runtime_mutex);
        if (runtime_snapshot.lease_count > 0)
            --runtime_snapshot.lease_count;
        active_ = false;
        if (runtime_snapshot.lease_count == 0 && runtime_snapshot.mode == RuntimeMode::HostOwned) {
            runtime_snapshot = {};
        }
    }

    RuntimeSnapshot RuntimeLease::snapshot() noexcept {
        std::lock_guard<std::mutex> lock(runtime_mutex);
        return runtime_snapshot;
    }
    BoundaryStatus RuntimeLease::request_finalize() noexcept {
        std::lock_guard<std::mutex> lock(runtime_mutex);
        if (runtime_snapshot.mode != RuntimeMode::CATChemOwned || runtime_snapshot.lease_count != 0)
            return BoundaryStatus::InvalidState;
        finalize_backend();
        runtime_snapshot.state = RuntimeState::Finalized;
        return BoundaryStatus::Success;
    }
} // namespace catchem
