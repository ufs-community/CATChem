#pragma once

#include "catchem_error.hpp"
#include <cstddef>
#include <memory>

namespace catchem {

    enum class RuntimeMode { Disabled, HostOwned, CATChemOwned };
    enum class RuntimeState { Unavailable, Initialized, FinalizePending, Finalized };

    struct RuntimeSnapshot {
        RuntimeMode mode = RuntimeMode::Disabled;
        RuntimeState state = RuntimeState::Unavailable;
        std::size_t lease_count = 0;
        bool initialized_by_catchem = false;
    };

    class RuntimeLease {
    public:
        RuntimeLease() noexcept = default;
        explicit RuntimeLease(RuntimeMode mode);
        RuntimeLease(const RuntimeLease&) = delete;
        RuntimeLease& operator=(const RuntimeLease&) = delete;
        RuntimeLease(RuntimeLease&& other) noexcept;
        RuntimeLease& operator=(RuntimeLease&& other) noexcept;
        ~RuntimeLease();

        bool active() const noexcept { return active_; }
        RuntimeMode mode() const noexcept { return mode_; }
        void release() noexcept;

        static RuntimeSnapshot snapshot() noexcept;
        static BoundaryStatus request_finalize() noexcept;

    private:
        RuntimeMode mode_ = RuntimeMode::Disabled;
        bool active_ = false;
    };

} // namespace catchem
