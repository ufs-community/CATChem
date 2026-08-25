#pragma once
#include "catchem_process_interface.hpp"
#include <Kokkos_Core.hpp>
#include <functional>

namespace catchem {

    class SettlingProcess : public ProcessInterface {
    private:
        std::string active_scheme;
        std::function<void(void*)> fortran_callback;

        // Device Views for aerosol properties
        Kokkos::View<int*, Kokkos::DefaultExecutionSpace::memory_space> dev_aero_indices;
        Kokkos::View<double*, Kokkos::DefaultExecutionSpace::memory_space> dev_radius_dry;
        Kokkos::View<double*, Kokkos::DefaultExecutionSpace::memory_space> dev_rhop_dry;

    public:
        SettlingProcess();
        std::string get_name() const override { return "settling"; }
        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override;

        // For legacy tests
        void set_fortran_bridge_callback(std::function<void(void*)> cb);
    };

} // namespace catchem
