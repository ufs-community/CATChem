#include "catchem_process_settling.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_settling_physics.hpp"
#include <iostream>

namespace catchem {

    SettlingProcess::SettlingProcess() : active_scheme("c++_kokkos"), fortran_callback(nullptr) {}

    void SettlingProcess::init(std::shared_ptr<StateManager> state) {
        int num_aerosols = state->chem.aerosol_indices.size();
        if (num_aerosols > 0) {
            dev_aero_indices =
                Kokkos::View<int*, Kokkos::DefaultExecutionSpace::memory_space>("dev_aero_indices", num_aerosols);
            dev_radius_dry =
                Kokkos::View<double*, Kokkos::DefaultExecutionSpace::memory_space>("dev_radius_dry", num_aerosols);
            dev_rhop_dry =
                Kokkos::View<double*, Kokkos::DefaultExecutionSpace::memory_space>("dev_rhop_dry", num_aerosols);

            auto host_aero_indices = Kokkos::create_mirror_view(dev_aero_indices);
            auto host_radius_dry = Kokkos::create_mirror_view(dev_radius_dry);
            auto host_rhop_dry = Kokkos::create_mirror_view(dev_rhop_dry);

            for (int i = 0; i < num_aerosols; ++i) {
                int ispec = state->chem.aerosol_indices[i];
                host_aero_indices(i) = ispec;
                host_radius_dry(i) = state->chem.species_list[ispec].radius * 1e-6; // Convert microns to meters
                host_rhop_dry(i) = state->chem.species_list[ispec].density;
                std::cout << "DEBUG INIT: Aerosol " << i << ": species index=" << ispec
                          << ", name=" << state->chem.species_list[ispec].short_name
                          << ", radius_dry=" << host_radius_dry(i) << ", density=" << host_rhop_dry(i) << std::endl;
            }

            Kokkos::deep_copy(dev_aero_indices, host_aero_indices);
            Kokkos::deep_copy(dev_radius_dry, host_radius_dry);
            Kokkos::deep_copy(dev_rhop_dry, host_rhop_dry);
        }
    }

    void SettlingProcess::set_fortran_bridge_callback(std::function<void(void*)> cb) {
        fortran_callback = cb;
    }

    void SettlingProcess::run(std::shared_ptr<StateManager> state) {
        if (fortran_callback) {
            // Fallback for tests explicitly requesting the Fortran bridge
            state->sync_to_host();
            fortran_callback(static_cast<void*>(state.get()));
            state->sync_to_device();
            return;
        }

        int num_aerosols = state->chem.aerosol_indices.size();
        if (num_aerosols == 0)
            return;

        if (!state->met.T || !state->met.AIRDEN || !state->met.PEDGE || !state->met.BXHEIGHT || !state->chem.conc) {
            std::cerr << "SettlingProcess: Missing required views!\n";
            return;
        }

        std::cout << "SettlingProcess: Views exist. Launching kernel for " << num_aerosols << " aerosols.\n";

        settling::SettlingFunctor functor;
        functor.conc = state->chem.conc->view();
        functor.t = state->met.T->view();
        functor.airden = state->met.AIRDEN->view();
        functor.pedge = state->met.PEDGE->view();
        functor.dz = state->met.BXHEIGHT->view();

        functor.aerosol_indices = dev_aero_indices;
        functor.aerosol_radius = dev_radius_dry;
        functor.aerosol_density = dev_rhop_dry;

        functor.cdt = state->time.timestep;
        functor.n_levels = state->n_levels;

        Kokkos::parallel_for("settling_compute_c++",
                             Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>(
                                 {0, 0}, {state->n_cols, num_aerosols}),
                             functor);
        Kokkos::fence();
        std::cout << "SettlingProcess: Kernel complete.\n";
    }

    void SettlingProcess::finalize() {
        // Kokkos views will be deallocated automatically when their reference count goes to zero.
    }

} // namespace catchem

extern "C" {
void catchem_register_settling_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "settling", []() { return std::make_shared<catchem::SettlingProcess>(); });
}
}
