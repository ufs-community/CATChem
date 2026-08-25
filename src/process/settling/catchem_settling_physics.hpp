#pragma once

#include "catchem_constants.hpp"
#include "catchem_state_manager.hpp"
#include <Kokkos_Core.hpp>

namespace catchem {
    namespace settling {

        struct SettlingFunctor {
            using ExecutionSpace = Kokkos::DefaultExecutionSpace;
            using MemorySpace = ExecutionSpace::memory_space;

            Kokkos::View<double***, Kokkos::LayoutLeft, MemorySpace> conc;
            Kokkos::View<double***, Kokkos::LayoutLeft, MemorySpace> t;
            Kokkos::View<double***, Kokkos::LayoutLeft, MemorySpace> airden;
            Kokkos::View<double***, Kokkos::LayoutLeft, MemorySpace> pedge;
            Kokkos::View<double***, Kokkos::LayoutLeft, MemorySpace> dz;

            Kokkos::View<int*, MemorySpace> aerosol_indices;
            Kokkos::View<double*, MemorySpace> aerosol_radius;
            Kokkos::View<double*, MemorySpace> aerosol_density;

            double cdt;
            int n_levels;

            KOKKOS_INLINE_FUNCTION
            void operator()(const int icol, const int iaero) const {
                int ispec = aerosol_indices(iaero);
                double r = aerosol_radius(iaero);
                double rhop = aerosol_density(iaero);

                if (r <= 0.0)
                    return;

                // Arrays for tau for this column/species
                // Max levels expected is ~256
                const int MAX_LEVELS = 256;
                double tau[MAX_LEVELS];

                double max_tau = 0.0;

                for (int k = 0; k < n_levels; ++k) {
                    double tmpu = t(icol, k, 0);
                    double rhoa = airden(icol, k, 0);

                    // Sutherland's Equation
                    double rmu = 1.8325e-5 * (416.16 / (tmpu + 120.0)) * std::pow(tmpu / 296.16, 1.5);

                    // Thermal velocity
                    double f_vt = 8.0 * 1.3807e-23 / constants::PI / 4.8096e-26;
                    double vt = std::sqrt(tmpu * f_vt);

                    // Mean free path
                    double rmfp = 2.0 * rmu / (rhoa * vt);

                    // Knudsen number
                    double rkn = rmfp / r;

                    // Cunningham slip correction
                    double bpm = 1.0 + 1.246 * rkn;

                    // Fall speed
                    double vsettle = (2.0 / 9.0) * rhop * r * r * constants::G0 * bpm / rmu;

                    // Reynolds number
                    double re = 2.0 * rhoa * r * vsettle / rmu;

                    if (re > 0.01) {
                        double x = std::log(24.0 * re / bpm);
                        double y =
                            -3.18657 +
                            x * (0.992696 +
                                 x * (-1.53193e-3 +
                                      x * (-9.870593e-4 + x * (-5.78878e-4 + x * (8.55176e-5 + -3.27815e-6 * x)))));
                        re = std::exp(y) * bpm;
                        vsettle = 0.5 * rmu * re / (rhoa * r);
                    }

                    double dz_val = dz(icol, k, 0);
                    if (dz_val > 0.0) {
                        tau[k] = vsettle / dz_val;
                    } else {
                        tau[k] = 0.0;
                    }

                    if (tau[k] > max_tau) {
                        max_tau = tau[k];
                    }
                }

                int nSubSteps = 0;
                double dt = cdt;

                if (max_tau > 0.0) {
                    double dt_cfl = 1.0 / max_tau;
                    if (dt_cfl > cdt) {
                        nSubSteps = 1;
                        dt = cdt;
                    } else {
                        nSubSteps = std::ceil(cdt / dt_cfl);
                        dt = cdt / static_cast<double>(nSubSteps);
                    }
                }

                // Local copy of concentration for sub-stepping
                double qa[MAX_LEVELS];
                for (int k = 0; k < n_levels; ++k) {
                    qa[k] = conc(icol, k, ispec);
                }

                for (int iit = 0; iit < nSubSteps; ++iit) {
                    int top = n_levels - 1;
                    qa[top] = qa[top] * (1.0 - dt * tau[top]);

                    for (int k = n_levels - 2; k >= 0; --k) {
                        double delp_k1 = pedge(icol, k + 1, 0) - pedge(icol, k + 2, 0);
                        double delp_k = pedge(icol, k, 0) - pedge(icol, k + 1, 0);
                        // Guard against non-positive layer pressure thickness
                        // (which can occur for degenerate/clamped pressure
                        // profiles) to avoid division by zero producing NaN/Inf.
                        double p_ratio = (delp_k > 0.0) ? (delp_k1 / delp_k) : 0.0;
                        qa[k] = qa[k] + p_ratio * dt * tau[k + 1] * qa[k + 1] - dt * tau[k] * qa[k];
                    }
                }

                for (int k = 0; k < n_levels; ++k) {
                    conc(icol, k, ispec) = qa[k];
                }
            }
        };

    } // namespace settling
} // namespace catchem
