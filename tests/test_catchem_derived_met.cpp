#include "catchem_core.hpp"
#include "catchem_kokkos_compat.hpp"
#include "catchem_state_manager.hpp"
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace {

    // CLDFRC parity: the upstream metstate_mod derives the column cloud
    // fraction as CLDFRC(:,:) = SUM(CLDF, DIM=3) -- a vertical sum over all
    // layers, not just the surface layer.

    enum FieldFlag : unsigned {
        F_CLDF = 1u << 0,
        F_CLDFRC_HOST = 1u << 1,
    };

    // Flat index into a bound (n_cols, n_levels) column-major buffer.
    std::size_t flat(int n_cols, int column, int level) {
        return static_cast<std::size_t>(column) + static_cast<std::size_t>(level) * n_cols;
    }

} // namespace

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    int failures = 0;
    auto check = [&](bool condition, const std::string& label) {
        std::cout << (condition ? "  PASS: " : "  FAIL: ") << label << '\n';
        if (!condition)
            ++failures;
    };

    {
        std::cout << "==========================================" << std::endl;
        std::cout << "RUNNING TEST: Derived MET (CLDFRC) Unit Test" << std::endl;
        std::cout << "==========================================" << std::endl;

        const int n_cols = 4;
        const int n_levels = 5;
        const int n_species = 1;

        // --- CLDFRC == vertical sum of CLDF per column ------------------------
        {
            auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
            auto state = core->get_state_manager();
            std::vector<double> cldf(n_cols * n_levels, 0.0);
            for (int c = 0; c < n_cols; ++c)
                for (int lev = 0; lev < n_levels; ++lev)
                    cldf[flat(n_cols, c, lev)] = 0.1 + 0.01 * lev;
            state->bind_met_field_3d("CLDF", cldf.data());

            state->derive_column_cloud_fraction();
            const double* cldfrc = state->read_field<2>("CLDFRC");
            check(cldfrc != nullptr, "CLDFRC derived and readable");
            // sum over lev=0..4 of (0.1 + 0.01*lev) = 0.5 + 0.1 = 0.6
            bool sum_ok = cldfrc != nullptr;
            for (int c = 0; c < n_cols && sum_ok; ++c)
                sum_ok = std::abs(cldfrc[c] - 0.6) < 1.0e-12;
            check(sum_ok, "CLDFRC equals the column-summed CLDF (0.6)");
        }

        // --- Layers aloft contribute to the column total ----------------------
        {
            auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
            auto state = core->get_state_manager();
            std::vector<double> cldf(n_cols * n_levels, 0.5); // every layer 0.5
            for (int c = 0; c < n_cols; ++c)
                cldf[flat(n_cols, c, 0)] = 0.2; // distinct surface value
            state->bind_met_field_3d("CLDF", cldf.data());

            state->derive_column_cloud_fraction();
            const double* cldfrc = state->read_field<2>("CLDFRC");
            // 0.2 (surface) + 4 x 0.5 (aloft) = 2.2; upstream applies no clamp.
            bool column_sum = cldfrc != nullptr;
            for (int c = 0; c < n_cols && column_sum; ++c)
                column_sum = std::abs(cldfrc[c] - 2.2) < 1.0e-12;
            check(column_sum, "CLDFRC sums the surface layer with layers aloft (2.2)");
        }

        // --- Host-provided CLDFRC is preserved untouched ---------------------
        {
            auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
            auto state = core->get_state_manager();
            std::vector<double> cldf(n_cols * n_levels, 0.5);
            std::vector<double> host_cldfrc(n_cols, 0.99);
            state->bind_met_field_3d("CLDF", cldf.data());
            state->bind_met_field_2d("CLDFRC", host_cldfrc.data());

            state->derive_column_cloud_fraction();
            const double* cldfrc = state->read_field<2>("CLDFRC");
            bool preserved = cldfrc != nullptr;
            for (int c = 0; c < n_cols && preserved; ++c)
                preserved = std::abs(cldfrc[c] - 0.99) < 1.0e-12;
            check(preserved, "host-provided CLDFRC is not overwritten by derivation");
        }

        std::cout << (failures == 0 ? "SUCCESS: all derived-met assertions passed.\n"
                                    : "FAILURE: " + std::to_string(failures) + " derived-met assertion(s) failed.\n");
    }
    Kokkos::finalize();
    return failures == 0 ? 0 : 1;
}
