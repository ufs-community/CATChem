#include "catchem_core.hpp"
#include "catchem_interop_field.hpp"
#include "catchem_kokkos_compat.hpp"
#include "catchem_state_manager.hpp"
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>

// Indexing-contract regression lock (spec 011 US5).  The CAP -> core ->
// science-bridge handoff passes flat column-major buffers that the Fortran
// side reinterprets with c_f_pointer(arr, [n_cols, n_lev(, n_species)]).  This
// test proves the C++ mdspan view agrees with that flat ordering:
//   element(i, j, k) == flat[(k-1)*n_lev*n_cols + (j-1)*n_cols + (i-1)]
// (1-based Fortran form; 0-based here), and that interface fields (Z, PEDGE)
// carry n_levels+1 in the vertical axis.

namespace {

    int failures = 0;
    void check(bool condition, const std::string& label) {
        std::cout << (condition ? "  PASS: " : "  FAIL: ") << label << '\n';
        if (!condition)
            ++failures;
    }

    // Encode a unique value per (column, level) so any index transposition is
    // detectable: value = 1000*level + column.
    double tag(int column, int level) {
        return 1000.0 * level + column;
    }

} // namespace

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "==========================================" << std::endl;
        std::cout << "RUNNING TEST: Interop Indexing Contract" << std::endl;
        std::cout << "==========================================" << std::endl;

        const int n_cols = 3;
        const int n_levels = 4;
        const int n_species = 5;

        auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
        auto state = core->get_state_manager();

        // --- 3D met field: level axis is the stride-n_cols axis --------------
        {
            std::vector<double> buffer(static_cast<size_t>(n_cols) * n_levels, 0.0);
            for (int c = 0; c < n_cols; ++c)
                for (int lev = 0; lev < n_levels; ++lev)
                    buffer[static_cast<size_t>(lev) * n_cols + c] = tag(c, lev);
            state->bind_met_field_3d("T", buffer.data());

            auto field = state->find_field<3>("T");
            check(field != nullptr, "3D met field resolvable");
            if (field) {
                auto view = field->mdspan();
                bool ok = true;
                for (int c = 0; c < n_cols && ok; ++c)
                    for (int lev = 0; lev < n_levels && ok; ++lev) {
                        std::size_t flat = static_cast<std::size_t>(lev) * n_cols + c;
                        ok = (view(c, lev, 0) == buffer[flat]) && (view(c, lev, 0) == tag(c, lev));
                    }
                check(ok, "3D mdspan(col,lev) == flat[lev*n_cols+col] (column-major, bottom-to-top)");
                check(field->extent(0) == static_cast<std::size_t>(n_cols) &&
                          field->extent(1) == static_cast<std::size_t>(n_levels),
                      "3D met extents are (n_cols, n_levels)");
            }
        }

        // --- Interface field (Z): vertical axis is n_levels+1 ----------------
        {
            const int n_interface = n_levels + 1;
            std::vector<double> zbuffer(static_cast<size_t>(n_cols) * n_interface, 0.0);
            for (int c = 0; c < n_cols; ++c)
                for (int lev = 0; lev < n_interface; ++lev)
                    zbuffer[static_cast<size_t>(lev) * n_cols + c] = tag(c, lev);
            state->bind_met_field_3d("Z", zbuffer.data());

            auto field = state->find_field<3>("Z");
            check(field != nullptr, "interface field resolvable");
            if (field) {
                auto view = field->mdspan();
                bool ok = true;
                for (int c = 0; c < n_cols && ok; ++c)
                    for (int lev = 0; lev < n_interface && ok; ++lev) {
                        std::size_t flat = static_cast<std::size_t>(lev) * n_cols + c;
                        ok = (view(c, lev, 0) == zbuffer[flat]);
                    }
                check(ok, "interface mdspan(col,lev) == flat[lev*n_cols+col]");
                check(field->extent(1) == static_cast<std::size_t>(n_levels + 1),
                      "interface field vertical extent is n_levels+1");
            }
        }

        // --- Species slab: species axis is the outermost stride --------------
        {
            std::vector<double> chem(static_cast<size_t>(n_cols) * n_levels * n_species, 0.0);
            for (int c = 0; c < n_cols; ++c)
                for (int lev = 0; lev < n_levels; ++lev)
                    for (int s = 0; s < n_species; ++s) {
                        std::size_t flat = static_cast<std::size_t>(s) * n_levels * n_cols +
                                           static_cast<std::size_t>(lev) * n_cols + c;
                        chem[flat] = 1000000.0 * s + tag(c, lev);
                    }
            state->bind_unified_chemistry(chem.data());

            auto conc = state->chemistry().conc;
            check(conc != nullptr, "chemistry concentration field resolvable");
            if (conc) {
                auto view = conc->mdspan();
                bool ok = true;
                for (int c = 0; c < n_cols && ok; ++c)
                    for (int lev = 0; lev < n_levels && ok; ++lev)
                        for (int s = 0; s < n_species && ok; ++s) {
                            std::size_t flat = static_cast<std::size_t>(s) * n_levels * n_cols +
                                               static_cast<std::size_t>(lev) * n_cols + c;
                            ok = (view(c, lev, s) == chem[flat]);
                        }
                check(ok, "species mdspan(col,lev,sp) == flat[sp*n_lev*n_cols + lev*n_cols + col]");
                check(conc->extent(2) == static_cast<std::size_t>(n_species), "species axis extent is n_species");
            }
        }

        std::cout << (failures == 0 ? "SUCCESS: indexing contract holds.\n"
                                    : "FAILURE: " + std::to_string(failures) + " index assertion(s) failed.\n");
    }
    Kokkos::finalize();
    return failures == 0 ? 0 : 1;
}
