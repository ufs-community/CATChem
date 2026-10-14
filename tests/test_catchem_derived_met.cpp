#include "catchem_core.hpp"
#include "catchem_kokkos_compat.hpp"
#include "catchem_state_manager.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
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

    {
        std::cout << "==========================================" << std::endl;
        std::cout << "RUNNING TEST: Derived MET (RH/AIRDEN/OBK) Unit Test" << std::endl;
        std::cout << "==========================================" << std::endl;

        const int n_cols = 3;
        const int n_levels = 4;
        const int n_species = 1;
        const int size_3d = n_cols * n_levels;

        // Reference formulas from contracts/derived-met-definitions.md.
        // Magnus/Alduchov-Eskridge saturation vapor pressure [Pa].
        auto es_ref = [](double T) {
            const double tc = T - 273.15;
            return 610.94 * std::exp(17.625 * tc / (tc + 243.04));
        };
        auto clip = [](double x, double lo, double hi) { return x < lo ? lo : (x > hi ? hi : x); };
        // Derived values are computed in catchem::fp (float unless USE_REAL8),
        // so compare against the double reference with a relative tolerance the
        // same way the property tests do.
        const double rel_tol = 32.0 * static_cast<double>(std::numeric_limits<catchem::fp>::epsilon());
        auto close = [&](double actual, double expected) {
            return std::abs(actual - expected) <= rel_tol * std::abs(expected);
        };
        // Moist air density factor AIRMW/H2OMW - 1 from the authoritative constants.
        const double moist_factor =
            static_cast<double>(catchem::constants::AIR_MW) / static_cast<double>(catchem::constants::H2O_MW) - 1.0;
        auto airden_moist_ref = [&](double p, double T, double qv) {
            return p / (static_cast<double>(catchem::constants::RD) * T * (1.0 + moist_factor * qv));
        };

        // --- (a) RH equals e/es with the Magnus form, clamped [0.005,0.99] ----
        // M-1..M-3. qv=0 hits the floor, saturated hits the ceiling, a mid-range
        // value matches the reference exactly.
        {
            auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
            auto state = core->get_state_manager();
            std::vector<double> t(size_3d), qv(size_3d), pmid(size_3d);
            for (int c = 0; c < n_cols; ++c) {
                for (int lev = 0; lev < n_levels; ++lev) {
                    const std::size_t i = flat(n_cols, c, lev);
                    t[i] = 290.0 + 2.0 * lev;
                    pmid[i] = 95000.0 - 15000.0 * lev;
                }
            }
            // column 0: dry (qv=0) -> floor; column 1: mid; column 2: saturated -> ceiling
            for (int lev = 0; lev < n_levels; ++lev)
                qv[flat(n_cols, 0, lev)] = 0.0;
            for (int lev = 0; lev < n_levels; ++lev)
                qv[flat(n_cols, 1, lev)] = 0.008;
            // qv large enough that e/es > 0.99 at every level, including the
            // coolest-highest-ratio top layer -> the whole column clamps to 0.99.
            for (int lev = 0; lev < n_levels; ++lev)
                qv[flat(n_cols, 2, lev)] = 0.050;
            state->bind_met_field_3d("T", t.data());
            state->bind_met_field_3d("QV", qv.data());
            state->bind_met_field_3d("PMID", pmid.data());

            state->derive_relative_humidity();
            const double* rh = state->read_field<3>("RH");
            check(rh != nullptr, "RH derived and readable");
            if (rh) {
                bool floor_ok = true;
                for (int lev = 0; lev < n_levels; ++lev)
                    floor_ok = floor_ok && close(rh[flat(n_cols, 0, lev)], 0.005);
                check(floor_ok, "qv=0 clamps RH to the 0.005 floor (M-2)");

                bool ceiling_ok = true;
                for (int lev = 0; lev < n_levels; ++lev)
                    ceiling_ok = ceiling_ok && close(rh[flat(n_cols, 2, lev)], 0.99);
                check(ceiling_ok, "supersaturated clamps RH to the 0.99 ceiling (M-2)");

                bool mid_ok = true;
                for (int lev = 0; lev < n_levels && mid_ok; ++lev) {
                    const std::size_t i = flat(n_cols, 1, lev);
                    const double e = pmid[i] * qv[i] / (0.622 + qv[i]);
                    const double expected = clip(e / es_ref(t[i]), 0.005, 0.99);
                    mid_ok = close(rh[i], expected);
                }
                check(mid_ok, "mid-range RH equals e/es with Magnus es (M-1)");
            }
        }

        // --- (b) moist AIRDEN, MAIRDEN resolves equal, differs from dry -------
        // M-5..M-7.
        {
            auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
            auto state = core->get_state_manager();
            std::vector<double> t(size_3d), qv(size_3d), pmid(size_3d);
            for (int c = 0; c < n_cols; ++c)
                for (int lev = 0; lev < n_levels; ++lev) {
                    const std::size_t i = flat(n_cols, c, lev);
                    t[i] = 288.0 + lev;
                    pmid[i] = 100000.0 - 10000.0 * lev;
                    qv[i] = 0.01; // nonzero moisture everywhere
                }
            state->bind_met_field_3d("T", t.data());
            state->bind_met_field_3d("QV", qv.data());
            state->bind_met_field_3d("PMID", pmid.data());

            state->derive_airden();
            const double* airden = state->read_field<3>("AIRDEN");
            check(airden != nullptr, "AIRDEN derived and readable");
            if (airden) {
                bool moist_ok = true;
                for (int i = 0; i < size_3d && moist_ok; ++i)
                    moist_ok = close(airden[i], airden_moist_ref(pmid[i], t[i], qv[i]));
                check(moist_ok, "AIRDEN uses the moist formula P/(RD*T*(1+factor*qv)) (M-5)");

                // MAIRDEN resolves to the same storage as AIRDEN.
                const double* mairden = state->read_field<3>("MAIRDEN");
                bool same = mairden != nullptr;
                for (int i = 0; i < size_3d && same; ++i)
                    same = mairden[i] == airden[i];
                check(same, "MAIRDEN resolves equal to moist AIRDEN (M-7)");

                // AIRDEN_DRY stays on its own (dry, humidity-corrected) formula
                // and differs from the moist value whenever qv>0 (M-6).
                state->derive_airden_dry();
                const double* dry = state->read_field<3>("AIRDEN_DRY");
                bool dry_differs = dry != nullptr;
                for (int i = 0; i < size_3d && dry_differs; ++i)
                    dry_differs = dry[i] != airden[i];
                check(dry_differs, "AIRDEN_DRY stays dry and differs from moist AIRDEN when qv>0 (M-6)");
            }
        }

        // --- (c) OBK uses lowest-layer T and moist density -------------------
        // M-8: OBK changes when T(:,:,1) (level 0) differs from TS, and the
        // density argument is the moist AIRDEN rather than the dry P/(RD*T).
        {
            auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
            auto state = core->get_state_manager();
            std::vector<double> t(size_3d), qv(size_3d), pmid(size_3d);
            std::vector<double> ustar(n_cols), hflux(n_cols), ts(n_cols);
            for (int c = 0; c < n_cols; ++c) {
                ustar[c] = 0.35;
                hflux[c] = 25.0; // positive -> unstable -> negative L
                ts[c] = 300.0;   // skin temperature distinct from lowest layer
                for (int lev = 0; lev < n_levels; ++lev) {
                    const std::size_t i = flat(n_cols, c, lev);
                    t[i] = 293.0 + lev; // lowest layer (lev 0) = 293 K, not 300 K
                    pmid[i] = 101000.0 - 9000.0 * lev;
                    qv[i] = 0.012;
                }
            }
            state->bind_met_field_3d("T", t.data());
            state->bind_met_field_3d("QV", qv.data());
            state->bind_met_field_3d("PMID", pmid.data());
            state->bind_met_field_2d("USTAR", ustar.data());
            state->bind_met_field_2d("HFLUX", hflux.data());
            state->bind_met_field_2d("TS", ts.data());

            state->derive_airden();
            state->derive_obk();
            const double* obk = state->read_field<2>("OBK");
            check(obk != nullptr, "OBK derived and readable");
            if (obk) {
                const double rho = airden_moist_ref(pmid[0], t[0], qv[0]);
                const double cp = static_cast<double>(catchem::constants::CP);
                const double g0 = static_cast<double>(catchem::constants::G0);
                const double expected = -(ustar[0] * ustar[0] * ustar[0] * rho * cp * t[0]) / (0.41 * g0 * hflux[0]);
                // Uses lowest-layer T (293 K) and moist density, not TS (300 K).
                // The TS value differs by ~2.4 % (7/293), far outside the float
                // tolerance, so matching `expected` and rejecting `wrong_ts`
                // together pin both the temperature source and the density.
                const double wrong_ts = -(ustar[0] * ustar[0] * ustar[0] * rho * cp * ts[0]) / (0.41 * g0 * hflux[0]);
                bool ok = close(obk[0], expected) && !close(obk[0], wrong_ts);
                check(ok, "OBK uses lowest-layer T and moist AIRDEN (M-8)");
            }
        }

        // --- (d) host precedence: current RH/AIRDEN/OBK are left untouched ----
        // FR-014, M-4/M-9.
        {
            auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
            auto state = core->get_state_manager();
            std::vector<double> t(size_3d, 288.0), qv(size_3d, 0.01), pmid(size_3d, 100000.0);
            std::vector<double> host_rh(size_3d, 0.55);
            std::vector<double> host_airden(size_3d, 1.111);
            std::vector<double> host_obk(n_cols, -42.0);
            std::vector<double> ustar(n_cols, 0.35), hflux(n_cols, 25.0), ts(n_cols, 300.0);
            state->bind_met_field_3d("T", t.data());
            state->bind_met_field_3d("QV", qv.data());
            state->bind_met_field_3d("PMID", pmid.data());
            state->bind_met_field_3d("RH", host_rh.data());
            state->bind_met_field_3d("AIRDEN", host_airden.data());
            state->bind_met_field_2d("OBK", host_obk.data());
            state->bind_met_field_2d("USTAR", ustar.data());
            state->bind_met_field_2d("HFLUX", hflux.data());
            state->bind_met_field_2d("TS", ts.data());

            state->derive_relative_humidity();
            state->derive_airden();
            state->derive_obk();
            const double* rh = state->read_field<3>("RH");
            const double* airden = state->read_field<3>("AIRDEN");
            const double* obk = state->read_field<2>("OBK");
            bool rh_kept = rh != nullptr;
            for (int i = 0; i < size_3d && rh_kept; ++i)
                rh_kept = rh[i] == 0.55;
            check(rh_kept, "host RH is preserved untouched by the derive (M-4)");
            bool airden_kept = airden != nullptr;
            for (int i = 0; i < size_3d && airden_kept; ++i)
                airden_kept = airden[i] == 1.111;
            check(airden_kept, "host AIRDEN is preserved untouched by the derive (FR-014)");
            bool obk_kept = obk != nullptr;
            for (int c = 0; c < n_cols && obk_kept; ++c)
                obk_kept = obk[c] == -42.0;
            check(obk_kept, "host OBK is preserved untouched by the derive (M-9)");
        }

        // --- (e) missing lowest-layer T reports the missing input -------------
        // OBK must not silently substitute TS when the profile T is unavailable.
        {
            auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
            auto state = core->get_state_manager();
            std::vector<double> ustar(n_cols, 0.35), hflux(n_cols, 25.0), ts(n_cols, 300.0);
            state->bind_met_field_2d("USTAR", ustar.data());
            state->bind_met_field_2d("HFLUX", hflux.data());
            state->bind_met_field_2d("TS", ts.data());
            // PMID current but T never bound -> derive_obk must name the missing T.
            std::vector<double> pmid(size_3d, 100000.0);
            state->bind_met_field_3d("PMID", pmid.data());

            bool threw = false;
            std::string message;
            try {
                state->derive_obk();
            } catch (const std::runtime_error& ex) {
                threw = true;
                message = ex.what();
            }
            check(threw, "OBK derive throws when lowest-layer T is missing");
            check(threw && message.find("T") != std::string::npos,
                  "OBK error names the missing T rather than substituting TS silently");
        }

        std::cout << (failures == 0 ? "SUCCESS: all derived-met assertions passed.\n"
                                    : "FAILURE: " + std::to_string(failures) + " derived-met assertion(s) failed.\n");
    }
    Kokkos::finalize();
    return failures == 0 ? 0 : 1;
}
