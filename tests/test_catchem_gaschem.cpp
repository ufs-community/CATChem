#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_process_gaschem.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_state_manager.hpp"
#include "catchem_test_config.hpp"
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

extern "C" {
void catchem_register_photolysis_cpp();
void catchem_register_gaschem_cpp();
}

bool file_exists(const std::string& name) {
    return std::filesystem::exists(name);
}

// Runs the full coupled photolysis+gaschem scenario once for a given gaschem
// mechanism configuration (a v0 CAMP directory, a v1 single file, or a v1
// directory relying on config.yaml auto-detection).
static void run_coupled_case(const std::string& label, const std::string& config_key, const std::string& config_value) {
    std::cout << "\n------------------------------------------" << std::endl;
    std::cout << "CASE: gaschem mechanism via " << config_key << " (" << label << ")" << std::endl;
    std::cout << "------------------------------------------\n" << std::endl;

    // 2. Set up core and states
    int n_cols = 1;
    int n_levels = 3;
    int n_species = 4; // matches tests/fixtures/mechanisms/chapman.yml

    auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
    auto state = core->get_state_manager();

    // Time setup (Noon, summer)
    state->clock().year = 2026;
    state->clock().month = 7;
    state->clock().day = 13;
    state->clock().hour = 12;
    state->clock().minute = 0;
    state->clock().second = 0;
    state->clock().calculate_derived_fields();

    // Meteorological arrays
    std::vector<double> lat(n_cols, 40.0);
    std::vector<double> lon(n_cols, -105.0);
    std::vector<double> temperature(n_cols * n_levels, 280.0);
    std::vector<double> airden(n_cols * n_levels, 1.2);
    std::vector<double> pedge(n_cols * (n_levels + 1), 101300.0);
    std::vector<double> bxheight(n_cols * n_levels, 100.0);

    for (int i = 0; i < n_levels; ++i) {
        temperature[i] = 280.0 - 0.5 * i;
        airden[i] = 1.2 * std::exp(-i / 10.0);
        bxheight[i] = 1000.0;
        pedge[i] = 101300.0 * std::exp(-i / 10.0);
    }
    pedge[n_levels] = 101300.0 * std::exp(-n_levels / 10.0);

    state->bind_met_field_2d("LAT", lat.data());
    state->bind_met_field_2d("LON", lon.data());
    state->bind_met_field_3d("T", temperature.data());
    state->bind_met_field_3d("AIRDEN", airden.data());
    state->bind_met_field_3d("AIRDEN_DRY", airden.data());
    state->bind_met_field_3d("PEDGE", pedge.data());
    state->bind_met_field_3d("PMID", pedge.data()); // PMID maps to PMID in tests
    state->bind_met_field_3d("BXHEIGHT", bxheight.data());

    // Load a mechanism whose species mirror chapman v0 (so every CATChem
    // species maps into the MICM variable map) and that carries the
    // photolysis.ozone role the photolysis contract requires.
    std::string species_config = std::string(catchem::test::TEST_DIR) + "/fixtures/mechanisms/chapman.yml";
    assert(file_exists(species_config) && "ERROR: Could not find the chapman mechanism fixture!");
    state->load_species_config(species_config);
    std::vector<double> conc_data(n_cols * n_levels * n_species, 1.0); // 1.0 ppmv initially
    state->bind_unified_chemistry(conc_data.data());

    // 3. Resolve configs explicitly using configured header
    std::string photolysis_config =
        std::string(catchem::test::SOURCE_DIR) + "/src/external/musica/configs/tuvx/from_host/config.json";
    assert(file_exists(photolysis_config) && "ERROR: Could not find photolysis config explicitly!");

    // The gaschem mechanism under test is supplied by the caller: either a v0
    // CAMP directory (config_dir) or a v1 single-file mechanism (config_file).
    std::string temp_main_coupled_config = "test_main_coupled_config.yml";
    std::ofstream main_conf_writer(temp_main_coupled_config);
    main_conf_writer << "process:\n";
    main_conf_writer << "  photolysis:\n";
    main_conf_writer << "    active: true\n";
    main_conf_writer << "    config_file: \"" << photolysis_config << "\"\n";
    main_conf_writer << "  gaschem:\n";
    main_conf_writer << "    active: true\n";
    main_conf_writer << "    " << config_key << ": \"" << config_value << "\"\n";
    main_conf_writer.close();

    state->set_configuration_path(temp_main_coupled_config);
    if (state->config_manager()) {
        state->config_manager()->load_from_file(temp_main_coupled_config);
    }

    // 4. Create and add processes to core pipeline
    auto photolysis_proc = catchem::ProcessRegistry::get_instance().create("photolysis");
    photolysis_proc->init(state);
    core->add_process(photolysis_proc);

    auto gaschem_proc = catchem::ProcessRegistry::get_instance().create("gaschem");
    gaschem_proc->init(state);
    core->add_process(gaschem_proc);

    // 5. Execute timestep
    core->run_timestep(3600.0);
    std::cout << "SUCCESS: Timestep run completed without error." << std::endl;

    // 6. Verify photolysis rate diagnostics calculation
    auto diag_mgr = core->get_diagnostic_manager();
    assert(diag_mgr->has_field("photolysis_rate_jfoo") && "TUV-x jfoo diagnostic must exist!");

    double* jfoo_rates = static_cast<double*>(diag_mgr->get_host_pointer("photolysis_rate_jfoo"));
    assert(jfoo_rates != nullptr);
    assert(jfoo_rates[0] > 0.0 && "Photolysis rates must be non-zero at noon!");

    std::cout << "SUCCESS: Dynamic J-rates calculated and coupled successfully!" << std::endl;

    std::remove(temp_main_coupled_config.c_str());
    std::cout << "CASE PASSED: " << label << std::endl;
}

// Runs a multi-column (n_cols = 6, n_levels = 3, n_species = 4) coupled
// photolysis + gaschem scenario for a given optional columns_per_batch setting,
// asserting that GasChemProcess sizes its MICM state to
// expected_batch_cols * n_levels and returns the final concentration buffer.
static std::vector<double> run_multi_column_batch_case(const std::string& label, const std::string& v1_file,
                                                       std::optional<int> configured_columns_per_batch,
                                                       std::size_t expected_batch_cols) {
    std::cout << "\n------------------------------------------" << std::endl;
    std::cout << "BATCH CASE: " << label << " (expected columns_per_batch=" << expected_batch_cols << ")" << std::endl;
    std::cout << "------------------------------------------\n" << std::endl;

    constexpr int n_cols = 6;
    constexpr int n_levels = 3;
    constexpr int n_species = 4; // matches tests/fixtures/mechanisms/chapman.yml

    auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
    auto state = core->get_state_manager();

    state->clock().year = 2026;
    state->clock().month = 7;
    state->clock().day = 13;
    state->clock().hour = 12;
    state->clock().minute = 0;
    state->clock().second = 0;
    state->clock().calculate_derived_fields();

    // Columns 1 and 4 are given identical physical profiles so they land in
    // different batches (including the partial final batch when
    // columns_per_batch = 4) while the remaining columns vary spatially.
    auto profile_slot = [](int icol) { return (icol == 4) ? 1 : icol; };

    std::vector<double> lat(n_cols, 40.0);
    std::vector<double> lon(n_cols, -105.0);
    std::vector<double> temperature(n_cols * n_levels, 280.0);
    std::vector<double> airden(n_cols * n_levels, 1.2);
    std::vector<double> pmid(n_cols * n_levels, 101300.0);
    std::vector<double> pedge(n_cols * (n_levels + 1), 101300.0);
    std::vector<double> bxheight(n_cols * n_levels, 1000.0);

    for (int icol = 0; icol < n_cols; ++icol) {
        const int slot = profile_slot(icol);
        lat[icol] = 35.0 + 2.0 * slot;
        lon[icol] = -105.0 + 1.5 * slot;
        for (int ilev = 0; ilev < n_levels; ++ilev) {
            const std::size_t idx = static_cast<std::size_t>(icol + ilev * n_cols);
            temperature[idx] = 280.0 + 1.5 * slot - 0.5 * ilev;
            airden[idx] = (1.2 - 0.02 * slot) * std::exp(-ilev / 10.0);
            pmid[idx] = (101300.0 - 500.0 * slot) * std::exp(-ilev / 10.0);
            pedge[idx] = pmid[idx];
            bxheight[idx] = 1000.0;
        }
        pedge[static_cast<std::size_t>(icol + n_levels * n_cols)] =
            (101300.0 - 500.0 * slot) * std::exp(-n_levels / 10.0);
    }

    state->bind_met_field_2d("LAT", lat.data());
    state->bind_met_field_2d("LON", lon.data());
    state->bind_met_field_3d("T", temperature.data());
    state->bind_met_field_3d("AIRDEN", airden.data());
    state->bind_met_field_3d("AIRDEN_DRY", airden.data());
    state->bind_met_field_3d("PEDGE", pedge.data());
    state->bind_met_field_3d("PMID", pmid.data());
    state->bind_met_field_3d("BXHEIGHT", bxheight.data());

    std::string species_config = std::string(catchem::test::TEST_DIR) + "/fixtures/mechanisms/chapman.yml";
    assert(file_exists(species_config) && "ERROR: Could not find the chapman mechanism fixture!");
    state->load_species_config(species_config);

    std::vector<double> initial_conc(n_cols * n_levels * n_species, 0.0);
    for (int ilev = 0; ilev < n_levels; ++ilev) {
        for (int icol = 0; icol < n_cols; ++icol) {
            const int slot = profile_slot(icol);
            const auto cell_idx = [&](int ispec) {
                return static_cast<std::size_t>(icol + ilev * n_cols + ispec * n_cols * n_levels);
            };
            // chapman.yml species order: 0=o2, 1=o, 2=o1d, 3=o3.
            // Keep O(1D) at 0 initially (its quenching lifetime with M is ~0.5 ns)
            // so O, O2, and O3 evolve smoothly on a multi-second timescale.
            initial_conc[cell_idx(0)] = 100.0 + 5.0 * slot + 2.0 * ilev; // o2 [ppmv]
            initial_conc[cell_idx(1)] = 1.0 + 0.1 * slot + 0.05 * ilev;  // o  [ppmv]
            initial_conc[cell_idx(2)] = 0.0;                             // o1d [ppmv]
            initial_conc[cell_idx(3)] = 10.0 + 0.5 * slot + 0.2 * ilev;  // o3 [ppmv]
        }
    }
    std::vector<double> conc_data = initial_conc;
    state->bind_unified_chemistry(conc_data.data());

    std::string photolysis_config =
        std::string(catchem::test::SOURCE_DIR) + "/src/external/musica/configs/tuvx/from_host/config.json";
    assert(file_exists(photolysis_config) && "ERROR: Could not find photolysis config explicitly!");

    std::string temp_config = "test_gaschem_batch_config.yml";
    std::ofstream conf_writer(temp_config);
    conf_writer << "processes:\n";
    conf_writer << "  photolysis:\n";
    conf_writer << "    activate: true\n";
    conf_writer << "    config_file: \"" << photolysis_config << "\"\n";
    conf_writer << "  gaschem:\n";
    conf_writer << "    activate: true\n";
    conf_writer << "    config_dir: \"" << v1_file << "\"\n";
    if (configured_columns_per_batch.has_value()) {
        conf_writer << "    columns_per_batch: " << *configured_columns_per_batch << "\n";
    }
    conf_writer.close();

    state->set_configuration_path(temp_config);
    if (state->config_manager()) {
        state->config_manager()->load_from_file(temp_config);
    }

    auto photolysis_proc = catchem::ProcessRegistry::get_instance().create("photolysis");
    photolysis_proc->init(state);
    core->add_process(photolysis_proc);

    auto gaschem_proc = catchem::ProcessRegistry::get_instance().create("gaschem");
    gaschem_proc->init(state);
    core->add_process(gaschem_proc);

    auto* gaschem_typed = dynamic_cast<catchem::GasChemProcess*>(gaschem_proc.get());
    assert(gaschem_typed != nullptr && "gaschem process must be a GasChemProcess instance");
    assert(gaschem_typed->columns_per_batch() == expected_batch_cols &&
           "GasChemProcess::columns_per_batch() must match configured batch size");
    assert(gaschem_typed->micm_grid_cell_count() == expected_batch_cols * static_cast<std::size_t>(n_levels) &&
           "MICM state cell count must equal columns_per_batch * n_levels (not n_cols * n_levels)");

    core->run_timestep(3600.0);
    state->sync_to_host();

    auto conc_view = state->chemistry().conc->host_view;
    bool any_species_changed = false;
    for (int icol = 0; icol < n_cols; ++icol) {
        for (int ilev = 0; ilev < n_levels; ++ilev) {
            for (int ispec = 0; ispec < n_species; ++ispec) {
                const double val = conc_view(icol, ilev, ispec);
                assert(std::isfinite(val) && val >= 0.0 && "Solved concentration must be finite and non-negative");
                const std::size_t idx = static_cast<std::size_t>(icol + ilev * n_cols + ispec * n_cols * n_levels);
                if (std::abs(val - initial_conc[idx]) > 1.0e-12) {
                    any_species_changed = true;
                }
            }
        }
    }
    assert(any_species_changed && "GasChem solve must evolve concentrations over a 3600s step");

    // Distinct columns (0 vs 1) must not have identical output concentrations.
    bool col0_differs_from_col1 = false;
    for (int ilev = 0; ilev < n_levels; ++ilev) {
        for (int ispec = 0; ispec < n_species; ++ispec) {
            if (std::abs(conc_view(0, ilev, ispec) - conc_view(1, ilev, ispec)) > 1.0e-6) {
                col0_differs_from_col1 = true;
            }
        }
    }
    assert(col0_differs_from_col1 && "Columns with distinct inputs must produce distinct outputs");

    std::remove(temp_config.c_str());
    std::cout << "BATCH CASE PASSED: " << label << std::endl;
    return conc_data;
}

static void test_invalid_batch_size(const std::string& v1_file) {
    auto core = std::make_shared<catchem::Core>(2, 3, 4);
    auto state = core->get_state_manager();

    std::string species_config = std::string(catchem::test::TEST_DIR) + "/fixtures/mechanisms/chapman.yml";
    state->load_species_config(species_config);

    std::string temp_config = "test_gaschem_invalid_batch.yml";
    std::ofstream conf_writer(temp_config);
    conf_writer << "processes:\n";
    conf_writer << "  gaschem:\n";
    conf_writer << "    activate: true\n";
    conf_writer << "    config_dir: \"" << v1_file << "\"\n";
    conf_writer << "    columns_per_batch: 0\n";
    conf_writer.close();

    state->set_configuration_path(temp_config);
    if (state->config_manager()) {
        state->config_manager()->load_from_file(temp_config);
    }

    auto gaschem_proc = catchem::ProcessRegistry::get_instance().create("gaschem");
    bool threw = false;
    try {
        gaschem_proc->init(state);
    } catch (const std::exception&) {
        threw = true;
    }
    std::remove(temp_config.c_str());
    assert(threw && "columns_per_batch < 1 must fail initialization");
    std::cout << "SUCCESS: Invalid columns_per_batch=0 rejected at init." << std::endl;
}

static void run_multi_column_batch_suite(const std::string& v1_file) {
    constexpr int n_cols = 6;
    constexpr int n_levels = 3;
    constexpr int n_species = 4;

    // 1. Default (omitted columns_per_batch -> 1)
    const auto conc_default = run_multi_column_batch_case("default (omitted -> 1)", v1_file, std::nullopt, 1);

    // 2. Explicit batch size 1 (6 batches of 1 column)
    const auto conc_b1 = run_multi_column_batch_case("columns_per_batch=1", v1_file, 1, 1);
    assert(conc_default == conc_b1 && "Omitted columns_per_batch must be bit-for-bit identical to columns_per_batch=1");

    // Verify identical-input columns (1 and 4) match bit-for-bit when solved in single-column batches.
    for (int ilev = 0; ilev < n_levels; ++ilev) {
        for (int ispec = 0; ispec < n_species; ++ispec) {
            const std::size_t idx1 = static_cast<std::size_t>(1 + ilev * n_cols + ispec * n_cols * n_levels);
            const std::size_t idx4 = static_cast<std::size_t>(4 + ilev * n_cols + ispec * n_cols * n_levels);
            assert(conc_b1[idx1] == conc_b1[idx4] &&
                   "Columns 1 and 4 with identical inputs must match bit-for-bit at columns_per_batch=1");
        }
    }

    // 3. Batch size 2 (even divisor: 3 batches of 2 columns)
    const auto conc_b2 = run_multi_column_batch_case("columns_per_batch=2 (even divisor)", v1_file, 2, 2);

    // 4. Batch size 4 (non-divisor partial final batch: batch 0 has 4 cols, batch 1 has 2 valid + 2 padded cols)
    const auto conc_b4 = run_multi_column_batch_case("columns_per_batch=4 (partial final batch)", v1_file, 4, 4);

    // 5. Batch size 6 (whole domain in 1 batch)
    const auto conc_b6 = run_multi_column_batch_case("columns_per_batch=6 (full domain)", v1_file, 6, 6);

    // Compare all multi-column batch sizes against the columns_per_batch=1 reference within Rosenbrock tolerance
    // (atol + rtol * |ref|, where v0/chapman species have atol = 1e-12 mol/m3 ~ 2.5e-8 ppmv and rtol = 1e-4).
    for (const auto& [batch_label, candidate] : std::vector<std::pair<std::string, std::vector<double>>>{
             {"columns_per_batch=2", conc_b2}, {"columns_per_batch=4", conc_b4}, {"columns_per_batch=6", conc_b6}}) {
        double max_abs_diff = 0.0;
        double max_rel_diff = 0.0;
        for (std::size_t i = 0; i < conc_b1.size(); ++i) {
            const double ref = conc_b1[i];
            const double val = candidate[i];
            const double abs_diff = std::abs(val - ref);
            const double tol = 1.0e-6 + 1.0e-2 * std::abs(ref);
            max_abs_diff = std::max(max_abs_diff, abs_diff);
            if (std::abs(ref) > 1.0e-6) {
                max_rel_diff = std::max(max_rel_diff, abs_diff / std::abs(ref));
            }
            if (abs_diff > tol) {
                std::cerr << "Tolerance exceeded in " << batch_label << " at index " << i << ": ref=" << ref
                          << " val=" << val << " abs_diff=" << abs_diff << " tol=" << tol << std::endl;
            }
            assert(abs_diff <= tol &&
                   "Batched solve concentrations must agree across batch sizes within solver tolerance");
        }
        std::cout << "SUCCESS: " << batch_label << " agrees with columns_per_batch=1 (max_abs_diff=" << max_abs_diff
                  << ", max_rel_diff=" << max_rel_diff << ")." << std::endl;
    }

    // 6. Reject non-positive columns_per_batch.
    test_invalid_batch_size(v1_file);
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "RUNNING INTEGRATION TEST: Photolysis and GasChem Coupling" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        // 1. Dynamic registrations
        catchem_register_photolysis_cpp();
        catchem_register_gaschem_cpp();
        assert(catchem::ProcessRegistry::get_instance().has_process("photolysis"));
        assert(catchem::ProcessRegistry::get_instance().has_process("gaschem"));

        // Case A: legacy v0 CAMP mechanism directory (multi-file layout).
        std::string v0_dir = std::string(catchem::test::SOURCE_DIR) + "/src/external/musica/configs/v0/chapman/";
        assert((file_exists(v0_dir + "config.json") || file_exists(v0_dir + "config.yaml")) &&
               "ERROR: Could not find gaschem v0 config dir explicitly!");
        run_coupled_case("v0 CAMP directory", "config_dir", v0_dir);

        // Case B: MUSICA v1 single-file mechanism passed via the new config_file key.
        std::string v1_file =
            std::string(catchem::test::SOURCE_DIR) + "/src/external/musica/configs/v1/chapman/config.yaml";
        assert(file_exists(v1_file) && "ERROR: Could not find gaschem v1 config file explicitly!");
        run_coupled_case("v1 single file (config_file)", "config_file", v1_file);

        // Case C: v1 mechanism directory exercising config.yaml auto-detection.
        std::string v1_dir = std::string(catchem::test::SOURCE_DIR) + "/src/external/musica/configs/v1/chapman/";
        assert(file_exists(v1_dir + "config.yaml") && "ERROR: Could not find gaschem v1 config dir explicitly!");
        run_coupled_case("v1 directory (config.yaml auto-detect)", "config_dir", v1_dir);

        // Case D: Multi-column batching suite across batch sizes (1, 2, 4, 6) and invalid batch rejection.
        run_multi_column_batch_suite(v0_dir);

        std::cout << "\n==========================================" << std::endl;
        std::cout << "COUPLED INTEGRATION TESTS PASSED (v0 dir, v1 file, v1 auto-detect, batching)!" << std::endl;
        std::cout << "==========================================\n" << std::endl;
    }
    Kokkos::finalize();
    return 0;
}
