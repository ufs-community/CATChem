#include "catchem_api.hpp"
#include "catchem_config_manager.hpp"
#include "catchem_constants.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_logger.hpp"
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
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <yaml-cpp/yaml.h>

extern "C" {
void catchem_register_photolysis_cpp();
void catchem_register_gaschem_cpp();
}

bool file_exists(const std::string& name) {
    return std::filesystem::exists(name);
}

// Runs the full coupled photolysis+gaschem scenario once for a given gaschem
// mechanism configuration (a v0 CAMP directory, a v1 single file, or a v1
// directory relying on config.yaml auto-detection). The species fixture must
// list every MICM variable of that mechanism.
static void run_coupled_case(const std::string& label, const std::string& config_key, const std::string& config_value,
                             const std::string& species_fixture = "chapman.yml", int n_species = 4) {
    std::cout << "\n------------------------------------------" << std::endl;
    std::cout << "CASE: gaschem mechanism via " << config_key << " (" << label << ")" << std::endl;
    std::cout << "------------------------------------------\n" << std::endl;

    // 2. Set up core and states
    int n_cols = 1;
    int n_levels = 3;

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

    // Load a species list that covers the mechanism under test and carries the
    // photolysis.ozone role the photolysis contract requires.
    std::string species_config = std::string(catchem::test::TEST_DIR) + "/fixtures/mechanisms/" + species_fixture;
    assert(file_exists(species_config) && "ERROR: Could not find the chapman mechanism fixture!");
    state->load_species_config(species_config);
    std::vector<double> conc_data(n_cols * n_levels * n_species, 1.0); // 1.0 ppmv initially
    // Molecular nitrogen, where the mechanism solves for it, starts at its
    // atmospheric mixing ratio rather than 1 ppmv.
    const auto n2 = state->chemistry().species_name_to_index.find("N2");
    if (n2 != state->chemistry().species_name_to_index.end()) {
        std::fill_n(conc_data.begin() + n2->second * n_cols * n_levels, n_cols * n_levels, 7.8084e5);
    }
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
    main_conf_writer << "    activate: true\n";
    main_conf_writer << "    config_file: \"" << photolysis_config << "\"\n";
    main_conf_writer << "  gaschem:\n";
    main_conf_writer << "    activate: true\n";
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

// One cell of a (column, level[, species]) field and the value to put there.
struct CellValue {
    int icol;
    int ilev;
    int ispec; // ignored for meteorology
    double value;
};

// What a scripted case can reach after its first timestep.
struct BatchCaseRun {
    catchem::Core& core;
    catchem::StateManager& state;
    catchem::GasChemProcess& gaschem;
};

// Departures from the baseline multi-column case, used by the guard and
// negative_concentration_policy tests.
struct BatchCaseOptions {
    // Species file under tests/fixtures/mechanisms. Alternatives must keep the
    // four chapman species in the same order.
    std::string species_fixture = "chapman.yml";
    // Written as processes/gaschem/negative_concentration_policy when not empty.
    std::string negative_concentration_policy;
    // Top-level key of the process blocks: "processes", or the legacy "process".
    std::string process_block = "processes";
    std::optional<CellValue> conc_override;
    std::optional<CellValue> temperature_override;
    // Photolysis runs before gaschem and would see a corrupted temperature first.
    bool with_photolysis = true;
    // The timestep must throw, and the error must contain every fragment.
    bool expect_failure = false;
    std::vector<std::string> expected_message_parts;
    // Builds the gaschem process; the registry builds it when this is empty. A
    // case uses it to run a GasChemProcess subclass that overrides the
    // solver-output test seam.
    std::function<std::shared_ptr<catchem::GasChemProcess>()> make_gaschem;
    // Called after the first timestep and before the final checks, to inspect
    // the process, change inputs, and run further timesteps.
    std::function<void(BatchCaseRun&)> after_first_step;
};

// What GasChemProcess reported after a case: the cells it floored to zero and
// the relative tolerance its solver ran with.
struct CaseStats {
    std::size_t inputs = 0;
    std::size_t outputs = 0;
    double min_input = 0.0;
    double min_output = 0.0;
    double solver_relative_tolerance = 0.0;
};

// Dry air density [kg/m3] of the multi-column case for a profile slot and level.
static double batch_case_air_density(int slot, int ilev) {
    return (1.2 - 0.02 * slot) * std::exp(-ilev / 10.0);
}

// The largest "absolute tolerance" [mol/m3] a v0 CAMP mechanism declares for
// its species, read from species.json in the mechanism directory.
static double largest_declared_absolute_tolerance(const std::string& mechanism_dir) {
    const YAML::Node species =
        YAML::LoadFile((std::filesystem::path(mechanism_dir) / "species.json").string())["camp-data"];
    assert(species && species.IsSequence() && "species.json must hold a camp-data list");
    double largest = 0.0;
    for (const auto& entry : species) {
        if (entry["absolute tolerance"]) {
            largest = std::max(largest, entry["absolute tolerance"].as<double>());
        }
    }
    assert(largest > 0.0 && "The mechanism must declare an absolute tolerance for a species");
    return largest;
}

// Runs a multi-column (n_cols = 6, n_levels = 3, n_species = 4) coupled
// photolysis + gaschem scenario for a given optional columns_per_batch setting,
// asserting that GasChemProcess sizes its MICM state to
// expected_batch_cols * n_levels and returns the final concentration buffer.
// When options.expect_failure is set, the timestep must throw instead and an
// empty buffer is returned.
static std::vector<double> run_multi_column_batch_case(const std::string& label, const std::string& mechanism_dir,
                                                       std::optional<int> configured_columns_per_batch,
                                                       std::size_t expected_batch_cols,
                                                       const BatchCaseOptions& options = {},
                                                       CaseStats* stats = nullptr) {
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
            airden[idx] = batch_case_air_density(slot, ilev);
            pmid[idx] = (101300.0 - 500.0 * slot) * std::exp(-ilev / 10.0);
            pedge[idx] = pmid[idx];
            bxheight[idx] = 1000.0;
        }
        pedge[static_cast<std::size_t>(icol + n_levels * n_cols)] =
            (101300.0 - 500.0 * slot) * std::exp(-n_levels / 10.0);
    }
    if (options.temperature_override.has_value()) {
        const auto& cell = *options.temperature_override;
        temperature[static_cast<std::size_t>(cell.icol + cell.ilev * n_cols)] = cell.value;
    }

    state->bind_met_field_2d("LAT", lat.data());
    state->bind_met_field_2d("LON", lon.data());
    state->bind_met_field_3d("T", temperature.data());
    state->bind_met_field_3d("AIRDEN", airden.data());
    state->bind_met_field_3d("AIRDEN_DRY", airden.data());
    state->bind_met_field_3d("PEDGE", pedge.data());
    state->bind_met_field_3d("PMID", pmid.data());
    state->bind_met_field_3d("BXHEIGHT", bxheight.data());

    std::string species_config =
        std::string(catchem::test::TEST_DIR) + "/fixtures/mechanisms/" + options.species_fixture;
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
    if (options.conc_override.has_value()) {
        const auto& cell = *options.conc_override;
        initial_conc[static_cast<std::size_t>(cell.icol + cell.ilev * n_cols + cell.ispec * n_cols * n_levels)] =
            cell.value;
    }
    std::vector<double> conc_data = initial_conc;
    state->bind_unified_chemistry(conc_data.data());

    std::string photolysis_config =
        std::string(catchem::test::SOURCE_DIR) + "/src/external/musica/configs/tuvx/from_host/config.json";
    assert(file_exists(photolysis_config) && "ERROR: Could not find photolysis config explicitly!");

    std::string temp_config = "test_gaschem_batch_config.yml";
    std::ofstream conf_writer(temp_config);
    conf_writer << options.process_block << ":\n";
    if (options.with_photolysis) {
        conf_writer << "  photolysis:\n";
        conf_writer << "    activate: true\n";
        conf_writer << "    config_file: \"" << photolysis_config << "\"\n";
    }
    conf_writer << "  gaschem:\n";
    conf_writer << "    activate: true\n";
    conf_writer << "    config_dir: \"" << mechanism_dir << "\"\n";
    if (configured_columns_per_batch.has_value()) {
        conf_writer << "    columns_per_batch: " << *configured_columns_per_batch << "\n";
    }
    if (!options.negative_concentration_policy.empty()) {
        conf_writer << "    negative_concentration_policy: " << options.negative_concentration_policy << "\n";
    }
    conf_writer.close();

    state->set_configuration_path(temp_config);
    if (state->config_manager()) {
        state->config_manager()->load_from_file(temp_config);
    }

    if (options.with_photolysis) {
        auto photolysis_proc = catchem::ProcessRegistry::get_instance().create("photolysis");
        photolysis_proc->init(state);
        core->add_process(photolysis_proc);
    }

    std::shared_ptr<catchem::ProcessInterface> gaschem_proc;
    if (options.make_gaschem) {
        gaschem_proc = options.make_gaschem();
    } else {
        gaschem_proc = catchem::ProcessRegistry::get_instance().create("gaschem");
    }
    gaschem_proc->init(state);
    core->add_process(gaschem_proc);

    auto* gaschem_typed = dynamic_cast<catchem::GasChemProcess*>(gaschem_proc.get());
    assert(gaschem_typed != nullptr && "gaschem process must be a GasChemProcess instance");
    assert(gaschem_typed->columns_per_batch() == expected_batch_cols &&
           "GasChemProcess::columns_per_batch() must match configured batch size");
    assert(gaschem_typed->micm_grid_cell_count() == expected_batch_cols * static_cast<std::size_t>(n_levels) &&
           "MICM state cell count must equal columns_per_batch * n_levels (not n_cols * n_levels)");

    if (options.expect_failure) {
        std::string message;
        try {
            core->run_timestep(3600.0);
        } catch (const std::domain_error& error) {
            message = error.what();
        }
        std::remove(temp_config.c_str());
        assert(!message.empty() && "A rejected gaschem concentration must raise std::domain_error");
        for (const auto& part : options.expected_message_parts) {
            const bool found = message.find(part) != std::string::npos;
            if (!found) {
                std::cerr << "Guard message is missing \"" << part << "\": " << message << std::endl;
            }
            assert(found && "GasChem guard message must contain every expected fragment");
        }
        std::cout << "GUARD CASE PASSED: " << label << " -> " << message << std::endl;
        return {};
    }

    core->run_timestep(3600.0);
    if (options.after_first_step) {
        BatchCaseRun run{*core, *state, *gaschem_typed};
        options.after_first_step(run);
    }
    state->sync_to_host();
    if (stats != nullptr) {
        stats->inputs = gaschem_typed->negative_inputs_clamped();
        stats->outputs = gaschem_typed->negative_outputs_clamped();
        stats->min_input = gaschem_typed->min_negative_input();
        stats->min_output = gaschem_typed->min_negative_output();
        stats->solver_relative_tolerance = gaschem_typed->solver_relative_tolerance();
    }

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

// GasChemProcess::init must reject the given processes/gaschem option line.
static void test_invalid_gaschem_option(const std::string& mechanism_dir, const std::string& option_line) {
    auto core = std::make_shared<catchem::Core>(2, 3, 4);
    auto state = core->get_state_manager();

    std::string species_config = std::string(catchem::test::TEST_DIR) + "/fixtures/mechanisms/chapman.yml";
    state->load_species_config(species_config);

    std::string temp_config = "test_gaschem_invalid_option.yml";
    std::ofstream conf_writer(temp_config);
    conf_writer << "processes:\n";
    conf_writer << "  gaschem:\n";
    conf_writer << "    activate: true\n";
    conf_writer << "    config_dir: \"" << mechanism_dir << "\"\n";
    conf_writer << "    " << option_line << "\n";
    conf_writer.close();

    state->set_configuration_path(temp_config);
    if (state->config_manager()) {
        state->config_manager()->load_from_file(temp_config);
    }

    auto gaschem_proc = catchem::ProcessRegistry::get_instance().create("gaschem");
    bool threw = false;
    try {
        gaschem_proc->init(state);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    std::remove(temp_config.c_str());
    if (!threw) {
        std::cerr << "Not rejected at init: " << option_line << std::endl;
    }
    assert(threw && "An invalid gaschem option must fail initialization with std::invalid_argument");
    std::cout << "SUCCESS: Invalid option rejected at init: " << option_line << std::endl;
}

// Initializes gaschem alone for one mechanism and species fixture and returns
// the std::invalid_argument message it raised, or an empty string if
// initialization succeeded.
static std::string gaschem_init_error(const std::string& config_key, const std::string& config_value,
                                      const std::string& species_fixture, int n_species) {
    auto core = std::make_shared<catchem::Core>(2, 3, n_species);
    auto state = core->get_state_manager();
    state->load_species_config(std::string(catchem::test::TEST_DIR) + "/fixtures/mechanisms/" + species_fixture);
    const std::string temp_config = "test_gaschem_init_error.yml";
    std::ofstream conf_writer(temp_config);
    conf_writer << "processes:\n";
    conf_writer << "  gaschem:\n";
    conf_writer << "    activate: true\n";
    conf_writer << "    " << config_key << ": \"" << config_value << "\"\n";
    conf_writer.close();
    state->set_configuration_path(temp_config);
    if (state->config_manager()) {
        state->config_manager()->load_from_file(temp_config);
    }
    auto gaschem_proc = catchem::ProcessRegistry::get_instance().create("gaschem");
    std::string message;
    try {
        gaschem_proc->init(state);
    } catch (const std::invalid_argument& error) {
        message = error.what();
    }
    std::remove(temp_config.c_str());
    return message;
}

// Asserts that a configuration error was raised and contains every fragment.
// stage says where the error was expected, for the log.
static void expect_init_error(const std::string& what, const std::string& message,
                              const std::vector<std::string>& fragments, const std::string& stage = "at init") {
    if (message.empty()) {
        std::cerr << "Not rejected " << stage << ": " << what << std::endl;
    }
    assert(!message.empty() && "the gaschem configuration must be refused");
    for (const auto& fragment : fragments) {
        const bool found = message.find(fragment) != std::string::npos;
        if (!found) {
            std::cerr << what << ": message is missing \"" << fragment << "\": " << message << std::endl;
        }
        assert(found && "gaschem configuration error must contain every expected fragment");
    }
    std::cout << "SUCCESS: " << what << " rejected " << stage << " -> " << message << std::endl;
}

// MICM variables are resolved to CATChem species through the mechanism
// definition: case-insensitively, through aliases, at most one MICM variable
// per CATChem species, and with no MICM variable left without one.
static void run_species_resolution_suite(const std::string& mechanism_dir, const std::string& v1_mechanism_file) {
    // 1. A species list in mixed case, with one species reachable only through an
    // alias, maps the same four species and so solves bit-for-bit like chapman.yml.
    const auto conc_reference = run_multi_column_batch_case("species resolution reference", mechanism_dir, 4, 4);
    BatchCaseOptions renamed;
    renamed.species_fixture = "chapman_alias_mixed_case.yml";
    const auto conc_renamed =
        run_multi_column_batch_case("mixed-case names and an alias", mechanism_dir, 4, 4, renamed);
    assert(conc_renamed == conc_reference &&
           "Species resolved by case-insensitive name or alias must solve bit-for-bit like exact names");

    // 2. Two MICM variables resolving to one CATChem species is a configuration error.
    expect_init_error("Duplicate species mapping",
                      gaschem_init_error("config_dir", mechanism_dir, "chapman_duplicate_alias.yml", 3),
                      {"O and O1D", "both resolve to CATChem species o"});

    // 3. Every MICM variable must be in the species list. MICM solves for all of
    // them, and with a state reused across column batches an unsupplied one would
    // start from another column's result. The v0 mechanism solves for O.
    expect_init_error("Species list without atomic oxygen (v0 mechanism)",
                      gaschem_init_error("config_dir", mechanism_dir, "chapman_without_o.yml", 3),
                      {"missing from the CATChem species list: O.", "is_advected: false"});

    // 4. The v1 mechanism also solves for N2, which chapman.yml does not list.
    expect_init_error("Species list without N2 (v1 mechanism)",
                      gaschem_init_error("config_file", v1_mechanism_file, "chapman.yml", 4),
                      {"missing from the CATChem species list: N2.", "is_advected: false"});
}

// processes/gaschem/negative_concentration_policy: reject (default) throws on a
// negative concentration, clamp floors it to zero. Neither relaxes the checks on
// non-finite concentrations or on meteorology. Atomic oxygen (species 1) carries
// the bad value because photolysis runs first and reads ozone.
static void run_negative_concentration_policy_suite(const std::string& mechanism_dir) {
    constexpr int o_atom = 1;
    constexpr double tiny_negative = -1.0e-12;
    const std::string input_guard = "GasChem received a negative or non-finite concentration";

    // 1. Key omitted and explicit reject: the message shows the magnitude and
    // where it is. Column 4 opens the short final batch when columns_per_batch = 4,
    // so the reported column must be the PET-local index, not the index in the batch.
    for (const char* policy : {"", "reject"}) {
        BatchCaseOptions rejected;
        rejected.negative_concentration_policy = policy;
        rejected.conc_override = CellValue{4, 2, o_atom, tiny_negative};
        rejected.expect_failure = true;
        rejected.expected_message_parts = {input_guard, "(-1.000000e-12)", " for O ", "column=4", "level=2", "batch=1"};
        run_multi_column_batch_case(std::string("negative input rejected (policy='") + policy + "')", mechanism_dir, 4,
                                    4, rejected);
    }

    // 2. clamp: the negative cell is solved exactly as if it had been zero.
    BatchCaseOptions zero_cell;
    zero_cell.negative_concentration_policy = "clamp";
    zero_cell.conc_override = CellValue{4, 2, o_atom, 0.0};
    CaseStats zero_stats;
    const auto conc_zero =
        run_multi_column_batch_case("clamp, cell at zero", mechanism_dir, 4, 4, zero_cell, &zero_stats);
    assert(zero_stats.inputs == 0 && zero_stats.outputs == 0 && "A zero concentration must not be counted as floored");

    BatchCaseOptions negative_cell = zero_cell;
    negative_cell.conc_override = CellValue{4, 2, o_atom, tiny_negative};
    CaseStats negative_stats;
    const auto conc_negative =
        run_multi_column_batch_case("clamp, cell tiny negative", mechanism_dir, 4, 4, negative_cell, &negative_stats);
    assert(conc_negative == conc_zero && "A floored negative input must solve bit-for-bit like a zero input");
    assert(negative_stats.inputs == 1 && "Exactly one input cell must be counted as floored");
    assert(negative_stats.min_input == tiny_negative && "The most negative floored input must be recorded");
    assert(negative_stats.outputs == 0 && "No solver output was negative in this case");

    // 3. clamp does not accept a non-finite concentration.
    BatchCaseOptions nan_cell = negative_cell;
    nan_cell.conc_override = CellValue{4, 2, o_atom, std::numeric_limits<double>::quiet_NaN()};
    nan_cell.expect_failure = true;
    nan_cell.expected_message_parts = {input_guard, " for O ", "column=4", "level=2"};
    run_multi_column_batch_case("clamp, NaN input rejected", mechanism_dir, 4, 4, nan_cell);

    // 4. clamp counts real cells only. Column 5 is the last valid column of the
    // short final batch, so its values are copied into the two padding columns.
    BatchCaseOptions padded_cell = negative_cell;
    padded_cell.conc_override = CellValue{5, 0, o_atom, tiny_negative};
    CaseStats padded_stats;
    run_multi_column_batch_case("clamp, negative in padded column", mechanism_dir, 4, 4, padded_cell, &padded_stats);
    assert(padded_stats.inputs == 1 && "Padding copies of a floored cell must not be counted");

    // 5. clamp does not relax gaschem's own meteorology check. The state manager
    // is left at its lenient legacy default here, so only gaschem can stop the run.
    BatchCaseOptions bad_temperature;
    bad_temperature.negative_concentration_policy = "clamp";
    bad_temperature.temperature_override = CellValue{1, 1, 0, -1.0};
    bad_temperature.with_photolysis = false;
    bad_temperature.expect_failure = true;
    bad_temperature.expected_message_parts = {"non-physical T, PMID, or AIRDEN_DRY"};
    run_multi_column_batch_case("clamp, negative temperature rejected", mechanism_dir, 4, 4, bad_temperature);

    // 6. Only reject and clamp are accepted.
    test_invalid_gaschem_option(mechanism_dir, "negative_concentration_policy: bogus");
}

// GasChemProcess whose solver "returns" a negative concentration for every
// cell of one batch on one run() call, through the test seam. The Chapman
// mechanism never produces a negative on its own.
class NegativeOutputGasChem : public catchem::GasChemProcess {
public:
    static constexpr double forced_mol_m3 = -1.0e-20;

    // on_run counts run() calls from 1; in_batch is zero-based.
    NegativeOutputGasChem(int on_run, int in_batch) : on_run_(on_run), in_batch_(in_batch) {}

    void run(std::shared_ptr<catchem::StateManager> state) override {
        ++run_count_;
        catchem::GasChemProcess::run(std::move(state));
    }

protected:
    void adjust_solver_output_for_testing(std::vector<double>& concentrations, int batch) override {
        if (run_count_ == on_run_ && batch == in_batch_) {
            std::fill(concentrations.begin(), concentrations.end(), forced_mol_m3);
        }
    }

private:
    int on_run_;
    int in_batch_;
    int run_count_ = 0;
};

// Holds everything written to std::clog, where Logger sends its debug, info and
// warn lines, and passes it on when destroyed so the test log stays complete.
class ClogCapture {
public:
    ClogCapture() : previous_(std::clog.rdbuf(buffer_.rdbuf())) {}
    ~ClogCapture() {
        std::clog.rdbuf(previous_);
        std::clog << buffer_.str() << std::flush;
    }
    ClogCapture(const ClogCapture&) = delete;
    ClogCapture& operator=(const ClogCapture&) = delete;

    // Number of "floored negative concentrations to zero" lines at one level
    // ("WARN" or "DEBUG") that also contain the fragment.
    std::size_t floored_lines(const std::string& level, const std::string& fragment = "") const {
        std::istringstream lines(buffer_.str());
        std::string line;
        std::size_t count = 0;
        while (std::getline(lines, line)) {
            if (line.find("floored negative concentrations to zero") != std::string::npos &&
                line.find(level) != std::string::npos && line.find(fragment) != std::string::npos) {
                ++count;
            }
        }
        return count;
    }

private:
    std::ostringstream buffer_;
    std::streambuf* previous_;
};

// What one run of the three-step clamp sequence left behind.
struct ClampSequenceResult {
    CaseStats stats;
    std::vector<double> conc;
    std::size_t warnings = 0;
    std::size_t input_warnings = 0;
    std::size_t output_warnings = 0;
    std::size_t debug_lines = 0;
};

// Three timesteps under clamp at the given log level: a negative input on step
// 1, a negative from the solver for all of batch 1 on step 2, and another
// negative input on step 3. With columns_per_batch = 4, batch 1 holds columns 4
// and 5 plus two padding columns.
static ClampSequenceResult run_clamp_warning_sequence(const std::string& mechanism_dir, catchem::Logger::Level level,
                                                      const std::string& label) {
    constexpr int n_cols = 6;
    constexpr int n_levels = 3;
    constexpr int n_species = 4;
    constexpr int o2 = 0;
    constexpr int o_atom = 1;
    constexpr double tiny_negative = -1.0e-12;
    // The real cells of batch 1: columns 4 and 5, every level and species.
    constexpr std::size_t cells_in_final_batch = 2 * n_levels * n_species;

    BatchCaseOptions sequence;
    sequence.negative_concentration_policy = "clamp";
    sequence.conc_override = CellValue{4, 2, o_atom, tiny_negative};
    sequence.make_gaschem = [] { return std::make_shared<NegativeOutputGasChem>(2, 1); };
    sequence.after_first_step = [&](BatchCaseRun& run) {
        // Step 1 floored one input and nothing from the solver.
        assert(run.gaschem.input_clamp_warned() && "The first floored input must be warned about");
        assert(!run.gaschem.output_clamp_warned() && "A floored input must not count as a solver-output warning");
        assert(run.gaschem.negative_inputs_clamped() == 1 && run.gaschem.negative_outputs_clamped() == 0);

        // Step 2: the solver returns negatives. They are floored, counted once
        // per real cell, and written back as zero.
        run.core.run_timestep(3600.0);
        run.state.sync_to_host();
        assert(run.gaschem.output_clamp_warned() && "The first floored solver output must be warned about");
        assert(run.gaschem.negative_inputs_clamped() == 1 && "Step 2 had no negative input");
        if (run.gaschem.negative_outputs_clamped() != cells_in_final_batch) {
            std::cerr << "outputs floored: " << run.gaschem.negative_outputs_clamped() << std::endl;
        }
        assert(run.gaschem.negative_outputs_clamped() == cells_in_final_batch &&
               "Every real cell of the batch must be counted once, and no padding cell");
        assert(run.gaschem.min_negative_output() < 0.0 && run.gaschem.min_negative_output() > -1.0e-12 &&
               "The most negative floored solver output must be recorded in ppmv");
        auto conc = run.state.chemistry().conc->host_view;
        for (int icol = 0; icol < n_cols; ++icol) {
            for (int ilev = 0; ilev < n_levels; ++ilev) {
                if (icol < 4) {
                    assert(conc(icol, ilev, o2) > 0.0 && "Batch 0 must keep its own solver results");
                    continue;
                }
                for (int ispec = 0; ispec < n_species; ++ispec) {
                    assert(conc(icol, ilev, ispec) == 0.0 && "A floored solver output must be written as zero");
                }
            }
        }

        // Step 3: a floored input in batch 0. Both warnings are already out.
        run.state.chemistry().conc->host_write();
        conc(0, 1, o_atom) = tiny_negative;
        run.core.run_timestep(3600.0);
        assert(run.gaschem.input_clamp_warned() && run.gaschem.output_clamp_warned());
        assert(run.gaschem.negative_inputs_clamped() == 2 && "The step 3 input must be counted");
        assert(run.gaschem.negative_outputs_clamped() == cells_in_final_batch && "Step 3 had no negative output");
    };

    ClampSequenceResult result;
    catchem::Logger::set_level(level);
    {
        ClogCapture log;
        result.conc = run_multi_column_batch_case(label, mechanism_dir, 4, 4, sequence, &result.stats);
        result.warnings = log.floored_lines("WARN");
        result.input_warnings = log.floored_lines("WARN", "first_occurrence=inputs inputs_floored=1 ");
        result.output_warnings = log.floored_lines("WARN", "first_occurrence=solver_outputs inputs_floored=0 ");
        result.debug_lines = log.floored_lines("DEBUG", "first_occurrence=none ");
    }
    catchem::Logger::clear_level();
    return result;
}

// Under clamp a PET warns once about floored inputs and once about floored
// solver outputs, and what it counts does not depend on the log level. Under
// reject a negative from the solver stops the run.
static void run_clamp_warning_suite(const std::string& mechanism_dir) {
    // 1. An input warning on step 1 does not silence the solver-output warning
    // on step 2, and step 3 logs at debug only.
    const auto with_debug =
        run_clamp_warning_sequence(mechanism_dir, catchem::Logger::Level::Debug, "clamp warnings, debug logging on");
    if (with_debug.warnings != 2 || with_debug.input_warnings != 1 || with_debug.output_warnings != 1 ||
        with_debug.debug_lines != 1) {
        std::cerr << "warnings=" << with_debug.warnings << " input_warnings=" << with_debug.input_warnings
                  << " output_warnings=" << with_debug.output_warnings << " debug_lines=" << with_debug.debug_lines
                  << std::endl;
    }
    assert(with_debug.warnings == 2 && "Exactly two warnings: one for inputs, one for solver outputs");
    assert(with_debug.input_warnings == 1 && "Step 1 must warn about the floored input");
    assert(with_debug.output_warnings == 1 && "Step 2 must warn about the floored solver outputs");
    assert(with_debug.debug_lines == 1 && "Step 3 must log at debug, not warn again");

    // 2. With debug logging off the same warnings appear, no debug line is
    // written, and the counters, minima and concentrations are unchanged.
    const auto without_debug =
        run_clamp_warning_sequence(mechanism_dir, catchem::Logger::Level::Warn, "clamp warnings, debug logging off");
    assert(without_debug.warnings == 2 && without_debug.input_warnings == 1 && without_debug.output_warnings == 1);
    assert(without_debug.debug_lines == 0 && "No debug line may be written at warn level");
    assert(without_debug.stats.inputs == with_debug.stats.inputs &&
           without_debug.stats.outputs == with_debug.stats.outputs &&
           without_debug.stats.min_input == with_debug.stats.min_input &&
           without_debug.stats.min_output == with_debug.stats.min_output &&
           "Clamp counters and minima must not depend on the log level");
    assert(without_debug.conc == with_debug.conc && "Concentrations must not depend on the log level");
    assert(with_debug.stats.inputs == 2 && with_debug.stats.min_input == -1.0e-12);

    // 3. Inputs and solver outputs floored for the first time on the same step
    // share one warning that names both.
    BatchCaseOptions both_at_once;
    both_at_once.negative_concentration_policy = "clamp";
    both_at_once.conc_override = CellValue{4, 2, 1, -1.0e-12};
    both_at_once.make_gaschem = [] { return std::make_shared<NegativeOutputGasChem>(1, 1); };
    both_at_once.after_first_step = [](BatchCaseRun& run) {
        assert(run.gaschem.input_clamp_warned() && run.gaschem.output_clamp_warned());
    };
    {
        ClogCapture log;
        run_multi_column_batch_case("clamp, input and solver output floored on one step", mechanism_dir, 4, 4,
                                    both_at_once);
        assert(log.floored_lines("WARN") == 1 &&
               log.floored_lines("WARN", "first_occurrence=inputs,solver_outputs inputs_floored=1 ") == 1 &&
               "One warning must name both inputs and solver outputs");
    }

    // 4. reject: a negative from the solver throws, with the value in ppmv and
    // its location. Batch 1 starts at column 4; O2 is the first mapped species.
    BatchCaseOptions rejected_output;
    rejected_output.make_gaschem = [] { return std::make_shared<NegativeOutputGasChem>(1, 1); };
    rejected_output.expect_failure = true;
    rejected_output.expected_message_parts = {"GasChem solver returned a negative or non-finite concentration",
                                              "(-2.",
                                              "e-16)",
                                              " for O2 ",
                                              "column=4",
                                              "level=0",
                                              "batch=1"};
    run_multi_column_batch_case("negative solver output rejected", mechanism_dir, 4, 4, rejected_output);
}

// Loads a runtime YAML whose gaschem block holds the given option lines and
// returns the parsed gaschem ProcessConfig. block is "processes" or the legacy
// "process".
static catchem::ProcessConfig gaschem_process_config(const std::vector<std::string>& option_lines,
                                                     const std::string& block = "processes") {
    const std::string temp_config = "test_gaschem_settings.yml";
    std::ofstream conf_writer(temp_config);
    conf_writer << block << ":\n";
    conf_writer << "  gaschem:\n";
    conf_writer << "    activate: true\n";
    for (const auto& line : option_lines) {
        conf_writer << "    " << line << "\n";
    }
    conf_writer.close();
    catchem::ConfigManager config;
    config.load_from_file(temp_config);
    std::remove(temp_config.c_str());
    return config.data.processes.at("gaschem");
}

// The settings validator registered for gaschem must reject these option lines.
static void expect_settings_rejected(const std::string& what, const std::vector<std::string>& option_lines,
                                     const std::vector<std::string>& fragments) {
    std::string message;
    try {
        catchem::ProcessRegistry::get_instance().validate_settings("gaschem", gaschem_process_config(option_lines));
    } catch (const std::invalid_argument& error) {
        message = error.what();
    }
    expect_init_error(what, message, fragments, "by the settings validator");
}

// Builds a Core from a runtime YAML that schedules gaschem with the given
// option lines and a mechanism directory that does not exist, and returns the
// message of the error that stops it.
static std::string configured_core_error(const std::vector<std::string>& option_lines) {
    const std::string temp_config = "test_gaschem_configured_core.yml";
    std::ofstream conf_writer(temp_config);
    conf_writer << "simulation:\n";
    conf_writer << "  name: gaschem-settings\n";
    conf_writer << "  species_filename: \"" << catchem::test::SOURCE_DIR
                << "/tests/fixtures/platform_integrity_species.yml\"\n";
    conf_writer << "  nx: 2\n  ny: 1\n  nz: 3\n  timestep: 60\n  nsteps: 1\n";
    conf_writer << "grid:\n  number_of_levels: 3\n";
    conf_writer << "diagnostics:\n  output:\n    enabled: false\n  collection:\n    enabled: false\n";
    conf_writer << "run_phases:\n  chemistry:\n    processes: [gaschem]\n";
    conf_writer << "processes:\n";
    conf_writer << "  gaschem:\n";
    conf_writer << "    activate: true\n";
    conf_writer << "    config_dir: \"/nonexistent/mechanism/\"\n";
    for (const auto& line : option_lines) {
        conf_writer << "    " << line << "\n";
    }
    conf_writer.close();
    std::string message;
    try {
        catchem::Core core(temp_config);
    } catch (const std::exception& error) {
        message = error.what();
    }
    std::remove(temp_config.c_str());
    return message;
}

// gaschem options are read by GasChemProcess::parse_settings. It is registered
// as the settings validator, which the core runs before init(), and init()
// calls it again for a process that is initialized directly.
static void run_settings_suite(const std::string& mechanism_dir) {
    using Policy = catchem::GasChemSettings::NegativeConcentrationPolicy;
    auto& registry = catchem::ProcessRegistry::get_instance();

    // 1. Defaults, then every option set next to framework keys, which are not options.
    const auto defaults = catchem::GasChemProcess::parse_settings(gaschem_process_config({}));
    assert(defaults.config_file.empty() && defaults.config_dir.empty() && defaults.columns_per_batch == 1 &&
           defaults.negative_concentration_policy == Policy::Reject && "gaschem option defaults");
    const auto all_set = gaschem_process_config({"diagnostics: true", "diag_species: [o3]",
                                                 "config_file: mechanism.yaml", "config_dir: mechanism/",
                                                 "columns_per_batch: 8", "negative_concentration_policy: clamp"});
    registry.validate_settings("gaschem", all_set);
    const auto parsed = catchem::GasChemProcess::parse_settings(all_set);
    assert(parsed.config_file == "mechanism.yaml" && parsed.config_dir == "mechanism/" &&
           parsed.columns_per_batch == 8 && parsed.negative_concentration_policy == Policy::Clamp &&
           "every gaschem option must be read from its key");

    // 2. The registered validator rejects a bad value, a misspelt key and a nested block.
    expect_settings_rejected("columns_per_batch: 0", {"columns_per_batch: 0"},
                             {"columns_per_batch must be >= 1, got 0"});
    expect_settings_rejected("columns_per_batch: many", {"columns_per_batch: many"},
                             {"not an integer", "columns_per_batch = many"});
    expect_settings_rejected("negative_concentration_policy: bogus", {"negative_concentration_policy: bogus"},
                             {"negative_concentration_policy must be reject or clamp, got bogus"});
    expect_settings_rejected(
        "misspelt column_per_batch", {"column_per_batch: 4"},
        {"Unknown option(s) for process 'gaschem': column_per_batch.",
         "Accepted options: columns_per_batch, config_dir, config_file, negative_concentration_policy."});
    expect_settings_rejected("nested solver block", {"solver:", "  rtol: 1.0e-6"},
                             {"Unknown option(s) for process 'gaschem': solver/rtol."});

    // 3. A process initialized directly, without the core, rejects an unknown key too.
    test_invalid_gaschem_option(mechanism_dir, "column_per_batch: 4");

    // 4. Both options take effect from a legacy top-level "process:" block: the
    // harness asserts the batch size, and one floored input shows clamp is on.
    BatchCaseOptions modern;
    modern.negative_concentration_policy = "clamp";
    modern.conc_override = CellValue{4, 2, 1, -1.0e-12};
    BatchCaseOptions legacy = modern;
    legacy.process_block = "process";
    CaseStats modern_stats;
    CaseStats legacy_stats;
    const auto conc_modern =
        run_multi_column_batch_case("processes: block", mechanism_dir, 4, 4, modern, &modern_stats);
    const auto conc_legacy =
        run_multi_column_batch_case("legacy process: block", mechanism_dir, 4, 4, legacy, &legacy_stats);
    assert(modern_stats.inputs == 1 && legacy_stats.inputs == 1 && "clamp must be in effect from either block");
    assert(conc_legacy == conc_modern && "The legacy process: block must solve bit-for-bit like processes:");

    // 5. A Core built from the runtime YAML stops on the bad option. The
    // mechanism directory does not exist, so the same configuration without the
    // misspelt key gets as far as MICM and fails there instead.
    expect_init_error("Core with a misspelt gaschem option", configured_core_error({"column_per_batch: 4"}),
                      {"Unknown option(s) for process 'gaschem': column_per_batch."}, "before MICM is constructed");
    expect_init_error("Core with a missing mechanism", configured_core_error({"columns_per_batch: 4"}),
                      {"failed to initialize MICM"}, "by MICM");
}

static void run_multi_column_batch_suite(const std::string& mechanism_dir) {
    constexpr int n_cols = 6;
    constexpr int n_levels = 3;
    constexpr int n_species = 4;

    // 1. Default (omitted columns_per_batch -> 1)
    const auto conc_default = run_multi_column_batch_case("default (omitted -> 1)", mechanism_dir, std::nullopt, 1);

    // 2. Explicit batch size 1 (6 batches of 1 column)
    CaseStats reference_stats;
    const auto conc_b1 = run_multi_column_batch_case("columns_per_batch=1", mechanism_dir, 1, 1, {}, &reference_stats);
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
    const auto conc_b2 = run_multi_column_batch_case("columns_per_batch=2 (even divisor)", mechanism_dir, 2, 2);

    // 4. Batch size 4 (non-divisor partial final batch: batch 0 has 4 cols, batch 1 has 2 valid + 2 padded cols)
    const auto conc_b4 = run_multi_column_batch_case("columns_per_batch=4 (partial final batch)", mechanism_dir, 4, 4);

    // 5. Batch size 6 (whole domain in 1 batch)
    const auto conc_b6 = run_multi_column_batch_case("columns_per_batch=6 (full domain)", mechanism_dir, 6, 6);

    // Compare every batch size against the columns_per_batch=1 reference. MICM controls its step size
    // with one error norm over all cells of a batch, so two batch sizes may step the same cell differently,
    // each to within the solver's local error tolerance of atol + rtol * |y|. The bound is ten times that
    // tolerance; the factor of ten is the allowance for comparing two independently stepped solutions over a
    // whole timestep. rtol is the value the solver ran with, which is MICM's default of 1e-6 because the
    // mechanism and the process set none. atol is the largest absolute tolerance the mechanism declares
    // (1e-12 mol/m3), converted to ppmv at the lowest air density in the case, where it is largest. O2
    // declares none, but oxygen conservation ties its error to the other three, so the same bound applies.
    const double rtol = reference_stats.solver_relative_tolerance;
    assert(rtol > 0.0 && rtol <= 1.0e-4 &&
           "The solver must report a relative tolerance tight enough for this comparison to mean something");
    constexpr double air_mw_kg_mol = catchem::constants::AIR_MW * 1.0e-3;
    const double lowest_air_density_mol = batch_case_air_density(n_cols - 1, n_levels - 1) / air_mw_kg_mol;
    const double atol_ppmv = largest_declared_absolute_tolerance(mechanism_dir) / lowest_air_density_mol * 1.0e6;
    constexpr double independent_solutions_factor = 10.0;
    std::cout << "Batch comparison bound: " << independent_solutions_factor << " * (" << atol_ppmv << " ppmv + " << rtol
              << " * |ref|)" << std::endl;
    for (const auto& [batch_label, candidate] : std::vector<std::pair<std::string, std::vector<double>>>{
             {"columns_per_batch=2", conc_b2}, {"columns_per_batch=4", conc_b4}, {"columns_per_batch=6", conc_b6}}) {
        double max_abs_diff = 0.0;
        double max_rel_diff = 0.0;
        for (std::size_t i = 0; i < conc_b1.size(); ++i) {
            const double ref = conc_b1[i];
            const double val = candidate[i];
            const double abs_diff = std::abs(val - ref);
            const double tol = independent_solutions_factor * (atol_ppmv + rtol * std::abs(ref));
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
    test_invalid_gaschem_option(mechanism_dir, "columns_per_batch: 0");
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

        // Case B: MUSICA v1 single-file mechanism passed via the new config_file key. The v1
        // mechanism solves for N2 as well, so cases B and C use a species list that includes it.
        std::string v1_file =
            std::string(catchem::test::SOURCE_DIR) + "/src/external/musica/configs/v1/chapman/config.yaml";
        assert(file_exists(v1_file) && "ERROR: Could not find gaschem v1 config file explicitly!");
        run_coupled_case("v1 single file (config_file)", "config_file", v1_file, "chapman_v1.yml", 5);

        // Case C: v1 mechanism directory exercising config.yaml auto-detection.
        std::string v1_dir = std::string(catchem::test::SOURCE_DIR) + "/src/external/musica/configs/v1/chapman/";
        assert(file_exists(v1_dir + "config.yaml") && "ERROR: Could not find gaschem v1 config dir explicitly!");
        run_coupled_case("v1 directory (config.yaml auto-detect)", "config_dir", v1_dir, "chapman_v1.yml", 5);

        // Case D: Multi-column batching suite across batch sizes (1, 2, 4, 6) and invalid batch rejection.
        run_multi_column_batch_suite(v0_dir);

        // Case E: negative_concentration_policy (reject, clamp) and the checks it must not relax.
        run_negative_concentration_policy_suite(v0_dir);

        // Case F: species resolution through the mechanism definition (case, aliases, duplicates).
        run_species_resolution_suite(v0_dir, v1_file);

        // Case G: separate clamp warnings for inputs and solver outputs, and the solver-output guard.
        run_clamp_warning_suite(v0_dir);

        // Case H: option parsing and validation, unknown keys, and the legacy process: block.
        run_settings_suite(v0_dir);

        std::cout << "\n==========================================" << std::endl;
        std::cout << "COUPLED INTEGRATION TESTS PASSED (v0 dir, v1 file, v1 auto-detect, batching, negative policy, "
                     "species resolution, clamp warnings, settings)!"
                  << std::endl;
        std::cout << "==========================================\n" << std::endl;
    }
    Kokkos::finalize();
    return 0;
}
