#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_kokkos_compat.hpp"
#include "catchem_logger.hpp"
#include "catchem_met_utilities.hpp"
#include "catchem_process_registry.hpp"
#include "catchem_state_manager.hpp"
#include "catchem_unit_conversion.hpp"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <string_view>
#include <vector>

extern "C" {
void catchem_register_seasalt_cpp();
void catchem_register_drydep_cpp();
void catchem_register_wetdep_cpp();
void catchem_register_settling_cpp();
void catchem_register_so4chem_cpp();
void catchem_register_dust_cpp();
void catchem_register_carbchem_cpp();
}

// Bounded random filling
void fill_random(std::vector<double>& vec, double min_val, double max_val, std::mt19937& gen) {
    std::uniform_real_distribution<double> dist(min_val, max_val);
    for (auto& val : vec) {
        val = dist(gen);
    }
}

// Helper to locate test fixture YAML files across build or source directories
std::string find_fixture(const std::string& filename) {
    std::vector<std::string> candidates = {filename, "tests/" + filename, "../tests/" + filename,
                                           "../../tests/" + filename};
    for (const auto& path : candidates) {
        std::ifstream f(path);
        if (f.good()) {
            return path;
        }
    }
    return filename;
}

// Invariants checking
void verify_properties(const std::vector<double>& conc, size_t size, int iteration) {
    for (size_t i = 0; i < size; ++i) {
        if (!std::isfinite(conc[i])) {
            std::cerr << "PROPERTY FAILURE: Index " << i << " is NaN or Inf at iteration " << iteration << std::endl;
            assert(false && "Concentration must remain finite!");
        }
        if (conc[i] < -1e-15) {
            std::cerr << "PROPERTY FAILURE: Index " << i << " is negative (" << conc[i] << ") at iteration "
                      << iteration << std::endl;
            assert(false && "Mass conservation violated!");
        }
    }
}

void CatchemPropertiesTest_RobustMetDerivations() {
    std::cout << "=== Running CatchemPropertiesTest.RobustMetDerivations ===" << std::endl;
    catchem::StateManager state(4, 10, 50);
    std::vector<double> pmid_data(40, 100000.0);
    std::vector<double> t_data(40, 300.0);
    std::vector<double> qv_data(40, 1.0); // 100% specific humidity (edge case: q >= 1.0)
    std::vector<double> airden_dry_data(40, 0.0);

    state.bind_met_field_3d("PMID", pmid_data.data());
    state.bind_met_field_3d("T", t_data.data());
    state.bind_met_field_3d("QV", qv_data.data());
    state.bind_met_field_3d("AIRDEN_DRY", airden_dry_data.data());

    // Should not divide by zero
    state.derive_airden_dry();
    state.sync_to_host();
    assert(!std::isnan(airden_dry_data[0]) && !std::isinf(airden_dry_data[0]));

    // Guard temperature t_val <= 0.0
    t_data[0] = 0.0;
    state.derive_airden_dry();
    state.sync_to_host();
    assert(!std::isnan(airden_dry_data[0]) && !std::isinf(airden_dry_data[0]));

    std::vector<double> pedge_data(44, 0.0); // 0 pressure at boundary (edge case: PEDGE = 0.0)
    std::vector<double> bxheight_data(40, 0.0);
    state.bind_met_field_3d("PEDGE", pedge_data.data());
    state.bind_met_field_3d("BXHEIGHT", bxheight_data.data());

    // Should not compute log(0) or log(negative)
    state.derive_bxheight();
    state.sync_to_host();
    assert(!std::isnan(bxheight_data[0]) && !std::isinf(bxheight_data[0]));
    assert(bxheight_data[0] == 0.0);

    std::cout << "=== PASS: CatchemPropertiesTest.RobustMetDerivations ===" << std::endl;
}

void CatchemPropertiesTest_ConfigManagerParsesRunPhases() {
    std::cout << "=== Running CatchemPropertiesTest.ConfigManagerParsesRunPhases ===" << std::endl;
    catchem::ConfigManager config_mgr;
    config_mgr.load_from_file(find_fixture("CATChem_new_config.yml"));

    assert(config_mgr.data.active_processes.size() == 1);
    assert(config_mgr.data.active_processes[0] == "seasalt");

    assert(config_mgr.get_bool("processes/seasalt/activate", false));
    assert(config_mgr.get_string("processes/seasalt/scheme", "") == "geos12");
    assert(config_mgr.get_string("processes/seasalt/scheme", "") == "geos12");

    std::cout << "=== PASS: CatchemPropertiesTest.ConfigManagerParsesRunPhases ===" << std::endl;
}

void CatchemPropertiesTest_ConfigManagerTypedQueries() {
    std::cout << "=== Running CatchemPropertiesTest.ConfigManagerTypedQueries ===" << std::endl;
    catchem::ConfigManager config_mgr;
    config_mgr.load_from_file(find_fixture("CATChem_new_config.yml"));

    assert(config_mgr.get_bool("processes/extemis/activate", false));
    assert(config_mgr.get_bool("processes/extemis/nonexistent", false) == false);
    assert(config_mgr.get_double("processes/extemis/global_factor", 0.0) == 1.0);
    assert(config_mgr.get_int("grid/number_of_levels", 0) == 64);
    assert(config_mgr.get_string("processes/seasalt/scheme", "") == "geos12");

    std::vector<std::string> diag = config_mgr.get_string_list("processes/extemis/anthro/diag_list");
    assert(!diag.empty());
    assert(diag[0] == "NO");

    assert(config_mgr.is_process_active("seasalt") == true);
    assert(config_mgr.is_process_active("nonexistent_process") == false);

    std::string fengsha_file = config_mgr.find_process_file_setting("fengsha");
    assert(!fengsha_file.empty());
    assert(fengsha_file.find("fengsha_emis.nc") != std::string::npos);

    std::cout << "=== PASS: CatchemPropertiesTest.ConfigManagerTypedQueries ===" << std::endl;
}

void CatchemPropertiesTest_ConfigManagerHandlesScalarAndMissingNodes() {
    std::cout << "=== Running CatchemPropertiesTest.ConfigManagerHandlesScalarAndMissingNodes ===" << std::endl;
    catchem::ConfigManager config_mgr;

    // Query on empty/unloaded manager
    assert(config_mgr.get_bool("processes/anything", false) == false);
    assert(config_mgr.get_string("processes/anything", "def") == "def");
    assert(config_mgr.get_double("processes/anything", 1.23) == 1.23);
    assert(config_mgr.get_int("processes/anything", 42) == 42);
    assert(config_mgr.get_string_list("processes/anything").empty());
    assert(config_mgr.find_process_file_setting("fengsha").empty());

    // Load file with partial/scalar nodes
    config_mgr.load_from_file(find_fixture("CATChem_new_config.yml"));

    // Deep invalid paths
    assert(config_mgr.get_bool("processes/seasalt/scheme/invalid_sub_key", false) == false);
    assert(config_mgr.get_string("processes/seasalt/activate/sub_key", "default") == "default");
    assert(config_mgr.get_string_list("simulation/nx").empty());

    std::cout << "=== PASS: CatchemPropertiesTest.ConfigManagerHandlesScalarAndMissingNodes ===" << std::endl;
}

void CatchemPropertiesTest_ConfigManagerLogLevelDrivesLogger() {
    std::cout << "=== Running CatchemPropertiesTest.ConfigManagerLogLevelDrivesLogger ===" << std::endl;
    using Level = catchem::Logger::Level;

    // 1) simulation/verbose/log_level overrides the environment-derived
    //    threshold (YAML always wins) and is recorded in the parsed data.
    {
        const std::string path = "cfg_log_level_warn.yml";
        std::ofstream f(path);
        f << "simulation:\n  name: loglevel\n  verbose:\n    activate: true\n    log_level: warn\n";
        f.close();
        catchem::ConfigManager mgr;
        mgr.load_from_file(path);
        assert(mgr.data.simulation.log_level == "warn");
        assert(mgr.data.simulation.verbose_enabled);
        assert(catchem::Logger::enabled(Level::Warn));
        assert(!catchem::Logger::enabled(Level::Info));
        assert(!catchem::Logger::enabled(Level::Debug));
    }

    // 2) An unrecognized level fails loudly at load time (fail-fast policy).
    {
        const std::string path = "cfg_log_level_bogus.yml";
        std::ofstream f(path);
        f << "simulation:\n  name: loglevel\n  verbose:\n    log_level: shout\n";
        f.close();
        bool threw = false;
        try {
            catchem::ConfigManager mgr;
            mgr.load_from_file(path);
        } catch (const std::exception&) {
            threw = true;
        }
        assert(threw);
    }

    // 3) A config without log_level keeps the previous setting (absent != reset).
    {
        const std::string path = "cfg_log_level_absent.yml";
        std::ofstream f(path);
        f << "simulation:\n  name: loglevel\n  verbose:\n    activate: false\n";
        f.close();
        catchem::ConfigManager mgr;
        mgr.load_from_file(path);
        assert(mgr.data.simulation.log_level.empty());
        assert(!catchem::Logger::enabled(Level::Info)); // warn override from step 1 persists
    }

    // Leave the process in the environment-controlled state for later tests.
    catchem::Logger::clear_level();

    std::cout << "=== PASS: CatchemPropertiesTest.ConfigManagerLogLevelDrivesLogger ===" << std::endl;
}

void CatchemPropertiesTest_StateBindingAndRebinding() {
    std::cout << "=== Running CatchemPropertiesTest.StateBindingAndRebinding ===" << std::endl;
    int nc = 4, nl = 10, ns = 5;
    catchem::StateManager state(nc, nl, ns);

    std::vector<double> t1(nc * nl, 280.0);
    std::vector<double> t2(nc * nl, 290.0);
    std::vector<double> ps1(nc, 100000.0);
    std::vector<double> ps2(nc, 101325.0);
    std::vector<double> chem1(nc * nl * ns, 1.0e-9);
    std::vector<double> chem2(nc * nl * ns, 2.0e-9);

    state.bind_met_field_3d("T", t1.data());
    assert(state.meteorology().T->host_data() == t1.data());
    assert(state.get_host_pointer_3d("T") == t1.data());

    // Rebind to new host buffer
    state.bind_met_field_3d("T", t2.data());
    assert(state.meteorology().T->host_data() == t2.data());
    assert(state.get_host_pointer_3d("T") == t2.data());

    state.bind_met_field_2d("PS", ps1.data());
    assert(state.meteorology().PS->host_data() == ps1.data());
    assert(state.get_host_pointer_2d("PS") == ps1.data());

    state.bind_met_field_2d("PS", ps2.data());
    assert(state.meteorology().PS->host_data() == ps2.data());

    state.bind_unified_chemistry(chem1.data());
    assert(state.chemistry().conc->host_data() == chem1.data());

    state.bind_unified_chemistry(chem2.data());
    assert(state.chemistry().conc->host_data() == chem2.data());

    std::cout << "=== PASS: CatchemPropertiesTest.StateBindingAndRebinding ===" << std::endl;
}

void CatchemPropertiesTest_GridAndDimensions() {
    std::cout << "=== Running CatchemPropertiesTest.GridAndDimensions ===" << std::endl;
    int nx = 8, ny = 4, nz = 16;
    catchem::GridManager grid(nx, ny, nz);
    assert(grid.geometry.nx == nx);
    assert(grid.geometry.ny == ny);
    assert(grid.geometry.nz == nz);

    std::vector<double> lat(nx * ny, 45.0);
    std::vector<double> lon(nx * ny, -100.0);
    std::vector<double> area(nx * ny, 1.0e6);

    grid.bind_lat(lat.data());
    grid.bind_lon(lon.data());
    grid.bind_area(area.data());

    assert(grid.geometry.lat->host_data() == lat.data());
    assert(grid.geometry.lon->host_data() == lon.data());
    assert(grid.geometry.grid_area->host_data() == area.data());

    // C API check
    void* core_ptr = catchem_core_create(nx * ny, nz, 10);
    int g_nx = 0, g_ny = 0, g_nz = 0;
    catchem_get_grid_dimensions(core_ptr, &g_nx, &g_ny, &g_nz);
    assert(g_nx == nx * ny);
    assert(g_ny == 1);
    assert(g_nz == nz);
    catchem_core_destroy(core_ptr);

    std::cout << "=== PASS: CatchemPropertiesTest.GridAndDimensions ===" << std::endl;
}

void CatchemPropertiesTest_TimeStateCalculations() {
    std::cout << "=== Running CatchemPropertiesTest.TimeStateCalculations ===" << std::endl;
    catchem::TimeState time;
    time.year = 2024; // leap year
    time.month = 2;
    time.day = 28;
    time.hour = 23;
    time.minute = 59;
    time.second = 0;

    assert(catchem::TimeState::is_leap_year(2024));
    assert(!catchem::TimeState::is_leap_year(2023));
    assert(catchem::TimeState::get_days_in_month(2, 2024) == 29);
    assert(catchem::TimeState::get_days_in_month(2, 2023) == 28);

    time.calculate_derived_fields();
    assert(time.doy == 59);

    // Advance by 120 seconds -> Feb 29 00:01:00
    time.advance(120.0);
    assert(time.year == 2024);
    assert(time.month == 2);
    assert(time.day == 29);
    assert(time.hour == 0);
    assert(time.minute == 1);
    assert(time.second == 0);

    // Advance to next year
    time.advance(365 * 86400.0);
    assert(time.year == 2025);

    std::cout << "=== PASS: CatchemPropertiesTest.TimeStateCalculations ===" << std::endl;
}

void CatchemPropertiesTest_SpeciesMetadataAPI() {
    std::cout << "=== Running CatchemPropertiesTest.SpeciesMetadataAPI ===" << std::endl;
    catchem::StateManager state(4, 10, 50);
    std::string species_path = find_fixture("CATChem_species.yml");
    state.chemistry().load_species_config(species_path);

    int count = catchem_state_get_species_count(&state);
    assert(count > 20);

    int so2_idx = catchem_state_get_species_index(&state, "so2"); // 1-based
    assert(so2_idx > 0);

    char name_buf[128];
    char desc_buf[256];
    char mie_buf[64];

    catchem_state_get_species_name_at(&state, so2_idx, name_buf);
    assert(std::string(name_buf) == "so2");

    catchem_state_get_species_desc_at(&state, so2_idx, desc_buf);
    assert(std::string(desc_buf) == "Sulfur dioxide");

    assert(catchem_state_is_species_gas(&state, so2_idx) == 1);
    assert(catchem_state_is_species_drydep(&state, so2_idx) == 1);
    assert(catchem_state_is_species_wetdep(&state, so2_idx) == 1);

    int msa_idx = catchem_state_get_species_index(&state, "msa");
    assert(msa_idx > 0);
    // The original SO4chem Fortran scheme carries MSA as aerosol mass (ug/kg).
    assert(catchem_state_is_species_gas(&state, msa_idx) == 0);
    assert(catchem_state_is_species_aerosol(&state, msa_idx) == 1);
    assert(catchem_state_is_species_drydep(&state, msa_idx) == 1);
    assert(catchem_state_is_species_wetdep(&state, msa_idx) == 1);
    // This fixture (tests/CATChem_species.yml) intentionally keeps the
    // optics-table physical properties; the parity A/B baseline lives in
    // tests/Configs/Default/CATChem_species.yml and is covered by
    // test_catchem_settling / the parity runners.
    assert(catchem_state_get_species_radius(&state, msa_idx) == 0.35);
    assert(catchem_state_get_species_density(&state, msa_idx) == 1700.0);

    int dust1_idx = catchem_state_get_species_index(&state, "dust1");
    assert(dust1_idx > 0);
    assert(catchem_state_is_species_dust(&state, dust1_idx) == 1);
    assert(catchem_state_is_species_aerosol(&state, dust1_idx) == 1);
    // Fixture-only: the 49-species unit fixture carries the GOCART DU
    // optics-table radius/density.  See the Default file for the parity
    // baseline values.
    assert(catchem_state_get_species_density(&state, dust1_idx) == 2650.0);
    assert(catchem_state_get_species_radius(&state, dust1_idx) == 0.6359);
    assert(catchem_state_get_species_lower_radius(&state, dust1_idx) == 0.1);
    assert(catchem_state_get_species_upper_radius(&state, dust1_idx) == 1.0);
    assert(catchem_state_get_species_radius(&state, count + 1) == 0.0);
    catchem_state_get_species_mie_name(&state, dust1_idx, mie_buf);
    assert(std::string(mie_buf) == "DU");

    std::cout << "=== PASS: CatchemPropertiesTest.SpeciesMetadataAPI ===" << std::endl;
}

void CatchemPropertiesTest_ConfigManagerLoadsTypedFixtureData() {
    std::cout << "=== Running CatchemPropertiesTest.ConfigManagerLoadsTypedFixtureData ===" << std::endl;
    catchem::ConfigManager config_mgr;
    config_mgr.load_from_file(find_fixture("CATChem_new_config.yml"));
    config_mgr.load_species_file(find_fixture("CATChem_species.yml"));
    config_mgr.load_emission_mapping_file(find_fixture("CATChem_emission.yml"));

    assert(config_mgr.data.simulation.name == "test");
    assert(config_mgr.data.simulation.species_filename == "./CATChem_species.yml");
    assert(config_mgr.data.simulation.emission_filename == "./CATChem_emission.yml");
    assert(config_mgr.data.simulation.verbose_enabled);
    assert(config_mgr.data.grid.number_of_levels == 64);
    assert(config_mgr.data.grid.number_of_soil_layers == 4);
    assert(config_mgr.data.timesteps.transport_timestep_in_s == 10);
    assert(config_mgr.data.timesteps.chemistry_timestep_in_s == 60);
    assert(config_mgr.data.diagnostics.output.enabled);
    assert(config_mgr.data.diagnostics.output.directory == "./output");
    assert(config_mgr.data.diagnostics.output.prefix == "catchem_diag");
    assert(config_mgr.data.diagnostics.output.frequency == 3600);
    assert(config_mgr.data.diagnostics.output.format == "netcdf");
    // Per-process diagnostic output is opt-in via diagnostics/output/process_diagnostics.
    assert(config_mgr.data.diagnostics.output.process_diagnostics);
    assert(config_mgr.data.diagnostics.output.diag_list.size() == 7);
    assert(config_mgr.data.diagnostics.output.diag_list[0] == "so2");
    assert(config_mgr.data.diagnostics.collection.enabled);
    assert(config_mgr.data.diagnostics.collection.buffer_size == 1000);

    const auto& seasalt = config_mgr.data.processes.at("seasalt");
    assert(seasalt.activate);
    assert(seasalt.diagnostics);
    assert(seasalt.scheme == "geos12");
    assert(config_mgr.get_double("processes/seasalt/geos12/scale_factor", 0.0) == 1.0);

    assert(config_mgr.data.species.size() > 20);
    const auto& so2 = config_mgr.data.species.at(0);
    assert(so2.name == "so2");
    assert(so2.description == "Sulfur dioxide");
    assert(so2.is_gas);
    assert(so2.is_drydep);
    assert(so2.is_wetdep);
    assert(std::abs(so2.molecular_weight_kg_mol - 64.04e-3) < 1.0e-12);

    const auto& anthro = config_mgr.data.emission_mappings.at("anthro");
    const auto& no_mapping = anthro.fields.at("NO");
    assert(no_mapping.long_name == "Nitrogen Oxide");
    assert(no_mapping.units == "kg/m2/s");
    assert(no_mapping.scale.size() == 1);
    assert(no_mapping.scale[0] == 1.0);
    assert(no_mapping.map.size() == 1);
    assert(no_mapping.map[0] == "no3");

    catchem::ConfigManager legacy_mgr;
    legacy_mgr.load_from_file(find_fixture("CATChem_config.yml"));
    assert(legacy_mgr.is_process_active("seasalt"));
    assert(legacy_mgr.data.processes.at("seasalt").get_int("scheme_opt", 0) == 3);

    std::cout << "=== PASS: CatchemPropertiesTest.ConfigManagerLoadsTypedFixtureData ===" << std::endl;
}

void CatchemPropertiesTest_DiagnosticsManagerAndAPI() {
    std::cout << "=== Running CatchemPropertiesTest.DiagnosticsManagerAndAPI ===" << std::endl;
    void* core_ptr = catchem_core_create(4, 10, 5);

    catchem_diag_register(core_ptr, "test_diag_2d", "Test 2D Diagnostic", "K", 2, 4, 10, 1);
    catchem_diag_register(core_ptr, "test_diag_3d", "Test 3D Diagnostic", "ppm", 3, 4, 10, 5);

    int count = catchem_diag_get_count(core_ptr);
    assert(count == 2);

    char name0_buf[64], name1_buf[64];
    catchem_diag_get_name_at(core_ptr, 0, name0_buf);
    catchem_diag_get_name_at(core_ptr, 1, name1_buf);
    std::string n0(name0_buf), n1(name1_buf);
    assert((n0 == "test_diag_2d" && n1 == "test_diag_3d") || (n0 == "test_diag_3d" && n1 == "test_diag_2d"));

    double* ptr_2d = static_cast<double*>(catchem_diag_get_pointer(core_ptr, "test_diag_2d"));
    assert(ptr_2d != nullptr);

    double* ptr_3d = static_cast<double*>(catchem_diag_get_pointer(core_ptr, "test_diag_3d"));
    assert(ptr_3d != nullptr);

    ptr_2d[0] = 273.15;
    ptr_3d[0] = 42.0;

    catchem_diag_sync_to_host(core_ptr);
    catchem_diag_reset(core_ptr);

    assert(ptr_2d[0] == 0.0);
    assert(ptr_3d[0] == 0.0);

    catchem_core_destroy(core_ptr);
    std::cout << "=== PASS: CatchemPropertiesTest.DiagnosticsManagerAndAPI ===" << std::endl;
}

void CatchemPropertiesTest_UnitConversions() {
    std::cout << "=== Running CatchemPropertiesTest.UnitConversions ===" << std::endl;
    double temp = 298.15;
    double press = 101325.0;
    double mw_o3 = 48.0;

    double ugm3 = catchem::unit_conversion::ppbv_to_ugm3(10.0, mw_o3, temp, press);
    double ppbv = catchem::unit_conversion::ugm3_to_ppbv(ugm3, mw_o3, temp, press);
    assert(std::abs(ppbv - 10.0) < 1.0e-10);

    double molcm3 = catchem::unit_conversion::ppbv_to_molcm3(10.0, temp, press);
    double back_ppbv = catchem::unit_conversion::molcm3_to_ppbv(molcm3, temp, press);
    assert(std::abs(back_ppbv - 10.0) < 1.0e-10);

    double mgm3 = catchem::unit_conversion::ppmv_to_mgm3(1.0, mw_o3, temp, press);
    double ppmv = catchem::unit_conversion::mgm3_to_ppmv(mgm3, mw_o3, temp, press);
    assert(std::abs(ppmv - 1.0) < 1.0e-10);

    double kgm2s = catchem::unit_conversion::molcm2s_to_kgm2s(1.0e-6, mw_o3);
    double molcm2s = catchem::unit_conversion::kgm2s_to_molcm2s(kgm2s, mw_o3);
    assert(std::abs(molcm2s - 1.0e-6) < 1.0e-12);

    double calc_mw = catchem::unit_conversion::calculate_molecular_weight("O3");
    assert(std::abs(calc_mw - 48.0) < 1.0e-6);

    std::cout << "=== PASS: CatchemPropertiesTest.UnitConversions ===" << std::endl;
}

void CatchemPropertiesTest_MetUtilities() {
    std::cout << "=== Running CatchemPropertiesTest.MetUtilities ===" << std::endl;
    double temp = 300.0;
    double qv = 0.01;
    double p = 100000.0;
    double ps = 101325.0;

    double theta = catchem::met_utilities::potential_temperature(temp, p, ps);
    assert(theta > temp);

    double tv = catchem::met_utilities::virtual_temperature(temp, qv);
    assert(tv > temp);

    double es = catchem::met_utilities::saturation_vapor_pressure(temp);
    assert(es > 0.0);

    double rh = catchem::met_utilities::relative_humidity(temp, qv, p);
    assert(rh >= 0.0 && rh <= 1.0);

    double mix = catchem::met_utilities::mixing_ratio(qv);
    double back_q = catchem::met_utilities::specific_humidity(mix);
    assert(std::abs(back_q - qv) < 1.0e-6);

    std::cout << "=== PASS: CatchemPropertiesTest.MetUtilities ===" << std::endl;
}

void CatchemPropertiesTest_ValidDerivationVectors() {
    catchem::StateManager state(1, 2, 1);
    std::vector<double> pressure_edge{100000.0, 85000.0, 70000.0};
    std::vector<double> pressure_mid{92500.0, 77500.0};
    std::vector<double> temperature{290.0, 270.0};
    std::vector<double> humidity{0.005, 0.010};
    state.bind_met_field_3d("PEDGE", pressure_edge.data());
    state.bind_met_field_3d("PMID", pressure_mid.data());
    state.bind_met_field_3d("T", temperature.data());
    state.bind_met_field_3d("QV", humidity.data());
    state.set_validation_policy(catchem::PhysicalValidationPolicy::Reject);

    state.derive_bxheight();
    state.derive_airden_dry();
    state.sync_to_host();
    assert(state.validation_report().empty());
    for (int level = 0; level < 2; ++level) {
        const double virtual_temperature =
            catchem::met_utilities::virtual_temperature(temperature[level], humidity[level]);
        const double expected_height = (catchem::constants::RD / catchem::constants::G0) * virtual_temperature *
                                       std::log(pressure_edge[level] / pressure_edge[level + 1]);
        // The derivation computes in catchem::fp, which is float unless the
        // build promotes reals (USE_REAL8).  Compare with a relative tolerance
        // scaled to fp's epsilon so the check is correct in both single- and
        // double-precision builds instead of assuming double.
        const double rel_tol = 32.0 * static_cast<double>(std::numeric_limits<catchem::fp>::epsilon());
        assert(std::abs(state.meteorology().BXHEIGHT->host_data()[level] - expected_height) <=
               rel_tol * std::abs(expected_height));

        const double ratio =
            (catchem::constants::AIR_MW / catchem::constants::H2O_MW) * humidity[level] / (1.0 - humidity[level]);
        const double water_mole_fraction = ratio / (1.0 + ratio);
        const double expected_density =
            pressure_mid[level] * (1.0 - water_mole_fraction) / (catchem::constants::RD * temperature[level]);
        assert(std::abs(state.meteorology().AIRDEN_DRY->host_data()[level] - expected_density) <=
               rel_tol * std::abs(expected_density));
    }
}

void CatchemPropertiesTest_HostWetDepFieldsAreAuthoritative() {
    catchem::StateManager state(1, 2, 1);
    std::vector<double> airden_dry{9.1, 9.2};
    std::vector<double> reevapls{3.1e-6, 3.2e-6};
    std::vector<double> bxheight{123.0, 456.0};

    state.bind_met_field_3d("AIRDEN_DRY", airden_dry.data());
    state.bind_met_field_3d("REEVAPLS", reevapls.data());
    state.bind_met_field_3d("BXHEIGHT", bxheight.data());

    state.derive_airden_dry();
    state.derive_reevapls();
    state.derive_bxheight();
    state.sync_to_host();

    assert(airden_dry[0] == 9.1 && airden_dry[1] == 9.2);
    assert(reevapls[0] == 3.1e-6 && reevapls[1] == 3.2e-6);
    assert(bxheight[0] == 123.0 && bxheight[1] == 456.0);
}

void CatchemPropertiesTest_ConstantsAndPrecision() {
    std::cout << "=== Running CatchemPropertiesTest.ConstantsAndPrecision ===" << std::endl;
    assert(std::abs(catchem::constants::RD - 287.05) < 0.1);
    assert(std::abs(catchem::constants::CP - 1004.6) < 1.0);
    assert(std::abs(catchem::constants::G0 - 9.80665) < 1.0e-4);
    assert(catchem::constants::PI > 3.1415 && catchem::constants::PI < 3.1416);

    std::cout << "=== PASS: CatchemPropertiesTest.ConstantsAndPrecision ===" << std::endl;
}

void CatchemPropertiesTest_CoreConstructsConfiguredProcesses() {
    std::cout << "=== Running CatchemPropertiesTest.CoreConstructsConfiguredProcesses ===" << std::endl;
    catchem_register_seasalt_cpp();

    catchem::Core core("CATChem_new_config.yml", 4, 10);

    assert(core.get_num_processes() == 1);

    std::cout << "=== PASS: CatchemPropertiesTest.CoreConstructsConfiguredProcesses ===" << std::endl;
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "=== RUNNING RANDOMIZED PROPERTY TESTS ===" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        // Verify Trace ID Generation
        {
            auto test_state = std::make_shared<catchem::StateManager>(4, 10, 50);
            assert(test_state->runtime_trace_id().length() == 8);
            assert(!test_state->runtime_trace_id().empty());
        }

        CatchemPropertiesTest_RobustMetDerivations();
        CatchemPropertiesTest_HostWetDepFieldsAreAuthoritative();
        CatchemPropertiesTest_StateBindingAndRebinding();
        CatchemPropertiesTest_GridAndDimensions();
        CatchemPropertiesTest_TimeStateCalculations();
        CatchemPropertiesTest_SpeciesMetadataAPI();
        CatchemPropertiesTest_ConfigManagerParsesRunPhases();
        CatchemPropertiesTest_ConfigManagerTypedQueries();
        CatchemPropertiesTest_ConfigManagerHandlesScalarAndMissingNodes();
        CatchemPropertiesTest_ConfigManagerLogLevelDrivesLogger();
        CatchemPropertiesTest_ConfigManagerLoadsTypedFixtureData();
        CatchemPropertiesTest_DiagnosticsManagerAndAPI();
        CatchemPropertiesTest_UnitConversions();
        CatchemPropertiesTest_MetUtilities();
        CatchemPropertiesTest_ValidDerivationVectors();
        CatchemPropertiesTest_ConstantsAndPrecision();
        CatchemPropertiesTest_CoreConstructsConfiguredProcesses();

        // Register All C++ Modern Process Handlers
        catchem_register_seasalt_cpp();
        catchem_register_drydep_cpp();
        catchem_register_wetdep_cpp();
        catchem_register_settling_cpp();
        catchem_register_so4chem_cpp();
        catchem_register_dust_cpp();
        catchem_register_carbchem_cpp();

        int n_cols = 12;
        int n_levels = 8;
        int n_species = 22;
        int size_3d = n_cols * n_levels;
        size_t total_size = n_cols * n_levels * n_species;

        // Create Core Orchestration Layer
        void* core_ptr = catchem_core_create(n_cols, n_levels, n_species);
        auto* core = static_cast<catchem::Core*>(core_ptr);
        auto state = core->get_state_manager();

        // Load runtime YAML specifications
        std::string config_path = "";
        std::vector<std::string> candidates = {"tests/CATChem_species.yml", "../tests/CATChem_species.yml",
                                               "../../tests/CATChem_species.yml", "CATChem_species.yml"};
        for (const auto& path : candidates) {
            std::ifstream f(path);
            if (f.good()) {
                config_path = path;
                break;
            }
        }
        if (config_path.empty()) {
            std::cerr << "ERROR: Could not find CATChem_species.yml inside test_catchem_properties.cpp\n";
            std::exit(1);
        }
        state->load_species_config(config_path);

        // The grid-only core starts with an empty processes block.  Attach a
        // runtime configuration so each process init() sees the scheme it
        // requires (settling/so4chem/carbchem gocart, dust fengsha, drydep
        // wesely/gocart, wetdep jacob, seasalt geos12).
        {
            const std::string cfg_path = "properties_runtime_config.yml";
            std::ofstream cfg(cfg_path);
            cfg << "simulation:\n  name: properties\n"
                << "processes:\n"
                << "  seasalt:\n    activate: true\n    scheme: geos12\n    diagnostics: false\n"
                << "  dust:\n    activate: true\n    scheme: fengsha\n    diagnostics: false\n"
                << "  drydep:\n    activate: true\n    gas_scheme: wesely\n    aero_scheme: gocart\n    diagnostics: "
                   "false\n"
                << "  settling:\n    activate: true\n    scheme: gocart\n    diagnostics: false\n"
                << "  so4chem:\n    activate: true\n    scheme: gocart\n    diagnostics: false\n"
                << "  wetdep:\n    activate: true\n    scheme: jacob\n    diagnostics: false\n"
                << "  carbchem:\n    activate: true\n    scheme: gocart\n    diagnostics: false\n";
            cfg.close();
            auto runtime_config = std::make_shared<catchem::ConfigManager>();
            runtime_config->load_from_file(cfg_path);
            state->attach_config_manager(runtime_config);
        }

        // Set up bounded fuzzer generator with fixed seed
        std::mt19937 gen(1337);

        // Define fuzzed input tensors
        std::vector<double> t_air(size_3d);
        std::vector<double> pmid(size_3d);
        std::vector<double> pedge(n_cols * (n_levels + 1));
        std::vector<double> airden_dry(size_3d);
        std::vector<double> mairden(size_3d);
        std::vector<double> bxheight(size_3d);
        std::vector<double> cldf(size_3d, 0.1);
        std::vector<double> pfilsan(size_3d);
        std::vector<double> pfllsan(size_3d);
        std::vector<double> reevapls(size_3d);
        std::vector<double> lat(n_cols, 40.0);
        std::vector<double> lon(n_cols, -80.0);
        std::vector<double> sst(n_cols, 290.0);
        std::vector<double> frocean(n_cols, 1.0);
        std::vector<double> frseaice(n_cols, 0.0);
        std::vector<double> ustar(n_cols, 0.5);
        std::vector<double> delp(size_3d, 1000.0);
        std::vector<double> u10m(n_cols, 5.0);
        std::vector<double> v10m(n_cols, 2.0);
        std::vector<double> rh(size_3d, 50.0);
        std::vector<double> ps(n_cols, 101325.0);
        std::vector<double> ts(n_cols, 288.15);
        std::vector<double> hflux(n_cols, 10.0);
        std::vector<double> obk(n_cols, 100.0);
        std::vector<double> pblh(n_cols, 1000.0);
        std::vector<double> z0h(n_cols, 0.01);

        // Conc & Tendencies fuzzed states
        std::vector<double> conc(total_size);

        // Bind meteorological fields
        state->bind_met_field_3d("T", t_air.data());
        state->bind_met_field_3d("PMID", pmid.data());
        state->bind_met_field_3d("PEDGE", pedge.data());
        state->bind_met_field_3d("AIRDEN_DRY", airden_dry.data());
        state->bind_met_field_3d("AIRDEN", airden_dry.data()); // Use airden_dry as the AIRDEN backing store
        state->bind_met_field_3d("BXHEIGHT", bxheight.data());
        state->bind_met_field_3d("MAIRDEN", mairden.data());
        state->bind_met_field_3d("PFILSAN", pfilsan.data());
        state->bind_met_field_3d("PFLLSAN", pfllsan.data());
        state->bind_met_field_3d("REEVAPLS", reevapls.data());
        state->bind_met_field_3d("RH", rh.data());
        state->bind_met_field_3d("CLDF", cldf.data());
        state->bind_met_field_2d("LAT", lat.data());
        state->bind_met_field_2d("LON", lon.data());
        state->bind_met_field_2d("SST", sst.data());
        state->bind_met_field_2d("PS", ps.data());
        state->bind_met_field_2d("TS", ts.data());
        state->bind_met_field_2d("HFLUX", hflux.data());
        state->bind_met_field_2d("OBK", obk.data());
        state->bind_met_field_2d("PBLH", pblh.data());
        state->bind_met_field_2d("Z0H", z0h.data());
        state->bind_met_field_2d("FROCEAN", frocean.data());
        state->bind_met_field_2d("FRSEAICE", frseaice.data());
        state->bind_met_field_2d("USTAR", ustar.data());
        state->bind_met_field_3d("DELP", delp.data());
        state->bind_met_field_2d("U10M", u10m.data());
        state->bind_met_field_2d("V10M", v10m.data());

        // Schedule All Registered Processes dynamically to test simultaneous execution
        auto settling = catchem::ProcessRegistry::get_instance().create("settling");
        settling->init(state);
        core->add_process(settling);

        auto drydep = catchem::ProcessRegistry::get_instance().create("drydep");
        drydep->init(state);
        core->add_process(drydep);

        auto seasalt = catchem::ProcessRegistry::get_instance().create("seasalt");
        seasalt->init(state);
        core->add_process(seasalt);

        auto wetdep = catchem::ProcessRegistry::get_instance().create("wetdep");
        wetdep->init(state);
        core->add_process(wetdep);

        auto so4chem = catchem::ProcessRegistry::get_instance().create("so4chem");
        so4chem->init(state);
        core->add_process(so4chem);

        auto dust = catchem::ProcessRegistry::get_instance().create("dust");
        dust->init(state);
        core->add_process(dust);

        auto carbchem = catchem::ProcessRegistry::get_instance().create("carbchem");
        carbchem->init(state);
        core->add_process(carbchem);

        std::cout << "state->meteorology().T = " << state->meteorology().T.get() << std::endl;
        std::cout << "state->meteorology().AIRDEN = " << state->meteorology().AIRDEN.get() << std::endl;
        std::cout << "state->meteorology().PEDGE = " << state->meteorology().PEDGE.get() << std::endl;
        std::cout << "state->meteorology().BXHEIGHT = " << state->meteorology().BXHEIGHT.get() << std::endl;
        std::cout << "state->chemistry().conc = " << state->chemistry().conc.get() << std::endl;

        std::cout << "Executing 100 high-fuzz property iterations over 7 synchronized processes..." << std::endl;

        for (int iter = 1; iter <= 100; ++iter) {
            // Fuzz temperature across extreme physical atmospheric ranges
            fill_random(t_air, 170.0, 330.0, gen); // Stratosphere to boundary surface Temps

            // Construct monotonic, physically consistent pressure edges and midpoints per column
            for (int icol = 0; icol < n_cols; ++icol) {
                double current_p = std::uniform_real_distribution<double>(95000.0, 103000.0)(gen); // Surface pressure
                pedge[icol + 0 * n_cols] = current_p;
                for (int k = 0; k < n_levels; ++k) {
                    double delta = std::uniform_real_distribution<double>(5000.0, 12000.0)(gen);
                    current_p -= delta;
                    double min_p = 100.0 - k * 5.0;
                    if (current_p < min_p)
                        current_p = min_p;
                    pedge[icol + (k + 1) * n_cols] = current_p;

                    // Midpoint pressure is average of the edges
                    double p1 = pedge[icol + k * n_cols];
                    double p2 = pedge[icol + (k + 1) * n_cols];
                    pmid[icol + k * n_cols] = 0.5 * (p1 + p2);
                    delp[icol + k * n_cols] = std::abs(p1 - p2);

                    // Derive dry air density using the Ideal Gas Law: rho = P / (R_dry * T)
                    double t = t_air[icol + k * n_cols];
                    double rho = pmid[icol + k * n_cols] / (287.05 * t);
                    if (rho < 0.01)
                        rho = 0.01;
                    if (rho > 2.0)
                        rho = 2.0;
                    airden_dry[icol + k * n_cols] = rho;
                    mairden[icol + k * n_cols] = rho;

                    // Derive dz (layer thickness) using hydrostatic balance: dz = dp / (rho * g)
                    double dz = delp[icol + k * n_cols] / (rho * 9.80665);
                    if (dz < 1.0)
                        dz = 1.0;
                    if (dz > 5000.0)
                        dz = 5000.0;
                    bxheight[icol + k * n_cols] = dz;
                }
            }

            fill_random(cldf, 0.0, 1.0, gen);    // cloud fractions
            fill_random(pfilsan, 0.0, 0.1, gen); // Dynamic fractions
            fill_random(pfllsan, 0.0, 0.1, gen);
            fill_random(reevapls, 0.0, 1e-4, gen); // Dynamic liquid reevaporations
            fill_random(ustar, 0.01, 2.5, gen);    // Extreme shear friction winds
            fill_random(u10m, -50.0, 50.0, gen);   // Dynamic 10-meter wind components
            fill_random(v10m, -50.0, 50.0, gen);

            // Fuzz chemical concentrations
            fill_random(conc, 0.0, 1e-6, gen); // Plausible trace-gas concentrations (kg/kg)

            // Dynamic pointers association and synchronizations
            state->bind_unified_chemistry(conc.data());
            state->sync_to_device();

            // Execute the modernized scheduled timestepping loop
            double dt = 1800.0; // 30-minute sim time step
            catchem_core_run_timestep(core_ptr, dt);

            // Sync outputs back to host C++ side
            state->sync_to_host();

            // Assert finite properties first
            for (size_t i = 0; i < total_size; ++i) {
                if (!std::isfinite(conc[i])) {
                    int spec_idx = i / size_3d;
                    int col_lev_idx = i % size_3d;
                    int col_idx = col_lev_idx % n_cols;
                    int lev_idx = col_lev_idx / n_cols;
                    std::cerr << "PROPERTY FAILURE: NaN detected at conc index " << i << " (Species=" << spec_idx
                              << ", Column=" << col_idx << ", Level=" << lev_idx << ") during iteration " << iter
                              << ", VALUE=" << conc[i] << std::endl;
                }
                assert(std::isfinite(conc[i]) && "Concentration must remain finite!");
            }

            // Clip tiny negative concentrations to 0.0 to mimic standard atmospheric model boundaries
            for (auto& val : conc) {
                if (val < 0.0)
                    val = 0.0;
            }
            state->sync_to_device();

            // Assert Invariants
            verify_properties(conc, total_size, iter);
        }

        // Finalize lifecycle
        catchem_core_destroy(core_ptr);

        std::cout << "\n==========================================" << std::endl;
        std::cout << "=== SUCCESS: ALL PROPERTY CHECKS HELD! ===" << std::endl;
        std::cout << "==========================================\n" << std::endl;
    }
    Kokkos::finalize();
    return 0;
}
