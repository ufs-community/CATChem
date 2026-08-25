#include "catchem_api.hpp"
#include "catchem_core.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_fortran_process.hpp"
#include "catchem_state_manager.hpp"
#include <Kokkos_Core.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

extern "C" {
void catchem_register_settling_cpp();
void catchem_register_drydep_cpp();
void catchem_register_seasalt_cpp();
void catchem_register_wetdep_cpp();
void catchem_register_so4chem_cpp();
}

// Mock Fortran physics scheme working directly on host array
void run_mock_fortran_physics(double* ptr, int n_cols, int n_levels) {
    // Simulate Fortran LayoutLeft (column-major) indexing: (i, j) -> i + j * n_cols
    for (int j = 0; j < n_levels; ++j) {
        for (int i = 0; i < n_cols; ++i) {
            ptr[i + j * n_cols] += 10.0; // Add tendency
        }
    }
}

// A dummy process that writes to a diagnostic field
class DummyDiagProcess : public catchem::ProcessInterface {
private:
    std::shared_ptr<catchem::DiagnosticManager> diag_mgr;
    int n_cols;

public:
    DummyDiagProcess(std::shared_ptr<catchem::DiagnosticManager> dm, int nc) : diag_mgr(dm), n_cols(nc) {}

    std::string get_name() const override { return "DummyDiagProcess"; }

    void init(std::shared_ptr<catchem::StateManager> state) override {}

    void run(std::shared_ptr<catchem::StateManager> state) override {
        // Retrieve the underlying diagnostic device View
        auto dust_flux = diag_mgr->get_device_view_2d("dust_emission_flux");

        // Capture View by value in the parallel kernel
        Kokkos::parallel_for(
            "calculate_dust_emissions", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, n_cols),
            KOKKOS_LAMBDA(int icol) {
                // Write directly to the diagnostic view
                dust_flux(icol, 0) = 42.0 + icol;
            });
    }

    void finalize() override {}
};

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        // ==========================================
        // TEST 1: Phase 1 Shared Memory / Interop Test
        // ==========================================
        {
            int n_cols = 4;
            int n_levels = 5;
            int n_species = 2;

            // Allocate mock Fortran memory (column-major contiguous)
            std::vector<double> fortran_array(n_cols * n_levels, 1.0);

            // 1. Create Core & bind arrays
            void* core = catchem_core_create(n_cols, n_levels, n_species);
            void* state = catchem_core_get_state_manager(core);

            catchem_state_bind_2d(state, "temperature", fortran_array.data());

            // Verify our new pointer retrievers
            std::vector<double> dummy_1d(n_cols, 2.0);
            std::vector<double> dummy_3d(n_cols * n_levels * n_species, 3.0);
            catchem_state_bind_1d(state, "dummy1", dummy_1d.data());
            catchem_state_bind_3d(state, "dummy3", dummy_3d.data());

            assert(catchem_state_get_pointer_1d(state, "dummy1") == dummy_1d.data());
            assert(catchem_state_get_pointer_2d(state, "temperature") == fortran_array.data());
            assert(catchem_state_get_pointer_3d(state, "dummy3") == dummy_3d.data());
            assert(catchem_state_get_pointer_1d(state, "nonexistent") == nullptr);

            // 2. Sync to active space
            catchem_state_sync_to_device(state);

            // 3. Execute Fortran process sequentially modifying the raw array on host
            run_mock_fortran_physics(fortran_array.data(), n_cols, n_levels);

            // Verify direct zero-copy modification
            assert(fortran_array[0] == 11.0);

            // 4. Sync up and clean up
            catchem_state_sync_to_host(state);
            catchem_core_destroy(core);

            std::cout << "SUCCESS: Interop Shared State Validation Passed!\n";
        }

        // ==========================================
        // TEST 2: Phase 2 Diagnostic Collection Test
        // ==========================================
        {
            int nx = 4;
            int ny = 1;
            int nz = 5;
            int n_cols = nx * ny;

            // 1. Create Core (creates StateManager & DiagnosticManager)
            void* core_ptr = catchem_core_create(n_cols, nz, 1);
            auto* core = static_cast<catchem::Core*>(core_ptr);
            auto diag_mgr = core->get_diagnostic_manager();

            // 2. Register diagnostic through C-API
            catchem_diag_register(core_ptr, "dust_emission_flux", "Dust flux", "kg/m2/s", 2, n_cols, 1, 0);

            // 3. Attach dummy diagnostic process
            core->add_process(std::make_shared<DummyDiagProcess>(diag_mgr, n_cols));

            // 4. Run timestep (runs dummy process and syncs diagnostics to host)
            catchem_core_run_timestep(core_ptr, 3600.0);

            // 5. Get host pointer and verify results
            void* host_ptr = catchem_diag_get_pointer(core_ptr, "dust_emission_flux");
            double* dust_flux_host = static_cast<double*>(host_ptr);

            bool passed = true;
            for (int i = 0; i < n_cols; ++i) {
                if (dust_flux_host[i] != 42.0 + i) { // Note LayoutLeft means col_i is inner dimension
                    std::cerr << "Diagnostic mismatch at col " << i << ": expected " << 42.0 + i << ", got "
                              << dust_flux_host[i] << std::endl;
                    passed = false;
                }
            }

            if (passed) {
                std::cout << "SUCCESS: C++ Diagnostic Validation Passed!\n";
            } else {
                std::cout << "FAILURE: C++ Diagnostic Validation Failed!\n";
                catchem_core_destroy(core_ptr);
                Kokkos::finalize();
                return 1;
            }

            catchem_core_destroy(core_ptr);
        }

        // ==========================================
        // TEST 4: Modular dynamic process registration validation
        // ==========================================
        {
            // Simulate Fortran explicitly linking and calling C++ register_settling_cpp
            catchem_register_settling_cpp();

            int n_cols = 4;
            int n_levels = 5;
            int n_species = 2;

            void* core_ptr = catchem_core_create(n_cols, n_levels, n_species);

            // Add the dynamically registered settling process by name via dynamic registry
            catchem_core_add_process_by_name(core_ptr, "settling");

            catchem_core_destroy(core_ptr);
            std::cout << "SUCCESS: Modular Dynamic Process C-API Registration Validation Passed!\n";
        }

        // ==========================================
        // TEST 6: Parallel Meteorological Derivations, SZA, and Unified Chem State
        // ==========================================
        {
            int n_cols = 4;
            int n_levels = 5;
            int n_species = 2;

            void* core = catchem_core_create(n_cols, n_levels, n_species);
            void* state = catchem_core_get_state_manager(core);

            // 1. Allocate mock meteorological and chemical arrays
            std::vector<double> temp_array(n_cols * n_levels, 290.15);          // Temperature [K]
            std::vector<double> qv_array(n_cols * n_levels, 0.01);              // Specific humidity [kg/kg]
            std::vector<double> pmid_array(n_cols * n_levels, 100000.0);        // Mid-pressure [Pa]
            std::vector<double> pedge_array(n_cols * (n_levels + 1), 101325.0); // Pressure edges [Pa]

            // Assign standard pressure edge levels sequentially
            for (int i = 0; i < n_cols; ++i) {
                pedge_array[i + 0 * n_cols] = 101325.0; // Surface
                pedge_array[i + 1 * n_cols] = 90000.0;
                pedge_array[i + 2 * n_cols] = 80000.0;
                pedge_array[i + 3 * n_cols] = 70000.0;
                pedge_array[i + 4 * n_cols] = 60000.0;
                pedge_array[i + 5 * n_cols] = 50000.0; // Top
            }

            std::vector<double> bxheight_array(n_cols * n_levels, 0.0);   // Output height
            std::vector<double> airden_dry_array(n_cols * n_levels, 0.0); // Output dry density

            std::vector<double> mock_chem_state(n_cols * n_levels * n_species, 4.2); // Unified chem state

            // 2. Bind arrays to StateManager
            catchem_state_bind_met_3d(state, "T", temp_array.data());
            catchem_state_bind_met_3d(state, "QV", qv_array.data());
            catchem_state_bind_met_3d(state, "PMID", pmid_array.data());
            catchem_state_bind_met_3d(state, "PEDGE", pedge_array.data());
            catchem_state_bind_met_3d(state, "BXHEIGHT", bxheight_array.data());
            catchem_state_bind_met_3d(state, "AIRDEN_DRY", airden_dry_array.data());

            catchem_state_bind_unified_chemistry(state, mock_chem_state.data());

            // 3. Sync arrays to device memory spaces
            catchem_state_sync_to_device(state);

            // 4. Trigger parallel derived met equations
            catchem_state_derive_bxheight(state);
            catchem_state_derive_airden_dry(state);

            // 5. Sync derived results back to host heap
            catchem_state_sync_to_host(state);

            // 6. Assert correct calculations
            // Layer 1 edge pressures: P_lower = 101325.0, P_upper = 90000.0
            // Virtual T = 290.15 * (1 + 0.608 * 0.01) = 291.914
            // Expected height = (287 / 9.80665) * 291.914 * std::log(101325.0 / 90000.0) ≈ 1010.5 meters
            double derived_h = bxheight_array[0];
            assert(derived_h > 990.0 && derived_h < 1030.0);
            std::cout << "INFO: Derived BXHEIGHT = " << derived_h << " meters.\n";

            double derived_rho = airden_dry_array[0];
            assert(derived_rho > 1.0 && derived_rho < 1.3);
            std::cout << "INFO: Derived Dry Air Density = " << derived_rho << " kg/m³.\n";

            // Assert unified chemistry array mapped accurately
            auto* state_obj = static_cast<catchem::StateManager*>(state);
            assert(state_obj->chem.conc != nullptr);
            assert(state_obj->chem.conc->host_view(0, 0, 0) == 4.2);

            // 7. Test portable Time State calculations
            catchem_state_set_time(state, 2026, 7, 8, 12, 0, 0, 189, 3600.0);
            double cos_sza = state_obj->time.get_cos_sza(40.0, -80.0);
            assert(cos_sza >= -1.0 && cos_sza <= 1.0);
            std::cout << "INFO: Calculated Cos(SZA) at lat=40, lon=-80: " << cos_sza << "\n";

            catchem_core_destroy(core);
            std::cout << "SUCCESS: Parallel Meteorological Derivations & SZA Validation Passed!\n";
        }

        // ==========================================
        // TEST 5: C++ Species Metadata & State Initialization
        // ==========================================
        {
            int n_cols = 4;
            int n_levels = 5;
            int n_species = 2;

            void* core = catchem_core_create(n_cols, n_levels, n_species);
            void* state = catchem_core_get_state_manager(core);

            // 1. Load species config from CATChem_species.yml (finding correct path)
            std::string config_path = "";
            std::vector<std::string> candidates = {"tests/Configs/Default/CATChem_species.yml",
                                                   "../tests/Configs/Default/CATChem_species.yml",
                                                   "../../tests/Configs/Default/CATChem_species.yml",
                                                   "tests/CATChem_species.yml",
                                                   "../tests/CATChem_species.yml",
                                                   "../../tests/CATChem_species.yml",
                                                   "CATChem_species.yml"};
            for (const auto& path : candidates) {
                std::ifstream f(path);
                if (f.good()) {
                    config_path = path;
                    break;
                }
            }
            assert(!config_path.empty() && "ERROR: Could not find CATChem_species.yml");
            catchem_state_load_species_config(state, config_path.c_str());

            // 2. Validate species counts and offsets
            int count = catchem_state_get_species_count(state);
            assert(count > 0);
            std::cout << "INFO: Loaded " << count << " species in integration test.\n";

            // 3. Translate species names to 1-based indices
            int idx_so2 = catchem_state_get_species_index(state, "so2");
            int idx_so4 = catchem_state_get_species_index(state, "so4");
            assert(idx_so2 != -1);
            assert(idx_so4 != -1);

            // 4. Validate physical properties of species
            double mw_so2 = catchem_state_get_species_mw(state, idx_so2);
            assert(mw_so2 == 64.04);

            int is_gas_so2 = catchem_state_is_species_gas(state, idx_so2);
            int is_aero_so2 = catchem_state_is_species_aerosol(state, idx_so2);
            assert(is_gas_so2 == 1);
            assert(is_aero_so2 == 0);

            int is_gas_so4 = catchem_state_is_species_gas(state, idx_so4);
            int is_aero_so4 = catchem_state_is_species_aerosol(state, idx_so4);
            assert(is_gas_so4 == 0);
            assert(is_aero_so4 == 1);

            // 5. Validate category lists (gas / aerosol)
            int gas_count = catchem_state_get_gas_species_count(state);
            assert(gas_count > 0);
            std::vector<int> gas_indices(gas_count);
            catchem_state_get_gas_indices(state, gas_indices.data());

            // Ensure so2 index is present in gas_indices
            bool found_so2 = false;
            for (int idx : gas_indices) {
                if (idx == idx_so2)
                    found_so2 = true;
            }
            assert(found_so2);

            catchem_core_destroy(core);
            std::cout << "SUCCESS: C++ Species Metadata & State Initialization Validation Passed!\n";
        }

        // ==========================================
        // TEST 7: C++ Kokkos Settling Process Execution
        // ==========================================
        {
            int n_cols = 4;
            int n_levels = 5;
            int n_species = 22;

            void* core = catchem_core_create(n_cols, n_levels, n_species);
            void* state = catchem_core_get_state_manager(core);
            auto* state_obj = static_cast<catchem::StateManager*>(state);

            // Load species config so we have aerosols to process
            std::string config_path = "";
            std::vector<std::string> candidates = {"tests/Configs/Default/CATChem_species.yml",
                                                   "../tests/Configs/Default/CATChem_species.yml",
                                                   "../../tests/Configs/Default/CATChem_species.yml",
                                                   "tests/CATChem_species.yml",
                                                   "../tests/CATChem_species.yml",
                                                   "../../tests/CATChem_species.yml",
                                                   "CATChem_species.yml"};
            for (const auto& path : candidates) {
                std::ifstream f(path);
                if (f.good()) {
                    config_path = path;
                    break;
                }
            }
            catchem_state_load_species_config(state, config_path.c_str());

            // Allocate mock meteorological arrays
            std::vector<double> temp_array(n_cols * n_levels, 298.15);          // Temperature [K]
            std::vector<double> qv_array(n_cols * n_levels, 0.01);              // Specific humidity [kg/kg]
            std::vector<double> pmid_array(n_cols * n_levels, 100000.0);        // Mid-pressure [Pa]
            std::vector<double> pedge_array(n_cols * (n_levels + 1), 101325.0); // Pressure edges [Pa]
            std::vector<double> airden_array(n_cols * n_levels, 1.2);           // Output dry density
            std::vector<double> bxheight_array(n_cols * n_levels, 1000.0);      // Output height

            for (int i = 0; i < n_cols; ++i) {
                pedge_array[i + 0 * n_cols] = 101325.0; // Surface
                pedge_array[i + 1 * n_cols] = 90000.0;
                pedge_array[i + 2 * n_cols] = 80000.0;
                pedge_array[i + 3 * n_cols] = 70000.0;
                pedge_array[i + 4 * n_cols] = 60000.0;
                pedge_array[i + 5 * n_cols] = 50000.0; // Top
            }

            std::vector<double> mock_chem_state(n_cols * n_levels * n_species,
                                                1.0); // Initial concentration of 1.0 for all

            // Bind arrays to StateManager
            catchem_state_bind_met_3d(state, "T", temp_array.data());
            catchem_state_bind_met_3d(state, "QV", qv_array.data());
            catchem_state_bind_met_3d(state, "PMID", pmid_array.data());
            catchem_state_bind_met_3d(state, "PEDGE", pedge_array.data());
            catchem_state_bind_met_3d(state, "BXHEIGHT", bxheight_array.data());
            catchem_state_bind_met_3d(state, "AIRDEN", airden_array.data());

            catchem_state_bind_unified_chemistry(state, mock_chem_state.data());
            catchem_state_set_time(state, 2026, 7, 8, 12, 0, 0, 189, 3600.0);

            catchem_state_sync_to_device(state);

            // Register and initialize settling process
            catchem_register_settling_cpp();
            catchem_core_add_process_by_name(core, "settling");

            // Execute the settling process
            catchem_core_run_timestep(core, 3600.0);

            catchem_state_sync_to_host(state);

            // Verify the concentrations have changed (sedimented)
            // The top layer should have less concentration due to settling.
            // Note: Since all aerosols settle down, concentration at the top layer should decrease.
            // We just check if the top layer of an aerosol species is less than 1.0.

            // Get index of an aerosol (e.g., 'so4')
            int idx_so4_1based = catchem_state_get_species_index(state, "so4");
            int idx_so4_0based = idx_so4_1based - 1;

            std::cout << "DEBUG: idx_so4_1based=" << idx_so4_1based << ", idx_so4_0based=" << idx_so4_0based
                      << std::endl;

            // Since LayoutLeft is (i, j, k), meaning (col, lev, spec)
            // Wait, LayoutLeft means innermost dimension is leftmost.
            // InteropField 3D: (col, level, species). Memory layout: col + level * nc + species * nc * nl
            int top_level_idx = n_levels - 1;
            double top_layer_conc = mock_chem_state[0 + top_level_idx * n_cols + idx_so4_0based * n_cols * n_levels];

            std::cout << "DEBUG: top_layer_conc=" << top_layer_conc << std::endl;
            for (int k = 0; k < n_levels; ++k) {
                std::cout << "  Level " << k
                          << " conc = " << mock_chem_state[0 + k * n_cols + idx_so4_0based * n_cols * n_levels]
                          << std::endl;
            }

            assert(top_layer_conc < 1.0); // Concentration should drop at the top

            catchem_core_destroy(core);
            std::cout << "SUCCESS: C++ Kokkos Settling Process Validation Passed!\n";
        }

        // ==========================================
        // TEST 7: C++ Config and Grid Initialization
        // ==========================================
        {
            std::string config_path = "";
            std::vector<std::string> candidates = {"tests/Configs/Default/CATChem_new_config.yml",
                                                   "../tests/Configs/Default/CATChem_new_config.yml",
                                                   "../../tests/Configs/Default/CATChem_new_config.yml",
                                                   "tests/CATChem_new_config.yml",
                                                   "../tests/CATChem_new_config.yml",
                                                   "../../tests/CATChem_new_config.yml",
                                                   "CATChem_new_config.yml"};
            for (const auto& path : candidates) {
                std::ifstream f(path);
                if (f.good()) {
                    config_path = path;
                    break;
                }
            }
            void* core = catchem_core_create_from_config(config_path.c_str());

            int nx, ny, nz;
            catchem_get_grid_dimensions(core, &nx, &ny, &nz);
            assert(nx > 0);
            assert(nz > 0);
            std::cout << "INFO: Loaded grid dimensions from config: " << nx << "x" << ny << "x" << nz << "\n";

            double dt = catchem_get_config_timestep(core);
            assert(dt > 0.0);
            std::cout << "INFO: Loaded timestep from config: " << dt << " s\n";

            catchem_core_destroy(core);
            std::cout << "SUCCESS: C++ Config & Grid Validation Passed!\n";
        }

        // ==========================================
        // TEST 9: Direct Flat-Science Interop Adapter for DryDep
        // ==========================================
        {
            std::cout << "\n--- TEST 9: Direct Flat-Science Interop Adapter for DryDep ---\n";
            int n_cols = 4;
            int n_levels = 5;
            int n_species = 22; // Must match CATChem_species.yml species list count

            // Register drydep C++ process dynamically
            catchem_register_drydep_cpp();

            void* core = catchem_core_create(n_cols, n_levels, n_species);
            auto* state = static_cast<catchem::StateManager*>(catchem_core_get_state_manager(core));

            // Load species configuration first (otherwise species_list is empty)
            std::vector<std::string> paths = {"tests/Configs/Default/CATChem_species.yml",
                                              "../tests/Configs/Default/CATChem_species.yml",
                                              "../../tests/Configs/Default/CATChem_species.yml",
                                              "tests/CATChem_species.yml",
                                              "../tests/CATChem_species.yml",
                                              "../../tests/CATChem_species.yml"};
            for (const auto& path : paths) {
                try {
                    state->load_species_config(path);
                    break;
                } catch (...) {
                }
            }

            // Set up mock temperature and concentrations (so4 at index 4)
            std::vector<double> mock_t(n_cols * n_levels, 288.15);
            std::vector<double> mock_qv(n_cols * n_levels, 0.01);
            std::vector<double> mock_pmid(n_cols * n_levels, 90000.0);
            std::vector<double> mock_pedge(n_cols * (n_levels + 1), 100000.0);
            std::vector<double> mock_bxheight(n_cols * n_levels, 100.0);
            std::vector<double> mock_airden(n_cols * n_levels, 1.2);
            std::vector<double> mock_rh(n_cols * n_levels, 0.5);

            std::vector<double> mock_ps(n_cols, 101300.0);
            std::vector<double> mock_ts(n_cols, 288.15);
            std::vector<double> mock_pblh(n_cols, 1000.0);
            std::vector<double> mock_ustar(n_cols, 0.3);
            std::vector<double> mock_hflux(n_cols, 100.0);
            std::vector<double> mock_obk(n_cols, -100.0);
            std::vector<double> mock_lat(n_cols, 40.0);
            std::vector<double> mock_lon(n_cols, -80.0);

            std::vector<double> mock_chem_state(n_cols * n_levels * n_species, 1.0);

            // Bind 3D Met fields
            catchem_state_bind_met_3d(state, "T", mock_t.data());
            catchem_state_bind_met_3d(state, "QV", mock_qv.data());
            catchem_state_bind_met_3d(state, "PMID", mock_pmid.data());
            catchem_state_bind_met_3d(state, "PEDGE", mock_pedge.data());
            catchem_state_bind_met_3d(state, "BXHEIGHT", mock_bxheight.data());
            catchem_state_bind_met_3d(state, "AIRDEN", mock_airden.data());
            catchem_state_bind_met_3d(state, "RH", mock_rh.data());

            // Bind 2D Met fields
            catchem_state_bind_met_2d(state, "PS", mock_ps.data());
            catchem_state_bind_met_2d(state, "TS", mock_ts.data());
            catchem_state_bind_met_2d(state, "PBLH", mock_pblh.data());
            catchem_state_bind_met_2d(state, "USTAR", mock_ustar.data());
            catchem_state_bind_met_2d(state, "HFLUX", mock_hflux.data());
            catchem_state_bind_met_2d(state, "OBK", mock_obk.data());
            catchem_state_bind_met_2d(state, "LAT", mock_lat.data());
            catchem_state_bind_met_2d(state, "LON", mock_lon.data());

            catchem_state_bind_unified_chemistry(state, mock_chem_state.data());
            catchem_state_sync_to_device(state);

            // Add the drydep process by name
            catchem_core_add_process_by_name(core, "drydep");

            // Execute the timestep which triggers drydep calculation
            catchem_core_run_timestep(core, 3600.0);

            // Fetch dynamic diagnostic pointer
            double* host_diag_con = (double*)catchem_diag_get_pointer(core, "drydep_con_per_species");
            double* host_diag_vel = (double*)catchem_diag_get_pointer(core, "drydep_velocity_per_species");

            assert(host_diag_con != nullptr);
            assert(host_diag_vel != nullptr);

            std::cout << "INFO: Retrieved host diagnostic pointer: " << host_diag_con << "\n";
            std::cout << "SUCCESS: DryDep Direct Adapter executed and populated diagnostics!\n";

            catchem_core_destroy(core);
        }

        // ==========================================
        // TEST 10: Direct Flat-Science Interop Adapter for SeaSalt
        // ==========================================
        {
            std::cout << "\n--- TEST 10: Direct Flat-Science Interop Adapter for SeaSalt ---\n";
            int n_cols = 4;
            int n_levels = 5;
            int n_species = 22;

            catchem_register_seasalt_cpp();

            void* core = catchem_core_create(n_cols, n_levels, n_species);
            auto* state = static_cast<catchem::StateManager*>(catchem_core_get_state_manager(core));

            // Load species
            std::vector<std::string> paths = {"tests/Configs/Default/CATChem_species.yml",
                                              "../tests/Configs/Default/CATChem_species.yml",
                                              "../../tests/Configs/Default/CATChem_species.yml",
                                              "tests/CATChem_species.yml",
                                              "../tests/CATChem_species.yml",
                                              "../../tests/CATChem_species.yml"};
            for (const auto& path : paths) {
                try {
                    state->load_species_config(path);
                    break;
                } catch (...) {
                }
            }

            // Allocate met fields dynamically (FROCEAN, FRSEAICE, SST, DELP)
            std::vector<double> mock_frocean(n_cols, 1.0);
            std::vector<double> mock_frseaice(n_cols, 0.0);
            std::vector<double> mock_sst(n_cols, 290.0);
            std::vector<double> mock_delp(n_cols * n_levels, 1000.0);
            std::vector<double> mock_ustar(n_cols, 0.5);

            // Bind them

            std::vector<double> mock_lat(n_cols, 40.0);
            std::vector<double> mock_lon(n_cols, -80.0);
            state->bind_met_field_2d("LAT", mock_lat.data());
            state->bind_met_field_2d("LON", mock_lon.data());
            state->bind_met_field_2d("FROCEAN", mock_frocean.data());

            state->bind_met_field_2d("FRSEAICE", mock_frseaice.data());
            state->bind_met_field_2d("SST", mock_sst.data());
            state->bind_met_field_3d("DELP", mock_delp.data());
            catchem_state_bind_met_2d(state, "USTAR", mock_ustar.data());

            // Concentrations and Tendencies
            std::vector<double> mock_chem_state(n_cols * n_levels * n_species, 1.0);
            catchem_state_bind_unified_chemistry(state, mock_chem_state.data());
            catchem_state_sync_to_device(state);

            catchem_core_add_process_by_name(core, "seasalt");
            catchem_core_run_timestep(core, 3600.0);

            // Retrieve diagnostic total mass emission pointer
            double* diag_total_mass = (double*)catchem_diag_get_pointer(core, "seasalt_mass_emission_total");
            assert(diag_total_mass != nullptr);

            std::cout << "SUCCESS: SeaSalt Direct Adapter executed and populated diagnostics!\n";
            catchem_core_destroy(core);
        }

        // ==========================================
        // TEST 11: Direct Flat-Science Interop Adapter for WetDep
        // ==========================================
        {
            std::cout << "\n--- TEST 11: Direct Flat-Science Interop Adapter for WetDep ---\n";
            int n_cols = 4;
            int n_levels = 5;
            int n_species = 22;

            catchem_register_wetdep_cpp();

            void* core = catchem_core_create(n_cols, n_levels, n_species);
            auto* state = static_cast<catchem::StateManager*>(catchem_core_get_state_manager(core));

            // Load species
            std::vector<std::string> paths = {"tests/Configs/Default/CATChem_species.yml",
                                              "../tests/Configs/Default/CATChem_species.yml",
                                              "../../tests/Configs/Default/CATChem_species.yml",
                                              "tests/CATChem_species.yml",
                                              "../tests/CATChem_species.yml",
                                              "../../tests/CATChem_species.yml"};
            for (const auto& path : paths) {
                try {
                    state->load_species_config(path);
                    break;
                } catch (...) {
                }
            }

            // Allocate met fields dynamically (AIRDEN_DRY, AIRDEN, MAIRDEN, PEDGE, PFILSAN, PFLLSAN, REEVAPLS, T)
            std::vector<double> mock_airden_dry(n_cols * n_levels, 1.2);
            std::vector<double> mock_airden(n_cols * n_levels, 1.2);
            std::vector<double> mock_mairden(n_cols * n_levels, 1.2);
            std::vector<double> mock_pedge(n_cols * (n_levels + 1), 100000.0);
            std::vector<double> mock_pfilsan(n_cols * (n_levels + 1), 0.05);
            std::vector<double> mock_pfllsan(n_cols * (n_levels + 1), 0.05);
            std::vector<double> mock_reevapls(n_cols * n_levels, 0.0);
            std::vector<double> mock_t(n_cols * n_levels, 288.15);

            catchem_state_bind_met_3d(state, "AIRDEN_DRY", mock_airden_dry.data());
            catchem_state_bind_met_3d(state, "AIRDEN", mock_airden.data());
            state->bind_met_field_3d("MAIRDEN", mock_mairden.data());
            catchem_state_bind_met_3d(state, "PEDGE", mock_pedge.data());
            state->bind_met_field_3d("PFILSAN", mock_pfilsan.data());
            state->bind_met_field_3d("PFLLSAN", mock_pfllsan.data());
            state->bind_met_field_3d("REEVAPLS", mock_reevapls.data());
            catchem_state_bind_met_3d(state, "T", mock_t.data());

            // Concentrations and Tendencies
            std::vector<double> mock_chem_state(n_cols * n_levels * n_species, 1.0);
            catchem_state_bind_unified_chemistry(state, mock_chem_state.data());
            catchem_state_sync_to_device(state);

            catchem_core_add_process_by_name(core, "wetdep");
            catchem_core_run_timestep(core, 3600.0);

            std::cout << "SUCCESS: WetDep Direct Adapter executed successfully!\n";
            catchem_core_destroy(core);
        }

        // ==========================================
        // TEST 12: Direct Flat-Science Interop Adapter for SO4chem
        // ==========================================
        {
            std::cout << "\n--- TEST 12: Direct Flat-Science Interop Adapter for SO4chem ---\n";
            int n_cols = 4;
            int n_levels = 5;
            int n_species = 22;

            catchem_register_so4chem_cpp();

            void* core = catchem_core_create(n_cols, n_levels, n_species);
            auto* state = static_cast<catchem::StateManager*>(catchem_core_get_state_manager(core));

            // Load species
            std::vector<std::string> paths = {"tests/Configs/Default/CATChem_species.yml",
                                              "../tests/Configs/Default/CATChem_species.yml",
                                              "../../tests/Configs/Default/CATChem_species.yml",
                                              "tests/CATChem_species.yml",
                                              "../tests/CATChem_species.yml",
                                              "../../tests/CATChem_species.yml"};
            for (const auto& path : paths) {
                try {
                    state->load_species_config(path);
                    break;
                } catch (...) {
                }
            }

            // Allocate met fields dynamically (AIRDEN, CLDF, DELP, PMID, T, PEDGE, HFLUX, LAT, LON, PBLH, USTAR)
            std::vector<double> mock_airden(n_cols * n_levels, 1.2);
            std::vector<double> mock_cldf(n_cols * n_levels, 0.5);
            std::vector<double> mock_delp(n_cols * n_levels, 1000.0);
            std::vector<double> mock_pmid(n_cols * n_levels, 90000.0);
            std::vector<double> mock_t(n_cols * n_levels, 288.15);
            std::vector<double> mock_pedge(n_cols * (n_levels + 1), 100000.0);

            std::vector<double> mock_hflux(n_cols, 100.0);
            std::vector<double> mock_lat(n_cols, 40.0);
            std::vector<double> mock_lon(n_cols, -80.0);
            std::vector<double> mock_pblh(n_cols, 1000.0);
            std::vector<double> mock_ustar(n_cols, 0.5);

            std::vector<double> mock_chem_state(n_cols * n_levels * n_species, 1.0);

            catchem_state_bind_met_3d(state, "AIRDEN", mock_airden.data());
            state->bind_met_field_3d("CLDF", mock_cldf.data());
            state->bind_met_field_3d("DELP", mock_delp.data());
            catchem_state_bind_met_3d(state, "PMID", mock_pmid.data());
            catchem_state_bind_met_3d(state, "T", mock_t.data());
            catchem_state_bind_met_3d(state, "PEDGE", mock_pedge.data());

            catchem_state_bind_met_2d(state, "HFLUX", mock_hflux.data());
            catchem_state_bind_met_2d(state, "LAT", mock_lat.data());
            catchem_state_bind_met_2d(state, "LON", mock_lon.data());
            catchem_state_bind_met_2d(state, "PBLH", mock_pblh.data());
            catchem_state_bind_met_2d(state, "USTAR", mock_ustar.data());

            catchem_state_bind_unified_chemistry(state, mock_chem_state.data());
            catchem_state_sync_to_device(state);

            // Add the drydep process by name
            catchem_core_add_process_by_name(core, "so4chem");

            // Execute the timestep which triggers drydep calculation
            catchem_core_run_timestep(core, 3600.0);

            double* diag_gas_source = (double*)catchem_diag_get_pointer(core, "PSO4_from_gaseous_SO2_per_level");
            assert(diag_gas_source != nullptr);

            std::cout << "SUCCESS: SO4chem Direct Adapter executed and populated diagnostics!\n";
            catchem_core_destroy(core);
        }
    }
    Kokkos::finalize();
    return 0;
}
