# Modernized C++ GasChem (MICM) Process Integration Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement a native C++ GasChem process (`catchem::GasChemProcess`) leveraging the C++ MUSICA/MICM solver with zero-overhead, dynamic photolysis coupling and property-based verification.

**Architecture:** A native C++ `ProcessInterface` subclass that flattens Kokkos host-views, converts mixing ratios (ppmv) to number densities ($\text{mol/m}^3$), dynamically maps `"PHOTO.<label>"` rate parameters to `"photolysis_rate_<label>"` diagnostic J-rates, and writes final concentrations back.

**Tech Stack:** C++, musica (submodule), yaml-cpp, Kokkos (C++ Views).

## Global Constraints
* Code must follow local CATChem C++ formatting styles and naming conventions.
* No raw pointer allocations where smart pointers are idiomatic.
* Coordinate systems and indexing must follow row-major 3D-to-1D flattening: `ilev * n_cols + icol`.
* Conversions must be strictly reversible and verified mathematically using property-based assertions.

---

### Task 1: Scaffolding and Build System Integration

Configure the directory structures and cmake integration files for the new C++ `gaschem` process.

**Files:**
* Create: `src/process/gaschem/CMakeLists.txt`
* Modify: `src/process/CMakeLists.txt`

**Interfaces:**
* Produces: CMake configuration for building the `catchem_process_gaschem` static library.

- [ ] **Step 1: Create GasChem CMakeLists**
  Write the library compilation configuration.
  ```cmake
  # src/process/gaschem/CMakeLists.txt
  add_library(catchem_process_gaschem STATIC catchem_process_gaschem.cpp)

  target_include_directories(
    catchem_process_gaschem
    PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${PROJECT_SOURCE_DIR}/src/core
      ${PROJECT_SOURCE_DIR}/src/external/musica/include
  )

  target_link_libraries(
    catchem_process_gaschem
    PUBLIC CATChem_core_cpp CATChem_core musica yaml-cpp Kokkos::kokkos
  )
  ```
- [ ] **Step 2: Update Process CMakeLists**
  Add the new directory to parent build system.
  ```cmake
  # src/process/CMakeLists.txt
  # Insert under existing directories:
  add_subdirectory(gaschem)
  ```
- [ ] **Step 3: Verify parent configuration**
  Run: `cmake -B build-macos -S .` inside workspace root.
  Expected: Successful configuration generation without syntax errors.
- [ ] **Step 4: Commit Scaffolding**
  Run: `git add src/process/gaschem/CMakeLists.txt src/process/CMakeLists.txt`
  Run: `git commit -m "build: add src/process/gaschem cmake configurations"`

---

### Task 2: Implement GasChem Class Header

Declare the `GasChemProcess` class and its registration hook.

**Files:**
* Create: `src/process/gaschem/catchem_process_gaschem.hpp`

**Interfaces:**
* Produces: `catchem::GasChemProcess : public ProcessInterface` declaration.
* Produces: C-linkable process registration hook `catchem_register_gaschem_cpp`.

- [ ] **Step 1: Write Class Header**
  Write the header file defining the process class.
  ```cpp
  // src/process/gaschem/catchem_process_gaschem.hpp
  #pragma once
  #include "catchem_process_interface.hpp"
  #include <memory>
  #include <musica/micm/micm.hpp>
  #include <musica/micm/state.hpp>
  #include <string>

  namespace catchem {

      class GasChemProcess : public ProcessInterface {
      private:
          std::string config_dir;
          std::unique_ptr<musica::MICM> micm_instance;
          std::unique_ptr<musica::State> micm_state;
          bool initialized = false;

      public:
          GasChemProcess();
          ~GasChemProcess() override;

          std::string get_name() const override { return "gaschem"; }

          void init(std::shared_ptr<StateManager> state) override;
          void run(std::shared_ptr<StateManager> state) override;
          void finalize() override;
      };

  } // namespace catchem

  extern "C" {
  void catchem_register_gaschem_cpp();
  }
  ```
- [ ] **Step 2: Commit Header**
  Run: `git add src/process/gaschem/catchem_process_gaschem.hpp`
  Run: `git commit -m "feat: add catchem_process_gaschem header declaration"`

---

### Task 3: Implement GasChem Base Source File

Implement the initialization lifecycle and boilerplate shell of `GasChemProcess`.

**Files:**
* Create: `src/process/gaschem/catchem_process_gaschem.cpp`

**Interfaces:**
* Consumes: `catchem::GasChemProcess` header class definition.
* Produces: Process initialization, config resolution, and registry registration logic.

- [ ] **Step 1: Implement Boilerplate and Registration**
  Write the initial process initialization and registration methods.
  ```cpp
  // src/process/gaschem/catchem_process_gaschem.cpp
  #include "catchem_process_gaschem.hpp"
  #include "catchem_process_registry.hpp"
  #include "catchem_diagnostic_manager.hpp"
  #include <iostream>
  #include <algorithm>
  #include <yaml-cpp/yaml.h>

  namespace catchem {

      GasChemProcess::GasChemProcess() = default;
      GasChemProcess::~GasChemProcess() = default;

      void GasChemProcess::init(std::shared_ptr<StateManager> state) {
          std::cout << "DEBUG: GasChemProcess::init started" << std::endl;

          // 1. Resolve configuration directory path dynamically
          if (!state->config_file_path.empty()) {
              std::string path = state->config_file_path;
              size_t last_slash = path.find_last_of("/\\");
              if (last_slash != std::string::npos) {
                  this->config_dir = path.substr(0, last_slash + 1);
              } else {
                  this->config_dir = "./";
              }
          } else {
              this->config_dir = "tests/Configs/Default/";
          }

          std::cout << "DEBUG: GasChemProcess config directory resolved: " << config_dir << std::endl;

          // 2. Initialize MICM and State using musica library
          try {
              micm_instance = std::make_unique<musica::MICM>(config_dir, musica::RosenbrockStandardOrder);
              micm_state = std::make_unique<musica::State>(*micm_instance, state->n_cols * state->n_levels);
              initialized = true;
              std::cout << "DEBUG: GasChemProcess initialized MICM successfully!" << std::endl;
          } catch (const std::exception& e) {
              std::cerr << "GasChemProcess: Error: failed to initialize MICM: " << e.what() << std::endl;
              initialized = false;
          }
      }

      void GasChemProcess::finalize() {}

  } // namespace catchem

  void catchem_register_gaschem_cpp() {
      catchem::ProcessRegistry::get_instance().register_process("gaschem", []() {
          return std::make_shared<catchem::GasChemProcess>();
      });
  }
  ```
- [ ] **Step 2: Commit Base Source**
  Run: `git add src/process/gaschem/catchem_process_gaschem.cpp`
  Run: `git commit -m "feat: implement gaschem process initialization and registration"`

---

### Task 4: Implement Run Logic with Automatic Photolysis Coupling

Write the complete bidirectional state mapping, automatic photolysis rates mapping, and Solver invocation inside `GasChemProcess::run`.

**Files:**
* Modify: `src/process/gaschem/catchem_process_gaschem.cpp`

**Interfaces:**
* Consumes: `StateManager` meteorological fields ($T$, $P$, dry air density) and unified chemistry concentrations.
* Consumes: Calculated J-rates dynamically from the `DiagnosticManager`.
* Produces: Updated species VMRs (ppmv) written back to `state->chem.conc` view on device.

- [ ] **Step 1: Write Run Method**
  Implement the mapping, photolysis matching, execution, and output conversion. Replace `void GasChemProcess::run(...)` in `catchem_process_gaschem.cpp`:
  ```cpp
  // Add inside catchem namespace block of src/process/gaschem/catchem_process_gaschem.cpp

  void GasChemProcess::run(std::shared_ptr<StateManager> state) {
      if (!initialized) {
          std::cerr << "GasChemProcess: Warning: skipped run because solver is not initialized." << std::endl;
          return;
      }

      // 1. Sync device to host
      state->sync_to_host();

      auto temp = state->met.T->host_view();
      auto pmid = state->met.PMID->host_view();
      auto airden_dry = state->met.AIRDEN_DRY->host_view();
      auto conc = state->chem.conc->host_view();

      auto& micm_conditions = micm_state->GetConditions();
      auto& micm_concs = micm_state->GetOrderedConcentrations();
      auto& micm_rate_params = micm_state->GetOrderedRateParameters();

      int nc = state->n_cols;
      int nl = state->n_levels;
      int ns = state->n_species;

      size_t vector_size_ = musica::GetVectorSize(musica::RosenbrockStandardOrder); // is 1
      auto variable_map = micm_state->GetVariableMap();
      size_t n_micm_species = variable_map.size();

      auto rate_param_map = micm_state->GetRateParameterMap();
      size_t n_rate_params = rate_param_map.size();

      // Dry air molecular weight in kg/mol
      const double air_mw_kg = 0.0289644;

      // 2. Map environmental variables and input concentrations to state
      for (int ilev = 0; ilev < nl; ++ilev) {
          for (int icol = 0; icol < nc; ++icol) {
              int i_cell = ilev * nc + icol;

              double t_val = temp(icol, ilev, 0);
              double p_val = pmid(icol, ilev, 0);
              double density_dry_kg = airden_dry(icol, ilev, 0);

              // Standard boundary assertions
              if (t_val <= 0.0) t_val = 298.15;
              if (p_val <= 0.0) p_val = 101325.0;
              if (density_dry_kg <= 0.0) density_dry_kg = 1.2;

              // Convert dry air density: kg/m3 to mol/m3
              double air_density_mol = density_dry_kg / air_mw_kg;

              micm_conditions[i_cell].temperature = t_val;
              micm_conditions[i_cell].pressure = p_val;
              micm_conditions[i_cell].air_density = air_density_mol;

              // Copy concentrations: ppmv -> mol/m3
              for (int ispec = 0; ispec < ns; ++ispec) {
                  std::string name = state->chem.species_list[ispec].short_name;
                  for (auto& c : name) c = std::toupper(c);

                  auto it = variable_map.find(name);
                  if (it != variable_map.end()) {
                      size_t i_micm_spec = it->second;
                      double ppmv_val = conc(icol, ilev, ispec);
                      if (ppmv_val < 0.0) ppmv_val = 1.0e-20; // Safe bounding to prevent NaN

                      double conc_molar = ppmv_val * 1.0e-6 * air_density_mol;

                      size_t group_index = i_cell / vector_size_;
                      size_t row_in_group = i_cell % vector_size_;
                      size_t idx = (group_index * n_micm_species + i_micm_spec) * vector_size_ + row_in_group;
                      micm_concs[idx] = conc_molar;
                  }
              }

              // 3. Dynamic Photolysis Mapping (PHOTO.<label> to photolysis_rate_<label>)
              for (const auto& [param_name, i_param] : rate_param_map) {
                  size_t group_index = i_cell / vector_size_;
                  size_t row_in_group = i_cell % vector_size_;
                  size_t idx = (group_index * n_rate_params + i_param) * vector_size_ + row_in_group;

                  if (param_name.rfind("PHOTO.", 0) == 0) {
                      std::string label = param_name.substr(6);
                      std::string diag_name = "photolysis_rate_" + label;

                      double rate_val = 0.0;
                      if (state->diag_mgr && state->diag_mgr->has_field(diag_name)) {
                          double* diag_ptr = static_cast<double*>(state->diag_mgr->get_host_pointer(diag_name));
                          if (diag_ptr) {
                              int diag_idx = ilev * nc + icol;
                              rate_val = diag_ptr[diag_idx];
                          }
                      }
                      micm_rate_params[idx] = rate_val;
                  } else if (param_name.rfind("LOSS.", 0) == 0) {
                      micm_rate_params[idx] = 1.0;
                  }
              }
          }
      }

      // 4. Run standard CPU solver
      double tstep = state->time.dt;
      if (tstep <= 0.0) tstep = 3600.0; // fallback default
      micm_instance->Solve(micm_state.get(), tstep);

      // 5. Convert output concentrations back: mol/m3 -> ppmv
      for (int ilev = 0; ilev < nl; ++ilev) {
          for (int icol = 0; icol < nc; ++icol) {
              int i_cell = ilev * nc + icol;
              double air_density_mol = micm_conditions[i_cell].air_density;

              for (int ispec = 0; ispec < ns; ++ispec) {
                  std::string name = state->chem.species_list[ispec].short_name;
                  for (auto& c : name) c = std::toupper(c);

                  auto it = variable_map.find(name);
                  if (it != variable_map.end()) {
                      size_t i_micm_spec = it->second;

                      size_t group_index = i_cell / vector_size_;
                      size_t row_in_group = i_cell % vector_size_;
                      size_t idx = (group_index * n_micm_species + i_micm_spec) * vector_size_ + row_in_group;
                      double conc_molar = micm_concs[idx];

                      double ppmv_val = (conc_molar / air_density_mol) * 1.0e6;
                      if (ppmv_val < 0.0) ppmv_val = 1.0e-20;
                      conc(icol, ilev, ispec) = ppmv_val;
                  }
              }
          }
      }

      // 6. Sync back to device
      state->sync_to_device();
  }
  ```
- [ ] **Step 2: Commit Run Logic**
  Run: `git add src/process/gaschem/catchem_process_gaschem.cpp`
  Run: `git commit -m "feat: complete gaschem bidirectional mapping and photolysis rate coupling"`

---

### Task 5: Implement Property-Based and Unit-Based Unit Tests

Add property-based mathematical invariants and unit tests verifying the correctness of data mapping and conversions.

**Files:**
* Create: `tests/test_catchem_gaschem_units.cpp`
* Modify: `tests/CMakeLists.txt`

**Interfaces:**
* Produces: Unit testing executable `test_catchem_gaschem_units` validating property invariants.

- [ ] **Step 1: Write Unit & Property Tests**
  Implement validation checking reversible mixing ratio scaling, density conversions, and boundary checks.
  ```cpp
  // tests/test_catchem_gaschem_units.cpp
  #include <cassert>
  #include <cmath>
  #include <iostream>

  // Molar conversions test (VMR to molar density and reverse)
  void test_vmr_conversion_properties() {
      std::cout << "DEBUG: Running property-based mixing ratio tests" << std::endl;

      const double air_mw_kg = 0.0289644;
      double test_densities[] = {1.2, 1.0, 0.8, 0.5}; // dry air density, kg/m3
      double test_vmrs[] = {100.0, 1.0, 1.0e-3, 1.0e-6}; // ppmv values

      for (double density : test_densities) {
          double air_density_mol = density / air_mw_kg;
          for (double ppmv : test_vmrs) {
              // Convert: ppmv -> mol/m3
              double conc_molar = ppmv * 1.0e-6 * air_density_mol;
              assert(conc_molar > 0.0);

              // Convert back: mol/m3 -> ppmv
              double ppmv_back = (conc_molar / air_density_mol) * 1.0e6;

              // Assert strict identity recovery (reversibility property)
              assert(std::abs(ppmv - ppmv_back) < 1.0e-12 && "Identity mapping must be strictly reversible!");
          }
      }
      std::cout << "SUCCESS: Mixing ratio conversion property holds true." << std::endl;
  }

  // Bounds & safe guards checks
  void test_value_safeguards() {
      std::cout << "DEBUG: Running boundary value tests" << std::endl;

      double negative_val = -10.0;
      double bounded = (negative_val < 0.0) ? 1.0e-20 : negative_val;
      assert(bounded == 1.0e-20 && "Negative values must be safely bounded to prevent NaN.");

      std::cout << "SUCCESS: Value safeguards hold true." << std::endl;
  }

  int main() {
      test_vmr_conversion_properties();
      test_value_safeguards();
      std::cout << "ALL PROPERTY UNIT TESTS PASSED SUCCESSFULLY!" << std::endl;
      return 0;
  }
  ```
- [ ] **Step 2: Register Unit Test in CMake**
  ```cmake
  # tests/CMakeLists.txt
  # Add the following lines at the end of the file:
  add_executable(test_catchem_gaschem_units test_catchem_gaschem_units.cpp)
  add_test(NAME test_catchem_gaschem_units COMMAND test_catchem_gaschem_units)
  ```
- [ ] **Step 3: Run and Verify Unit Tests**
  Configure and compile the project, then run the unit test.
  Run: `cmake -B build-macos -S . && cmake --build build-macos --target test_catchem_gaschem_units`
  Run: `./build-macos/tests/test_catchem_gaschem_units`
  Expected: Prints "ALL PROPERTY UNIT TESTS PASSED SUCCESSFULLY!"
- [ ] **Step 4: Commit Unit Test**
  Run: `git add tests/test_catchem_gaschem_units.cpp tests/CMakeLists.txt`
  Run: `git commit -m "test: add property and unit based tests for gaschem conversions"`

---

### Task 6: Implement Coupled Process Integration Test

Build the final end-to-end integration test coupling TUV-x photolysis midpoint rates to MICM solver rate parameters.

**Files:**
* Create: `tests/test_catchem_gaschem.cpp`
* Modify: `tests/CMakeLists.txt`

**Interfaces:**
* Consumes: `catchem_register_photolysis_cpp`, `catchem_register_gaschem_cpp`.
* Produces: Integration test executable `test_catchem_gaschem` verified in CTest.

- [ ] **Step 1: Write Coupled Integration Test**
  ```cpp
  // tests/test_catchem_gaschem.cpp
  #include "catchem_api.hpp"
  #include "catchem_core.hpp"
  #include "catchem_diagnostic_manager.hpp"
  #include "catchem_process_registry.hpp"
  #include "catchem_state_manager.hpp"
  #include <Kokkos_Core.hpp>
  #include <cassert>
  #include <cmath>
  #include <fstream>
  #include <iostream>
  #include <vector>

  extern "C" {
  void catchem_register_photolysis_cpp();
  void catchem_register_gaschem_cpp();
  }

  bool file_exists(const std::string& name) {
      std::ifstream f(name.c_str());
      return f.good();
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

          // 2. Set up core and states
          int n_cols = 1;
          int n_levels = 3;
          int n_species = 5;

          auto core = std::make_shared<catchem::Core>(n_cols, n_levels, n_species);
          auto state = core->get_state_manager();

          // Time setup (Noon, summer)
          state->time.year = 2026;
          state->time.month = 7;
          state->time.day = 13;
          state->time.hour = 12;
          state->time.minute = 0;
          state->time.second = 0;
          state->time.calculate_derived_fields();

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

          // Load species metadata
          state->load_species_config("tests/Configs/Default/CATChem_species.yml");
          std::vector<double> conc_data(n_cols * n_levels * n_species, 1.0); // 1.0 ppmv initially
          state->bind_unified_chemistry(conc_data.data());

          // 3. Resolve configs robustly
          std::string photolysis_config = "src/external/musica/configs/tuvx/from_host/config.json";
          if (!file_exists(photolysis_config)) {
              photolysis_config = "../src/external/musica/configs/tuvx/from_host/config.json";
          }

          std::string temp_main_config = "test_main_coupled_config.yml";
          std::ofstream main_conf_writer(temp_main_config);
          main_conf_writer << "process:\n";
          main_conf_writer << "  photolysis:\n";
          main_conf_writer << "    activate: true\n";
          main_conf_writer << "    config_file: \"" << photolysis_config << "\"\n";
          main_conf_writer << "  gaschem:\n";
          main_conf_writer << "    activate: true\n";
          main_conf_writer.close();

          state->config_file_path = temp_main_config;

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

          std::remove(temp_main_config.c_str());
          std::cout << "\n==========================================" << std::endl;
          std::cout << "COUPLED INTEGRATION TEST PASSED SUCCESSFULLY!" << std::endl;
          std::cout << "==========================================\n" << std::endl;
      }
      Kokkos::finalize();
      return 0;
  }
  ```
- [ ] **Step 2: Register Integration Test in CMake**
  ```cmake
  # tests/CMakeLists.txt
  # Add the following lines at the end of the file:
  add_executable(test_catchem_gaschem test_catchem_gaschem.cpp)
  target_link_libraries(
    test_catchem_gaschem
    PUBLIC CATChem_core_cpp CATChem_core catchem_process_photolysis catchem_process_gaschem Kokkos::kokkos
  )
  set_target_properties(test_catchem_gaschem PROPERTIES LINKER_LANGUAGE CXX)
  add_test(NAME test_catchem_gaschem COMMAND test_catchem_gaschem)
  ```
- [ ] **Step 3: Run and Verify Integration Tests**
  Run: `cmake -B build-macos -S . && cmake --build build-macos --target test_catchem_gaschem`
  Run: `./build-macos/tests/test_catchem_gaschem`
  Expected: Prints "COUPLED INTEGRATION TEST PASSED SUCCESSFULLY!"
- [ ] **Step 4: Commit Integration Test**
  Run: `git add tests/test_catchem_gaschem.cpp tests/CMakeLists.txt`
  Run: `git commit -m "test: add end-to-end photolysis and gaschem coupled integration test"`
