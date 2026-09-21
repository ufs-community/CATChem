# MUSICA TUV-x Photolysis Process Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement a native C++ physical process class `catchem::PhotolysisProcess` extending `catchem::ProcessInterface`, which directly interfaces with the MUSICA TUV-x photolysis library on the CPU host and maps midpoint-interpolated $J$-values to registered diagnostics.

**Architecture:** Native C++ Process implementing the `ProcessInterface`. It synchronizes GPU memory to host on execution, extracts 1D column meteorological variables (altitude/thickness, air density, SZA, O3 density), executes the TUV-x solver column-wise on the host CPU, interpolates the resulting edge-level photolysis rates to layer midpoints, and writes the output back to diagnostics.

**Tech Stack:** C++20, Kokkos (Views), yaml-cpp, and MUSICA TUV-x (Fortran/C++).

## Global Constraints
*   All code must reside in the specified C++ namespaces and follow existing CATChem architectural conventions (such as `catchem::ProcessInterface`).
*   No changes may revert existing modernization work (e.g. Kokkos InteropField mappings).
*   Code must be compile-ready on GCC 14 with C++20 enabled.

---

### Task 1: Extend StateManager and Core for Config Propagation

**Files:**
- Modify: `src/core/catchem_state_manager.hpp`
- Modify: `src/core/catchem_core.cpp`

**Interfaces:**
- Produces: `std::string config_file_path` inside `catchem::StateManager`.

- [ ] **Step 1: Modify catchem_state_manager.hpp**
  Add a `config_file_path` public string member variable to the `StateManager` class declaration to hold the main YAML config file path.

  ```cpp
  // Add this inside the public section of class StateManager in src/core/catchem_state_manager.hpp
  std::string config_file_path;
  ```

- [ ] **Step 2: Modify catchem_core.cpp**
  Propagate the `config_file` parameter to `state_mgr->config_file_path` inside the `Core` constructor.

  ```cpp
  // Inside Core::Core(const std::string& config_file) in src/core/catchem_core.cpp, set:
  state_mgr->config_file_path = config_file;
  ```

- [ ] **Step 3: Verification**
  Run CMake configure to verify the changes don't break existing compilation.

  Run: `docker run --rm -v $(pwd):/opt/catchem/src ufschem-spack-base-ubuntu-gcc-13-dev:latest /bin/bash -c "source /opt/ufschem/spack-stack/setup.sh && spack env activate ufschem && cd /opt/catchem/src && mkdir -p build-test && cd build-test && cmake .. -DCATCHEM_BUILD_TESTING=ON && make -j2 catchem_core"`
  Expected: Successful compilation of the CATChem core library.

- [ ] **Step 4: Commit**
  ```bash
  git add src/core/catchem_state_manager.hpp src/core/catchem_core.cpp
  git commit -m "feat(core): propagate config path to state manager for process initialization"
  ```

---

### Task 2: Create Photolysis Process Header and Registration

**Files:**
- Create: `src/process/photolysis/catchem_process_photolysis.hpp`
- Create: `src/process/photolysis/catchem_process_photolysis.cpp` (Stub registration)

**Interfaces:**
- Produces: Class `catchem::PhotolysisProcess` extending `catchem::ProcessInterface`.
- Produces: `extern "C" void catchem_register_photolysis_cpp()` registration function.

- [ ] **Step 1: Write catchem_process_photolysis.hpp**
  Declare the photolysis process interface class with its standard lifecycle hooks.

  ```cpp
  // src/process/photolysis/catchem_process_photolysis.hpp
  #pragma once
  #include "catchem_process_interface.hpp"
  #include <musica/tuvx/tuvx.hpp>
  #include <memory>
  #include <string>

  namespace catchem {

      class PhotolysisProcess : public ProcessInterface {
      private:
          std::string config_path;
          std::unique_ptr<musica::TUVX> tuvx_instance;
          musica::Mappings photo_mappings;

      public:
          PhotolysisProcess();
          ~PhotolysisProcess() override;

          std::string get_name() const override { return "photolysis"; }

          void init(std::shared_ptr<StateManager> state) override;
          void run(std::shared_ptr<StateManager> state) override;
          void finalize() override;
      };

  } // namespace catchem
  ```

- [ ] **Step 2: Write stub catchem_process_photolysis.cpp**
  Define constructor, destructor, stub hooks, and the `extern "C"` registration entry.

  ```cpp
  // src/process/photolysis/catchem_process_photolysis.cpp
  #include "catchem_process_photolysis.hpp"
  #include "catchem_process_registry.hpp"
  #include <iostream>

  namespace catchem {

      PhotolysisProcess::PhotolysisProcess() : config_path("") {}
      PhotolysisProcess::~PhotolysisProcess() = default;

      void PhotolysisProcess::init(std::shared_ptr<StateManager> state) {
          std::cout << "DEBUG: PhotolysisProcess init" << std::endl;
      }

      void PhotolysisProcess::run(std::shared_ptr<StateManager> state) {
          std::cout << "DEBUG: PhotolysisProcess run" << std::endl;
      }

      void PhotolysisProcess::finalize() {
          std::cout << "DEBUG: PhotolysisProcess finalize" << std::endl;
      }

  } // namespace catchem

  extern "C" {
  void catchem_register_photolysis_cpp() {
      catchem::ProcessRegistry::get_instance().register_process(
          "photolysis", []() { return std::make_shared<catchem::PhotolysisProcess>(); });
  }
  }
  ```

- [ ] **Step 3: Verification**
  Verify the syntax of newly created header and source files.

  Run: `g++ -std=c++20 -Isrc/core -Isrc/external/musica/include -c src/process/photolysis/catchem_process_photolysis.cpp -o /tmp/photo_stub.o && rm /tmp/photo_stub.o`
  Expected: Successful compilation without syntax errors.

- [ ] **Step 4: Commit**
  ```bash
  git add src/process/photolysis/catchem_process_photolysis.hpp src/process/photolysis/catchem_process_photolysis.cpp
  git commit -m "feat(photolysis): add photolysis process class definition and stub registration"
  ```

---

### Task 3: Implement Process Initialization and Dynamic Registration

**Files:**
- Modify: `src/process/photolysis/catchem_process_photolysis.cpp`

**Interfaces:**
- Consumes: `state->config_file_path`
- Produces: Initialized `tuvx_instance` and dynamically registered diagnostics `photolysis_rate_<reaction>` in `DiagnosticManager`.

- [ ] **Step 1: Add yaml-cpp parsing & TUV-x creation to catchem_process_photolysis.cpp**
  Update the `init` method to read `state->config_file_path` and construct the `musica::TUVX` instance.

  Replace `init` block in `src/process/photolysis/catchem_process_photolysis.cpp` with:
  ```cpp
  #include "catchem_diagnostic_manager.hpp"
  #include <yaml-cpp/yaml.h>

  void PhotolysisProcess::init(std::shared_ptr<StateManager> state) {
      if (!state->config_file_path.empty()) {
          try {
              YAML::Node main_config = YAML::LoadFile(state->config_file_path);
              if (main_config["process"] && main_config["process"]["photolysis"]) {
                  auto photo_node = main_config["process"]["photolysis"];
                  if (photo_node["config_file"]) {
                      this->config_path = photo_node["config_file"].as<std::string>();
                  }
              }
          } catch (const std::exception& e) {
              std::cerr << "PhotolysisProcess: Warning: failed to parse main config: " << e.what() << std::endl;
          }
      }

      if (this->config_path.empty()) {
          this->config_path = "src/external/musica/configs/tuvx/tuv_5_4.yml";
      }

      musica::Error err;
      std::unique_ptr<musica::GridMap> grids(musica::CreateGridMap(&err));
      std::unique_ptr<musica::ProfileMap> profiles(musica::CreateProfileMap(&err));
      std::unique_ptr<musica::RadiatorMap> radiators(musica::CreateRadiatorMap(&err));

      tuvx_instance = std::make_unique<musica::TUVX>();
      tuvx_instance->Create(config_path.c_str(), grids.get(), profiles.get(), radiators.get(), &err);

      if (err.status_ != 0) {
          std::cerr << "PhotolysisProcess: Error: Failed to initialize TUV-x! " << err.message_ << std::endl;
          return;
      }

      tuvx_instance->GetPhotolysisRateConstantsOrdering(&photo_mappings, &err);

      if (state->diag_mgr) {
          std::vector<int> dims_2d = {state->n_cols, state->n_levels};
          for (size_t i = 0; i < photo_mappings.size_; ++i) {
              std::string rx_name = photo_mappings.mappings_[i].name_;
              state->diag_mgr->register_field("photolysis_rate_" + rx_name,
                                              "Photolysis rate for " + rx_name,
                                              "s-1", DiagType::FIELD_2D, dims_2d);
          }
      }
  }
  ```

- [ ] **Step 2: Commit**
  ```bash
  git add src/process/photolysis/catchem_process_photolysis.cpp
  git commit -m "feat(photolysis): implement config parsing and dynamic diagnostics registration in init"
  ```

---

### Task 4: Implement Solver Column-Wise Execution Loop

**Files:**
- Modify: `src/process/photolysis/catchem_process_photolysis.cpp`

**Interfaces:**
- Consumes: Meteorological and chemical profiles from `StateManager`
- Produces: Midpoint-interpolated $J$-rates stored in registered diagnostics.

- [ ] **Step 1: Write solver run loop inside catchem_process_photolysis.cpp**
  Populate the `run` method to execute TUV-x column-by-column, map heights and densities, execute the solver, interpolate edge rates to midpoints, and populate diagnostic arrays.

  Replace the `run` block in `src/process/photolysis/catchem_process_photolysis.cpp` with:
  ```cpp
  #include <cmath>
  #include <algorithm>

  void PhotolysisProcess::run(std::shared_ptr<StateManager> state) {
      if (!tuvx_instance) return;

      state->sync_to_host();

      int i_o3 = -1;
      for (size_t i = 0; i < state->chem.species_list.size(); ++i) {
          if (state->chem.species_list[i].short_name == "O3") {
              i_o3 = i;
              break;
          }
      }

      musica::Error err;
      int num_reactions = tuvx_instance->GetPhotolysisRateConstantCount();

      std::unique_ptr<musica::ProfileMap> profiles(tuvx_instance->GetProfileMap(&err));
      if (err.status_ != 0) {
          std::cerr << "PhotolysisProcess: Error getting ProfileMap: " << err.message_ << std::endl;
          return;
      }

      musica::Profile* profile_air = profiles->GetProfile("air", "molecule cm-3", &err);
      musica::Profile* profile_o2  = profiles->GetProfile("O2", "molecule cm-3", &err);
      musica::Profile* profile_o3  = profiles->GetProfile("O3", "molecule cm-3", &err);

      std::vector<double> air_profile(state->n_levels, 0.0);
      std::vector<double> o2_profile(state->n_levels, 0.0);
      std::vector<double> o3_profile(state->n_levels, 0.0);

      for (int i_col = 0; i_col < state->n_cols; ++i_col) {
          double lat_deg = state->met.LAT ? state->met.LAT->host_view(i_col, 0) : 0.0;
          double lon_deg = state->met.LON ? state->met.LON->host_view(i_col, 0) : 0.0;
          double cos_sza = state->time.get_cos_sza(lat_deg, lon_deg, true);
          double sza_rad = std::acos(std::max(-1.0, std::min(1.0, cos_sza)));

          for (int i_lvl = 0; i_lvl < state->n_levels; ++i_lvl) {
              double airden_kg_m3 = state->met.AIRDEN ? state->met.AIRDEN->host_view(i_col, i_lvl, 0) : 1.2;
              air_profile[i_lvl] = airden_kg_m3 * 2.079153e19;
              o2_profile[i_lvl] = air_profile[i_lvl] * 0.2095;

              if (i_o3 >= 0 && state->chem.conc) {
                  o3_profile[i_lvl] = state->chem.conc->host_view(i_col, i_lvl, i_o3);
              } else {
                  o3_profile[i_lvl] = air_profile[i_lvl] * 3e-7;
              }
          }

          profile_air->SetMidpointValues(air_profile.data(), state->n_levels, &err);
          profile_o2->SetMidpointValues(o2_profile.data(), state->n_levels, &err);
          profile_o3->SetMidpointValues(o3_profile.data(), state->n_levels, &err);

          std::vector<double> edge_photolysis_rates((state->n_levels + 1) * num_reactions, 0.0);
          std::vector<double> edge_heating_rates((state->n_levels + 1) * tuvx_instance->GetHeatingRateCount(), 0.0);

          tuvx_instance->Run(
              sza_rad,
              1.0,
              edge_photolysis_rates.data(),
              edge_heating_rates.data(),
              nullptr,
              nullptr,
              nullptr,
              &err);

          if (err.status_ != 0) {
              std::cerr << "PhotolysisProcess: Solver error in column " << i_col << ": " << err.message_ << std::endl;
              continue;
          }

          if (state->diag_mgr) {
              for (size_t rx_idx = 0; rx_idx < photo_mappings.size_; ++rx_idx) {
                  std::string rx_name = photo_mappings.mappings_[rx_idx].name_;
                  std::string diag_name = "photolysis_rate_" + rx_name;
                  double* diag_ptr = static_cast<double*>(state->diag_mgr->get_host_pointer(diag_name));

                  if (diag_ptr) {
                      for (int i_lvl = 0; i_lvl < state->n_levels; ++i_lvl) {
                          int idx_edge1 = rx_idx * (state->n_levels + 1) + i_lvl;
                          int idx_edge2 = rx_idx * (state->n_levels + 1) + (i_lvl + 1);

                          double rate_midpoint = 0.5 * (edge_photolysis_rates[idx_edge1] + edge_photolysis_rates[idx_edge2]);

                          int diag_idx = i_lvl * state->n_cols + i_col;
                          diag_ptr[diag_idx] = rate_midpoint;
                      }
                  }
              }
          }
      }

      state->sync_to_device();
      if (state->diag_mgr) {
          state->diag_mgr->sync_to_device();
      }
  }
  ```

- [ ] **Step 2: Commit**
  ```bash
  git add src/process/photolysis/catchem_process_photolysis.cpp
  git commit -m "feat(photolysis): implement column-wise SZA calculation and TUV-x Run step"
  ```

---

### Task 5: Integrate Build System

**Files:**
- Create: `src/process/photolysis/CMakeLists.txt`
- Modify: `src/process/CMakeLists.txt`

**Interfaces:**
- Consumes: Target `musica_tuvx`
- Produces: CMake library target `catchem_process_photolysis` and adds it to the list of build-subdirectories.

- [ ] **Step 1: Create src/process/photolysis/CMakeLists.txt**
  Configure the build targets.

  ```cmake
  # src/process/photolysis/CMakeLists.txt
  add_library(catchem_process_photolysis STATIC
    catchem_process_photolysis.cpp
  )

  target_include_directories(catchem_process_photolysis PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src/core
    ${PROJECT_SOURCE_DIR}/src/external/musica/include
  )

  target_link_libraries(catchem_process_photolysis PUBLIC
    catchem_core
    musica_tuvx
    yaml-cpp
  )
  ```

- [ ] **Step 2: Modify src/process/CMakeLists.txt**
  Include the new photolysis folder in CMake.

  ```cmake
  # Add this to src/process/CMakeLists.txt
  add_subdirectory(photolysis)
  ```

- [ ] **Step 3: Verification**
  Reconfigure CMake and run target compilation using Docker.

  Run: `docker run --rm -v $(pwd):/opt/catchem/src ufschem-spack-base-ubuntu-gcc-13-dev:latest /bin/bash -c "source /opt/ufschem/spack-stack/setup.sh && spack env activate ufschem && cd /opt/catchem/src && mkdir -p build-test && cd build-test && cmake .. -DCATCHEM_BUILD_TESTING=ON && make -j2 catchem_process_photolysis"`
  Expected: Successful compilation of `catchem_process_photolysis` library without errors.

- [ ] **Step 4: Commit**
  ```bash
  git add src/process/CMakeLists.txt src/process/photolysis/CMakeLists.txt
  git commit -m "build(cmake): integrate photolysis process into build system targets"
  ```

---

### Task 6: Unit and Integration Testing

**Files:**
- Create: `tests/test_catchem_photolysis.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `catchem_register_photolysis_cpp`, `catchem_core_add_process_by_name`.
- Produces: Compiled executable test `test_catchem_photolysis` running photolysis and validating rates.

- [ ] **Step 1: Create test_catchem_photolysis.cpp**
  Set up the C++ integration and unit test.

  ```cpp
  // tests/test_catchem_photolysis.cpp
  #include <gtest/gtest.h>
  #include "catchem_core.hpp"
  #include "catchem_api.hpp"
  #include "catchem_process_registry.hpp"
  #include <memory>
  #include <vector>

  extern "C" {
  void catchem_register_photolysis_cpp();
  }

  TEST(PhotolysisTest, Registration) {
      catchem_register_photolysis_cpp();
      EXPECT_TRUE(catchem::ProcessRegistry::get_instance().has_process("photolysis"));
  }

  TEST(PhotolysisTest, InitAndRun) {
      catchem_register_photolysis_cpp();

      // Instantiate Core
      auto core = std::make_shared<catchem::Core>(1, 64, 5); // 1 column, 64 levels, 5 species

      // Load standard Met State
      double temperature[64];
      double airden[64];
      double pedge[65];
      double bxheight[64];
      double lat[1] = {40.0};
      double lon[1] = {-105.0};

      for (int i = 0; i < 64; ++i) {
          temperature[i] = 280.0 - 0.5 * i;
          airden[i] = 1.2 * std::exp(-i / 10.0);
          bxheight[i] = 100.0;
          pedge[i] = 101300.0 * std::exp(-i / 10.0);
      }
      pedge[64] = 101300.0 * std::exp(-64 / 10.0);

      auto state = core->get_state_manager();
      state->bind_met_field_2d("LAT", lat);
      state->bind_met_field_2d("LON", lon);
      state->bind_met_field_3d("T", temperature);
      state->bind_met_field_3d("AIRDEN", airden);
      state->bind_met_field_3d("PEDGE", pedge);
      state->bind_met_field_3d("BXHEIGHT", bxheight);

      // Add and run photolysis process
      auto process = catchem::ProcessRegistry::get_instance().create("photolysis");
      process->init(state);
      core->add_process(process);

      // Run timestep
      core->run_timestep(60.0);

      // Validate that photolysis diagnostic rates are initialized and populated (non-zero)
      double* o3_photo_rate = static_cast<double*>(core->get_diagnostic_manager()->get_host_pointer("photolysis_rate_O3+hv->O2+O(1D)"));
      EXPECT_NE(o3_photo_rate, nullptr);

      // Ensure we have non-zero solar photolysis rates
      double sum_rates = 0.0;
      for (int i = 0; i < 64; ++i) {
          sum_rates += o3_photo_rate[i];
      }
      EXPECT_GT(sum_rates, 0.0);
  }

  int main(int argc, char **argv) {
      ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
  }
  ```

- [ ] **Step 2: Modify tests/CMakeLists.txt**
  Add our new test suite to CMake.

  ```cmake
  # Add this to tests/CMakeLists.txt
  add_executable(test_catchem_photolysis test_catchem_photolysis.cpp)
  target_link_libraries(test_catchem_photolysis PRIVATE
    catchem_core
    catchem_process_photolysis
    gtest
    gtest_main
  )
  add_test(NAME test_catchem_photolysis COMMAND test_catchem_photolysis)
  ```

- [ ] **Step 3: Verification**
  Run our compilation and execution tests.

  Run: `docker run --rm -v $(pwd):/opt/catchem/src ufschem-spack-base-ubuntu-gcc-13-dev:latest /bin/bash -c "source /opt/ufschem/spack-stack/setup.sh && spack env activate ufschem && cd /opt/catchem/src && mkdir -p build-test && cd build-test && cmake .. -DCATCHEM_BUILD_TESTING=ON && make -j2 test_catchem_photolysis && ctest -R test_catchem_photolysis --output-on-failure"`
  Expected: Compiled and passed the `test_catchem_photolysis` test.

- [ ] **Step 4: Commit**
  ```bash
  git add tests/CMakeLists.txt tests/test_catchem_photolysis.cpp
  git commit -m "test(photolysis): add integration tests for photolysis process and verify J-rate computation"
  ```
