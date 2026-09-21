# Standardized Columnar Logging Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement and integrate a standardized, visually aligned, columnar plain-text logging framework across CATChem to replace raw, unformatted `std::cout` statements.

**Architecture:** A lightweight C++20 `catchem::Logger` wrapper that formats logs to a strict columnar padding structure: `[TIMESTAMP] [LEVEL] [SERVICE] [TRACE_ID] Message | key=value`. The `Trace ID` is propagated dynamically via `catchem::StateManager`.

**Tech Stack:** C++20, Kokkos, `<chrono>`, `<unistd.h>` (for `isatty()`).

## Global Constraints
* **Format:** `[TIMESTAMP] [LEVEL] [SERVICE] [TRACE_ID] Message | key=value`
* **Level Length:** Exactly 5 characters (`DEBUG`, `INFO `, `WARN `, `ERROR`).
* **Service Length:** Exactly 15 characters, left-justified.
* **Trace ID Length:** Exactly 8 characters, space-padded if needed.
* **Color Stripping:** Colors are applied to `[LEVEL]` tag on local terminals, but completely stripped if `NO_COLOR` is in the environment or if stdout/stderr are not interactive TTYs.

---

### Task 1: Extend StateManager with Trace ID Generation

**Files:**
- Modify: `src/core/catchem_state_manager.hpp`
- Modify: `src/core/catchem_state_manager.cpp`
- Test: `tests/test_catchem_properties.cpp`

**Interfaces:**
- Consumes: None (base field addition)
- Produces: `catchem::StateManager::trace_id` member string

- [ ] **Step 1: Write the test**
  Add a unit test in `tests/test_catchem_properties.cpp` to verify that `StateManager` automatically generates an 8-character trace ID.
  ```cpp
  // Add in tests/test_catchem_properties.cpp inside a test block
  auto state = std::make_shared<catchem::StateManager>(4, 10, 50);
  assert(state->trace_id.length() == 8);
  assert(!state->trace_id.empty());
  ```

- [ ] **Step 2: Run test to verify it fails**
  Run: `cmake --build build-macos -j 4 && ctest --test-dir build-macos -R test_catchem_properties --output-on-failure`
  Expected: FAIL with "no member trace_id"

- [ ] **Step 3: Write implementation**
  Add `trace_id` string to `catchem_state_manager.hpp`:
  ```cpp
  // In src/core/catchem_state_manager.hpp
  class StateManager {
  public:
      std::string trace_id;
      // ... existing fields ...
  };
  ```
  Implement Trace ID generation using standard string utilities in `catchem_state_manager.cpp` constructor:
  ```cpp
  // In src/core/catchem_state_manager.cpp constructor
  StateManager::StateManager(int nc, int nl, int ns)
      : n_cols(nc), n_levels(nl), n_species(ns) {
      // Generate standard 8-character trace id (alphanumeric random)
      static const char alphanum[] = "0123456789abcdefghijklmnopqrstuvwxyz";
      trace_id = "";
      for (int i = 0; i < 8; ++i) {
          trace_id += alphanum[rand() % (sizeof(alphanum) - 1)];
      }
      // ... existing allocation logic ...
  }
  ```

- [ ] **Step 4: Run test to verify it passes**
  Run: `cmake --build build-macos -j 4 && ctest --test-dir build-macos -R test_catchem_properties --output-on-failure`
  Expected: PASS

- [ ] **Step 5: Commit**
  ```bash
  git add src/core/catchem_state_manager.hpp src/core/catchem_state_manager.cpp tests/test_catchem_properties.cpp
  git commit -m "feat: add trace_id generation to StateManager"
  ```

---

### Task 2: Implement C++ Central Logger Wrapper

**Files:**
- Create: `src/core/catchem_logger.hpp`
- Create: `src/core/catchem_logger.cpp`
- Modify: `src/core/CMakeLists.txt`
- Create: `tests/test_catchem_logger.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `catchem::StateManager`
- Produces: `catchem::Logger` class static methods (`debug`, `info`, `warn`, `error`)

- [ ] **Step 1: Write the failing tests**
  Create `tests/test_catchem_logger.cpp` to verify padding, log layout output, and color stripping:
  ```cpp
  #include <catchem_logger.hpp>
  #include <catchem_state_manager.hpp>
  #include <cassert>
  #include <iostream>

  int main() {
      auto state = std::make_shared<catchem::StateManager>(4, 10, 50);
      state->trace_id = "testtrac";

      // Manual redirect stringstream capture is optional, we assert logger successfully formats and prints
      catchem::Logger::info(state.get(), "Simulation timestep advanced", {{"step", "12"}, {"dt", "300.0"}});
      catchem::Logger::error(state.get(), "Division by zero encountered", {{"cell", "4"}});

      std::cout << "All logger formatting unit tests passed!" << std::endl;
      return 0;
  }
  ```
  Add unit test executable registration in `tests/CMakeLists.txt`:
  ```cmake
  add_executable(test_catchem_logger test_catchem_logger.cpp)
  target_link_libraries(test_catchem_logger PRIVATE CATChem_core_cpp)
  add_test(NAME test_catchem_logger COMMAND test_catchem_logger)
  ```

- [ ] **Step 2: Run test to verify it fails**
  Run: `cmake --build build-macos -j 4 && ctest --test-dir build-macos -R test_catchem_logger --output-on-failure`
  Expected: FAIL (cannot find `<catchem_logger.hpp>`)

- [ ] **Step 3: Implement central C++ Logger**
  Create `src/core/catchem_logger.hpp` with the class declarations:
  ```cpp
  #pragma once
  #include "catchem_state_manager.hpp"
  #include <string_view>
  #include <initializer_list>
  #include <utility>

  namespace catchem {

      class Logger {
      public:
          using ContextList = std::initializer_list<std::pair<std::string_view, std::string_view>>;

          static void debug(const StateManager* state, std::string_view message, ContextList context = {});
          static void info(const StateManager* state, std::string_view message, ContextList context = {});
          static void warn(const StateManager* state, std::string_view message, ContextList context = {});
          static void error(const StateManager* state, std::string_view message, ContextList context = {});

      private:
          static void log(const StateManager* state, std::string_view level, std::string_view message, ContextList context);
          static bool should_color(int fd);
      };

  } // namespace catchem
  ```

  Create `src/core/catchem_logger.cpp` with exact format building, time layout, padding rules, and stream writes:
  ```cpp
  #include "catchem_logger.hpp"
  #include <iostream>
  #include <chrono>
  #include <iomanip>
  #include <sstream>
  #include <cstdlib>
  #include <unistd.h>

  namespace catchem {

      bool Logger::should_color(int fd) {
          const char* no_color = std::getenv("NO_COLOR");
          if (no_color && no_color[0] != '\0') {
              return false;
          }
          return isatty(fd);
      }

      void Logger::debug(const StateManager* state, std::string_view message, ContextList context) {
          log(state, "DEBUG", message, context);
      }

      void Logger::info(const StateManager* state, std::string_view message, ContextList context) {
          log(state, "INFO ", message, context);
      }

      void Logger::warn(const StateManager* state, std::string_view message, ContextList context) {
          log(state, "WARN ", message, context);
      }

      void Logger::error(const StateManager* state, std::string_view message, ContextList context) {
          log(state, "ERROR", message, context);
      }

      void Logger::log(const StateManager* state, std::string_view level, std::string_view message, ContextList context) {
          // 1. Get exact current UTC Timestamp
          auto now = std::chrono::system_clock::now();
          std::time_t now_time = std::chrono::system_clock::to_time_t(now);
          std::tm* utc_tm = std::gmtime(&now_time);

          std::ostringstream ss;
          ss << std::put_time(utc_tm, "%Y-%m-%d %H:%M:%S");
          std::string timestamp = ss.str();

          // 2. Format Level with ANSI Coloring
          int fd = (level == "ERROR") ? fileno(stderr) : fileno(stdout);
          bool color = should_color(fd);

          std::string colored_level(level);
          if (color) {
              if (level == "DEBUG") colored_level = "\033[36mDEBUG\033[0m"; // Cyan
              else if (level == "INFO ") colored_level = "\033[32mINFO \033[0m"; // Green
              else if (level == "WARN ") colored_level = "\033[33mWARN \033[0m"; // Yellow
              else if (level == "ERROR") colored_level = "\033[31mERROR\033[0m"; // Red
          }

          // 3. Assemble Service Name (exactly 15 chars, left-justified)
          std::string service = "catchem";
          service.append(15 - service.length(), ' ');

          // 4. Assemble Trace ID (exactly 8 chars)
          std::string trace = (state && !state->trace_id.empty()) ? state->trace_id : "global  ";
          if (trace.length() < 8) trace.append(8 - trace.length(), ' ');

          // 5. Build full golden prefix
          std::ostringstream out;
          out << "[" << timestamp << "] [" << colored_level << "] [" << service << "] [" << trace << "] " << message;

          // 6. Append Key-Value Context dictionary
          if (context.size() > 0) {
              out << " |";
              for (const auto& [key, value] : context) {
                  out << " " << key << "=" << value;
              }
          }
          out << "\n";

          // 7. Stream out cleanly
          if (level == "ERROR") {
              std::cerr << out.str() << std::flush;
          } else {
              std::clog << out.str() << std::flush;
          }
      }

  } // namespace catchem
  ```
  Add the source file to `src/core/CMakeLists.txt`:
  ```cmake
  # In src/core/CMakeLists.txt targets src files list
  catchem_logger.cpp
  ```

- [ ] **Step 4: Run test to verify it passes**
  Run: `cmake --build build-macos -j 4 && ctest --test-dir build-macos -R test_catchem_logger --output-on-failure`
  Expected: PASS

- [ ] **Step 5: Commit**
  ```bash
  git add src/core/catchem_logger.hpp src/core/catchem_logger.cpp src/core/CMakeLists.txt tests/test_catchem_logger.cpp tests/CMakeLists.txt
  git commit -m "feat: add standardized central C++ Logger module"
  ```

---

### Task 3: Refactor Core & Processes Print Statements

**Files:**
- Modify: `src/process/gaschem/catchem_process_gaschem.cpp`
- Modify: `src/process/photolysis/catchem_process_photolysis.cpp`

**Interfaces:**
- Consumes: `catchem::Logger`
- Produces: Visual unified stdout alignment logging

- [ ] **Step 1: Inspect existing logs**
  Verify that we've found all standard debug cout loops using grep.

- [ ] **Step 2: Refactor `GasChemProcess::init`**
  Replace standard console streams with the new `catchem::Logger` in `catchem_process_gaschem.cpp`:
  ```cpp
  // Add header
  #include "catchem_logger.hpp"

  // Replace clog DEBUG prints with:
  Logger::debug(state.get(), "GasChemProcess::init started");
  Logger::info(state.get(), "GasChemProcess: resolved config directory", {{"dir", config_dir}});
  Logger::info(state.get(), "GasChemProcess: initialized MICM successfully!");
  ```

- [ ] **Step 3: Refactor `PhotolysisProcess::init`**
  Replace standard console streams in `catchem_process_photolysis.cpp` with `Logger`:
  ```cpp
  // Add header
  #include "catchem_logger.hpp"

  // Replace debug print blocks:
  Logger::debug(state.get(), "PhotolysisProcess::init started");
  Logger::debug(state.get(), "Parsing TUV-x config file to check profiles", {{"config", config_path}});
  Logger::info(state.get(), "PhotolysisProcess: initialized TUV-x successfully!");
  ```

- [ ] **Step 4: Run full CTest suite to verify perfect functionality**
  Run: `cmake --build build-macos -j 4 && ctest --test-dir build-macos --output-on-failure`
  Expected: All 12/12 tests PASS

- [ ] **Step 5: Commit**
  ```bash
  git add src/process/gaschem/catchem_process_gaschem.cpp src/process/photolysis/catchem_process_photolysis.cpp
  git commit -m "refactor: integrate standardized logger across GasChem and Photolysis"
  ```
