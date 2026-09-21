# Plan: Modular C++ Process Layout & Linker-Safe Registration

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Re-architect physical process layout by moving C++ process controllers from `@src/core` to their respective `@src/process/<name>` subdirectories, using CMake mixed-language targets, and introducing a linker-safe Dynamic Process Registration mechanism.

**Architecture:** Moving concrete physics out of core. `@src/core` will retain the abstract `catchem::ProcessInterface` and `ProcessRegistry`, while concrete process C++ implementations reside under `@src/process/<name>/**` and register their instantiators to the Core at runtime via standard C dynamic bridges.

**Tech Stack:** C++20, Fortran 2018, CMake, Kokkos, ISO_C_BINDING.

## Global Constraints
- **C++ Standard:** Target C++20.
- **Layout Alignment:** Fortran Column-Major Layout (`Kokkos::LayoutLeft`) across pointer boundaries.
- **Modularity:** `@src/core` must have zero knowledge of individual physical processes.
- **Linker Safety:** Standard static library dynamic registration relies on side-effects which are optimized out by the linker; explicit `extern "C"` registration routines must be declared and called directly from Fortran initializers.

---

## Tasks

### Task 1: Relocate Process C++ Sources & Update Core CMakeLists.txt

**Files:**
- Create: `src/process/settling/catchem_process_settling.hpp`, `src/process/settling/catchem_process_settling.cpp`
- Create: `src/process/seasalt/catchem_process_seasalt.hpp`, `src/process/seasalt/catchem_process_seasalt.cpp`
- Create: `src/process/drydep/catchem_process_drydep.hpp`, `src/process/drydep/catchem_process_drydep.cpp`
- Create: `src/process/wetdep/catchem_process_wetdep.hpp`, `src/process/wetdep/catchem_process_wetdep.cpp`
- Create: `src/process/so4chem/catchem_process_so4chem.hpp`, `src/process/so4chem/catchem_process_so4chem.cpp`
- Delete: `src/core/catchem_process_settling.*`, `src/core/catchem_process_seasalt.*`, `src/core/catchem_process_drydep.*`, `src/core/catchem_process_wetdep.*`, `src/core/catchem_process_so4chem.*`
- Modify: `src/core/CMakeLists.txt`

**Interfaces:**
- Consumes: None
- Produces: Cleaned `CATChem_core_cpp` that builds without concrete physics definitions.

- [ ] **Step 1: Move files using standard git command**
Run:
```bash
git mv src/core/catchem_process_settling.hpp src/process/settling/catchem_process_settling.hpp
git mv src/core/catchem_process_settling.cpp src/process/settling/catchem_process_settling.cpp

git mv src/core/catchem_process_seasalt.hpp src/process/seasalt/catchem_process_seasalt.hpp
git mv src/core/catchem_process_seasalt.cpp src/process/seasalt/catchem_process_seasalt.cpp

git mv src/core/catchem_process_drydep.hpp src/process/drydep/catchem_process_drydep.hpp
git mv src/core/catchem_process_drydep.cpp src/process/drydep/catchem_process_drydep.cpp

git mv src/core/catchem_process_wetdep.hpp src/process/wetdep/catchem_process_wetdep.hpp
git mv src/core/catchem_process_wetdep.cpp src/process/wetdep/catchem_process_wetdep.cpp

git mv src/core/catchem_process_so4chem.hpp src/process/so4chem/catchem_process_so4chem.hpp
git mv src/core/catchem_process_so4chem.cpp src/process/so4chem/catchem_process_so4chem.cpp
```

- [ ] **Step 2: Update `src/core/CMakeLists.txt` to remove concrete sources**
Modify `src/core/CMakeLists.txt` around lines 59-70:
```cmake
# Before:
  set(
    _cpp_core_srcs
    catchem_core.cpp
    catchem_api.cpp
    catchem_diagnostic.cpp
    catchem_process_settling.cpp
    catchem_process_seasalt.cpp
    catchem_process_drydep.cpp
    catchem_process_wetdep.cpp
    catchem_process_so4chem.cpp
  )

# After:
  set(
    _cpp_core_srcs
    catchem_core.cpp
    catchem_api.cpp
    catchem_diagnostic.cpp
  )
```

- [ ] **Step 3: Commit relocation**
```bash
git add src/core/CMakeLists.txt src/process/
git commit -m "refactor(core): relocate physical process sources from core to process directories"
```

---

### Task 2: Reconfigure Process CMakeLists.txt files for Mixed C++/Fortran Target Compilations

**Files:**
- Modify: `src/process/settling/CMakeLists.txt`
- Modify: `src/process/seasalt/CMakeLists.txt`
- Modify: `src/process/drydep/CMakeLists.txt`
- Modify: `src/process/wetdep/CMakeLists.txt`
- Modify: `src/process/so4chem/CMakeLists.txt`

**Interfaces:**
- Consumes: C++20 process controllers relocated in Task 1.
- Produces: Dynamic library compilation supporting both Fortran and C++ files when `ENABLE_KOKKOS` is active.

- [ ] **Step 1: Modify `src/process/settling/CMakeLists.txt`**
Update to conditionally include `catchem_process_settling.cpp` and link with `Kokkos::kokkos` when `ENABLE_KOKKOS` is enabled.
```cmake
# Around Settling All Sources
set(SETTLING_ALL_SOURCES ${SETTLING_PROCESS_SOURCES} ${SETTLING_SCHEME_SOURCES})
if(ENABLE_KOKKOS)
  list(APPEND SETTLING_ALL_SOURCES catchem_process_settling.cpp)
endif()

set(_lib CATChem_process_settling)
add_library(${_lib} ${SETTLING_ALL_SOURCES})

target_link_libraries(${_lib} PUBLIC CATChem_core)
if(ENABLE_KOKKOS)
  target_link_libraries(${_lib} PUBLIC Kokkos::kokkos)
  target_compile_definitions(${_lib} PRIVATE ENABLE_KOKKOS)
endif()
```

- [ ] **Step 2: Modify `src/process/seasalt/CMakeLists.txt`**
Apply same structure:
```cmake
set(SEASALT_ALL_SOURCES ${SEASALT_PROCESS_SOURCES} ${SEASALT_SCHEME_SOURCES})
if(ENABLE_KOKKOS)
  list(APPEND SEASALT_ALL_SOURCES catchem_process_seasalt.cpp)
endif()

set(_lib CATChem_process_seasalt)
add_library(${_lib} ${SEASALT_ALL_SOURCES})

target_link_libraries(${_lib} PUBLIC CATChem_core)
if(ENABLE_KOKKOS)
  target_link_libraries(${_lib} PUBLIC Kokkos::kokkos)
  target_compile_definitions(${_lib} PRIVATE ENABLE_KOKKOS)
endif()
```

- [ ] **Step 3: Modify `src/process/drydep/CMakeLists.txt`**
Apply same structure:
```cmake
set(DRYDEP_ALL_SOURCES ${DRYDEP_PROCESS_SOURCES} ${DRYDEP_SCHEME_SOURCES})
if(ENABLE_KOKKOS)
  list(APPEND DRYDEP_ALL_SOURCES catchem_process_drydep.cpp)
endif()

set(_lib CATChem_process_drydep)
add_library(${_lib} ${DRYDEP_ALL_SOURCES})

target_link_libraries(${_lib} PUBLIC CATChem_core)
if(ENABLE_KOKKOS)
  target_link_libraries(${_lib} PUBLIC Kokkos::kokkos)
  target_compile_definitions(${_lib} PRIVATE ENABLE_KOKKOS)
endif()
```

- [ ] **Step 4: Modify `src/process/wetdep/CMakeLists.txt`**
Apply same structure:
```cmake
set(WETDEP_ALL_SOURCES ${WETDEP_PROCESS_SOURCES} ${WETDEP_SCHEME_SOURCES})
if(ENABLE_KOKKOS)
  list(APPEND WETDEP_ALL_SOURCES catchem_process_wetdep.cpp)
endif()

set(_lib CATChem_process_wetdep)
add_library(${_lib} ${WETDEP_ALL_SOURCES})

target_link_libraries(${_lib} PUBLIC CATChem_core)
if(ENABLE_KOKKOS)
  target_link_libraries(${_lib} PUBLIC Kokkos::kokkos)
  target_compile_definitions(${_lib} PRIVATE ENABLE_KOKKOS)
endif()
```

- [ ] **Step 5: Modify `src/process/so4chem/CMakeLists.txt`**
Apply same structure:
```cmake
set(SO4CHEM_ALL_SOURCES ${SO4CHEM_PROCESS_SOURCES} ${SO4CHEM_SCHEME_SOURCES})
if(ENABLE_KOKKOS)
  list(APPEND SO4CHEM_ALL_SOURCES catchem_process_so4chem.cpp)
endif()

set(_lib CATChem_process_so4chem)
add_library(${_lib} ${SO4CHEM_ALL_SOURCES})

target_link_libraries(${_lib} PUBLIC CATChem_core)
if(ENABLE_KOKKOS)
  target_link_libraries(${_lib} PUBLIC Kokkos::kokkos)
  target_compile_definitions(${_lib} PRIVATE ENABLE_KOKKOS)
endif()
```

- [ ] **Step 6: Commit build changes**
```bash
git add src/process/**/CMakeLists.txt
git commit -m "build(process): enable mixed C++/Fortran compilation conditionally on ENABLE_KOKKOS"
```

---

### Task 3: Implement C-Linkage Dynamic Registration & Fortran Activation Calls

**Files:**
- Modify: `src/process/settling/catchem_process_settling.cpp` & `ProcessSettlingInterface_Mod.F90`
- Modify: `src/process/seasalt/catchem_process_seasalt.cpp` & `ProcessSeaSaltInterface_Mod.F90`
- Modify: `src/process/drydep/catchem_process_drydep.cpp` & `ProcessDryDepInterface_Mod.F90`
- Modify: `src/process/wetdep/catchem_process_wetdep.cpp` & `ProcessWetDepInterface_Mod.F90`
- Modify: `src/process/so4chem/catchem_process_so4chem.cpp` & `ProcessSO4chemInterface_Mod.F90`

**Interfaces:**
- Consumes: C++ global `ProcessRegistry::get_instance()`.
- Produces: Linker-safe invocation that binds each process dynamically at startup.

- [ ] **Step 1: Modify `src/process/settling/catchem_process_settling.cpp`**
Add C linkage helper at the end of the file:
```cpp
#include "catchem_process_registry.hpp"

extern "C" {
void catchem_register_settling_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "settling",
        []() { return std::make_shared<catchem::SettlingProcess>(); }
    );
}
}
```

- [ ] **Step 2: Modify `src/process/settling/ProcessSettlingInterface_Mod.F90`**
Under `process_init`, declare and invoke this callback:
```fortran
   subroutine process_init(this, container, rc)
      ...
#ifdef ENABLE_KOKKOS
      block
         interface
            subroutine catchem_register_settling_cpp() bind(c, name="catchem_register_settling_cpp")
            end subroutine catchem_register_settling_cpp
         end interface
         call catchem_register_settling_cpp()
      end block
#endif
```

- [ ] **Step 3: Modify `src/process/seasalt/catchem_process_seasalt.cpp` and `ProcessSeaSaltInterface_Mod.F90`**
Apply the same pattern.
C++:
```cpp
#include "catchem_process_registry.hpp"

extern "C" {
void catchem_register_seasalt_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "seasalt",
        []() { return std::make_shared<catchem::SeaSaltProcess>(); }
    );
}
}
```
Fortran in `process_init`:
```fortran
#ifdef ENABLE_KOKKOS
      block
         interface
            subroutine catchem_register_seasalt_cpp() bind(c, name="catchem_register_seasalt_cpp")
            end subroutine catchem_register_seasalt_cpp
         end interface
         call catchem_register_seasalt_cpp()
      end block
#endif
```

- [ ] **Step 4: Modify `src/process/drydep/catchem_process_drydep.cpp` and `ProcessDryDepInterface_Mod.F90`**
Apply the same pattern.
C++:
```cpp
#include "catchem_process_registry.hpp"

extern "C" {
void catchem_register_drydep_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "drydep",
        []() { return std::make_shared<catchem::DryDepProcess>(); }
    );
}
}
```
Fortran in `process_init`:
```fortran
#ifdef ENABLE_KOKKOS
      block
         interface
            subroutine catchem_register_drydep_cpp() bind(c, name="catchem_register_drydep_cpp")
            end subroutine catchem_register_drydep_cpp
         end interface
         call catchem_register_drydep_cpp()
      end block
#endif
```

- [ ] **Step 5: Modify `src/process/wetdep/catchem_process_wetdep.cpp` and `ProcessWetDepInterface_Mod.F90`**
Apply the same pattern.
C++:
```cpp
#include "catchem_process_registry.hpp"

extern "C" {
void catchem_register_wetdep_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "wetdep",
        []() { return std::make_shared<catchem::WetDepProcess>(); }
    );
}
}
```
Fortran in `process_init`:
```fortran
#ifdef ENABLE_KOKKOS
      block
         interface
            subroutine catchem_register_wetdep_cpp() bind(c, name="catchem_register_wetdep_cpp")
            end subroutine catchem_register_wetdep_cpp
         end interface
         call catchem_register_wetdep_cpp()
      end block
#endif
```

- [ ] **Step 6: Modify `src/process/so4chem/catchem_process_so4chem.cpp` and `ProcessSO4chemInterface_Mod.F90`**
Apply the same pattern.
C++:
```cpp
#include "catchem_process_registry.hpp"

extern "C" {
void catchem_register_so4chem_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "so4chem",
        []() { return std::make_shared<catchem::SO4chemProcess>(); }
    );
}
}
```
Fortran in `process_init`:
```fortran
#ifdef ENABLE_KOKKOS
      block
         interface
            subroutine catchem_register_so4chem_cpp() bind(c, name="catchem_register_so4chem_cpp")
            end subroutine catchem_register_so4chem_cpp
         end interface
         call catchem_register_so4chem_cpp()
      end block
#endif
```

- [ ] **Step 7: Commit Registration Pattern**
```bash
git add src/process/
git commit -m "feat(process): implement C registration wrappers and call from Fortran initialization"
```

---

### Task 4: Extend Core C-API for Dynamic Process Instantiation

**Files:**
- Modify: `src/core/catchem_api.hpp`
- Modify: `src/core/catchem_api.cpp`

**Interfaces:**
- Consumes: `catchem::ProcessRegistry`
- Produces: `void catchem_core_add_process_by_name(void* core_ptr, const char* name)` C symbol.

- [ ] **Step 1: Modify `src/core/catchem_api.hpp`**
Add symbol declaration inside `extern "C" {}`:
```cpp
void catchem_core_add_process_by_name(void* core_ptr, const char* name);
```

- [ ] **Step 2: Modify `src/core/catchem_api.cpp`**
Implement the registry extraction:
```cpp
#include "catchem_process_registry.hpp"

void catchem_core_add_process_by_name(void* core_ptr, const char* name) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    core->add_process(catchem::ProcessRegistry::get_instance().create(name));
}
```

- [ ] **Step 3: Commit C-API Extensions**
```bash
git add src/core/catchem_api.hpp src/core/catchem_api.cpp
git commit -m "feat(core): expose catchem_core_add_process_by_name to bind dynamically registered C++ processes"
```

---

### Task 5: Update Interop Tests and Run Full Verification Suite

**Files:**
- Modify: `tests/test_catchem_interop.cpp`

**Interfaces:**
- Consumes: `catchem_core_add_process_by_name`
- Produces: Green test suite in Docker container.

- [ ] **Step 1: Modify `tests/test_catchem_interop.cpp`**
Add dynamic registration verification section (e.g. Test 4):
```cpp
        // ==========================================
        // TEST 4: Modular dynamic process registration validation
        // ==========================================
        {
            // Simulate Fortran explicitly linking and calling C++ register_settling_cpp
            extern "C" void catchem_register_settling_cpp();
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
```

- [ ] **Step 2: Run verification in Docker cece-dev**
Confirm the library compiles and the new test execution succeeds.
```bash
mkdir -p build && cd build && cmake -DENABLE_KOKKOS=ON .. && make -j$(nproc) && ./tests/test_catchem_interop
```

- [ ] **Step 3: Commit and finalize**
```bash
git add tests/test_catchem_interop.cpp
git commit -m "test(interop): verify modular dynamic C-API process registration in test suite"
```
