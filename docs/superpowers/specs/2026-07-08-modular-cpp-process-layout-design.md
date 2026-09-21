# Spec: Modular C++ Process Layout and Linker-Safe Registration

*   **Date:** July 8, 2026
*   **Status:** APPROVED
*   **Author:** Gemini CLI

---

## 1. Executive Summary
The CATChem C++ core architecture is designed to orchestrate physical and chemical processes using a centralized, high-performance C++20 layer. Initially, concrete process implementations (e.g. `SettlingProcess`, `SeaSaltProcess`, etc.) were placed under `@src/core`. To achieve absolute separation of concerns, improve modularity, and conform to the project's original directories structure, we are re-architecting the system to place all physical processes under `@src/process/<name>/**`.

This design establishes **mixed-language (Fortran/C++) process libraries** under `@src/process/` and introduces a **linker-safe, ISO_C_BINDING registration pattern** that prevents C++ symbols from being optimized out by the linker when compiled inside process-specific static libraries.

---

## 2. Directory Layout & Source Relocation
All concrete C++ physical process controllers will be relocated from `@src/core/` to their respective subdirectories within `@src/process/`.

### 2.1 File Map

| Process | Source File in `src/core/` (Old) | New Location in `src/process/` (New) |
|---|---|---|
| **Settling** | `catchem_process_settling.hpp` / `.cpp` | `src/process/settling/catchem_process_settling.hpp` / `.cpp` |
| **SeaSalt** | `catchem_process_seasalt.hpp` / `.cpp` | `src/process/seasalt/catchem_process_seasalt.hpp` / `.cpp` |
| **DryDep** | `catchem_process_drydep.hpp` / `.cpp` | `src/process/drydep/catchem_process_drydep.hpp` / `.cpp` |
| **WetDep** | `catchem_process_wetdep.hpp` / `.cpp` | `src/process/wetdep/catchem_process_wetdep.hpp` / `.cpp` |
| **SO4chem** | `catchem_process_so4chem.hpp` / `.cpp` | `src/process/so4chem/catchem_process_so4chem.hpp` / `.cpp` |

---

## 3. CMake Target Reconfiguration
Each process directory builds its own static library (e.g., `CATChem_process_settling`). These libraries will be converted to mixed-language (Fortran + C++) targets when `ENABLE_KOKKOS` is active.

### 3.1 Template `CMakeLists.txt` Integration
In each of the five process `CMakeLists.txt` files (e.g. `src/process/settling/CMakeLists.txt`), C++ source compilation and linking will be configured as follows:

```cmake
# Define settling process sources
set(
  SETTLING_PROCESS_SOURCES
  SettlingCommon_Mod.F90
  ProcessSettlingInterface_Mod.F90
  SettlingProcessCreator_Mod.F90
)

# Define settling scheme sources
set(
  SETTLING_SCHEME_SOURCES
  schemes/SettlingScheme_GOCART_Mod.F90
  schemes/SettlingPhysics_Mod.F90
)

# Combine all sources
set(SETTLING_ALL_SOURCES ${SETTLING_PROCESS_SOURCES} ${SETTLING_SCHEME_SOURCES})

# Conditional C++ component adding
if(ENABLE_KOKKOS)
  list(APPEND SETTLING_ALL_SOURCES catchem_process_settling.cpp)
endif()

# Create the settling process library
set(_lib CATChem_process_settling)
add_library(${_lib} ${SETTLING_ALL_SOURCES})

# Link with required libraries
target_link_libraries(${_lib} PUBLIC CATChem_core)

# Conditional Kokkos linking
if(ENABLE_KOKKOS)
  target_link_libraries(${_lib} PUBLIC Kokkos::kokkos)
  target_compile_definitions(${_lib} PRIVATE ENABLE_KOKKOS)
endif()
```

---

## 4. Linker-Safe Process Registration
When concrete process implementations are compiled into static libraries, linker optimizations ("dead-code stripping") often remove unused object files if they are not explicitly referenced from `main`. To guarantee that C++ creator lambdas are always registered in the dynamic `ProcessRegistry` at initialization, we implement an explicit, linker-safe callback mechanism across the C/Fortran boundary.

### 4.1 C-Linkage Registration Routine
Each C++ process controller exports a standard C-linkage registration routine:

```cpp
// src/process/settling/catchem_process_settling.cpp
#include "catchem_process_settling.hpp"
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

### 4.2 Fortran Invocation
The process Fortran interface (e.g., `ProcessSettlingInterface_Mod.F90` under `process_init`) declares and calls this registration routine. Because the Fortran module is always linked by the driver, this explicitly forces the linker to preserve and bind the associated C++ object files and register the C++ class:

```fortran
! src/process/settling/ProcessSettlingInterface_Mod.F90
subroutine process_init(this, container, rc)
   ...
#ifdef ENABLE_KOKKOS
   interface
      subroutine catchem_register_settling_cpp() bind(c, name="catchem_register_settling_cpp")
      end subroutine catchem_register_settling_cpp
   end interface

   ! Register C++ dynamic process in the global ProcessRegistry
   call catchem_register_settling_cpp()
#endif
   ...
end subroutine process_init
```

---

## 5. C-API Core Orchestration Extension
We extend the C++ Core and C-API layers with standard methods for process instantiation by string name.

### 5.1 C-API Header additions (`src/core/catchem_api.hpp`):
```cpp
void catchem_core_add_process_by_name(void* core_ptr, const char* name);
```

### 5.2 C-API Implementation additions (`src/core/catchem_api.cpp`):
```cpp
void catchem_core_add_process_by_name(void* core_ptr, const char* name) {
    auto* core = static_cast<catchem::Core*>(core_ptr);
    core->add_process(catchem::ProcessRegistry::get_instance().create(name));
}
```

---

## 6. Testing & Validation Plan
*   **Verification Target:** `tests/test_catchem_interop.cpp` will be updated to:
    1. Register mock or concrete physical processes via their registration callbacks.
    2. Instantiate processes using `catchem_core_add_process_by_name`.
    3. Run simulation cycles verifying exact state updates, diagnostic updates, and execution correctness.
*   **Compilation Check:** Compile and run inside the target Docker container, checking that there are no linker errors or unused symbol stripped-out behavior.
