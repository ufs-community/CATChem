# Creating a New Column-Based Process

This guide provides a step-by-step tutorial for creating a new atmospheric process in CATChem using the modern C++ Core and Fortran ScienceBridge architecture. We will use a simplified version of the `seasalt` process as an example.

## 1. Directory and File Structure

Create a directory for your new process within `src/process`. The modern package structure looks like this:

```
src/process/newprocess/
├── catchem_process_newprocess.hpp    # C++ Process class header (extends catchem::ProcessInterface)
├── catchem_process_newprocess.cpp    # C++ Process class implementation & ProcessRegistry registration
├── NewProcessScienceBridge.F90      # BIND(C) Fortran Science Bridge
├── NewProcessCommon_Mod.F90         # Process configuration and common data structures
├── schemes/
│   ├── NewProcessScheme_DEFAULT_Mod.F90 # Pure Fortran science scheme implementation
│   └── CMakeLists.txt               # Schemes CMake target
└── CMakeLists.txt                   # Process library CMake configuration
```

-   **`catchem_process_newprocess.hpp / .cpp`**: C++ process class wrapper inheriting from `catchem::ProcessInterface`. Manages process lifecycle (`init`, `run`, `finalize`), state binding, and registers with `catchem::ProcessRegistry`.
-   **`NewProcessScienceBridge.F90`**: C-interoperable Fortran `BIND(C)` module that unpacks Kokkos host view pointers into Fortran arrays, handles column iteration, and delegates to pure Fortran science schemes.
-   **`NewProcessCommon_Mod.F90`**: Fortran configuration module for scheme parameters.
-   **`schemes/`**: Pure Fortran science algorithms with no infrastructure or C-API dependencies.

---

## 2. The C++ Process Class

The C++ process class inherits from `catchem::ProcessInterface`.

`catchem_process_newprocess.hpp`:
```cpp
#pragma once
#include "catchem_process_interface.hpp"
#include <string>
#include <vector>

namespace catchem {

    class NewProcessProcess : public ProcessInterface {
    public:
        std::string active_scheme;
        bool diagnostics_enabled;

        NewProcessProcess();

        std::string get_name() const override { return "newprocess"; }
        void init(std::shared_ptr<StateManager> state) override;
        void run(std::shared_ptr<StateManager> state) override;
        void finalize() override;
    };

} // namespace catchem
```

`catchem_process_newprocess.cpp`:
```cpp
#include "catchem_process_newprocess.hpp"
#include "catchem_diagnostic_manager.hpp"
#include "catchem_process_registry.hpp"
#include <iostream>

extern "C" {
void run_newprocess_science_bridge(
    int n_cols, int n_levels, int n_species, int n_soil, double dt,
    const char* active_scheme, int diagnostics,
    double* u10m, double* v10m, double* sst,
    double* species_density, double* species_radius,
    double* species_lower_radius, double* species_upper_radius,
    double* conc, double* tendency,
    double* diag_emission_total,
    const int* diagnostic_species_id, int n_diag_species
);
}

namespace catchem {

    NewProcessProcess::NewProcessProcess()
        : active_scheme("default"), diagnostics_enabled(true) {}

    void NewProcessProcess::init(std::shared_ptr<StateManager> state) {
        // Register diagnostic fields with C++ DiagnosticManager
        std::vector<int> dims_2d = {state->n_cols, 1};
        state->diag_mgr->register_field("newprocess_emission_total", "Total emission flux",
                                         "kg/m2/s", DiagType::FIELD_2D, dims_2d);
    }

    void NewProcessProcess::run(std::shared_ptr<StateManager> state) {
        state->sync_to_host();

        // Retrieve Meteorological state pointers from StateManager
        double* u10m_ptr = state->met.fields_2d["U10M"] ? state->met.fields_2d["U10M"]->host_view.data() : nullptr;
        double* v10m_ptr = state->met.fields_2d["V10M"] ? state->met.fields_2d["V10M"]->host_view.data() : nullptr;
        double* sst_ptr = state->met.fields_2d["SST"] ? state->met.fields_2d["SST"]->host_view.data() : nullptr;

        double* conc_ptr = state->chem.conc ? state->chem.conc->host_view.data() : nullptr;
        std::vector<double> local_tendency(state->n_cols * state->n_levels * state->n_species, 0.0);

        // Invoke Fortran Science Bridge
        run_newprocess_science_bridge(
            state->n_cols, state->n_levels, state->n_species, 4, state->time.timestep,
            active_scheme.c_str(), diagnostics_enabled ? 1 : 0,
            u10m_ptr, v10m_ptr, sst_ptr,
            nullptr, nullptr, nullptr, nullptr,
            conc_ptr, local_tendency.data(),
            nullptr,
            nullptr, 0
        );

        // Update state
        if (conc_ptr) {
            for (size_t i = 0; i < local_tendency.size(); ++i) {
                conc_ptr[i] += state->time.timestep * local_tendency[i];
            }
        }

        state->sync_to_device();
    }

    void NewProcessProcess::finalize() {}

    CATCHEM_REGISTER_PROCESS(NewProcessProcess, "newprocess")

} // namespace catchem

extern "C" {
void catchem_register_newprocess_cpp() {
    catchem::ProcessRegistry::get_instance().register_process(
        "newprocess", []() { return std::make_shared<catchem::NewProcessProcess>(); }
    );
}
}
```

---

## 3. The Fortran Science Bridge

The science bridge provides C-interoperability (`BIND(C)`) and delegates column calculations to pure Fortran schemes.

`NewProcessScienceBridge.F90`:
```fortran
module NewProcessScienceBridge_Mod
   use iso_c_binding
   use precision_mod, only: fp
   use NewProcessCommon_Mod
   use NewProcessScheme_DEFAULT_Mod, only: compute_default

   implicit none
   private

   public :: run_newprocess_science_bridge

contains

   subroutine run_newprocess_science_bridge( &
      n_cols, n_levels, n_species, n_soil, dt, &
      active_scheme, diagnostics, &
      u10m, v10m, sst, &
      species_density, species_radius, species_lower_radius, species_upper_radius, &
      conc, tendency, &
      diag_emission_total, &
      diagnostic_species_id, n_diag_species &
      ) bind(C, name="run_newprocess_science_bridge")

      integer(c_int), value :: n_cols, n_levels, n_species, n_soil
      real(c_double), value :: dt
      character(kind=c_char), intent(in) :: active_scheme(*)
      integer(c_int), value :: diagnostics

      type(c_ptr), value :: u10m, v10m, sst
      type(c_ptr), value :: species_density, species_radius, species_lower_radius, species_upper_radius
      type(c_ptr), value :: conc, tendency, diag_emission_total
      integer(c_int), value :: n_diag_species
      integer(c_int), intent(in) :: diagnostic_species_id(n_diag_species)

      real(c_double), pointer :: f_u10m(:), f_v10m(:), f_sst(:)
      real(c_double), pointer :: f_conc(:,:,:), f_tendency(:,:,:)

      integer :: icol

      if (c_associated(u10m)) call c_f_pointer(u10m, f_u10m, [n_cols])
      if (c_associated(v10m)) call c_f_pointer(v10m, f_v10m, [n_cols])
      if (c_associated(sst)) call c_f_pointer(sst, f_sst, [n_cols])
      if (c_associated(conc)) call c_f_pointer(conc, f_conc, [n_cols, n_levels, n_species])
      if (c_associated(tendency)) call c_f_pointer(tendency, f_tendency, [n_cols, n_levels, n_species])

      do icol = 1, n_cols
         call compute_default( &
            n_levels, n_species, real(dt, fp), &
            f_u10m(icol), f_v10m(icol), f_sst(icol), &
            f_conc(icol, :, :), f_tendency(icol, :, :) &
         )
      end do

   end subroutine run_newprocess_science_bridge

end module NewProcessScienceBridge_Mod
```

---

## 4. Pure Science Scheme

The scheme contains pure Fortran algorithms with no C-API or framework dependencies.

`schemes/NewProcessScheme_DEFAULT_Mod.F90`:
```fortran
module NewProcessScheme_DEFAULT_Mod
   use precision_mod, only: fp
   implicit none
   private
   public :: compute_default

contains

   subroutine compute_default(n_levels, n_species, dt, u10m, v10m, sst, conc, tendency)
      integer, intent(in) :: n_levels, n_species
      real(fp), intent(in) :: dt, u10m, v10m, sst
      real(fp), intent(in) :: conc(n_levels, n_species)
      real(fp), intent(out) :: tendency(n_levels, n_species)

      integer :: k, ispec

      do ispec = 1, n_species
         do k = 1, n_levels
            tendency(k, ispec) = (u10m**2 + v10m**2) * sst * 1.0e-12_fp
         end do
      end do
   end subroutine compute_default

end module NewProcessScheme_DEFAULT_Mod
```

---

## 5. CMake Integration

`src/process/newprocess/CMakeLists.txt`:
```cmake
set(
  NEWPROCESS_PROCESS_SOURCES
  NewProcessCommon_Mod.F90
  NewProcessScienceBridge.F90
)

set(
  NEWPROCESS_SCHEME_SOURCES
  schemes/NewProcessScheme_DEFAULT_Mod.F90
)

set(NEWPROCESS_ALL_SOURCES ${NEWPROCESS_PROCESS_SOURCES} ${NEWPROCESS_SCHEME_SOURCES})

if(CATCHEM_ENABLE_KOKKOS)
  list(APPEND NEWPROCESS_ALL_SOURCES catchem_process_newprocess.cpp)
endif()

set(_lib CATChem_process_newprocess)
add_library(${_lib} ${NEWPROCESS_ALL_SOURCES})

target_link_libraries(${_lib} PUBLIC CATChem_core)

if(CATCHEM_ENABLE_KOKKOS)
  target_link_libraries(${_lib} PUBLIC Kokkos::kokkos CATChem_core_cpp)
  target_compile_definitions(${_lib} PRIVATE CATCHEM_ENABLE_KOKKOS)
endif()

set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
```

---

## 6. Testing

Generated process packages automatically include a science unit test (`tests/test_newprocess_science.f90`) and can also be tested using the C++ property test harness in `tests/test_catchem_properties.cpp`.
