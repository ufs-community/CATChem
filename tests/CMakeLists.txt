foreach(opt IN LISTS _catchem_compiler_options)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${opt}>")
endforeach()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/catchem_test_config.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/catchem_test_config.hpp"
)

include_directories("${CMAKE_CURRENT_BINARY_DIR}")

# Automatically link the musica configs folder so that CTest can find relative NetCDF files on GHA CI
if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/configs")
  file(
    CREATE_LINK
      "${PROJECT_SOURCE_DIR}/src/external/musica/configs"
      "${CMAKE_CURRENT_BINARY_DIR}/configs"
    SYMBOLIC
  )
endif()

# Restore compilation of foundational Fortran testing helper library
add_library(testing testing_mod.f90)
target_link_libraries(testing PRIVATE CATChem_core)
set_target_properties(
  testing
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
target_include_directories(testing PRIVATE ${CMAKE_BINARY_DIR}/include)

# Restored Low-Level Foundational Unit Tests
macro(add_fortran_test name_test)
  add_executable(${name_test} ${name_test}.f90)
  target_link_libraries(${name_test} PRIVATE CATChem_core testing)
  set_target_properties(
    ${name_test}
    PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
      LINKER_LANGUAGE Fortran
  )
  target_include_directories(${name_test} PRIVATE ${CMAKE_BINARY_DIR}/include)
  add_dependencies(${name_test} CATChem_core)
  add_test(
    NAME ${name_test}
    COMMAND ${name_test}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endmacro()

add_fortran_test(test_Precision)
add_fortran_test(test_Error)
add_fortran_test(test_Constants)
add_fortran_test(test_GridGeometry)
add_fortran_test(test_TimeState)
add_fortran_test(test_MetState)
add_fortran_test(test_ChemSpeciesUtils)
add_fortran_test(test_DiagnosticInterface)

# Process Science Unit Tests
macro(add_process_science_test name_test process_lib)
  add_executable(${name_test} ${name_test}.f90)
  target_link_libraries(
    ${name_test}
    PRIVATE ${process_lib} CATChem_core testing
  )
  set_target_properties(
    ${name_test}
    PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
      LINKER_LANGUAGE Fortran
  )
  target_include_directories(${name_test} PRIVATE ${CMAKE_BINARY_DIR}/include)
  add_dependencies(${name_test} ${process_lib} CATChem_core)
  add_test(
    NAME ${name_test}
    COMMAND ${name_test}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endmacro()

add_process_science_test(test_settling_science CATChem_process_settling)
add_process_science_test(test_seasalt_science CATChem_process_seasalt)
add_process_science_test(test_dust_science CATChem_process_dust)
add_process_science_test(test_drydep_science CATChem_process_drydep)
add_process_science_test(test_wetdep_science CATChem_process_wetdep)
add_process_science_test(test_so4chem_science CATChem_process_so4chem)
add_process_science_test(test_carbchem_science CATChem_process_carbchem)

# UnitConversion needs specific include search dirs
add_executable(test_UnitConversion test_UnitConversion.f90)
target_link_libraries(test_UnitConversion PRIVATE CATChem_core testing)
set_target_properties(
  test_UnitConversion
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
target_include_directories(
  test_UnitConversion
  PRIVATE ${CMAKE_BINARY_DIR}/include
)
add_dependencies(test_UnitConversion CATChem_core)
add_test(
  NAME test_UnitConversion
  COMMAND test_UnitConversion
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# New C++ Interop and Integration Suite
if(CATCHEM_ENABLE_KOKKOS)
  add_executable(test_catchem_logger test_catchem_logger.cpp)
  target_link_libraries(
    test_catchem_logger
    PRIVATE CATChem_core_cpp Kokkos::kokkos
  )
  add_test(NAME test_catchem_logger COMMAND test_catchem_logger)

  add_executable(
    test_catchem_process_registry
    test_catchem_process_registry.cpp
  )
  target_link_libraries(
    test_catchem_process_registry
    PRIVATE
      CATChem_core_cpp
      CATChem_core
      Kokkos::kokkos
      CATChem_process_settling
      CATChem_process_seasalt
      CATChem_process_drydep
      CATChem_process_wetdep
      CATChem_process_so4chem
      CATChem_process_dust
      CATChem_process_carbchem
  )
  set_target_properties(
    test_catchem_process_registry
    PROPERTIES LINKER_LANGUAGE CXX
  )
  add_test(
    NAME test_catchem_process_registry
    COMMAND test_catchem_process_registry
  )

  add_executable(test_catchem_interop test_catchem_interop.cpp)
  target_link_libraries(
    test_catchem_interop
    PRIVATE
      CATChem_core_cpp
      CATChem_core
      Kokkos::kokkos
      CATChem_process_settling
      CATChem_process_seasalt
      CATChem_process_drydep
      CATChem_process_wetdep
      CATChem_process_so4chem
      CATChem_process_dust
      CATChem_process_carbchem
  )
  set_target_properties(test_catchem_interop PROPERTIES LINKER_LANGUAGE CXX)
  add_test(NAME test_catchem_interop COMMAND test_catchem_interop)

  # Modern Randomized Property-Based Testing
  add_executable(test_catchem_properties test_catchem_properties.cpp)
  target_link_libraries(
    test_catchem_properties
    PRIVATE
      CATChem_core_cpp
      CATChem_core
      Kokkos::kokkos
      CATChem_process_settling
      CATChem_process_seasalt
      CATChem_process_drydep
      CATChem_process_wetdep
      CATChem_process_so4chem
      CATChem_process_dust
      CATChem_process_carbchem
  )
  set_target_properties(test_catchem_properties PROPERTIES LINKER_LANGUAGE CXX)
  add_test(NAME test_catchem_properties COMMAND test_catchem_properties)

  if(CATCHEM_BUILD_MUSICA)
    # Photolysis Process Unit/Integration Test
    add_executable(test_catchem_photolysis test_catchem_photolysis.cpp)
    target_link_libraries(
      test_catchem_photolysis
      PRIVATE
        CATChem_core_cpp
        CATChem_core
        Kokkos::kokkos
        catchem_process_photolysis
        yaml-cpp::yaml-cpp
    )
    set_target_properties(
      test_catchem_photolysis
      PROPERTIES LINKER_LANGUAGE CXX
    )
    add_test(NAME test_catchem_photolysis COMMAND test_catchem_photolysis)

    # Coupled Process Integration Test
    add_executable(test_catchem_gaschem test_catchem_gaschem.cpp)
    target_link_libraries(
      test_catchem_gaschem
      PUBLIC
        CATChem_core_cpp
        CATChem_core
        catchem_process_photolysis
        catchem_process_gaschem
        Kokkos::kokkos
    )
    set_target_properties(test_catchem_gaschem PROPERTIES LINKER_LANGUAGE CXX)
    add_test(NAME test_catchem_gaschem COMMAND test_catchem_gaschem)
  endif()
endif()

add_executable(test_catchem_gaschem_units test_catchem_gaschem_units.cpp)
add_test(NAME test_catchem_gaschem_units COMMAND test_catchem_gaschem_units)
