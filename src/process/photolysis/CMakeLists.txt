# src/process/photolysis/CMakeLists.txt
add_library(catchem_process_photolysis STATIC catchem_process_photolysis.cpp)

target_include_directories(
  catchem_process_photolysis
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src/core
    ${PROJECT_SOURCE_DIR}/src/external/musica/include
)

target_link_libraries(
  catchem_process_photolysis
  PUBLIC CATChem_core_cpp CATChem_core musica yaml-cpp::yaml-cpp Kokkos::kokkos
)
