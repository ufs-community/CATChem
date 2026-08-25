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
  PUBLIC CATChem_core_cpp CATChem_core musica yaml-cpp::yaml-cpp Kokkos::kokkos
)

# Robust target-level linking to propagate fetched include directories automatically
if(TARGET musica::micm)
  target_link_libraries(catchem_process_gaschem PUBLIC musica::micm)
elseif(TARGET micm::micm)
  target_link_libraries(catchem_process_gaschem PUBLIC micm::micm)
elseif(TARGET micm)
  target_link_libraries(catchem_process_gaschem PUBLIC micm)
endif()

if(TARGET musica::mechanism_configuration)
  target_link_libraries(
    catchem_process_gaschem
    PUBLIC musica::mechanism_configuration
  )
elseif(TARGET mechanism_configuration)
  target_link_libraries(catchem_process_gaschem PUBLIC mechanism_configuration)
endif()
