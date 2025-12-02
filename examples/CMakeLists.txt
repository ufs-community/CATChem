cmake_minimum_required(VERSION 3.12)
project(CATChemNetCDFMPIExample LANGUAGES Fortran)

# Find required packages
find_package(PkgConfig REQUIRED)

# Find NetCDF
pkg_check_modules(NETCDF REQUIRED netcdf-fortran)

# Find MPI (optional)
find_package(MPI COMPONENTS Fortran)

# Set compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O2")

# Add preprocessor definitions
add_definitions(-DNETCDF_ENABLED)

if(MPI_FOUND)
  add_definitions(-DMPI_ENABLED)
  message(STATUS "MPI support enabled")
else()
  message(STATUS "MPI support disabled (MPI not found)")
endif()

# Check for parallel NetCDF
if(NETCDF_PARALLEL)
  add_definitions(-DNETCDF_PARALLEL)
  message(STATUS "Parallel NetCDF support enabled")
else()
  message(STATUS "Parallel NetCDF support disabled")
endif()

# Include directories
include_directories(${NETCDF_INCLUDE_DIRS})
include_directories(../src/core)
include_directories(../src/api)

if(MPI_FOUND)
  include_directories(${MPI_Fortran_INCLUDE_DIRS})
endif()

# Add executable
add_executable(
  netcdf_mpi_example
  netcdf_mpi_example.F90
  ../src/core/precision_mod.F90
  ../src/core/error_mod.F90
  ../src/api/CATChemNetCDF_Mod.F90
)

# Link libraries
target_link_libraries(netcdf_mpi_example ${NETCDF_LIBRARIES})

if(MPI_FOUND)
  target_link_libraries(netcdf_mpi_example ${MPI_Fortran_LIBRARIES})
  target_compile_options(
    netcdf_mpi_example
    PRIVATE ${MPI_Fortran_COMPILE_FLAGS}
  )
endif()

# Add custom target to run the example
add_custom_target(
  run_netcdf_mpi_example
  COMMAND ./netcdf_mpi_example
  DEPENDS netcdf_mpi_example
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Running NetCDF MPI compatibility example"
)

# Add MPI run target if MPI is available
if(MPI_FOUND)
  add_custom_target(
    run_netcdf_mpi_example_parallel
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ./netcdf_mpi_example
    DEPENDS netcdf_mpi_example
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running NetCDF MPI compatibility example with 4 MPI ranks"
  )
endif()
