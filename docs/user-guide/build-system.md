# Build System

CATChem uses CMake as its primary build system, providing a flexible and robust framework for building the model across different platforms and compilers.

## Overview

The CATChem build system is designed with the following principles:

- **Modern CMake**: Uses CMake 3.18+ features for better maintainability
- **Modular Structure**: Each process and component has its own CMakeLists.txt
- **Cross-Platform**: Supports Linux, macOS, and Windows (via WSL)
- **Dependency Management**: Automatic detection and configuration of external libraries
- **Testing Integration**: Built-in support for unit tests and validation

## Build Requirements

### System Requirements

- **CMake**: Version 3.18 or later
- **Fortran Compiler**:
  - GNU Fortran (gfortran) 9.0 or later (recommended)
  - Intel Fortran Compiler 2021.1 or later
  - NVIDIA HPC SDK 21.3 or later

### Required Dependencies

- **NetCDF-Fortran**: For I/O operations
- **MPI**: For parallel execution (optional but recommended)
- **YAML-Fortran**: For configuration file parsing

### Optional Dependencies

- **HDF5**: Required for NetCDF-4 support
- **PnetCDF**: For parallel NetCDF I/O
- **ESMF**: For Earth System Modeling Framework integration

## Quick Build

```bash
# Clone the repository
git clone https://github.com/UFS-Community/CATChem.git
cd CATChem

# Create build directory
mkdir build && cd build

# Configure and build
cmake ..
make -j4

# Run a test case
make test
```

## Detailed Build Configuration

### Basic Configuration

```bash
# Create and enter build directory
mkdir build && cd build

# Configure with default options
cmake ..
```

### Advanced Configuration Options

```bash
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DENABLE_MPI=ON \
    -DENABLE_OPENMP=ON \
    -DENABLE_TESTING=ON
```

#### Common CMake Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `CMAKE_BUILD_TYPE` | Build type (Debug, Release, RelWithDebInfo) | `Release` |
| `CMAKE_Fortran_COMPILER` | Fortran compiler path | Auto-detected |
| `CMAKE_INSTALL_PREFIX` | Installation directory | `/usr/local` |
| `ENABLE_MPI` | Enable MPI parallelization | `ON` |
| `ENABLE_OPENMP` | Enable OpenMP threading | `ON` |
| `ENABLE_TESTING` | Build unit tests | `ON` |
| `ENABLE_DOCS` | Build documentation | `OFF` |

#### Dependency Configuration

```bash
# Specify NetCDF location
cmake .. -DNetCDF_ROOT=/path/to/netcdf

# Specify MPI implementation
cmake .. -DMPI_Fortran_COMPILER=mpif90

# Use custom YAML-Fortran
cmake .. -DYAML_ROOT=/path/to/yaml-fortran
```

## Build Structure

The CATChem source tree is organized as follows:

```
cc_restructure/
├── CMakeLists.txt              # Root CMake configuration
├── cmake/                      # CMake modules and utilities
│   └── FindNetCDF.cmake       # NetCDF detection module
├── src/                        # Source code
│   ├── CMakeLists.txt         # Core library configuration
│   ├── api/                   # Public API
│   ├── core/                  # Core modules
│   └── process/               # Process modules
│       ├── CMakeLists.txt     # Process library configuration
│       └── settling/          # Example process
│           └── CMakeLists.txt # Process-specific build
├── drivers/                    # Model drivers
│   ├── ccpp/                  # CCPP integration
│   └── nuopc/                 # NUOPC integration
└── tests/                      # Test suite
    └── CMakeLists.txt         # Test configuration
```

## Process Module Integration

Each process module follows a standard CMakeLists.txt structure:

```cmake
# src/process/settling/CMakeLists.txt
set(PROCESS_NAME settling)

# Source files
set(PROCESS_SOURCES
    settlingProcess_Mod.F90
    settlingCommon_Mod.F90
)

# Create process library
add_library(${PROCESS_NAME}_process ${PROCESS_SOURCES})

# Link dependencies
target_link_libraries(${PROCESS_NAME}_process
    PRIVATE
        catchem_core
        catchem_process_interface
)

# Include schemes subdirectory
add_subdirectory(schemes)
```

### Adding New Processes

To add a new process to the build system:

1. Create process directory structure:
   ```bash
   mkdir src/process/myprocess
   mkdir src/process/myprocess/schemes
   ```

2. Create `src/process/myprocess/CMakeLists.txt`:
   ```cmake
   set(PROCESS_NAME myprocess)

   set(PROCESS_SOURCES
       myprocessProcess_Mod.F90
       myprocessCommon_Mod.F90
   )

   add_library(${PROCESS_NAME}_process ${PROCESS_SOURCES})
   target_link_libraries(${PROCESS_NAME}_process
       PRIVATE catchem_core catchem_process_interface
   )

   add_subdirectory(schemes)
   ```

3. Add the new process to `src/process/CMakeLists.txt`:
   ```cmake
   # Add new process
   add_subdirectory(myprocess)
   ```

## Compiler-Specific Configuration

### GNU Fortran (gfortran)

```bash
cmake .. \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_Fortran_FLAGS="-O3 -ffast-math -funroll-loops"
```

### Intel Fortran

```bash
cmake .. \
    -DCMAKE_Fortran_COMPILER=ifort \
    -DCMAKE_Fortran_FLAGS="-O3 -xHost -ipo"
```

### NVIDIA HPC SDK

```bash
cmake .. \
    -DCMAKE_Fortran_COMPILER=nvfortran \
    -DCMAKE_Fortran_FLAGS="-O3 -fast"
```

## Platform-Specific Notes

### HPC Systems

For HPC systems, use the appropriate compiler wrappers:

```bash
# Example for systems with modules
module load cmake netcdf-fortran mpi

cmake .. \
    -DCMAKE_Fortran_COMPILER=ftn \
    -DENABLE_MPI=ON
```

### macOS with Homebrew

```bash
# Install dependencies
brew install cmake netcdf gfortran

# Configure build
cmake .. \
    -DCMAKE_Fortran_COMPILER=gfortran-12 \
    -DNetCDF_ROOT=$(brew --prefix netcdf)
```

## Testing and Validation

### Running Tests

```bash
# Build and run all tests
make test

# Run specific test categories
ctest -L unit      # Unit tests only
ctest -L integration  # Integration tests only
ctest -L validation   # Validation tests only
```

### Debug Builds

For debugging and development:

```bash
cmake .. \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_Fortran_FLAGS="-g -O0 -fcheck=all -Wall"

make -j4
```

## Installation

```bash
# Install to configured prefix
make install

# Create package
make package
```

## Troubleshooting

### Common Issues

1. **NetCDF not found**:
   ```bash
   cmake .. -DNetCDF_ROOT=/path/to/netcdf
   ```

2. **MPI compiler wrapper issues**:
   ```bash
   cmake .. -DMPI_Fortran_COMPILER=mpif90
   ```

3. **Fortran module path issues**:
   ```bash
   # Clean build directory
   rm -rf build/*
   cmake ..
   ```

### Getting Build Information

```bash
# View detailed build configuration
cmake .. -LAH

# Show compile commands
make VERBOSE=1
```

## Best Practices

1. **Always use out-of-source builds**
2. **Clean build directory when changing major configuration**
3. **Use compiler-specific optimization flags for production builds**
4. **Enable testing during development**
5. **Use appropriate build types (Debug vs Release)**

For more advanced build configuration and platform-specific instructions, see the [HPC Installation Guide](hpc-installation.md).
