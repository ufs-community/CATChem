# CATChem

**CATChem** (Configurable ATmospheric Chemistry) is a modelling component that includes all chemical and aerosol processes needed to perform atmospheric chemistry and composition simulations within a model through a flexible, easy to modify, and well-documented infrastructure.

## Features

- **NUOPC-compliant interface** for integration with Earth system models
- **ESMF I/O integration** for efficient parallel file operations
- **Flexible configuration system** with YAML-based setup
- **Automatic regridding** with weight file management
- **CF-compliant input/output** supporting Climate and Forecast conventions
- **Comprehensive utility tools** for configuration and weight generation

## Documentation

Comprehensive documentation is available in the `docs/` directory:

- **[Documentation Home](docs/index.md)** - Complete documentation portal
- **[Quick Start Guide](docs/quick-start/index.md)** - Get started quickly
- **[User Guide](docs/user-guide/index.md)** - Complete user documentation
- **[Developer Guide](docs/developer-guide/index.md)** - Technical documentation for developers
- **[API Reference](docs/api/index.md)** - Complete API documentation

## Quick Start

### Container

The official UFS-Chem and CATChem images are recommended for developers and interested users.

1. Base image for developers:

    ```bash
    # Mounts the CATChem root directory to /opt/project inside the container
    docker run -it --rm --platform linux/amd64 -v <CATChem root directory>:/opt/project noaaepic/ufschem-spack-base-ubuntu-gcc-13:latest
    # Now inside the container...
    cd /opt/project # <-- changes made in this directory will be reflected in the local file system
    ```

2. Running a pre-installed and tested CATChem:

    ```bash
    # Mounts the CATChem root directory to /opt/project inside the container
    docker run -it --rm --platform linux/amd64 -v <CATChem root directory>:/opt/project noaaepic/catchem-ubuntu-gcc-13:latest
    # Now inside the container...
    cd /opt/catchem
    ```

### NUOPC Interface

```bash
# Build the NUOPC interface
cd drivers/nuopc
make

# Set up utilities and weight files
../../util/generate_esmf_weights.sh setup

# Validate configuration
python ../../util/validate_weight_config.py catchem_input_config.yml

# Run standalone test
./catchem_nuopc_driver
```

### Utility Scripts

The `util/` directory contains helpful scripts:

- **`generate_esmf_weights.sh`** - ESMF weight file generation and management
- **`manage_weights.sh`** - Weight file optimization and maintenance
- **`validate_weight_config.py`** - Configuration validation and checking

## Development

### pre-commit

If you don't have the command-line tool `pre-commit` available, [install it](https://pre-commit.com/#install).

Install the pre-commit hooks with

```
pre-commit install --install-hooks
```

Now, some checks and auto-formatting will run automatically when you commit.

> [!NOTE]
> For source files with their own formatting that we don't intend to modify (or only modify slightly),
> e.g. "vendored" modules, we generally don't want `findent` to be applied.
> For such files, update the `exclude` section in `.pre-commit-config.yaml` accordingly.

> [!TIP]
> ```
> pre-commit run --all-files
> ```
> can be used to check the pre-commit config
> (e.g. to make sure that the `findent` `exclude` section is working as you intended)
> or to catch up if commits were made without pre-commit or the pre-commit config was updated.

### Build

To test the build locally, first configure, using `FC` to specify the compiler you want to use, then build. For example:

```
FC=gfortran-12 cmake -B build
```
```
cmake --build build -j
```

To clean, you can use

```
cmake --build build --target clean
```

or remove the build directory (`./build`).

### Test

To run the tests, after building, use

```
ctest --test-dir build/tests
```

There [are options](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html#testing-using-ctest) for selecting specific tests.

Edit `tests/CMakelists.txt` to add new tests.
