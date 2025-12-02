# CATChem Process Generator

The CATChem Process Generator is a powerful tool that automates the creation of new atmospheric processes following established architecture patterns. This guide covers everything from basic usage to advanced customization.

## Overview

The process generator is a command-line tool that creates all the necessary boilerplate for a new process, including:

-   Fortran source files for the process interface, schemes, and common utilities.
-   CMake build files.
-   Unit tests.
-   Comprehensive documentation.

By using the generator, you can focus on the scientific implementation of your process, while ensuring that it adheres to CATChem's architectural standards.

## Quick Start

The simplest way to create a new process is to use a YAML configuration file.

1.  **Create a configuration file** (`my_process.yaml`):

    ```yaml
    process:
      name: MyProcess
      type: transformation
      description: "A new example process."
      schemes: [scheme1, scheme2]
      species: [O3, NO2]
    ```

2.  **Run the generator:**

    ```bash
    cd /path/to/catchem/tools/process_generator
    python process_generator.py generate --config configs/my_process.yaml
    ```

This will generate a complete process structure in `src/process/myprocess`, along with tests and documentation.

## Command-Line Usage

You can also use command-line arguments to configure the generator.

```bash
python process_generator.py generate [OPTIONS]
```

**Required:**

*   `--name NAME`: Process name (e.g., `DryDeposition`).
*   `--type TYPE`: Process type (`emission`, `transformation`, `loss`, etc.).

**Optional:**

*   `--schemes SCHEMES`: Comma-separated list of scheme names.
*   `--species SPECIES`: Comma-separated list of species.
*   `--description TEXT`: Process description.
*   `--config FILE`: Use a YAML configuration file.

## Tutorial: Creating a Sea Salt Emission Process

Let's walk through creating a `seasalt` emission process.

1.  **Create the configuration file** (`seasalt.yaml`):

    ```yaml
    process:
      name: SeaSalt
      type: emission
      description: "Sea salt aerosol emissions."
      schemes: [gong97, gong03]
      species: [SS_A01, SS_A02]
      required_met_fields: [U10M, V10M, SST, FRSEAICE, FROCEAN]
    ```

2.  **Run the generator:**

    ```bash
    python process_generator.py generate --config configs/seasalt.yaml
    ```

3.  **Review the generated files:**

    The generator creates the following structure:

    ```
    src/process/seasalt/
    ├── ProcessSeaSaltInterface_Mod.F90
    ├── SeaSaltCommon_Mod.F90
    ├── SeaSaltCreator_Mod.F90
    ├── schemes/
    │   ├── SeaSalt_gong97_Scheme_Mod.F90
    │   └── SeaSalt_gong03_Scheme_Mod.F90
    └── CMakeLists.txt
    ```

4.  **Implement the science:**

    Open the scheme files (e.g., `SeaSalt_gong97_Scheme_Mod.F90`) and add your scientific algorithm to the `compute_scheme` subroutine. The generator has already set up the necessary inputs and outputs.

    ```fortran
    subroutine compute_gong97_scheme(u10m, v10m, sst, tendency, rc)
      real(fp), intent(in) :: u10m, v10m, sst
      real(fp), intent(out) :: tendency
      integer, intent(out) :: rc

      ! Your science goes here
      tendency = calculate_seasalt_flux(u10m, v10m, sst)
      rc = 0
    end subroutine
    ```

5.  **Build and test:**

    The generator also creates unit tests. Build and run them to verify your implementation.

    ```bash
    cd /path/to/catchem/build
    make
    ctest -R test_seasalt
    ```

## Advanced Configuration

The process generator supports a hierarchical YAML structure for more complex configurations.

```yaml
processes:
  my_process:
    name: "MyProcess"
    version: "1.0.0"
    active: true
    scheme: "scheme1"

    species:
      - name: "O3"
        active: true

    scheme_config:
      scale_factor: 1.0
      # other scheme-specific parameters

    diagnostics:
      rate:
        active: true
        units: "mol/m3/s"
```

The generator uses this structure to create a unified configuration type in the `...Common_Mod.F90` file, which is then used by the process and its schemes.

## Scheme State Management

By default, schemes are assumed to be stateless. However, for schemes that need to maintain state between time steps (e.g., for iterative solvers or time-dependent calculations), you can specify this in the configuration.

### When to Use Scheme State

-   **Working Arrays:** For complex calculations that require temporary arrays that are expensive to allocate on every time step.
-   **Persistent Variables:** To remember values from previous time steps.
-   **Scheme-Specific Diagnostics:** To generate diagnostics that are internal to a specific scheme.

### Example Configuration with State

```yaml
schemes:
  - name: "my_stateful_scheme"
    has_persistent_state: true
    working_arrays:
      - name: "size_distribution"
        dimensions: ["n_bins", "n_columns"]
    persistent_variables:
      - name: "last_time"
        type: "real(fp)"
        default: "0.0"
```

When `has_persistent_state` is true, the generator will create a `...SchemeState` type in the scheme module to hold the specified arrays and variables.

## Generated Documentation

The generator creates a complete set of documentation for your process in `docs/processes/`. This includes:

-   `index.md`: An overview of the process.
-   `schemes.md`: Detailed descriptions of each scheme.
-   `examples.md`: Usage examples.

You should edit these files to add scientific details, references, and validation results.

## Best Practices

-   **Use configuration files:** For any non-trivial process, a YAML configuration file is easier to manage than command-line arguments.
-   **Specify required met fields:** Always list the meteorological fields your process needs in the `required_met_fields` section of your configuration. This allows the framework to provide them to your process efficiently.
-   **Start with the science:** The generator handles the boilerplate, so you can focus on implementing the scientific algorithm in the scheme modules.
-   **Test thoroughly:** The generated unit tests are a starting point. Add more tests to cover all aspects of your implementation.
-   **Document well:** The generated documentation provides a template. Fill it in with the scientific details of your process.
-   **Default to stateless schemes:** Only add persistent state to a scheme when it is truly necessary. This will result in cleaner and more efficient code.
