# CATChem Process Generator

A powerful Python tool for generating standardized atmospheric chemistry and physics process implementations matching the CATChem C++ core architecture (`catchem::Core`, `catchem::StateManager`, `catchem::ProcessRegistry`).

## Features

- 🎯 **C++ Core Architecture Integration** - Generates C++ process wrapper headers/sources (`catchem_process_<name>.hpp/cpp`), flat C-interoperable Fortran science bridges (`<ClassName>ScienceBridge.F90`), and pure Fortran science schemes
- ⚙️ **YAML Configuration** - Simple YAML configuration files define process parameters, schemes, species filters, and diagnostics
- 🔍 **Validation** - Built-in configuration validation (`validate` subcommand) checking schema syntax and meteorological field compatibility
- 📦 **Kokkos & CMake Ready** - Generates CMake build files supporting `CATCHEM_ENABLE_KOKKOS` and standalone CTest science unit tests (`test_<name>_science.f90`)
- 🧪 **Template Command** - Generates starter YAML configs for `emission`, `chemistry`, `deposition`, `settling`, and `plume_rise` process types

## Usage

### 1. Generate a Process Package

```bash
# Generate C++/Fortran process package from YAML config
python process_generator.py generate --config configs/dust_emission.yaml --output-dir src/process/dust_test --force
```

### 2. Validate Process Configuration

```bash
# Validate YAML configuration syntax and field completeness
python process_generator.py validate --config configs/dust_emission.yaml
```

### 3. Generate Starter YAML Configuration

```bash
# Generate starter YAML template for chemistry process
python process_generator.py template --type chemistry --output my_chem_template.yaml
```

## Generated Process Layout

```text
src/process/<process_name>/
├── catchem_process_<name>.hpp       # C++ Process Interface Class Header
├── catchem_process_<name>.cpp       # C++ Process Implementation & C-Registration
├── <ClassName>ScienceBridge.F90      # Fortran BIND(C) Science Bridge
├── <ClassName>Common_Mod.F90         # Fortran Common Types & Config
├── schemes/                          # Pure Fortran Science Schemes
│   ├── <ClassName>Scheme_<SCHEME>_Mod.F90
│   └── CMakeLists.txt
├── CMakeLists.txt                    # Process CMake Build Script
└── README.md                         # Process Documentation
```

## Configuration Format

Process configurations use YAML format:

```yaml
process:
  name: "MyChemicalProcess"
  description: "Custom atmospheric chemistry process"
  version: "1.0.0"

variables:
  - name: "temperature"
    type: "real"
    units: "K"
    description: "Temperature"

  - name: "concentration_o3"
    type: "real"
    units: "molecules/cm3"
    description: "Ozone concentration"

chemistry:
  reactions:
    - name: "R1"
      equation: "O + O2 + M -> O3 + M"
      rate_constant: "6.0e-34 * (300/T)**2.4"
```

## Dependencies

The process generator supports two installation approaches:

### Minimal Installation (requirements.txt)
For users who just want to run the tool:
```bash
pip install -r requirements.txt
```

**Core Dependencies:**
- **PyYAML** (≥6.0) - YAML configuration parsing
- **Jinja2** (≥3.1.0) - Template engine
- **click** (≥8.0.0) - Command-line interface
- **pathlib-mate** (≥1.0.0) - Enhanced path operations
- **colorama** (≥0.4.4) - Cross-platform colored terminal output
- **rich** (≥13.0.0) - Rich terminal formatting

### Full Development Installation (pyproject.toml)
For contributors and developers:
```bash
pip install -e ".[dev,test,docs,validation]"
```

**Additional Development Dependencies:**
- **pytest** - Testing framework
- **black** - Code formatting
- **isort** - Import sorting
- **flake8** - Linting
- **mypy** - Type checking
- **pre-commit** - Git hooks

**Optional Extensions:**
- **sphinx** - Documentation generation
- **jsonschema** - Configuration validation
- **yamllint** - YAML file validation
- **fortls** - Fortran language server support

## Project Structure

```
process_generator/
├── process_generator.py     # Main generator script
├── templates/              # Jinja2 templates
│   ├── process_interface.F90.j2
│   ├── process_creator.F90.j2
│   └── ...
├── configs/               # Example configurations
│   ├── example_process.yaml
│   └── ...
├── pyproject.toml        # Package configuration
├── setup_dev_env.sh      # Development setup script
└── README.md            # This file
```

## Development Workflow

### Testing
```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=process_generator

# Run specific test categories
pytest -m unit
pytest -m integration
pytest -m template
```

### Code Quality
```bash
# Format code
black .
isort .

# Check linting
flake8

# Type checking
mypy process_generator.py

# Pre-commit hooks
pre-commit install
pre-commit run --all-files
```

### Documentation
```bash
# Build documentation
sphinx-build docs docs/_build

# Live documentation server
sphinx-autobuild docs docs/_build
```

## Templates

The generator includes several Jinja2 templates:

- `process_interface.F90.j2` - Main process implementation
- `process_creator.F90.j2` - Factory pattern creator
- `process_config.yaml.j2` - Configuration template
- `process_test.F90.j2` - Unit test template

### Custom Templates

You can create custom templates by:

1. Adding `.j2` files to the `templates/` directory
2. Using the CATChem template variables
3. Following the established naming conventions

## Command Line Interface

```bash
# Basic usage
python process_generator.py COMMAND [OPTIONS]

# Main commands
python process_generator.py generate --config CONFIG_FILE
python process_generator.py validate --config CONFIG_FILE
python process_generator.py template --type TYPE --output FILE
python process_generator.py fields --type TYPE [--verbose]

# Generate options
  --config CONFIG_FILE    YAML configuration file (required)
  --output DIR           Output directory for generated files
  --templates DIR        Custom template directory
  --metstate FILE        Path to MetState file for field discovery
  --verbose              Enable verbose output

# Validate options
  --config CONFIG_FILE    YAML configuration file (required)
  --metstate FILE        Path to MetState file for field discovery
  --verbose              Enable verbose output

# Template options
  --type TYPE            Template type: emission, chemistry, transport
  --output FILE          Output file for template (required)

# Fields options
  --type TYPE            Field type: 2d, 3d, categorical, all
  --metstate FILE        Path to MetState file for field discovery
  --verbose              Show detailed field lists

# Examples
python process_generator.py generate --config seasalt.yaml
python process_generator.py generate --config seasalt.yaml --metstate ../../src/core/metstate_mod.F90
python process_generator.py validate --config seasalt.yaml --verbose
python process_generator.py fields --type all --verbose
python process_generator.py template --type emission --output template.yaml
```

## Integration with CATChem Build System

Generated processes automatically integrate with the CATChem CMake build system:

```cmake
# Generated processes include proper CMake targets
# No manual build system modification required
```

## Troubleshooting

### Common Issues

1. **Import errors**: Ensure all dependencies are installed
   ```bash
   pip install -e ".[dev]"
   ```

2. **MetState file not found**: Use explicit path or check auto-discovery
   ```bash
   # Check if MetState file exists
   ls -la ../../src/core/metstate_mod.F90

   # Use explicit path
   python process_generator.py generate --config config.yaml --metstate /path/to/metstate_mod.F90
   ```

3. **Field not discovered**: Check field lists and MetState definition
   ```bash
   # Check discovered fields
   python process_generator.py fields --type all --verbose

   # Verify field exists in MetState
   grep -n "FIELD_NAME" ../../src/core/metstate_mod.F90
   ```

4. **Template not found**: Check template directory path
   ```bash
   python process_generator.py generate --config config.yaml --templates ./templates
   ```

5. **YAML parsing errors**: Validate YAML syntax
   ```bash
   # Test YAML parsing
   python -c "import yaml; yaml.safe_load(open('config.yaml'))"
   ```

### Getting Help

- **Command help**: `python process_generator.py --help`
- **Field discovery**: `python process_generator.py fields --type all --verbose`
- **Verbose output**: Use `--verbose` flag for detailed information
- **Examples**: Review configurations in `configs/`

## Contributing

1. Set up development environment: `./setup_dev_env.sh`
2. Make changes with tests: `pytest`
3. Ensure code quality: `pre-commit run --all-files`
4. Update documentation as needed

## License

Apache License 2.0 - See the main CATChem repository for details.
