# CATChem Process Generator

A powerful Python tool for generating standardized atmospheric chemistry process implementations for the CATChem framework.

## Features

- 🎯 **Template-driven generation** - Uses Jinja2 templates for consistent code structure
- ⚙️ **YAML configuration** - Simple configuration files define process parameters
- 🔍 **Validation** - Built-in validation for process configurations
- 🚀 **Modern Fortran** - Generates clean, modern Fortran 2008+ code
- 🔧 **CATChem integration** - Seamlessly integrates with CATChem infrastructure
- 🌡️ **MetState Field Detection** - Automatically discovers and classifies meteorological fields

## MetState Field Detection

The process generator includes automatic meteorological field discovery and classification from the CATChem MetState definition. This ensures generated processes are synchronized with the actual MetState implementation.

### Field Classification

The generator automatically categorizes fields into three types:

- **2D Surface Fields**: Surface-level scalar fields (e.g., `FROCEAN`, `SST`, `USTAR`, `PBLH`)
- **3D Atmospheric Fields**: Vertical profile arrays (e.g., `T`, `QV`, `P`, `U`, `V`)
- **Categorical Fields**: Special-purpose arrays with non-vertical dimensions (e.g., `SOILM`, `FRLANDUSE`)

### Field Discovery Commands

```bash
# Show all discovered fields (summary)
python process_generator.py fields --type all

# Show detailed list of 2D surface fields
python process_generator.py fields --type 2d --verbose

# Show 3D atmospheric fields
python process_generator.py fields --type 3d --verbose

# Show categorical fields
python process_generator.py fields --type categorical --verbose

# Use specific MetState file
python process_generator.py fields --metstate /path/to/metstate_mod.F90 --type all
```

### Generated VirtualMet Pattern

The generator produces code that uses the modern VirtualMet pattern:

```fortran
! Generated process implementation
type(VirtualMetType), pointer :: met => null()

! Get meteorological data pointer
met => column%get_met()

! Access fields based on automatic classification
frocean(1) = met%FROCEAN    ! 2D surface field - scalar access
sst(1) = met%SST            ! 2D surface field - scalar access

! For 3D fields (when affects_full_column = true):
! temperature(:) = met%T     ! 3D atmospheric field - array access
```

### MetState File Discovery

The generator automatically finds the MetState file using these search paths:

1. Automatic discovery: `../../src/core/metstate_mod.F90` (relative to generator)
2. Explicit path: `--metstate /path/to/metstate_mod.F90`
3. Fallback: Uses hardcoded field lists if MetState file not found

### Benefits

- **Always Up-to-Date**: New MetState fields are immediately available
- **Type Safety**: Based on actual Fortran type definitions
- **Reduced Maintenance**: No manual field list updates required
- **Accurate Classification**: Respects actual field dimensions and usage patterns

## Quick Start

### 1. Set up development environment

```bash
# Quick setup - Full development environment (recommended for contributors)
./setup_dev_env.sh my-env

# Minimal setup - Just run the tool (for users)
./setup_dev_env.sh my-env 3.10 minimal

# Manual setup options:

# Option A: Full development setup with pyproject.toml
conda create -n catchem-process-gen python=3.10
conda activate catchem-process-gen
pip install -e ".[dev,test]"

# Option B: Minimal setup with requirements.txt
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Option C: Manual venv with pyproject.toml
python3 -m venv venv
source venv/bin/activate
pip install -e .
```

### 2. Generate a process

```bash
# Basic usage with automatic MetState field discovery
python process_generator.py generate --config configs/example_process.yaml

# With explicit MetState file path
python process_generator.py generate --config configs/example_process.yaml --metstate ../../src/core/metstate_mod.F90

# With custom output directory
python process_generator.py generate --config configs/example_process.yaml --output /path/to/output

# Validate configuration only
python process_generator.py validate --config configs/example_process.yaml

# Inspect discovered meteorological fields
python process_generator.py fields --type all --verbose
```

### 3. Integration with CATChem

The generated processes automatically integrate with the CATChem framework:

```fortran
! The generated process includes:
! - Proper ProcessInterface implementation
! - VirtualMet pattern with automatic field classification
! - Error handling integration
! - State variable management
! - Diagnostic capability
! - Registry-compatible factory pattern
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
