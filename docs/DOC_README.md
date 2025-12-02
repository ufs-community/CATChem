# CATChem Documentation

This directory contains the source files for CATChem's documentation website, built with [MkDocs](https://www.mkdocs.org/) and the [Material theme](https://squidfunk.github.io/mkdocs-material/).

## Quick Start

### Prerequisites

```bash
# Install Python 3.8+
python --version

# Install documentation dependencies
pip install -r docs/requirements.txt
```

### Development Server

```bash
# Start live development server
mkdocs serve

# Open browser to http://localhost:8000
```

### Building Documentation

```bash
# Build static site
mkdocs build

# Build with API documentation
mkdocs build --config-file mkdocs.yml
```

## Structure

```
docs/
├── index.md                  # Main landing page
├── user-guide/               # End-user documentation
├── processes/                # Process-specific documentation used in the user-guide
├── developer-guide/          # Technical documentation for developers
├── api/                      # API reference
├── ufschem/                  # How CATChem is integrated into the UFS to create UFS-Chem
├── community/                # How to get involved?
├── evaluation/               # How CATChem and UFS-Chem are evaluated
├── assets/                   # Images, stylesheets, scripts
└── ../mkdocs.yml             # Documentation configuration
```

## Features

### Material Theme with NOAA Colors
- Custom NOAA color palette (blue, teal, navy)
- Responsive design for all devices
- Dark/light mode toggle
- Search functionality

### MkDoxy Integration
- Automatic API documentation from Fortran source code
- Doxygen-style comment parsing
- Cross-references between code and docs

### Enhanced Markdown
- Mermaid diagrams for architecture visualization
- Code syntax highlighting with Fortran support
- Admonitions for notes, warnings, tips
- Tabbed content sections
- Mathematical equations with MathJax

## Contributing

Update documentation alongside code changes and test locally with `mkdocs serve`.

For detailed information, see the full documentation setup guide in the built documentation.
