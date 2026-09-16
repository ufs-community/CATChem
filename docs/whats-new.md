# What's New

## Recent Updates

### Version 2.0.0 (Latest)

CATChem version 2.0 was released in July 2026 (**[v2.0.0](https://github.com/ufs-community/CATChem/releases/tag/v2.0.0)**) and is integrated into the **[UFS Weather Model](https://github.com/ufs-community/ufs-weather-model)** as a NUOPC component. See **[UFS-Chem](ufschem/index.md)** for details.

#### 🚀 Major Features
- **Modernized Process Infrastructure** - Complete overhaul of process and scheme architecture
- **Stokes Settling Scheme** - Advanced physics-based settling with slip correction
- **Improved Documentation** - Comprehensive user and developer guides with API references
- **YAML Configuration** - Flexible, modern configuration system
- **StateContainer** - All core states use a common interface and are accessible through a common "container"
- **ConfigManager** - Hierarchical YAML configuration loading with file inheritance, environment variable substitution, and validation

#### 🔧 Process Updates
- **Settling Process** - Renamed from `settlingvelocity` to `settling` with new Stokes scheme
- **Modular Schemes** - Support for multiple schemes per process with subdirectories
- **Enhanced Diagnostics** - Process-level diagnostic output and monitoring
- **Template Generator** - Improved code generation for new processes

#### 📚 Documentation
- **MkDocs + Material Theme** - Modern documentation with NOAA color palette
- **API Documentation** - Auto-generated API docs with MkDoxy integration
- **User Guides** - Comprehensive process documentation and tutorials
- **Developer Resources** - Architecture guides and contribution guidelines

#### 🛠️ Technical Improvements
- **CMake Modernization** - Updated build system for modular components
- **Error Handling** - Robust error management and validation
- **Code Standards** - Consistent coding practices and style guidelines

## Previous Versions

### Version 1.x
- Legacy settling velocity implementations
- Basic process framework
- Original CCPP and NUOPC integration

---

For detailed release notes and migration guides, visit our **[GitHub releases](https://github.com/UFS-Community/CATChem/releases)**.
