# CATChem Process Generator Templates Inventory

This directory contains the Jinja2 templates used by `tools/process_generator/process_generator.py` to generate process packages.

## Active Templates

* `catchem_process.hpp.j2` - C++ Process class header template (`catchem_process_<name>.hpp`)
* `catchem_process.cpp.j2` - C++ Process class source template (`catchem_process_<name>.cpp`)
* `science_bridge.F90.j2` - Fortran `BIND(C)` ScienceBridge module template (`<Name>ScienceBridge.F90`)
* `process_common.F90.j2` - Fortran process common module & config types (`<Name>Common_Mod.F90`)
* `scheme_module.F90.j2` - Pure Fortran science scheme module template (`schemes/<Name>Scheme_<SCHEME>_Mod.F90`)
* `CMakeLists.txt.j2` - Process CMake build script template (`CMakeLists.txt`)
* `schemes_CMakeLists.txt.j2` - Schemes CMake build script template (`schemes/CMakeLists.txt`)
* `test_science.f90.j2` - Standalone Fortran CTest science test template (`tests/test_<name>_science.f90`)
* `test_CMakeLists.txt.j2` - Test CMake build script template (`tests/process/<name>/CMakeLists.txt`)
* `process_documentation.md.j2` - Process README documentation template (`README.md`)

## Legacy Templates (Removed in v2.0 C++ Core Architecture)

* `process_interface.F90.j2` (Removed - replaced by `catchem_process.hpp.j2` / `cpp.j2` and `science_bridge.F90.j2`)
* `process_creator.F90.j2`, `process_creator_new.F90.j2`, `process_creator_old.F90.j2` (Removed - process registration handled by C++ `ProcessRegistry`)
* `integration_test.F90.j2`, `unit_test.F90.j2` (Removed - replaced by standalone science unit test `test_science.f90.j2`)
