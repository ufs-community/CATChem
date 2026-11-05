

# File ConfigManager\_Mod.F90



[**FileList**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**ConfigManager\_Mod.F90**](_config_manager___mod_8_f90.md)

[Go to the source code of this file](_config_manager___mod_8_f90_source.md)

_Enhanced configuration management for CATChem._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**configmanager\_mod**](namespaceconfigmanager__mod.md) <br> |




















































## Detailed Description




**Author:**

CATChem Development Team 




**Date:**

2025 




**Version:**

2.0


This module provides a modern, flexible configuration management system for CATChem, supporting modular configuration loading, validation, and runtime configuration updates.


The ConfigManager consolidates and modernizes CATChem configuration management by replacing two previous modules:


**Replaces config\_opt\_mod.F90:**
* Legacy config support removed, now using modern YAML-based ConfigDataType
* All configuration parameters and options
* Process-specific configuration flags




**Replaces config\_mod.F90:**
* Read\_Input\_File subroutine -&gt; ConfigManagerload\_from\_file() method
* YAML parsing logic with enhanced error handling
* Configuration validation and initialization




**New features:**
* Hierarchical configuration loading with includes and inheritance
* Schema-based validation with detailed error reporting
* Plugin-based process configuration loading
* Runtime configuration updates and hot-reloading
* Configuration templates and presets
* Environment variable and command-line override support




## Usage Example




```Fortran
use configmanager_mod
type(ConfigManagerType) :: config_mgr
type(StateContainer) :: container
integer :: rc

call config_mgr%init(rc)
call config_mgr%load_from_file('CATChem_config.yml', rc)
call config_mgr%apply_to_container(container, rc)
```
 



    

------------------------------
The documentation for this class was generated from the following file `src/core/ConfigManager_Mod.F90`

