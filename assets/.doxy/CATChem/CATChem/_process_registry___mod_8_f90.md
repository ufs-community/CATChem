

# File ProcessRegistry\_Mod.F90



[**FileList**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**ProcessRegistry\_Mod.F90**](_process_registry___mod_8_f90.md)

[Go to the source code of this file](_process_registry___mod_8_f90_source.md)

_Process registration system for dynamic process discovery._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**processregistry\_mod**](namespaceprocessregistry__mod.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| interface | [**processcreatorinterface**](interfaceprocessregistry__mod_1_1processcreatorinterface.md) <br>_Function pointer interface for process creators._  |


















































## Detailed Description




**Author:**

CATChem Development Team 




**Date:**

2025 




**Version:**

2.0


This module provides a registry system for atmospheric processes, enabling dynamic discovery and instantiation of processes at runtime. It supports the plugin architecture outlined in the refactoring plan.


**Key Features:**
* Dynamic process registration and discovery
* Category-based process organization
* Runtime metadata querying
* Process availability checking
* Plugin-style architecture support




**Usage Example:** 
```Fortran
use processregistry_mod
type(ProcessRegistryType) :: registry
class(ProcessInterface), allocatable :: process

call registry%init()
call registry%register_process('dust_fengsha', create_dust_fengsha)
call registry%create_process('dust_fengsha', process, rc)
```
 


    

------------------------------
The documentation for this class was generated from the following file `src/core/ProcessRegistry_Mod.F90`

