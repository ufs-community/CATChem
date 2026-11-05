

# File DiagnosticManager\_Mod.F90



[**FileList**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**DiagnosticManager\_Mod.F90**](_diagnostic_manager___mod_8_f90.md)

[Go to the source code of this file](_diagnostic_manager___mod_8_f90_source.md)

_Central diagnostic manager integrating with CATChem framework._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**diagnosticmanager\_mod**](namespacediagnosticmanager__mod.md) <br> |




















































## Detailed Description




**Author:**

CATChem Development Team 




**Date:**

2025 




**Version:**

2.0


This module provides a central diagnostic manager that integrates the dynamic diagnostic system with the existing CATChem framework (StateContainer, ProcessManager, ConfigManager, ErrorManager).


The DiagnosticManager:
* Manages diagnostic registries for all processes
* Integrates with StateContainer for state management
* Uses ErrorManager for consistent error handling
* Supports ConfigManager for output configuration
* Provides centralized diagnostic collection and output
* Replaces static diagstate\_mod with dynamic system




## Usage Example




```Fortran
use diagnosticmanager_mod
type(DiagnosticManagerType) :: diag_mgr
type(StateContainerType) :: container
integer :: rc

! Initialize diagnostic manager
call diag_mgr%init(container, rc)

! Processes register their diagnostics via process manager
call process_mgr%run_all(container, rc)

! Collect and write diagnostics
call diag_mgr%collect_all_diagnostics(container, rc)
call diag_mgr%write_output(container, rc)
```
 



    

------------------------------
The documentation for this class was generated from the following file `src/core/DiagnosticManager_Mod.F90`

