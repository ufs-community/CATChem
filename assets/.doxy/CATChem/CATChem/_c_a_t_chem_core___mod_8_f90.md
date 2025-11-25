

# File CATChemCore\_Mod.F90



[**FileList**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**CATChemCore\_Mod.F90**](_c_a_t_chem_core___mod_8_f90.md)

[Go to the source code of this file](_c_a_t_chem_core___mod_8_f90_source.md)

_Central CATChem core framework with unified component management._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**catchemcore\_mod**](namespacecatchemcore__mod.md) <br> |




















































## Detailed Description




**Author:**

CATChem Development Team 




**Date:**

2025 




**Version:**

2.0


This module provides the central CATChem core that owns and manages all major components (StateManager, GridManager, DiagnosticManager, etc.) and provides controlled access between them. This eliminates circular dependencies and provides a clean, centralized architecture.


Key Features:
* Centralized component lifecycle management
* Controlled inter-component communication
* Simplified initialization and error handling
* Better testability and modularity




## Usage Example




```Fortran
use catchemcore_mod
type(CATChemCoreType) :: core
integer :: rc

! Initialize the core with configuration
call core%init('config.yaml', rc)

! Run simulation timesteps
do timestep = 1, n_timesteps
   call core%run_timestep(timestep, dt, rc)
enddo

! Clean shutdown
call core%finalize(rc)
```
 



    

------------------------------
The documentation for this class was generated from the following file `src/core/CATChemCore_Mod.F90`

