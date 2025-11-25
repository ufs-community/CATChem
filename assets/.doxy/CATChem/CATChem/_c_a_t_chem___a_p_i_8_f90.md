

# File CATChem\_API.F90



[**FileList**](files.md) **>** [**api**](dir_da61e3e9a357748887e3ca8d7c5a0c16.md) **>** [**CATChem\_API.F90**](_c_a_t_chem___a_p_i_8_f90.md)

[Go to the source code of this file](_c_a_t_chem___a_p_i_8_f90_source.md)

_Streamlined CATChem API for host model integration._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**catchem\_api**](namespacecatchem__api.md) <br> |




















































## Detailed Description




**Author:**

CATChem Development Team 




**Date:**

2025 




**Version:**

2.1


This module provides a streamlined, lightweight API for integrating CATChem into different modeling architectures. It leverages the existing core architecture without duplicating functionality, providing clean interfaces for the most common integration patterns.


Key design principles:
* Lightweight wrapper around existing core components
* Support for multiple processes and run phases
* Streamlined data exchange with host models
* Clear error handling and status reporting
* No duplication of existing types (ConfigManager, StateManager, etc.)




Usage pattern:
* Initialize with configuration file
* Setup grid geometry
* Add processes as needed
* Configure run phases (optional)
* Execute timesteps or phases
* Exchange data with host model
* Retrieve diagnostics
* Finalize 




    

------------------------------
The documentation for this class was generated from the following file `src/api/CATChem_API.F90`

