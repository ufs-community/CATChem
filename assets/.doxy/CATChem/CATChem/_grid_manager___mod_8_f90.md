

# File GridManager\_Mod.F90



[**FileList**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**GridManager\_Mod.F90**](_grid_manager___mod_8_f90.md)

[Go to the source code of this file](_grid_manager___mod_8_f90_source.md)

_Advanced grid management with column virtualization support._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**gridmanager\_mod**](namespacegridmanager__mod.md) <br> |




















































## Detailed Description




**Author:**

CATChem Development Team 




**Date:**

2025 




**Version:**

1.0


This module provides a comprehensive grid management system that enables full column virtualization at the process level while maintaining 3D spatial awareness. All processes operate on virtual columns but exist in 3D space through the grid manager.


Key Features:
* Grid virtualization: processes see 1D columns, grid manager handles 3D
* Unified grid interface for 1D, 2D, and 3D models
* Zero-copy column access through smart pointers
* Thread-safe parallelization support (OpenMP/OpenACC)
* Dynamic grid configuration and decomposition
* Grid metrics and geometric calculations
* Coordinate transformations and projections 




    

------------------------------
The documentation for this class was generated from the following file `src/core/GridManager_Mod.F90`

