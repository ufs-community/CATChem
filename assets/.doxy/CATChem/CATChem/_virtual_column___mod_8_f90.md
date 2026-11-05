

# File VirtualColumn\_Mod.F90



[**FileList**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**VirtualColumn\_Mod.F90**](_virtual_column___mod_8_f90.md)

[Go to the source code of this file](_virtual_column___mod_8_f90_source.md)

_Virtual column data container for CATChem processes with macro-generated meteorological fields._ [More...](#detailed-description)

* `#include "virtualmet_type.inc"`
* `#include "virtualmet_cleanup.inc"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**virtualcolumn\_mod**](namespacevirtualcolumn__mod.md) <br> |




















































## Detailed Description




**Author:**

CATChem Development Team 




**Date:**

2025 




**Version:**

3.0


This module provides a column-based data container that allows processes to work with vertical column data extracted from 3D grid fields. The VirtualColumn contains a VirtualMetType with direct pointers to meteorological fields, eliminating data copying overhead. The VirtualMetType definition is now generated automatically from MetState field definitions using macros. 


    

------------------------------
The documentation for this class was generated from the following file `src/core/VirtualColumn_Mod.F90`

