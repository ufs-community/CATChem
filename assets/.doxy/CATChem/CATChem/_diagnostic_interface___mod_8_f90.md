

# File DiagnosticInterface\_Mod.F90



[**FileList**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**DiagnosticInterface\_Mod.F90**](_diagnostic_interface___mod_8_f90.md)

[Go to the source code of this file](_diagnostic_interface___mod_8_f90_source.md)

_Dynamic diagnostic system interfaces and types._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**diagnosticinterface\_mod**](namespacediagnosticinterface__mod.md) <br> |




















































## Detailed Description




**Author:**

CATChem Development Team 




**Date:**

2025 




**Version:**

2.0


This module provides interfaces and types for the dynamic diagnostic system, allowing processes to register and manage their own diagnostic outputs at runtime.


The diagnostic system supports:
* Multiple data types (scalar, 1D, 2D, 3D arrays)
* Flexible metadata (units, description, output frequency)
* Process-specific diagnostic registration
* Runtime diagnostic query and collection
* Optional diagnostic output control




## Usage Example




```Fortran
use diagnosticinterface_mod
type(DiagnosticFieldType) :: diag_field
integer :: rc

call diag_field%create('dust_flux', 'Total dust emission flux', &
                       'kg m-2 s-1', diag_real_2d, rc)
call diag_mgr%register_diagnostic('dust_process', diag_field, rc)
```
 



    

------------------------------
The documentation for this class was generated from the following file `src/core/DiagnosticInterface_Mod.F90`

