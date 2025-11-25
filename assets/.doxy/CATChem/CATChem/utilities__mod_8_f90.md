

# File utilities\_mod.F90



[**FileList**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**utilities\_mod.F90**](utilities__mod_8_f90.md)

[Go to the source code of this file](utilities__mod_8_f90_source.md)

_General utility functions for CATChem._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**utilities\_mod**](namespaceutilities__mod.md) <br> |




















































## Detailed Description




**Author:**

CATChem Development Team 




**Date:**

2025 




**Version:**

2.0


This module provides general utility functions for unit conversions, calculations, and validation that are used throughout CATChem.


The utilities module includes:
* Unit conversion functions (pressure, temperature)
* Atmospheric calculation utilities
* Validation functions for physical consistency
* Common mathematical operations
* Error checking and validation helpers




## Usage Example




```Fortran
use utilities_mod
real(fp) :: pressure_hpa, pressure_pa
integer :: rc

pressure_pa = convert_pressure_units(pressure_hpa, 'hPa', 'Pa', rc)
call validate_atmospheric_constants(rc)
```
 



    

------------------------------
The documentation for this class was generated from the following file `src/core/utilities_mod.F90`

