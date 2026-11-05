

# File constants.F90



[**FileList**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**constants.F90**](constants_8_f90.md)

[Go to the source code of this file](constants_8_f90_source.md)

_Physical and mathematical constants for CATChem._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**constants**](namespaceconstants.md) <br> |




















































## Detailed Description




**Author:**

CATChem Development Team 




**Date:**

2025 




**Version:**

2.0


This module provides fundamental physical and mathematical constants used throughout the CATChem atmospheric chemistry modeling system.


The constants module defines all physical constants, conversion factors, and mathematical constants used in atmospheric chemistry calculations. All values are given in SI units unless otherwise specified.


## Usage Example




```Fortran
use constants
real(fp) :: air_density
air_density = pressure / (rd * temperature)
```
 



    

------------------------------
The documentation for this class was generated from the following file `src/core/constants.F90`

