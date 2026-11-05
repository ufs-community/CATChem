

# File met\_utilities\_mod.F90



[**FileList**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**met\_utilities\_mod.F90**](met__utilities__mod_8_f90.md)

[Go to the source code of this file](met__utilities__mod_8_f90_source.md)

_Meteorological utility functions for CATChem._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**met\_utilities\_mod**](namespacemet__utilities__mod.md) <br> |




















































## Detailed Description




**Author:**

CATChem Development Team 




**Date:**

2025 




**Version:**

1.0


This module provides meteorological and atmospheric utility functions commonly used in atmospheric chemistry and physics, including calculations for potential temperature, virtual temperature, dew point, relative humidity, saturation vapor pressure, and more.


The met\_utilities module includes:
* Potential temperature calculation
* Virtual temperature calculation
* Dew point calculation
* Relative humidity calculation
* Saturation vapor pressure (Clausius-Clapeyron)
* Mixing ratio and specific humidity conversions
* Lapse rate calculations




## Usage Example




```Fortran
use met_utilities_mod
real(fp) :: T, p, theta, Tv, rh, Td, es
theta = potential_temperature(t, p, p0)
tv = virtual_temperature(t, qv)
td = dew_point(t, rh)
es = saturation_vapor_pressure(t)
```
 



    

------------------------------
The documentation for this class was generated from the following file `src/core/met_utilities_mod.F90`

