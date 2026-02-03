

# File SettlingScheme\_GOCART\_Mod.F90



[**FileList**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**settling**](dir_1a0bba2ffdf6e6637fcb76856471cb75.md) **>** [**schemes**](dir_34df91cc26d24067840a7381fe21b817.md) **>** [**SettlingScheme\_GOCART\_Mod.F90**](_settling_scheme___g_o_c_a_r_t___mod_8_f90.md)

[Go to the source code of this file](_settling_scheme___g_o_c_a_r_t___mod_8_f90_source.md)

_GOCART gravitational settling scheme._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**settlingscheme\_gocart\_mod**](namespacesettlingscheme__gocart__mod.md) <br> |




















































## Detailed Description


Pure science kernel for gocart scheme in settling process. This module contains ONLY the computational algorithm with NO infrastructure dependencies. Uses only basic Fortran types for maximum portability and reusability.


SCIENCE CUSTOMIZATION GUIDE:
* Modify the algorithm in compute\_gocart (search for "TODO")
* Add scheme-specific helper subroutines as needed
* Update physical constants for your scheme
* Customize the environmental response functions




INFRASTRUCTURE RESPONSIBILITIES (handled by host model):
* Parameter initialization and validation
* Input array validation and error handling
* Memory management and array allocation
* Integration with host model time stepping




Generated on: 2025-12-17T15:27:52.203209 Author: Wei Li Reference: GOCART2G process library Chem\_SettlingSimple function 


    

------------------------------
The documentation for this class was generated from the following file `src/process/settling/schemes/SettlingScheme_GOCART_Mod.F90`

