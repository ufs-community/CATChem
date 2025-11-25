

# File DryDepScheme\_WESELY\_Mod.F90



[**FileList**](files.md) **>** [**drydep**](dir_57fb5aa14ddb2cd518a6d90b65ffd000.md) **>** [**schemes**](dir_5a3c86e36f17958630366ebc2b7ca21b.md) **>** [**DryDepScheme\_WESELY\_Mod.F90**](_dry_dep_scheme___w_e_s_e_l_y___mod_8_f90.md)

[Go to the source code of this file](_dry_dep_scheme___w_e_s_e_l_y___mod_8_f90_source.md)

_Wesely 1989 gas dry deposition scheme._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**drydepscheme\_wesely\_mod**](namespacedrydepscheme__wesely__mod.md) <br> |




















































## Detailed Description


Pure science kernel for wesely scheme in drydep process. This module contains ONLY the computational algorithm with NO infrastructure dependencies. Uses only basic Fortran types for maximum portability and reusability.


SCIENCE CUSTOMIZATION GUIDE:
* Modify the algorithm in compute\_wesely (search for "TODO")
* Add scheme-specific helper subroutines as needed
* Update physical constants for your scheme
* Customize the environmental response functions




INFRASTRUCTURE RESPONSIBILITIES (handled by host model):
* Parameter initialization and validation
* Input array validation and error handling
* Memory management and array allocation
* Integration with host model time stepping




Generated on: 2025-11-13T14:35:43.237148 Author: Wei Li Reference: Wesely, M. L. [1989] Parameterization of surface resistances to gaseous dry deposition... 


    

------------------------------
The documentation for this class was generated from the following file `src/process/drydep/schemes/DryDepScheme_WESELY_Mod.F90`

