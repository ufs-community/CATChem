

# File SeaSaltScheme\_GONG97\_Mod.F90



[**FileList**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**seasalt**](dir_3e6fd2e121e43ca7d4114b6c0b4e05b6.md) **>** [**schemes**](dir_ec083b49fedbd640552af85049fd7226.md) **>** [**SeaSaltScheme\_GONG97\_Mod.F90**](_sea_salt_scheme___g_o_n_g97___mod_8_f90.md)

[Go to the source code of this file](_sea_salt_scheme___g_o_n_g97___mod_8_f90_source.md)

_Gong 1997 sea salt emission scheme._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**seasaltscheme\_gong97\_mod**](namespaceseasaltscheme__gong97__mod.md) <br> |




















































## Detailed Description


Pure science kernel for gong97 scheme in seasalt process. This module contains ONLY the computational algorithm with NO infrastructure dependencies. Uses only basic Fortran types for maximum portability and reusability.


SCIENCE CUSTOMIZATION GUIDE:
* Modify the algorithm in compute\_gong97 (search for "TODO")
* Add scheme-specific helper subroutines as needed
* Update physical constants for your scheme
* Customize the environmental response functions




INFRASTRUCTURE RESPONSIBILITIES (handled by host model):
* Parameter initialization and validation
* Input array validation and error handling
* Memory management and array allocation
* Integration with host model time stepping




Generated on: 2025-09-16T00:40:10.206610 Author: Barry Baker Reference: Gong et al. [1997] 


    

------------------------------
The documentation for this class was generated from the following file `src/process/seasalt/schemes/SeaSaltScheme_GONG97_Mod.F90`

