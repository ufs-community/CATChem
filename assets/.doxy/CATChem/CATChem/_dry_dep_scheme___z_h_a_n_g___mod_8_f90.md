

# File DryDepScheme\_ZHANG\_Mod.F90



[**FileList**](files.md) **>** [**drydep**](dir_57fb5aa14ddb2cd518a6d90b65ffd000.md) **>** [**schemes**](dir_5a3c86e36f17958630366ebc2b7ca21b.md) **>** [**DryDepScheme\_ZHANG\_Mod.F90**](_dry_dep_scheme___z_h_a_n_g___mod_8_f90.md)

[Go to the source code of this file](_dry_dep_scheme___z_h_a_n_g___mod_8_f90_source.md)

_Zhang et al._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**drydepscheme\_zhang\_mod**](namespacedrydepscheme__zhang__mod.md) <br> |




















































## Detailed Description


[2001] scheme with Emerson et al. [2020] updates. The Ra and Rb are still from Wesely (1989) for now.


Pure science kernel for zhang scheme in drydep process. This module contains ONLY the computational algorithm with NO infrastructure dependencies. Uses only basic Fortran types for maximum portability and reusability. Reference: (1) Wesely, M. L. (1989). Parameterization of surface resistances to gaseous dry deposition in regional-scale numerical models. Atmospheric Environment. (2) Zhang, L., Gong, S., Padro, J., & Barrie, L. (2001). A size-segregated particle dry deposition scheme for an atmospheric aerosol module. Atmospheric environment. (3) Emerson, E. W., et al. (2020). Revisiting particle dry deposition and its role in radiative effect estimates. PNAS, 117(42), 26076-26082. (4) Most of the codes are adopted from GEOS-Chem drydep\_mod.F90 module. [https://github.com/geoschem/geos-chem](https://github.com/geoschem/geos-chem)


SCIENCE CUSTOMIZATION GUIDE:
* Modify the algorithm in compute\_zhang (search for "TODO")
* Add scheme-specific helper subroutines as needed
* Update physical constants for your scheme
* Customize the environmental response functions




INFRASTRUCTURE RESPONSIBILITIES (handled by host model):
* Parameter initialization and validation
* Input array validation and error handling
* Memory management and array allocation
* Integration with host model time stepping




Generated on: 2025-11-13T17:12:59.281598 Author: Wei Li Reference: Zhang et al., 2001; Emerson et al., 2020 


    

------------------------------
The documentation for this class was generated from the following file `src/process/drydep/schemes/DryDepScheme_ZHANG_Mod.F90`

