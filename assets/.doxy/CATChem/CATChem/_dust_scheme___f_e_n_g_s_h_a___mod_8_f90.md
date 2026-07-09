

# File DustScheme\_FENGSHA\_Mod.F90



[**FileList**](files.md) **>** [**dust**](dir_1c14dfbaca1e3f4c2e26e74290119ebd.md) **>** [**schemes**](dir_11b2254edcf6ee5df673de29b129f986.md) **>** [**DustScheme\_FENGSHA\_Mod.F90**](_dust_scheme___f_e_n_g_s_h_a___mod_8_f90.md)

[Go to the source code of this file](_dust_scheme___f_e_n_g_s_h_a___mod_8_f90_source.md)

_Fengsha Dust emission scheme developed at NOAA ARL for use at NOAA NWS._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**dustscheme\_fengsha\_mod**](namespacedustscheme__fengsha__mod.md) <br> |




















































## Detailed Description


Pure science kernel for fengsha scheme in dust process. This module contains ONLY the computational algorithm with NO infrastructure dependencies. Uses only basic Fortran types for maximum portability and reusability.


SCIENCE CUSTOMIZATION GUIDE:
* Modify the algorithm in compute\_fengsha (search for "TODO")
* Add scheme-specific helper subroutines as needed
* Update physical constants for your scheme
* Customize the environmental response functions




INFRASTRUCTURE RESPONSIBILITIES (handled by host model):
* Parameter initialization and validation
* Input array validation and error handling
* Memory management and array allocation
* Integration with host model time stepping




Generated on: 2026-04-17T13:57:10.251187 Author: Barry Baker & Wei Li Reference: Zhang et al. 2022 


    

------------------------------
The documentation for this class was generated from the following file `src/process/dust/schemes/DustScheme_FENGSHA_Mod.F90`

