

# File WetDepScheme\_JACOB\_Mod.F90



[**FileList**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**wetdep**](dir_8b9a0ce556ea4a65f6920dfb49dcd69d.md) **>** [**schemes**](dir_8ca87c5e2f5cf830ab1a41055168a46b.md) **>** [**WetDepScheme\_JACOB\_Mod.F90**](_wet_dep_scheme___j_a_c_o_b___mod_8_f90.md)

[Go to the source code of this file](_wet_dep_scheme___j_a_c_o_b___mod_8_f90_source.md)

_Jacob et al._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**wetdepscheme\_jacob\_mod**](namespacewetdepscheme__jacob__mod.md) <br> |




















































## Detailed Description


[2000] wet deposition scheme


Pure science kernel for jacob scheme in wetdep process. This module contains ONLY the computational algorithm with NO infrastructure dependencies. Uses only basic Fortran types for maximum portability and reusability.


SCIENCE CUSTOMIZATION GUIDE:
* Modify the algorithm in compute\_jacob (search for "TODO")
* Add scheme-specific helper subroutines as needed
* Update physical constants for your scheme
* Customize the environmental response functions




INFRASTRUCTURE RESPONSIBILITIES (handled by host model):
* Parameter initialization and validation
* Input array validation and error handling
* Memory management and array allocation
* Integration with host model time stepping




Generated on: 2025-11-25T14:47:14.566559 Author: Wei Li (1) Jacob, D. J., Liu, H., Mari, C., and Yantosca, B. M., Harvard wet deposition scheme for GMI, available at: [http://acmg.seas.harvard.edu/geos/wiki\_docs/deposition/wetdep.jacob\_etal\_2000.pdf](http://acmg.seas.harvard.edu/geos/wiki_docs/deposition/wetdep.jacob_etal_2000.pdf) (2) GEOS-Chem's source codes in the module file of wetscav\_mod.F90 and reference therein. ([https://github.com/geoschem/geos-chem/blob/main/GeosCore/wetscav\_mod.F90](https://github.com/geoschem/geos-chem/blob/main/GeosCore/wetscav_mod.F90)) (3) The above scheme was also adopted in GOCART2G\_process.F90 for aerosols, which is shorter and cleaner. [https://github.com/GEOS-ESM/GOCART/blob/develop/Process\_Library/GOCART2G\_Process.F90#L3525-L4115](https://github.com/GEOS-ESM/GOCART/blob/develop/Process_Library/GOCART2G_Process.F90#L3525-L4115) 


    

------------------------------
The documentation for this class was generated from the following file `src/process/wetdep/schemes/WetDepScheme_JACOB_Mod.F90`

