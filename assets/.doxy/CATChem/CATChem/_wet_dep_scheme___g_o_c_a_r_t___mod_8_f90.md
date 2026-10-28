

# File WetDepScheme\_GOCART\_Mod.F90



[**FileList**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**wetdep**](dir_8b9a0ce556ea4a65f6920dfb49dcd69d.md) **>** [**schemes**](dir_8ca87c5e2f5cf830ab1a41055168a46b.md) **>** [**WetDepScheme\_GOCART\_Mod.F90**](_wet_dep_scheme___g_o_c_a_r_t___mod_8_f90.md)

[Go to the source code of this file](_wet_dep_scheme___g_o_c_a_r_t___mod_8_f90_source.md)

_GOCART2G wet removal scheme: SU\_Wet\_Removal for sulfate species (DMS/SO2/SO4/MSA) and WetRemovalUFS for all other species._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**wetdepscheme\_gocart\_mod**](namespacewetdepscheme__gocart__mod.md) <br> |




















































## Detailed Description


Pure science kernel for gocart scheme in wetdep process. This module wraps the GOCART2G process-library wet-removal routines. Rather than re-implementing the science, it prepares the column meteorology in the layout the GOCART routines expect (top-to-bottom ordering, (1,1,k) shaped pointers) and calls them directly:
* SU\_Wet\_Removal : the sulfate group (DMS/SO2/SO4/MSA) in a single call.
* WetRemovalUFS : every other aerosol species, one call per species.




Ported from ufs-community/CATChem PR #206 (feature/standalone) for the C++ core. Two deliberate adaptations relative to the PR source:
* The PR tree ran on a pre-filtered wetdep species list, so participation was implicit. The C++ core passes the FULL mechanism species list with an explicit species\_is\_wetdep flag (same contract as the Jacob kernel), so this kernel gates every removal route on that flag: non-participants are echoed bit-exactly into species\_tendencies.
* Module imports use the rework's bridge precision/constants utilities instead of the legacy precision\_mod/Constants.




Author: Wei Li Reference: GOCART2G Process Library (GEOS-ESM/GOCART): SU\_Wet\_Removal and WetRemovalUFS (Jacob et al. [2000] Harvard scheme) 


    

------------------------------
The documentation for this class was generated from the following file `src/process/wetdep/schemes/WetDepScheme_GOCART_Mod.F90`

