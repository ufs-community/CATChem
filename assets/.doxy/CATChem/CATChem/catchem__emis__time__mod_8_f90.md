

# File catchem\_emis\_time\_mod.F90



[**FileList**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**catchem\_emis\_time\_mod.F90**](catchem__emis__time__mod_8_f90.md)

[Go to the source code of this file](catchem__emis__time__mod_8_f90_source.md)

_Pure, ESMF-free temporal-interpolation helpers for external emissions._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**catchem\_emis\_time**](namespacecatchem__emis__time.md) <br> |




















































## Detailed Description


These routines reproduce the MAPL ExtData / GOCART2G emission time handling (feature 014, US4) so the interpolation math can be unit-tested without pulling in the ESMF framework. The NUOPC emission driver (drivers/nuopc/catchem\_emis\_mod.F90) converts an ESMF\_Time into the plain (day-of-year, integer datetime-key) arguments below and delegates here.


Two behaviours live here:
* catchem\_emis\_time\_weight - weight between two record valid-times on a fixed 365-day day-of-year axis with cyclic year-end wrap (climatological anchoring).
* catchem\_emis\_file\_time\_bracket - select the lower bracketing record and the upper-record weight for a monthly climatology anchored at the file's own record timestamps ('file' anchoring), spanning vs cyclic regimes.
* catchem\_emis\_dt\_weight - daily-mean linear weight on 12Z knots, with the MAPL 'daily\_hold' piecewise-constant variant.




All routines are pure and depend only on the bridge precision kind. 


    

------------------------------
The documentation for this class was generated from the following file `src/process/catchem_emis_time_mod.F90`

