

# File CATChem\_API.F90



[**FileList**](files.md) **>** [**api**](dir_da61e3e9a357748887e3ca8d7c5a0c16.md) **>** [**CATChem\_API.F90**](_c_a_t_chem___a_p_i_8_f90.md)

[Go to the source code of this file](_c_a_t_chem___a_p_i_8_f90_source.md)

_Modernized high-level BIND(C) API for host model integration._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**catchem\_api**](namespacecatchem__api.md) <br> |
| namespace | [**iso\_c\_binding**](namespaceiso__c__binding.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| interface | [**catchem\_diag\_get\_contract**](interfacecatchem__api_1_1catchem__diag__get__contract.md) <br> |
| interface | [**catchem\_diag\_register\_contract\_checked**](interfacecatchem__api_1_1catchem__diag__register__contract__checked.md) <br> |


















































## Detailed Description


This module provides the standard, backward-compatible CATChem\_Model derived type used by Earth System drivers (like NUOPC and UFS). It delegates all memory management, synchronization, and process scheduling directly to the modernized top-down C++ Core via standard BIND(C) bridges. 


    

------------------------------
The documentation for this class was generated from the following file `src/api/CATChem_API.F90`

