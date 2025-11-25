

# File yaml\_interface\_mod.F90



[**FileList**](files.md) **>** [**external**](dir_805a0af995e93a362739e98abd740eb2.md) **>** [**yaml\_interface**](dir_d0b1a67acd809cff502adc02c61e9ebd.md) **>** [**yaml\_interface\_mod.F90**](yaml__interface__mod_8_f90.md)

[Go to the source code of this file](yaml__interface__mod_8_f90_source.md)

_High-level Fortran interface for yaml-cpp._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**iso\_c\_binding**](namespaceiso__c__binding.md) <br> |
| namespace | [**yaml\_interface\_mod**](namespaceyaml__interface__mod.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| interface | [**yaml\_get**](interfaceyaml__interface__mod_1_1yaml__get.md) <br>_Generic interface for getting values from YAML This allows uniform syntax: call yaml\_get(node, key, value, rc) for any supported data type._  |
| interface | [**yaml\_get\_array**](interfaceyaml__interface__mod_1_1yaml__get__array.md) <br>_Generic interface for getting arrays from YAML This allows uniform syntax: call yaml\_get\_array(node, key, values, rc) for any supported array type._  |
| interface | [**yaml\_set**](interfaceyaml__interface__mod_1_1yaml__set.md) <br>_Generic interface for setting values in YAML This allows uniform syntax: call yaml\_set(node, key, value, rc) for any supported data type._  |


















































## Detailed Description




**Author:**

CATChem Development Team 




**Date:**

2025


This module provides a high-level Fortran interface to yaml-cpp with generic procedures that can handle all data types seamlessly. 


    

------------------------------
The documentation for this class was generated from the following file `src/external/yaml_interface/yaml_interface_mod.F90`

