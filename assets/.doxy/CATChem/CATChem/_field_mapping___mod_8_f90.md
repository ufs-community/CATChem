

# File FieldMapping\_Mod.F90



[**FileList**](files.md) **>** [**api**](dir_da61e3e9a357748887e3ca8d7c5a0c16.md) **>** [**FieldMapping\_Mod.F90**](_field_mapping___mod_8_f90.md)

[Go to the source code of this file](_field_mapping___mod_8_f90_source.md)

_Field mapping system for CATChem high-level API._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**fieldmapping\_mod**](namespacefieldmapping__mod.md) <br> |




















































## Detailed Description




**Author:**

CATChem Development Team 




**Date:**

2025 




**Version:**

2.0


This module provides a flexible field mapping system that allows host models to map their field names to CATChem field names. This enables seamless integration without requiring host models to use CATChem's internal naming conventions.


Example usage: 
```Fortran
type(FieldMappingType) :: mapper
call mapper%init()

! Map host model fields to CATChem fields
call mapper%add_mapping('host_temp', 'temperature', 'meteo', rc)
call mapper%add_mapping('host_pres', 'pressure', 'meteo', rc)
call mapper%add_mapping('host_o3', 'O3', 'chemistry', rc)

! Use mapping to set data
call catchem%set_field_by_mapping(mapper, 'host_temp', temp_data, rc)
```
 


    

------------------------------
The documentation for this class was generated from the following file `src/api/FieldMapping_Mod.F90`

