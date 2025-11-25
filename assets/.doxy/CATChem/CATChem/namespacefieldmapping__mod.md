

# Namespace fieldmapping\_mod



[**Namespace List**](namespaces.md) **>** [**fieldmapping\_mod**](namespacefieldmapping__mod.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  integer, parameter, public | [**mapping\_failure**](#variable-mapping_failure)   = `-1`<br> |
|  integer, parameter, public | [**mapping\_field\_not\_found**](#variable-mapping_field_not_found)   = `-2`<br> |
|  integer, parameter, public | [**mapping\_invalid\_category**](#variable-mapping_invalid_category)   = `-3`<br> |
|  integer, parameter, public | [**mapping\_success**](#variable-mapping_success)   = `0`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**mapping\_init**](#function-mapping_init) (class([**fieldmappingtype**](namespacefieldmapping__mod.md#none-fieldmappingtype)), intent(inout) this) <br>_Initialize the field mapping system._  |




























## Public Attributes Documentation




### variable mapping\_failure 

```Fortran
integer, parameter, public fieldmapping_mod::mapping_failure;
```




<hr>



### variable mapping\_field\_not\_found 

```Fortran
integer, parameter, public fieldmapping_mod::mapping_field_not_found;
```




<hr>



### variable mapping\_invalid\_category 

```Fortran
integer, parameter, public fieldmapping_mod::mapping_invalid_category;
```




<hr>



### variable mapping\_success 

```Fortran
integer, parameter, public fieldmapping_mod::mapping_success;
```




<hr>
## Public Functions Documentation




### function mapping\_init 

_Initialize the field mapping system._ 
```Fortran
subroutine fieldmapping_mod::mapping_init (
    class( fieldmappingtype ), intent(inout) this
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/api/FieldMapping_Mod.F90`

