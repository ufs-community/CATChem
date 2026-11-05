

# Interface yaml\_interface\_mod::yaml\_get



[**ClassList**](annotated.md) **>** [**yaml\_interface\_mod**](namespaceyaml__interface__mod.md) **>** [**yaml\_get**](interfaceyaml__interface__mod_1_1yaml__get.md)



_Generic interface for getting values from YAML This allows uniform syntax: call yaml\_get(node, key, value, rc) for any supported data type._ 






































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**yaml\_get\_integer\_generic**](#function-yaml_get_integer_generic) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, integer, intent(out) value, integer, intent(out), optional rc, integer, intent(in), optional default\_value) <br>_Generic integer getter with optional default and return code._  |
|  subroutine | [**yaml\_get\_logical\_generic**](#function-yaml_get_logical_generic) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, logical, intent(out) value, integer, intent(out), optional rc, logical, intent(in), optional default\_value) <br>_Generic logical getter with optional default and return code._  |
|  subroutine | [**yaml\_get\_real\_dp\_generic**](#function-yaml_get_real_dp_generic) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, real(kind=real64), intent(out) value, integer, intent(out), optional rc, real(kind=real64), intent(in), optional default\_value) <br>_Generic double precision real getter._  |
|  subroutine | [**yaml\_get\_real\_sp\_generic**](#function-yaml_get_real_sp_generic) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, real(kind=real32), intent(out) value, integer, intent(out), optional rc, real(kind=real32), intent(in), optional default\_value) <br>_Generic single precision real getter._  |
|  subroutine | [**yaml\_get\_string\_generic**](#function-yaml_get_string_generic) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, character(len=\*), intent(out) value, integer, intent(out), optional rc, character(len=\*), intent(in), optional default\_value) <br>_Generic string getter with optional default and return code._  |




























## Public Functions Documentation




### function yaml\_get\_integer\_generic 

_Generic integer getter with optional default and return code._ 
```Fortran
subroutine yaml_interface_mod::yaml_get::yaml_get_integer_generic (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    integer, intent(out) value,
    integer, intent(out), optional rc,
    integer, intent(in), optional default_value
) 
```




<hr>



### function yaml\_get\_logical\_generic 

_Generic logical getter with optional default and return code._ 
```Fortran
subroutine yaml_interface_mod::yaml_get::yaml_get_logical_generic (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    logical, intent(out) value,
    integer, intent(out), optional rc,
    logical, intent(in), optional default_value
) 
```




<hr>



### function yaml\_get\_real\_dp\_generic 

_Generic double precision real getter._ 
```Fortran
subroutine yaml_interface_mod::yaml_get::yaml_get_real_dp_generic (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    real(kind=real64), intent(out) value,
    integer, intent(out), optional rc,
    real(kind=real64), intent(in), optional default_value
) 
```




<hr>



### function yaml\_get\_real\_sp\_generic 

_Generic single precision real getter._ 
```Fortran
subroutine yaml_interface_mod::yaml_get::yaml_get_real_sp_generic (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    real(kind=real32), intent(out) value,
    integer, intent(out), optional rc,
    real(kind=real32), intent(in), optional default_value
) 
```




<hr>



### function yaml\_get\_string\_generic 

_Generic string getter with optional default and return code._ 
```Fortran
subroutine yaml_interface_mod::yaml_get::yaml_get_string_generic (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    character(len=*), intent(out) value,
    integer, intent(out), optional rc,
    character(len=*), intent(in), optional default_value
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/external/yaml_interface/yaml_interface_mod.F90`

