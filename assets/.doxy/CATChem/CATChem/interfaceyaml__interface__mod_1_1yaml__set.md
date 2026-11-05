

# Interface yaml\_interface\_mod::yaml\_set



[**ClassList**](annotated.md) **>** [**yaml\_interface\_mod**](namespaceyaml__interface__mod.md) **>** [**yaml\_set**](interfaceyaml__interface__mod_1_1yaml__set.md)



_Generic interface for setting values in YAML This allows uniform syntax: call yaml\_set(node, key, value, rc) for any supported data type._ 






































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**yaml\_set\_integer\_generic**](#function-yaml_set_integer_generic) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, integer, intent(in) value, integer, intent(out), optional rc) <br>_Generic integer setter with return code._  |
|  subroutine | [**yaml\_set\_logical\_generic**](#function-yaml_set_logical_generic) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, logical, intent(in) value, integer, intent(out), optional rc) <br>_Generic logical setter with return code._  |
|  subroutine | [**yaml\_set\_real\_dp\_generic**](#function-yaml_set_real_dp_generic) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, real(kind=fp), intent(in) value, integer, intent(out), optional rc) <br>_Generic double precision real setter._  |
|  subroutine | [**yaml\_set\_string\_generic**](#function-yaml_set_string_generic) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, character(len=\*), intent(in) value, integer, intent(out), optional rc) <br>_Generic string setter with return code._  |




























## Public Functions Documentation




### function yaml\_set\_integer\_generic 

_Generic integer setter with return code._ 
```Fortran
subroutine yaml_interface_mod::yaml_set::yaml_set_integer_generic (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    integer, intent(in) value,
    integer, intent(out), optional rc
) 
```




<hr>



### function yaml\_set\_logical\_generic 

_Generic logical setter with return code._ 
```Fortran
subroutine yaml_interface_mod::yaml_set::yaml_set_logical_generic (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    logical, intent(in) value,
    integer, intent(out), optional rc
) 
```




<hr>



### function yaml\_set\_real\_dp\_generic 

_Generic double precision real setter._ 
```Fortran
subroutine yaml_interface_mod::yaml_set::yaml_set_real_dp_generic (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    real(kind=fp), intent(in) value,
    integer, intent(out), optional rc
) 
```




<hr>



### function yaml\_set\_string\_generic 

_Generic string setter with return code._ 
```Fortran
subroutine yaml_interface_mod::yaml_set::yaml_set_string_generic (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    character(len=*), intent(in) value,
    integer, intent(out), optional rc
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/external/yaml_interface/yaml_interface_mod.F90`

