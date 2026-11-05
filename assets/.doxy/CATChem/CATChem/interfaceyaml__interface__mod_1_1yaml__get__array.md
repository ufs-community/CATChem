

# Interface yaml\_interface\_mod::yaml\_get\_array



[**ClassList**](annotated.md) **>** [**yaml\_interface\_mod**](namespaceyaml__interface__mod.md) **>** [**yaml\_get\_array**](interfaceyaml__interface__mod_1_1yaml__get__array.md)



_Generic interface for getting arrays from YAML This allows uniform syntax: call yaml\_get\_array(node, key, values, rc) for any supported array type._ 






































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**yaml\_get\_integer\_array\_generic**](#function-yaml_get_integer_array_generic) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, integer, dimension(:), intent(out) values, integer, intent(out), optional rc, integer, intent(out), optional actual\_size) <br>_Generic integer array getter._  |
|  subroutine | [**yaml\_get\_real\_dp\_array\_generic**](#function-yaml_get_real_dp_array_generic) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, real(kind=fp), dimension(:), intent(out) values, integer, intent(out), optional rc, integer, intent(out), optional actual\_size) <br>_Generic double precision real array getter._  |
|  subroutine | [**yaml\_get\_string\_array\_generic**](#function-yaml_get_string_array_generic) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, character(len=\*), dimension(:), intent(out) values, integer, intent(out), optional rc, integer, intent(out), optional actual\_size) <br>_Generic string array getter._  |




























## Public Functions Documentation




### function yaml\_get\_integer\_array\_generic 

_Generic integer array getter._ 
```Fortran
subroutine yaml_interface_mod::yaml_get_array::yaml_get_integer_array_generic (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    integer, dimension(:), intent(out) values,
    integer, intent(out), optional rc,
    integer, intent(out), optional actual_size
) 
```




<hr>



### function yaml\_get\_real\_dp\_array\_generic 

_Generic double precision real array getter._ 
```Fortran
subroutine yaml_interface_mod::yaml_get_array::yaml_get_real_dp_array_generic (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    real(kind=fp), dimension(:), intent(out) values,
    integer, intent(out), optional rc,
    integer, intent(out), optional actual_size
) 
```




<hr>



### function yaml\_get\_string\_array\_generic 

_Generic string array getter._ 
```Fortran
subroutine yaml_interface_mod::yaml_get_array::yaml_get_string_array_generic (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    character(len=*), dimension(:), intent(out) values,
    integer, intent(out), optional rc,
    integer, intent(out), optional actual_size
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/external/yaml_interface/yaml_interface_mod.F90`

