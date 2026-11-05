

# Interface yaml\_interface\_mod::c\_yaml\_get\_integer\_array



[**ClassList**](annotated.md) **>** [**c\_yaml\_get\_integer\_array**](interfaceyaml__interface__mod_1_1c__yaml__get__integer__array.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  logical(c\_bool) function | [**c\_yaml\_get\_integer\_array**](#function-c_yaml_get_integer_array) (type(c\_ptr), value node, character(kind=c\_char), dimension(\*), intent(in) key, integer(c\_int), dimension(\*), intent(out) values, integer(c\_int), value max\_size, integer(c\_int), intent(out) actual\_size) <br> |




























## Public Functions Documentation




### function c\_yaml\_get\_integer\_array 

```Fortran
logical(c_bool) function c_yaml_get_integer_array::c_yaml_get_integer_array (
    type(c_ptr), value node,
    character(kind=c_char), dimension(*), intent(in) key,
    integer(c_int), dimension(*), intent(out) values,
    integer(c_int), value max_size,
    integer(c_int), intent(out) actual_size
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/external/yaml_interface/yaml_interface_mod.F90`

