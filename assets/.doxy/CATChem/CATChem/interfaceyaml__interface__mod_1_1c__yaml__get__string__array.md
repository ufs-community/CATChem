

# Interface yaml\_interface\_mod::c\_yaml\_get\_string\_array



[**ClassList**](annotated.md) **>** [**c\_yaml\_get\_string\_array**](interfaceyaml__interface__mod_1_1c__yaml__get__string__array.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  logical(c\_bool) function | [**c\_yaml\_get\_string\_array**](#function-c_yaml_get_string_array) (type(c\_ptr), value node, character(kind=c\_char), dimension(\*), intent(in) key, character(kind=c\_char), dimension(\*), intent(out) values, integer(c\_int), value max\_strings, integer(c\_int), value max\_len, integer(c\_int), intent(out) actual\_size) <br> |




























## Public Functions Documentation




### function c\_yaml\_get\_string\_array 

```Fortran
logical(c_bool) function c_yaml_get_string_array::c_yaml_get_string_array (
    type(c_ptr), value node,
    character(kind=c_char), dimension(*), intent(in) key,
    character(kind=c_char), dimension(*), intent(out) values,
    integer(c_int), value max_strings,
    integer(c_int), value max_len,
    integer(c_int), intent(out) actual_size
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/external/yaml_interface/yaml_interface_mod.F90`

