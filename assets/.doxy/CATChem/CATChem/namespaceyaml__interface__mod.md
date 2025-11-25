

# Namespace yaml\_interface\_mod



[**Namespace List**](namespaces.md) **>** [**yaml\_interface\_mod**](namespaceyaml__interface__mod.md)




















## Classes

| Type | Name |
| ---: | :--- |
| interface | [**yaml\_get**](interfaceyaml__interface__mod_1_1yaml__get.md) <br>_Generic interface for getting values from YAML This allows uniform syntax: call yaml\_get(node, key, value, rc) for any supported data type._  |
| interface | [**yaml\_get\_array**](interfaceyaml__interface__mod_1_1yaml__get__array.md) <br>_Generic interface for getting arrays from YAML This allows uniform syntax: call yaml\_get\_array(node, key, values, rc) for any supported array type._  |
| interface | [**yaml\_set**](interfaceyaml__interface__mod_1_1yaml__set.md) <br>_Generic interface for setting values in YAML This allows uniform syntax: call yaml\_set(node, key, value, rc) for any supported data type._  |






















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**safe\_yaml\_get\_integer**](#function-safe_yaml_get_integer) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) yaml\_root, character(len=\*), intent(in) key, integer, intent(out) value, integer, intent(out) rc) <br>_Safe YAML integer reader Reads value as string first, then converts to avoid yaml-cpp conversion errors._  |
|  subroutine, public | [**safe\_yaml\_get\_logical**](#function-safe_yaml_get_logical) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) yaml\_root, character(len=\*), intent(in) key, logical, intent(out) value, integer, intent(out) rc) <br>_Safe YAML logical reader Reads value as string first, then converts to avoid yaml-cpp conversion errors._  |
|  subroutine, public | [**safe\_yaml\_get\_real**](#function-safe_yaml_get_real) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) yaml\_root, character(len=\*), intent(in) key, real(fp), intent(out) value, integer, intent(out) rc) <br>_Safe YAML real number reader Reads value as string first, then converts to avoid yaml-cpp conversion errors._  |
|  subroutine, public | [**yaml\_destroy\_node**](#function-yaml_destroy_node) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(inout) node) <br>_Destroy YAML node._  |
|  logical function, public | [**yaml\_get\_all\_keys**](#function-yaml_get_all_keys) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), dimension(:), intent(out) keys, integer, intent(out) actual\_count) <br>_Get all keys from a YAML map._  |
|  logical function, public | [**yaml\_get\_integer**](#function-yaml_get_integer) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, integer, intent(out) value) <br>_Get integer value._  |
|  logical function, public | [**yaml\_get\_integer\_array**](#function-yaml_get_integer_array) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, integer, dimension(:), intent(out) values, integer, intent(out) actual\_size) <br>_Get integer array._  |
|  logical function, public | [**yaml\_get\_logical**](#function-yaml_get_logical) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, logical, intent(out) value) <br>_Get logical value._  |
|  logical function, public | [**yaml\_get\_real**](#function-yaml_get_real) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, real(fp), intent(out) value) <br>_Get real value._  |
|  logical function, public | [**yaml\_get\_real\_array**](#function-yaml_get_real_array) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, real(fp), dimension(:), intent(out) values, integer, intent(out) actual\_size) <br>_Get real array._  |
|  integer function, public | [**yaml\_get\_size**](#function-yaml_get_size) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node) <br>_Get size of node._  |
|  logical function, public | [**yaml\_get\_string**](#function-yaml_get_string) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, character(len=\*), intent(out) value) <br>_Get string value._  |
|  logical function, public | [**yaml\_get\_string\_array**](#function-yaml_get_string_array) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, character(len=\*), dimension(:), intent(out) values, integer, intent(out) actual\_size) <br>_Get string array._  |
|  logical function, public | [**yaml\_has\_key**](#function-yaml_has_key) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key) <br>_Check if key exists._  |
|  logical function, public | [**yaml\_is\_map**](#function-yaml_is_map) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node) <br>_Check if node is map._  |
|  logical function, public | [**yaml\_is\_sequence**](#function-yaml_is_sequence) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node) <br>_Check if node is sequence._  |
|  type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)) function, public | [**yaml\_load\_file**](#function-yaml_load_file) (character(len=\*), intent(in) filename) <br>_Load YAML from file._  |
|  type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)) function, public | [**yaml\_load\_string**](#function-yaml_load_string) (character(len=\*), intent(in) yaml\_string) <br>_Load YAML from string._  |
|  logical function, public | [**yaml\_save\_file**](#function-yaml_save_file) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) filename) <br>_Save YAML to file._  |
|  logical function, public | [**yaml\_set\_integer**](#function-yaml_set_integer) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, integer, intent(in) value) <br>_Set integer value._  |
|  logical function, public | [**yaml\_set\_logical**](#function-yaml_set_logical) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, logical, intent(in) value) <br>_Set logical value._  |
|  logical function, public | [**yaml\_set\_real**](#function-yaml_set_real) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, real(fp), intent(in) value) <br>_Set real value._  |
|  logical function, public | [**yaml\_set\_string**](#function-yaml_set_string) (type([**yaml\_node\_t**](namespaceyaml__interface__mod.md#none-yaml_node_t)), intent(in) node, character(len=\*), intent(in) key, character(len=\*), intent(in) value) <br>_Set string value._  |




























## Public Functions Documentation




### function safe\_yaml\_get\_integer 

_Safe YAML integer reader Reads value as string first, then converts to avoid yaml-cpp conversion errors._ 
```Fortran
subroutine, public yaml_interface_mod::safe_yaml_get_integer (
    type( yaml_node_t ), intent(in) yaml_root,
    character(len=*), intent(in) key,
    integer, intent(out) value,
    integer, intent(out) rc
) 
```




<hr>



### function safe\_yaml\_get\_logical 

_Safe YAML logical reader Reads value as string first, then converts to avoid yaml-cpp conversion errors._ 
```Fortran
subroutine, public yaml_interface_mod::safe_yaml_get_logical (
    type( yaml_node_t ), intent(in) yaml_root,
    character(len=*), intent(in) key,
    logical, intent(out) value,
    integer, intent(out) rc
) 
```




<hr>



### function safe\_yaml\_get\_real 

_Safe YAML real number reader Reads value as string first, then converts to avoid yaml-cpp conversion errors._ 
```Fortran
subroutine, public yaml_interface_mod::safe_yaml_get_real (
    type( yaml_node_t ), intent(in) yaml_root,
    character(len=*), intent(in) key,
    real(fp), intent(out) value,
    integer, intent(out) rc
) 
```




<hr>



### function yaml\_destroy\_node 

_Destroy YAML node._ 
```Fortran
subroutine, public yaml_interface_mod::yaml_destroy_node (
    type( yaml_node_t ), intent(inout) node
) 
```




<hr>



### function yaml\_get\_all\_keys 

_Get all keys from a YAML map._ 
```Fortran
logical function, public yaml_interface_mod::yaml_get_all_keys (
    type( yaml_node_t ), intent(in) node,
    character(len=*), dimension(:), intent(out) keys,
    integer, intent(out) actual_count
) 
```




<hr>



### function yaml\_get\_integer 

_Get integer value._ 
```Fortran
logical function, public yaml_interface_mod::yaml_get_integer (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    integer, intent(out) value
) 
```




<hr>



### function yaml\_get\_integer\_array 

_Get integer array._ 
```Fortran
logical function, public yaml_interface_mod::yaml_get_integer_array (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    integer, dimension(:), intent(out) values,
    integer, intent(out) actual_size
) 
```




<hr>



### function yaml\_get\_logical 

_Get logical value._ 
```Fortran
logical function, public yaml_interface_mod::yaml_get_logical (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    logical, intent(out) value
) 
```




<hr>



### function yaml\_get\_real 

_Get real value._ 
```Fortran
logical function, public yaml_interface_mod::yaml_get_real (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    real(fp), intent(out) value
) 
```




<hr>



### function yaml\_get\_real\_array 

_Get real array._ 
```Fortran
logical function, public yaml_interface_mod::yaml_get_real_array (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    real(fp), dimension(:), intent(out) values,
    integer, intent(out) actual_size
) 
```




<hr>



### function yaml\_get\_size 

_Get size of node._ 
```Fortran
integer function, public yaml_interface_mod::yaml_get_size (
    type( yaml_node_t ), intent(in) node
) 
```




<hr>



### function yaml\_get\_string 

_Get string value._ 
```Fortran
logical function, public yaml_interface_mod::yaml_get_string (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    character(len=*), intent(out) value
) 
```




<hr>



### function yaml\_get\_string\_array 

_Get string array._ 
```Fortran
logical function, public yaml_interface_mod::yaml_get_string_array (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    character(len=*), dimension(:), intent(out) values,
    integer, intent(out) actual_size
) 
```




<hr>



### function yaml\_has\_key 

_Check if key exists._ 
```Fortran
logical function, public yaml_interface_mod::yaml_has_key (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key
) 
```




<hr>



### function yaml\_is\_map 

_Check if node is map._ 
```Fortran
logical function, public yaml_interface_mod::yaml_is_map (
    type( yaml_node_t ), intent(in) node
) 
```




<hr>



### function yaml\_is\_sequence 

_Check if node is sequence._ 
```Fortran
logical function, public yaml_interface_mod::yaml_is_sequence (
    type( yaml_node_t ), intent(in) node
) 
```




<hr>



### function yaml\_load\_file 

_Load YAML from file._ 
```Fortran
type( yaml_node_t ) function, public yaml_interface_mod::yaml_load_file (
    character(len=*), intent(in) filename
) 
```




<hr>



### function yaml\_load\_string 

_Load YAML from string._ 
```Fortran
type( yaml_node_t ) function, public yaml_interface_mod::yaml_load_string (
    character(len=*), intent(in) yaml_string
) 
```




<hr>



### function yaml\_save\_file 

_Save YAML to file._ 
```Fortran
logical function, public yaml_interface_mod::yaml_save_file (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) filename
) 
```




<hr>



### function yaml\_set\_integer 

_Set integer value._ 
```Fortran
logical function, public yaml_interface_mod::yaml_set_integer (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    integer, intent(in) value
) 
```




<hr>



### function yaml\_set\_logical 

_Set logical value._ 
```Fortran
logical function, public yaml_interface_mod::yaml_set_logical (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    logical, intent(in) value
) 
```




<hr>



### function yaml\_set\_real 

_Set real value._ 
```Fortran
logical function, public yaml_interface_mod::yaml_set_real (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    real(fp), intent(in) value
) 
```




<hr>



### function yaml\_set\_string 

_Set string value._ 
```Fortran
logical function, public yaml_interface_mod::yaml_set_string (
    type( yaml_node_t ), intent(in) node,
    character(len=*), intent(in) key,
    character(len=*), intent(in) value
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/external/yaml_interface/yaml_interface_mod.F90`

