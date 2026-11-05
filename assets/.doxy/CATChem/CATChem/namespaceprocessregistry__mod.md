

# Namespace processregistry\_mod



[**Namespace List**](namespaces.md) **>** [**processregistry\_mod**](namespaceprocessregistry__mod.md)




















## Classes

| Type | Name |
| ---: | :--- |
| interface | [**processcreatorinterface**](interfaceprocessregistry__mod_1_1processcreatorinterface.md) <br>_Function pointer interface for process creators._  |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  type([**processregistrytype**](namespaceprocessregistry__mod.md#none-processregistrytype)), target, save | [**global\_registry**](#variable-global_registry)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  type([**processregistrytype**](namespaceprocessregistry__mod.md#none-processregistrytype)) function, pointer, public | [**get\_global\_registry**](#function-get_global_registry) () <br>_Get the global process registry instance._  |




























## Public Attributes Documentation




### variable global\_registry 

```Fortran
type(processregistrytype), target, save processregistry_mod::global_registry;
```




<hr>
## Public Functions Documentation




### function get\_global\_registry 

_Get the global process registry instance._ 
```Fortran
type( processregistrytype ) function, pointer, public processregistry_mod::get_global_registry () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/ProcessRegistry_Mod.F90`

