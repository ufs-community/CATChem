

# Interface catchem\_api::catchem\_diag\_get\_contract



[**ClassList**](annotated.md) **>** [**catchem\_api**](namespacecatchem__api.md) **>** [**catchem\_diag\_get\_contract**](interfacecatchem__api_1_1catchem__diag__get__contract.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  integer(c\_int) function | [**catchem\_diag\_get\_contract**](#function-catchem_diag_get_contract) (type(c\_ptr), value core\_ptr, character(kind=c\_char), dimension(\*), intent(in) name, integer(c\_int), intent(out) generation, integer(c\_int), intent(out) availability, integer(c\_int), intent(out) latest\_writer, integer(c\_int), intent(out) policy) <br> |




























## Public Functions Documentation




### function catchem\_diag\_get\_contract 

```Fortran
integer(c_int) function catchem_api::catchem_diag_get_contract::catchem_diag_get_contract (
    type(c_ptr), value core_ptr,
    character(kind=c_char), dimension(*), intent(in) name,
    integer(c_int), intent(out) generation,
    integer(c_int), intent(out) availability,
    integer(c_int), intent(out) latest_writer,
    integer(c_int), intent(out) policy
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/api/CATChem_API.F90`

