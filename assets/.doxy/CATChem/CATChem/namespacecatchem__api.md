

# Namespace catchem\_api



[**Namespace List**](namespaces.md) **>** [**catchem\_api**](namespacecatchem__api.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**model\_initialize**](#function-model_initialize) (class([**catchem\_model**](namespacecatchem__api.md#none-catchem_model)), intent(inout) this, character(len=\*), intent(in) config\_file, integer, intent(in) nx, integer, intent(in) ny, integer, intent(in) nz, integer, intent(in), optional nsoil, integer, intent(in), optional nsoiltype, integer, intent(in), optional nsurftype, integer, intent(out) rc) <br>_Initialize the CATChem model with configuration file and grid dimensions This method sets up the core CATChem infrastructure using the builder pattern, loads configuration from the specified file, and sets up the grid geometry._  |




























## Public Functions Documentation




### function model\_initialize 

_Initialize the CATChem model with configuration file and grid dimensions This method sets up the core CATChem infrastructure using the builder pattern, loads configuration from the specified file, and sets up the grid geometry._ 
```Fortran
subroutine catchem_api::model_initialize (
    class( catchem_model ), intent(inout) this,
    character(len=*), intent(in) config_file,
    integer, intent(in) nx,
    integer, intent(in) ny,
    integer, intent(in) nz,
    integer, intent(in), optional nsoil,
    integer, intent(in), optional nsoiltype,
    integer, intent(in), optional nsurftype,
    integer, intent(out) rc
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/api/CATChem_API.F90`

