

# Namespace catchem\_highlevel\_api



[**Namespace List**](namespaces.md) **>** [**catchem\_highlevel\_api**](namespacecatchem__highlevel__api.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**catchem\_add\_process**](#function-catchem_add_process) (type([**catchemmodeltype**](namespacecatchem__highlevel__api.md#none-catchemmodeltype)), intent(inout) model, character(len=\*), intent(in) process\_name, character(len=\*), intent(in), optional process\_config, integer, intent(out) rc) <br>_Add a process to CATChem._  |
|  subroutine, public | [**catchem\_finalize**](#function-catchem_finalize) (type([**catchemmodeltype**](namespacecatchem__highlevel__api.md#none-catchemmodeltype)), intent(inout) model, integer, intent(out) rc) <br>_Finalize CATChem._  |
|  subroutine, public | [**catchem\_get\_concentrations**](#function-catchem_get_concentrations) (type([**catchemmodeltype**](namespacecatchem__highlevel__api.md#none-catchemmodeltype)), intent(in) model, real(fp), dimension(:,:,:,:), intent(out), allocatable concentrations, character(len=32), dimension(:), intent(out), allocatable species\_names, integer, intent(out) rc) <br>_Get chemical concentrations._  |
|  subroutine, public | [**catchem\_get\_diagnostics**](#function-catchem_get_diagnostics) (type([**catchemmodeltype**](namespacecatchem__highlevel__api.md#none-catchemmodeltype)), intent(in) model, type([**catchemdiagnosticstype**](namespacecatchem__highlevel__api.md#none-catchemdiagnosticstype)), intent(out) diagnostics, integer, intent(out) rc) <br>_Get diagnostic data._  |
|  subroutine, public | [**catchem\_init**](#function-catchem_init) (type([**catchemmodeltype**](namespacecatchem__highlevel__api.md#none-catchemmodeltype)), intent(out) model, type([**catchemconfigtype**](namespacecatchem__highlevel__api.md#none-catchemconfigtype)), intent(in) config, integer, intent(out) rc) <br>_Initialize CATChem with simple configuration._  |
|  subroutine, public | [**catchem\_run\_timestep**](#function-catchem_run_timestep) (type([**catchemmodeltype**](namespacecatchem__highlevel__api.md#none-catchemmodeltype)), intent(inout) model, type([**catchemdatatype**](namespacecatchem__highlevel__api.md#none-catchemdatatype)), intent(in) met\_data, type([**catchemdatatype**](namespacecatchem__highlevel__api.md#none-catchemdatatype)), intent(in) emis\_data, integer, intent(out) rc) <br>_Run a single CATChem timestep._  |
|  subroutine, public | [**catchem\_set\_emissions**](#function-catchem_set_emissions) (type([**catchemmodeltype**](namespacecatchem__highlevel__api.md#none-catchemmodeltype)), intent(inout) model, type([**catchemdatatype**](namespacecatchem__highlevel__api.md#none-catchemdatatype)), intent(in) emis\_data, integer, intent(out) rc) <br>_Set emission data._  |
|  subroutine, public | [**catchem\_set\_meteorology**](#function-catchem_set_meteorology) (type([**catchemmodeltype**](namespacecatchem__highlevel__api.md#none-catchemmodeltype)), intent(inout) model, type([**catchemdatatype**](namespacecatchem__highlevel__api.md#none-catchemdatatype)), intent(in) met\_data, integer, intent(out) rc) <br>_Set meteorological data._  |




























## Public Functions Documentation




### function catchem\_add\_process 

_Add a process to CATChem._ 
```Fortran
subroutine, public catchem_highlevel_api::catchem_add_process (
    type( catchemmodeltype ), intent(inout) model,
    character(len=*), intent(in) process_name,
    character(len=*), intent(in), optional process_config,
    integer, intent(out) rc
) 
```




<hr>



### function catchem\_finalize 

_Finalize CATChem._ 
```Fortran
subroutine, public catchem_highlevel_api::catchem_finalize (
    type( catchemmodeltype ), intent(inout) model,
    integer, intent(out) rc
) 
```




<hr>



### function catchem\_get\_concentrations 

_Get chemical concentrations._ 
```Fortran
subroutine, public catchem_highlevel_api::catchem_get_concentrations (
    type( catchemmodeltype ), intent(in) model,
    real(fp), dimension(:,:,:,:), intent(out), allocatable concentrations,
    character(len=32), dimension(:), intent(out), allocatable species_names,
    integer, intent(out) rc
) 
```




<hr>



### function catchem\_get\_diagnostics 

_Get diagnostic data._ 
```Fortran
subroutine, public catchem_highlevel_api::catchem_get_diagnostics (
    type( catchemmodeltype ), intent(in) model,
    type( catchemdiagnosticstype ), intent(out) diagnostics,
    integer, intent(out) rc
) 
```




<hr>



### function catchem\_init 

_Initialize CATChem with simple configuration._ 
```Fortran
subroutine, public catchem_highlevel_api::catchem_init (
    type( catchemmodeltype ), intent(out) model,
    type( catchemconfigtype ), intent(in) config,
    integer, intent(out) rc
) 
```




<hr>



### function catchem\_run\_timestep 

_Run a single CATChem timestep._ 
```Fortran
subroutine, public catchem_highlevel_api::catchem_run_timestep (
    type( catchemmodeltype ), intent(inout) model,
    type( catchemdatatype ), intent(in) met_data,
    type( catchemdatatype ), intent(in) emis_data,
    integer, intent(out) rc
) 
```




<hr>



### function catchem\_set\_emissions 

_Set emission data._ 
```Fortran
subroutine, public catchem_highlevel_api::catchem_set_emissions (
    type( catchemmodeltype ), intent(inout) model,
    type( catchemdatatype ), intent(in) emis_data,
    integer, intent(out) rc
) 
```




<hr>



### function catchem\_set\_meteorology 

_Set meteorological data._ 
```Fortran
subroutine, public catchem_highlevel_api::catchem_set_meteorology (
    type( catchemmodeltype ), intent(inout) model,
    type( catchemdatatype ), intent(in) met_data,
    integer, intent(out) rc
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/api/CATChem_HighLevel_API.F90`

