

# Namespace statemanager\_mod



[**Namespace List**](namespaces.md) **>** [**statemanager\_mod**](namespacestatemanager__mod.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  integer, parameter, public | [**state\_status\_error**](#variable-state_status_error)   = `-1`<br> |
|  integer, parameter, public | [**state\_status\_initialized**](#variable-state_status_initialized)   = `1`<br> |
|  integer, parameter, public | [**state\_status\_uninitialized**](#variable-state_status_uninitialized)   = `0`<br>_State status enumeration._  |
|  integer, parameter, public | [**state\_status\_valid**](#variable-state_status_valid)   = `2`<br> |
|  integer, parameter, public | [**state\_type\_chem**](#variable-state_type_chem)   = `2`<br> |
|  integer, parameter, public | [**state\_type\_diag**](#variable-state_type_diag)   = `4`<br> |
|  integer, parameter, public | [**state\_type\_emis**](#variable-state_type_emis)   = `3`<br> |
|  integer, parameter, public | [**state\_type\_met**](#variable-state_type_met)   = `1`<br>_State type enumeration._  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**allocate\_met\_field**](#function-allocate_met_field) (class([**metstatetype**](namespacemetstate__mod.md#none-metstatetype)), intent(inout) met\_state, character(len=\*), intent(in) field\_name, integer, intent(out) rc) <br>_Allocate a specific field in a MetStateType object._  |
|  character(len=32) function, public | [**get\_state\_type\_name**](#function-get_state_type_name) (integer, intent(in) state\_type) <br>_Get state type name from enumeration._  |
|  subroutine | [**manager\_init**](#function-manager_init) (class([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) this, character(len=\*), intent(in), optional name, integer, intent(out) rc) <br>_Initialize the state manager (called by CATChemCore)_  |




























## Public Attributes Documentation




### variable state\_status\_error 

```Fortran
integer, parameter, public statemanager_mod::state_status_error;
```




<hr>



### variable state\_status\_initialized 

```Fortran
integer, parameter, public statemanager_mod::state_status_initialized;
```




<hr>



### variable state\_status\_uninitialized 

_State status enumeration._ 
```Fortran
integer, parameter, public statemanager_mod::state_status_uninitialized;
```




<hr>



### variable state\_status\_valid 

```Fortran
integer, parameter, public statemanager_mod::state_status_valid;
```




<hr>



### variable state\_type\_chem 

```Fortran
integer, parameter, public statemanager_mod::state_type_chem;
```




<hr>



### variable state\_type\_diag 

```Fortran
integer, parameter, public statemanager_mod::state_type_diag;
```




<hr>



### variable state\_type\_emis 

```Fortran
integer, parameter, public statemanager_mod::state_type_emis;
```




<hr>



### variable state\_type\_met 

_State type enumeration._ 
```Fortran
integer, parameter, public statemanager_mod::state_type_met;
```




<hr>
## Public Functions Documentation




### function allocate\_met\_field 

_Allocate a specific field in a MetStateType object._ 
```Fortran
subroutine, public statemanager_mod::allocate_met_field (
    class( metstatetype ), intent(inout) met_state,
    character(len=*), intent(in) field_name,
    integer, intent(out) rc
) 
```



Direct utility function for field allocation in MetStateType. For more complex state management, use CATChemCore\_Mod.




**Parameters:**


* `met_state` MetStateType object 
* `field_name` Name of the field to allocate (or 'ALL') 
* `rc` Return code (CC\_SUCCESS or error code) 




        

<hr>



### function get\_state\_type\_name 

_Get state type name from enumeration._ 
```Fortran
character(len=32) function, public statemanager_mod::get_state_type_name (
    integer, intent(in) state_type
) 
```




<hr>



### function manager\_init 

_Initialize the state manager (called by CATChemCore)_ 
```Fortran
subroutine statemanager_mod::manager_init (
    class( statemanagertype ), intent(inout) this,
    character(len=*), intent(in), optional name,
    integer, intent(out) rc
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/StateManager_Mod.F90`

