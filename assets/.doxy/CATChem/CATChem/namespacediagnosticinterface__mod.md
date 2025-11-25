

# Namespace diagnosticinterface\_mod



[**Namespace List**](namespaces.md) **>** [**diagnosticinterface\_mod**](namespacediagnosticinterface__mod.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  integer, parameter, public | [**diag\_freq\_custom**](#variable-diag_freq_custom)   = `99`<br>_Custom frequency._  |
|  integer, parameter, public | [**diag\_freq\_daily**](#variable-diag_freq_daily)   = `3`<br>_Daily output._  |
|  integer, parameter, public | [**diag\_freq\_hourly**](#variable-diag_freq_hourly)   = `2`<br>_Hourly output._  |
|  integer, parameter, public | [**diag\_freq\_never**](#variable-diag_freq_never)   = `0`<br>_Never output._  |
|  integer, parameter, public | [**diag\_freq\_timestep**](#variable-diag_freq_timestep)   = `1`<br>_Every timestep._  |
|  integer, parameter, public | [**diag\_integer\_1d**](#variable-diag_integer_1d)   = `12`<br> |
|  integer, parameter, public | [**diag\_integer\_2d**](#variable-diag_integer_2d)   = `13`<br> |
|  integer, parameter, public | [**diag\_integer\_3d**](#variable-diag_integer_3d)   = `14`<br> |
|  integer, parameter, public | [**diag\_integer\_scalar**](#variable-diag_integer_scalar)   = `11`<br> |
|  integer, parameter, public | [**diag\_logical\_1d**](#variable-diag_logical_1d)   = `22`<br> |
|  integer, parameter, public | [**diag\_logical\_2d**](#variable-diag_logical_2d)   = `23`<br> |
|  integer, parameter, public | [**diag\_logical\_3d**](#variable-diag_logical_3d)   = `24`<br> |
|  integer, parameter, public | [**diag\_logical\_scalar**](#variable-diag_logical_scalar)   = `21`<br> |
|  integer, parameter, public | [**diag\_real\_1d**](#variable-diag_real_1d)   = `2`<br> |
|  integer, parameter, public | [**diag\_real\_2d**](#variable-diag_real_2d)   = `3`<br> |
|  integer, parameter, public | [**diag\_real\_3d**](#variable-diag_real_3d)   = `4`<br> |
|  integer, parameter, public | [**diag\_real\_scalar**](#variable-diag_real_scalar)   = `1`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**diag\_data\_allocate**](#function-diag_data_allocate) (class([**diagnosticdatatype**](namespacediagnosticinterface__mod.md#none-diagnosticdatatype)), intent(inout) this, integer, intent(in) data\_type, integer, dimension(:), intent(in), optional dims, integer, intent(out) rc) <br>_Allocate diagnostic data storage._  |




























## Public Attributes Documentation




### variable diag\_freq\_custom 

_Custom frequency._ 
```Fortran
integer, parameter, public diagnosticinterface_mod::diag_freq_custom;
```




<hr>



### variable diag\_freq\_daily 

_Daily output._ 
```Fortran
integer, parameter, public diagnosticinterface_mod::diag_freq_daily;
```




<hr>



### variable diag\_freq\_hourly 

_Hourly output._ 
```Fortran
integer, parameter, public diagnosticinterface_mod::diag_freq_hourly;
```




<hr>



### variable diag\_freq\_never 

_Never output._ 
```Fortran
integer, parameter, public diagnosticinterface_mod::diag_freq_never;
```




<hr>



### variable diag\_freq\_timestep 

_Every timestep._ 
```Fortran
integer, parameter, public diagnosticinterface_mod::diag_freq_timestep;
```




<hr>



### variable diag\_integer\_1d 

```Fortran
integer, parameter, public diagnosticinterface_mod::diag_integer_1d;
```




<hr>



### variable diag\_integer\_2d 

```Fortran
integer, parameter, public diagnosticinterface_mod::diag_integer_2d;
```




<hr>



### variable diag\_integer\_3d 

```Fortran
integer, parameter, public diagnosticinterface_mod::diag_integer_3d;
```




<hr>



### variable diag\_integer\_scalar 

```Fortran
integer, parameter, public diagnosticinterface_mod::diag_integer_scalar;
```




<hr>



### variable diag\_logical\_1d 

```Fortran
integer, parameter, public diagnosticinterface_mod::diag_logical_1d;
```




<hr>



### variable diag\_logical\_2d 

```Fortran
integer, parameter, public diagnosticinterface_mod::diag_logical_2d;
```




<hr>



### variable diag\_logical\_3d 

```Fortran
integer, parameter, public diagnosticinterface_mod::diag_logical_3d;
```




<hr>



### variable diag\_logical\_scalar 

```Fortran
integer, parameter, public diagnosticinterface_mod::diag_logical_scalar;
```




<hr>



### variable diag\_real\_1d 

```Fortran
integer, parameter, public diagnosticinterface_mod::diag_real_1d;
```




<hr>



### variable diag\_real\_2d 

```Fortran
integer, parameter, public diagnosticinterface_mod::diag_real_2d;
```




<hr>



### variable diag\_real\_3d 

```Fortran
integer, parameter, public diagnosticinterface_mod::diag_real_3d;
```




<hr>



### variable diag\_real\_scalar 

```Fortran
integer, parameter, public diagnosticinterface_mod::diag_real_scalar;
```




<hr>
## Public Functions Documentation




### function diag\_data\_allocate 

_Allocate diagnostic data storage._ 
```Fortran
subroutine diagnosticinterface_mod::diag_data_allocate (
    class( diagnosticdatatype ), intent(inout) this,
    integer, intent(in) data_type,
    integer, dimension(:), intent(in), optional dims,
    integer, intent(out) rc
) 
```





**Parameters:**


* `this` DiagnosticDataType instance 
* `data_type` Type of data to allocate 
* `dims` Dimensions for multi-dimensional arrays 
* `rc` Return code 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/DiagnosticInterface_Mod.F90`

