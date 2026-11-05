

# Namespace utilities\_mod



[**Namespace List**](namespaces.md) **>** [**utilities\_mod**](namespaceutilities__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**calculate\_air\_density**](#function-calculate_air_density) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) pressure, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) temperature, integer, intent(out) rc) <br>_Calculate air density using ideal gas law._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**calculate\_geopotential\_height**](#function-calculate_geopotential_height) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) p1, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) p2, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) tv\_mean, integer, intent(out) rc) <br>_Calculate geopotential height difference between two pressure levels._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**calculate\_scale\_height**](#function-calculate_scale_height) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) temperature, integer, intent(out) rc, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in), optional z) <br>_Calculate atmospheric scale height using geopotential height._  |
|  subroutine, public | [**check\_array\_bounds**](#function-check_array_bounds) (integer, intent(in) index, integer, intent(in) array\_size, character(len=\*), intent(in) array\_name, integer, intent(out) rc) <br>_Check array bounds safely._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**convert\_pressure\_units**](#function-convert_pressure_units) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) pressure\_in, character(len=\*), intent(in) unit\_in, character(len=\*), intent(in) unit\_out, integer, intent(out) rc) <br>_Convert pressure between different units._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**convert\_temperature\_units**](#function-convert_temperature_units) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) temp\_in, character(len=\*), intent(in) unit\_in, character(len=\*), intent(in) unit\_out, integer, intent(out) rc) <br>_Convert temperature between different units._  |
|  logical function, public | [**is\_valid\_pressure**](#function-is_valid_pressure) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) pressure) <br>_Check if pressure is in valid range._  |
|  logical function, public | [**is\_valid\_temperature**](#function-is_valid_temperature) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) temperature) <br>_Check if temperature is in valid range._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**safe\_divide**](#function-safe_divide) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) numerator, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) denominator, integer, intent(out) rc) <br>_Safe division with zero check._  |
|  subroutine, public | [**validate\_atmospheric\_constants**](#function-validate_atmospheric_constants) (integer, intent(out) rc, type([**errormanagertype**](namespaceerror__mod.md#none-errormanagertype)), intent(inout), optional error\_mgr) <br>_Validate atmospheric constants for physical consistency._  |




























## Public Functions Documentation




### function calculate\_air\_density 

_Calculate air density using ideal gas law._ 
```Fortran
real( fp ) function, public utilities_mod::calculate_air_density (
    real( fp ), intent(in) pressure,
    real( fp ), intent(in) temperature,
    integer, intent(out) rc
) 
```





**Parameters:**


* `pressure` Pressure [Pa] 
* `temperature` Temperature [K] 
* `rc` Return code 




        

<hr>



### function calculate\_geopotential\_height 

_Calculate geopotential height difference between two pressure levels._ 
```Fortran
real( fp ) function, public utilities_mod::calculate_geopotential_height (
    real( fp ), intent(in) p1,
    real( fp ), intent(in) p2,
    real( fp ), intent(in) tv_mean,
    integer, intent(out) rc
) 
```





**Parameters:**


* `p1` Pressure at lower level [Pa] 
* `p2` Pressure at upper level [Pa] 
* `Tv_mean` Mean virtual temperature between levels [K] 
* `rc` Return code 




        

<hr>



### function calculate\_scale\_height 

_Calculate atmospheric scale height using geopotential height._ 
```Fortran
real( fp ) function, public utilities_mod::calculate_scale_height (
    real( fp ), intent(in) temperature,
    integer, intent(out) rc,
    real( fp ), intent(in), optional z
) 
```





**Parameters:**


* `temperature` Temperature [K] 
* `rc` Return code 
* `z` (optional) Geometric height above surface [m] 




        

<hr>



### function check\_array\_bounds 

_Check array bounds safely._ 
```Fortran
subroutine, public utilities_mod::check_array_bounds (
    integer, intent(in) index,
    integer, intent(in) array_size,
    character(len=*), intent(in) array_name,
    integer, intent(out) rc
) 
```





**Parameters:**


* `index` Index to check 
* `array_size` Size of array 
* `array_name` Name of array for error reporting 
* `rc` Return code 




        

<hr>



### function convert\_pressure\_units 

_Convert pressure between different units._ 
```Fortran
real( fp ) function, public utilities_mod::convert_pressure_units (
    real( fp ), intent(in) pressure_in,
    character(len=*), intent(in) unit_in,
    character(len=*), intent(in) unit_out,
    integer, intent(out) rc
) 
```



This function converts pressure values between common atmospheric units.




**Parameters:**


* `pressure_in` Input pressure value 
* `unit_in` Input unit ('Pa', 'hPa', 'mbar', 'atm', 'mmHg', 'torr') 
* `unit_out` Output unit 
* `rc` Return code 




        

<hr>



### function convert\_temperature\_units 

_Convert temperature between different units._ 
```Fortran
real( fp ) function, public utilities_mod::convert_temperature_units (
    real( fp ), intent(in) temp_in,
    character(len=*), intent(in) unit_in,
    character(len=*), intent(in) unit_out,
    integer, intent(out) rc
) 
```



This function converts temperature values between Kelvin, Celsius, and Fahrenheit.




**Parameters:**


* `temp_in` Input temperature value 
* `unit_in` Input unit ('K', 'C', 'F') 
* `unit_out` Output unit 
* `rc` Return code 




        

<hr>



### function is\_valid\_pressure 

_Check if pressure is in valid range._ 
```Fortran
logical function, public utilities_mod::is_valid_pressure (
    real( fp ), intent(in) pressure
) 
```





**Parameters:**


* `pressure` Pressure [Pa] 




        

<hr>



### function is\_valid\_temperature 

_Check if temperature is in valid range._ 
```Fortran
logical function, public utilities_mod::is_valid_temperature (
    real( fp ), intent(in) temperature
) 
```





**Parameters:**


* `temperature` Temperature [K] 




        

<hr>



### function safe\_divide 

_Safe division with zero check._ 
```Fortran
real( fp ) function, public utilities_mod::safe_divide (
    real( fp ), intent(in) numerator,
    real( fp ), intent(in) denominator,
    integer, intent(out) rc
) 
```





**Parameters:**


* `numerator` Numerator 
* `denominator` Denominator 
* `rc` Return code 




        

<hr>



### function validate\_atmospheric\_constants 

_Validate atmospheric constants for physical consistency._ 
```Fortran
subroutine, public utilities_mod::validate_atmospheric_constants (
    integer, intent(out) rc,
    type( errormanagertype ), intent(inout), optional error_mgr
) 
```



This subroutine performs runtime validation of atmospheric constants to ensure physical consistency and catch any compilation issues.




**Parameters:**


* `rc` Return code (CC\_SUCCESS if all constants are valid) 
* `error_mgr` Optional error manager for enhanced reporting 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/utilities_mod.F90`

