

# Namespace unitconversion\_mod



[**Namespace List**](namespaces.md) **>** [**unitconversion\_mod**](namespaceunitconversion__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  real(fp) function, public | [**calculate\_air\_density**](#function-calculate_air_density) (real(fp), intent(in) temperature, real(fp), intent(in) pressure, real(fp), intent(in), optional humidity) <br>_Calculate air density._  |
|  real(fp) function, public | [**calculate\_molecular\_weight**](#function-calculate_molecular_weight) (character(len=\*), intent(in) formula) <br>_Calculate molecular weight from formula._  |
|  subroutine, public | [**convert\_concentration**](#function-convert_concentration) (real(fp), intent(in) input\_value, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, real(fp), intent(in) molecular\_weight, real(fp), intent(in) temperature, real(fp), intent(in) pressure, real(fp), intent(out) output\_value, integer, intent(out) rc) <br>_Convert concentration units between different systems._  |
|  real(fp) function, public | [**convert\_flux**](#function-convert_flux) (real(fp), intent(in) flux\_in, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, real(fp), intent(in) molecular\_weight, integer, intent(out) rc) <br>_Convert flux units._  |
|  real(fp) function, public | [**convert\_imperial\_area**](#function-convert_imperial_area) (real(fp), intent(in) area\_in, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, integer, intent(out) rc) <br>_Convert imperial area units to metric._  |
|  real(fp) function, public | [**convert\_imperial\_energy**](#function-convert_imperial_energy) (real(fp), intent(in) energy\_in, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, integer, intent(out) rc) <br>_Convert imperial energy units to metric._  |
|  real(fp) function, public | [**convert\_imperial\_force**](#function-convert_imperial_force) (real(fp), intent(in) force\_in, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, integer, intent(out) rc) <br>_Convert imperial force units to metric._  |
|  real(fp) function, public | [**convert\_imperial\_length**](#function-convert_imperial_length) (real(fp), intent(in) length\_in, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, integer, intent(out) rc) <br>_Convert imperial length units to metric._  |
|  real(fp) function, public | [**convert\_imperial\_mass**](#function-convert_imperial_mass) (real(fp), intent(in) mass\_in, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, integer, intent(out) rc) <br>_Convert imperial mass units to metric._  |
|  real(fp) function, public | [**convert\_imperial\_pressure**](#function-convert_imperial_pressure) (real(fp), intent(in) pressure\_in, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, integer, intent(out) rc) <br>_Convert imperial pressure units to metric._  |
|  real(fp) function, public | [**convert\_imperial\_speed**](#function-convert_imperial_speed) (real(fp), intent(in) speed\_in, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, integer, intent(out) rc) <br>_Convert imperial speed units to metric._  |
|  real(fp) function, public | [**convert\_imperial\_temperature**](#function-convert_imperial_temperature) (real(fp), intent(in) temp\_in, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, integer, intent(out) rc) <br>_Convert imperial temperature units to metric._  |
|  real(fp) function, public | [**convert\_imperial\_volume**](#function-convert_imperial_volume) (real(fp), intent(in) volume\_in, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, integer, intent(out) rc) <br>_Convert imperial volume units to metric._  |
|  real(fp) function, public | [**convert\_mass\_units**](#function-convert_mass_units) (real(fp), intent(in) mass\_in, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, integer, intent(out) rc) <br>_Convert mass units._  |
|  real(fp) function, public | [**convert\_pressure**](#function-convert_pressure) (real(fp), intent(in) pressure\_in, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, integer, intent(out) rc) <br>_Convert pressure units._  |
|  real(fp) function, public | [**convert\_rate\_constant**](#function-convert_rate_constant) (real(fp), intent(in) rate\_in, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, integer, intent(out) rc) <br>_Convert rate constant units._  |
|  real(fp) function, public | [**convert\_temperature**](#function-convert_temperature) (real(fp), intent(in) temp\_in, character(len=\*), intent(in) input\_units, character(len=\*), intent(in) output\_units, integer, intent(out) rc) <br>_Convert temperature units._  |




























## Public Functions Documentation




### function calculate\_air\_density 

_Calculate air density._ 
```Fortran
real(fp) function, public unitconversion_mod::calculate_air_density (
    real(fp), intent(in) temperature,
    real(fp), intent(in) pressure,
    real(fp), intent(in), optional humidity
) 
```





**Parameters:**


* `temperature` [K] 
* `pressure` [Pa] 
* `humidity` relative humidity [0-1] 



**Returns:**

[kg/m³] 





        

<hr>



### function calculate\_molecular\_weight 

_Calculate molecular weight from formula._ 
```Fortran
real(fp) function, public unitconversion_mod::calculate_molecular_weight (
    character(len=*), intent(in) formula
) 
```




<hr>



### function convert\_concentration 

_Convert concentration units between different systems._ 
```Fortran
subroutine, public unitconversion_mod::convert_concentration (
    real(fp), intent(in) input_value,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    real(fp), intent(in) molecular_weight,
    real(fp), intent(in) temperature,
    real(fp), intent(in) pressure,
    real(fp), intent(out) output_value,
    integer, intent(out) rc
) 
```





**Parameters:**


* `molecular_weight` [g/mol] 
* `temperature` [K] 
* `pressure` [Pa] 




        

<hr>



### function convert\_flux 

_Convert flux units._ 
```Fortran
real(fp) function, public unitconversion_mod::convert_flux (
    real(fp), intent(in) flux_in,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    real(fp), intent(in) molecular_weight,
    integer, intent(out) rc
) 
```





**Parameters:**


* `molecular_weight` [g/mol] 




        

<hr>



### function convert\_imperial\_area 

_Convert imperial area units to metric._ 
```Fortran
real(fp) function, public unitconversion_mod::convert_imperial_area (
    real(fp), intent(in) area_in,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    integer, intent(out) rc
) 
```




<hr>



### function convert\_imperial\_energy 

_Convert imperial energy units to metric._ 
```Fortran
real(fp) function, public unitconversion_mod::convert_imperial_energy (
    real(fp), intent(in) energy_in,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    integer, intent(out) rc
) 
```




<hr>



### function convert\_imperial\_force 

_Convert imperial force units to metric._ 
```Fortran
real(fp) function, public unitconversion_mod::convert_imperial_force (
    real(fp), intent(in) force_in,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    integer, intent(out) rc
) 
```




<hr>



### function convert\_imperial\_length 

_Convert imperial length units to metric._ 
```Fortran
real(fp) function, public unitconversion_mod::convert_imperial_length (
    real(fp), intent(in) length_in,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    integer, intent(out) rc
) 
```




<hr>



### function convert\_imperial\_mass 

_Convert imperial mass units to metric._ 
```Fortran
real(fp) function, public unitconversion_mod::convert_imperial_mass (
    real(fp), intent(in) mass_in,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    integer, intent(out) rc
) 
```




<hr>



### function convert\_imperial\_pressure 

_Convert imperial pressure units to metric._ 
```Fortran
real(fp) function, public unitconversion_mod::convert_imperial_pressure (
    real(fp), intent(in) pressure_in,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    integer, intent(out) rc
) 
```




<hr>



### function convert\_imperial\_speed 

_Convert imperial speed units to metric._ 
```Fortran
real(fp) function, public unitconversion_mod::convert_imperial_speed (
    real(fp), intent(in) speed_in,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    integer, intent(out) rc
) 
```




<hr>



### function convert\_imperial\_temperature 

_Convert imperial temperature units to metric._ 
```Fortran
real(fp) function, public unitconversion_mod::convert_imperial_temperature (
    real(fp), intent(in) temp_in,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    integer, intent(out) rc
) 
```




<hr>



### function convert\_imperial\_volume 

_Convert imperial volume units to metric._ 
```Fortran
real(fp) function, public unitconversion_mod::convert_imperial_volume (
    real(fp), intent(in) volume_in,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    integer, intent(out) rc
) 
```




<hr>



### function convert\_mass\_units 

_Convert mass units._ 
```Fortran
real(fp) function, public unitconversion_mod::convert_mass_units (
    real(fp), intent(in) mass_in,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    integer, intent(out) rc
) 
```




<hr>



### function convert\_pressure 

_Convert pressure units._ 
```Fortran
real(fp) function, public unitconversion_mod::convert_pressure (
    real(fp), intent(in) pressure_in,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    integer, intent(out) rc
) 
```




<hr>



### function convert\_rate\_constant 

_Convert rate constant units._ 
```Fortran
real(fp) function, public unitconversion_mod::convert_rate_constant (
    real(fp), intent(in) rate_in,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    integer, intent(out) rc
) 
```




<hr>



### function convert\_temperature 

_Convert temperature units._ 
```Fortran
real(fp) function, public unitconversion_mod::convert_temperature (
    real(fp), intent(in) temp_in,
    character(len=*), intent(in) input_units,
    character(len=*), intent(in) output_units,
    integer, intent(out) rc
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/UnitConversion_Mod.F90`

