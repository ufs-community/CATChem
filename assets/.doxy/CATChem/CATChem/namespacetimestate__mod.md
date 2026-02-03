

# Namespace timestate\_mod



[**Namespace List**](namespaces.md) **>** [**timestate\_mod**](namespacetimestate__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  real function | [**get\_cos\_sza**](#function-get_cos_sza) (class([**timestatetype**](namespacetimestate__mod.md#none-timestatetype)), intent(in) this, real, intent(in) lat, real, intent(in) lon, logical, intent(in), optional mid\_timestep) <br>_Compute solar zenith angle (degrees) using latitude, longitude, and time of day._  |
|  logical function, public | [**is\_global\_holiday**](#function-is_global_holiday) (integer, intent(in) month, integer, intent(in) day) <br>_Check if a date is a global holiday._  |
|  logical function, public | [**is\_us\_holiday**](#function-is_us_holiday) (integer, intent(in) month, integer, intent(in) day) <br>_Check if a date is a U.S._  |




























## Public Functions Documentation




### function get\_cos\_sza 

_Compute solar zenith angle (degrees) using latitude, longitude, and time of day._ 
```Fortran
real function timestate_mod::get_cos_sza (
    class( timestatetype ), intent(in) this,
    real, intent(in) lat,
    real, intent(in) lon,
    logical, intent(in), optional mid_timestep
) 
```




<hr>



### function is\_global\_holiday 

_Check if a date is a global holiday._ 
```Fortran
logical function, public timestate_mod::is_global_holiday (
    integer, intent(in) month,
    integer, intent(in) day
) 
```




<hr>



### function is\_us\_holiday 

_Check if a date is a U.S._ 
```Fortran
logical function, public timestate_mod::is_us_holiday (
    integer, intent(in) month,
    integer, intent(in) day
) 
```



holiday 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/TimeState_Mod.F90`

