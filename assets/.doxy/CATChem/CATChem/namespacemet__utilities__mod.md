

# Namespace met\_utilities\_mod



[**Namespace List**](namespaces.md) **>** [**met\_utilities\_mod**](namespacemet__utilities__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**arrhenius\_rate**](#function-arrhenius_rate) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) a, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) ea, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t) <br>_Calculate Arrhenius rate constant._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**brunt\_vaisala\_frequency**](#function-brunt_vaisala_frequency) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t0, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) dtdz) <br>_Calculate Brunt–Väisälä frequency squared (N^2)_  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**bulk\_richardson\_number**](#function-bulk_richardson_number) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t0, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) tz, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) u, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) z) <br>_Calculate the bulk Richardson number._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**cunningham\_correction\_factor**](#function-cunningham_correction_factor) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) dp, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) lambda) <br>_Calculate Cunningham correction factor._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**dew\_point**](#function-dew_point) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) rh) <br>_Calculate dew point temperature._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**dry\_adiabatic\_lapse\_rate**](#function-dry_adiabatic_lapse_rate) () <br>_Calculate dry adiabatic lapse rate._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**friction\_velocity**](#function-friction_velocity) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) tau, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) rho) <br>_Calculate friction velocity (u\*)_  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**henrys\_law\_constant**](#function-henrys_law_constant) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) h0, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) dh, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t0) <br>_Calculate Henry's Law constant (temperature dependent)_  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**latent\_heat\_vaporization**](#function-latent_heat_vaporization) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t) <br>_Calculate latent heat of vaporization (temperature dependent)_  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**mean\_free\_path\_air**](#function-mean_free_path_air) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) p) <br>_Calculate the mean free path of air molecules._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**mixing\_ratio**](#function-mixing_ratio) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) q) <br>_Calculate mixing ratio from specific humidity._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**monin\_obukhov\_length**](#function-monin_obukhov_length) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) ustar, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t0, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) h, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) rho) <br>_Calculate the Monin-Obukhov length._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**nuclear\_decay**](#function-nuclear_decay) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) n0, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) lambda, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t) <br>_Calculate nuclear decay (first-order)_  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**photolysis\_rate\_scaling**](#function-photolysis_rate_scaling) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) j0, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) sza) <br>_Scale photolysis rate for solar zenith angle._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**potential\_temperature**](#function-potential_temperature) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) p, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) p0) <br>_Calculate potential temperature (theta)_  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**ppm\_to\_ugm3**](#function-ppm_to_ugm3) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) ppm, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) m, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) p) <br>_Convert ppm to ug/m3._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**psi\_h\_businger**](#function-psi_h_businger) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) zeta) <br>_Businger-Dyer stability correction for heat._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**psi\_m\_businger**](#function-psi_m_businger) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) zeta) <br>_Businger-Dyer stability correction for momentum._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**psychrometric\_constant**](#function-psychrometric_constant) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) p, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) lv) <br>_Calculate the psychrometric constant._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**relative\_humidity**](#function-relative_humidity) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) qv, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) p) <br>_Calculate relative humidity._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**saturation\_mixing\_ratio**](#function-saturation_mixing_ratio) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) p, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t) <br>_Calculate saturation mixing ratio._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**saturation\_vapor\_pressure**](#function-saturation_vapor_pressure) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t) <br>_Calculate saturation vapor pressure (Clausius-Clapeyron)_  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**specific\_humidity**](#function-specific_humidity) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) r) <br>_Calculate specific humidity from mixing ratio._  |
|  integer function, public | [**stability\_classification**](#function-stability_classification) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) l) <br>_Classify atmospheric stability based on Monin-Obukhov length._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**stokes\_number**](#function-stokes_number) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) rho\_p, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) d\_p, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) u, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) mu, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) l) <br>_Calculate Stokes number from base state variables._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**stokes\_settling\_velocity**](#function-stokes_settling_velocity) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) dp, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) rho\_p, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) rho\_a, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) mu, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) cc) <br>_Calculate Stokes settling velocity for a particle._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**ugm3\_to\_ppm**](#function-ugm3_to_ppm) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) ugm3, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) m, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) p) <br>_Convert ug/m3 to ppm._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**virtual\_temperature**](#function-virtual_temperature) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) t, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) qv) <br>_Calculate virtual temperature._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)) function, public | [**wind\_profile\_loglaw**](#function-wind_profile_loglaw) (real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) ustar, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) z, real([**fp**](namespaceprecision__mod.md#variable-fp)), intent(in) z0) <br>_Calculate wind speed at height z using the log-law._  |




























## Public Functions Documentation




### function arrhenius\_rate 

_Calculate Arrhenius rate constant._ 
```Fortran
real( fp ) function, public met_utilities_mod::arrhenius_rate (
    real( fp ), intent(in) a,
    real( fp ), intent(in) ea,
    real( fp ), intent(in) t
) 
```





**Parameters:**


* `A` Pre-exponential factor [units vary] 
* `Ea` Activation energy [J/mol] 
* `T` Temperature [K] 



**Returns:**

Rate constant [units of A] SeinfeldPandis2016 





        

<hr>



### function brunt\_vaisala\_frequency 

_Calculate Brunt–Väisälä frequency squared (N^2)_ 
```Fortran
real( fp ) function, public met_utilities_mod::brunt_vaisala_frequency (
    real( fp ), intent(in) t0,
    real( fp ), intent(in) dtdz
) 
```





**Parameters:**


* `T0` Reference temperature [K] 
* `dTdz` Vertical temperature gradient [K/m] 



**Returns:**

Brunt–Väisälä frequency squared [1/s^2] WallaceHobbs2006 





        

<hr>



### function bulk\_richardson\_number 

_Calculate the bulk Richardson number._ 
```Fortran
real( fp ) function, public met_utilities_mod::bulk_richardson_number (
    real( fp ), intent(in) t0,
    real( fp ), intent(in) tz,
    real( fp ), intent(in) u,
    real( fp ), intent(in) z
) 
```





**Parameters:**


* `T0` Surface temperature [K] 
* `Tz` Temperature at height z [K] 
* `u` Wind speed at height z [m/s] 
* `z` Height above ground [m] 



**Returns:**

Bulk Richardson number (dimensionless) 





        

<hr>



### function cunningham\_correction\_factor 

_Calculate Cunningham correction factor._ 
```Fortran
real( fp ) function, public met_utilities_mod::cunningham_correction_factor (
    real( fp ), intent(in) dp,
    real( fp ), intent(in) lambda
) 
```





**Parameters:**


* `dp` Particle diameter [m] 
* `lambda` Mean free path of air [m] 



**Returns:**

Cunningham correction factor (dimensionless) 





        

<hr>



### function dew\_point 

_Calculate dew point temperature._ 
```Fortran
real( fp ) function, public met_utilities_mod::dew_point (
    real( fp ), intent(in) t,
    real( fp ), intent(in) rh
) 
```





**Parameters:**


* `T` Temperature [K] 
* `rh` Relative humidity [0-1] 



**Returns:**

Dew point temperature [K] Bolton1980 





        

<hr>



### function dry\_adiabatic\_lapse\_rate 

_Calculate dry adiabatic lapse rate._ 
```Fortran
real( fp ) function, public met_utilities_mod::dry_adiabatic_lapse_rate () 
```





**Returns:**

Dry adiabatic lapse rate [K/m] 





        

<hr>



### function friction\_velocity 

_Calculate friction velocity (u\*)_ 
```Fortran
real( fp ) function, public met_utilities_mod::friction_velocity (
    real( fp ), intent(in) tau,
    real( fp ), intent(in) rho
) 
```





**Parameters:**


* `tau` Surface shear stress [N/m^2] 
* `rho` Air density [kg/m^3] 



**Returns:**

Friction velocity [m/s] 





        

<hr>



### function henrys\_law\_constant 

_Calculate Henry's Law constant (temperature dependent)_ 
```Fortran
real( fp ) function, public met_utilities_mod::henrys_law_constant (
    real( fp ), intent(in) h0,
    real( fp ), intent(in) dh,
    real( fp ), intent(in) t,
    real( fp ), intent(in) t0
) 
```





**Parameters:**


* `H0` Reference Henry's constant [mol/(m^3\*Pa)] 
* `dH` Enthalpy of solution [J/mol] 
* `T` Temperature [K] 
* `T0` Reference temperature [K] 



**Returns:**

Henry's Law constant at T [mol/(m^3\*Pa)] Sander2015 





        

<hr>



### function latent\_heat\_vaporization 

_Calculate latent heat of vaporization (temperature dependent)_ 
```Fortran
real( fp ) function, public met_utilities_mod::latent_heat_vaporization (
    real( fp ), intent(in) t
) 
```





**Parameters:**


* `T` Temperature [K] 



**Returns:**

Latent heat of vaporization [J/kg] 





        

<hr>



### function mean\_free\_path\_air 

_Calculate the mean free path of air molecules._ 
```Fortran
real( fp ) function, public met_utilities_mod::mean_free_path_air (
    real( fp ), intent(in) t,
    real( fp ), intent(in) p
) 
```





**Parameters:**


* `T` Temperature [K] 
* `p` Pressure [Pa] 



**Returns:**

Mean free path [m] SeinfeldPandis2016 





        

<hr>



### function mixing\_ratio 

_Calculate mixing ratio from specific humidity._ 
```Fortran
real( fp ) function, public met_utilities_mod::mixing_ratio (
    real( fp ), intent(in) q
) 
```





**Parameters:**


* `q` Specific humidity [kg/kg] 



**Returns:**

Mixing ratio [kg/kg] 





        

<hr>



### function monin\_obukhov\_length 

_Calculate the Monin-Obukhov length._ 
```Fortran
real( fp ) function, public met_utilities_mod::monin_obukhov_length (
    real( fp ), intent(in) ustar,
    real( fp ), intent(in) t0,
    real( fp ), intent(in) h,
    real( fp ), intent(in) rho
) 
```





**Parameters:**


* `ustar` Friction velocity [m/s] 
* `T0` Surface temperature [K] 
* `H` Sensible heat flux [W/m^2] 
* `rho` Air density [kg/m^3] 



**Returns:**

Monin-Obukhov length [m] 





        

<hr>



### function nuclear\_decay 

_Calculate nuclear decay (first-order)_ 
```Fortran
real( fp ) function, public met_utilities_mod::nuclear_decay (
    real( fp ), intent(in) n0,
    real( fp ), intent(in) lambda,
    real( fp ), intent(in) t
) 
```





**Parameters:**


* `N0` Initial quantity 
* `lambda` Decay constant [1/s] 
* `t` Time [s] 



**Returns:**

Remaining quantity after time t 





        

<hr>



### function photolysis\_rate\_scaling 

_Scale photolysis rate for solar zenith angle._ 
```Fortran
real( fp ) function, public met_utilities_mod::photolysis_rate_scaling (
    real( fp ), intent(in) j0,
    real( fp ), intent(in) sza
) 
```





**Parameters:**


* `J0` Base photolysis rate [1/s] 
* `sza` Solar zenith angle [degrees] 



**Returns:**

Scaled photolysis rate [1/s] 





        

<hr>



### function potential\_temperature 

_Calculate potential temperature (theta)_ 
```Fortran
real( fp ) function, public met_utilities_mod::potential_temperature (
    real( fp ), intent(in) t,
    real( fp ), intent(in) p,
    real( fp ), intent(in) p0
) 
```





**Parameters:**


* `T` Temperature [K] 
* `p` Pressure [Pa] 
* `p0` Surface pressure [Pa] 



**Returns:**

Potential temperature [K] WallaceHobbs2006 





        

<hr>



### function ppm\_to\_ugm3 

_Convert ppm to ug/m3._ 
```Fortran
real( fp ) function, public met_utilities_mod::ppm_to_ugm3 (
    real( fp ), intent(in) ppm,
    real( fp ), intent(in) m,
    real( fp ), intent(in) t,
    real( fp ), intent(in) p
) 
```





**Parameters:**


* `ppm` Concentration [ppm] 
* `M` Molar mass [g/mol] 
* `T` Temperature [K] 
* `p` Pressure [Pa] 



**Returns:**

Concentration [ug/m3] 





        

<hr>



### function psi\_h\_businger 

_Businger-Dyer stability correction for heat._ 
```Fortran
real( fp ) function, public met_utilities_mod::psi_h_businger (
    real( fp ), intent(in) zeta
) 
```





**Parameters:**


* `zeta` z/L (dimensionless stability parameter) 



**Returns:**

Psi\_h (stability correction for heat) 





        

<hr>



### function psi\_m\_businger 

_Businger-Dyer stability correction for momentum._ 
```Fortran
real( fp ) function, public met_utilities_mod::psi_m_businger (
    real( fp ), intent(in) zeta
) 
```





**Parameters:**


* `zeta` z/L (dimensionless stability parameter) 



**Returns:**

Psi\_m (stability correction for momentum) 





        

<hr>



### function psychrometric\_constant 

_Calculate the psychrometric constant._ 
```Fortran
real( fp ) function, public met_utilities_mod::psychrometric_constant (
    real( fp ), intent(in) p,
    real( fp ), intent(in) lv
) 
```





**Parameters:**


* `p` Pressure [Pa] 
* `Lv` Latent heat of vaporization [J/kg] 



**Returns:**

Psychrometric constant [Pa/K] 





        

<hr>



### function relative\_humidity 

_Calculate relative humidity._ 
```Fortran
real( fp ) function, public met_utilities_mod::relative_humidity (
    real( fp ), intent(in) t,
    real( fp ), intent(in) qv,
    real( fp ), intent(in) p
) 
```





**Parameters:**


* `T` Temperature [K] 
* `qv` Water vapor mixing ratio [kg/kg] 
* `p` Pressure [Pa] 



**Returns:**

Relative humidity [0-1] WallaceHobbs2006 





        

<hr>



### function saturation\_mixing\_ratio 

_Calculate saturation mixing ratio._ 
```Fortran
real( fp ) function, public met_utilities_mod::saturation_mixing_ratio (
    real( fp ), intent(in) p,
    real( fp ), intent(in) t
) 
```





**Parameters:**


* `p` Pressure [Pa] 
* `T` Temperature [K] 



**Returns:**

Saturation mixing ratio [kg/kg] 





        

<hr>



### function saturation\_vapor\_pressure 

_Calculate saturation vapor pressure (Clausius-Clapeyron)_ 
```Fortran
real( fp ) function, public met_utilities_mod::saturation_vapor_pressure (
    real( fp ), intent(in) t
) 
```





**Parameters:**


* `T` Temperature [K] 



**Returns:**

Saturation vapor pressure [Pa] Bolton1980 





        

<hr>



### function specific\_humidity 

_Calculate specific humidity from mixing ratio._ 
```Fortran
real( fp ) function, public met_utilities_mod::specific_humidity (
    real( fp ), intent(in) r
) 
```





**Parameters:**


* `r` Mixing ratio [kg/kg] 



**Returns:**

Specific humidity [kg/kg] 





        

<hr>



### function stability\_classification 

_Classify atmospheric stability based on Monin-Obukhov length._ 
```Fortran
integer function, public met_utilities_mod::stability_classification (
    real( fp ), intent(in) l
) 
```





**Parameters:**


* `L` Monin-Obukhov length [m] 



**Returns:**

Stability class: -1 (unstable), 0 (neutral), 1 (stable) 





        

<hr>



### function stokes\_number 

_Calculate Stokes number from base state variables._ 
```Fortran
real( fp ) function, public met_utilities_mod::stokes_number (
    real( fp ), intent(in) rho_p,
    real( fp ), intent(in) d_p,
    real( fp ), intent(in) u,
    real( fp ), intent(in) mu,
    real( fp ), intent(in) l
) 
```





**Parameters:**


* `rho_p` Particle density [kg/m^3] 
* `d_p` Particle diameter [m] 
* `U` Characteristic velocity [m/s] 
* `mu` Dynamic viscosity [kg/m/s] 
* `L` Characteristic length scale [m] 



**Returns:**

Stokes number (dimensionless) 





        

<hr>



### function stokes\_settling\_velocity 

_Calculate Stokes settling velocity for a particle._ 
```Fortran
real( fp ) function, public met_utilities_mod::stokes_settling_velocity (
    real( fp ), intent(in) dp,
    real( fp ), intent(in) rho_p,
    real( fp ), intent(in) rho_a,
    real( fp ), intent(in) mu,
    real( fp ), intent(in) cc
) 
```





**Parameters:**


* `dp` Particle diameter [m] 
* `rho_p` Particle density [kg/m3] 
* `rho_a` Air density [kg/m3] 
* `mu` Air dynamic viscosity [kg/m/s] 
* `Cc` Cunningham correction factor 



**Returns:**

Settling velocity [m/s] 





        

<hr>



### function ugm3\_to\_ppm 

_Convert ug/m3 to ppm._ 
```Fortran
real( fp ) function, public met_utilities_mod::ugm3_to_ppm (
    real( fp ), intent(in) ugm3,
    real( fp ), intent(in) m,
    real( fp ), intent(in) t,
    real( fp ), intent(in) p
) 
```





**Parameters:**


* `ugm3` Concentration [ug/m3] 
* `M` Molar mass [g/mol] 
* `T` Temperature [K] 
* `p` Pressure [Pa] 



**Returns:**

Concentration [ppm] 





        

<hr>



### function virtual\_temperature 

_Calculate virtual temperature._ 
```Fortran
real( fp ) function, public met_utilities_mod::virtual_temperature (
    real( fp ), intent(in) t,
    real( fp ), intent(in) qv
) 
```





**Parameters:**


* `T` Temperature [K] 
* `qv` Water vapor mixing ratio [kg/kg] 



**Returns:**

Virtual temperature [K] WallaceHobbs2006 





        

<hr>



### function wind\_profile\_loglaw 

_Calculate wind speed at height z using the log-law._ 
```Fortran
real( fp ) function, public met_utilities_mod::wind_profile_loglaw (
    real( fp ), intent(in) ustar,
    real( fp ), intent(in) z,
    real( fp ), intent(in) z0
) 
```





**Parameters:**


* `ustar` Friction velocity [m/s] 
* `z` Height above ground [m] 
* `z0` Surface roughness length [m] 



**Returns:**

Wind speed at height z [m/s] 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/met_utilities_mod.F90`

