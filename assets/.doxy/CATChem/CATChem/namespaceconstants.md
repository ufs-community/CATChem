

# Namespace constants



[**Namespace List**](namespaces.md) **>** [**constants**](namespaceconstants.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**airmw**](#variable-airmw)   = `28.9644\_fp`<br>_Average molecular weight of dry air [g/mol]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**atm**](#variable-atm)   = `1.01325[**e**](namespaceconstants.md#variable-e)+5\_fp`<br>_Standard atmospheric pressure [Pa]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**avo**](#variable-avo)   = `6.022140857[**e**](namespaceconstants.md#variable-e)+23\_fp`<br>_Avogadro's number [particles/mol]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**boltz**](#variable-boltz)   = `1.38064852[**e**](namespaceconstants.md#variable-e)-23\_fp`<br>_Boltzmann's constant [J/K]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**cconst**](#variable-cconst)   = `2.99792458[**e**](namespaceconstants.md#variable-e)+8\_fp`<br>_Speed of light in vacuum [m/s]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**consvap**](#variable-consvap)   = `6.1078[**e**](namespaceconstants.md#variable-e)+03\_fp / ( BOLTZ \* 1[**e**](namespaceconstants.md#variable-e)+7\_fp )`<br>_Condensation vapor pressure factor._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**cp**](#variable-cp)   = `1.0046[**e**](namespaceconstants.md#variable-e)+3\_fp`<br>_Specific heat of dry air at constant pressure [J/kg/K]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**cv**](#variable-cv)   = `7.1760e+2\_fp`<br>_Specific heat of dry air at constant volume [J/kg/K]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**e**](#variable-e)   = `2.718281828459045235360287471352\_fp`<br>_Euler's number (dimensionless)_  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**g0**](#variable-g0)   = `9.80665[**e**](namespaceconstants.md#variable-e)+0\_fp`<br>_Standard gravity acceleration [m/s^2]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**g0\_100**](#variable-g0_100)   = `100.0\_fp / [**g0**](namespaceconstants.md#variable-g0)`<br>_100 divided by standard gravity_  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**h2omw**](#variable-h2omw)   = `18.016\_fp`<br>_Molecular weight of water [g/mol]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**mwcarb**](#variable-mwcarb)   = `12.01[**e**](namespaceconstants.md#variable-e)-3\_fp`<br>_Molecular weight of carbon [kg/mol]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**pi**](#variable-pi)   = `3.14159265358979323\_fp`<br>_Pi (dimensionless)_  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**pi\_180**](#variable-pi_180)   = `PI / 180.0\_fp`<br>_Radians per degree conversion factor._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**planck**](#variable-planck)   = `6.62606957[**e**](namespaceconstants.md#variable-e)-34\_fp`<br>_Planck's constant [J⋅s]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**rd**](#variable-rd)   = `287.0\_fp`<br>_Gas constant for dry air [J/K/kg]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**rdg0**](#variable-rdg0)   = `Rd / [**g0**](namespaceconstants.md#variable-g0)`<br>_Gas constant for dry air divided by gravity._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**re**](#variable-re)   = `6.3710072[**e**](namespaceconstants.md#variable-e)+6\_fp`<br>_Earth's radius [m]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**rgaslatm**](#variable-rgaslatm)   = `8.2057[**e**](namespaceconstants.md#variable-e)-2\_fp`<br>_Gas constant in L⋅atm/(K⋅mol)_  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**rstarg**](#variable-rstarg)   = `8.3144598\_fp`<br>_Universal gas constant [J/K/mol]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**rv**](#variable-rv)   = `461.00\_fp`<br>_Gas constant for water vapor [J/K/kg]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**scale\_height**](#variable-scale_height)   = `7600.0\_fp`<br>_Atmospheric scale height [m]._  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**von\_karman**](#variable-von_karman)   = `0.41\_fp`<br>_Von Karman's constant (dimensionless)_  |
|  real([**fp**](namespaceprecision__mod.md#variable-fp)), parameter, public | [**xnumolair**](#variable-xnumolair)   = `AVO / ( AIRMW \* 1.[**e**](namespaceconstants.md#variable-e)-3\_fp )`<br>_Molecules of dry air per kg dry air._  |












































## Public Attributes Documentation




### variable airmw 

_Average molecular weight of dry air [g/mol]._ 
```Fortran
real(fp), parameter, public constants::airmw;
```




<hr>



### variable atm 

_Standard atmospheric pressure [Pa]._ 
```Fortran
real(fp), parameter, public constants::atm;
```




<hr>



### variable avo 

_Avogadro's number [particles/mol]._ 
```Fortran
real(fp), parameter, public constants::avo;
```




<hr>



### variable boltz 

_Boltzmann's constant [J/K]._ 
```Fortran
real(fp), parameter, public constants::boltz;
```




<hr>



### variable cconst 

_Speed of light in vacuum [m/s]._ 
```Fortran
real(fp), parameter, public constants::cconst;
```




<hr>



### variable consvap 

_Condensation vapor pressure factor._ 
```Fortran
real(fp), parameter, public constants::consvap;
```




<hr>



### variable cp 

_Specific heat of dry air at constant pressure [J/kg/K]._ 
```Fortran
real(fp), parameter, public constants::cp;
```




<hr>



### variable cv 

_Specific heat of dry air at constant volume [J/kg/K]._ 
```Fortran
real(fp), parameter, public constants::cv;
```




<hr>



### variable e 

_Euler's number (dimensionless)_ 
```Fortran
real(fp), parameter, public constants::e;
```




<hr>



### variable g0 

_Standard gravity acceleration [m/s^2]._ 
```Fortran
real(fp), parameter, public constants::g0;
```




<hr>



### variable g0\_100 

_100 divided by standard gravity_ 
```Fortran
real(fp), parameter, public constants::g0_100;
```




<hr>



### variable h2omw 

_Molecular weight of water [g/mol]._ 
```Fortran
real(fp), parameter, public constants::h2omw;
```




<hr>



### variable mwcarb 

_Molecular weight of carbon [kg/mol]._ 
```Fortran
real(fp), parameter, public constants::mwcarb;
```




<hr>



### variable pi 

_Pi (dimensionless)_ 
```Fortran
real(fp), parameter, public constants::pi;
```




<hr>



### variable pi\_180 

_Radians per degree conversion factor._ 
```Fortran
real(fp), parameter, public constants::pi_180;
```




<hr>



### variable planck 

_Planck's constant [J⋅s]._ 
```Fortran
real(fp), parameter, public constants::planck;
```




<hr>



### variable rd 

_Gas constant for dry air [J/K/kg]._ 
```Fortran
real(fp), parameter, public constants::rd;
```




<hr>



### variable rdg0 

_Gas constant for dry air divided by gravity._ 
```Fortran
real(fp), parameter, public constants::rdg0;
```




<hr>



### variable re 

_Earth's radius [m]._ 
```Fortran
real(fp), parameter, public constants::re;
```




<hr>



### variable rgaslatm 

_Gas constant in L⋅atm/(K⋅mol)_ 
```Fortran
real(fp), parameter, public constants::rgaslatm;
```




<hr>



### variable rstarg 

_Universal gas constant [J/K/mol]._ 
```Fortran
real(fp), parameter, public constants::rstarg;
```




<hr>



### variable rv 

_Gas constant for water vapor [J/K/kg]._ 
```Fortran
real(fp), parameter, public constants::rv;
```




<hr>



### variable scale\_height 

_Atmospheric scale height [m]._ 
```Fortran
real(fp), parameter, public constants::scale_height;
```




<hr>



### variable von\_karman 

_Von Karman's constant (dimensionless)_ 
```Fortran
real(fp), parameter, public constants::von_karman;
```




<hr>



### variable xnumolair 

_Molecules of dry air per kg dry air._ 
```Fortran
real(fp), parameter, public constants::xnumolair;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/constants.F90`

