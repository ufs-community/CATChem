

# Namespace catchem\_bridge\_precision



[**Namespace List**](namespaces.md) **>** [**catchem\_bridge\_precision**](namespacecatchem__bridge__precision.md)




















## Classes

| Type | Name |
| ---: | :--- |
| interface | [**exact\_equal**](interfacecatchem__bridge__precision_1_1exact__equal.md) <br>_Exact IEEE equality, spelled without == so gfortran's -Wcompare-reals stays quiet at deliberate exact comparisons (sentinels, literal zeros)._  |
| interface | [**is\_exact\_zero**](interfacecatchem__bridge__precision_1_1is__exact__zero.md) <br>_Exact test for zero: exact\_equal(x, 0)._  |
| interface | [**rae**](interfacecatchem__bridge__precision_1_1rae.md) <br> |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  integer, parameter, public | [**default\_real\_kind**](#variable-default_real_kind)   = `KIND(0.0)`<br> |
|  integer, parameter, public | [**f4**](#variable-f4)   = `KIND( 0.0\_4 )`<br>_KIND parameter for 4-byte precision._  |
|  integer, parameter, public | [**f8**](#variable-f8)   = `KIND( 0.0\_8 )`<br>_KIND parameter for 8-byte precision._  |
|  integer, parameter, public | [**fp**](#variable-fp)   = `[**f8**](namespacecatchem__bridge__precision.md#variable-f8)`<br>_8-byte (default)_  |
|  real([**fp**](namespacecatchem__bridge__precision.md#variable-fp)), parameter, public | [**missing**](#variable-missing)   = `-999.0\_fp`<br>_Missing real value (kind=fp)_  |
|  logical, parameter, public | [**missing\_bool**](#variable-missing_bool)   = `.FALSE.`<br>_Missing boolean value._  |
|  real([**f8**](namespacecatchem__bridge__precision.md#variable-f8)), parameter, public | [**missing\_dble**](#variable-missing_dble)   = `-999.0\_f8`<br>_Missing real value (kind=f8)_  |
|  integer, parameter, public | [**missing\_int**](#variable-missing_int)   = `-999`<br>_Missing integer value._  |
|  real([**f4**](namespacecatchem__bridge__precision.md#variable-f4)), parameter, public | [**missing\_real**](#variable-missing_real)   = `-999.0\_f4`<br>_Missing real value (kind=f4)_  |
|  character(len=7), parameter, public | [**missing\_str**](#variable-missing_str)   = `"UNKNOWN"`<br>_Missing string._  |
|  real([**fp**](namespacecatchem__bridge__precision.md#variable-fp)), parameter, public | [**one**](#variable-one)   = `1.0\_fp`<br>_One value (kind=fp)_  |
|  real([**f8**](namespacecatchem__bridge__precision.md#variable-f8)), parameter, public | [**one\_dble**](#variable-one_dble)   = `1.0\_f8`<br>_One value (kind=f8)_  |
|  real([**f4**](namespacecatchem__bridge__precision.md#variable-f4)), parameter, public | [**one\_real**](#variable-one_real)   = `1.0\_f4`<br>_One value (kind=f4)_  |
|  real([**fp**](namespacecatchem__bridge__precision.md#variable-fp)), parameter, public | [**tiny\_**](#variable-tiny_)   = `TINY\_REAL`<br> |
|  real([**f8**](namespacecatchem__bridge__precision.md#variable-f8)), parameter, public | [**tiny\_dble**](#variable-tiny_dble)   = `1.0e-31\_f8`<br>_A small value (kind=f8)_  |
|  real([**f4**](namespacecatchem__bridge__precision.md#variable-f4)), parameter, public | [**tiny\_real**](#variable-tiny_real)   = `1.0e-16\_f4`<br>_A small value (kind=f4)_  |
|  real([**fp**](namespacecatchem__bridge__precision.md#variable-fp)), parameter, public | [**zero**](#variable-zero)   = `0.0\_fp`<br>_Zero value (kind=fp)_  |
|  real([**f8**](namespacecatchem__bridge__precision.md#variable-f8)), parameter, public | [**zero\_dble**](#variable-zero_dble)   = `0.0\_f8`<br>_Zero value (kind=f8)_  |
|  real([**f4**](namespacecatchem__bridge__precision.md#variable-f4)), parameter, public | [**zero\_real**](#variable-zero_real)   = `0.0\_f4`<br>_Zero value (kind=f4)_  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  elemental logical function | [**exact\_equal\_f4**](#function-exact_equal_f4) (real([**f4**](namespacecatchem__bridge__precision.md#variable-f4)), intent(in) a, real([**f4**](namespacecatchem__bridge__precision.md#variable-f4)), intent(in) b) <br>_Exact equality (see the_ `exact_equal` _interface)_ |
|  elemental logical function | [**exact\_equal\_f8**](#function-exact_equal_f8) (real([**f8**](namespacecatchem__bridge__precision.md#variable-f8)), intent(in) a, real([**f8**](namespacecatchem__bridge__precision.md#variable-f8)), intent(in) b) <br>_Exact equality (see the_ `exact_equal` _interface)_ |
|  elemental logical function | [**is\_exact\_zero\_f4**](#function-is_exact_zero_f4) (real([**f4**](namespacecatchem__bridge__precision.md#variable-f4)), intent(in) x) <br>_Exact zero test (see the_ `is_exact_zero` _interface)_ |
|  elemental logical function | [**is\_exact\_zero\_f8**](#function-is_exact_zero_f8) (real([**f8**](namespacecatchem__bridge__precision.md#variable-f8)), intent(in) x) <br>_Exact zero test (see the_ `is_exact_zero` _interface)_ |
|  logical function | [**rae\_f4**](#function-rae_f4) (real([**f4**](namespacecatchem__bridge__precision.md#variable-f4)), intent(in) a, real([**f4**](namespacecatchem__bridge__precision.md#variable-f4)), intent(in) b) <br>_Real approximately equal:_ `abs(a - b) < tiny(a)` __ |
|  logical function | [**rae\_f8**](#function-rae_f8) (real([**f8**](namespacecatchem__bridge__precision.md#variable-f8)), intent(in) a, real([**f8**](namespacecatchem__bridge__precision.md#variable-f8)), intent(in) b) <br>_Real approximately equal:_ `abs(a - b) < tiny(a)` __ |




























## Public Attributes Documentation




### variable default\_real\_kind 

```Fortran
integer, parameter, public catchem_bridge_precision::default_real_kind;
```




<hr>



### variable f4 

_KIND parameter for 4-byte precision._ 
```Fortran
integer, parameter, public catchem_bridge_precision::f4;
```




<hr>



### variable f8 

_KIND parameter for 8-byte precision._ 
```Fortran
integer, parameter, public catchem_bridge_precision::f8;
```




<hr>



### variable fp 

_8-byte (default)_ 
```Fortran
integer, parameter, public catchem_bridge_precision::fp;
```




<hr>



### variable missing 

_Missing real value (kind=fp)_ 
```Fortran
real(fp), parameter, public catchem_bridge_precision::missing;
```




<hr>



### variable missing\_bool 

_Missing boolean value._ 
```Fortran
logical, parameter, public catchem_bridge_precision::missing_bool;
```




<hr>



### variable missing\_dble 

_Missing real value (kind=f8)_ 
```Fortran
real(f8), parameter, public catchem_bridge_precision::missing_dble;
```




<hr>



### variable missing\_int 

_Missing integer value._ 
```Fortran
integer, parameter, public catchem_bridge_precision::missing_int;
```




<hr>



### variable missing\_real 

_Missing real value (kind=f4)_ 
```Fortran
real(f4), parameter, public catchem_bridge_precision::missing_real;
```




<hr>



### variable missing\_str 

_Missing string._ 
```Fortran
character(len=7), parameter, public catchem_bridge_precision::missing_str;
```




<hr>



### variable one 

_One value (kind=fp)_ 
```Fortran
real(fp), parameter, public catchem_bridge_precision::one;
```




<hr>



### variable one\_dble 

_One value (kind=f8)_ 
```Fortran
real(f8), parameter, public catchem_bridge_precision::one_dble;
```




<hr>



### variable one\_real 

_One value (kind=f4)_ 
```Fortran
real(f4), parameter, public catchem_bridge_precision::one_real;
```




<hr>



### variable tiny\_ 

```Fortran
real(fp), parameter, public catchem_bridge_precision::tiny_;
```




<hr>



### variable tiny\_dble 

_A small value (kind=f8)_ 
```Fortran
real(f8), parameter, public catchem_bridge_precision::tiny_dble;
```




<hr>



### variable tiny\_real 

_A small value (kind=f4)_ 
```Fortran
real(f4), parameter, public catchem_bridge_precision::tiny_real;
```




<hr>



### variable zero 

_Zero value (kind=fp)_ 
```Fortran
real(fp), parameter, public catchem_bridge_precision::zero;
```




<hr>



### variable zero\_dble 

_Zero value (kind=f8)_ 
```Fortran
real(f8), parameter, public catchem_bridge_precision::zero_dble;
```




<hr>



### variable zero\_real 

_Zero value (kind=f4)_ 
```Fortran
real(f4), parameter, public catchem_bridge_precision::zero_real;
```




<hr>
## Public Functions Documentation




### function exact\_equal\_f4 

_Exact equality (see the_ `exact_equal` _interface)_
```Fortran
elemental logical function catchem_bridge_precision::exact_equal_f4 (
    real( f4 ), intent(in) a,
    real( f4 ), intent(in) b
) 
```




<hr>



### function exact\_equal\_f8 

_Exact equality (see the_ `exact_equal` _interface)_
```Fortran
elemental logical function catchem_bridge_precision::exact_equal_f8 (
    real( f8 ), intent(in) a,
    real( f8 ), intent(in) b
) 
```




<hr>



### function is\_exact\_zero\_f4 

_Exact zero test (see the_ `is_exact_zero` _interface)_
```Fortran
elemental logical function catchem_bridge_precision::is_exact_zero_f4 (
    real( f4 ), intent(in) x
) 
```




<hr>



### function is\_exact\_zero\_f8 

_Exact zero test (see the_ `is_exact_zero` _interface)_
```Fortran
elemental logical function catchem_bridge_precision::is_exact_zero_f8 (
    real( f8 ), intent(in) x
) 
```




<hr>



### function rae\_f4 

_Real approximately equal:_ `abs(a - b) < tiny(a)` __
```Fortran
logical function catchem_bridge_precision::rae_f4 (
    real( f4 ), intent(in) a,
    real( f4 ), intent(in) b
) 
```




<hr>



### function rae\_f8 

_Real approximately equal:_ `abs(a - b) < tiny(a)` __
```Fortran
logical function catchem_bridge_precision::rae_f8 (
    real( f8 ), intent(in) a,
    real( f8 ), intent(in) b
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/catchem_bridge_utilities.F90`

