

# Interface catchem\_bridge\_precision::is\_exact\_zero



[**ClassList**](annotated.md) **>** [**catchem\_bridge\_precision**](namespacecatchem__bridge__precision.md) **>** [**is\_exact\_zero**](interfacecatchem__bridge__precision_1_1is__exact__zero.md)



_Exact test for zero: exact\_equal(x, 0)._ [More...](#detailed-description)






































## Public Functions

| Type | Name |
| ---: | :--- |
|  elemental logical function | [**is\_exact\_zero\_f4**](#function-is_exact_zero_f4) (real([**f4**](namespacecatchem__bridge__precision.md#variable-f4)), intent(in) x) <br>_Exact zero test (see the_ `is_exact_zero` _interface)_ |
|  elemental logical function | [**is\_exact\_zero\_f8**](#function-is_exact_zero_f8) (real([**f8**](namespacecatchem__bridge__precision.md#variable-f8)), intent(in) x) <br>_Exact zero test (see the_ `is_exact_zero` _interface)_ |




























## Detailed Description


True for +0 and -0 only. 


    
## Public Functions Documentation




### function is\_exact\_zero\_f4 

_Exact zero test (see the_ `is_exact_zero` _interface)_
```Fortran
elemental logical function catchem_bridge_precision::is_exact_zero::is_exact_zero_f4 (
    real( f4 ), intent(in) x
) 
```




<hr>



### function is\_exact\_zero\_f8 

_Exact zero test (see the_ `is_exact_zero` _interface)_
```Fortran
elemental logical function catchem_bridge_precision::is_exact_zero::is_exact_zero_f8 (
    real( f8 ), intent(in) x
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/catchem_bridge_utilities.F90`

