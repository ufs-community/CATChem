

# Namespace catchem\_emis\_time



[**Namespace List**](namespaces.md) **>** [**catchem\_emis\_time**](namespacecatchem__emis__time.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  integer(8), parameter, public | [**catchem\_emis\_no\_key**](#variable-catchem_emis_no_key)   = `-1\_8`<br>_Sentinel datetime key meaning "no current-time key available"; forces the bracket helper into the cyclic-climatology (day-of-year) regime._  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  pure real(fp) function, public | [**catchem\_emis\_doy\_fraction**](#function-catchem_emis_doy_fraction) (integer, intent(in) month, integer, intent(in) day, integer, intent(in) secs) <br>_Fractional day-of-year (0-based) for (month, day, seconds-of-day)._  |
|  pure real(fp) function, public | [**catchem\_emis\_dt\_weight**](#function-catchem_emis_dt_weight) (integer, intent(in) hour, integer, intent(in) minute, integer, intent(in) second, logical, intent(in) daily\_hold) <br>_Daily-mean linear interpolation weight on 12Z knots._  |
|  pure subroutine, public | [**catchem\_emis\_file\_time\_bracket**](#function-catchem_emis_file_time_bracket) (integer(8), intent(in) curr\_key, real(fp), intent(in) curr\_dn, integer, intent(in) n, integer, dimension(:), intent(in) tc\_dates, integer, dimension(:), intent(in) tc\_secs, integer, intent(out) lo\_rec, real(fp), intent(out) w\_next) <br>_Monthly bracket anchored at the file's actual record timestamps._  |
|  pure subroutine, public | [**catchem\_emis\_time\_weight**](#function-catchem_emis_time_weight) (real(fp), intent(in) curr\_dn, integer, intent(in) d1, integer, intent(in) s1, integer, intent(in) d2, integer, intent(in) s2, real(fp), intent(out) w\_next) <br>_Linear interpolation weight between two record valid-times._  |




























## Public Attributes Documentation




### variable catchem\_emis\_no\_key 

_Sentinel datetime key meaning "no current-time key available"; forces the bracket helper into the cyclic-climatology (day-of-year) regime._ 
```Fortran
integer(8), parameter, public catchem_emis_time::catchem_emis_no_key;
```



The key is 8-byte because yyyymmdd\*1e5 exceeds the 32-bit integer range. 


        

<hr>
## Public Functions Documentation




### function catchem\_emis\_doy\_fraction 

_Fractional day-of-year (0-based) for (month, day, seconds-of-day)._ 
```Fortran
pure real(fp) function, public catchem_emis_time::catchem_emis_doy_fraction (
    integer, intent(in) month,
    integer, intent(in) day,
    integer, intent(in) secs
) 
```



Month outside 1..12 is clamped to January. A non-leap 365-day axis is used, consistent with the climatological weighting convention. 


        

<hr>



### function catchem\_emis\_dt\_weight 

_Daily-mean linear interpolation weight on 12Z knots._ 
```Fortran
pure real(fp) function, public catchem_emis_time::catchem_emis_dt_weight (
    integer, intent(in) hour,
    integer, intent(in) minute,
    integer, intent(in) second,
    logical, intent(in) daily_hold
) 
```



Daily-mean files are valid at 12Z (GEOS/ExtData), so the two bracketing knots are consecutive-day 12Z means. The fraction-of-day is shifted by half a day and wrapped into [0,1): the morning blends [yesterday, today] and the afternoon [today, tomorrow], matching MAPL ExtData rather than a 00Z ramp. With `daily_hold` (MAPL ExtData refresh-cadence match) the [D-1, D] bracket is held constant for the whole day at the 00Z blend (weight 0.5), giving a piecewise-constant daily value 0.5\*(emis(D-1) + emis(D)).




**Parameters:**


* `hour` current clock time (UTC) 
* `daily_hold` hold the bracket piecewise-constant 



**Returns:**

weight of the upper (later) knot, in [0,1) 





        

<hr>



### function catchem\_emis\_file\_time\_bracket 

_Monthly bracket anchored at the file's actual record timestamps._ 
```Fortran
pure subroutine, public catchem_emis_time::catchem_emis_file_time_bracket (
    integer(8), intent(in) curr_key,
    real(fp), intent(in) curr_dn,
    integer, intent(in) n,
    integer, dimension(:), intent(in) tc_dates,
    integer, dimension(:), intent(in) tc_secs,
    integer, intent(out) lo_rec,
    real(fp), intent(out) w_next
) 
```



Selects the lower bracketing record (1-based) and the upper-record weight from a multi-record climatology, honouring the records' own time coordinate. Two regimes:
* spanning file: the current datetime key lies within [first,last] record. Bracket by the records' real datetime (integer key = date\*1e5 + secs, so it is strictly ordered). Needed for files that repeat months, e.g. the padded 14-record GMI oxidant file.
* cyclic climatology: the current time is OUTSIDE the file's years (a 12-record climatology used in a later run), or no current key is supplied (CATCHEM\_EMIS\_NO\_KEY). Bracket by day-of-year, wrapping Dec-&gt;Jan. The upper record used downstream is lo\_rec+1 (wrapping to 1 for the cyclic Dec-&gt;Jan case). The weight itself is delegated to catchem\_emis\_time\_weight.






**Parameters:**


* `curr_key` current datetime key (yyyymmdd\*1e5 + secs) or NO\_KEY 
* `curr_dn` current fractional day-of-year 
* `n` number of records 
* `tc_dates` yyyymmdd for each record 
* `tc_secs` seconds-of-day for each record 
* `lo_rec` 1-based lower bracketing record index 
* `w_next` weight of the upper record, clamped to [0,1] 




        

<hr>



### function catchem\_emis\_time\_weight 

_Linear interpolation weight between two record valid-times._ 
```Fortran
pure subroutine, public catchem_emis_time::catchem_emis_time_weight (
    real(fp), intent(in) curr_dn,
    integer, intent(in) d1,
    integer, intent(in) s1,
    integer, intent(in) d2,
    integer, intent(in) s2,
    real(fp), intent(out) w_next
) 
```



Storage-independent: given two bracketing record timestamps (yyyymmdd, seconds-of-day) and the current fractional day-of-year, return the weight `w_next` (0..1) of the UPPER record so value = (1 - w\_next) \* record1 + w\_next \* record2. Every time is mapped onto a fixed 365-day day-of-year axis and wrapped cyclically across the year end, so it works for single-file records, multi-file bracket times, and any record stamping. For ~monthly record spacing the day-of-year fraction equals the true elapsed-time fraction.




**Parameters:**


* `curr_dn` current fractional day-of-year (catchem\_emis\_doy\_fraction) 
* `d1` lower record date (yyyymmdd) and seconds-of-day 
* `d2` upper record date (yyyymmdd) and seconds-of-day 
* `w_next` weight of the upper record, clamped to [0,1] 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/catchem_emis_time_mod.F90`

