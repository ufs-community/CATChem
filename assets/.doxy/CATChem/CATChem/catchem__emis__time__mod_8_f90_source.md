

# File catchem\_emis\_time\_mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**catchem\_emis\_time\_mod.F90**](catchem__emis__time__mod_8_f90.md)

[Go to the documentation of this file](catchem__emis__time__mod_8_f90.md)


```Fortran

module catchem_emis_time
   use catchem_bridge_precision, only: fp
   implicit none
   private

   public :: catchem_emis_doy_fraction
   public :: catchem_emis_time_weight
   public :: catchem_emis_file_time_bracket
   public :: catchem_emis_dt_weight
   public :: catchem_emis_no_key

   integer(8), parameter :: CATCHEM_EMIS_NO_KEY = -1_8

contains

   pure real(fp) function catchem_emis_doy_fraction(month, day, secs) result(dn)
      integer, intent(in) :: month
      integer, intent(in) :: day
      integer, intent(in) :: secs
      integer, parameter :: cum(12) = (/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/)
      integer :: m

      m = month
      if (m < 1 .or. m > 12) m = 1
      dn = real(cum(m), fp) + real(day - 1, fp) + real(secs, fp) / 86400.0_fp
   end function catchem_emis_doy_fraction

   pure subroutine catchem_emis_time_weight(curr_dn, d1, s1, d2, s2, w_next)
      real(fp), intent(in)  :: curr_dn
      integer,  intent(in)  :: d1
      integer,  intent(in)  :: s1
      integer,  intent(in)  :: d2
      integer,  intent(in)  :: s2
      real(fp), intent(out) :: w_next
      real(fp), parameter :: YEARLEN = 365.0_fp
      real(fp) :: dn, d1f, d2f

      d1f = catchem_emis_doy_fraction(mod(d1 / 100, 100), mod(d1, 100), s1)
      d2f = catchem_emis_doy_fraction(mod(d2 / 100, 100), mod(d2, 100), s2)
      dn = curr_dn
      if (d2f <= d1f) d2f = d2f + yearlen   ! upper record wraps past year end
      if (dn < d1f) dn = dn + yearlen       ! current time before the lower knot
      w_next = 0.0_fp
      if (d2f > d1f) w_next = (dn - d1f) / (d2f - d1f)
      if (w_next < 0.0_fp) w_next = 0.0_fp
      if (w_next > 1.0_fp) w_next = 1.0_fp
   end subroutine catchem_emis_time_weight

   pure subroutine catchem_emis_file_time_bracket(curr_key, curr_dn, n, tc_dates, tc_secs, &
      lo_rec, w_next)
      integer(8),  intent(in)  :: curr_key
      real(fp),    intent(in)  :: curr_dn
      integer,     intent(in)  :: n
      integer,     intent(in)  :: tc_dates(:)
      integer,     intent(in)  :: tc_secs(:)
      integer,     intent(out) :: lo_rec
      real(fp),    intent(out) :: w_next
      integer(8) :: curr_key8, key_i, key_lo, key_hi
      real(fp) :: doy_i
      integer :: i, up
      logical :: spanning

      lo_rec = 1
      w_next = 0.0_fp
      if (n < 2) return

      curr_key8 = curr_key
      key_lo = int(tc_dates(1), 8) * 100000_8 + int(tc_secs(1), 8)
      key_hi = int(tc_dates(n), 8) * 100000_8 + int(tc_secs(n), 8)
      spanning = (curr_key8 /= int(catchem_emis_no_key, 8)) .and. &
         (curr_key8 >= key_lo .and. curr_key8 <= key_hi)

      if (spanning) then
         ! Bracket by actual record datetime (padded/spanning file).
         lo_rec = 1
         do i = 1, n
            key_i = int(tc_dates(i), 8) * 100000_8 + int(tc_secs(i), 8)
            if (key_i <= curr_key8) then
               lo_rec = i
            else
               exit
            end if
         end do
      else
         ! Cyclic climatology: day-of-year bracket (records in calendar order).
         lo_rec = 0
         do i = 1, n
            doy_i = catchem_emis_doy_fraction(mod(tc_dates(i) / 100, 100), &
               mod(tc_dates(i), 100), tc_secs(i))
            if (doy_i <= curr_dn) then
               lo_rec = i
            else
               exit
            end if
         end do
         if (lo_rec == 0) lo_rec = n   ! before the first record -> wrap to last
      end if

      up = lo_rec + 1
      if (up > n) up = 1
      call catchem_emis_time_weight(curr_dn, tc_dates(lo_rec), tc_secs(lo_rec), &
         tc_dates(up), tc_secs(up), w_next)
   end subroutine catchem_emis_file_time_bracket

   pure real(fp) function catchem_emis_dt_weight(hour, minute, second, daily_hold) result(w_next)
      integer, intent(in) :: hour
      integer, intent(in) :: minute
      integer, intent(in) :: second
      logical, intent(in) :: daily_hold
      real(fp) :: w_local

      if (daily_hold) then
         w_next = 0.5_fp
         return
      end if
      w_local = (real(hour, fp) + real(minute, fp) / 60.0_fp + &
         real(second, fp) / 3600.0_fp) / 24.0_fp - 0.5_fp
      if (w_local < 0.0_fp) w_local = w_local + 1.0_fp
      w_next = w_local
   end function catchem_emis_dt_weight

end module catchem_emis_time
```


