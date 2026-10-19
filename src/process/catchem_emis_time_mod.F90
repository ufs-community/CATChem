!> \file catchem_emis_time_mod.F90
!! \brief Pure, ESMF-free temporal-interpolation helpers for external emissions.
!!
!! These routines reproduce the MAPL ExtData / GOCART2G emission time handling
!! (feature 014, US4) so the interpolation math can be unit-tested without
!! pulling in the ESMF framework.  The NUOPC emission driver
!! (drivers/nuopc/catchem_emis_mod.F90) converts an ESMF_Time into the plain
!! (day-of-year, integer datetime-key) arguments below and delegates here.
!!
!! Two behaviours live here:
!!   * catchem_emis_time_weight      - weight between two record valid-times on a
!!                                     fixed 365-day day-of-year axis with cyclic
!!                                     year-end wrap (climatological anchoring).
!!   * catchem_emis_file_time_bracket - select the lower bracketing record and the
!!                                     upper-record weight for a monthly climatology
!!                                     anchored at the file's own record timestamps
!!                                     ('file' anchoring), spanning vs cyclic regimes.
!!   * catchem_emis_dt_weight        - daily-mean linear weight on 12Z knots, with the
!!                                     MAPL 'daily_hold' piecewise-constant variant.
!!
!! All routines are pure and depend only on the bridge precision kind.
module catchem_emis_time
   use catchem_bridge_precision, only: fp
   implicit none
   private

   public :: catchem_emis_doy_fraction
   public :: catchem_emis_time_weight
   public :: catchem_emis_file_time_bracket
   public :: catchem_emis_dt_weight
   public :: CATCHEM_EMIS_NO_KEY

   !> Sentinel datetime key meaning "no current-time key available"; forces the
   !! bracket helper into the cyclic-climatology (day-of-year) regime.  The key is
   !! 8-byte because yyyymmdd*1e5 exceeds the 32-bit integer range.
   integer(8), parameter :: CATCHEM_EMIS_NO_KEY = -1_8

contains

   !> \brief Fractional day-of-year (0-based) for (month, day, seconds-of-day).
   !! Month outside 1..12 is clamped to January.  A non-leap 365-day axis is used,
   !! consistent with the climatological weighting convention.
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

   !> \brief Linear interpolation weight between two record valid-times.
   !!
   !! Storage-independent: given two bracketing record timestamps (yyyymmdd,
   !! seconds-of-day) and the current fractional day-of-year, return the weight
   !! `w_next` (0..1) of the UPPER record so
   !!   value = (1 - w_next) * record1 + w_next * record2.
   !! Every time is mapped onto a fixed 365-day day-of-year axis and wrapped
   !! cyclically across the year end, so it works for single-file records,
   !! multi-file bracket times, and any record stamping.  For ~monthly record
   !! spacing the day-of-year fraction equals the true elapsed-time fraction.
   !!
   !! \param[in]  curr_dn  current fractional day-of-year (catchem_emis_doy_fraction)
   !! \param[in]  d1, s1   lower record date (yyyymmdd) and seconds-of-day
   !! \param[in]  d2, s2   upper record date (yyyymmdd) and seconds-of-day
   !! \param[out] w_next   weight of the upper record, clamped to [0,1]
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
      if (d2f <= d1f) d2f = d2f + YEARLEN   ! upper record wraps past year end
      if (dn < d1f) dn = dn + YEARLEN       ! current time before the lower knot
      w_next = 0.0_fp
      if (d2f > d1f) w_next = (dn - d1f) / (d2f - d1f)
      if (w_next < 0.0_fp) w_next = 0.0_fp
      if (w_next > 1.0_fp) w_next = 1.0_fp
   end subroutine catchem_emis_time_weight

   !> \brief Monthly bracket anchored at the file's actual record timestamps.
   !!
   !! Selects the lower bracketing record (1-based) and the upper-record weight
   !! from a multi-record climatology, honouring the records' own time coordinate.
   !! Two regimes:
   !!   * spanning file: the current datetime key lies within [first,last] record.
   !!     Bracket by the records' real datetime (integer key = date*1e5 + secs, so
   !!     it is strictly ordered).  Needed for files that repeat months, e.g. the
   !!     padded 14-record GMI oxidant file.
   !!   * cyclic climatology: the current time is OUTSIDE the file's years (a
   !!     12-record climatology used in a later run), or no current key is supplied
   !!     (CATCHEM_EMIS_NO_KEY).  Bracket by day-of-year, wrapping Dec->Jan.
   !! The upper record used downstream is lo_rec+1 (wrapping to 1 for the cyclic
   !! Dec->Jan case).  The weight itself is delegated to catchem_emis_time_weight.
   !!
   !! \param[in]  curr_key  current datetime key (yyyymmdd*1e5 + secs) or NO_KEY
   !! \param[in]  curr_dn   current fractional day-of-year
   !! \param[in]  n         number of records
   !! \param[in]  tc_dates  yyyymmdd for each record
   !! \param[in]  tc_secs   seconds-of-day for each record
   !! \param[out] lo_rec    1-based lower bracketing record index
   !! \param[out] w_next    weight of the upper record, clamped to [0,1]
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
      spanning = (curr_key8 /= int(CATCHEM_EMIS_NO_KEY, 8)) .and. &
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

   !> \brief Daily-mean linear interpolation weight on 12Z knots.
   !!
   !! Daily-mean files are valid at 12Z (GEOS/ExtData), so the two bracketing knots
   !! are consecutive-day 12Z means.  The fraction-of-day is shifted by half a day
   !! and wrapped into [0,1): the morning blends [yesterday, today] and the
   !! afternoon [today, tomorrow], matching MAPL ExtData rather than a 00Z ramp.
   !! With `daily_hold` (MAPL ExtData refresh-cadence match) the [D-1, D] bracket is
   !! held constant for the whole day at the 00Z blend (weight 0.5), giving a
   !! piecewise-constant daily value 0.5*(emis(D-1) + emis(D)).
   !!
   !! \param[in]  hour, minute, second  current clock time (UTC)
   !! \param[in]  daily_hold            hold the bracket piecewise-constant
   !! \returns    weight of the upper (later) knot, in [0,1)
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
