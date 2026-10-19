!> \file test_emis_time_weight.f90
!! \brief Unit tests for the pure emission time-interpolation helpers (feature 014, US4).
!!
!! These exercise the ESMF-free math in catchem_emis_time against hand-computed
!! reference weights, mirroring the MAPL ExtData / GOCART2G conventions:
!!   * catchem_emis_time_weight       - DOY axis with 365-day cyclic wrap
!!   * catchem_emis_file_time_bracket - record-timestamp anchoring (spanning vs cyclic)
!!   * catchem_emis_dt_weight         - daily-mean 12Z knots + daily_hold piecewise constant
program test_emis_time_weight
   use testing_mod, only: assert, assert_close
   use catchem_bridge_precision, only: fp
   use catchem_emis_time, only: catchem_emis_doy_fraction, catchem_emis_time_weight, &
      catchem_emis_file_time_bracket, catchem_emis_dt_weight, &
      CATCHEM_EMIS_NO_KEY

   implicit none

   real(fp) :: w
   integer :: lo
   integer :: i
   integer :: tc_dates(12), tc_secs(12)
   integer :: tc_dates14(14), tc_secs14(14)
   real(fp) :: dn
   real(fp), parameter :: tol = 1.0e-3_fp

   write(*,*) 'Testing emission time-weight helpers...'
   write(*,*) ''

   ! --- Test 1: day-of-year fraction ---------------------------------------
   write(*,*) 'Test 1: day-of-year fraction'
   ! Jan 1 00Z -> 0
   call assert_close(catchem_emis_doy_fraction(1, 1, 0), 0.0_fp, 1.0e-6_fp, &
      'Jan 1 00Z is day-of-year 0')
   ! Jan 20 00Z -> 19
   call assert_close(catchem_emis_doy_fraction(1, 20, 0), 19.0_fp, 1.0e-6_fp, &
      'Jan 20 is day-of-year 19')
   ! Feb 1 00Z -> 31 (non-leap axis)
   call assert_close(catchem_emis_doy_fraction(2, 1, 0), 31.0_fp, 1.0e-6_fp, &
      'Feb 1 is day-of-year 31')
   ! Dec 31 12Z -> 334 + 30 + 0.5 = 364.5
   dn = catchem_emis_doy_fraction(12, 31, 12 * 3600)
   call assert_close(dn, 364.5_fp, 1.0e-6_fp, 'Dec 31 12Z is day-of-year 364.5')
   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! --- Test 2: time weight, month-start records (MEGAN-style) --------------
   write(*,*) 'Test 2: time weight with month-start records'
   ! d1 = Jan 1 (20210101) 00Z, d2 = Feb 1 (20210201) 00Z, now = Jan 20 00Z.
   ! Expected: (19 - 0) / (31 - 0) = 19/31 ~= 0.6129.
   call catchem_emis_time_weight(catchem_emis_doy_fraction(1, 20, 0), &
      20210101, 0, 20210201, 0, w)
   call assert_close(w, 19.0_fp / 31.0_fp, tol, 'Jan 20 blends 19/31 toward Feb 1')
   ! Exactly at the lower knot -> 0.
   call catchem_emis_time_weight(catchem_emis_doy_fraction(1, 1, 0), &
      20210101, 0, 20210201, 0, w)
   call assert_close(w, 0.0_fp, tol, 'weight is 0 at the lower knot')
   ! Exactly at the upper knot -> 1.
   call catchem_emis_time_weight(catchem_emis_doy_fraction(2, 1, 0), &
      20210101, 0, 20210201, 0, w)
   call assert_close(w, 1.0_fp, tol, 'weight is 1 at the upper knot')
   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! --- Test 3: time weight, mid-month records (DMS-style) ------------------
   write(*,*) 'Test 3: time weight with mid-month records'
   ! d1 = 20111214 @ 12Z, d2 = 20120114 @ 12Z, now = Dec 20 00Z.
   ! doy(Dec14,12Z) = 334 + 13 + 0.5 = 347.5
   ! doy(Dec20,00Z) = 334 + 19      = 353
   ! doy(Jan14,12Z) = 0 + 13 + 0.5  = 13.5 -> wraps to 378.5
   ! w = (353 - 347.5) / (378.5 - 347.5) = 5.5 / 31 ~= 0.1774
   call catchem_emis_time_weight(catchem_emis_doy_fraction(12, 20, 0), &
      20111214, 12 * 3600, 20120114, 12 * 3600, w)
   call assert_close(w, 5.5_fp / 31.0_fp, tol, 'Dec 20 blends ~0.177 across the mid-month bracket')
   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! --- Test 4: time weight, year-end wrap ----------------------------------
   write(*,*) 'Test 4: time weight wraps across the year end'
   ! Dec 15 -> Jan 15 bracket, now = Dec 31.
   ! doy(Dec15) = 334 + 14 = 348, doy(Jan15) = 14 -> wraps to 379
   ! doy(Dec31) = 334 + 30 = 364
   ! w = (364 - 348) / (379 - 348) = 16/31 ~= 0.5161
   call catchem_emis_time_weight(catchem_emis_doy_fraction(12, 31, 0), &
      20201215, 0, 20210115, 0, w)
   call assert_close(w, 16.0_fp / 31.0_fp, tol, 'Dec 31 blends ~0.516 across the Dec->Jan wrap')
   ! Current time before the lower knot wraps too: now = Jan 5 with a
   ! Dec 15 -> Jan 15 bracket. doy(Jan5) = 4 < 348 -> +365 = 369.
   ! w = (369 - 348) / (379 - 348) = 21/31 ~= 0.6774
   call catchem_emis_time_weight(catchem_emis_doy_fraction(1, 5, 0), &
      20201215, 0, 20210115, 0, w)
   call assert_close(w, 21.0_fp / 31.0_fp, tol, 'Jan 5 wraps ahead of the Dec 15 knot')
   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! --- Test 5: cyclic climatology bracket (month-start, 12 records) --------
   write(*,*) 'Test 5: file-time bracket, cyclic climatology'
   ! A 12-record month-start climatology stamped 2018-01-01 .. 2018-12-01.
   do i = 1, 12
      tc_dates(i) = 20180101 + (i - 1) * 100
      tc_secs(i) = 0
   end do
   ! now = Dec 20 -> outside the file's years (2021 run) -> cyclic regime.
   ! doy(Dec20) = 353; doy(Dec 1) = 334 <= 353, doy(Jan 1)=0 is not > 353 in
   ! calendar order, so lo_rec = 12 (Dec), upper = 1 (Jan) with wrap.
   ! w = (353 - 334) / (365 + 0 - 334) = 19/31 ~= 0.6129
   dn = catchem_emis_doy_fraction(12, 20, 0)
   call catchem_emis_file_time_bracket(CATCHEM_EMIS_NO_KEY, dn, 12, tc_dates, tc_secs, lo, w)
   call assert(lo == 12, 'cyclic Dec bracket selects record 12 (Dec)')
   call assert_close(w, 19.0_fp / 31.0_fp, tol, 'cyclic Dec 20 weight is 19/31 toward Jan')
   ! now = Feb 10 -> doy(41): lo_rec = 2 (Feb 1, doy 31), upper = Mar 1 (doy 59)
   ! w = (40 - 31) / (59 - 31) = 9/28 ~= 0.321
   dn = catchem_emis_doy_fraction(2, 10, 0)
   call catchem_emis_file_time_bracket(CATCHEM_EMIS_NO_KEY, dn, 12, tc_dates, tc_secs, lo, w)
   call assert(lo == 2, 'cyclic Feb bracket selects record 2 (Feb)')
   call assert_close(w, 9.0_fp / 28.0_fp, tol, 'cyclic Feb 10 weight is 9/28 toward Mar')
   ! now = Jan 5 -> doy(4) < doy(record 1)=0? No: 4 > 0, so lo_rec = 1, upper = 2.
   ! w = (4 - 0) / (31 - 0) = 4/31
   dn = catchem_emis_doy_fraction(1, 5, 0)
   call catchem_emis_file_time_bracket(CATCHEM_EMIS_NO_KEY, dn, 12, tc_dates, tc_secs, lo, w)
   call assert(lo == 1, 'cyclic Jan bracket selects record 1 (Jan)')
   call assert_close(w, 4.0_fp / 31.0_fp, tol, 'cyclic Jan 5 weight is 4/31 toward Feb')
   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   ! --- Test 6: spanning file bracket (padded 14-record GMI-style) ----------
   write(*,*) 'Test 6: file-time bracket, spanning regime'
   ! 14 records: Dec2020, Jan2021 ... Dec2021, Jan2022 (month-start stamps).
   tc_dates14(1) = 20201201
   tc_secs14(1) = 0
   do i = 2, 13
      tc_dates14(i) = 20210101 + (i - 2) * 100
      tc_secs14(i) = 0
   end do
   tc_dates14(14) = 20220101
   tc_secs14(14) = 0
   ! now = Dec 20 2021 -> within [Dec2020, Jan2022] -> spanning regime.
   ! The last record <= now is record 13 (Dec 1 2021); upper = 14 (Jan 2022).
   ! Weight via the day-of-year axis: doy(Dec20)=353, doy(Dec1)=334, doy(Jan1)=0->365
   ! w = (353 - 334) / (365 - 334) = 19/31 ~= 0.6129
   dn = catchem_emis_doy_fraction(12, 20, 0)
   call catchem_emis_file_time_bracket(20211220_8 * 100000_8, dn, 14, tc_dates14, tc_secs14, lo, w)
   call assert(lo == 13, 'spanning Dec 20 2021 selects record 13 (Dec 2021)')
   call assert_close(w, 19.0_fp / 31.0_fp, tol, 'spanning Dec 20 weight is 19/31 toward Jan 2022')
   ! now = Feb 10 2021 -> spanning; last record <= now is record 3 (Feb 1 2021).
   ! doy(Feb10)=40, doy(Feb1)=31, doy(Mar1)=59 -> w = 9/28
   dn = catchem_emis_doy_fraction(2, 10, 0)
   call catchem_emis_file_time_bracket(20210210_8 * 100000_8, dn, 14, tc_dates14, tc_secs14, lo, w)
   call assert(lo == 3, 'spanning Feb 10 2021 selects record 3 (Feb 2021)')
   call assert_close(w, 9.0_fp / 28.0_fp, tol, 'spanning Feb 10 weight is 9/28 toward Mar')
   write(*,*) 'Test 6 passed!'
   write(*,*) ''

   ! --- Test 7: daily-mean 12Z knot weight + daily_hold ---------------------
   write(*,*) 'Test 7: daily weight and daily_hold'
   ! 06Z -> (6/24) - 0.5 = -0.25 -> wraps to 0.75: morning blends [yesterday,today].
   call assert_close(catchem_emis_dt_weight(6, 0, 0, .false.), 0.75_fp, 1.0e-6_fp, &
      '06Z daily weight is 0.75')
   ! 12Z -> 0.5 - 0.5 = 0.0
   call assert_close(catchem_emis_dt_weight(12, 0, 0, .false.), 0.0_fp, 1.0e-6_fp, &
      '12Z daily weight is 0')
   ! 18Z -> 0.75 - 0.5 = 0.25
   call assert_close(catchem_emis_dt_weight(18, 0, 0, .false.), 0.25_fp, 1.0e-6_fp, &
      '18Z daily weight is 0.25')
   ! 00Z -> 0 - 0.5 = -0.5 -> wraps to 0.5
   call assert_close(catchem_emis_dt_weight(0, 0, 0, .false.), 0.5_fp, 1.0e-6_fp, &
      '00Z daily weight is 0.5')
   ! daily_hold -> constant 0.5 regardless of hour
   call assert_close(catchem_emis_dt_weight(3, 0, 0, .true.), 0.5_fp, 1.0e-6_fp, &
      'daily_hold weight is 0.5 at 03Z')
   call assert_close(catchem_emis_dt_weight(23, 30, 0, .true.), 0.5_fp, 1.0e-6_fp, &
      'daily_hold weight is 0.5 at 23:30Z')
   write(*,*) 'Test 7 passed!'
   write(*,*) ''

   write(*,*) 'All emission time-weight tests passed!'
end program test_emis_time_weight
