!> \file test_TimeState.f90
!! \brief Test program for TimeState module
!!
!!!>
program test_TimeState
   use testing_mod, only: assert, assert_close
   use TimeState_Mod
   use Error_Mod, only: ErrorManagerType, CC_SUCCESS

   implicit none

   type(TimeStateType) :: time_state
   type(ErrorManagerType), pointer :: error_mgr
   integer :: year, month, day
   real :: jd
   integer :: doy
   real :: sza
   integer :: rc

   write(*,*) 'Testing TimeState module...'
   write(*,*) ''

   ! Test 1: Error manager initialization
   write(*,*) 'Test 1: Error manager initialization'
   allocate(error_mgr)
   call error_mgr%init()

   ! Test 2: TimeState initialization (default values)
   write(*,*) 'Test 2: TimeState initialization (default values)'
   call time_state%init(error_mgr=error_mgr, rc=rc)
   call assert(rc == CC_SUCCESS, "TimeState initialization should succeed")

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: Default values
   write(*,*) 'Test 3: Default values'
   call time_state%get_current_date(year, month, day)
   jd = time_state%get_julian_date()
   doy = time_state%get_doy()

   call assert(year == 2000, "Default year should be 2000")
   call assert(month == 1, "Default month should be 1")
   call assert(day == 1, "Default day should be 1")
   call assert(doy == 1, "Default day of year should be 1")

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! Test 4: Solar zenith angle calculation
   write(*,*) 'Test 4: Solar zenith angle calculation'
   ! Calculate SZA for noon at equator on day 1
   sza = time_state%get_sza(0.0, 0.0)  ! lat=0, lon=0
   call assert(sza >= 0.0 .and. sza <= 90.0, "SZA should be between 0 and 90 degrees")

   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! Test 5: Timestep
   write(*,*) 'Test 5: Timestep'
   block
      real :: dt
      dt = time_state%get_timestep()
      call assert_close(dt, 3600.0, 1e-6, "Default timestep should be 3600 seconds")
   end block

   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   ! Test 6: Time string functions
   write(*,*) 'Test 6: Time string functions'
   block
      character(len=25) :: iso_time
      character(len=25) :: human_time
      character(len=16) :: compact_time

      iso_time = time_state%get_time_iso8601()
      human_time = time_state%get_time_human()
      compact_time = time_state%get_time_compact()

      ! Just check that they return non-empty strings
      call assert(len_trim(iso_time) > 0, "ISO time should not be empty")
      call assert(len_trim(human_time) > 0, "Human time should not be empty")
      call assert(len_trim(compact_time) > 0, "Compact time should not be empty")
   end block

   write(*,*) 'Test 6 passed!'
   write(*,*) ''

   ! Test 7: Custom initialization
   write(*,*) 'Test 7: Custom initialization'
   call time_state%init(year=2023, month=6, day=15, hour=12, minute=30, second=45, timestep=1800.0, error_mgr=error_mgr, rc=rc)
   call assert(rc == CC_SUCCESS, "Custom TimeState initialization should succeed")

   call time_state%get_current_date(year, month, day)
   call assert(year == 2023, "Custom year should be 2023")
   call assert(month == 6, "Custom month should be 6")
   call assert(day == 15, "Custom day should be 15")
   call assert_close(time_state%get_timestep(), 1800.0, 1e-6, "Custom timestep should be 1800 seconds")

   write(*,*) 'Test 7 passed!'
   write(*,*) ''

   ! Test 8: Validation
   write(*,*) 'Test 8: Validation'
   call time_state%validate(error_mgr, rc)
   call assert(rc == CC_SUCCESS, "TimeState validation should succeed")

   write(*,*) 'Test 8 passed!'
   write(*,*) ''

   ! Test 9: Reset functionality
   write(*,*) 'Test 9: Reset functionality'
   call time_state%reset(error_mgr, rc)
   call assert(rc == CC_SUCCESS, "TimeState reset should succeed")

   call time_state%get_current_date(year, month, day)
   call assert(year == 2000, "Reset year should be 2000")
   call assert(month == 1, "Reset month should be 1")
   call assert(day == 1, "Reset day should be 1")

   write(*,*) 'Test 9 passed!'
   write(*,*) ''

   ! Test 10: State status
   write(*,*) 'Test 10: State status'
   block
      logical :: is_ready
      is_ready = time_state%is_ready()
      call assert(is_ready, "TimeState should be ready after initialization")
   end block

   write(*,*) 'Test 10 passed!'
   write(*,*) ''

   ! Test 11: Memory usage
   write(*,*) 'Test 11: Memory usage'
   block
      integer(8) :: memory_bytes
      memory_bytes = time_state%get_memory_usage()
      call assert(memory_bytes > 0, "Memory usage should be positive")
   end block

   write(*,*) 'Test 11 passed!'
   write(*,*) ''

   ! Test 12: Timezone offset
   write(*,*) 'Test 12: Timezone offset'
   block
      integer :: tz_offset
      tz_offset = time_state%get_timezone_offset(0.0)   ! UTC
      call assert(tz_offset == 0, "UTC timezone offset should be 0")

      tz_offset = time_state%get_timezone_offset(-75.0)  ! Eastern US
      call assert(tz_offset == -5, "Eastern US timezone offset should be -5")

      tz_offset = time_state%get_timezone_offset(120.0)  ! China
      call assert(tz_offset == 8, "China timezone offset should be 8")
   end block

   write(*,*) 'Test 12 passed!'
   write(*,*) ''

   ! Test 13: Cleanup
   write(*,*) 'Test 13: Cleanup'
   call time_state%cleanup(error_mgr, rc)
   call assert(rc == CC_SUCCESS, "TimeState cleanup should succeed")

   ! Verify cleanup worked
   block
      logical :: is_ready
      is_ready = time_state%is_ready()
      call assert(.not. is_ready, "TimeState should not be ready after cleanup")
   end block

   write(*,*) 'Test 13 passed!'
   write(*,*) ''

   write(*,*) 'All TimeState tests passed!'

end program test_TimeState
