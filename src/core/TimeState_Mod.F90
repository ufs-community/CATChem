!> \file TimeState_Mod.F90
!! \brief Time state and common time/solar functions for atmospheric chemistry
!!
!! Provides timekeeping, solar zenith angle, and calendar utilities for CATChem.
!!
module TimeState_Mod
   use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE
   use constants, only: PI, PI_180
   implicit none
   private
   public :: TimeStateType, is_global_holiday, is_us_holiday

   ! Local status constants to avoid circular dependency
   integer, parameter :: STATE_STATUS_UNINITIALIZED = 0
   integer, parameter :: STATE_STATUS_INITIALIZED = 1

   !> \brief Time state for model
   type :: TimeStateType
      integer :: year = 2000
      integer :: month = 1
      integer :: day = 1
      integer :: hour = 0
      integer :: minute = 0
      integer :: second = 0
      real    :: timestep = 3600.0 !< seconds
      real    :: julian_date = 0.0
      integer :: doy = 1
   contains
      procedure :: get_sza
      procedure :: get_cos_sza
      procedure :: get_timestep
      procedure :: get_current_date
      procedure :: get_julian_date
      procedure :: get_doy
      procedure :: init => timestate_init
      procedure :: validate => timestate_validate
      procedure :: cleanup => timestate_cleanup
      procedure :: reset => timestate_reset
      procedure :: get_status => timestate_get_status
      procedure :: get_memory_usage => timestate_get_memory_usage
      procedure :: print_info => timestate_print_info
      procedure :: is_ready => timestate_is_ready
      procedure :: get_time_iso8601
      procedure :: get_time_human
      procedure :: get_time_compact
      procedure :: get_timezone_offset
      procedure, private :: calculate_derived_fields
   end type TimeStateType

contains

   !> \brief Compute solar zenith angle (degrees) using latitude, longitude, and time of day
   real function get_cos_sza(this, lat, lon, mid_timestep) result(cos_sza_val)
      class(TimeStateType), intent(in) :: this
      real, intent(in) :: lat, lon
      logical, intent(in), optional :: mid_timestep
      ! Accurate solar zenith angle calculation
      ! Inputs: lat, lon in degrees; time from this%hour, this%minute, this%second; day of year from this%doy
      real :: lat_rad, lon_rad, decl_rad, ha_rad
      real :: decl, eqtime, time_offset, tst, ha
      !real :: cos_sza_val
      real :: fractional_hour, gamma

      ! Convert latitude and longitude to radians
      lat_rad = lat * PI_180
      lon_rad = lon * PI_180

      ! Calculate fractional hour of the day (UTC)
      if (present(mid_timestep)) then
         if (mid_timestep) then
            fractional_hour = real(this%hour) + real(this%minute)/60.0 + real(this%second)/3600.0 + (this%timestep/2.0)/3600.0
         end if
      else
         fractional_hour = real(this%hour) + real(this%minute)/60.0 + real(this%second)/3600.0
      end if

      ! Calculate day angle (in radians)
      gamma = 2.0 * PI * (real(this%doy) - 1.0) / 365.0

      ! Solar declination (in degrees, then radians)
      !decl = 23.44 * sin(2.0 * PI * (real(this%doy) - 81.0) / 365.0)
      !decl_rad = decl * PI_180

      !use a more accurate formula for declination
      decl = 0.006918 - 0.399912*cos(gamma) + 0.070257*sin(gamma) &
         - 0.006758*cos(2.0*gamma) + 0.000907*sin(2.0*gamma) &
         - 0.002697*cos(3.0*gamma) + 0.00148*sin(3.0*gamma)
      decl_rad = decl

      ! Equation of time (in minutes).
      eqtime = 229.18 * (0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) \
      - 0.014615 * cos(2.0*gamma) - 0.040849 * sin(2.0*gamma))

      ! Time offset (in minutes). Note here we assume longitude between -180 and 180 degrees
      time_offset = eqtime + 4.0 * lon

      ! True solar time (in minutes)
      tst = fractional_hour * 60.0 + time_offset

      ! Hour angle (in degrees, then radians)
      ha = (tst / 4.0) - 180.0
      ha_rad = ha * PI_180

      ! Solar zenith angle calculation
      cos_sza_val = sin(lat_rad) * sin(decl_rad) + cos(lat_rad) * cos(decl_rad) * cos(ha_rad)
      cos_sza_val = max(-1.0, min(1.0, cos_sza_val)) ! Clamp for safety
   end function get_cos_sza

   !> \brief Compute cosine of solar zenith angle
   real function get_sza(this, lat, lon) result(sza)
      class(TimeStateType), intent(in) :: this
      real, intent(in) :: lat, lon
      real :: cos_sza_val
      cos_sza_val = this%get_cos_sza(lat, lon)
      sza = acos(cos_sza_val) / PI_180
      sza = min(max(sza, 0.0), 90.0) ! clamp to [0, 90] (daylight) degrees
   end function get_sza

   !> \brief Get model timestep (seconds)
   real function get_timestep(this) result(dt)
      class(TimeStateType), intent(in) :: this
      dt = this%timestep
   end function get_timestep

   !> \brief Get current date as (year, month, day)
   subroutine get_current_date(this, year, month, day)
      class(TimeStateType), intent(in) :: this
      integer, intent(out) :: year, month, day
      year = this%year
      month = this%month
      day = this%day
   end subroutine get_current_date

   !> \brief Get Julian date
   real function get_julian_date(this) result(jd)
      class(TimeStateType), intent(in) :: this
      jd = this%julian_date
   end function get_julian_date

   !> \brief Get day of year
   integer function get_doy(this) result(doy)
      class(TimeStateType), intent(in) :: this
      doy = this%doy
   end function get_doy

   !> \brief Initialize TimeStateType
   !! Sets default values, calculates derived quantities (Julian date, DOY),
   !! and validates the initial time configuration.
   subroutine timestate_init(this, year, month, day, hour, minute, second, timestep, error_mgr, rc)
      use error_mod, only: ERROR_INVALID_INPUT

      class(TimeStateType), intent(inout) :: this
      integer, optional, intent(in) :: year, month, day, hour, minute, second
      real, optional, intent(in) :: timestep
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      character(len=256) :: thisLoc

      thisLoc = 'timestate_init (in core/TimeState_Mod.F90)'
      call error_mgr%push_context('timestate_init', 'initializing time state')

      rc = CC_SUCCESS

      ! Set time components with defaults or provided values
      this%year = 2000
      this%month = 1
      this%day = 1
      this%hour = 0
      this%minute = 0
      this%second = 0
      this%timestep = 3600.0  ! 1 hour default

      if (present(year)) this%year = year
      if (present(month)) this%month = month
      if (present(day)) this%day = day
      if (present(hour)) this%hour = hour
      if (present(minute)) this%minute = minute
      if (present(second)) this%second = second
      if (present(timestep)) this%timestep = timestep

      ! Calculate derived quantities
      call this%calculate_derived_fields(error_mgr, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%pop_context()
         return
      endif

      ! Validate the initialized state
      call this%validate(error_mgr, rc)

      call error_mgr%pop_context()
   end subroutine timestate_init

   !> \brief Validate TimeStateType
   !! Checks that all time components are within valid ranges and
   !! that derived quantities are consistent.
   subroutine timestate_validate(this, error_mgr, rc)
      use error_mod, only: ERROR_INVALID_INPUT

      class(TimeStateType), intent(in) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      character(len=256) :: thisLoc
      integer :: days_in_month

      thisLoc = 'timestate_validate (in core/TimeState_Mod.F90)'
      call error_mgr%push_context('timestate_validate', 'validating time state')

      rc = CC_SUCCESS

      ! Validate year (reasonable range)
      if (this%year < 1900 .or. this%year > 2200) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
            'Year must be between 1900 and 2200', rc, &
            thisLoc, 'Adjust year to reasonable range')
         call error_mgr%pop_context()
         return
      endif

      ! Validate month
      if (this%month < 1 .or. this%month > 12) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
            'Month must be between 1 and 12', rc, &
            thisLoc, 'Set month to valid range [1-12]')
         call error_mgr%pop_context()
         return
      endif

      ! Validate day based on month and leap year
      days_in_month = get_days_in_month(this%month, this%year)
      if (this%day < 1 .or. this%day > days_in_month) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
            'Day is invalid for the given month/year', rc, &
            thisLoc, 'Set day to valid range for the month')
         call error_mgr%pop_context()
         return
      endif

      ! Validate hour
      if (this%hour < 0 .or. this%hour > 23) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
            'Hour must be between 0 and 23', rc, &
            thisLoc, 'Set hour to valid range [0-23]')
         call error_mgr%pop_context()
         return
      endif

      ! Validate minute
      if (this%minute < 0 .or. this%minute > 59) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
            'Minute must be between 0 and 59', rc, &
            thisLoc, 'Set minute to valid range [0-59]')
         call error_mgr%pop_context()
         return
      endif

      ! Validate second
      if (this%second < 0 .or. this%second > 59) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
            'Second must be between 0 and 59', rc, &
            thisLoc, 'Set second to valid range [0-59]')
         call error_mgr%pop_context()
         return
      endif

      ! Validate timestep
      if (this%timestep <= 0.0) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
            'Timestep must be positive', rc, &
            thisLoc, 'Set timestep to positive value in seconds')
         call error_mgr%pop_context()
         return
      endif

      ! Validate day of year
      if (this%doy < 1 .or. this%doy > 366) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
            'Day of year must be between 1 and 366', rc, &
            thisLoc, 'Check DOY calculation')
         call error_mgr%pop_context()
         return
      endif

      call error_mgr%pop_context()
   end subroutine timestate_validate

   !> \brief Cleanup TimeStateType
   !! Resets all time state variables to uninitialized values.
   !! This prepares the object for deallocation or reinitialization.
   subroutine timestate_cleanup(this, error_mgr, rc)
      class(TimeStateType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      character(len=256) :: thisLoc

      thisLoc = 'timestate_cleanup (in core/TimeState_Mod.F90)'
      call error_mgr%push_context('timestate_cleanup', 'cleaning up time state')

      rc = CC_SUCCESS

      ! Reset all time components to invalid/uninitialized values
      this%year = -1
      this%month = -1
      this%day = -1
      this%hour = -1
      this%minute = -1
      this%second = -1
      this%timestep = -1.0
      this%julian_date = -1.0
      this%doy = -1

      call error_mgr%pop_context()
   end subroutine timestate_cleanup

   !> \brief Reset TimeStateType
   !! Resets time state to default values (Jan 1, 2000, 00:00:00)
   !! and recalculates derived quantities.
   subroutine timestate_reset(this, error_mgr, rc)
      class(TimeStateType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      character(len=256) :: thisLoc

      thisLoc = 'timestate_reset (in core/TimeState_Mod.F90)'
      call error_mgr%push_context('timestate_reset', 'resetting time state')

      rc = CC_SUCCESS

      ! Reset to default values
      this%year = 2000
      this%month = 1
      this%day = 1
      this%hour = 0
      this%minute = 0
      this%second = 0
      this%timestep = 3600.0  ! 1 hour

      ! Recalculate derived quantities
      call this%calculate_derived_fields(error_mgr, rc)

      call error_mgr%pop_context()
   end subroutine timestate_reset

   !> \brief Get status of TimeStateType
   !! Returns initialized status if all time components are valid
   function timestate_get_status(this) result(status)
      class(TimeStateType), intent(in) :: this
      integer :: status

      ! Check if time state has been properly initialized
      if (this%year > 0 .and. this%month > 0 .and. this%day > 0 .and. &
         this%hour >= 0 .and. this%minute >= 0 .and. this%second >= 0 .and. &
         this%timestep > 0.0 .and. this%doy > 0) then
         status = STATE_STATUS_INITIALIZED
      else
         status = STATE_STATUS_UNINITIALIZED
      endif
   end function timestate_get_status

   !> \brief Get memory usage of TimeStateType
   !! Returns estimated memory usage in bytes for the time state object
   function timestate_get_memory_usage(this) result(memory_bytes)
      class(TimeStateType), intent(in) :: this
      integer(8) :: memory_bytes

      ! Estimate memory usage:
      ! 6 integers (year, month, day, hour, minute, second, doy) = 6 * 4 = 24 bytes
      ! 2 reals (timestep, julian_date) = 2 * 4 = 8 bytes (assuming real is 4 bytes)
      ! Total ≈ 32 bytes
      memory_bytes = 32_8
   end function timestate_get_memory_usage

   !> \brief Print info for TimeStateType
   subroutine timestate_print_info(this, unit)
      class(TimeStateType), intent(in) :: this
      integer, optional, intent(in) :: unit
      if (present(unit)) then
         write(unit,*) 'TimeStateType: ', this%year, this%month, this%day, this%hour, this%minute, this%second
      else
         print *, 'TimeStateType: ', this%year, this%month, this%day, this%hour, this%minute, this%second
      end if
   end subroutine timestate_print_info

   !> \brief Is TimeStateType ready?
   !! Returns true if the time state is properly initialized and valid
   function timestate_is_ready(this) result(ready)
      class(TimeStateType), intent(in) :: this
      logical :: ready

      ready = (this%get_status() == STATE_STATUS_INITIALIZED)
   end function timestate_is_ready

   !> \brief Check if a date is a global holiday
   logical function is_global_holiday(month, day)
      integer, intent(in) :: month, day
      is_global_holiday = ( (month==1 .and. day==1) .or. (month==12 .and. day==25) )
   end function is_global_holiday

   !> \brief Check if a date is a U.S. holiday
   logical function is_us_holiday(month, day)
      integer, intent(in) :: month, day
      is_us_holiday = ( (month==7 .and. day==4) .or. (month==11 .and. day>=22 .and. day<=28) )
   end function is_us_holiday

   !> \brief Get current simulation time as ISO 8601 string (YYYY-MM-DDTHH:MM:SS)
   pure function get_time_iso8601(this) result(timestr)
      class(TimeStateType), intent(in) :: this
      character(len=25) :: timestr
      write(timestr, '(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2)') &
         this%year, this%month, this%day, this%hour, this%minute, this%second
   end function get_time_iso8601

   !> \brief Get current simulation time as a human-readable string (e.g., 2025-06-25 14:30:00)
   pure function get_time_human(this) result(timestr)
      class(TimeStateType), intent(in) :: this
      character(len=25) :: timestr
      write(timestr, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
         this%year, this%month, this%day, this%hour, this%minute, this%second
   end function get_time_human

   !> \brief Get current simulation time as a compact string (YYYYMMDD_HHMMSS)
   pure function get_time_compact(this) result(timestr)
      class(TimeStateType), intent(in) :: this
      character(len=16) :: timestr
      write(timestr, '(I4.4,I2.2,I2.2,"_",I2.2,I2.2,I2.2)') &
         this%year, this%month, this%day, this%hour, this%minute, this%second
   end function get_time_compact

   !> \brief Calculate local timezone offset (hours) from longitude (robust, clamped)
   pure integer function get_timezone_offset(this, lon) result(tz_offset)
      class(TimeStateType), intent(in) :: this
      real, intent(in) :: lon
      ! Truncate toward zero, clamp to [-12, 14] (real-world timezones)
      tz_offset = int(lon / 15.0)
      tz_offset = max(-12, min(14, tz_offset))
   end function get_timezone_offset

   !> \brief Calculate derived fields (Julian date, DOY) from date components
   !! \param[inout] this       TimeStateType object
   !! \param[inout] error_mgr  Error manager
   !! \param[out]   rc         Return code
   subroutine calculate_derived_fields(this, error_mgr, rc)
      use error_mod, only: ERROR_INVALID_INPUT

      class(TimeStateType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      integer :: a, y, m, jdn
      character(len=256) :: thisLoc

      thisLoc = 'calculate_derived_fields (in core/TimeState_Mod.F90)'
      rc = CC_SUCCESS

      ! Calculate Julian Date Number using standard algorithm
      if (this%month > 2) then
         a = 0
         y = this%year
         m = this%month
      else
         a = 1
         y = this%year - 1
         m = this%month + 12
      endif

      ! Julian Day Number (integer part)
      jdn = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + this%day - 1524 - a

      ! Julian Date (with fractional day)
      this%julian_date = real(jdn) + (this%hour + this%minute/60.0 + this%second/3600.0) / 24.0

      ! Calculate day of year
      this%doy = calculate_day_of_year(this%year, this%month, this%day)

   end subroutine calculate_derived_fields

   !> \brief Get number of days in a given month, accounting for leap years
   !! \param[in] month  Month (1-12)
   !! \param[in] year   Year
   !! \return           Number of days in the month
   pure integer function get_days_in_month(month, year) result(days)
      integer, intent(in) :: month, year

      integer, parameter :: days_per_month(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

      days = days_per_month(month)

      ! Adjust February for leap years
      if (month == 2 .and. is_leap_year(year)) then
         days = 29
      endif

   end function get_days_in_month

   !> \brief Check if a year is a leap year
   !! \param[in] year  Year to check
   !! \return          True if leap year, false otherwise
   pure logical function is_leap_year(year) result(is_leap)
      integer, intent(in) :: year

      is_leap = (mod(year, 4) == 0 .and. mod(year, 100) /= 0) .or. (mod(year, 400) == 0)

   end function is_leap_year

   !> \brief Calculate day of year from month and day
   !! \param[in] year   Year
   !! \param[in] month  Month (1-12)
   !! \param[in] day    Day of month
   !! \return           Day of year (1-366)
   pure integer function calculate_day_of_year(year, month, day) result(doy)
      integer, intent(in) :: year, month, day

      integer :: i
      integer, parameter :: days_per_month(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

      doy = day

      ! Add days from previous months
      do i = 1, month - 1
         doy = doy + days_per_month(i)
         ! Add extra day for February in leap years
         if (i == 2 .and. is_leap_year(year)) then
            doy = doy + 1
         endif
      enddo

   end function calculate_day_of_year

end module TimeState_Mod
