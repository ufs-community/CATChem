

# File TimeState\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**TimeState\_Mod.F90**](_time_state___mod_8_f90.md)

[Go to the documentation of this file](_time_state___mod_8_f90.md)


```Fortran

module timestate_mod
   use statemanager_mod, only: state_status_uninitialized, state_status_initialized
   use error_mod, only: errormanagertype, cc_success, cc_failure
   use constants, only: pi, pi_180
   implicit none
   private
   public :: timestatetype, is_global_holiday, is_us_holiday

   type :: timestatetype
      integer :: year = 2000
      integer :: month = 1
      integer :: day = 1
      integer :: hour = 0
      integer :: minute = 0
      integer :: second = 0
      real    :: timestep = 3600.0
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
   end type timestatetype

contains

   real function get_sza(this, lat, lon) result(sza)
      class(TimeStateType), intent(in) :: this
      real, intent(in) :: lat, lon
      ! Accurate solar zenith angle calculation
      ! Inputs: lat, lon in degrees; time from this%hour, this%minute, this%second; day of year from this%doy
      real :: lat_rad, lon_rad, decl_rad, ha_rad
      real :: decl, eqtime, time_offset, tst, ha
      real :: cos_sza_val
      real :: fractional_hour, gamma

      ! Convert latitude and longitude to radians
      lat_rad = lat * pi_180
      lon_rad = lon * pi_180

      ! Calculate fractional hour of the day (UTC)
      fractional_hour = real(this%hour) + real(this%minute)/60.0 + real(this%second)/3600.0

      ! Calculate day angle (in radians)
      gamma = 2.0 * pi * (real(this%doy) - 1.0) / 365.0

      ! Solar declination (in degrees, then radians)
      decl = 23.44 * sin(2.0 * pi * (real(this%doy) - 81.0) / 365.0)
      decl_rad = decl * pi_180

      ! Equation of time (in minutes)
      eqtime = 229.18 * (0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) \
      - 0.014615 * cos(2.0*gamma) - 0.040849 * sin(2.0*gamma))

      ! Time offset (in minutes)
      time_offset = eqtime + 4.0 * lon

      ! True solar time (in minutes)
      tst = fractional_hour * 60.0 + time_offset

      ! Hour angle (in degrees, then radians)
      ha = (tst / 4.0) - 180.0
      ha_rad = ha * pi_180

      ! Solar zenith angle calculation
      cos_sza_val = sin(lat_rad) * sin(decl_rad) + cos(lat_rad) * cos(decl_rad) * cos(ha_rad)
      cos_sza_val = max(-1.0, min(1.0, cos_sza_val)) ! Clamp for safety
      sza = acos(cos_sza_val) / pi_180
      sza = min(max(sza, 0.0), 90.0) ! Clamp to [0, 90] degrees
   end function get_sza

   real function get_cos_sza(this, lat, lon) result(cos_sza)
      class(TimeStateType), intent(in) :: this
      real, intent(in) :: lat, lon
      cos_sza = cos(this%get_sza(lat, lon) * pi_180)
   end function get_cos_sza

   real function get_timestep(this) result(dt)
      class(TimeStateType), intent(in) :: this
      dt = this%timestep
   end function get_timestep

   subroutine get_current_date(this, year, month, day)
      class(TimeStateType), intent(in) :: this
      integer, intent(out) :: year, month, day
      year = this%year
      month = this%month
      day = this%day
   end subroutine get_current_date

   real function get_julian_date(this) result(jd)
      class(TimeStateType), intent(in) :: this
      jd = this%julian_date
   end function get_julian_date

   integer function get_doy(this) result(doy)
      class(TimeStateType), intent(in) :: this
      doy = this%doy
   end function get_doy

   subroutine timestate_init(this, error_mgr, rc)
      class(TimeStateType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Suppress unused parameter warning
      if (associated(error_mgr)) continue
      rc = cc_success
   end subroutine timestate_init

   subroutine timestate_validate(this, error_mgr, rc)
      class(TimeStateType), intent(in) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Suppress unused parameter warning
      if (associated(error_mgr)) continue
      rc = cc_success
   end subroutine timestate_validate

   subroutine timestate_cleanup(this, error_mgr, rc)
      class(TimeStateType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Suppress unused parameter warning
      if (associated(error_mgr)) continue
      rc = cc_success
   end subroutine timestate_cleanup

   subroutine timestate_reset(this, error_mgr, rc)
      class(TimeStateType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      ! Suppress unused parameter warning
      if (associated(error_mgr)) continue
      rc = cc_success
   end subroutine timestate_reset

   function timestate_get_status(this) result(status)
      class(TimeStateType), intent(in) :: this
      integer :: status
      status = state_status_initialized
   end function timestate_get_status

   function timestate_get_memory_usage(this) result(memory_bytes)
      class(TimeStateType), intent(in) :: this
      integer(8) :: memory_bytes
      memory_bytes = 0_8
   end function timestate_get_memory_usage

   subroutine timestate_print_info(this, unit)
      class(TimeStateType), intent(in) :: this
      integer, optional, intent(in) :: unit
      if (present(unit)) then
         write(unit,*) 'TimeStateType: ', this%year, this%month, this%day, this%hour, this%minute, this%second
      else
         print *, 'TimeStateType: ', this%year, this%month, this%day, this%hour, this%minute, this%second
      end if
   end subroutine timestate_print_info

   function timestate_is_ready(this) result(ready)
      class(TimeStateType), intent(in) :: this
      logical :: ready
      ready = .true.
   end function timestate_is_ready

   logical function is_global_holiday(month, day)
      integer, intent(in) :: month, day
      is_global_holiday = ( (month==1 .and. day==1) .or. (month==12 .and. day==25) )
   end function is_global_holiday

   logical function is_us_holiday(month, day)
      integer, intent(in) :: month, day
      is_us_holiday = ( (month==7 .and. day==4) .or. (month==11 .and. day>=22 .and. day<=28) )
   end function is_us_holiday

   pure function get_time_iso8601(this) result(timestr)
      class(TimeStateType), intent(in) :: this
      character(len=25) :: timestr
      write(timestr, '(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2)') &
         this%year, this%month, this%day, this%hour, this%minute, this%second
   end function get_time_iso8601

   pure function get_time_human(this) result(timestr)
      class(TimeStateType), intent(in) :: this
      character(len=25) :: timestr
      write(timestr, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
         this%year, this%month, this%day, this%hour, this%minute, this%second
   end function get_time_human

   pure function get_time_compact(this) result(timestr)
      class(TimeStateType), intent(in) :: this
      character(len=16) :: timestr
      write(timestr, '(I4.4,I2.2,I2.2,"_",I2.2,I2.2,I2.2)') &
         this%year, this%month, this%day, this%hour, this%minute, this%second
   end function get_time_compact

   pure integer function get_timezone_offset(this, lon) result(tz_offset)
      class(TimeStateType), intent(in) :: this
      real, intent(in) :: lon
      ! Truncate toward zero, clamp to [-12, 14] (real-world timezones)
      tz_offset = int(lon / 15.0)
      tz_offset = max(-12, min(14, tz_offset))
   end function get_timezone_offset

end module timestate_mod
```


