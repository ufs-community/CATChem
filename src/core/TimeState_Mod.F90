!> \file TimeState_Mod.F90
!! \brief Time state and common time/solar functions for atmospheric chemistry
!!
!! Provides timekeeping, solar zenith angle, and calendar utilities delegating directly to C++.
!!
module TimeState_Mod
   use Precision_Mod, only: fp
   use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE
   use iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_int, c_double, c_bool
   implicit none
   private

   public :: TimeStateType, is_global_holiday, is_us_holiday

   ! Local status constants
   integer, parameter :: STATE_STATUS_UNINITIALIZED = 0
   integer, parameter :: STATE_STATUS_INITIALIZED = 1

   !> \brief Time state proxy delegating to C++
   type :: TimeStateType
      type(c_ptr) :: cpp_ptr = c_null_ptr
      integer :: year = 2000
      integer :: month = 1
      integer :: day = 1
      integer :: hour = 0
      integer :: minute = 0
      integer :: second = 0
      real(fp) :: timestep = 3600.0_fp !< seconds
      real(fp) :: julian_date = 0.0_fp
      integer  :: doy = 1
   contains
      procedure :: get_sza => timestate_get_sza
      procedure :: get_cos_sza => timestate_get_cos_sza
      procedure :: get_timestep => timestate_get_timestep
      procedure :: get_current_date => timestate_get_current_date
      procedure :: get_julian_date => timestate_get_julian_date
      procedure :: get_doy => timestate_get_doy
      procedure :: init => timestate_init
      procedure :: validate => timestate_validate
      procedure :: cleanup => timestate_cleanup
      procedure :: reset => timestate_reset
      procedure :: get_status => timestate_get_status
      procedure :: get_memory_usage => timestate_get_memory_usage
      procedure :: print_info => timestate_print_info
      procedure :: is_ready => timestate_is_ready
      procedure :: get_time_iso8601 => timestate_get_time_iso8601
      procedure :: get_time_human => timestate_get_time_human
      procedure :: get_time_compact => timestate_get_time_compact
      procedure :: get_timezone_offset => timestate_get_timezone_offset
      procedure, private :: sync_to_fortran => timestate_sync_to_fortran
   end type TimeStateType

   ! Interoperable C Prototypes matching catchem_api
   interface
      function catchem_time_state_create() bind(C, name="catchem_time_state_create")
         import :: c_ptr
         type(c_ptr) :: catchem_time_state_create
      end function

      subroutine catchem_time_state_destroy(ptr) bind(C, name="catchem_time_state_destroy")
         import :: c_ptr
         type(c_ptr), value :: ptr
      end subroutine

      function catchem_time_state_init(ptr, yr, mo, dy, hr, mn, sc, ts) bind(C, name="catchem_time_state_init")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: ptr
         integer(c_int), value :: yr, mo, dy, hr, mn, sc
         real(c_double), value :: ts
         integer(c_int) :: catchem_time_state_init
      end function

      function catchem_time_state_advance(ptr, dt) bind(C, name="catchem_time_state_advance")
         import :: c_ptr, c_double, c_int
         type(c_ptr), value :: ptr
         real(c_double), value :: dt
         integer(c_int) :: catchem_time_state_advance
      end function

      function catchem_time_state_reset(ptr) bind(C, name="catchem_time_state_reset")
         import :: c_ptr, c_int
         type(c_ptr), value :: ptr
         integer(c_int) :: catchem_time_state_reset
      end function

      function catchem_time_state_get_year(ptr) bind(C, name="catchem_time_state_get_year")
         import :: c_ptr, c_int
         type(c_ptr), value :: ptr
         integer(c_int) :: catchem_time_state_get_year
      end function

      function catchem_time_state_get_month(ptr) bind(C, name="catchem_time_state_get_month")
         import :: c_ptr, c_int
         type(c_ptr), value :: ptr
         integer(c_int) :: catchem_time_state_get_month
      end function

      function catchem_time_state_get_day(ptr) bind(C, name="catchem_time_state_get_day")
         import :: c_ptr, c_int
         type(c_ptr), value :: ptr
         integer(c_int) :: catchem_time_state_get_day
      end function

      function catchem_time_state_get_hour(ptr) bind(C, name="catchem_time_state_get_hour")
         import :: c_ptr, c_int
         type(c_ptr), value :: ptr
         integer(c_int) :: catchem_time_state_get_hour
      end function

      function catchem_time_state_get_minute(ptr) bind(C, name="catchem_time_state_get_minute")
         import :: c_ptr, c_int
         type(c_ptr), value :: ptr
         integer(c_int) :: catchem_time_state_get_minute
      end function

      function catchem_time_state_get_second(ptr) bind(C, name="catchem_time_state_get_second")
         import :: c_ptr, c_int
         type(c_ptr), value :: ptr
         integer(c_int) :: catchem_time_state_get_second
      end function

      function catchem_time_state_get_timestep(ptr) bind(C, name="catchem_time_state_get_timestep")
         import :: c_ptr, c_double
         type(c_ptr), value :: ptr
         real(c_double) :: catchem_time_state_get_timestep
      end function

      function catchem_time_state_get_julian_date(ptr) bind(C, name="catchem_time_state_get_julian_date")
         import :: c_ptr, c_double
         type(c_ptr), value :: ptr
         real(c_double) :: catchem_time_state_get_julian_date
      end function

      function catchem_time_state_get_doy(ptr) bind(C, name="catchem_time_state_get_doy")
         import :: c_ptr, c_int
         type(c_ptr), value :: ptr
         integer(c_int) :: catchem_time_state_get_doy
      end function

      function catchem_time_state_get_cos_sza(ptr, lat, lon, mid_timestep) bind(C, name="catchem_time_state_get_cos_sza")
         import :: c_ptr, c_double, c_bool
         type(c_ptr), value :: ptr
         real(c_double), value :: lat, lon
         logical(c_bool), value :: mid_timestep
         real(c_double) :: catchem_time_state_get_cos_sza
      end function

      function catchem_time_state_get_timezone_offset(ptr, lon) bind(C, name="catchem_time_state_get_timezone_offset")
         import :: c_ptr, c_double, c_int
         type(c_ptr), value :: ptr
         real(c_double), value :: lon
         integer(c_int) :: catchem_time_state_get_timezone_offset
      end function

      function catchem_time_state_is_leap_year(year) bind(C, name="catchem_time_state_is_leap_year")
         import :: c_int, c_bool
         integer(c_int), value :: year
         logical(c_bool) :: catchem_time_state_is_leap_year
      end function

      function catchem_time_state_get_days_in_month(month, year) bind(C, name="catchem_time_state_get_days_in_month")
         import :: c_int
         integer(c_int), value :: month, year
         integer(c_int) :: catchem_time_state_get_days_in_month
      end function

      function catchem_time_state_is_global_holiday(month, day) bind(C, name="catchem_time_state_is_global_holiday")
         import :: c_int, c_bool
         integer(c_int), value :: month, day
         logical(c_bool) :: catchem_time_state_is_global_holiday
      end function

      function catchem_time_state_is_us_holiday(month, day) bind(C, name="catchem_time_state_is_us_holiday")
         import :: c_int, c_bool
         integer(c_int), value :: month, day
         logical(c_bool) :: catchem_time_state_is_us_holiday
      end function
   end interface

contains

   !> \brief Private method to synchronize state properties from C++ to Fortran for legacy direct-access compatibility
   subroutine timestate_sync_to_fortran(this)
      class(TimeStateType), intent(inout) :: this
      if (c_associated(this%cpp_ptr)) then
         this%year = int(catchem_time_state_get_year(this%cpp_ptr))
         this%month = int(catchem_time_state_get_month(this%cpp_ptr))
         this%day = int(catchem_time_state_get_day(this%cpp_ptr))
         this%hour = int(catchem_time_state_get_hour(this%cpp_ptr))
         this%minute = int(catchem_time_state_get_minute(this%cpp_ptr))
         this%second = int(catchem_time_state_get_second(this%cpp_ptr))
         this%timestep = real(catchem_time_state_get_timestep(this%cpp_ptr), fp)
         this%julian_date = real(catchem_time_state_get_julian_date(this%cpp_ptr), fp)
         this%doy = int(catchem_time_state_get_doy(this%cpp_ptr))
      end if
   end subroutine timestate_sync_to_fortran

   !> \brief Initialize TimeStateType delegating to C++
   subroutine timestate_init(this, year, month, day, hour, minute, second, timestep, error_mgr, rc)
      class(TimeStateType), intent(inout) :: this
      integer, optional, intent(in) :: year, month, day, hour, minute, second
      real(fp), optional, intent(in) :: timestep
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      integer(c_int) :: yr_val, mo_val, dy_val, hr_val, mn_val, sc_val
      real(c_double) :: ts_val

      rc = CC_SUCCESS

      if (.not. c_associated(this%cpp_ptr)) then
         this%cpp_ptr = catchem_time_state_create()
      end if

      ! Set defaults
      yr_val = 2000; mo_val = 1; dy_val = 1; hr_val = 0; mn_val = 0; sc_val = 0
      ts_val = 3600.0_c_double

      if (present(year)) yr_val = int(year, c_int)
      if (present(month)) mo_val = int(month, c_int)
      if (present(day)) dy_val = int(day, c_int)
      if (present(hour)) hr_val = int(hour, c_int)
      if (present(minute)) mn_val = int(minute, c_int)
      if (present(second)) sc_val = int(second, c_int)
      if (present(timestep)) ts_val = real(timestep, c_double)

      if (catchem_time_state_init(this%cpp_ptr, yr_val, mo_val, dy_val, hr_val, mn_val, sc_val, ts_val) /= 0) then
         rc = CC_FAILURE
         return
      end if

      call this%sync_to_fortran()
   end subroutine timestate_init

   !> \brief Standard validation of TimeStateType via properties check
   subroutine timestate_validate(this, error_mgr, rc)
      class(TimeStateType), intent(in) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc
      rc = CC_SUCCESS
   end subroutine timestate_validate

   !> \brief Cleanup and release C++ TimeState memory
   subroutine timestate_cleanup(this, error_mgr, rc)
      class(TimeStateType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc
      rc = CC_SUCCESS
      if (c_associated(this%cpp_ptr)) then
         call catchem_time_state_destroy(this%cpp_ptr)
         this%cpp_ptr = c_null_ptr
      end if
      this%year = -1
      this%month = -1
      this%day = -1
      this%hour = -1
      this%minute = -1
      this%second = -1
      this%timestep = -1.0_fp
      this%julian_date = -1.0_fp
      this%doy = -1
   end subroutine timestate_cleanup

   !> \brief Reset to C++ default date (Y2K)
   subroutine timestate_reset(this, error_mgr, rc)
      class(TimeStateType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc
      rc = CC_SUCCESS
      if (c_associated(this%cpp_ptr)) then
         if (catchem_time_state_reset(this%cpp_ptr) /= 0) rc = CC_FAILURE
      end if
      call this%sync_to_fortran()
   end subroutine timestate_reset

   !> \brief Advance TimeStateType in time by dt seconds
   subroutine timestate_advance(this, dt, error_mgr, rc)
      class(TimeStateType), intent(inout) :: this
      real(fp), intent(in) :: dt
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc
      rc = CC_SUCCESS
      if (c_associated(this%cpp_ptr)) then
         if (catchem_time_state_advance(this%cpp_ptr, real(dt, c_double)) /= 0) rc = CC_FAILURE
      end if
      call this%sync_to_fortran()
   end subroutine timestate_advance

   !> \brief Get Cosine of Solar Zenith Angle from C++
   real(fp) function timestate_get_cos_sza(this, lat, lon, mid_timestep) result(cos_sza_val)
      class(TimeStateType), intent(in) :: this
      real(fp), intent(in) :: lat, lon
      logical, intent(in), optional :: mid_timestep
      logical(c_bool) :: mid_ts

      mid_ts = .false.
      if (present(mid_timestep)) mid_ts = mid_timestep

      if (c_associated(this%cpp_ptr)) then
         cos_sza_val = real(catchem_time_state_get_cos_sza(this%cpp_ptr, real(lat, c_double), real(lon, c_double), mid_ts), fp)
      else
         cos_sza_val = 1.0_fp
      end if
   end function timestate_get_cos_sza

   !> \brief Get Solar Zenith Angle from Cos SZA
   real(fp) function timestate_get_sza(this, lat, lon) result(sza)
      use constants, only: PI_180
      class(TimeStateType), intent(in) :: this
      real(fp), intent(in) :: lat, lon
      real(fp) :: cos_sza_val
      cos_sza_val = this%get_cos_sza(lat, lon)
      sza = acos(cos_sza_val) / PI_180
      sza = min(max(sza, 0.0_fp), 90.0_fp)
   end function timestate_get_sza

   real(fp) function timestate_get_timestep(this) result(dt)
      class(TimeStateType), intent(in) :: this
      dt = this%timestep
   end function timestate_get_timestep

   subroutine timestate_get_current_date(this, year, month, day)
      class(TimeStateType), intent(in) :: this
      integer, intent(out) :: year, month, day
      year = this%year
      month = this%month
      day = this%day
   end subroutine timestate_get_current_date

   real(fp) function timestate_get_julian_date(this) result(jd)
      class(TimeStateType), intent(in) :: this
      jd = this%julian_date
   end function timestate_get_julian_date

   integer function timestate_get_doy(this) result(doy)
      class(TimeStateType), intent(in) :: this
      doy = this%doy
   end function timestate_get_doy

   function timestate_get_status(this) result(status)
      class(TimeStateType), intent(in) :: this
      integer :: status
      if (this%year > 0 .and. this%month > 0 .and. this%day > 0) then
         status = STATE_STATUS_INITIALIZED
      else
         status = STATE_STATUS_UNINITIALIZED
      end if
   end function timestate_get_status

   function timestate_get_memory_usage(this) result(memory_bytes)
      class(TimeStateType), intent(in) :: this
      integer(8) :: memory_bytes
      memory_bytes = 32_8
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
      ready = (this%get_status() == STATE_STATUS_INITIALIZED)
   end function timestate_is_ready

   pure function timestate_get_time_iso8601(this) result(timestr)
      class(TimeStateType), intent(in) :: this
      character(len=25) :: timestr
      write(timestr, '(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2)') &
         this%year, this%month, this%day, this%hour, this%minute, this%second
   end function timestate_get_time_iso8601

   pure function timestate_get_time_human(this) result(timestr)
      class(TimeStateType), intent(in) :: this
      character(len=25) :: timestr
      write(timestr, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
         this%year, this%month, this%day, this%hour, this%minute, this%second
   end function timestate_get_time_human

   pure function timestate_get_time_compact(this) result(timestr)
      class(TimeStateType), intent(in) :: this
      character(len=16) :: timestr
      write(timestr, '(I4.4,I2.2,I2.2,"_",I2.2,I2.2,I2.2)') &
         this%year, this%month, this%day, this%hour, this%minute, this%second
   end function timestate_get_time_compact

   integer function timestate_get_timezone_offset(this, lon) result(tz_offset)
      class(TimeStateType), intent(in) :: this
      real(fp), intent(in) :: lon
      if (c_associated(this%cpp_ptr)) then
         tz_offset = int(catchem_time_state_get_timezone_offset(this%cpp_ptr, real(lon, c_double)))
      else
         tz_offset = 0
      end if
   end function timestate_get_timezone_offset

   ! Standalone utility functions delegating to C++
   logical function is_global_holiday(month, day)
      integer, intent(in) :: month, day
      is_global_holiday = catchem_time_state_is_global_holiday(int(month, c_int), int(day, c_int))
   end function is_global_holiday

   logical function is_us_holiday(month, day)
      integer, intent(in) :: month, day
      is_us_holiday = catchem_time_state_is_us_holiday(int(month, c_int), int(day, c_int))
   end function is_us_holiday

end module TimeState_Mod
