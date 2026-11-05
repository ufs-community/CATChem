

# File utilities\_mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**utilities\_mod.F90**](utilities__mod_8_f90.md)

[Go to the documentation of this file](utilities__mod_8_f90.md)


```Fortran

module utilities_mod
   use precision_mod
   use constants
   use error_mod, only: cc_success, cc_failure, error_numerical_instability, error_invalid_input, &
      errormanagertype

   implicit none
   private

   ! Public utility functions
   public :: validate_atmospheric_constants
   public :: convert_pressure_units
   public :: convert_temperature_units
   public :: calculate_air_density
   public :: calculate_scale_height
   public :: is_valid_temperature
   public :: is_valid_pressure
   public :: check_array_bounds
   public :: safe_divide
   public :: calculate_geopotential_height

contains

   subroutine validate_atmospheric_constants(rc, error_mgr)
      implicit none
      integer, intent(out) :: rc
      type(ErrorManagerType), intent(inout), optional :: error_mgr

      rc = cc_success

      if (present(error_mgr)) then
         call error_mgr%push_context('validate_atmospheric_constants', 'checking physical constants')
      endif

      ! Check fundamental constants
      if (g0 <= 0.0_fp) then
         if (present(error_mgr)) then
            call error_mgr%report_error(error_numerical_instability, &
               'Invalid gravitational acceleration', rc)
         else
            rc = cc_failure
         endif
         return
      endif

      if (rd <= 0.0_fp) then
         if (present(error_mgr)) then
            call error_mgr%report_error(error_numerical_instability, &
               'Invalid dry air gas constant', rc)
         else
            rc = cc_failure
         endif
         return
      endif

      ! Check derived relationships
      if (abs(cp - cv - rd) > 1.0e-6_fp) then
         if (present(error_mgr)) then
            call error_mgr%report_error(error_numerical_instability, &
               ≠'Inconsistent thermodynamic constants: Cp - Cv  Rd', rc)
         else
            rc = cc_failure
         endif
         return
      endif

      if (present(error_mgr)) then
         call error_mgr%pop_context()
      endif

   end subroutine validate_atmospheric_constants

   function convert_pressure_units(pressure_in, unit_in, unit_out, rc) result(pressure_out)
      implicit none
      real(fp), intent(in) :: pressure_in
      character(len=*), intent(in) :: unit_in, unit_out
      integer, intent(out) :: rc
      real(fp) :: pressure_out

      real(fp) :: pressure_pa  ! Intermediate value in Pascals

      rc = cc_success
      pressure_out = 0.0_fp

      ! Convert input to Pascals
      select case (trim(unit_in))
       case ('Pa')
         pressure_pa = pressure_in
       case ('hPa', 'mbar')
         pressure_pa = pressure_in * 100.0_fp
       case ('atm')
         pressure_pa = pressure_in * atm
       case ('mmHg', 'torr')
         pressure_pa = pressure_in * 133.322_fp
       case default
         rc = error_invalid_input
         return
      end select

      ! Convert from Pascals to output unit
      select case (trim(unit_out))
       case ('Pa')
         pressure_out = pressure_pa
       case ('hPa', 'mbar')
         pressure_out = pressure_pa / 100.0_fp
       case ('atm')
         pressure_out = pressure_pa / atm
       case ('mmHg', 'torr')
         pressure_out = pressure_pa / 133.322_fp
       case default
         rc = error_invalid_input
         return
      end select

   end function convert_pressure_units

   function convert_temperature_units(temp_in, unit_in, unit_out, rc) result(temp_out)
      implicit none
      real(fp), intent(in) :: temp_in
      character(len=*), intent(in) :: unit_in, unit_out
      integer, intent(out) :: rc
      real(fp) :: temp_out

      real(fp) :: temp_k  ! Intermediate value in Kelvin

      rc = cc_success
      temp_out = 0.0_fp

      ! Convert input to Kelvin
      select case (trim(unit_in))
       case ('K')
         temp_k = temp_in
       case ('C')
         temp_k = temp_in + 273.15_fp
       case ('F')
         temp_k = (temp_in - 32.0_fp) * 5.0_fp/9.0_fp + 273.15_fp
       case default
         rc = error_invalid_input
         return
      end select

      ! Convert from Kelvin to output unit
      select case (trim(unit_out))
       case ('K')
         temp_out = temp_k
       case ('C')
         temp_out = temp_k - 273.15_fp
       case ('F')
         temp_out = (temp_k - 273.15_fp) * 9.0_fp/5.0_fp + 32.0_fp
       case default
         rc = error_invalid_input
         return
      end select

   end function convert_temperature_units

   function calculate_air_density(pressure, temperature, rc) result(density)
      implicit none
      real(fp), intent(in) :: pressure, temperature
      integer, intent(out) :: rc
      real(fp) :: density

      rc = cc_success

      if (temperature <= 0.0_fp) then
         rc = error_invalid_input
         density = 0.0_fp
         return
      endif

      if (pressure < 0.0_fp) then
         rc = error_invalid_input
         density = 0.0_fp
         return
      endif

      density = pressure / (rd * temperature)

   end function calculate_air_density

   function calculate_scale_height(temperature, rc, z) result(scale_height)
      use constants, only: g0, rd, re
      implicit none
      real(fp), intent(in) :: temperature
      integer, intent(out) :: rc
      real(fp), intent(in), optional :: z
      real(fp) :: scale_height
      real(fp) :: g_local

      rc = cc_success

      if (temperature <= 0.0_fp) then
         rc = error_invalid_input
         scale_height = 0.0_fp
         return
      endif

      if (present(z)) then
         g_local = g0 * (re / (re + z))
      else
         g_local = g0
      endif

      scale_height = rd * temperature / g_local

   end function calculate_scale_height

   function is_valid_temperature(temperature) result(is_valid)
      implicit none
      real(fp), intent(in) :: temperature
      logical :: is_valid

      ! Valid range: 100K to 400K (typical atmospheric range)
      is_valid = (temperature >= 100.0_fp .and. temperature <= 400.0_fp)

   end function is_valid_temperature

   function is_valid_pressure(pressure) result(is_valid)
      implicit none
      real(fp), intent(in) :: pressure
      logical :: is_valid

      ! Valid range: 1 Pa to 110000 Pa (typical atmospheric range)
      is_valid = (pressure >= 1.0_fp .and. pressure <= 110000.0_fp)

   end function is_valid_pressure

   subroutine check_array_bounds(index, array_size, array_name, rc)
      implicit none
      integer, intent(in) :: index, array_size
      character(len=*), intent(in) :: array_name
      integer, intent(out) :: rc

      rc = cc_success

      if (index < 1 .or. index > array_size) then
         rc = error_invalid_input
      endif

   end subroutine check_array_bounds

   function safe_divide(numerator, denominator, rc) result(quotient)
      implicit none
      real(fp), intent(in) :: numerator, denominator
      integer, intent(out) :: rc
      real(fp) :: quotient

      rc = cc_success

      if (abs(denominator) < epsilon(denominator) * 1000.0_fp) then
         rc = error_numerical_instability
         quotient = 0.0_fp
      else
         quotient = numerator / denominator
      endif

   end function safe_divide

   function calculate_geopotential_height(p1, p2, Tv_mean, rc) result(z)
      use constants, only: g0, rd
      implicit none
      real(fp), intent(in) :: p1, p2, Tv_mean
      integer, intent(out) :: rc
      real(fp) :: z

      rc = cc_success
      z = 0.0_fp

      if (p1 <= 0.0_fp .or. p2 <= 0.0_fp .or. tv_mean <= 0.0_fp) then
         rc = error_invalid_input
         return
      endif

      if (p1 <= p2) then
         rc = error_invalid_input
         return
      endif

      z = (rd * tv_mean / g0) * log(p1 / p2)

   end function calculate_geopotential_height

end module utilities_mod
```


