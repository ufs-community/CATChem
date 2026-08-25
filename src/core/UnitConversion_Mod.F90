!> \file UnitConversion_Mod.F90
!! \brief Comprehensive unit conversion utilities for atmospheric chemistry
!! \ingroup core_modules
!!
!! This module delegates all unit conversion arithmetic and chemical constants
!! to the optimized C++20 core.
!!
module UnitConversion_Mod
   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE
   use Constants, only: ATM
   use iso_c_binding, only: c_double, c_int, c_char, c_null_char, c_bool

   implicit none
   private

   public :: UnitConverterType
   public :: convert_concentration
   public :: convert_pressure
   public :: convert_temperature
   public :: convert_flux
   public :: convert_rate_constant
   public :: convert_mass_units
   public :: convert_imperial_length
   public :: convert_imperial_area
   public :: convert_imperial_volume
   public :: convert_imperial_speed
   public :: convert_imperial_force
   public :: convert_imperial_pressure
   public :: convert_imperial_temperature
   public :: convert_imperial_mass
   public :: convert_imperial_energy
   public :: calculate_air_density
   public :: calculate_molecular_weight
   public :: convert_process_concentration_units
   public :: convert_process_flux_units

   ! Standard conditions using constants from Constants module
   real(fp), parameter :: STANDARD_TEMP = 273.15_fp  !< Standard temperature [K]
   real(fp), parameter :: STANDARD_PRESS = ATM       !< Standard pressure [Pa] from Constants

   !> \brief Unit converter type for managing conversions delegating to C++
   type :: UnitConverterType
      real(fp) :: temperature = STANDARD_TEMP    !< Reference temperature [K]
      real(fp) :: pressure = STANDARD_PRESS      !< Reference pressure [Pa]
      real(fp) :: air_density = 1.225_fp         !< Air density [kg/m³]
      logical :: use_standard_conditions = .true. !< Use STP conditions

   contains
      procedure :: init => converter_init
      procedure :: set_conditions => converter_set_conditions
      procedure :: get_air_density => converter_get_air_density
      procedure :: calculate_number_density => converter_calculate_number_density

      ! Concentration conversions
      procedure :: ppbv_to_ugm3 => converter_ppbv_to_ugm3
      procedure :: ugm3_to_ppbv => converter_ugm3_to_ppbv
      procedure :: molcm3_to_ppbv => converter_molcm3_to_ppbv
      procedure :: ppbv_to_molcm3 => converter_ppbv_to_molcm3
      procedure :: ppmv_to_mgm3 => converter_ppmv_to_mgm3
      procedure :: mgm3_to_ppmv => converter_mgm3_to_ppmv

      ! Column integrals
      procedure :: calculate_column_mass => converter_calculate_column_mass
      procedure :: calculate_dobson_units => converter_calculate_dobson_units

      ! Flux conversions
      procedure :: molcm2s_to_kgm2s => converter_molcm2s_to_kgm2s
      procedure :: kgm2s_to_molcm2s => converter_kgm2s_to_molcm2s

      ! Rate constant conversions
      procedure :: convert_rate_units => converter_convert_rate_units
   end type UnitConverterType

   ! Interoperable C Prototypes matching catchem_api
   interface
      function catchem_convert_concentration(val, from_u, to_u, mw, temp, press, rc) bind(C, name="catchem_convert_concentration")
         import :: c_double, c_char, c_int
         real(c_double), value :: val
         character(c_char), intent(in) :: from_u(*), to_u(*)
         real(c_double), value :: mw, temp, press
         integer(c_int), intent(out) :: rc
         real(c_double) :: catchem_convert_concentration
      end function

      function catchem_convert_pressure(val, from_u, to_u, rc) bind(C, name="catchem_convert_pressure")
         import :: c_double, c_char, c_int
         real(c_double), value :: val
         character(c_char), intent(in) :: from_u(*), to_u(*)
         integer(c_int), intent(out) :: rc
         real(c_double) :: catchem_convert_pressure
      end function

      function catchem_convert_temperature(val, from_u, to_u, rc) bind(C, name="catchem_convert_temperature")
         import :: c_double, c_char, c_int
         real(c_double), value :: val
         character(c_char), intent(in) :: from_u(*), to_u(*)
         integer(c_int), intent(out) :: rc
         real(c_double) :: catchem_convert_temperature
      end function

      function catchem_convert_flux(val, from_u, to_u, mw, rc) bind(C, name="catchem_convert_flux")
         import :: c_double, c_char, c_int
         real(c_double), value :: val
         character(c_char), intent(in) :: from_u(*), to_u(*)
         real(c_double), value :: mw
         integer(c_int), intent(out) :: rc
         real(c_double) :: catchem_convert_flux
      end function

      function catchem_convert_rate_constant(val, from_u, to_u, rc) bind(C, name="catchem_convert_rate_constant")
         import :: c_double, c_char, c_int
         real(c_double), value :: val
         character(c_char), intent(in) :: from_u(*), to_u(*)
         integer(c_int), intent(out) :: rc
         real(c_double) :: catchem_convert_rate_constant
      end function

      function catchem_convert_mass_units(val, from_u, to_u, rc) bind(C, name="catchem_convert_mass_units")
         import :: c_double, c_char, c_int
         real(c_double), value :: val
         character(c_char), intent(in) :: from_u(*), to_u(*)
         integer(c_int), intent(out) :: rc
         real(c_double) :: catchem_convert_mass_units
      end function

      function catchem_calculate_air_density(temp, press, humidity, use_humidity) bind(C, name="catchem_calculate_air_density")
         import :: c_double, c_bool
         real(c_double), value :: temp, press, humidity
         logical(c_bool), value :: use_humidity
         real(c_double) :: catchem_calculate_air_density
      end function

      function catchem_calculate_molecular_weight(formula) bind(C, name="catchem_calculate_molecular_weight")
         import :: c_double, c_char
         character(c_char), intent(in) :: formula(*)
         real(c_double) :: catchem_calculate_molecular_weight
      end function

      function catchem_convert_imperial(val, from_u, to_u, category, rc) bind(C, name="catchem_convert_imperial")
         import :: c_double, c_char, c_int
         real(c_double), value :: val
         character(c_char), intent(in) :: from_u(*), to_u(*), category(*)
         integer(c_int), intent(out) :: rc
         real(c_double) :: catchem_convert_imperial
      end function

      function catchem_convert_process_concentration_units(values, size, from_u, to_u, mw, temp, press) bind(C, name="catchem_convert_process_concentration_units")
         import :: fp, c_char, c_int
         real(fp), intent(inout) :: values(*)
         integer(c_int), value :: size
         character(c_char), intent(in) :: from_u(*), to_u(*)
         real(fp), value :: mw, temp, press
         integer(c_int) :: catchem_convert_process_concentration_units
      end function

      function catchem_convert_process_flux_units(values, size, from_u, to_u, mw) bind(C, name="catchem_convert_process_flux_units")
         import :: fp, c_char, c_int
         real(fp), intent(inout) :: values(*)
         integer(c_int), value :: size
         character(c_char), intent(in) :: from_u(*), to_u(*)
         real(fp), value :: mw
         integer(c_int) :: catchem_convert_process_flux_units
      end function
   end interface

contains

   !> \brief Convert concentration units between different systems
   subroutine convert_concentration(input_value, input_units, output_units, &
      molecular_weight, temperature, pressure, output_value, rc)
      real(fp), intent(in) :: input_value
      character(len=*), intent(in) :: input_units, output_units
      real(fp), intent(in) :: molecular_weight
      real(fp), intent(in) :: temperature
      real(fp), intent(in) :: pressure
      real(fp), intent(out) :: output_value
      integer, intent(out) :: rc

      integer(c_int) :: rc_c

      output_value = real(catchem_convert_concentration(real(input_value, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, &
         real(molecular_weight, c_double), real(temperature, c_double), real(pressure, c_double), rc_c), fp)
      rc = int(rc_c)
   end subroutine convert_concentration

   !> \brief Convert pressure units
   function convert_pressure(pressure_in, input_units, output_units, rc) result(pressure_out)
      real(fp), intent(in) :: pressure_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: pressure_out

      integer(c_int) :: rc_c

      pressure_out = real(catchem_convert_pressure(real(pressure_in, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, rc_c), fp)
      rc = int(rc_c)
   end function convert_pressure

   !> \brief Convert temperature units
   function convert_temperature(temp_in, input_units, output_units, rc) result(temp_out)
      real(fp), intent(in) :: temp_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: temp_out

      integer(c_int) :: rc_c

      temp_out = real(catchem_convert_temperature(real(temp_in, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, rc_c), fp)
      rc = int(rc_c)
   end function convert_temperature

   !> \brief Convert flux units
   function convert_flux(flux_in, input_units, output_units, molecular_weight, rc) result(flux_out)
      real(fp), intent(in) :: flux_in
      character(len=*), intent(in) :: input_units, output_units
      real(fp), intent(in) :: molecular_weight
      integer, intent(out) :: rc
      real(fp) :: flux_out

      integer(c_int) :: rc_c

      flux_out = real(catchem_convert_flux(real(flux_in, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, &
         real(molecular_weight, c_double), rc_c), fp)
      rc = int(rc_c)
   end function convert_flux

   !> \brief Convert rate constant units
   function convert_rate_constant(rate_in, input_units, output_units, rc) result(rate_out)
      real(fp), intent(in) :: rate_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: rate_out

      integer(c_int) :: rc_c

      rate_out = real(catchem_convert_rate_constant(real(rate_in, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, rc_c), fp)
      rc = int(rc_c)
   end function convert_rate_constant

   !> \brief Convert mass units
   function convert_mass_units(mass_in, input_units, output_units, rc) result(mass_out)
      real(fp), intent(in) :: mass_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: mass_out

      integer(c_int) :: rc_c

      mass_out = real(catchem_convert_mass_units(real(mass_in, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, rc_c), fp)
      rc = int(rc_c)
   end function convert_mass_units

   !> \brief Calculate air density
   function calculate_air_density(temperature, pressure, humidity) result(air_density)
      real(fp), intent(in) :: temperature
      real(fp), intent(in) :: pressure
      real(fp), intent(in), optional :: humidity
      real(fp) :: air_density

      real(c_double) :: hum_val
      logical(c_bool) :: use_hum

      hum_val = 0.0_c_double
      use_hum = .false.
      if (present(humidity)) then
         hum_val = real(humidity, c_double)
         use_hum = .true.
      end if

      air_density = real(catchem_calculate_air_density(real(temperature, c_double), &
         real(pressure, c_double), hum_val, use_hum), fp)
   end function calculate_air_density

   !> \brief Calculate molecular weight from formula
   function calculate_molecular_weight(formula) result(mw)
      character(len=*), intent(in) :: formula
      real(fp) :: mw

      mw = real(catchem_calculate_molecular_weight(trim(formula) // c_null_char), fp)
   end function calculate_molecular_weight

   ! Standard Imperial Conversion Proxies
   function convert_imperial_length(length_in, input_units, output_units, rc) result(length_out)
      real(fp), intent(in) :: length_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: length_out
      integer(c_int) :: rc_c
      length_out = real(catchem_convert_imperial(real(length_in, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, "LENGTH" // c_null_char, rc_c), fp)
      rc = int(rc_c)
   end function convert_imperial_length

   function convert_imperial_area(area_in, input_units, output_units, rc) result(area_out)
      real(fp), intent(in) :: area_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: area_out
      integer(c_int) :: rc_c
      area_out = real(catchem_convert_imperial(real(area_in, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, "AREA" // c_null_char, rc_c), fp)
      rc = int(rc_c)
   end function convert_imperial_area

   function convert_imperial_volume(volume_in, input_units, output_units, rc) result(volume_out)
      real(fp), intent(in) :: volume_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: volume_out
      integer(c_int) :: rc_c
      volume_out = real(catchem_convert_imperial(real(volume_in, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, "VOLUME" // c_null_char, rc_c), fp)
      rc = int(rc_c)
   end function convert_imperial_volume

   function convert_imperial_speed(speed_in, input_units, output_units, rc) result(speed_out)
      real(fp), intent(in) :: speed_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: speed_out
      integer(c_int) :: rc_c
      speed_out = real(catchem_convert_imperial(real(speed_in, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, "SPEED" // c_null_char, rc_c), fp)
      rc = int(rc_c)
   end function convert_imperial_speed

   function convert_imperial_force(force_in, input_units, output_units, rc) result(force_out)
      real(fp), intent(in) :: force_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: force_out
      integer(c_int) :: rc_c
      force_out = real(catchem_convert_imperial(real(force_in, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, "FORCE" // c_null_char, rc_c), fp)
      rc = int(rc_c)
   end function convert_imperial_force

   function convert_imperial_pressure(pressure_in, input_units, output_units, rc) result(pressure_out)
      real(fp), intent(in) :: pressure_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: pressure_out
      integer(c_int) :: rc_c
      pressure_out = real(catchem_convert_imperial(real(pressure_in, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, "PRESSURE" // c_null_char, rc_c), fp)
      rc = int(rc_c)
   end function convert_imperial_pressure

   function convert_imperial_temperature(temp_in, input_units, output_units, rc) result(temp_out)
      real(fp), intent(in) :: temp_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: temp_out
      integer(c_int) :: rc_c
      temp_out = real(catchem_convert_imperial(real(temp_in, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, "TEMPERATURE" // c_null_char, rc_c), fp)
      rc = int(rc_c)
   end function convert_imperial_temperature

   function convert_imperial_mass(mass_in, input_units, output_units, rc) result(mass_out)
      real(fp), intent(in) :: mass_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: mass_out
      integer(c_int) :: rc_c
      mass_out = real(catchem_convert_imperial(real(mass_in, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, "MASS" // c_null_char, rc_c), fp)
      rc = int(rc_c)
   end function convert_imperial_mass

   function convert_imperial_energy(energy_in, input_units, output_units, rc) result(energy_out)
      real(fp), intent(in) :: energy_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: energy_out
      integer(c_int) :: rc_c
      energy_out = real(catchem_convert_imperial(real(energy_in, c_double), &
         trim(input_units) // c_null_char, trim(output_units) // c_null_char, "ENERGY" // c_null_char, rc_c), fp)
      rc = int(rc_c)
   end function convert_imperial_energy

   subroutine convert_process_concentration_units(values, from_units, to_units, &
      molecular_weight, temperature, pressure, rc)
      real(fp), intent(inout) :: values(:)
      character(len=*), intent(in) :: from_units, to_units
      real(fp), intent(in), optional :: molecular_weight, temperature, pressure
      integer, intent(out) :: rc

      real(fp) :: mw, temp, press
      integer(c_int) :: rc_c

      mw = 28.9644_fp
      if (present(molecular_weight)) mw = molecular_weight

      temp = STANDARD_TEMP
      if (present(temperature)) temp = temperature

      press = STANDARD_PRESS
      if (present(pressure)) press = pressure

      rc_c = int(catchem_convert_process_concentration_units(values, size(values), &
         trim(from_units) // c_null_char, trim(to_units) // c_null_char, &
         mw, temp, press))
      rc = int(rc_c)
   end subroutine convert_process_concentration_units

   subroutine convert_process_flux_units(flux_values, from_units, to_units, molecular_weight, rc)
      real(fp), intent(inout) :: flux_values(:)
      character(len=*), intent(in) :: from_units, to_units
      real(fp), intent(in) :: molecular_weight
      integer, intent(out) :: rc

      integer(c_int) :: rc_c

      rc_c = int(catchem_convert_process_flux_units(flux_values, size(flux_values), &
         trim(from_units) // c_null_char, trim(to_units) // c_null_char, &
         molecular_weight))
      rc = int(rc_c)
   end subroutine convert_process_flux_units

   !========================================================================
   ! UnitConverterType Implementation
   ! =========================================================================

   subroutine converter_init(this, temperature, pressure)
      class(UnitConverterType), intent(inout) :: this
      real(fp), intent(in), optional :: temperature
      real(fp), intent(in), optional :: pressure

      if (present(temperature)) then
         this%temperature = temperature
      else
         this%temperature = STANDARD_TEMP
      end if

      if (present(pressure)) then
         this%pressure = pressure
      else
         this%pressure = STANDARD_PRESS
      end if

      this%air_density = calculate_air_density(this%temperature, this%pressure)
   end subroutine converter_init

   subroutine converter_set_conditions(this, temperature, pressure, humidity)
      class(UnitConverterType), intent(inout) :: this
      real(fp), intent(in) :: temperature
      real(fp), intent(in) :: pressure
      real(fp), intent(in), optional :: humidity

      this%temperature = temperature
      this%pressure = pressure
      if (present(humidity)) then
         this%air_density = calculate_air_density(temperature, pressure, humidity)
      else
         this%air_density = calculate_air_density(temperature, pressure)
      end if
   end subroutine converter_set_conditions

   function converter_get_air_density(this) result(air_density)
      class(UnitConverterType), intent(in) :: this
      real(fp) :: air_density
      air_density = this%air_density
   end function converter_get_air_density

   function converter_calculate_number_density(this) result(number_density)
      class(UnitConverterType), intent(in) :: this
      real(fp) :: number_density
      real(fp) :: boltz_val
      ! Boltz constant J/K
      boltz_val = 1.380649e-23_fp
      number_density = this%pressure / (boltz_val * this%temperature) * 1.0e-6_fp
   end function converter_calculate_number_density

   function converter_ppbv_to_ugm3(this, ppbv, molecular_weight) result(ugm3)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: ppbv
      real(fp), intent(in) :: molecular_weight
      real(fp) :: ugm3
      integer :: rc
      call convert_concentration(ppbv, "ppbv", "ug/m3", molecular_weight, this%temperature, this%pressure, ugm3, rc)
   end function converter_ppbv_to_ugm3

   function converter_ugm3_to_ppbv(this, ugm3, molecular_weight) result(ppbv)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: ugm3
      real(fp), intent(in) :: molecular_weight
      real(fp) :: ppbv
      integer :: rc
      call convert_concentration(ugm3, "ug/m3", "ppbv", molecular_weight, this%temperature, this%pressure, ppbv, rc)
   end function converter_ugm3_to_ppbv

   function converter_molcm3_to_ppbv(this, molcm3, temperature, pressure) result(ppbv)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: molcm3
      real(fp), intent(in) :: temperature
      real(fp), intent(in) :: pressure
      real(fp) :: ppbv
      integer :: rc
      call convert_concentration(molcm3, "molec/cm3", "ppbv", 1.0_fp, temperature, pressure, ppbv, rc)
   end function converter_molcm3_to_ppbv

   function converter_ppbv_to_molcm3(this, ppbv, temperature, pressure) result(molcm3)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: ppbv
      real(fp), intent(in) :: temperature
      real(fp), intent(in) :: pressure
      real(fp) :: molcm3
      integer :: rc
      call convert_concentration(ppbv, "ppbv", "molec/cm3", 1.0_fp, temperature, pressure, molcm3, rc)
   end function converter_ppbv_to_molcm3

   function converter_ppmv_to_mgm3(this, ppmv, molecular_weight) result(mgm3)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: ppmv
      real(fp), intent(in) :: molecular_weight
      real(fp) :: mgm3
      integer :: rc
      call convert_concentration(ppmv, "ppmv", "mg/m3", molecular_weight, this%temperature, this%pressure, mgm3, rc)
   end function converter_ppmv_to_mgm3

   function converter_mgm3_to_ppmv(this, mgm3, molecular_weight) result(ppmv)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: mgm3
      real(fp), intent(in) :: molecular_weight
      real(fp) :: ppmv
      integer :: rc
      call convert_concentration(mgm3, "mg/m3", "ppmv", molecular_weight, this%temperature, this%pressure, ppmv, rc)
   end function converter_mgm3_to_ppmv

   function converter_calculate_column_mass(this, concentrations, layer_heights, &
      molecular_weight) result(column_mass)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: concentrations(:)
      real(fp), intent(in) :: layer_heights(:)
      real(fp), intent(in) :: molecular_weight
      real(fp) :: column_mass
      integer :: k, n
      real(fp) :: mass_density

      column_mass = 0.0_fp
      n = min(size(concentrations), size(layer_heights))
      do k = 1, n
         ! Convert from ppbv to kg/m³
         ! concentrations(k) * molecular_weight * pressure / (RSTARG * temperature) * 1.0e-12
         mass_density = concentrations(k) * molecular_weight * this%pressure / &
            (8.314462618_fp * this%temperature) * 1.0e-12_fp
         column_mass = column_mass + mass_density * layer_heights(k)
      end do
   end function converter_calculate_column_mass

   function converter_calculate_dobson_units(this, concentrations, layer_heights) result(dobson_units)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: concentrations(:)
      real(fp), intent(in) :: layer_heights(:)
      real(fp) :: dobson_units
      integer :: k, n
      real(fp) :: num_density, column_density

      column_density = 0.0_fp
      n = min(size(concentrations), size(layer_heights))
      do k = 1, n
         ! Convert ppbv to molecules/cm³
         num_density = concentrations(k) * this%pressure / (1.380649e-23_fp * this%temperature) * 1.0e-15_fp
         ! Layer height is in meters, convert to cm (1 m = 100 cm)
         column_density = column_density + num_density * (layer_heights(k) * 100.0_fp)
      end do

      ! Convert to Dobson Units (1 DU = 2.687e16 molecules/cm²)
      dobson_units = column_density / 2.687e16_fp
   end function converter_calculate_dobson_units

   function converter_molcm2s_to_kgm2s(this, molcm2s, molecular_weight) result(kgm2s)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: molcm2s
      real(fp), intent(in) :: molecular_weight
      real(fp) :: kgm2s
      integer :: rc
      kgm2s = convert_flux(molcm2s, "molec/cm2/s", "kg/m2/s", molecular_weight, rc)
   end function converter_molcm2s_to_kgm2s

   function converter_kgm2s_to_molcm2s(this, kgm2s, molecular_weight) result(molcm2s)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: kgm2s
      real(fp), intent(in) :: molecular_weight
      real(fp) :: molcm2s
      integer :: rc
      molcm2s = convert_flux(kgm2s, "kg/m2/s", "molec/cm2/s", molecular_weight, rc)
   end function converter_kgm2s_to_molcm2s

   function converter_convert_rate_units(this, rate_in, input_units, output_units, rc) result(rate_out)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: rate_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: rate_out
      rate_out = convert_rate_constant(rate_in, input_units, output_units, rc)
   end function converter_convert_rate_units

end module UnitConversion_Mod
