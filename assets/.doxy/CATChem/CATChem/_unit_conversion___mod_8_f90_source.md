

# File UnitConversion\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**UnitConversion\_Mod.F90**](_unit_conversion___mod_8_f90.md)

[Go to the documentation of this file](_unit_conversion___mod_8_f90.md)


```Fortran

module unitconversion_mod
   use precision_mod, only: fp
   use error_mod, only: cc_success, cc_failure
   use constants, only: avo, rstarg, boltz, airmw, h2omw, atm

   implicit none
   private

   public :: unitconvertertype
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

   ! Standard conditions using constants from Constants module
   real(fp), parameter :: STANDARD_TEMP = 273.15_fp  
   real(fp), parameter :: STANDARD_PRESS = atm       

   ! Conversion factors
   real(fp), parameter :: PA_TO_HPA = 0.01_fp        
   real(fp), parameter :: PA_TO_ATM = 1.0_fp/atm     
   real(fp), parameter :: PA_TO_TORR = 0.007501_fp   
   real(fp), parameter :: PPB_TO_PPMV = 1.0e-3_fp    
   real(fp), parameter :: PPT_TO_PPBV = 1.0e-3_fp    

   ! Imperial unit conversion factors
   real(fp), parameter :: INCH_TO_M = 0.0254_fp      
   real(fp), parameter :: FOOT_TO_M = 0.3048_fp      
   real(fp), parameter :: YARD_TO_M = 0.9144_fp      
   real(fp), parameter :: MILE_TO_M = 1609.344_fp    
   real(fp), parameter :: SQFT_TO_M2 = 0.09290304_fp 
   real(fp), parameter :: ACRE_TO_M2 = 4046.856_fp   
   real(fp), parameter :: CUFT_TO_M3 = 0.02831685_fp 
   real(fp), parameter :: GALLON_TO_M3 = 0.003785412_fp 
   real(fp), parameter :: MPH_TO_MS = 0.44704_fp     
   real(fp), parameter :: FTS_TO_MS = 0.3048_fp      
   real(fp), parameter :: KNOT_TO_MS = 0.514444_fp   
   real(fp), parameter :: LBF_TO_N = 4.448222_fp     
   real(fp), parameter :: PSI_TO_PA = 6894.757_fp    
   real(fp), parameter :: INHG_TO_PA = 3386.389_fp   
   real(fp), parameter :: LB_TO_KG = 0.4535924_fp    
   real(fp), parameter :: OZ_TO_KG = 0.02834952_fp   
   real(fp), parameter :: TON_TO_KG = 907.1847_fp    
   real(fp), parameter :: BTU_TO_J = 1055.056_fp     
   real(fp), parameter :: CALORIE_TO_J = 4.184_fp    

   ! Imperial to metric conversion factors
   real(fp), parameter :: FT_TO_M = 0.3048_fp        
   real(fp), parameter :: INCH_TO_M = 0.0254_fp      
   real(fp), parameter :: MILE_TO_M = 1609.344_fp    
   real(fp), parameter :: YD_TO_M = 0.9144_fp        
   real(fp), parameter :: FT2_TO_M2 = 0.092903_fp    
   real(fp), parameter :: ACRE_TO_M2 = 4046.86_fp    
   real(fp), parameter :: FT3_TO_M3 = 0.028317_fp    
   real(fp), parameter :: GAL_TO_M3 = 0.003785_fp    
   real(fp), parameter :: MPH_TO_MS = 0.44704_fp     
   real(fp), parameter :: KNOT_TO_MS = 0.514444_fp   
   real(fp), parameter :: LB_TO_KG = 0.453592_fp     
   real(fp), parameter :: OZ_TO_KG = 0.0283495_fp    

   type :: unitconvertertype
      private
      real(fp) :: temperature = standard_temp    
      real(fp) :: pressure = standard_press      
      real(fp) :: air_density = 1.225_fp         
      logical :: use_standard_conditions = .true. 

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
   end type unitconvertertype

contains

   !========================================================================
   ! Standalone conversion functions
   !========================================================================

   subroutine convert_concentration(input_value, input_units, output_units, &
      molecular_weight, temperature, pressure, &
      output_value, rc)
      real(fp), intent(in) :: input_value
      character(len=*), intent(in) :: input_units, output_units
      real(fp), intent(in) :: molecular_weight
      real(fp), intent(in) :: temperature
      real(fp), intent(in) :: pressure
      real(fp), intent(out) :: output_value
      integer, intent(out) :: rc

      type(UnitConverterType) :: converter

      rc = cc_success

      call converter%init(temperature, pressure)

      select case (trim(input_units) // ' -> ' // trim(output_units))
       case ('ppbv -> ug/m3', 'ppb -> ug/m3')
         output_value = converter%ppbv_to_ugm3(input_value, molecular_weight)
       case ('ug/m3 -> ppbv', 'ug/m3 -> ppb')
         output_value = converter%ugm3_to_ppbv(input_value, molecular_weight)
       case ('ppmv -> mg/m3', 'ppm -> mg/m3')
         output_value = converter%ppmv_to_mgm3(input_value, molecular_weight)
       case ('mg/m3 -> ppmv', 'mg/m3 -> ppm')
         output_value = converter%mgm3_to_ppmv(input_value, molecular_weight)
       case ('molec/cm3 -> ppbv')
         output_value = converter%molcm3_to_ppbv(input_value, temperature, pressure)
       case ('ppbv -> molec/cm3')
         output_value = converter%ppbv_to_molcm3(input_value, temperature, pressure)
       case default
         rc = cc_failure
         output_value = input_value
      end select

   end subroutine convert_concentration

   function convert_pressure(pressure_in, input_units, output_units, rc) result(pressure_out)
      real(fp), intent(in) :: pressure_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: pressure_out

      real(fp) :: pressure_pa

      rc = cc_success

      ! First convert to Pa
      select case (trim(input_units))
       case ('Pa')
         pressure_pa = pressure_in
       case ('hPa', 'mb', 'mbar')
         pressure_pa = pressure_in / pa_to_hpa
       case ('atm')
         pressure_pa = pressure_in / pa_to_atm
       case ('Torr', 'mmHg')
         pressure_pa = pressure_in / pa_to_torr
       case ('psi')
         pressure_pa = pressure_in * 6894.76_fp
       case default
         rc = cc_failure
         pressure_out = pressure_in
         return
      end select

      ! Then convert from Pa to output units
      select case (trim(output_units))
       case ('Pa')
         pressure_out = pressure_pa
       case ('hPa', 'mb', 'mbar')
         pressure_out = pressure_pa * pa_to_hpa
       case ('atm')
         pressure_out = pressure_pa * pa_to_atm
       case ('Torr', 'mmHg')
         pressure_out = pressure_pa * pa_to_torr
       case ('psi')
         pressure_out = pressure_pa / 6894.76_fp
       case default
         rc = cc_failure
         pressure_out = pressure_in
      end select

   end function convert_pressure

   function convert_temperature(temp_in, input_units, output_units, rc) result(temp_out)
      real(fp), intent(in) :: temp_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: temp_out

      real(fp) :: temp_k

      rc = cc_success

      ! First convert to Kelvin
      select case (trim(input_units))
       case ('K', 'Kelvin')
         temp_k = temp_in
       case ('C', 'Celsius')
         temp_k = temp_in + 273.15_fp
       case ('F', 'Fahrenheit')
         temp_k = (temp_in - 32.0_fp) * 5.0_fp/9.0_fp + 273.15_fp
       case default
         rc = cc_failure
         temp_out = temp_in
         return
      end select

      ! Then convert from Kelvin to output units
      select case (trim(output_units))
       case ('K', 'Kelvin')
         temp_out = temp_k
       case ('C', 'Celsius')
         temp_out = temp_k - 273.15_fp
       case ('F', 'Fahrenheit')
         temp_out = (temp_k - 273.15_fp) * 9.0_fp/5.0_fp + 32.0_fp
       case default
         rc = cc_failure
         temp_out = temp_in
      end select

   end function convert_temperature

   function convert_flux(flux_in, input_units, output_units, molecular_weight, rc) result(flux_out)
      real(fp), intent(in) :: flux_in
      character(len=*), intent(in) :: input_units, output_units
      real(fp), intent(in) :: molecular_weight
      integer, intent(out) :: rc
      real(fp) :: flux_out

      real(fp) :: flux_kgm2s

      rc = cc_success

      ! First convert to kg/m²/s
      select case (trim(input_units))
       case ('kg/m2/s')
         flux_kgm2s = flux_in
       case ('g/m2/s')
         flux_kgm2s = flux_in * 1.0e-3_fp
       case ('mg/m2/s')
         flux_kgm2s = flux_in * 1.0e-6_fp
       case ('ug/m2/s')
         flux_kgm2s = flux_in * 1.0e-9_fp
       case ('mol/cm2/s')
         flux_kgm2s = flux_in * molecular_weight * 1.0e-3_fp * 1.0e4_fp
       case ('molec/cm2/s')
         flux_kgm2s = flux_in * molecular_weight / avo * 1.0e-3_fp * 1.0e4_fp
       case default
         rc = cc_failure
         flux_out = flux_in
         return
      end select

      ! Then convert from kg/m²/s to output units
      select case (trim(output_units))
       case ('kg/m2/s')
         flux_out = flux_kgm2s
       case ('g/m2/s')
         flux_out = flux_kgm2s * 1.0e3_fp
       case ('mg/m2/s')
         flux_out = flux_kgm2s * 1.0e6_fp
       case ('ug/m2/s')
         flux_out = flux_kgm2s * 1.0e9_fp
       case ('mol/cm2/s')
         flux_out = flux_kgm2s / molecular_weight * 1.0e3_fp * 1.0e-4_fp
       case ('molec/cm2/s')
         flux_out = flux_kgm2s / molecular_weight * avo * 1.0e3_fp * 1.0e-4_fp
       case default
         rc = cc_failure
         flux_out = flux_in
      end select

   end function convert_flux

   function convert_rate_constant(rate_in, input_units, output_units, rc) result(rate_out)
      real(fp), intent(in) :: rate_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: rate_out

      real(fp) :: rate_s

      rc = cc_success

      ! First convert to s⁻¹
      select case (trim(input_units))
       case ('s-1', '/s', '1/s')
         rate_s = rate_in
       case ('min-1', '/min', '1/min')
         rate_s = rate_in / 60.0_fp
       case ('hr-1', 'h-1', '/hr', '/h', '1/hr', '1/h')
         rate_s = rate_in / 3600.0_fp
       case ('day-1', 'd-1', '/day', '/d', '1/day', '1/d')
         rate_s = rate_in / 86400.0_fp
       case default
         rc = cc_failure
         rate_out = rate_in
         return
      end select

      ! Then convert from s⁻¹ to output units
      select case (trim(output_units))
       case ('s-1', '/s', '1/s')
         rate_out = rate_s
       case ('min-1', '/min', '1/min')
         rate_out = rate_s * 60.0_fp
       case ('hr-1', 'h-1', '/hr', '/h', '1/hr', '1/h')
         rate_out = rate_s * 3600.0_fp
       case ('day-1', 'd-1', '/day', '/d', '1/day', '1/d')
         rate_out = rate_s * 86400.0_fp
       case default
         rc = cc_failure
         rate_out = rate_in
      end select

   end function convert_rate_constant

   function convert_mass_units(mass_in, input_units, output_units, rc) result(mass_out)
      real(fp), intent(in) :: mass_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: mass_out

      real(fp) :: mass_kg

      rc = cc_success

      ! First convert to kg
      select case (trim(input_units))
       case ('kg')
         mass_kg = mass_in
       case ('g')
         mass_kg = mass_in * 1.0e-3_fp
       case ('mg')
         mass_kg = mass_in * 1.0e-6_fp
       case ('ug')
         mass_kg = mass_in * 1.0e-9_fp
       case ('ng')
         mass_kg = mass_in * 1.0e-12_fp
       case ('Tg')
         mass_kg = mass_in * 1.0e9_fp
       case ('Gg')
         mass_kg = mass_in * 1.0e6_fp
       case default
         rc = cc_failure
         mass_out = mass_in
         return
      end select

      ! Then convert from kg to output units
      select case (trim(output_units))
       case ('kg')
         mass_out = mass_kg
       case ('g')
         mass_out = mass_kg * 1.0e3_fp
       case ('mg')
         mass_out = mass_kg * 1.0e6_fp
       case ('ug')
         mass_out = mass_kg * 1.0e9_fp
       case ('ng')
         mass_out = mass_kg * 1.0e12_fp
       case ('Tg')
         mass_out = mass_kg * 1.0e-9_fp
       case ('Gg')
         mass_out = mass_kg * 1.0e-6_fp
       case default
         rc = cc_failure
         mass_out = mass_in
      end select

   end function convert_mass_units

   function calculate_air_density(temperature, pressure, humidity) result(air_density)
      real(fp), intent(in) :: temperature
      real(fp), intent(in) :: pressure
      real(fp), intent(in), optional :: humidity
      real(fp) :: air_density

      real(fp) :: rh, p_sat, p_dry
      real(fp), parameter :: MW_H2O = 18.015_fp  

      if (present(humidity)) then
         rh = humidity
         ! Calculate saturation vapor pressure (simplified Antoine equation)
         p_sat = 610.78_fp * exp(17.27_fp * (temperature - 273.15_fp) / (temperature - 35.86_fp))
         p_dry = pressure - rh * p_sat
         ! Calculate density accounting for humidity
         air_density = (p_dry * airmw + rh * p_sat * h2omw) / (rstarg * temperature) * 1.0e-3_fp
      else
         ! Dry air density
         air_density = pressure * airmw / (rstarg * temperature) * 1.0e-3_fp
      endif

   end function calculate_air_density

   function calculate_molecular_weight(formula) result(mw)
      character(len=*), intent(in) :: formula
      real(fp) :: mw

      ! Simplified implementation - actual implementation would parse chemical formula
      ! For now, return common molecular weights based on formula
      select case (trim(formula))
       case ('O3')
         mw = 48.0_fp
       case ('NO2')
         mw = 46.0_fp
       case ('NO')
         mw = 30.0_fp
       case ('CO')
         mw = 28.0_fp
       case ('SO2')
         mw = 64.1_fp
       case ('NH3')
         mw = 17.0_fp
       case ('CH4')
         mw = 16.0_fp
       case ('H2O')
         mw = 18.0_fp
       case ('CO2')
         mw = 44.0_fp
       case default
         mw = airmw  ! Default to air molecular weight
      end select

   end function calculate_molecular_weight

   !========================================================================
   ! UnitConverterType Implementation
   !========================================================================

   subroutine converter_init(this, temperature, pressure)
      class(UnitConverterType), intent(inout) :: this
      real(fp), intent(in), optional :: temperature
      real(fp), intent(in), optional :: pressure

      if (present(temperature)) then
         this%temperature = temperature
         this%use_standard_conditions = .false.
      else
         this%temperature = standard_temp
      endif

      if (present(pressure)) then
         this%pressure = pressure
         this%use_standard_conditions = .false.
      else
         this%pressure = standard_press
      endif

      ! Calculate air density
      this%air_density = calculate_air_density(this%temperature, this%pressure)

   end subroutine converter_init

   subroutine converter_set_conditions(this, temperature, pressure, humidity)
      class(UnitConverterType), intent(inout) :: this
      real(fp), intent(in) :: temperature
      real(fp), intent(in) :: pressure
      real(fp), intent(in), optional :: humidity

      this%temperature = temperature
      this%pressure = pressure
      this%use_standard_conditions = .false.

      if (present(humidity)) then
         this%air_density = calculate_air_density(temperature, pressure, humidity)
      else
         this%air_density = calculate_air_density(temperature, pressure)
      endif

   end subroutine converter_set_conditions

   function converter_get_air_density(this) result(air_density)
      class(UnitConverterType), intent(in) :: this
      real(fp) :: air_density

      air_density = this%air_density

   end function converter_get_air_density

   function converter_calculate_number_density(this) result(number_density)
      class(UnitConverterType), intent(in) :: this
      real(fp) :: number_density

      number_density = this%pressure / (boltz * this%temperature) * 1.0e-6_fp

   end function converter_calculate_number_density

   function converter_ppbv_to_ugm3(this, ppbv, molecular_weight) result(ugm3)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: ppbv
      real(fp), intent(in) :: molecular_weight
      real(fp) :: ugm3

      ugm3 = ppbv * molecular_weight * this%pressure / (rstarg * this%temperature) * 1.0e-3_fp

   end function converter_ppbv_to_ugm3

   function converter_ugm3_to_ppbv(this, ugm3, molecular_weight) result(ppbv)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: ugm3
      real(fp), intent(in) :: molecular_weight
      real(fp) :: ppbv

      ppbv = ugm3 * rstarg * this%temperature / (molecular_weight * this%pressure) * 1.0e3_fp

   end function converter_ugm3_to_ppbv

   function converter_molcm3_to_ppbv(this, molcm3, temperature, pressure) result(ppbv)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: molcm3
      real(fp), intent(in) :: temperature
      real(fp), intent(in) :: pressure
      real(fp) :: ppbv

      real(fp) :: number_density

      number_density = pressure / (boltzmann * temperature) * 1.0e-6_fp
      ppbv = molcm3 / number_density * 1.0e9_fp

   end function converter_molcm3_to_ppbv

   function converter_ppbv_to_molcm3(this, ppbv, temperature, pressure) result(molcm3)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: ppbv
      real(fp), intent(in) :: temperature
      real(fp), intent(in) :: pressure
      real(fp) :: molcm3

      real(fp) :: number_density

      number_density = pressure / (boltzmann * temperature) * 1.0e-6_fp
      molcm3 = ppbv * number_density * 1.0e-9_fp

   end function converter_ppbv_to_molcm3

   function converter_ppmv_to_mgm3(this, ppmv, molecular_weight) result(mgm3)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: ppmv
      real(fp), intent(in) :: molecular_weight
      real(fp) :: mgm3

      mgm3 = ppmv * molecular_weight * this%pressure / (r_gas * this%temperature)

   end function converter_ppmv_to_mgm3

   function converter_mgm3_to_ppmv(this, mgm3, molecular_weight) result(ppmv)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: mgm3
      real(fp), intent(in) :: molecular_weight
      real(fp) :: ppmv

      ppmv = mgm3 * r_gas * this%temperature / (molecular_weight * this%pressure)

   end function converter_mgm3_to_ppmv

   function converter_calculate_column_mass(this, concentrations, layer_heights, &
      molecular_weight) result(column_mass)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: concentrations(:)
      real(fp), intent(in) :: layer_heights(:)
      real(fp), intent(in) :: molecular_weight
      real(fp) :: column_mass

      integer :: k
      real(fp) :: mass_density

      column_mass = 0.0_fp
      do k = 1, size(concentrations)
         ! Convert ppbv to kg/m³
         mass_density = concentrations(k) * molecular_weight * this%pressure / &
            (r_gas * this%temperature) * 1.0e-12_fp
         column_mass = column_mass + mass_density * layer_heights(k)
      end do

   end function converter_calculate_column_mass

   function converter_calculate_dobson_units(this, concentrations, layer_heights) result(dobson_units)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: concentrations(:)
      real(fp), intent(in) :: layer_heights(:)
      real(fp) :: dobson_units

      integer :: k
      real(fp) :: number_density, column_density

      column_density = 0.0_fp
      do k = 1, size(concentrations)
         number_density = this%pressure / (boltzmann * this%temperature) * 1.0e-6_fp
         column_density = column_density + concentrations(k) * number_density * &
            layer_heights(k) * 1.0e-9_fp * 1.0e2_fp
      end do

      ! Convert to Dobson Units (1 DU = 2.687 × 10¹⁶ molecules/cm²)
      dobson_units = column_density / 2.687e16_fp

   end function converter_calculate_dobson_units

   function converter_molcm2s_to_kgm2s(this, molcm2s, molecular_weight) result(kgm2s)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: molcm2s
      real(fp), intent(in) :: molecular_weight
      real(fp) :: kgm2s

      kgm2s = molcm2s * molecular_weight * 1.0e-3_fp * 1.0e4_fp

   end function converter_molcm2s_to_kgm2s

   function converter_kgm2s_to_molcm2s(this, kgm2s, molecular_weight) result(molcm2s)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: kgm2s
      real(fp), intent(in) :: molecular_weight
      real(fp) :: molcm2s

      molcm2s = kgm2s / molecular_weight * 1.0e3_fp * 1.0e-4_fp

   end function converter_kgm2s_to_molcm2s

   function converter_convert_rate_units(this, rate_in, input_units, output_units, rc) result(rate_out)
      class(UnitConverterType), intent(in) :: this
      real(fp), intent(in) :: rate_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: rate_out

      rate_out = convert_rate_constant(rate_in, input_units, output_units, rc)

   end function converter_convert_rate_units

   !========================================================================
   ! Imperial Unit Conversion Functions
   !========================================================================

   function convert_imperial_length(length_in, input_units, output_units, rc) result(length_out)
      real(fp), intent(in) :: length_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: length_out

      real(fp) :: length_m

      rc = cc_success

      ! First convert to meters
      select case (trim(input_units))
       case ('m', 'meter', 'metre')
         length_m = length_in
       case ('cm', 'centimeter')
         length_m = length_in * 0.01_fp
       case ('mm', 'millimeter')
         length_m = length_in * 0.001_fp
       case ('km', 'kilometer')
         length_m = length_in * 1000.0_fp
       case ('in', 'inch', 'inches')
         length_m = length_in * inch_to_m
       case ('ft', 'foot', 'feet')
         length_m = length_in * foot_to_m
       case ('yd', 'yard', 'yards')
         length_m = length_in * yard_to_m
       case ('mi', 'mile', 'miles')
         length_m = length_in * mile_to_m
       case ('nmi', 'nautical_mile')
         length_m = length_in * 1852.0_fp
       case default
         rc = cc_failure
         length_out = length_in
         return
      end select

      ! Then convert from meters to output units
      select case (trim(output_units))
       case ('m', 'meter', 'metre')
         length_out = length_m
       case ('cm', 'centimeter')
         length_out = length_m / 0.01_fp
       case ('mm', 'millimeter')
         length_out = length_m / 0.001_fp
       case ('km', 'kilometer')
         length_out = length_m / 1000.0_fp
       case ('in', 'inch', 'inches')
         length_out = length_m / inch_to_m
       case ('ft', 'foot', 'feet')
         length_out = length_m / foot_to_m
       case ('yd', 'yard', 'yards')
         length_out = length_m / yard_to_m
       case ('mi', 'mile', 'miles')
         length_out = length_m / mile_to_m
       case ('nmi', 'nautical_mile')
         length_out = length_m / 1852.0_fp
       case default
         rc = cc_failure
         length_out = length_in
      end select

   end function convert_imperial_length

   function convert_imperial_area(area_in, input_units, output_units, rc) result(area_out)
      real(fp), intent(in) :: area_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: area_out

      real(fp) :: area_m2

      rc = cc_success

      ! First convert to m²
      select case (trim(input_units))
       case ('m2', 'm^2', 'sq_m')
         area_m2 = area_in
       case ('cm2', 'cm^2', 'sq_cm')
         area_m2 = area_in * 1.0e-4_fp
       case ('km2', 'km^2', 'sq_km')
         area_m2 = area_in * 1.0e6_fp
       case ('in2', 'in^2', 'sq_in')
         area_m2 = area_in * (inch_to_m**2)
       case ('ft2', 'ft^2', 'sq_ft')
         area_m2 = area_in * sqft_to_m2
       case ('yd2', 'yd^2', 'sq_yd')
         area_m2 = area_in * (yard_to_m**2)
       case ('acre', 'acres')
         area_m2 = area_in * acre_to_m2
       case ('mi2', 'mi^2', 'sq_mi')
         area_m2 = area_in * (mile_to_m**2)
       case default
         rc = cc_failure
         area_out = area_in
         return
      end select

      ! Then convert from m² to output units
      select case (trim(output_units))
       case ('m2', 'm^2', 'sq_m')
         area_out = area_m2
       case ('cm2', 'cm^2', 'sq_cm')
         area_out = area_m2 / 1.0e-4_fp
       case ('km2', 'km^2', 'sq_km')
         area_out = area_m2 / 1.0e6_fp
       case ('in2', 'in^2', 'sq_in')
         area_out = area_m2 / (inch_to_m**2)
       case ('ft2', 'ft^2', 'sq_ft')
         area_out = area_m2 / sqft_to_m2
       case ('yd2', 'yd^2', 'sq_yd')
         area_out = area_m2 / (yard_to_m**2)
       case ('acre', 'acres')
         area_out = area_m2 / acre_to_m2
       case ('mi2', 'mi^2', 'sq_mi')
         area_out = area_m2 / (mile_to_m**2)
       case default
         rc = cc_failure
         area_out = area_in
      end select

   end function convert_imperial_area

   function convert_imperial_volume(volume_in, input_units, output_units, rc) result(volume_out)
      real(fp), intent(in) :: volume_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: volume_out

      real(fp) :: volume_m3

      rc = cc_success

      ! First convert to m³
      select case (trim(input_units))
       case ('m3', 'm^3', 'cu_m')
         volume_m3 = volume_in
       case ('cm3', 'cm^3', 'cu_cm')
         volume_m3 = volume_in * 1.0e-6_fp
       case ('L', 'liter', 'litre')
         volume_m3 = volume_in * 1.0e-3_fp
       case ('mL', 'ml', 'milliliter')
         volume_m3 = volume_in * 1.0e-6_fp
       case ('in3', 'in^3', 'cu_in')
         volume_m3 = volume_in * (inch_to_m**3)
       case ('ft3', 'ft^3', 'cu_ft')
         volume_m3 = volume_in * cuft_to_m3
       case ('yd3', 'yd^3', 'cu_yd')
         volume_m3 = volume_in * (yard_to_m**3)
       case ('gal', 'gallon', 'us_gal')
         volume_m3 = volume_in * gallon_to_m3
       case ('qt', 'quart')
         volume_m3 = volume_in * gallon_to_m3 / 4.0_fp
       case ('pt', 'pint')
         volume_m3 = volume_in * gallon_to_m3 / 8.0_fp
       case ('fl_oz', 'fluid_ounce')
         volume_m3 = volume_in * gallon_to_m3 / 128.0_fp
       case default
         rc = cc_failure
         volume_out = volume_in
         return
      end select

      ! Then convert from m³ to output units
      select case (trim(output_units))
       case ('m3', 'm^3', 'cu_m')
         volume_out = volume_m3
       case ('cm3', 'cm^3', 'cu_cm')
         volume_out = volume_m3 / 1.0e-6_fp
       case ('L', 'liter', 'litre')
         volume_out = volume_m3 / 1.0e-3_fp
       case ('mL', 'ml', 'milliliter')
         volume_out = volume_m3 / 1.0e-6_fp
       case ('in3', 'in^3', 'cu_in')
         volume_out = volume_m3 / (inch_to_m**3)
       case ('ft3', 'ft^3', 'cu_ft')
         volume_out = volume_m3 / cuft_to_m3
       case ('yd3', 'yd^3', 'cu_yd')
         volume_out = volume_m3 / (yard_to_m**3)
       case ('gal', 'gallon', 'us_gal')
         volume_out = volume_m3 / gallon_to_m3
       case ('qt', 'quart')
         volume_out = volume_m3 / gallon_to_m3 * 4.0_fp
       case ('pt', 'pint')
         volume_out = volume_m3 / gallon_to_m3 * 8.0_fp
       case ('fl_oz', 'fluid_ounce')
         volume_out = volume_m3 / gallon_to_m3 * 128.0_fp
       case default
         rc = cc_failure
         volume_out = volume_in
      end select

   end function convert_imperial_volume

   function convert_imperial_speed(speed_in, input_units, output_units, rc) result(speed_out)
      real(fp), intent(in) :: speed_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: speed_out

      real(fp) :: speed_ms

      rc = cc_success

      ! First convert to m/s
      select case (trim(input_units))
       case ('m/s', 'ms', 'mps')
         speed_ms = speed_in
       case ('km/h', 'kmh', 'kph')
         speed_ms = speed_in / 3.6_fp
       case ('cm/s', 'cms')
         speed_ms = speed_in * 0.01_fp
       case ('ft/s', 'fts', 'fps')
         speed_ms = speed_in * fts_to_ms
       case ('mph', 'mi/h')
         speed_ms = speed_in * mph_to_ms
       case ('knot', 'kn', 'kt')
         speed_ms = speed_in * knot_to_ms
       case ('in/s', 'ips')
         speed_ms = speed_in * inch_to_m
       case default
         rc = cc_failure
         speed_out = speed_in
         return
      end select

      ! Then convert from m/s to output units
      select case (trim(output_units))
       case ('m/s', 'ms', 'mps')
         speed_out = speed_ms
       case ('km/h', 'kmh', 'kph')
         speed_out = speed_ms * 3.6_fp
       case ('cm/s', 'cms')
         speed_out = speed_ms / 0.01_fp
       case ('ft/s', 'fts', 'fps')
         speed_out = speed_ms / fts_to_ms
       case ('mph', 'mi/h')
         speed_out = speed_ms / mph_to_ms
       case ('knot', 'kn', 'kt')
         speed_out = speed_ms / knot_to_ms
       case ('in/s', 'ips')
         speed_out = speed_ms / inch_to_m
       case default
         rc = cc_failure
         speed_out = speed_in
      end select

   end function convert_imperial_speed

   function convert_imperial_force(force_in, input_units, output_units, rc) result(force_out)
      real(fp), intent(in) :: force_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: force_out

      real(fp) :: force_n

      rc = cc_success

      ! First convert to Newtons
      select case (trim(input_units))
       case ('N', 'newton', 'newtons')
         force_n = force_in
       case ('kN', 'kilonewton')
         force_n = force_in * 1000.0_fp
       case ('dyne', 'dynes')
         force_n = force_in * 1.0e-5_fp
       case ('lbf', 'lb', 'pound_force')
         force_n = force_in * lbf_to_n
       case ('ozf', 'oz', 'ounce_force')
         force_n = force_in * lbf_to_n / 16.0_fp
       case ('kip', 'kips')
         force_n = force_in * lbf_to_n * 1000.0_fp
       case default
         rc = cc_failure
         force_out = force_in
         return
      end select

      ! Then convert from Newtons to output units
      select case (trim(output_units))
       case ('N', 'newton', 'newtons')
         force_out = force_n
       case ('kN', 'kilonewton')
         force_out = force_n / 1000.0_fp
       case ('dyne', 'dynes')
         force_out = force_n / 1.0e-5_fp
       case ('lbf', 'lb', 'pound_force')
         force_out = force_n / lbf_to_n
       case ('ozf', 'oz', 'ounce_force')
         force_out = force_n / lbf_to_n * 16.0_fp
       case ('kip', 'kips')
         force_out = force_n / (lbf_to_n * 1000.0_fp)
       case default
         rc = cc_failure
         force_out = force_in
      end select

   end function convert_imperial_force

   function convert_imperial_pressure(pressure_in, input_units, output_units, rc) result(pressure_out)
      real(fp), intent(in) :: pressure_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: pressure_out

      real(fp) :: pressure_pa

      rc = cc_success

      ! First convert to Pascals
      select case (trim(input_units))
       case ('Pa', 'pascal')
         pressure_pa = pressure_in
       case ('kPa', 'kilopascal')
         pressure_pa = pressure_in * 1000.0_fp
       case ('MPa', 'megapascal')
         pressure_pa = pressure_in * 1.0e6_fp
       case ('hPa', 'hectopascal', 'mb', 'mbar', 'millibar')
         pressure_pa = pressure_in * 100.0_fp
       case ('bar')
         pressure_pa = pressure_in * 1.0e5_fp
       case ('atm', 'atmosphere')
         pressure_pa = pressure_in * atm
       case ('Torr', 'mmHg')
         pressure_pa = pressure_in * 133.3224_fp
       case ('psi', 'lb/in2')
         pressure_pa = pressure_in * psi_to_pa
       case ('psf', 'lb/ft2')
         pressure_pa = pressure_in * psi_to_pa / 144.0_fp
       case ('inHg', 'in_hg')
         pressure_pa = pressure_in * inhg_to_pa
       case ('inH2O', 'in_h2o')
         pressure_pa = pressure_in * 248.84_fp
       case default
         rc = cc_failure
         pressure_out = pressure_in
         return
      end select

      ! Then convert from Pascals to output units
      select case (trim(output_units))
       case ('Pa', 'pascal')
         pressure_out = pressure_pa
       case ('kPa', 'kilopascal')
         pressure_out = pressure_pa / 1000.0_fp
       case ('MPa', 'megapascal')
         pressure_out = pressure_pa / 1.0e6_fp
       case ('hPa', 'hectopascal', 'mb', 'mbar', 'millibar')
         pressure_out = pressure_pa / 100.0_fp
       case ('bar')
         pressure_out = pressure_pa / 1.0e5_fp
       case ('atm', 'atmosphere')
         pressure_out = pressure_pa / atm
       case ('Torr', 'mmHg')
         pressure_out = pressure_pa / 133.3224_fp
       case ('psi', 'lb/in2')
         pressure_out = pressure_pa / psi_to_pa
       case ('psf', 'lb/ft2')
         pressure_out = pressure_pa / psi_to_pa * 144.0_fp
       case ('inHg', 'in_hg')
         pressure_out = pressure_pa / inhg_to_pa
       case ('inH2O', 'in_h2o')
         pressure_out = pressure_pa / 248.84_fp
       case default
         rc = cc_failure
         pressure_out = pressure_in
      end select

   end function convert_imperial_pressure

   function convert_imperial_temperature(temp_in, input_units, output_units, rc) result(temp_out)
      real(fp), intent(in) :: temp_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: temp_out

      real(fp) :: temp_k

      rc = cc_success

      ! First convert to Kelvin
      select case (trim(input_units))
       case ('K', 'Kelvin', 'kelvin')
         temp_k = temp_in
       case ('C', 'Celsius', 'celsius')
         temp_k = temp_in + 273.15_fp
       case ('F', 'Fahrenheit', 'fahrenheit')
         temp_k = (temp_in - 32.0_fp) * 5.0_fp/9.0_fp + 273.15_fp
       case ('R', 'Rankine', 'rankine')
         temp_k = temp_in * 5.0_fp/9.0_fp
       case default
         rc = cc_failure
         temp_out = temp_in
         return
      end select

      ! Then convert from Kelvin to output units
      select case (trim(output_units))
       case ('K', 'Kelvin', 'kelvin')
         temp_out = temp_k
       case ('C', 'Celsius', 'celsius')
         temp_out = temp_k - 273.15_fp
       case ('F', 'Fahrenheit', 'fahrenheit')
         temp_out = (temp_k - 273.15_fp) * 9.0_fp/5.0_fp + 32.0_fp
       case ('R', 'Rankine', 'rankine')
         temp_out = temp_k * 9.0_fp/5.0_fp
       case default
         rc = cc_failure
         temp_out = temp_in
      end select

   end function convert_imperial_temperature

   function convert_imperial_mass(mass_in, input_units, output_units, rc) result(mass_out)
      real(fp), intent(in) :: mass_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: mass_out

      real(fp) :: mass_kg

      rc = cc_success

      ! First convert to kg
      select case (trim(input_units))
       case ('kg', 'kilogram')
         mass_kg = mass_in
       case ('g', 'gram')
         mass_kg = mass_in * 1.0e-3_fp
       case ('mg', 'milligram')
         mass_kg = mass_in * 1.0e-6_fp
       case ('ug', 'microgram')
         mass_kg = mass_in * 1.0e-9_fp
       case ('lb', 'pound', 'lbs')
         mass_kg = mass_in * lb_to_kg
       case ('oz', 'ounce')
         mass_kg = mass_in * oz_to_kg
       case ('ton', 'short_ton', 'us_ton')
         mass_kg = mass_in * ton_to_kg
       case ('long_ton', 'uk_ton')
         mass_kg = mass_in * 1016.047_fp
       case ('stone')
         mass_kg = mass_in * 6.350293_fp
       case ('grain')
         mass_kg = mass_in * 6.479891e-5_fp
       case default
         rc = cc_failure
         mass_out = mass_in
         return
      end select

      ! Then convert from kg to output units
      select case (trim(output_units))
       case ('kg', 'kilogram')
         mass_out = mass_kg
       case ('g', 'gram')
         mass_out = mass_kg / 1.0e-3_fp
       case ('mg', 'milligram')
         mass_out = mass_kg / 1.0e-6_fp
       case ('ug', 'microgram')
         mass_out = mass_kg / 1.0e-9_fp
       case ('lb', 'pound', 'lbs')
         mass_out = mass_kg / lb_to_kg
       case ('oz', 'ounce')
         mass_out = mass_kg / oz_to_kg
       case ('ton', 'short_ton', 'us_ton')
         mass_out = mass_kg / ton_to_kg
       case ('long_ton', 'uk_ton')
         mass_out = mass_kg / 1016.047_fp
       case ('stone')
         mass_out = mass_kg / 6.350293_fp
       case ('grain')
         mass_out = mass_kg / 6.479891e-5_fp
       case default
         rc = cc_failure
         mass_out = mass_in
      end select

   end function convert_imperial_mass

   function convert_imperial_energy(energy_in, input_units, output_units, rc) result(energy_out)
      real(fp), intent(in) :: energy_in
      character(len=*), intent(in) :: input_units, output_units
      integer, intent(out) :: rc
      real(fp) :: energy_out

      real(fp) :: energy_j

      rc = cc_success

      ! First convert to Joules
      select case (trim(input_units))
       case ('J', 'joule', 'joules')
         energy_j = energy_in
       case ('kJ', 'kilojoule')
         energy_j = energy_in * 1000.0_fp
       case ('MJ', 'megajoule')
         energy_j = energy_in * 1.0e6_fp
       case ('cal', 'calorie', 'calories')
         energy_j = energy_in * calorie_to_j
       case ('kcal', 'kilocalorie', 'Cal')
         energy_j = energy_in * calorie_to_j * 1000.0_fp
       case ('BTU', 'btu', 'british_thermal_unit')
         energy_j = energy_in * btu_to_j
       case ('therm', 'therms')
         energy_j = energy_in * btu_to_j * 100000.0_fp
       case ('kWh', 'kilowatt_hour')
         energy_j = energy_in * 3.6e6_fp
       case ('eV', 'electron_volt')
         energy_j = energy_in * 1.602176e-19_fp
       case ('ft_lb', 'foot_pound')
         energy_j = energy_in * 1.355818_fp
       case default
         rc = cc_failure
         energy_out = energy_in
         return
      end select

      ! Then convert from Joules to output units
      select case (trim(output_units))
       case ('J', 'joule', 'joules')
         energy_out = energy_j
       case ('kJ', 'kilojoule')
         energy_out = energy_j / 1000.0_fp
       case ('MJ', 'megajoule')
         energy_out = energy_j / 1.0e6_fp
       case ('cal', 'calorie', 'calories')
         energy_out = energy_j / calorie_to_j
       case ('kcal', 'kilocalorie', 'Cal')
         energy_out = energy_j / (calorie_to_j * 1000.0_fp)
       case ('BTU', 'btu', 'british_thermal_unit')
         energy_out = energy_j / btu_to_j
       case ('therm', 'therms')
         energy_out = energy_j / (btu_to_j * 100000.0_fp)
       case ('kWh', 'kilowatt_hour')
         energy_out = energy_j / 3.6e6_fp
       case ('eV', 'electron_volt')
         energy_out = energy_j / 1.602176e-19_fp
       case ('ft_lb', 'foot_pound')
         energy_out = energy_j / 1.355818_fp
       case default
         rc = cc_failure
         energy_out = energy_in
      end select

   end function convert_imperial_energy

end module unitconversion_mod
```


