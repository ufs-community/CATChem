!> \file constants.F90
!! \brief Physical and mathematical constants for CATChem
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module provides fundamental physical and mathematical constants
!! used throughout the CATChem atmospheric chemistry modeling system.
!!
!! \details
!! The constants module defines all physical constants, conversion factors,
!! and mathematical constants used in atmospheric chemistry calculations.
!! All values are given in SI units unless otherwise specified.
!!
!! \section constants_usage Usage Example
!! \code{.f90}
!! use constants
!! real(fp) :: air_density
!! air_density = pressure / (Rd * temperature)
!! \endcode
!!
module Constants
   use precision_mod

   implicit none
   private

   ! \name Fundamental Physical Constants (must be defined first for dependencies)
   !! \brief Universal physical constants
   !! \{
   REAL(fp), PARAMETER, PUBLIC :: AVO = 6.022140857e+23_fp         !< Avogadro's number [particles/mol]
   REAL(fp), PARAMETER, PUBLIC :: g0     = 9.80665e+0_fp           !< Standard gravity acceleration [m/s^2]
   REAL(fp), PARAMETER, PUBLIC :: g0_100 = 100.0_fp / g0           !< 100 divided by standard gravity
   REAL(fp), PARAMETER, PUBLIC :: Re = 6.3710072e+6_fp             !< Earth's radius [m]
   REAL(fp), PARAMETER, PUBLIC :: RSTARG = 8.3144598_fp            !< Universal gas constant [J/K/mol]
   REAL(fp), PARAMETER, PUBLIC :: BOLTZ = 1.38064852e-23_fp        !< Boltzmann's constant [J/K]
   REAL(fp), PARAMETER, PUBLIC :: PLANCK = 6.62606957e-34_fp       !< Planck's constant [J⋅s]
   REAL(fp), PARAMETER, PUBLIC :: CCONST = 2.99792458e+8_fp        !< Speed of light in vacuum [m/s]
   ! \}

   ! \name Atmospheric Properties
   !! \brief Constants related to atmospheric composition and properties
   !! \{
   REAL(fp), PARAMETER, PUBLIC :: Cp = 1.0046e+3_fp                !< Specific heat of dry air at constant pressure [J/kg/K]
   REAL(fp), PARAMETER, PUBLIC :: Cv = 7.1760e+2_fp                !< Specific heat of dry air at constant volume [J/kg/K]
   REAL(fp), PARAMETER, PUBLIC :: AIRMW = 28.9644_fp               !< Average molecular weight of dry air [g/mol]
   REAL(fp), PARAMETER, PUBLIC :: H2OMW = 18.016_fp                !< Molecular weight of water [g/mol]
   REAL(fp), PARAMETER, PUBLIC :: Rd   = 287.0_fp                  !< Gas constant for dry air [J/K/kg]
   REAL(fp), PARAMETER, PUBLIC :: Rdg0 = Rd / g0                   !< Gas constant for dry air divided by gravity
   REAL(fp), PARAMETER, PUBLIC :: Rv = 461.00_fp                   !< Gas constant for water vapor [J/K/kg]
   REAL(fp), PARAMETER, PUBLIC :: SCALE_HEIGHT = 7600.0_fp         !< Atmospheric scale height [m]
   REAL(fp), PARAMETER, PUBLIC :: VON_KARMAN = 0.41_fp             !< Von Karman's constant (dimensionless)
   REAL(fp), PARAMETER, PUBLIC :: ATM = 1.01325e+5_fp              !< Standard atmospheric pressure [Pa]
   REAL(fp), PARAMETER, PUBLIC :: XNUMOLAIR = AVO / ( AIRMW * 1.e-3_fp )  !< Molecules of dry air per kg dry air
   ! \}

   ! \name Mathematical Constants
   !! \brief Mathematical constants and conversion factors
   !! \{
   REAL(fp), PARAMETER, PUBLIC :: PI     = 3.14159265358979323_fp  !< Pi (dimensionless)
   REAL(fp), PARAMETER, PUBLIC :: PI_180 = PI / 180.0_fp           !< Radians per degree conversion factor
   REAL(fp), PARAMETER, PUBLIC :: E = 2.718281828459045235360287471352_fp  !< Euler's number (dimensionless)
   ! \}

   ! \name Chemistry-Specific Constants
   !! \brief Constants for atmospheric chemistry calculations
   !! \{
   REAL(fp), PARAMETER, PUBLIC :: CONSVAP = 6.1078e+03_fp / ( BOLTZ * 1e+7_fp ) !< Condensation vapor pressure factor
   REAL(fp), PARAMETER, PUBLIC :: RGASLATM = 8.2057e-2_fp          !< Gas constant in L⋅atm/(K⋅mol)
   REAL(fp), PARAMETER, PUBLIC :: MWCARB = 12.01e-3_fp             !< Molecular weight of carbon [kg/mol]
   ! \}

contains

   !> \brief Validate atmospheric constants for physical consistency
   !!
   !! This subroutine performs runtime validation of atmospheric constants
   !! to ensure physical consistency and catch any compilation issues.
   !!
   !! \param[out] rc Return code
   subroutine validate_atmospheric_constants(rc)
      use error_mod, only: CC_SUCCESS, ERROR_NUMERICAL_INSTABILITY
      implicit none
      integer, intent(out) :: rc

      real(fp) :: test_value
      real(fp), parameter :: TOLERANCE = 1.0e-12_fp

      rc = CC_SUCCESS

      ! Test fundamental relationships
      ! Ideal gas law consistency
      test_value = RSTARG / AIRMW * 1000.0_fp  ! Should equal Rd
      if (abs(test_value - Rd) > TOLERANCE) then
         rc = ERROR_NUMERICAL_INSTABILITY
         return
      endif

      ! Test that gravity is reasonable
      if (g0 < 9.0_fp .or. g0 > 10.0_fp) then
         rc = ERROR_NUMERICAL_INSTABILITY
         return
      endif

      ! Test Avogadro's number order of magnitude
      if (AVO < 6.0e23_fp .or. AVO > 7.0e23_fp) then
         rc = ERROR_NUMERICAL_INSTABILITY
         return
      endif

      ! Test derived constants
      test_value = 100.0_fp / g0
      if (abs(test_value - g0_100) > TOLERANCE) then
         rc = ERROR_NUMERICAL_INSTABILITY
         return
      endif

   end subroutine validate_atmospheric_constants

end module Constants
