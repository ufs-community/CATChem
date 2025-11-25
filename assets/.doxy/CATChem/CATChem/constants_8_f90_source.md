

# File constants.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**constants.F90**](constants_8_f90.md)

[Go to the documentation of this file](constants_8_f90.md)


```Fortran

module constants
   use precision_mod

   implicit none
   private

   ! \name Fundamental Physical Constants (must be defined first for dependencies)
   !! \brief Universal physical constants
   !! \{
   REAL(fp), PARAMETER, PUBLIC :: AVO = 6.022140857e+23_fp         
   REAL(fp), PARAMETER, PUBLIC :: g0     = 9.80665e+0_fp           
   REAL(fp), PARAMETER, PUBLIC :: g0_100 = 100.0_fp / g0           
   REAL(fp), PARAMETER, PUBLIC :: Re = 6.3710072e+6_fp             
   REAL(fp), PARAMETER, PUBLIC :: RSTARG = 8.3144598_fp            
   REAL(fp), PARAMETER, PUBLIC :: BOLTZ = 1.38064852e-23_fp        
   REAL(fp), PARAMETER, PUBLIC :: PLANCK = 6.62606957e-34_fp       
   REAL(fp), PARAMETER, PUBLIC :: CCONST = 2.99792458e+8_fp        
   ! \}

   ! \name Atmospheric Properties
   !! \brief Constants related to atmospheric composition and properties
   !! \{
   REAL(fp), PARAMETER, PUBLIC :: Cp = 1.0046e+3_fp                
   REAL(fp), PARAMETER, PUBLIC :: Cv = 7.1760e+2_fp                
   REAL(fp), PARAMETER, PUBLIC :: AIRMW = 28.9644_fp               
   REAL(fp), PARAMETER, PUBLIC :: H2OMW = 18.016_fp                
   REAL(fp), PARAMETER, PUBLIC :: Rd   = 287.0_fp                  
   REAL(fp), PARAMETER, PUBLIC :: Rdg0 = rd / g0                   
   REAL(fp), PARAMETER, PUBLIC :: Rv = 461.00_fp                   
   REAL(fp), PARAMETER, PUBLIC :: SCALE_HEIGHT = 7600.0_fp         
   REAL(fp), PARAMETER, PUBLIC :: VON_KARMAN = 0.41_fp             
   REAL(fp), PARAMETER, PUBLIC :: ATM = 1.01325e+5_fp              
   REAL(fp), PARAMETER, PUBLIC :: XNUMOLAIR = avo / ( airmw * 1.e-3_fp )  
   ! \}

   ! \name Mathematical Constants
   !! \brief Mathematical constants and conversion factors
   !! \{
   REAL(fp), PARAMETER, PUBLIC :: PI     = 3.14159265358979323_fp  
   REAL(fp), PARAMETER, PUBLIC :: PI_180 = pi / 180.0_fp           
   REAL(fp), PARAMETER, PUBLIC :: E = 2.718281828459045235360287471352_fp  
   ! \}

   ! \name Chemistry-Specific Constants
   !! \brief Constants for atmospheric chemistry calculations
   !! \{
   REAL(fp), PARAMETER, PUBLIC :: CONSVAP = 6.1078e+03_fp / ( boltz * 1e+7_fp ) 
   REAL(fp), PARAMETER, PUBLIC :: RGASLATM = 8.2057e-2_fp          
   REAL(fp), PARAMETER, PUBLIC :: MWCARB = 12.01e-3_fp             
   ! \}

contains

   subroutine validate_atmospheric_constants(rc)
      use error_mod, only: cc_success, error_numerical_instability
      implicit none
      integer, intent(out) :: rc

      real(fp) :: test_value
      real(fp), parameter :: TOLERANCE = 1.0e-12_fp

      rc = cc_success

      ! Test fundamental relationships
      ! Ideal gas law consistency
      test_value = rstarg / airmw * 1000.0_fp  ! Should equal Rd
      if (abs(test_value - rd) > tolerance) then
         rc = error_numerical_instability
         return
      endif

      ! Test that gravity is reasonable
      if (g0 < 9.0_fp .or. g0 > 10.0_fp) then
         rc = error_numerical_instability
         return
      endif

      ! Test Avogadro's number order of magnitude
      if (avo < 6.0e23_fp .or. avo > 7.0e23_fp) then
         rc = error_numerical_instability
         return
      endif

      ! Test derived constants
      test_value = 100.0_fp / g0
      if (abs(test_value - g0_100) > tolerance) then
         rc = error_numerical_instability
         return
      endif

   end subroutine validate_atmospheric_constants

end module constants
```


