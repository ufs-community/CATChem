! Jian.He@noaa.gov, 05/2023
! Move to parameters folder for CATChem

module catchem_constants

   implicit none

   public

   integer, parameter :: kind_chem = 8

   real(kind=kind_chem),parameter:: con_g      =9.80665e+0_kind_chem                !< gravity (\f$m/s^{2}\f$)
   real(kind=kind_chem),parameter:: con_rd     =2.8705e+2_kind_chem                 !< gas constant air (\f$J/kg/K\f$)
   real(kind=kind_chem),parameter:: con_rv     =4.6150e+2_kind_chem                 !< gas constant H2O (\f$J/kg/K\f$)
   real(kind=kind_chem),parameter:: con_cp     =1.0046e+3_kind_chem                 !< spec heat air at p (\f$J/kg/K\f$)
   real(kind=kind_chem),parameter:: con_cv     =7.1760e+2_kind_chem                 !< spec heat air at v (\f$J/kg/K\f$)
   real(kind=kind_chem),parameter:: con_pi     =4.0d0*atan(1.0d0)                   !< pi

end module catchem_constants
