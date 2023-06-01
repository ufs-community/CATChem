!Jian.He@noaa.gov, 05/2023
!Restructure for CATChem

module seas_data_mod

  use catchem_constants , only : kind_chem

  ! -- parameters from NGAC v2.4.0 (rev. d48932c)
  integer,                             parameter :: number_ss_bins  = 5
  ! -- lower/upper particle radii (um) for each bin
  real(kind=kind_chem), dimension(number_ss_bins), parameter :: ra = (/  0.03,   0.1,   0.5,   1.5,   5.0 /)
  real(kind=kind_chem), dimension(number_ss_bins), parameter :: rb = (/   0.1,   0.5,   1.5,   5.0,  10.0 /)
  ! -- global scaling factors for sea salt emissions (originally 0.875 in NGAC namelist)
  !real(kind=kind_chem), dimension(number_ss_bins), parameter :: emission_scale = (/ 0.100, 0.100, 0.100, 0.100, 0.100 /)
  real(kind=kind_chem), dimension(number_ss_bins), parameter :: emission_scale = (/ 1.0, 1.0, 1.0, 1.0, 1.0 /)
  ! -- sea salt density
  real(kind=kind_chem), dimension(number_ss_bins), parameter :: den_seas  = (/    2200.,    2200.,    2200.,    2200.,    2200. /)
  ! -- particle effective radius (m)
  real(kind=kind_chem), dimension(number_ss_bins), parameter :: reff_seas = (/ 0.079e-6, 0.316e-6, 1.119e-6, 2.818e-6, 7.772e-6 /)

  ! -- NGAC parameters
  integer, parameter :: emission_scheme = 3    ! GEOSS 2012
  real(kind=kind_chem), parameter    :: r80fac = 1.65     ! ratio of radius(RH=0.8)/radius(RH=0.) [Gerber]
  real(kind=kind_chem), parameter    :: rhop = 2200.      ! dry seasalt density [kg m-3]

end module seas_data_mod
