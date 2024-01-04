module ddep_data_mod

   use catchem_constants , only : kind_chem

   implicit none

!--------------------------------------------------
! many of these parameters will depend on the RADM mechanism!
! if you change it, lets talk about it and get it done!!!
!--------------------------------------------------

   INTEGER, PARAMETER :: dep_seasons = 5
   INTEGER, PARAMETER :: nlu = 25
   REAL, parameter    :: small_value = 1.e-36
   REAL, parameter    :: large_value = 1.e36

!--------------------------------------------------
! following currently hardwired to USGS
!--------------------------------------------------
   integer, parameter :: isice_temp   = 24
   integer, parameter :: iswater_temp = 16
   integer, parameter :: wrf2mz_lt_map(nlu) = (/ 1, 2, 2, 2, 2, &
      4, 3, 3, 3, 3, &
      4, 5, 4, 5, 6, &
      7, 9, 6, 8, 9, &
      6, 6, 8, 0, 0 /)
   real, parameter    :: wh2o = 18.0153
   real, parameter    :: wpan = 121.04793
   real, PARAMETER ::  KARMAN=0.4
   INTEGER,  parameter :: luse2usgs(21) = (/14,13,12,11,15,8,9,10,10,7, &
      17,4,1,5,24,19,16,21,22,23,16 /)
   character(len=4), parameter :: mminlu = 'USGS'

   INTEGER :: month = 0
   INTEGER :: ixxxlu(nlu)
!     include modis landuse
!--

   REAL    :: kpart(nlu)
   REAL    :: rac(nlu,dep_seasons), rclo(nlu,dep_seasons), rcls(nlu,dep_seasons)
   REAL    :: rgso(nlu,dep_seasons), rgss(nlu,dep_seasons)
   REAL    :: ri(nlu,dep_seasons), rlu(nlu,dep_seasons)
   REAL    :: ri_pan(5,11)
   real    :: c0_pan(11) = (/ 0.000, 0.006, 0.002, 0.009, 0.015, &
      0.006, 0.000, 0.000, 0.000, 0.002, 0.002 /)
   real    :: k_pan (11) = (/ 0.000, 0.010, 0.005, 0.004, 0.003, &
      0.005, 0.000, 0.000, 0.000, 0.075, 0.002 /)

   public

end module ddep_data_mod
