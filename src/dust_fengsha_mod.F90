module dust_fengsha_mod
!
!  This module developed by Barry Baker (NOAA ARL)
!  For serious questions contact barry.baker@noaa.gov
!
!  07/16/2019 - Adapted for NUOPC/GOCART, R. Montuoro
!  02/01/2020 - Adapted for FV3/CCPP, Haiqin Li
!  06/2023, Restructure for CATChem, Jian.He@noaa.gov
!  08/2023, Update to latest version Barry.Baker@noaa.gov
  use catchem_constants, only : kind_chem, g=>con_g, pi=>con_pi
  use catchem_config, only: num_chem,num_emis_dust,&
                            p_dust_1,p_dust_2,p_dust_3,p_dust_4,p_dust_5, &
                            p_edust1,p_edust2,p_edust3,p_edust4,p_edust5
  use dust_data_mod

  implicit none

  private

  public :: gocart_dust_fengsha_driver

contains

  subroutine gocart_dust_fengsha_driver(dt,           &
          chem_arr,rho_phy,dz8w,smois,         &
          delp,ssm,isltyp,vegfra,snowh,area,  &
          emis_dust,ust,znt,clay,sand, &
          rdrag,uthr,num_soil_layers,random_factor)

    IMPLICIT NONE

    ! Input Variables
    ! ---------------
    INTEGER, INTENT(IN) :: isltyp
    INTEGER, INTENT(IN) :: num_soil_layers ! number of soil layers 
    REAL(kind=kind_chem), INTENT(IN) :: dt ! timestep 
    REAL(kind=kind_chem), INTENT(IN) :: rho_phy ! air density
    REAL(kind=kind_chem), INTENT(IN) :: dz8w 
    REAL(kind=kind_chem), INTENT(IN) :: delp    ! del Pressure in first layer
    REAL(kind=kind_chem), INTENT(IN) :: ssm     ! sediment supply map 
    REAL(kind=kind_chem), INTENT(IN) :: vegfra  ! vegetation fraction 
    REAL(kind=kind_chem), INTENT(IN) :: snowh   ! snow height 
    REAL(kind=kind_chem), INTENT(IN) :: area    ! area of grid cell 
    REAL(kind=kind_chem), INTENT(IN) :: ust     ! ustar 
    REAL(kind=kind_chem), INTENT(IN) :: znt     ! surface roughness length 
    REAL(kind=kind_chem), INTENT(IN) :: clay    ! fractional clay content 
    REAL(kind=kind_chem), INTENT(IN) :: sand    ! fractional sand content 
    REAL(kind=kind_chem), INTENT(IN) :: rdrag   ! drag partition 
    REAL(kind=kind_chem), INTENT(IN) :: uthr    ! dry threshold friction velocity
    REAL(kind=kind_chem), INTENT(IN) :: random_factor
  
    REAL(kind=kind_chem), DIMENSION( num_chem ), INTENT(INOUT)                 :: chem_arr  ! chemical array 
    REAL(kind=kind_chem), DIMENSION( num_emis_dust ), OPTIONAL, INTENT(INOUT ) :: emis_dust ! final dust emission 
    REAL(kind=kind_chem), DIMENSION( num_soil_layers ), INTENT(IN)             :: smois     ! volumetric soil moisture at 0-5cm 

    ! Local variables
    ! ---------------
    integer :: do_dust ! 0 - no dust emission 1 - dust emission 
    integer :: n       ! looping variable 
    
    real(kind=kind_chem), DIMENSION (num_emis_dust) :: tc ! tracer concentration 
    REAL(kind=kind_chem), DIMENSION (num_emis_dust) :: bin_emis ! bin emisions 
    real(kind=kind_chem)                            :: emis ! total emission
    real(kind=kind_chem), DIMENSION (num_emis_dust) :: distribution ! fractional distribution for size bins 
    
    real(kind_chem), dimension (3)  :: massfrac      ! fractional mass of sand silt and clay 
    
    real(kind=kind_chem), PARAMETER :: conver=1.e-9 ! parameter to convert units   
    real(kind=kind_chem), PARAMETER :: converi=1.e9 ! parameter to convert units back 

    ! Total concentration at lowest model level. 
    ! ------------------------------------------
    do n=0,num_emis_dust-1 
       tc(n+1)=chem_arr(p_dust_1+n)*conver
    end do

    ! Air mass at lowest model level.
    airmas=area * delp / g

    ! ====================================
    ! Don't do dust over certain criteria 
    ! ====================================
    do_dust = 1
    
    ! limit where there is lots of vegetation
    ! redundent with rdrag but keep because it is updated in NRT right now
    if (vegfra .gt. .2) then 
       do_dust = 0
    endif

    ! limit where there is snow on the ground
    if (snowh .gt. 0) then
       do_dust = 0
    endif

    ! Do not allow areas with bedrock, lava, or land-ice to loft
    IF (isltyp.eq. 15 .or. isltyp .eq. 16. .or. &
          isltyp .eq. 18) then
      do_dust=0
    ENDIF
    
    ! do not allow dust over the ocean 
    IF (isltyp .eq. 0)then
      do_dust=0
    endif
    
    ! check ssm input valid range
    if (ssm .lt. 0.05 .and. ssm .gt. 1.0)  then ! ensure values are realistic
       do_dust = 0
    end if
    
    ! check drag partition valid range 
    if (rdrag .lt. 0.05 .and. rdrag .gt. 1.0) then 
       do_dust = 0
    endif
    
    if(do_dust == 1 ) return
    ! ====================================

    ! Call fengsha dust emission for total emission 
    call DustEmissionFENGSHA(slc, clay, sand, ssm, rdrag, airdens, ustar, uthrs, area, dust_alpha, dust_gamma, emis)

    ! call dust distribution 
    call DustAerosolDistributionKok(reff_dust, lo_dust, up_dust, distribution)

    ! Distribute emissions to bins and convert to mass flux (kg s-1)
    ! --------------------------------------------------------------
    bin_emis = distribution * emis * dt * random_factor

    ! now convert to tracer concentration 
    ! -----------------------------------
    DO n=1,nmx
       dsrc = bin_emis(n)
       IF (dsrc < 0.0) dsrc = 0.0

       ! Update dust mixing ratio at first model level.
       tc(n) = tc(n) + dsrc / airmas ! (kg/kg)
       emis_dust(n)= dsrc/(dxy*dt1) ! diagnostic (kg/m2/s)
    END DO
    
    do n = 0, nmx-1
       chem_arr(p_dust_1+n))=tc(n + 1)*converi ! (ug/kg)
    end do

  end subroutine gocart_dust_fengsha_driver

  subroutine DustEmissionFENGSHA(slc, clay, ssm, rdrag, airdens, ustar, uthrs, area, ,alpha, gamma, emissions)

      ! !USES:
      implicit NONE

      ! !INPUT PARAMETERS:
      real(kind=kind_chem), intent(in) :: slc      ! liquid water content of soil layer, volumetric fraction [1]
      real(kind=kind_chem), intent(in) :: clay     ! fractional clay content [1] - range: [0 1]
      real(kind=kind_chem), intent(in) :: sand     ! fractional clay content [1] - range: [0 1]
      real(kind=kind_chem), intent(in) :: ssm      ! erosion map [1] - range: [0 1]
      real(kind=kind_chem), intent(in) :: rdrag    ! drag partition [1/m] - range: [0 1]
      real(kind=kind_chem), intent(in) :: airdens  ! air density at lowest level [kg/m^3]
      real(kind=kind_chem), intent(in) :: ustar    ! friction velocity [m/sec]
      real(kind=kind_chem), intent(in) :: uthrs    ! threshold velocity [m/2]
      real(kind=kind_chem), intent(in) :: kvhmax   ! max. vertical to horizontal mass flux ratio [1]
      real(kind=kind_chem), intent(in) :: area     ! area of dust emission (can be fractional area of grid cell)

      ! this will need to be read in by a configuration file
      ! right now I just hard coded it in dust_data_mod.F90
      !=====================================================
      real(kind=kind_chem), intent(in) :: alpha    ! scaling factor [1]
      real(kind=kind_chem), intent(in) :: gamma    ! scaling factor [1]
      !=====================================================

      ! !OUTPUT PARAMETERS:
      REAL(kind=kind_chem), dimension(:), intent(inout) :: total_emissions ! binned surface emissions [kg/(m^2 sec)]

      ! !DESCRIPTION: Compute dust emissions using NOAA/ARL FENGSHA model
      !
      ! !REVISION HISTORY:
      !
      ! 22Feb2020 B.Baker/NOAA    - Original implementation
      ! 29Mar2021 R.Montuoro/NOAA - Refactored for process library
      ! 09Aug2022 B.Baker/NOAA    - Adapted for CCPP-Physics

      ! !Local Variables
      real(kind=kind_chem) :: alpha_grav
      real(kind=kind_chem) :: h
      real(kind=kind_chem) :: kvh
      real(kind=kind_chem) :: q
      real(kind=kind_chem) :: rustar
      real(kind=kind_chem) :: u_sum
      real(kind=kind_chem) :: u_thresh
      real(kind=kind_chem) :: emission
      real(kind=kind_chem) :: stotal
      real(kind=kind_chem) :: dmass
      real(kind=kind_chem), dimension(nsalt) :: dsurface
      real(kind=kind_chem), dimension(3) :: massfrac ! fractional soil content of sand(1) silt(2) and clay(3)

      real(kind=kind_chem), parameter:: clay_thresh = 0.2
      real(kind=kind_chem), parameter :: rhow = 1000.

      !EOP
      !-------------------------------------------------------------------------
      !  Begin

      !  Initialize emissions
      !  --------------------
      total_emissions = 0.

      !  Prepare scaling factor
      !  ----------------------
      alpha_grav = dust_alpha / con_g

      ! Compute vertical-to-horizontal mass flux ratio
      ! ----------------------------------------------
      if (clay > clay_thresh) then
         kvh = kvhmax
      else
         kvh = 10.0**(13.4*clay-6.0)
      end if

      ! Compute total emissions
      ! -----------------------
      emission = alpha_grav * (ssm ** dust_gamma) * airdens * kvh

      !  Compute threshold wind friction velocity using drag partition
      !  -------------------------------------------------------------
      rustar = rdrag * ustar

      !  Now compute size-dependent total emission flux
      !  ----------------------------------------------
      ! Fecan moisture correction
      ! -------------------------
      if (DUST_OPT_FENGSHA_FECAN .eq. .true.) then
         !  Saturated volumetric water content (sand-dependent) ! [m3 m-3]
         vsat = 0.489 - 0.00126 * ( 100. * sand )

         !  Gravimetric soil content
         grvsoilm = vsoil * rhow / (dust_den * (1. - vsat))

         !   Compute fecan dry limit
         drylimit = clay * (14.0 * clay + 17.0)

         !  Compute soil moisture correction
         h = sqrt(1.0 + 1.21 * max(0., grvsoilm - drylimit)**0.68)
      else
         ! Shao soil mositure
         !--------------------
         if (slc <= 0.03) then
            h = exp(22.7 * slc)
         else
            h = exp(93.5 * slc - 2.029)
         end if
      end if

      ! Adjust threshold
      ! ----------------
      u_thresh = uthrs * h

      u_sum = rustar + u_thresh

      ! Compute Horizontal Saltation Flux according to Eq (9) in Webb et al. (2020)
      ! ---------------------------------------------------------------------------
      q = max(0., rustar - u_thresh) * u_sum * u_sum

      ! Calculate total dust using the dust potential (q) 
      ! -------------------------------------------------
      emission = emission * q 
      
      ! Calculate saltation surface area distribution from sand, silt, and clay
      ! mass fractions and saltation bin fraction. Based on Eqn. (32) in 
      ! Marticorena & Bergametti, 1995 (hereon, MB95).
      ! ----------------------------------------------
      DO n=1,nsalt
         dmass=massfrac(spoint(n))*frac_salt(n)
         dsurface(n)=0.75*dmass/(den_salt(n)*reff_salt(n))
      ENDDO

      ! The following equation yields relative surface area fraction.  
      ! ------------------------------------------------------------
      stotal=SUM(dsurface(:))
      DO n=1,nsalt
         total_emissions = total_emissions + emission * dsurface(n)/stotal
      ENDDO
      
      ! Distribute emissions to bins and convert to mass flux (kg s-1)
      ! --------------------------------------------------------------
      !emissions = distribution * total_emissions * q


   end subroutine DustEmissionFENGSHA
   !-----------------------------------------------------------------
   subroutine DustAerosolDistributionKok ( radius, rLow, rUp, distribution )
     
     ! !USES:
     implicit NONE
     
     ! !INPUT PARAMETERS:
     real, dimension(:), intent(in)  :: radius      ! Dry particle bin effective radius [um]
     real, dimension(:), intent(in)  :: rLow, rUp   ! Dry particle bin edge radii [um]
     
     ! !OUTPUT PARAMETERS:
     real, dimension(:), intent(out) :: distribution    ! Normalized dust aerosol distribution [1]
     
     ! !DESCRIPTION: Computes lognormal aerosol size distribution for dust bins according to
     !               J.F.Kok, PNAS, Jan 2011, 108 (3) 1016-1021; doi:10.1073/pnas.1014798108
     !
     ! !REVISION HISTORY:
     !
     ! 22Feb2020 B.Baker/NOAA    - Original implementation
     ! 01Apr2021 R.Montuoro/NOAA - Refactored for GOCART process library
     !
     
     ! !Local Variables
     integer :: n, nbins
     real    :: diameter, dlam, dvol
     
     !   !CONSTANTS
     real, parameter    :: mmd    = 3.4          ! median mass diameter [um]
     real, parameter    :: stddev = 3.0          ! geometric standard deviation [1]
     real, parameter    :: lambda = 12.0         ! crack propagation length [um]
     real, parameter    :: factor = 1.e0 / (sqrt(2.e0) * log(stddev))  ! auxiliary constant
     
     character(len=*), parameter :: myname = 'DustAerosolDistributionKok'
     
     !EOP
     !-------------------------------------------------------------------------
     !  Begin...
     
     distribution = 0.
     
     !  Assume all arrays are dimensioned consistently
     nbins = size(radius)
     
     dvol = 0.
     do n = 1, nbins
        diameter = 2 * radius(n)
        dlam = diameter/lambda
        distribution(n) = diameter * (1. + erf(factor * log(diameter/mmd))) * exp(-dlam * dlam * dlam) * log(rUp(n)/rLow(n))
        dvol = dvol + distribution(n)
     end do
     
     !  Normalize distribution
     do n = 1, nbins
        distribution(n) = distribution(n) / dvol
     end do
     
   end subroutine DustAerosolDistributionKok
   
end module dust_fengsha_mod
