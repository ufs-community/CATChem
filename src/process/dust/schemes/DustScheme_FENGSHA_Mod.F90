!> \file DustScheme_FENGSHA_Mod.F90
!! \brief Fengsha Dust emission scheme developed at NOAA ARL for use at NOAA NWS
!!
!! Pure science kernel for fengsha scheme in dust process.
!! This module contains ONLY the computational algorithm with NO infrastructure dependencies.
!! Uses only basic Fortran types for maximum portability and reusability.
!!
!! SCIENCE CUSTOMIZATION GUIDE:
!! 1. Modify the algorithm in compute_fengsha (search for "TODO")
!! 2. Add scheme-specific helper subroutines as needed
!! 3. Update physical constants for your scheme
!! 4. Customize the environmental response functions
!!
!! INFRASTRUCTURE RESPONSIBILITIES (handled by host model):
!! - Parameter initialization and validation
!! - Input array validation and error handling
!! - Memory management and array allocation
!! - Integration with host model time stepping
!!
!! Generated on: 2026-04-17T13:57:10.251187
!! Author: Barry Baker & Wei Li
!! Reference: Zhang et al. 2022
module DustScheme_FENGSHA_Mod

   use precision_mod, only: fp
   use DustCommon_Mod, only: DustSchemeFENGSHAConfig

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_fengsha

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter    :: SSM_THRESH  = 1.0E-02_fp  ! Minimum erodibility threshold
   real(fp), parameter    :: VEG_THRESH  = 0.4_fp      ! Maximum vegetation threshold
   real(fp), parameter    :: SMALL       = 1.0E-10_fp  ! Small number for division protection
   real(fp), parameter    :: MAX_RDRAG   = 0.3_fp      ! Maximum drag partition ratio
   real(fp), parameter    :: CLAY_THRESH = 0.2_fp      ! clay fraction above which the maximum flux ratio is returned

contains

   !> Pure science computation for fengsha scheme
   !!
   !! This is a pure computational kernel implementing Fengsha Dust emission scheme developed at NOAA ARL for use at NOAA NWS.
   !! NO error checking, validation, or infrastructure concerns.
   !! Host model must ensure all inputs are valid before calling.
   !!
   !! @param[in]  num_layers     Number of vertical layers
   !! @param[in]  num_species    Number of chemical species
   !! @param[in]  params         Scheme parameters (pre-validated by host)
   !! @param[in]  g0    Required constant from Constants module
   !! @param[in]  airden    AIRDEN field [appropriate units]
   !! @param[in]  clayfrac    CLAYFRAC field [appropriate units]
   !! @param[in]  frlake    FRLAKE field [appropriate units]
   !! @param[in]  frsno    FRSNO field [appropriate units]
   !! @param[in]  gvf    GVF field [appropriate units]
   !! @param[in]  lai    LAI field [appropriate units]
   !! @param[in]  lwi    LWI field [appropriate units]
   !! @param[in]  rdrag    RDRAG field [appropriate units]
   !! @param[in]  sandfrac    SANDFRAC field [appropriate units]
   !! @param[in]  soilm    SOILM field [appropriate units]
   !! @param[in]  ssm    SSM field [appropriate units]
   !! @param[in]  tskin    TSKIN field [appropriate units]
   !! @param[in]  ustar    USTAR field [appropriate units]
   !! @param[in]  ustar_threshold    USTAR_THRESHOLD field [appropriate units]
   !! @param[in]  z0    Z0 field [appropriate units]
   !! @param[in]  species_radius    Species radius property
   !! @param[in]  species_lower_radius    Species lower_radius property
   !! @param[in]  species_upper_radius    Species upper_radius property
   !! @param[in]  species_conc   Species concentrations [ppm or ug/kg] (num_layers, num_species)
   !! @param[inout] species_tendencies  Species tendency terms [mol/mol/s] (num_layers, num_species)
   !! @param[inout] dust_emission_total    Total dust emissions for all bins [kg/m2/s]
   !! @param[inout] dust_emission_per_bin    Dust emission flux per bin [kg/m2/s] (num_species)
   !! @param[inout] dust_horizontal_flux    Total horizontal flux - Q [kg/m2/s]
   !! @param[inout] dust_moisture_correction    Moisture Correction - H [1.0]
   !! @param[inout] dust_effective_threshold    Effective Dust threshold friction velocity: u_thres * H / R [m/s]
   !! @param[in] diagnostic_species_id Indices mapping diagnostic species to species array (optional, for per-species diagnostics)
   subroutine compute_fengsha( &
      num_layers, &
      num_species, &
      params, &
      g0, &
      airden, &
      clayfrac, &
      frlake, &
      frsno, &
      gvf, &
      lai, &
      lwi, &
      rdrag, &
      sandfrac, &
      soilm, &
      ssm, &
      tskin, &
      ustar, &
      ustar_threshold, &
      z0, &
      species_radius, &
      species_lower_radius, &
      species_upper_radius, &
      species_conc, &
      species_tendencies, &
      dust_emission_total, &
      dust_emission_per_bin, &
      dust_horizontal_flux, &
      dust_moisture_correction, &
      dust_effective_threshold, &
      diagnostic_species_id &
      )

      ! Arguments
      integer, intent(in) :: num_layers
      integer, intent(in) :: num_species
      type(DustSchemeFENGSHAConfig), intent(in) :: params
      real(fp), intent(in) :: g0  ! Required constant from Constants module
      real(fp), intent(in) :: airden(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: clayfrac  ! Surface field - scalar
      real(fp), intent(in) :: frlake  ! Surface field - scalar
      real(fp), intent(in) :: frsno  ! Surface field - scalar
      real(fp), intent(in) :: gvf  ! Surface field - scalar
      real(fp), intent(in) :: lai  ! Surface field - scalar
      integer, intent(in) :: lwi  ! Surface field - scalar
      real(fp), intent(in) :: rdrag  ! Surface field - scalar
      real(fp), intent(in) :: sandfrac  ! Surface field - scalar
      real(fp), intent(in) :: soilm(:)  ! variable dimension array
      real(fp), intent(in) :: ssm  ! Surface field - scalar
      real(fp), intent(in) :: tskin  ! Surface field - scalar
      real(fp), intent(in) :: ustar  ! Surface field - scalar
      real(fp), intent(in) :: ustar_threshold  ! Surface field - scalar
      real(fp), intent(in) :: z0  ! Surface field - scalar
      real(fp), intent(in) :: species_radius(:)  ! Species radius property
      real(fp), intent(in) :: species_lower_radius(:)  ! Species lower_radius property
      real(fp), intent(in) :: species_upper_radius(:)  ! Species upper_radius property
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(inout) :: species_tendencies(num_layers, num_species)
      real(fp), intent(inout), optional :: dust_emission_total
      real(fp), intent(inout), optional :: dust_emission_per_bin(:)
      real(fp), intent(inout), optional :: dust_horizontal_flux
      real(fp), intent(inout), optional :: dust_moisture_correction
      real(fp), intent(inout), optional :: dust_effective_threshold
      integer, intent(in), optional :: diagnostic_species_id(:)  ! Indices mapping diagnostic species to species array

      ! Local variables
      logical :: skip
      integer :: k, species_idx
      integer :: diag_idx  ! For diagnostic species indexing
      real(fp) :: fracland, total_emissions
      real(fp) :: R                              !< Drag Partition [1]
      real(fp) :: H                              !< Soil Moisture Attenuation Factor
      real(fp) :: alpha_grav                     !< Alpha Parameter over Gravity
      real(fp) :: q                              !< Horizontal Mass Flux
      real(fp) :: FengshaScale                   !< Total Scaling Factor
      real(fp) :: h_to_v_ratio                   !< Horizontal to Vertical Mass Flux Ratio
      real(fp) :: distribution(num_species)      !< Distribution Weights

      !needs to reinitialize otherwise the skip condition below will cause weird maps.
      if (present(dust_effective_threshold)) dust_effective_threshold = 0.0_fp
      if (present(dust_horizontal_flux)) dust_horizontal_flux = 0.0_fp
      if (present(dust_moisture_correction)) dust_moisture_correction = 0.0_fp
      if (present(dust_emission_total)) dust_emission_total = 0.0_fp
      if (present(dust_emission_per_bin)) dust_emission_per_bin = 0.0_fp

      ! Precompute scaling factor
      alpha_grav = params%alpha / max(g0, SMALL)

      ! Skip criteria evaluation
      skip = .false.
      skip = (LWI /= 1)  !land = 1, water = 0, ice = 2

      select case(params%drag_option)
       case(2)  ! Darmenova scheme
         if (.not. skip) then
            skip = (gvf < 0.0_fp) .or. (gvf >= VEG_THRESH) .or. &
               (rdrag > MAX_RDRAG)
         endif
       case(3)  ! Leung scheme
         if (.not. skip) then
            skip = (gvf < 0.0_fp) .or. (lai >= VEG_THRESH)
         endif
       case default
         if (.not. skip) skip = (rdrag < 0.0_fp .or. rdrag > 1.0_fp)
      end select

      if (.not. skip) then
         skip = (clayfrac /= clayfrac) .or. (sandfrac /= sandfrac) ! check for NaNs
         if (skip) return !return here to avoid floating point checking below.
      endif

      if (.not. skip) then
         skip = (SSM < SSM_THRESH) .or. &
            (clayfrac < 0.0_fp) .or. (sandfrac < 0.0_fp) .or. &
            (clayfrac > 1.0_fp) .or. (sandfrac > 1.0_fp)
      endif

      ! Don't do dust over frozen soil
      !--------------------------------
      if (TSKIN <= 273.15_fp) then
         ! skip = .true.
      endif

      ! Skip computation if criteria not met
      if (skip) then
         return
      end if

      ! Calculate land fraction (TODO: I am using 1 - frlake - frsno, not following GOCART below)
      ! fracland = max(0.0_fp, min(1.0_fp, 1.0_fp - frlake)) * &
      !    max(0.0_fp, min(1.0_fp, 1.0_fp - frsno))

      fracland = max(0.0_fp, min(1.0_fp, 1.0_fp - frsno - frlake))  ! my calculation

      ! Compute vertical-to-horizontal mass flux ratio
      ! B.Marticorena, G.Bergametti, J.Geophys.Res., 1995
      ! doi:10.1029/95JD00690
      ! ----------------------------------------------
      if (CLAYFRAC > clay_thresh) then
         h_to_v_ratio = params%kvhmax
      else
         h_to_v_ratio = 10.0_fp**(13.4_fp*CLAYFRAC-6.0_fp)
      end if

      ! Compute the Drag Partition
      ! 1: Input Drag Partition
      ! 2: Darmenova 2009
      ! 3: Leung 2022
      ! 4: MB95 Drag Partition
      !----------------------------
      select case(params%drag_option)
       case(1)
         R = rdrag
       case(2)
         call DarmenovaDragPartition(rdrag, gvf, VEG_THRESH, R)
       case(3)
         call LeungDragPartition(rdrag, lai, gvf, VEG_THRESH, R)
       case(4)
         call MB95_DragPartition(z0, R)
      end select

      ! compute moisture correction factor
      select case(params%moist_option)
       case(1)
         call Fecan_SoilMoisture(clayfrac, sandfrac, soilm(1) * params%moist_correction_factor, params%drylimit_factor, h)
       case(2)
         call Zhao_SoilMoisture(soilm(1), h)
      end select

      ! Compute the Horizontal Mass Flux
      ! 1: White 1979 (in GOCART2G version of Fengsha)
      ! 2: Draxler 2001
      ! 3: Kawamura 1964 / Webb 2020
      !----------------------------------
      select case (params%horizflux_option)
       case(1)
         call White_HorizFlux(ustar, ustar_threshold, R, h, q)
       case(2)
         call Draxler_HorizFlux(ustar, ustar_threshold, R, h, q)
       case(3)
         call Kawamura_HorizFlux(ustar, ustar_threshold, R, h, q)
      end select

      ! Calculate total emissions potential
      FengshaScale = alpha_grav * fracland * (ssm ** params%gamma) * airden(1)
      total_emissions = FengshaScale * h_to_v_ratio * q

      !debug only
      ! if (total_emissions > 1.0e-5_fp) then
      !    write(*,'(A,F12.8)') 'Debug: Total Emissions = ', total_emissions
      !    write(*,'(A,F12.8)') 'Debug: Total Fengsha Scale = ', FengshaScale
      !    write(*,'(A,F12.8)') 'Debug: h_to_v_ratio = ', h_to_v_ratio
      !    write(*,'(A,F12.8)') 'Debug: q = ', q
      !    write(*,'(A,F12.8)') 'Debug: ustar = ', ustar
      !    write(*,'(A,F12.8)') 'Debug: ustar_threshold = ', ustar_threshold
      !    write(*,'(A,F12.8)') 'Debug: h = ', h
      !    write(*,'(A,F12.8)') 'Debug: R = ', R
      !    write(*,'(A,F12.8)') 'Debug: clayfrac = ', clayfrac
      !    write(*,'(A,F12.8)') 'Debug: sandfrac = ', sandfrac
      !    write(*,'(A,F12.8)') 'Debug: soilm = ', soilm(1)
      !    write(*,'(A,F12.8)') 'Debug: LAI = ', LAI
      !    write(*,'(A,F12.8)') 'Debug: fracland = ', fracland
      !    write(*,'(A,F12.8)') 'Debug: airden = ', airden(1)
      !    write(*,'(A,F12.8)') 'Debug: ssm = ', ssm
      !    write(*,'(A,F12.8)') 'Debug: alpha_grav = ', alpha_grav
      ! end if


      ! get distribution of dust and map total emissions to species bins
      !--------------------------------
      select case (params%distribution_option)
       case(1)
         call KokDistribution(species_radius, species_lower_radius, species_upper_radius, distribution)
         !case(2) !not implemented yet
         !   call MengDistribution(species_radius, species_lower_radius, species_upper_radius, distribution)
      end select

      ! Main computation loop - CUSTOMIZE THIS SECTION FOR YOUR SCHEME
      do k = 1, num_layers

         ! Apply to each species
         do species_idx = 1, num_species

            species_tendencies(k, species_idx) = total_emissions * distribution(species_idx)

            ! Ensure non-negative emissions
            species_tendencies(k, species_idx) = max(0.0_fp, species_tendencies(k, species_idx))

            ! TODO: Update diagnostic fields here based on your scheme's requirements
            ! Each process should implement custom diagnostic calculations
            ! Example patterns:
            if (present(dust_emission_total)) then
               ! Add your custom total dust emissions for all bins calculation
               dust_emission_total = dust_emission_total + species_tendencies(k, species_idx)
            end if
            ! Per-species diagnostic: only update for diagnostic species
            if (present(dust_emission_per_bin) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == species_idx) then
                     ! Add your custom dust emission flux per bin calculation
                     dust_emission_per_bin(diag_idx) = species_tendencies(k, species_idx)
                     exit
                  end if
               end do
            end if
         end do ! species loop
      end do  ! layer loop

      ! save other species independent diagnostics
      if (present(dust_horizontal_flux)) then
         ! Add your custom total horizontal flux - q calculation
         dust_horizontal_flux = q
      end if
      if (present(dust_moisture_correction)) then
         ! Add your custom moisture correction - h calculation
         dust_moisture_correction = H
      end if
      if (present(dust_effective_threshold)) then
         ! Add your custom effective dust threshold friction velocity: u_thres * h / r calculation
         dust_effective_threshold = ustar_threshold * H / R
      end if

   end subroutine compute_fengsha

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines
   ! Examples: environmental response functions, species-specific calculations, etc.

   !>
   !!
   !! Kok, J. F. (2011a), A scaling theory for the size distribution of emitted
   !! dust aerosols suggests climate models underestimate the size of the global
   !! dust cycle, Proc. Natl. Acad. Sci. U. S. A., 108(3), 1016–1021,
   !! doi:10.1073/pnas.1014798108.
   !!
   !! \param radius Radius
   !! \param rLow Lower radius
   !! \param rUp Upper radius
   !! \param dist Distribution
   !!
   !!!>
   subroutine KokDistribution(radius, rLow, rUp, dist)
      !use constants, only: pi
      IMPLICIT NONE
      ! Parameters
      real(fp), dimension(:), intent(in)  :: radius
      real(fp), dimension(:), intent(in)  :: rLow
      real(fp), dimension(:), intent(in)  :: rUp
      real(fp), dimension(:), intent(out) :: dist

      ! Local Variables
      integer :: n          !< looping variable
      integer :: nbins      !< number of bins
      real(fp) :: diameter  !< effective diameter of particle
      real(fp) :: dlam      !< diameter / lambda
      real(fp) :: dvol      !< volume of particle

      ! Constants
      real(fp), parameter :: mmd = 3.4_fp                               !< median mass diameter [microns]
      real(fp), parameter :: stddev = 3.0_fp                            !< standard deviation [microns]
      real(fp), parameter :: lambda = 12.0_fp                           !< crack propagation length [um]
      real(fp), parameter :: factor = 1.0_fp / ( sqrt(2.0_fp) * log(stddev)) !< auxiliary constant for the distribution

      ! Initialize
      dvol = 0.0_fp
      dist = 0.0_fp
      nbins = size(radius)

      do n = 1, nbins
         diameter = radius(n) * 2.0_fp
         dlam = diameter / lambda
         dist(n) = diameter * (1._fp + erf(factor * log(diameter/mmd))) * exp(-dlam * dlam * dlam) * log(rUp(n)/rLow(n))
         dvol = dvol + dist(n)
      end do

      ! Normalize Distribution
      do n = 1, nbins
         dist(n) = dist(n) / dvol
      end do

   end subroutine KokDistribution

   !>
   !! Computes the soil moisture attenuation factor for dust emission
   !!
   !! Fecan, F., Marticorena, B., and Bergametti, G.: Parametrization of the increase of the aeolian
   !! erosion threshold wind friction velocity due to soil moisture for arid and semi-arid areas,
   !! Ann. Geophys., 17, 149–157, https://doi.org/10.1007/s00585-999-0149-7, 1999.
   !!
   !! \param clay Fractional clay content
   !! \param sand Fractional sand content
   !! \param volumetric_soil_moisture Volumetric soil moisture
   !! \param H Soil moisture attenuation factor for dust emission
   !!
   !!!>
   subroutine Fecan_SoilMoisture( clay, sand, volumetric_soil_moisture, b, H)
      IMPLICIT NONE
      ! Parameters
      !-----------
      real(fp), intent(in)  :: clay                      !< Fractional Clay Content
      real(fp), intent(in)  :: sand                      !< Fractional Sand Content
      real(fp), intent(in)  :: volumetric_soil_moisture  !< volumetric soil moisture fraction [m3 / m3]
      real(fp), intent(in)  :: b                         ! drylimit factor from Zender 2003
      real(fp), intent(out) :: H                         !< Soil Moisture attenuation factor for dust emission [1]

      ! Local Variables
      !----------------
      real(fp) :: vsat                      !< Saturated volumetric water content (sand-dependent) [m3 m-3]

      real(fp) :: gravimetric_soil_moisture !< Gravimetric soil moisture [kg/kg]
      real(fp) :: DryLimit                  !< Dry limit of the soil moisture [kg/kg]

      !CONSTANTS:
      real(fp), parameter :: waterDensity = 1000.0_fp    ! density of water [kg m-3]
      real(fp), parameter :: particleDensity = 1700.0_fp ! density of soil particles [kg m-3]

      ! Initialize
      !-----------
      H = 0.0_fp

      !--------------------------------------------
      ! Compute Saturated Volumetric Water Content
      !--------------------------------------------
      vsat = 0.489_fp - 0.00126_fp * (100._fp * sand)

      !--------------------------------------------
      ! Compute Gravimetric Soil moisture
      !--------------------------------------------
      gravimetric_soil_moisture = 100.0_fp * volumetric_soil_moisture * waterDensity / (particleDensity * (1.0_fp - vsat))

      !--------------------------------------------
      ! Compute Dry Limit
      !--------------------------------------------
      DryLimit = b * clay * (14.0_fp * clay + 17.0_fp)

      !--------------------------------------------
      ! Compute attenuation factor
      !--------------------------------------------
      H = sqrt(1.0_fp + 1.21_fp * max(0._fp, gravimetric_soil_moisture - DryLimit)**0.68_fp)

   end subroutine Fecan_SoilMoisture

   !>
   !! \brief Computes the soil moisture attenuation factor for dust emission
   !!
   !! Zhao, T. L., S. L. Gong, X. Y. Zhang, A. Abdel-Mawgoud, and Y. P. Shao (2006),
   !! An assessment of dust emission schemes in modeling east Asian dust storms,
   !! J. Geophys. Res., 111, D05S90, doi:10.1029/2004JD005746.
   !!
   !! \param volumetric_soil_moisture Volumetric soil moisture
   !! \param H Soil moisture attenuation factor for dust emission
   !!
   !! \ingroup catchem_dust_process
   !!!>
   subroutine Zhao_SoilMoisture( volumetric_soil_moisture, H)
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)  :: volumetric_soil_moisture  !< Volumetric soil moisture [m3 m-3]
      real(fp), intent(out) :: H                         !< Soil Moisture attenuation factor for dust emission [1]

      ! Initialize
      H = 0.0_fp

      !--------------------------------------------
      ! Compute attenuation factor
      !--------------------------------------------
      if (volumetric_soil_moisture <= 0.03_fp) THEN
         H = exp(22.7_fp * volumetric_soil_moisture)
      else
         H = exp(93.5_fp * volumetric_soil_moisture - 2.029_fp)
      endif

      return

   end subroutine Zhao_SoilMoisture


   function calc_drag_partition(sig, m, Beta, Lc) result(feff)
      real(fp), intent(in) :: sig, m, Beta, Lc
      real(fp) :: feff
      real(fp) :: R1, R2

      R1 = 1.0_fp / sqrt(1.0_fp - sig * m * Lc)
      R2 = 1.0_fp / sqrt(1.0_fp + m * Beta * Lc)
      feff = R1 * R2
   end function calc_drag_partition

   !>
   !! Calculates the double drag partition  from Darmenova et al. 2009
   !! DESCRIPTION: Computes the drag partition according to
   !!              Darmenova, K. et al. 2009 Dust emission parameterization scheme
   !!              regions in Central and East Asia, JGR Atmospheres, 114, D14201
   !
   ! !REVISION HISTORY:
   ! 27Jun2024 B.Baker/NOAA    - Original implementation
   ! DD MMM YYYY Author  - Refactored for improved structure
   !
   subroutine DarmenovaDragPartition(Lc, vegfrac, thresh, dragpartition)

      ! !USES:
      implicit NONE

      ! !INPUT PARAMETERS:
      real(fp), intent(in) :: Lc       ! Roughness length
      real(fp), intent(in) :: vegfrac  ! Vegetative fraction [0-1]
      real(fp), intent(in) :: thresh   ! Threshold for vegetation fraction
      real(fp), intent(out) :: dragpartition  ! Output drag partition

      !-------------------------------------------------------------------------

      ! !CONSTANTS:
      real(fp), parameter :: DRAG_MIN = 1.0e-3_fp  ! Minimum allowable drag partition
      real(fp), parameter :: sigb = 1.0_fp         ! Bare surface sigma
      real(fp), parameter :: mb = 0.5_fp           ! Bare surface m
      real(fp), parameter :: Betab = 90.0_fp       ! Bare surface Beta
      real(fp), parameter :: sigv = 1.45_fp        ! Vegetation sigma
      real(fp), parameter :: mv = 0.16_fp          ! Vegetation m
      real(fp), parameter :: Betav = 202.0_fp      ! Vegetation Beta

      ! !LOCAL VARIABLES:
      real(fp) :: Lc_veg        ! Vegetation roughness length
      real(fp) :: Lc_bare       ! Bare surface roughness length
      real(fp) :: feff_bare     ! Bare surface drag partition
      real(fp) :: feff_veg      ! Vegetation drag partition
      real(fp) :: feff          ! Total drag partition
      logical  :: skip          ! Flag to skip calculations
      real(fp) :: tmpVal        ! Temp value for numerical check

      ! Skip conditions logic
      skip = .false.
      if (vegfrac < 0.0_fp .or. vegfrac >= thresh) skip = .true.

      if (.not. skip) then
         ! Calculate vegetation effect
         Lc_veg = -0.35_fp * LOG(1.0_fp - vegfrac)
         feff_veg = calc_drag_partition(sigv, mv, Betav, Lc_veg)
      else
         feff_veg = DRAG_MIN
      endif

      ! Calculate bare surface effect
      Lc_bare = Lc / (1.0_fp - vegfrac)  ! Avoid numerical issues at high Lc
      tmpVal = 1 - sigb * mb * Lc_bare

      skip = .false.
      if (vegfrac < 0.0_fp .or. vegfrac >= thresh) skip = .true.
      if (.not. skip) skip = (Lc > 0.2_fp) .or. (tmpVal <= 0.0_fp)

      if (.not. skip) then
         feff_bare = calc_drag_partition(sigb, mb, Betab, Lc_bare)
      else
         feff_bare = DRAG_MIN
      endif

      ! Calculate total drag partition
      feff = feff_veg * feff_bare

      ! Apply bounds
      if (feff > 1.0_fp .or. feff < 1.0e-5_fp) then
         dragpartition = DRAG_MIN
      else
         dragpartition = feff
      endif

   end subroutine DarmenovaDragPartition

   !>
   !! Calculates drag partition for mixed surfaces
   !!
   !! DESCRIPTION: Computes the drag partition coefficient for surfaces with both
   !!              vegetative and bare components based on Leung's formulation
   !!
   !! REVISION HISTORY:
   !! 15Aug2024 B.Baker/NOAA    - Original implementation
   !!
   subroutine LeungDragPartition(Lc, lai, gvf, thresh, dragpartition)

      ! !USES:
      implicit NONE

      ! !INPUT PARAMETERS:
      real(fp), intent(in) :: Lc     ! Canopy length scale [m]
      real(fp), intent(in) :: lai    ! Leaf Area Index [m²/m²]
      real(fp), intent(in) :: gvf    ! Green Vegetation Fraction [0-1]
      real(fp), intent(in) :: thresh ! LAI threshold value [m²/m²]

      ! !OUTPUT PARAMETERS:
      real(fp), intent(out) :: dragpartition ! Drag partition coefficient [-]

      ! !LOCAL VARIABLES:
      real(fp) :: frac_bare ! Fraction of bare surface [0-1]
      real(fp) :: K         ! Normalized gap length [-]
      real(fp) :: feff_bare ! Effective drag partition for bare surfaces [-]
      real(fp) :: feff_veg  ! Effective drag partition for vegetated surfaces [-]
      real(fp) :: Rbare1    ! Intermediate bare surface calculation [-]
      real(fp) :: Rbare2    ! Intermediate bare surface calculation [-]
      real(fp) :: Lc_bare   ! Bare surface canopy length scale [m]
      real(fp) :: feff      ! Final effective drag partition [-]
      real(fp) :: tmpVal    ! Temporary calculation value [-]

      ! !CONSTANTS:
      !real(fp), parameter :: LAI_THR = 0.33_fp  ! LAI threshold [-]
      real(fp), parameter :: C = 4.8_fp         ! Empirical constant [-]
      real(fp), parameter :: F0 = 0.32_fp       ! Base efficiency factor [-]
      real(fp), parameter :: SIGB = 1.0_fp      ! Roughness density parameter [-]
      real(fp), parameter :: MB = 0.5_fp        ! Empirical constant [-]
      real(fp), parameter :: BETAB = 90.0_fp    ! Empirical constant [-]
      real(fp), parameter :: MIN_FEFF = 1.0E-5_fp ! Minimum allowable drag partition [-]
      real(fp), parameter :: MAX_FEFF = 1.0_fp    ! Maximum allowable drag partition [-]
      real(fp), parameter :: SMALL = 1.0E-10_fp   ! Small number to prevent division by zero

      ! Initialize variables
      feff_bare = 0.0_fp
      feff_veg = 0.0_fp

      ! Calculate bare surface fraction with bounds checking
      frac_bare = MAX(MIN(1.0_fp - lai / thresh, 1.0_fp), SMALL)

      ! Calculate vegetative component
      if ((lai <= 0.0_fp) .or. (lai >= thresh)) then
         feff_veg = 0.0_fp
      else
         K = 2.0_fp * (1.0_fp / MAX(1.0_fp - lai, SMALL) - 1.0_fp)
         feff_veg = (K + F0 * C) / (K + C)
      endif

      ! Calculate bare surface component
      if ((Lc <= 0.2_fp) .and. (Lc > 0.0_fp) .and. (lai < thresh)) then
         Lc_bare = Lc / MAX(frac_bare, SMALL)
         tmpVal = 1.0_fp - SIGB * MB * Lc_bare

         if (tmpVal > SMALL) then
            Rbare1 = 1.0_fp / SQRT(MAX(1.0_fp - SIGB * MB * Lc_bare, SMALL))
            Rbare2 = 1.0_fp / SQRT(1.0_fp + BETAB * MB * Lc_bare)
            feff_bare = Rbare1 * Rbare2
         else
            feff_bare = 0.0_fp
         endif
      else
         feff_bare = 0.0_fp
      endif

      ! Calculate final effective drag partition
      feff = (gvf * feff_veg**3 + frac_bare * feff_bare**3) ** (1.0_fp/3.0_fp)

      ! Apply bounds
      if (feff > MAX_FEFF .or. feff < MIN_FEFF) then
         dragpartition = MIN_FEFF
      else
         dragpartition = feff
      endif

   end subroutine LeungDragPartition

   !>
   !! \brief Computes the Drag Partition from MB95
   !!
   !! Marticorena, B. and Bergametti, G.: Modeling the atmospheric dust cycle:
   !! 1. Design of a soil-derived dust emission scheme,
   !! J. Geophys. Res.-Atmos., 100, 16415–16430, https://doi.org/10.1029/95JD00690, 1995
   !!
   !! \param z0 roughness length
   !! \param R Drag partition
   !!
   !! \ingroup catchem_dust_process
   !!!>
   subroutine MB95_DragPartition(z0, R)
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)  :: z0   !< roughness length [m]
      real(fp), intent(out) :: R    !< Drag partition (0-1)

      ! Local Variables
      real(fp), parameter :: z0s = 0.0008467_fp !< ideal roughness length of soil

      ! Initialize
      R = 0.0_fp

      !--------------------------------------------
      ! MB95 Drag Partition
      !--------------------------------------------
      R = 1.0_fp - (log(z0 / z0s ) / log(0.7_fp * (0.1_fp / z0s) ** 0.8_fp))
      return

   end subroutine MB95_DragPartition

   !>
   !! \brief Computes Draxler Hoirizontal Flux
   !!
   !! Draxler, R.R, D.A. Gillette, J.S. Kirkpatrick, and J. Heller (2001),
   !! Estimating PM10 air concentrations from dust storms in Iraq, Kuwait,
   !! and Saudi Arabia, Atm. Environ, 35: 4315-4330.
   !! https://doi.org/10.1016/S1352-2310(01)00159-5
   !!
   !! \param ustar friction velocity
   !! \param ustar_threshold dry threshold friction velocity
   !! \param R Drag partition
   !! \param H Soil Moisture Attenuation Factor
   !! \param HorizFlux Horizontal Mass Flux
   !!
   !!
   !! \ingroup catchem_dust_process
   !!!>
   subroutine Draxler_HorizFlux(ustar, ustar_threshold, R, H, HorizFlux)
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)  :: ustar            !< friction velocity [m/s]
      real(fp), intent(in)  :: ustar_threshold  !< dry threshold friction velocity [m/s]
      real(fp), intent(in)  :: R                !< Drag partition (0-1)
      real(fp), intent(in)  :: H                !< Soil Moisture Attenuation Factor
      real(fp), intent(inout) :: HorizFlux      !< Horizontal Mass Flux [kg/m2/s]

      ! Local Variables
      !----------------
      real(fp) :: u_ts    !< Modified threshold friction velocity

      ! Initialize
      !-----------
      HorizFlux = 0.0_fp

      !--------------------------------------------
      ! Compute Draxler Horizontal Flux
      !--------------------------------------------
      u_ts = ustar_threshold * H / R

      if (ustar >= ustar_threshold) then
         HorizFlux = max(0._fp ,(ustar * R) ** 3.0_fp * (1.0_fp - ( u_ts / ustar ) ** 2.0_fp))
      endif

   end subroutine Draxler_HorizFlux

   !>
   !! \brief Computes Kawamura Hoirizontal Flux
   !!
   !! Kawamura, R., 1951. Study on sand movement by wind. Report, 5(3), pp.95-112.
   !!
   !! Webb, N., Chappell, A., LeGrand, S., Ziegler, N., Edwards, B. 2020.
   !! A note on the use of drag partition in aeolian transport models.
   !! Aeolian Research. 42:100560. https://doi.org/10.1016/j.aeolia.2019.100560.
   !!
   !! \param ustar friction velocity
   !! \param ustar_threshold dry threshold friction velocity
   !! \param R Drag partition
   !! \param H Soil Moisture Attenuation Factor
   !! \param HorizFlux Horizontal Mass Flux
   !!
   !! \ingroup catchem_dust_process
   !!!>
   subroutine Kawamura_HorizFlux(ustar, ustar_threshold, R, H, HorizFlux)
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)  :: ustar           !< friction velocity [m/s]
      real(fp), intent(in)  :: ustar_threshold !< dry threshold friction velocity [m/s]
      real(fp), intent(in)  :: R               !< Drag partition (0-1)
      real(fp), intent(in)  :: H               !< Soil Moisture Attenuation Factor
      real(fp), intent(inout) :: HorizFlux     !<

      ! Local Variables
      real(fp) :: u_ts !< Modified threshold friction velocity

      ! Initialize
      HorizFlux = 0.0_fp

      !--------------------------------------------
      ! Compute Kawamura Horizontal Flux
      !--------------------------------------------
      u_ts = ustar_threshold * H / R

      HorizFlux = MAX(0._fp, (ustar ** 3.0_fp * (1.0_fp - (u_ts / ustar) ** 2.0_fp) * (1.0_fp + (u_ts / ustar) ** 2.0_fp ) ) )

   end subroutine Kawamura_HorizFlux

   !>
   !! \brief Computes White Horizontal Flux (used in GOCART2G)
   !!
   !! White, B. R. (1979). Soil transport by winds on Mars and Earth.
   !! JGR: Solid Earth, 84(B8), 4643–4651. https://doi.org/10.1029/JB084iB08p04643
   !!
   !! \param ustar friction velocity
   !! \param ustar_threshold dry threshold friction velocity
   !! \param R Drag partition
   !! \param H Soil Moisture Attenuation Factor
   !! \param HorizFlux Horizontal Mass Flux
   !!
   !! \ingroup catchem_dust_process
   !!!>
   subroutine White_HorizFlux(ustar, ustar_threshold, R, H, HorizFlux)
      IMPLICIT NONE
      ! Parameters
      real(fp), intent(in)  :: ustar           !< friction velocity [m/s]
      real(fp), intent(in)  :: ustar_threshold !< dry threshold friction velocity [m/s]
      real(fp), intent(in)  :: R               !< Drag partition (0-1)
      real(fp), intent(in)  :: H               !< Soil Moisture Attenuation Factor
      real(fp), intent(inout) :: HorizFlux     !<

      ! Local Variables
      real(fp) :: rustar !< Modified friction velocity
      real(fp) :: u_thresh !< Modified threshold friction velocity
      real(fp) :: u_sum !< Sum of modified friction and threshold velocities

      ! Initialize
      HorizFlux = 0.0_fp

      !--------------------------------------------
      ! Compute White Horizontal Flux
      !--------------------------------------------
      rustar = R * ustar
      ! Calculate threshold velocity
      u_thresh = ustar_threshold * H
      u_sum = rustar + u_thresh

      ! Calculate horizontal saltation flux
      HorizFlux = max(0.0_fp, rustar - u_thresh) * u_sum * u_sum

   end subroutine White_HorizFlux


end module DustScheme_FENGSHA_Mod
