module DryDepCommon_Mod

   use catchem_bridge_precision, only: fp
   implicit none
   public


   type :: DryDepConfig

      ! Process settings
      character(len=32) :: gas_scheme = 'wesely'   ! Scheme for gas species
      character(len=32) :: aero_scheme = 'gocart'  ! Scheme for aerosol species
      logical :: is_active = .true.
      logical :: diagnostics = .false.  ! Diagnostic switch

      ! Diagnostic species configuration
      integer :: n_diagnostic_species = 0
      character(len=32), allocatable :: diagnostic_species(:)  ! User-defined species for diagnostics
      integer, allocatable :: diagnostic_species_id(:)  ! Indices mapping diagnostic_species to species_names
      real(fp) :: dt_min = 1.0_fp     ! Minimum time step (seconds)
      real(fp) :: dt_max = 3600.0_fp  ! Maximum time step (seconds)

      ! Species configuration
      integer :: n_species = 0
      character(len=32), allocatable :: species_names(:)
      integer, allocatable :: species_indices(:)  ! Indices of drydep species in ChemState
      logical, allocatable :: is_gas(:)           ! Is gas species (true for gas, false for aerosol)



      ! Species properties
      real(fp), allocatable :: species_dd_DvzAerSnow(:)      ! dd_DvzAerSnow for each species
      real(fp), allocatable :: species_dd_DvzMinVal_land(:)      ! dd_DvzMinVal_land for each species
      real(fp), allocatable :: species_dd_DvzMinVal_snow(:)      ! dd_DvzMinVal_snow for each species
      real(fp), allocatable :: species_dd_f0(:)      ! dd_f0 for each species
      real(fp), allocatable :: species_dd_hstar(:)      ! dd_hstar for each species
      real(fp), allocatable :: species_density(:)      ! density for each species
      logical, allocatable :: species_is_dust(:)      ! is_dust for each species
      logical, allocatable :: species_is_seasalt(:)      ! is_seasalt for each species
      real(fp), allocatable :: species_lower_radius(:)      ! lower_radius for each species
      real(fp), allocatable :: species_mw_g(:)      ! mw_g for each species
      real(fp), allocatable :: species_radius(:)      ! radius for each species
      real(fp), allocatable :: species_upper_radius(:)      ! upper_radius for each species

      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   end type


   type :: DryDepSchemeWESELYConfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'wesely'
      character(len=256) :: description = 'Wesely 1989 gas dry deposition scheme'
      character(len=64) :: author = 'Wei Li'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .false.  ! Surface-only processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! DryDep velocity scale factor
      logical :: co2_effect = .true.  ! Apply CO2 effect on stomatal conductance
      real(fp) :: co2_level = 600.0  ! Ambient CO2 level for stomatal conductance adjustment
      real(fp) :: co2_reference = 380.0  ! Reference CO2 level for stomatal conductance adjustment

      ! Required meteorological fields
      integer :: n_required_met_fields = 21
      character(len=32) :: required_met_fields(21)

   end type


   type :: DryDepSchemeGOCARTConfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'gocart'
      character(len=256) :: description = 'GOCART-2G aerosol dry deposition scheme'
      character(len=64) :: author = 'Wei Li & Lacey Holland'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .false.  ! Surface-only processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Dry deposition velocity scale factor
      logical :: resuspension = .false.  ! Apply resuspension for dry deposition
      logical :: dust_resuspension_only = .true.  ! If true, resuspension only applies to dust species

      ! Required meteorological fields
      integer :: n_required_met_fields = 13
      character(len=32) :: required_met_fields(13)

   end type


   type :: DryDepSchemeZHANGConfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'zhang'
      character(len=256) :: description = 'Zhang et al. [2001] scheme with Emerson et al. [2020] updates'
      character(len=64) :: author = 'Wei Li'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .false.  ! Surface-only processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Dry deposition velocity scale factor

      ! Required meteorological fields
      integer :: n_required_met_fields = 15
      character(len=32) :: required_met_fields(15)

   end type


   type :: DryDepProcessConfig

      ! Process metadata
      character(len=64) :: process_name = 'drydep'
      character(len=16) :: process_version = '1.0.0'
      logical :: is_active = .true.

      ! Process-specific configuration (delegate to DryDepConfig)
      type(DryDepConfig) :: drydep_config

      ! Scheme configurations
      ! Separate gas and aerosol scheme configurations
      type(DryDepSchemeWESELYConfig) :: wesely_config
      type(DryDepSchemeGOCARTConfig) :: gocart_config
      type(DryDepSchemeZHANGConfig) :: zhang_config

   end type

end module DryDepCommon_Mod
