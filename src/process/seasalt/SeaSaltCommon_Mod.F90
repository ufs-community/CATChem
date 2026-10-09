module SeaSaltCommon_Mod

   use catchem_bridge_precision, only: fp
   implicit none
   public


   type :: SeaSaltConfig

      ! Process settings
      character(len=32) :: scheme = 'gong97'
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
      integer, allocatable :: species_indices(:)  ! Indices of seasalt species in ChemState



      ! Species properties
      real(fp), allocatable :: species_density(:)      ! density for each species
      real(fp), allocatable :: species_lower_radius(:)      ! lower_radius for each species
      real(fp), allocatable :: species_radius(:)      ! radius for each species
      real(fp), allocatable :: species_upper_radius(:)      ! upper_radius for each species

      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   end type


   type :: SeaSaltSchemeGONG97Config

      ! Scheme metadata
      character(len=64) :: scheme_name = 'gong97'
      character(len=256) :: description = 'Gong 1997 sea salt emission scheme'
      character(len=64) :: author = 'Barry Baker'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .false.  ! Surface-only processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Emission scale factor
      logical :: weibull_flag = .false.  ! Apply Weibull distribution for particle size

      ! Required meteorological fields
      integer :: n_required_met_fields = 7
      character(len=32) :: required_met_fields(7)

   end type


   type :: SeaSaltSchemeGONG03Config

      ! Scheme metadata
      character(len=64) :: scheme_name = 'gong03'
      character(len=256) :: description = 'Gong 2003 sea salt emission scheme with improved sub- and super-micron treatment'
      character(len=64) :: author = 'Barry Baker'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .false.  ! Surface-only processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Emission scale factor
      logical :: weibull_flag = .false.  ! Apply Weibull distribution for particle size

      ! Required meteorological fields
      integer :: n_required_met_fields = 7
      character(len=32) :: required_met_fields(7)

   end type


   type :: SeaSaltSchemeGEOS12Config

      ! Scheme metadata
      character(len=64) :: scheme_name = 'geos12'
      character(len=256) :: description = 'GEOS-Chem 2012 sea salt emission scheme with observational constraints'
      character(len=64) :: author = 'Barry Baker'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .false.  ! Surface-only processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Emission scale factor
      logical :: weibull_flag = .false.  ! Apply Weibull distribution for particle size

      ! Required meteorological fields
      integer :: n_required_met_fields = 8
      character(len=32) :: required_met_fields(8)

   end type


   type :: SeaSaltProcessConfig

      ! Process metadata
      character(len=64) :: process_name = 'seasalt'
      character(len=16) :: process_version = '1.0.0'
      logical :: is_active = .true.

      ! Process-specific configuration (delegate to SeaSaltConfig)
      type(SeaSaltConfig) :: seasalt_config

      ! Scheme configurations
      type(SeaSaltSchemeGONG97Config) :: gong97_config
      type(SeaSaltSchemeGONG03Config) :: gong03_config
      type(SeaSaltSchemeGEOS12Config) :: geos12_config


   end type

end module SeaSaltCommon_Mod
