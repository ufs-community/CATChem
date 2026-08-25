module DustCommon_Mod

   use precision_mod, only: fp
   implicit none
   public


   type :: DustConfig

      ! Process settings
      character(len=32) :: scheme = 'fengsha'
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
      integer, allocatable :: species_indices(:)  ! Indices of dust species in ChemState



      ! Species properties
      real(fp), allocatable :: species_density(:)      ! density for each species
      real(fp), allocatable :: species_lower_radius(:)      ! lower_radius for each species
      real(fp), allocatable :: species_radius(:)      ! radius for each species
      real(fp), allocatable :: species_upper_radius(:)      ! upper_radius for each species

      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   end type


   type :: DustSchemeFENGSHAConfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'fengsha'
      character(len=256) :: description = 'Fengsha Dust emission scheme developed at NOAA ARL for use at NOAA NWS'
      character(len=64) :: author = 'Barry Baker & Wei Li'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .false.  ! Surface-only processing

      ! Scheme parameters
      real(fp) :: alpha = 0.2_fp  ! linear scaling factor
      real(fp) :: gamma = 1.0_fp  ! Exponential scaling factor on source parameter
      real(fp) :: drylimit_factor = 1.0_fp  ! Dry Limit factor modifying the Fecan dry limit following Zender 2003
      real(fp) :: moist_correction_factor = 1.0_fp  ! Moisture correction factor
      real(fp) :: kvhmax = 0.0002_fp  ! Maximum vertical to horizontal flux ratio
      integer :: drag_option = 1  ! Drag Partition Option: 1 - use input drag, 2 - Darmenova, 3 - Leung 2022, 4 - MB95
      integer :: horizflux_option = 1  ! Horizontal flux option: 1 - White (1979), 2 - Draxler (2001), 3 - Kawamura (1964)
      integer :: moist_option = 1  ! Moisture parameterization: 1 - Fecan, 2 - Zhao
      integer :: distribution_option = 1  ! Dust Distribution option: 1 - Kok 2011, 2 - Meng 2022 (not implemented yet)

      ! Required meteorological fields
      integer :: n_required_met_fields = 15
      character(len=32) :: required_met_fields(15)

   end type


   type :: DustSchemeGINOUXConfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'ginoux'
      character(len=256) :: description = 'Ginoux dust emission scheme'
      character(len=64) :: author = 'Barry Baker & Wei Li'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .false.  ! Surface-only processing

      ! Scheme parameters
      real(fp) :: Ch_DU(5) = (/ 1.0_fp, 1.0_fp, 1.0_fp, 1.0_fp, 1.0_fp /)  ! Dust tuning coefficient per species bin

      ! Required meteorological fields
      integer :: n_required_met_fields = 9
      character(len=32) :: required_met_fields(9)

   end type


   type :: DustProcessConfig

      ! Process metadata
      character(len=64) :: process_name = 'dust'
      character(len=16) :: process_version = '1.0.0'
      logical :: is_active = .true.

      ! Process-specific configuration (delegate to DustConfig)
      type(DustConfig) :: dust_config

      ! Scheme configurations
      type(DustSchemeFENGSHAConfig) :: fengsha_config
      type(DustSchemeGINOUXConfig) :: ginoux_config


   end type

end module DustCommon_Mod
