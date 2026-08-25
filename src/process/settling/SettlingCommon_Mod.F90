module SettlingCommon_Mod

   use precision_mod, only: fp
   implicit none
   public


   type :: SettlingConfig

      ! Process settings
      character(len=32) :: scheme = 'gocart'
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
      integer, allocatable :: species_indices(:)  ! Indices of settling species in ChemState



      ! Species properties
      real(fp), allocatable :: species_density(:)      ! density for each species
      real(fp), allocatable :: species_mie_map(:)      ! mie_map for each species
      real(fp), allocatable :: species_radius(:)      ! radius for each species

      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   end type


   type :: SettlingSchemeGOCARTConfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'gocart'
      character(len=256) :: description = 'GOCART gravitational settling scheme'
      character(len=64) :: author = 'Wei Li'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .true.  ! Full column processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! settling velocity factor
      logical :: simple_scheme = .false.  ! read in mie data for wet particles if true; otherwise calculate particles wet swelling internally
      integer :: swelling_method = 1  ! method for calculating particle swelling: 1 Fitzgerald 1975; 2 for Gerber 1985
      logical :: correction_maring = .false.  ! correct the settling velocity following Maring et al, 2003

      ! Required meteorological fields
      integer :: n_required_met_fields = 7
      character(len=32) :: required_met_fields(7)

   end type


   type :: SettlingProcessConfig

      ! Process metadata
      character(len=64) :: process_name = 'settling'
      character(len=16) :: process_version = '1.0.0'
      logical :: is_active = .true.

      ! Process-specific configuration (delegate to SettlingConfig)
      type(SettlingConfig) :: settling_config

      ! Scheme configurations
      type(SettlingSchemeGOCARTConfig) :: gocart_config

   end type

end module SettlingCommon_Mod
