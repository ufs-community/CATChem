module SO4chemCommon_Mod

   use precision_mod, only: fp
   implicit none
   public


   type :: SO4chemConfig

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
      integer, allocatable :: species_indices(:)  ! Indices of so4chem species in ChemState



      ! Species properties
      real(fp), allocatable :: species_mw_g(:)      ! mw_g for each species

      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   end type


   type :: SO4chemSchemeGOCARTConfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'gocart'
      character(len=256) :: description = 'GOCART SO2 to SO4 production scheme'
      character(len=64) :: author = 'Wei Li'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .true.  ! Full column processing

      ! Scheme parameters
      logical :: update_so2 = .true.  ! whether to update SO2 concentration based on chemical production/loss

      ! Required meteorological fields
      integer :: n_required_met_fields = 16
      character(len=32) :: required_met_fields(16)

   end type


   type :: SO4chemGOCARTPersistentState
      logical :: firsttime  ! flag for first time step
      integer :: nymd_last  ! last day of H2O2 update
      integer :: nhms_last_recycle  ! last time step of H2O2 recycle
      real(fp), allocatable :: xh2o2_init(:)  ! H2O2 column initialization

   end type


   type :: SO4chemProcessConfig

      ! Process metadata
      character(len=64) :: process_name = 'so4chem'
      character(len=16) :: process_version = '1.0.0'
      logical :: is_active = .true.

      ! Process-specific configuration (delegate to SO4chemConfig)
      type(SO4chemConfig) :: so4chem_config

      ! Scheme configurations
      type(SO4chemSchemeGOCARTConfig) :: gocart_config

      ! Persistent state arrays for column processing
      type(SO4chemGocartPersistentState), allocatable :: gocart_persistent_state(:)  ! Per-column state for gocart scheme
      integer :: total_columns = 0  ! Total number of columns for state allocation

   end type

end module SO4chemCommon_Mod
