module WetDepCommon_Mod

   use precision_mod, only: fp
   implicit none
   public


   type :: WetDepConfig

      ! Process settings
      character(len=32) :: scheme = 'jacob'
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
      integer, allocatable :: species_indices(:)  ! Indices of wetdep species in ChemState



      ! Species properties
      real(fp), allocatable :: species_henry_cr(:)      ! henry_cr for each species
      real(fp), allocatable :: species_henry_k0(:)      ! henry_k0 for each species
      real(fp), allocatable :: species_henry_pKa(:)      ! henry_pKa for each species
      logical, allocatable :: species_is_aerosol(:)      ! is_aerosol for each species
      real(fp), allocatable :: species_mw_g(:)      ! mw_g for each species
      real(fp), allocatable :: species_radius(:)      ! radius for each species
      logical, allocatable :: species_wd_LiqAndGas(:)      ! wd_LiqAndGas for each species
      real(fp), allocatable :: species_wd_convfacI2G(:)      ! wd_convfacI2G for each species
      real(fp), allocatable :: species_wd_rainouteff(:,:)      ! wd_rainouteff for each species
      real(fp), allocatable :: species_wd_retfactor(:)      ! wd_retfactor for each species
      real(fp), allocatable :: species_wd_reevap_frac(:)      ! wd_reevap_frac for each species

      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   end type


   type :: WetDepSchemeJACOBConfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'jacob'
      character(len=256) :: description = 'Jacob et al. [2000] wet deposition scheme'
      character(len=64) :: author = 'Wei Li'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .true.  ! Full column processing

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Washout tuning factor
      real(fp) :: radius_threshold = 1.0  ! Radius threshold for aerosol wet deposition (um)
      logical :: so4_gocart_resusp = .true.  ! Use GOCART SU_Wet_Removal resuspension (alpha) for sulfate (SO4/SO2) only
      real(fp) :: so4_washout_eff = 1.0  ! Sulfate-only below-cloud washout efficiency multiplier (SO4 column tuning; 1.0 = unchanged)

      ! Required meteorological fields
      integer :: n_required_met_fields = 8
      character(len=32) :: required_met_fields(8)

   end type


   type :: WetDepProcessConfig

      ! Process metadata
      character(len=64) :: process_name = 'wetdep'
      character(len=16) :: process_version = '1.0.0'
      logical :: is_active = .true.

      ! Process-specific configuration (delegate to WetDepConfig)
      type(WetDepConfig) :: wetdep_config

      ! Scheme configurations
      type(WetDepSchemeJACOBConfig) :: jacob_config

   end type

end module WetDepCommon_Mod
