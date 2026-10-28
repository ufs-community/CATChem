

# File CarbChemCommon\_Mod.F90

[**File List**](files.md) **>** [**carbchem**](dir_5dbdd03f815becc4c35f94e0692b4e09.md) **>** [**CarbChemCommon\_Mod.F90**](_carb_chem_common___mod_8_f90.md)

[Go to the documentation of this file](_carb_chem_common___mod_8_f90.md)


```Fortran
module carbchemcommon_mod

   use catchem_bridge_precision, only: fp
   implicit none
   public


   type :: carbchemconfig

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
      integer, allocatable :: species_indices(:)  ! Indices of carbchem species in ChemState



      ! Species properties
      real(fp), allocatable :: species_t_chem_loss(:)      ! t_chem_loss for each species

      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   end type


   type :: carbchemschemegocartconfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'gocart'
      character(len=256) :: description = 'GOCART carbon species chemical production and loss scheme'
      character(len=64) :: author = 'Wei Li'
      character(len=16) :: algorithm_type = 'explicit'

      ! Process configuration
      logical :: affects_full_column = .true.  ! Full column processing

      ! Scheme parameters
      real(fp) :: time_days_hydrophobic_to_hydrophilic = 2.5  ! Rate of conversion of hydrophobic to hydrophilic [days]

      ! Required meteorological fields
      integer :: n_required_met_fields = 4
      character(len=32) :: required_met_fields(4)

   end type


   type :: carbchemprocessconfig

      ! Process metadata
      character(len=64) :: process_name = 'carbchem'
      character(len=16) :: process_version = '1.0.0'
      logical :: is_active = .true.

      ! Process-specific configuration (delegate to CarbChemConfig)
      type(CarbChemConfig) :: carbchem_config

      ! Scheme configurations
      type(CarbChemSchemeGOCARTConfig) :: gocart_config


   end type

end module carbchemcommon_mod
```


