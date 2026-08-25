!> \file ConfigManager_Mod.F90
!! \brief Lightweight backward-compatible Fortran wrapper for configuration values
!!
module ConfigManager_Mod
   use Precision_Mod, only: fp
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE

   implicit none
   private

   public :: ConfigManagerType, ConfigDataType, EmissionCategoryMapping, EmisSpeciesMappingEntry, EmissionMappingConfig, &
      RuntimeConfig, FilePathConfig

   type :: EmisSpeciesMappingEntry
      character(len=128) :: emission_field = ""
      character(len=64) :: units = ""
      character(len=128) :: long_name = ""
      integer :: n_mappings = 0
      character(len=64), allocatable :: map(:)
      real(fp), allocatable :: scale(:)
      integer, allocatable :: index(:)
   end type EmisSpeciesMappingEntry

   type :: EmissionCategoryMapping
      character(len=128) :: category_name = ""
      logical :: is_active = .false.
      integer :: n_emission_species = 0
      type(EmisSpeciesMappingEntry), allocatable :: species_mappings(:)
   end type EmissionCategoryMapping

   type :: EmissionMappingConfig
      logical :: is_loaded = .false.
      integer :: n_categories = 0
      type(EmissionCategoryMapping), allocatable :: categories(:)
   end type EmissionMappingConfig

   type :: RuntimeConfig
      integer :: Output_Frequency = 1
      integer :: CompressLev = 0
      logical :: latlon_output = .false.
      logical :: DiagEnabled = .true.
      character(len=64), allocatable :: diag_species(:)
   end type RuntimeConfig

   type :: FilePathConfig
      character(len=256) :: Output_Directory = "./"
      character(len=256) :: Output_Prefix = "catchem_diag"
   end type FilePathConfig

   type :: ConfigDataType
      type(EmissionMappingConfig) :: emission_mapping
      type(RuntimeConfig) :: runtime
      type(FilePathConfig) :: file_paths
   end type ConfigDataType

   type :: ConfigManagerType
      type(ConfigDataType) :: config_data
   contains
      procedure :: get_logical => config_get_logical
      procedure :: get_string => config_get_string
      procedure :: get_real => config_get_real
      procedure :: get_array => config_get_array
   end type ConfigManagerType

contains

   subroutine config_get_logical(this, path, val, rc, default)
      class(ConfigManagerType), intent(in) :: this
      character(len=*), intent(in) :: path
      logical, intent(out) :: val
      integer, intent(out) :: rc
      logical, intent(in), optional :: default

      associate(unused1 => this, unused2 => path)
      end associate

      val = .false.
      if (present(default)) val = default
      rc = CC_SUCCESS
   end subroutine config_get_logical

   subroutine config_get_string(this, path, val, rc, default)
      class(ConfigManagerType), intent(in) :: this
      character(len=*), intent(in) :: path
      character(len=*), intent(out) :: val
      integer, intent(out) :: rc
      character(len=*), intent(in), optional :: default

      associate(unused1 => this, unused2 => path)
      end associate

      val = ""
      if (present(default)) val = default
      rc = CC_SUCCESS
   end subroutine config_get_string

   subroutine config_get_real(this, path, val, rc, default)
      class(ConfigManagerType), intent(in) :: this
      character(len=*), intent(in) :: path
      real(fp), intent(out) :: val
      integer, intent(out) :: rc
      real(fp), intent(in), optional :: default

      associate(unused1 => this, unused2 => path)
      end associate

      val = 0.0_fp
      if (present(default)) val = default
      rc = CC_SUCCESS
   end subroutine config_get_real

   subroutine config_get_array(this, path, val, rc, default_values)
      class(ConfigManagerType), intent(in) :: this
      character(len=*), intent(in) :: path
      character(len=*), allocatable, intent(out) :: val(:)
      integer, intent(out) :: rc
      character(len=*), intent(in), optional :: default_values(:)

      associate(unused1 => this, unused2 => path)
      end associate

      if (present(default_values)) then
         allocate(val(size(default_values)))
         val = default_values
      else
         allocate(val(0))
      end if
      rc = CC_SUCCESS
   end subroutine config_get_array

end module ConfigManager_Mod
