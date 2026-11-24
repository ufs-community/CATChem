!> \brief Fixed character lengths shared by the NUOPC driver modules.
!!
!! Every name, path, keyword, and message buffer in drivers/nuopc is declared
!! with one of these constants, and the same constant sizes the buffers handed
!! to the C API, so a Fortran string and the C buffer it is copied into cannot
!! drift apart.  Values are the lengths the driver has always used.
module catchem_nuopc_lengths
   implicit none
   public

   integer, parameter :: NUOPC_SHORT_LEN = 16       !< tile suffixes, apply methods
   integer, parameter :: NUOPC_KEYWORD_LEN = 32     !< units, method and mode keywords
   integer, parameter :: NUOPC_NAME_LEN = 64        !< field, category, species, and Mie-table names (C API name buffers)
   integer, parameter :: NUOPC_LONG_NAME_LEN = 128  !< standard/long names, variable and tracer names
   integer, parameter :: NUOPC_PATH_LEN = 256       !< file paths, descriptions, configuration lines
   integer, parameter :: NUOPC_MSG_LEN = 512        !< log and error messages
   integer, parameter :: NUOPC_WARN_MSG_LEN = 4096  !< long diagnostic-selector warnings
end module catchem_nuopc_lengths
