module CCPr_Chem_mod
   implicit none

   private
   public :: get_micm_version

contains

   function get_micm_version() result(res)
      use musica_util, only: string_t
      use musica_micm, only: get_micm_version_ => get_micm_version
      character(len=256) :: res

      type(string_t) :: micm_version_

      micm_version_ = get_micm_version_()
      res = micm_version_%get_char_array()

   end function get_micm_version

end module CCPr_Chem_mod
