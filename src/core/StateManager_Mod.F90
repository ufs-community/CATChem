!> \file StateManager_Mod.F90
!! \brief Lightweight backward-compatible Fortran wrapper for StateManager pointer associations
!!
module StateManager_Mod
   use Precision_Mod, only: fp
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE, ErrorManagerType
   use MetState_Mod, only: MetStateType
   use ChemState_Mod, only: ChemStateType
   use ConfigManager_Mod, only: ConfigManagerType
   use TimeState_Mod, only: TimeStateType
   use iso_c_binding, only: c_ptr, c_null_ptr, c_associated

   implicit none
   private

   public :: StateManagerType

   type :: StateManagerType
      type(c_ptr) :: cpp_ptr = c_null_ptr
      real(fp) :: tstep = 0.0_fp
      type(MetStateType), pointer :: met_state => null()
      type(ChemStateType), pointer :: chem_state => null()
      type(ConfigManagerType), pointer :: config_mgr => null()
      type(ErrorManagerType), pointer :: error_mgr => null()
      type(TimeStateType), pointer :: time_state => null()
   contains
      procedure :: get_met_state_ptr => state_mgr_get_met_state_ptr
      procedure :: get_chem_state_ptr => state_mgr_get_chem_state_ptr
      procedure :: get_config_ptr => state_mgr_get_config_ptr
      procedure :: get_error_manager => state_mgr_get_error_mgr
      procedure :: get_time_state_ptr => state_mgr_get_time_state_ptr
   end type StateManagerType

contains

   function state_mgr_get_met_state_ptr(this) result(ptr)
      use Interop_Mod, only: get_cpp_field
      class(StateManagerType), intent(in) :: this
      type(MetStateType), pointer :: ptr
      integer :: rc, nx, ny, nz

      ptr => this%met_state
      if (associated(ptr) .and. c_associated(this%cpp_ptr)) then
         ptr%cpp_ptr = this%cpp_ptr
         call ptr%geometry%get_dimensions(nx, ny, nz)

         ! Bind volumetric 3D arrays
         call get_cpp_field(this%cpp_ptr, "T", ptr%T, [nx, ny, nz], rc)
         call get_cpp_field(this%cpp_ptr, "QV", ptr%QV, [nx, ny, nz], rc)
         call get_cpp_field(this%cpp_ptr, "RH", ptr%RH, [nx, ny, nz], rc)
         call get_cpp_field(this%cpp_ptr, "PMID", ptr%PMID, [nx, ny, nz], rc)
         call get_cpp_field(this%cpp_ptr, "PEDGE", ptr%PEDGE, [nx, ny, nz+1], rc)
         call get_cpp_field(this%cpp_ptr, "AIRDEN", ptr%AIRDEN, [nx, ny, nz], rc)
         call get_cpp_field(this%cpp_ptr, "AIRDEN_DRY", ptr%AIRDEN_DRY, [nx, ny, nz], rc)
         call get_cpp_field(this%cpp_ptr, "BXHEIGHT", ptr%BXHEIGHT, [nx, ny, nz], rc)
         call get_cpp_field(this%cpp_ptr, "DELP", ptr%DELP, [nx, ny, nz], rc)
         call get_cpp_field(this%cpp_ptr, "DELP_DRY", ptr%DELP_DRY, [nx, ny, nz], rc)

         ! Bind surface 2D arrays
         call get_cpp_field(this%cpp_ptr, "PS", ptr%PS, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "TS", ptr%TS, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "PBLH", ptr%PBLH, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "USTAR", ptr%USTAR, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "HFLUX", ptr%HFLUX, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "OBK", ptr%OBK, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "LAT", ptr%LAT, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "LON", ptr%LON, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "FROCEAN", ptr%FROCEAN, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "FRSEAICE", ptr%FRSEAICE, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "SST", ptr%SST, [nx, ny], rc)
      end if
   end function state_mgr_get_met_state_ptr

   function state_mgr_get_chem_state_ptr(this) result(ptr)
      class(StateManagerType), intent(in) :: this
      type(ChemStateType), pointer :: ptr
      ptr => this%chem_state
   end function state_mgr_get_chem_state_ptr

   function state_mgr_get_config_ptr(this) result(ptr)
      class(StateManagerType), intent(in) :: this
      type(ConfigManagerType), pointer :: ptr
      ptr => this%config_mgr
   end function state_mgr_get_config_ptr

   function state_mgr_get_error_mgr(this) result(ptr)
      class(StateManagerType), intent(in) :: this
      type(ErrorManagerType), pointer :: ptr
      ptr => this%error_mgr
   end function state_mgr_get_error_mgr

   function state_mgr_get_time_state_ptr(this) result(ptr)
      class(StateManagerType), intent(in) :: this
      type(TimeStateType), pointer :: ptr
      ptr => this%time_state
   end function state_mgr_get_time_state_ptr

end module StateManager_Mod
