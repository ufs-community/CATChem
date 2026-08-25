!> \file DiagnosticManager_Mod.F90
!! \brief Lightweight backward-compatible Fortran wrapper for DiagnosticManager procedures
!!
module DiagnosticManager_Mod
   use Precision_Mod, only: fp
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE
   use DiagnosticInterface_Mod, only: DiagnosticRegistryType

   implicit none
   private

   public :: DiagnosticManagerType

   type :: DiagnosticManagerType
      character(len=64) :: dummy = ""
   contains
      procedure :: list_processes => diag_mgr_list_processes
      procedure :: get_process_registry => diag_mgr_get_process_registry
      procedure :: get_field_value => diag_mgr_get_field_value
      procedure :: register_process => diag_mgr_register_process
   end type DiagnosticManagerType

contains

   subroutine diag_mgr_list_processes(this, process_list, num_processes, rc)
      class(DiagnosticManagerType), intent(in) :: this
      character(len=*), allocatable, intent(out) :: process_list(:)
      integer, intent(out) :: num_processes
      integer, intent(out) :: rc

      associate(unused => this)
      end associate

      allocate(process_list(0))
      num_processes = 0
      rc = CC_SUCCESS
   end subroutine diag_mgr_list_processes

   subroutine diag_mgr_get_process_registry(this, process_name, registry, rc)
      class(DiagnosticManagerType), intent(in) :: this
      character(len=*), intent(in) :: process_name
      type(DiagnosticRegistryType), pointer, intent(out) :: registry
      integer, intent(out) :: rc

      associate(unused1 => this, unused2 => process_name)
      end associate

      nullify(registry)
      rc = CC_SUCCESS
   end subroutine diag_mgr_get_process_registry

   subroutine diag_mgr_get_field_value(this, process_name, field_name, &
      scalar_value, array_1d_ptr, array_2d_ptr, array_3d_ptr, &
      data_type, description, units, rc)
      class(DiagnosticManagerType), intent(in) :: this
      character(len=*), intent(in) :: process_name
      character(len=*), intent(in) :: field_name
      real(fp), intent(out) :: scalar_value
      real(fp), pointer, intent(out) :: array_1d_ptr(:)
      real(fp), pointer, intent(out) :: array_2d_ptr(:,:)
      real(fp), pointer, intent(out) :: array_3d_ptr(:,:,:)
      integer, intent(out) :: data_type
      character(len=*), intent(out) :: description
      character(len=*), intent(out) :: units
      integer, intent(out) :: rc

      associate(unused1 => this, unused2 => process_name, unused3 => field_name)
      end associate

      scalar_value = 0.0_fp
      nullify(array_1d_ptr, array_2d_ptr, array_3d_ptr)
      data_type = 0
      description = ""
      units = ""
      rc = CC_SUCCESS
   end subroutine diag_mgr_get_field_value

   subroutine diag_mgr_register_process(this, process_name, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      integer, intent(out) :: rc

      associate(unused1 => this, unused2 => process_name)
      end associate

      rc = CC_SUCCESS
   end subroutine diag_mgr_register_process

end module DiagnosticManager_Mod
