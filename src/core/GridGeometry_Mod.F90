module GridGeometry_Mod
   implicit none
   private
   public :: GridGeometryType

   type :: GridGeometryType
      integer :: nx = 1
      integer :: ny = 1
      integer :: nz = 1
   contains
      procedure :: get_nx
      procedure :: get_ny
      procedure :: get_nz
      procedure :: get_dimensions
      procedure :: set
   end type GridGeometryType

contains

   function get_nx(this) result(nx)
      class(GridGeometryType), intent(in) :: this
      integer :: nx
      nx = this%nx
   end function get_nx

   function get_ny(this) result(ny)
      class(GridGeometryType), intent(in) :: this
      integer :: ny
      ny = this%ny
   end function get_ny

   function get_nz(this) result(nz)
      class(GridGeometryType), intent(in) :: this
      integer :: nz
      nz = this%nz
   end function get_nz

   subroutine get_dimensions(this, nx, ny, nz)
      class(GridGeometryType), intent(in) :: this
      integer, intent(out) :: nx, ny, nz
      nx = this%nx
      ny = this%ny
      nz = this%nz
   end subroutine get_dimensions

   !> \brief Set the grid dimensions
   !!
   !! \param[in] this GridGeometryType object
   !! \param[in] nx   Number of grid points in x-direction
   !! \param[in] ny   Number of grid points in y-direction
   !! \param[in] nz   Number of grid points in z-direction
   subroutine set(this, nx, ny, nz)
      class(GridGeometryType), intent(inout) :: this
      integer, intent(in) :: nx, ny, nz
      this%nx = nx
      this%ny = ny
      this%nz = nz
   end subroutine set

end module GridGeometry_Mod
