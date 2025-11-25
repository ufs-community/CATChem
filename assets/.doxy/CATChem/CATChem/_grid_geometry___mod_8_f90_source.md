

# File GridGeometry\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**GridGeometry\_Mod.F90**](_grid_geometry___mod_8_f90.md)

[Go to the documentation of this file](_grid_geometry___mod_8_f90.md)


```Fortran
module gridgeometry_mod
   implicit none
   private
   public :: gridgeometrytype

   type :: gridgeometrytype
      integer :: nx = 1
      integer :: ny = 1
      integer :: nz = 1
   contains
      procedure :: get_nx
      procedure :: get_ny
      procedure :: get_nz
      procedure :: get_dimensions
      procedure :: set
   end type gridgeometrytype

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

   subroutine set(this, nx, ny, nz)
      class(GridGeometryType), intent(inout) :: this
      integer, intent(in) :: nx, ny, nz
      this%nx = nx
      this%ny = ny
      this%nz = nz
   end subroutine set

end module gridgeometry_mod
```


