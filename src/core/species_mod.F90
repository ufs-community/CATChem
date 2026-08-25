!> \file species_mod.F90
!! \brief Compatibility-preserving Fortran proxy mirroring the unified C++ species metadata database
!!
module species_mod
   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE
   use iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_int, c_double, c_char, c_null_char

   implicit none
   private

   public :: SpeciesType, SpeciesManagerType, validate_species, find_species_by_name, populate_species_from_cpp

   ! Missing value sentinel
   real(fp), parameter, public :: MISSING_VV = 1.0e-20_fp

   !> Derived type representing a chemical species, populated dynamically from C++ records
   type :: SpeciesType
      character(len=30) :: long_name = ""
      character(len=30) :: short_name = ""
      character(len=50) :: description = ""

      logical :: is_gas = .true.
      logical :: is_aerosol = .false.
      logical :: is_tracer = .false.
      logical :: is_advected = .true.
      logical :: is_drydep = .false.
      logical :: is_wetdep = .false.
      logical :: is_photolysis = .false.
      logical :: is_gocart_aero = .false.
      logical :: is_dust = .false.
      logical :: is_seasalt = .false.

      real(fp) :: mw_g = 0.0_fp
      real(fp) :: density = 1000.0_fp
      real(fp) :: radius = 1.0e-9_fp
      real(fp) :: lower_radius = 0.0_fp
      real(fp) :: upper_radius = 0.0_fp
      real(fp) :: viscosity = 0.0_fp

      ! Dry deposition
      real(fp) :: dd_f0 = 0.0_fp
      real(fp) :: dd_hstar = 0.0_fp
      real(fp) :: dd_DvzAerSnow = 0.0_fp
      real(fp) :: dd_DvzMinVal_snow = 0.0_fp
      real(fp) :: dd_DvzMinVal_land = 0.0_fp

      ! Wet deposition
      real(fp) :: henry_k0 = 0.0_fp
      real(fp) :: henry_cr = 0.0_fp
      real(fp) :: henry_pKa = 0.0_fp
      real(fp) :: wd_retfactor = 0.0_fp
      logical :: wd_LiqAndGas = .false.
      real(fp) :: wd_convfacI2G = 0.0_fp
      real(fp) :: wd_rainouteff(3) = 0.0_fp
      real(fp) :: wd_reevap_frac = 0.5_fp

      character(len=30) :: mie_name = ""
      real(fp) :: t_chem_loss = -1.0_fp
      real(fp) :: BackgroundVV = MISSING_VV

      integer :: species_index = -1
      integer :: drydep_index = -1
      integer :: photolysis_index = -1
      integer :: gocart_aero_index = -1

      real(fp), pointer :: conc(:,:,:) => null()
      logical :: is_valid = .false.
   contains
      procedure :: init => species_init
      procedure :: validate => species_validate
      procedure :: cleanup => species_cleanup
      procedure :: set_concentration => species_set_concentration
      procedure :: get_concentration => species_get_concentration
      procedure :: copy => species_copy
      procedure :: print_info => species_print_info
   end type SpeciesType

   !> Compatibility manager database
   type :: SpeciesManagerType
      type(SpeciesType), allocatable :: species_db(:)
      integer :: num_species = 0
      logical :: is_initialized = .false.
   contains
      procedure :: init => species_manager_init
      procedure :: add_species => species_manager_add_species
      procedure :: find_species => species_manager_find_species
      procedure :: validate_database => species_manager_validate_database
      procedure :: load_from_file => species_manager_load_from_cpp
      procedure :: cleanup => species_manager_cleanup
      procedure :: print_database => species_manager_print_database
   end type SpeciesManagerType

   ! C Interoperable Interface Definitions
   interface
      integer(c_int) function catchem_state_get_species_count(state_ptr) bind(C, name="catchem_state_get_species_count")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
      end function

      subroutine catchem_state_get_species_name_at(state_ptr, index, name_out) bind(C, name="catchem_state_get_species_name_at")
         import :: c_ptr, c_int, c_char
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
         character(kind=c_char), intent(out) :: name_out(*)
      end subroutine

      subroutine catchem_state_get_species_long_name_at(state_ptr, index, name_out) bind(C, name="catchem_state_get_species_long_name_at")
         import :: c_ptr, c_int, c_char
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
         character(kind=c_char), intent(out) :: name_out(*)
      end subroutine

      subroutine catchem_state_get_species_desc_at(state_ptr, index, name_out) bind(C, name="catchem_state_get_species_desc_at")
         import :: c_ptr, c_int, c_char
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
         character(kind=c_char), intent(out) :: name_out(*)
      end subroutine

      real(c_double) function catchem_state_get_species_density(state_ptr, index) bind(C, name="catchem_state_get_species_density")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_radius(state_ptr, index) bind(C, name="catchem_state_get_species_radius")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_lower_radius(state_ptr, index) bind(C, name="catchem_state_get_species_lower_radius")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_upper_radius(state_ptr, index) bind(C, name="catchem_state_get_species_upper_radius")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_viscosity(state_ptr, index) bind(C, name="catchem_state_get_species_viscosity")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      integer(c_int) function catchem_state_get_species_is_tracer(state_ptr, index) bind(C, name="catchem_state_get_species_is_tracer")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      integer(c_int) function catchem_state_get_species_is_advected(state_ptr, index) bind(C, name="catchem_state_get_species_is_advected")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      integer(c_int) function catchem_state_get_species_is_drydep(state_ptr, index) bind(C, name="catchem_state_get_species_is_drydep")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      integer(c_int) function catchem_state_get_species_is_wetdep(state_ptr, index) bind(C, name="catchem_state_get_species_is_wetdep")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      integer(c_int) function catchem_state_get_species_is_photolysis(state_ptr, index) bind(C, name="catchem_state_get_species_is_photolysis")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      integer(c_int) function catchem_state_get_species_is_dust(state_ptr, index) bind(C, name="catchem_state_get_species_is_dust")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      integer(c_int) function catchem_state_get_species_is_seasalt(state_ptr, index) bind(C, name="catchem_state_get_species_is_seasalt")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_dd_f0(state_ptr, index) bind(C, name="catchem_state_get_species_dd_f0")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_dd_hstar(state_ptr, index) bind(C, name="catchem_state_get_species_dd_hstar")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_dd_DvzAerSnow(state_ptr, index) bind(C, name="catchem_state_get_species_dd_DvzAerSnow")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_dd_DvzMinVal_snow(state_ptr, index) bind(C, name="catchem_state_get_species_dd_DvzMinVal_snow")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_dd_DvzMinVal_land(state_ptr, index) bind(C, name="catchem_state_get_species_dd_DvzMinVal_land")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_henry_k0(state_ptr, index) bind(C, name="catchem_state_get_species_henry_k0")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_henry_cr(state_ptr, index) bind(C, name="catchem_state_get_species_henry_cr")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_henry_pKa(state_ptr, index) bind(C, name="catchem_state_get_species_henry_pKa")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_wd_retfactor(state_ptr, index) bind(C, name="catchem_state_get_species_wd_retfactor")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      integer(c_int) function catchem_state_get_species_wd_LiqAndGas(state_ptr, index) bind(C, name="catchem_state_get_species_wd_LiqAndGas")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_wd_convfacI2G(state_ptr, index) bind(C, name="catchem_state_get_species_wd_convfacI2G")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      subroutine catchem_state_get_species_wd_rainouteff(state_ptr, index, eff_out) bind(C, name="catchem_state_get_species_wd_rainouteff")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
         real(c_double), intent(out) :: eff_out(3)
      end subroutine

      real(c_double) function catchem_state_get_species_wd_reevap_frac(state_ptr, index) bind(C, name="catchem_state_get_species_wd_reevap_frac")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_t_chem_loss(state_ptr, index) bind(C, name="catchem_state_get_species_t_chem_loss")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_BackgroundVV(state_ptr, index) bind(C, name="catchem_state_get_species_BackgroundVV")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      subroutine catchem_state_get_species_mie_name(state_ptr, index, name_out) bind(C, name="catchem_state_get_species_mie_name")
         import :: c_ptr, c_int, c_char
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
         character(kind=c_char), intent(out) :: name_out(*)
      end subroutine

      real(c_double) function catchem_state_get_species_mw(state_ptr, index) bind(C, name="catchem_state_get_species_mw")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      integer(c_int) function catchem_state_is_species_gas(state_ptr, index) bind(C, name="catchem_state_is_species_gas")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      integer(c_int) function catchem_state_is_species_aerosol(state_ptr, index) bind(C, name="catchem_state_is_species_aerosol")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function
   end interface

contains

   ! Helper to convert null-terminated C buffers back to Fortran fixed strings
   subroutine c_to_f_string(c_str, f_str)
      character(kind=c_char), intent(in) :: c_str(*)
      character(len=*), intent(out) :: f_str
      integer :: i
      f_str = ""
      do i = 1, len(f_str)
         if (c_str(i) == c_null_char) exit
         f_str(i:i) = c_str(i)
      end do
   end subroutine c_to_f_string

   subroutine populate_species_from_cpp(this, state_ptr, index, rc)
      type(SpeciesType), intent(inout) :: this
      type(c_ptr), intent(in) :: state_ptr
      integer, intent(in) :: index
      integer, intent(out) :: rc

      character(kind=c_char) :: c_buf(128)
      real(c_double) :: r_eff(3)

      rc = CC_SUCCESS

      call catchem_state_get_species_name_at(state_ptr, int(index, c_int), c_buf)
      call c_to_f_string(c_buf, this%short_name)

      call catchem_state_get_species_long_name_at(state_ptr, int(index, c_int), c_buf)
      call c_to_f_string(c_buf, this%long_name)

      call catchem_state_get_species_desc_at(state_ptr, int(index, c_int), c_buf)
      call c_to_f_string(c_buf, this%description)

      this%is_gas = (catchem_state_is_species_gas(state_ptr, int(index, c_int)) /= 0)
      this%is_aerosol = (catchem_state_is_species_aerosol(state_ptr, int(index, c_int)) /= 0)
      this%is_tracer = (catchem_state_get_species_is_tracer(state_ptr, int(index, c_int)) /= 0)
      this%is_advected = (catchem_state_get_species_is_advected(state_ptr, int(index, c_int)) /= 0)
      this%is_drydep = (catchem_state_get_species_is_drydep(state_ptr, int(index, c_int)) /= 0)
      this%is_wetdep = (catchem_state_get_species_is_wetdep(state_ptr, int(index, c_int)) /= 0)
      this%is_photolysis = (catchem_state_get_species_is_photolysis(state_ptr, int(index, c_int)) /= 0)
      this%is_dust = (catchem_state_get_species_is_dust(state_ptr, int(index, c_int)) /= 0)
      this%is_seasalt = (catchem_state_get_species_is_seasalt(state_ptr, int(index, c_int)) /= 0)

      this%mw_g = real(catchem_state_get_species_mw(state_ptr, int(index, c_int)), fp)
      this%density = real(catchem_state_get_species_density(state_ptr, int(index, c_int)), fp)
      this%radius = real(catchem_state_get_species_radius(state_ptr, int(index, c_int)), fp)
      this%lower_radius = real(catchem_state_get_species_lower_radius(state_ptr, int(index, c_int)), fp)
      this%upper_radius = real(catchem_state_get_species_upper_radius(state_ptr, int(index, c_int)), fp)
      this%viscosity = real(catchem_state_get_species_viscosity(state_ptr, int(index, c_int)), fp)

      this%dd_f0 = real(catchem_state_get_species_dd_f0(state_ptr, int(index, c_int)), fp)
      this%dd_hstar = real(catchem_state_get_species_dd_hstar(state_ptr, int(index, c_int)), fp)
      this%dd_DvzAerSnow = real(catchem_state_get_species_dd_DvzAerSnow(state_ptr, int(index, c_int)), fp)
      this%dd_DvzMinVal_snow = real(catchem_state_get_species_dd_DvzMinVal_snow(state_ptr, int(index, c_int)), fp)
      this%dd_DvzMinVal_land = real(catchem_state_get_species_dd_DvzMinVal_land(state_ptr, int(index, c_int)), fp)

      this%henry_k0 = real(catchem_state_get_species_henry_k0(state_ptr, int(index, c_int)), fp)
      this%henry_cr = real(catchem_state_get_species_henry_cr(state_ptr, int(index, c_int)), fp)
      this%henry_pKa = real(catchem_state_get_species_henry_pKa(state_ptr, int(index, c_int)), fp)
      this%wd_retfactor = real(catchem_state_get_species_wd_retfactor(state_ptr, int(index, c_int)), fp)
      this%wd_LiqAndGas = (catchem_state_get_species_wd_LiqAndGas(state_ptr, int(index, c_int)) /= 0)
      this%wd_convfacI2G = real(catchem_state_get_species_wd_convfacI2G(state_ptr, int(index, c_int)), fp)
      call catchem_state_get_species_wd_rainouteff(state_ptr, int(index, c_int), r_eff)
      this%wd_rainouteff = real(r_eff, fp)
      this%wd_reevap_frac = real(catchem_state_get_species_wd_reevap_frac(state_ptr, int(index, c_int)), fp)

      this%t_chem_loss = real(catchem_state_get_species_t_chem_loss(state_ptr, int(index, c_int)), fp)
      this%BackgroundVV = real(catchem_state_get_species_BackgroundVV(state_ptr, int(index, c_int)), fp)

      call catchem_state_get_species_mie_name(state_ptr, int(index, c_int), c_buf)
      call c_to_f_string(c_buf, this%mie_name)

      this%species_index = index
      this%is_valid = .true.
   end subroutine populate_species_from_cpp

   subroutine species_init(this, species_name, long_name, molecular_weight, rc)
      class(SpeciesType), intent(inout) :: this
      character(len=*), intent(in) :: species_name, long_name
      real(fp), intent(in) :: molecular_weight
      integer, intent(out) :: rc
      rc = CC_SUCCESS
      this%short_name = trim(species_name)
      this%long_name = trim(long_name)
      this%mw_g = molecular_weight
      this%is_valid = .true.
   end subroutine species_init

   function species_validate(this, rc) result(is_val)
      class(SpeciesType), intent(in) :: this
      integer, intent(out) :: rc
      logical :: is_val
      rc = CC_SUCCESS
      is_val = this%is_valid
   end function species_validate

   subroutine species_cleanup(this, rc)
      class(SpeciesType), intent(inout) :: this
      integer, intent(out) :: rc
      rc = CC_SUCCESS
      nullify(this%conc)
      this%is_valid = .false.
   end subroutine species_cleanup

   subroutine species_set_concentration(this, concentration, grid_index, rc)
      class(SpeciesType), intent(inout) :: this
      real(fp), intent(in) :: concentration
      integer, intent(in), dimension(3) :: grid_index
      integer, intent(out) :: rc
      rc = CC_SUCCESS
      if (associated(this%conc)) then
         this%conc(grid_index(1), grid_index(2), grid_index(3)) = concentration
      end if
   end subroutine species_set_concentration

   function species_get_concentration(this, grid_index, rc) result(concentration)
      class(SpeciesType), intent(in) :: this
      integer, intent(in), dimension(3) :: grid_index
      integer, intent(out) :: rc
      real(fp) :: concentration
      rc = CC_SUCCESS
      concentration = 0.0_fp
      if (associated(this%conc)) then
         concentration = this%conc(grid_index(1), grid_index(2), grid_index(3))
      end if
   end function species_get_concentration

   subroutine species_copy(this, source, rc)
      class(SpeciesType), intent(inout) :: this
      class(SpeciesType), intent(in) :: source
      integer, intent(out) :: rc
      rc = CC_SUCCESS
      this%long_name = source%long_name
      this%short_name = source%short_name
      this%description = source%description
      this%is_gas = source%is_gas
      this%is_aerosol = source%is_aerosol
      this%is_tracer = source%is_tracer
      this%is_advected = source%is_advected
      this%is_drydep = source%is_drydep
      this%is_wetdep = source%is_wetdep
      this%is_photolysis = source%is_photolysis
      this%is_gocart_aero = source%is_gocart_aero
      this%is_dust = source%is_dust
      this%is_seasalt = source%is_seasalt
      this%mw_g = source%mw_g
      this%density = source%density
      this%radius = source%radius
      this%lower_radius = source%lower_radius
      this%upper_radius = source%upper_radius
      this%viscosity = source%viscosity
      this%dd_f0 = source%dd_f0
      this%dd_hstar = source%dd_hstar
      this%dd_DvzAerSnow = source%dd_DvzAerSnow
      this%dd_DvzMinVal_snow = source%dd_DvzMinVal_snow
      this%dd_DvzMinVal_land = source%dd_DvzMinVal_land
      this%henry_k0 = source%henry_k0
      this%henry_cr = source%henry_cr
      this%henry_pKa = source%henry_pKa
      this%wd_retfactor = source%wd_retfactor
      this%wd_LiqAndGas = source%wd_LiqAndGas
      this%wd_convfacI2G = source%wd_convfacI2G
      this%wd_rainouteff = source%wd_rainouteff
      this%wd_reevap_frac = source%wd_reevap_frac
      this%mie_name = source%mie_name
      this%t_chem_loss = source%t_chem_loss
      this%BackgroundVV = source%BackgroundVV
      this%species_index = source%species_index
      this%is_valid = source%is_valid
      this%conc => source%conc
   end subroutine species_copy

   subroutine species_print_info(this)
      class(SpeciesType), intent(in) :: this
      print *, 'Species properties: ', trim(this%short_name), ' MW=', this%mw_g, ' density=', this%density
   end subroutine species_print_info

   !-------------------
   ! SpeciesManagerType Procedures
   !-------------------
   subroutine species_manager_init(this, max_species, rc)
      class(SpeciesManagerType), intent(inout) :: this
      integer, intent(in) :: max_species
      integer, intent(out) :: rc
      rc = CC_SUCCESS
      if (allocated(this%species_db)) deallocate(this%species_db)
      allocate(this%species_db(max_species))
      this%num_species = 0
      this%is_initialized = .true.
   end subroutine species_manager_init

   subroutine species_manager_add_species(this, species, rc)
      class(SpeciesManagerType), intent(inout) :: this
      type(SpeciesType), intent(in) :: species
      integer, intent(out) :: rc
      rc = CC_SUCCESS
      if (this%num_species < size(this%species_db)) then
         this%num_species = this%num_species + 1
         call this%species_db(this%num_species)%copy(species, rc)
         this%species_db(this%num_species)%species_index = this%num_species
      else
         rc = CC_FAILURE
      end if
   end subroutine species_manager_add_species

   subroutine species_manager_find_species(this, species_name, species_index, rc)
      class(SpeciesManagerType), intent(in) :: this
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: species_index
      integer, intent(out) :: rc
      integer :: i
      rc = CC_FAILURE
      species_index = -1
      do i = 1, this%num_species
         if (trim(this%species_db(i)%short_name) == trim(species_name)) then
            species_index = i
            rc = CC_SUCCESS
            return
         end if
      end do
   end subroutine species_manager_find_species

   subroutine species_manager_validate_database(this, rc)
      class(SpeciesManagerType), intent(inout) :: this
      integer, intent(out) :: rc
      rc = CC_SUCCESS
   end subroutine species_manager_validate_database

   subroutine species_manager_load_from_cpp(this, state_mgr_ptr, rc)
      class(SpeciesManagerType), intent(inout) :: this
      type(c_ptr), intent(in) :: state_mgr_ptr
      integer, intent(out) :: rc
      integer :: n_spec, i

      rc = CC_SUCCESS
      if (.not. c_associated(state_mgr_ptr)) then
         rc = CC_FAILURE
         return
      end if

      n_spec = int(catchem_state_get_species_count(state_mgr_ptr))
      if (allocated(this%species_db)) deallocate(this%species_db)
      allocate(this%species_db(n_spec))

      do i = 1, n_spec
         call populate_species_from_cpp(this%species_db(i), state_mgr_ptr, i, rc)
         if (rc /= CC_SUCCESS) return
      end do
      this%num_species = n_spec
      this%is_initialized = .true.
   end subroutine species_manager_load_from_cpp

   subroutine species_manager_cleanup(this, rc)
      class(SpeciesManagerType), intent(inout) :: this
      integer, intent(out) :: rc
      rc = CC_SUCCESS
      if (allocated(this%species_db)) deallocate(this%species_db)
      this%num_species = 0
      this%is_initialized = .false.
   end subroutine species_manager_cleanup

   subroutine species_manager_print_database(this)
      class(SpeciesManagerType), intent(in) :: this
      print *, 'Species database count: ', this%num_species
   end subroutine species_manager_print_database

   !-------------------
   ! Standalone Helpers
   !-------------------
   function validate_species(species, rc) result(is_val)
      type(SpeciesType), intent(in) :: species
      integer, intent(out) :: rc
      logical :: is_val
      is_val = species%validate(rc)
   end function validate_species

   function find_species_by_name(species_db, num_species, species_name, rc) result(species_index)
      type(SpeciesType), intent(in) :: species_db(:)
      integer, intent(in) :: num_species
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      integer :: species_index
      integer :: i
      rc = CC_FAILURE
      species_index = -1
      do i = 1, num_species
         if (trim(species_db(i)%short_name) == trim(species_name)) then
            species_index = i
            rc = CC_SUCCESS
            return
         end if
      end do
   end function find_species_by_name

   subroutine create_species_database(species_mgr, rc)
      type(SpeciesManagerType), intent(out) :: species_mgr
      integer, intent(out) :: rc
      rc = CC_SUCCESS
      call species_mgr%init(100, rc)
   end subroutine create_species_database

end module species_mod
