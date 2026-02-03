

# File species\_mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**species\_mod.F90**](species__mod_8_f90.md)

[Go to the documentation of this file](species__mod_8_f90.md)


```Fortran


module species_mod
   use precision_mod
   use error_mod
   use, intrinsic :: ieee_arithmetic
   implicit none
   private

   public :: speciestype
   public :: speciesmanagertype
   public :: validate_species
   public :: find_species_by_name
   public :: create_species_database

   type :: speciestype

      ! Names
      character(len=30) :: long_name
      character(len=30) :: short_name
      character(len=50) :: description

      ! Logical switches
      logical :: is_gas
      logical :: is_aerosol
      logical :: is_tracer
      logical :: is_advected
      logical :: is_drydep
      logical :: is_wetdep
      logical :: is_photolysis
      logical :: is_gocart_aero
      logical :: is_dust
      logical :: is_seasalt

      ! Numerical properties
      real(kind=fp) :: mw_g                 
      real(kind=fp) :: density              
      real(kind=fp) :: radius               
      real(kind=fp) :: lower_radius         
      real(kind=fp) :: upper_radius         
      real(kind=fp) :: viscosity            

      ! used for dry deposition
      real(kind=fp) :: dd_f0                
      real(kind=fp) :: dd_hstar             
      real(kind=fp) :: dd_dvzaersnow        
      real(kind=fp) :: dd_dvzminval_snow    
      real(kind=fp) :: dd_dvzminval_land    

      ! used for wet deposition
      !real(kind=fp) :: radius_wet           !< mean molecular diameter in meters for wet conditions (use the same radius for both dry and wet deposition for now)
      real(kind=fp) :: henry_k0             
      real(kind=fp) :: henry_cr             
      real(kind=fp) :: henry_pka            
      real(kind=fp) :: wd_retfactor         
      logical       :: wd_LiqAndGas
      real(kind=fp) :: wd_convfaci2g        
      real(kind=fp) :: wd_rainouteff(3)     

      !used for settling
      character(len=30) :: mie_name

      ! Default background concentration
      real(kind=fp) :: backgroundvv        

      ! Indices
      integer :: species_index
      integer :: drydep_index
      integer :: photolysis_index
      integer :: gocart_aero_index

      ! Concentration
      real(kind=fp), POINTER :: conc(:,:,:)             

      ! Validation and status
      logical :: is_valid = .false.                     

   contains
      ! Enhanced methods with error handling
      procedure :: init => species_init
      procedure :: validate => species_validate
      procedure :: cleanup => species_cleanup
      procedure :: set_concentration => species_set_concentration
      procedure :: get_concentration => species_get_concentration
      procedure :: copy => species_copy
      procedure :: print_info => species_print_info

   end type speciestype

   type :: speciesmanagertype
      private

      type(SpeciesType), allocatable :: species_db(:)
      integer :: num_species = 0
      type(ErrorManagerType) :: error_mgr
      logical :: is_initialized = .false.              

   contains
      procedure :: init => species_manager_init
      procedure :: add_species => species_manager_add_species
      procedure :: find_species => species_manager_find_species
      procedure :: validate_database => species_manager_validate_database
      procedure :: load_from_file => species_manager_load_from_file
      procedure :: cleanup => species_manager_cleanup
      procedure :: print_database => species_manager_print_database

   end type speciesmanagertype

   !
   ! !DEFINED PARAMETERS:
   !
   !=========================================================================
   ! Missing species concentration value if not in restart file and special
   ! background value not defined
   !=========================================================================
   REAL(fp), PARAMETER, PUBLIC :: MISSING_VV  = 1.0e-20_fp 

contains

   subroutine species_init(this, species_name, long_name, molecular_weight, rc)
      implicit none
      class(SpeciesType), intent(inout) :: this
      character(len=*), intent(in) :: species_name
      character(len=*), intent(in) :: long_name
      real(fp), intent(in) :: molecular_weight
      integer, intent(out) :: rc

      ! Validate inputs
      if (len_trim(species_name) == 0) then
         rc = error_invalid_input
         return
      endif

      if (molecular_weight <= 0.0_fp) then
         rc = error_invalid_input
         return
      endif

      ! Initialize species properties
      this%short_name = trim(species_name)
      this%long_name = trim(long_name)
      this%description = ''  ! Initialize description
      this%mw_g = molecular_weight

      ! Set defaults
      this%is_gas = .true.
      this%is_aerosol = .false.
      this%is_tracer = .false.
      this%is_advected = .true.
      this%is_drydep = .false.
      this%is_wetdep = .false.  ! Initialize wet deposition flag
      this%is_photolysis = .false.
      this%is_gocart_aero = .false.
      this%is_dust = .false.
      this%is_seasalt = .false.

      this%density = 1000.0_fp  ! Default density
      this%radius = 1.0e-9_fp   ! Default radius
      this%lower_radius = 0.0_fp
      this%upper_radius = 0.0_fp
      this%viscosity = 1.0e-5_fp

      ! Initialize dry deposition properties
      this%dd_f0 = 0.0_fp
      this%dd_hstar = 0.0_fp
      this%dd_DvzAerSnow = 0.0_fp
      this%dd_DvzMinVal_snow = 0.0_fp
      this%dd_DvzMinVal_land = 0.0_fp

      ! Initialize wet deposition properties
      this%henry_k0 = 0.0_fp
      this%henry_cr = 0.0_fp
      this%henry_pKa = 0.0_fp
      this%wd_retfactor = 0.0_fp
      this%wd_LiqAndGas = .false.
      this%wd_convfacI2G = 0.0_fp
      this%wd_rainouteff(:) = 0.0_fp

      this%BackgroundVV = missing_vv
      this%mie_name = ''  ! Initialize Mie name to empty

      this%species_index = -1
      this%drydep_index = -1
      this%photolysis_index = -1
      this%gocart_aero_index = -1

      this%is_valid = .true.
      rc = cc_success

   end subroutine species_init

   function species_validate(this, rc) result(is_valid)
      implicit none
      class(SpeciesType), intent(in) :: this
      integer, intent(out) :: rc
      logical :: is_valid

      is_valid = .true.
      rc = cc_success

      ! Check basic properties
      if (len_trim(this%short_name) == 0) then
         is_valid = .false.
         rc = error_invalid_input
         return
      endif

      if (this%mw_g <= 0.0_fp) then
         is_valid = .false.
         rc = error_invalid_input
         return
      endif

      ! Check logical consistency
      if (this%is_gas .and. this%is_aerosol) then
         is_valid = .false.
         rc = error_state_inconsistency
         return
      endif

      if (this%is_dust .and. .not. this%is_aerosol) then
         is_valid = .false.
         rc = error_state_inconsistency
         return
      endif

      if (this%is_seasalt .and. .not. this%is_aerosol) then
         is_valid = .false.
         rc = error_state_inconsistency
         return
      endif

      ! Check physical properties for aerosols
      if (this%is_aerosol) then
         if (this%density <= 0.0_fp) then
            is_valid = .false.
            rc = error_invalid_input
            return
         endif

         if (this%radius <= 0.0_fp) then
            is_valid = .false.
            rc = error_invalid_input
            return
         endif
      endif

      ! Check concentration array: all values must be positive and finite
      if (associated(this%conc)) then
         if (any(this%conc < 0.0_fp)) then
            is_valid = .false.
            rc = error_invalid_input
            return
         endif
         if (any(.not. ieee_is_finite(this%conc))) then
            is_valid = .false.
            rc = error_invalid_input
            return
         endif
      endif

   end function species_validate

   subroutine species_set_concentration(this, concentration, grid_index, rc)
      implicit none
      class(SpeciesType), intent(inout) :: this
      real(fp), intent(in) :: concentration
      integer, intent(in), dimension(3) :: grid_index
      integer, intent(out) :: rc

      ! Check if concentration array is allocated
      if (.not. associated(this%conc)) then
         rc = error_state_inconsistency
         return
      endif

      ! Check bounds
      if (any(grid_index < 1) .or. any(grid_index > shape(this%conc))) then
         rc = error_bounds_check
         return
      endif

      ! Check for valid concentration
      if (concentration < 0.0_fp) then
         rc = error_invalid_input
         return
      endif

      this%conc(grid_index(1), grid_index(2), grid_index(3)) = concentration
      rc = cc_success

   end subroutine species_set_concentration

   subroutine species_cleanup(this, rc)
      implicit none
      class(SpeciesType), intent(inout) :: this
      integer, intent(out) :: rc

      if (associated(this%conc)) deallocate(this%conc)

      this%description = ''  ! Clear description
      this%mie_name = ''  ! Clear Mie name
      this%is_valid = .false.
      rc = cc_success

   end subroutine species_cleanup

   subroutine species_manager_init(this, max_species, rc)
      implicit none
      class(SpeciesManagerType), intent(inout) :: this
      integer, intent(in) :: max_species
      integer, intent(out) :: rc

      if (max_species <= 0) then
         rc = error_invalid_input
         return
      endif

      allocate(this%species_db(max_species), stat=rc)
      if (rc /= 0) then
         rc = error_memory_allocation
         return
      endif

      call this%error_mgr%init(verbose=.true.)
      this%num_species = 0
      this%is_initialized = .true.
      rc = cc_success

   end subroutine species_manager_init

   subroutine species_manager_find_species(this, species_name, species_index, rc)
      implicit none
      class(SpeciesManagerType), intent(in) :: this
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: species_index
      integer, intent(out) :: rc

      integer :: i

      species_index = -1
      rc = cc_success

      if (.not. this%is_initialized) then
         rc = error_state_inconsistency
         return
      endif

      do i = 1, this%num_species
         if (trim(this%species_db(i)%short_name) == trim(species_name)) then
            species_index = i
            return
         endif
      enddo

      ! Species not found
      rc = error_invalid_input

   end subroutine species_manager_find_species

   function validate_species(species, rc) result(is_valid)
      implicit none
      type(SpeciesType), intent(in) :: species
      integer, intent(out) :: rc
      logical :: is_valid

      is_valid = species%validate(rc)

   end function validate_species

   function find_species_by_name(species_db, num_species, species_name, rc) result(species_index)
      implicit none
      type(SpeciesType), intent(in) :: species_db(:)
      integer, intent(in) :: num_species
      character(len=*), intent(in) :: species_name
      integer, intent(out) :: rc
      integer :: species_index

      integer :: i

      species_index = -1
      rc = cc_success

      do i = 1, num_species
         if (trim(species_db(i)%short_name) == trim(species_name)) then
            species_index = i
            return
         endif
      enddo

      rc = error_invalid_input

   end function find_species_by_name

   subroutine create_species_database(species_mgr, rc)
      implicit none
      type(SpeciesManagerType), intent(out) :: species_mgr
      integer, intent(out) :: rc

      call species_mgr%init(100, rc)  ! Support up to 100 species
      if (rc /= cc_success) return

      ! Add some common species (this would be expanded)
      ! This is a placeholder implementation
      rc = cc_success

   end subroutine create_species_database

   function species_get_concentration(this, grid_index, rc) result(concentration)
      implicit none
      class(SpeciesType), intent(in) :: this
      integer, intent(in), dimension(3) :: grid_index
      integer, intent(out) :: rc
      real(fp) :: concentration

      rc = cc_success
      concentration = 0.0_fp

      ! Check if concentration array is allocated
      if (.not. associated(this%conc)) then
         rc = error_memory_allocation
         return
      endif

      ! Check bounds
      if (any(grid_index < 1) .or. any(grid_index > shape(this%conc))) then
         rc = error_bounds_check
         return
      endif

      concentration = this%conc(grid_index(1), grid_index(2), grid_index(3))

   end function species_get_concentration

   subroutine species_copy(this, source, rc)
      implicit none
      class(SpeciesType), intent(inout) :: this
      class(SpeciesType), intent(in) :: source
      integer, intent(out) :: rc

      rc = cc_success

      ! Copy all scalar properties
      this%short_name = source%short_name
      this%long_name = source%long_name
      this%description = source%description

      ! Copy logical switches
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

      ! Copy numerical properties
      this%mw_g = source%mw_g
      this%density = source%density
      this%radius = source%radius
      this%lower_radius = source%lower_radius
      this%upper_radius = source%upper_radius
      this%viscosity = source%viscosity

      ! Copy dry deposition properties
      this%dd_f0 = source%dd_f0
      this%dd_hstar = source%dd_hstar
      this%dd_DvzAerSnow = source%dd_DvzAerSnow
      this%dd_DvzMinVal_snow = source%dd_DvzMinVal_snow
      this%dd_DvzMinVal_land = source%dd_DvzMinVal_land

      ! Copy wet deposition properties
      this%henry_k0 = source%henry_k0
      this%henry_cr = source%henry_cr
      this%henry_pKa = source%henry_pKa
      this%wd_retfactor = source%wd_retfactor
      this%wd_LiqAndGas = source%wd_LiqAndGas
      this%wd_convfacI2G = source%wd_convfacI2G
      this%wd_rainouteff = source%wd_rainouteff

      this%BackgroundVV = source%BackgroundVV
      this%mie_name = source%mie_name

      ! Copy indices
      this%species_index = source%species_index
      this%drydep_index = source%drydep_index
      this%photolysis_index = source%photolysis_index
      this%gocart_aero_index = source%gocart_aero_index

      ! Copy concentration array if allocated
      if (associated(source%conc)) then
         if (associated(this%conc)) deallocate(this%conc)
         allocate(this%conc(size(source%conc,1), size(source%conc,2), size(source%conc,3)), stat=rc)
         if (rc /= 0) then
            rc = error_memory_allocation
            return
         endif
         this%conc = source%conc
      endif

      this%is_valid = source%is_valid
      rc = cc_success

   end subroutine species_copy

   subroutine species_print_info(this)
      implicit none
      class(SpeciesType), intent(in) :: this

      write(*, '(A)') '=== Species Information ==='
      write(*, '(A,A)') 'Short name: ', trim(this%short_name)
      write(*, '(A,A)') 'Long name:  ', trim(this%long_name)
      write(*, '(A,A)') 'Description: ', trim(this%description)
      write(*, '(A,F12.6)') 'Molecular weight [g/mol]: ', this%mw_g
      write(*, '(A,F12.3)') ³'Density [kg/m]: ', this%density
      write(*, '(A,E12.5)') 'Radius [m]: ', this%radius
      write(*, '(A,L1)') 'Is gas: ', this%is_gas
      write(*, '(A,L1)') 'Is aerosol: ', this%is_aerosol
      write(*, '(A,L1)') 'Is tracer: ', this%is_tracer
      write(*, '(A,L1)') 'Is advected: ', this%is_advected
      write(*, '(A,L1)') 'Undergoes dry deposition: ', this%is_drydep
      write(*, '(A,L1)') 'Undergoes wet deposition: ', this%is_wetdep
      write(*, '(A,L1)') 'Undergoes photolysis: ', this%is_photolysis
      write(*, '(A,L1)') 'Is GOCART aerosol: ', this%is_gocart_aero
      write(*, '(A,L1)') 'Is dust: ', this%is_dust
      write(*, '(A,L1)') 'Is seasalt: ', this%is_seasalt
      write(*, '(A,A)') 'Mie data name: ', trim(this%mie_name)
      write(*, '(A,E12.5)') 'Background concentration [v/v]: ', this%BackgroundVV
      write(*, '(A,I0)') 'Species index: ', this%species_index
      if (associated(this%conc)) then
         write(*, '(A,I0)') 'Concentration grid size: ', size(this%conc)
      else
         write(*, '(A)') 'Concentration: Not allocated'
      endif
      write(*, '(A,L1)') 'Valid: ', this%is_valid
      write(*, '(A)') '=========================='

   end subroutine species_print_info

   subroutine species_manager_add_species(this, species, rc)
      implicit none
      class(SpeciesManagerType), intent(inout) :: this
      type(SpeciesType), intent(in) :: species
      integer, intent(out) :: rc

      rc = cc_success

      if (.not. this%is_initialized) then
         rc = error_process_initialization
         return
      endif

      if (.not. allocated(this%species_db)) then
         rc = error_memory_allocation
         return
      endif

      if (this%num_species >= size(this%species_db)) then
         rc = error_bounds_check
         return
      endif

      ! Add species and increment counter
      this%num_species = this%num_species + 1
      call this%species_db(this%num_species)%copy(species, rc)
      if (rc /= cc_success) return

      ! Set the species index
      this%species_db(this%num_species)%species_index = this%num_species

   end subroutine species_manager_add_species

   subroutine species_manager_validate_database(this, rc)
      implicit none
      class(SpeciesManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i, local_rc
      character(len=100) :: error_msg
      logical :: is_valid

      rc = cc_success

      if (.not. this%is_initialized) then
         rc = error_process_initialization
         return
      endif

      ! Validate each species
      do i = 1, this%num_species
         is_valid = this%species_db(i)%validate(local_rc)
         if (local_rc /= cc_success .or. .not. is_valid) then
            write(error_msg, '(A,I0)') 'Species validation failed for species index ', i
            call this%error_mgr%report_error(local_rc, error_msg, &
               rc, 'species_manager_validate_database')
            if (rc /= cc_success) return
         endif
      enddo

   end subroutine species_manager_validate_database

   subroutine species_manager_load_from_file(this, filename, rc)
      implicit none
      class(SpeciesManagerType), intent(inout) :: this
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      rc = cc_success

      if (.not. this%is_initialized) then
         rc = error_process_initialization
         return
      endif

      ! Placeholder implementation
      ! In a full implementation, this would:
      ! 1. Open and parse the configuration file
      ! 2. Create SpeciesType objects from file data
      ! 3. Add them to the database using add_species
      ! 4. Validate the loaded database

      write(*, '(A,A)') 'INFO: Loading species from file: ', trim(filename)
      write(*, '(A)') 'WARNING: species_manager_load_from_file is a placeholder implementation'

   end subroutine species_manager_load_from_file

   subroutine species_manager_cleanup(this, rc)
      implicit none
      class(SpeciesManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i, local_rc

      rc = cc_success

      ! Clean up each species
      if (allocated(this%species_db)) then
         do i = 1, this%num_species
            call this%species_db(i)%cleanup(local_rc)
            ! Continue even if individual cleanup fails
         enddo
         deallocate(this%species_db)
      endif

      this%num_species = 0
      this%is_initialized = .false.

   end subroutine species_manager_cleanup

   subroutine species_manager_print_database(this)
      implicit none
      class(SpeciesManagerType), intent(in) :: this

      integer :: i

      write(*, '(A)') '=== Species Database Summary ==='
      write(*, '(A,L1)') 'Initialized: ', this%is_initialized
      write(*, '(A,I0)') 'Number of species: ', this%num_species

      if (allocated(this%species_db)) then
         write(*, '(A,I0)') 'Database capacity: ', size(this%species_db)

         write(*, '(A)') 'Species list:'
         do i = 1, this%num_species
            write(*, '(I4,A,A,A,A)') i, ': ', trim(this%species_db(i)%short_name), &
               ' (', trim(this%species_db(i)%long_name), ')'
         enddo
      else
         write(*, '(A)') 'Database: Not allocated'
      endif

      write(*, '(A)') '==============================='

   end subroutine species_manager_print_database

   ! function get_name(this) result(species_name)
   !    character(len=30) :: species_name

   !    species_name = this%short_name
   ! end function get_name

   ! function get_atomic_number(this) result(atomic_num)
   !    class(Species), intent(in) :: this
   !    integer :: atomic_num

   !    atomic_num = this%atomic_number
   ! end function get_atomic_number

end module species_mod
```


