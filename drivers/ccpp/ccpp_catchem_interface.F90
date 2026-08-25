!> \file ccpp_catchem_interface.F90
!! \brief CCPP interface module for CATChem integration with dynamic constituents
!!
module ccpp_catchem_interface

  use iso_c_binding,             only: c_loc, c_null_char, c_double, c_ptr, c_associated, c_int, c_char
  use machine,                   only: kind_phys
  use CATChem_API,               only: CATChem_Model
  use Error_Mod,                 only: CC_SUCCESS, CC_FAILURE
  use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t, ccpp_constituent_prop_ptr_t
  use precision_mod,             only: fp

  implicit none

  private

  public :: ccpp_catchem_interface_register, &
            ccpp_catchem_interface_init, &
            ccpp_catchem_interface_run, &
            ccpp_catchem_interface_finalize

  ! Module-scope wrapper mapping directly to the C++ core
  type(CATChem_Model), save :: cc_model

  ! Saved dynamic index mapper array to store host index for each species
  integer, allocatable, save :: catchem_indices_constituent_props(:)

  ! Linkable BIND(C) external declarations to avoid circular dependencies
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

     real(c_double) function catchem_state_get_species_mw(state_ptr, index) bind(C, name="catchem_state_get_species_mw")
        import :: c_ptr, c_int, c_double
        type(c_ptr), value :: state_ptr
        integer(c_int), value :: index
     end function

     integer(c_int) function catchem_state_get_species_is_advected(state_ptr, index) bind(C, name="catchem_state_get_species_is_advected")
        import :: c_ptr, c_int
        type(c_ptr), value :: state_ptr
        integer(c_int), value :: index
     end function

     integer(c_int) function catchem_state_is_species_gas(state_ptr, index) bind(C, name="catchem_state_is_species_gas")
        import :: c_ptr, c_int
        type(c_ptr), value :: state_ptr
        integer(c_int), value :: index
     end function

     subroutine catchem_state_bind_met_3d(state_ptr, name, ptr) bind(C, name="catchem_state_bind_met_3d")
        import :: c_ptr, c_char
        type(c_ptr), value :: state_ptr
        character(kind=c_char), intent(in) :: name(*)
        type(c_ptr), value :: ptr
     end subroutine

     subroutine catchem_state_bind_met_2d(state_ptr, name, ptr) bind(C, name="catchem_state_bind_met_2d")
        import :: c_ptr, c_char
        type(c_ptr), value :: state_ptr
        character(kind=c_char), intent(in) :: name(*)
        type(c_ptr), value :: ptr
     end subroutine

     subroutine catchem_state_bind_unified_chemistry(state_ptr, ptr) bind(C, name="catchem_state_bind_unified_chemistry")
        import :: c_ptr
        type(c_ptr), value :: state_ptr
        type(c_ptr), value :: ptr
     end subroutine
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

   !> \brief Register dynamic constituent properties with CCPP
   !!
   !! \details This subroutine reads the chemical species configuration database using the path
   !! retrieved from the `CATCHEM_CONFIG` environment variable (falling back to standard path),
   !! queries the total species count from the modern C++ Core StateManager, and dynamically
   !! instantiates the CCPP framework metadata registry (`ccpp_constituent_properties_t`) for
   !! each active species (including name, molecular weight, and advection parameters).
   !!
   !! \param[out] constituent_props Array of dynamically registered CCPP constituent properties
   !! \param[out] errmsg Error message string set if registration fails
   !! \param[out] errflg Error code flag (0=success, non-zero=failure)
   subroutine ccpp_catchem_interface_register(constituent_props, errmsg, errflg)
      implicit none

      type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituent_props(:)
      character(len=*),                                 intent(out) :: errmsg
      integer,                                          intent(out) :: errflg

      ! Local variables
      character(len=512)     :: config_path
      character(kind=c_char) :: c_buf(128)
      character(len=64)      :: f_name
      real(kind_phys)        :: molar_mass_g_mol, molar_mass_kg_mol
      logical                :: is_advected, is_gas
      integer                :: n_spec, i
      type(c_ptr)            :: state_mgr

      errmsg = ''
      errflg = 0

      ! 1. Locate YAML configuration path via environment variable
      call get_environment_variable("CATCHEM_CONFIG", config_path)
      if (trim(config_path) == "") then
         config_path = "./tests/Configs/Default/CATChem_config.yml"
      end if

      ! 2. Lightweight core initialization to load species configuration (spatial bounds are dummy 1, 1, 1 during register)
      call cc_model%initialize(config_path, 1, 1, 1, rc=errflg)
      if (errflg /= CC_SUCCESS) then
         errmsg = 'CATChem Register: Failed to load species configuration.'
         return
      end if

      state_mgr = cc_model%get_state_manager()
      n_spec = int(catchem_state_get_species_count(state_mgr))

      ! 3. Allocate and populate constituent_props
      allocate(constituent_props(n_spec), stat=errflg)
      if (errflg /= 0) {
         errmsg = 'CATChem Register: Memory allocation failure.'
         return
      }

      do i = 1, n_spec
         call catchem_state_get_species_name_at(state_mgr, int(i, c_int), c_buf)
         call c_to_f_string(c_buf, f_name)

         molar_mass_g_mol = real(catchem_state_get_species_mw(state_mgr, int(i, c_int)), kind_phys)
         molar_mass_kg_mol = molar_mass_g_mol * 1.0e-3_kind_phys
         is_advected = (catchem_state_get_species_is_advected(state_mgr, int(i, c_int)) == 1)
         is_gas = (catchem_state_is_species_gas(state_mgr, int(i, c_int)) == 1)

         if (is_gas) then
            call constituent_props(i)%instantiate( &
               std_name      = trim(f_name), &
               long_name     = trim(f_name), &
               diag_name     = trim(f_name), &
               units         = 'kg kg-1', &
               vertical_dim  = 'vertical_layer_dimension', &
               default_value = 0.0_kind_phys, &
               min_value     = 0.0_kind_phys, &
               molar_mass    = molar_mass_kg_mol, &
               advected      = is_advected, &
               errcode       = errflg, &
               errmsg        = errmsg )
         else
            ! Aerosol species do not have a constant gas-phase molar mass, so we omit the molar_mass argument
            call constituent_props(i)%instantiate( &
               std_name      = trim(f_name), &
               long_name     = trim(f_name), &
               diag_name     = trim(f_name), &
               units         = 'kg kg-1', &
               vertical_dim  = 'vertical_layer_dimension', &
               default_value = 0.0_kind_phys, &
               min_value     = 0.0_kind_phys, &
               advected      = is_advected, &
               errcode       = errflg, &
               errmsg        = errmsg )
         end if
         if (errflg /= 0) return
      end do

      ! 4. Clean up lightweight setup so model can be fully initialized with real grid sizes during init phase
      call cc_model%finalize(errflg)

   end subroutine ccpp_catchem_interface_register

   !> \brief Initialize the CATChem CCPP interface and resolve host-mapped indices
   !!
   !! \details Fully configures the high-performance C++ Core manager dynamically using the
   !! host-provided physical, meteorological, and land/soil grid dimensions (such as horizontal columns,
   !! vertical layers, soil layers, soil categories, and vegetation categories). It then resolves the
   !! global host-model runtime index mapping for each registered chemical species using CCPP's standard
   !! index retrieval utility `ccpp_const_get_idx`.
   !!
   !! \param[in]  im Horizontal loop extent / dimension of the chunk
   !! \param[in]  kte Vertical layer dimension / number of layers
   !! \param[in]  nsoil Vertical dimension of soil layers
   !! \param[in]  nlndcat Number of vegetation categories
   !! \param[in]  nsoilcat Number of soil categories
   !! \param[in]  do_catchem Logical coupling flag controlling activation of CATChem
   !! \param[in]  catchem_configfile_in Path to the CATChem YAML configuration file
   !! \param[in]  constituent_props_ptr CCPP constituent properties pointer array
   !! \param[out] errmsg Error message string set if initialization fails
   !! \param[out] errflg Error code flag (0=success, non-zero=failure)
   subroutine ccpp_catchem_interface_init(im, kte, nsoil, nlndcat, nsoilcat, &
                                          do_catchem, catchem_configfile_in, &
                                          constituent_props_ptr, errmsg, errflg)
      use ccpp_const_utils, only: ccpp_const_get_idx
      implicit none

      integer,                           intent(in)  :: im
      integer,                           intent(in)  :: kte
      integer,                           intent(in)  :: nsoil
      integer,                           intent(in)  :: nlndcat
      integer,                           intent(in)  :: nsoilcat
      logical,                           intent(in)  :: do_catchem
      character(len=*),                  intent(in)  :: catchem_configfile_in
      type(ccpp_constituent_prop_ptr_t), intent(in)  :: constituent_props_ptr(:)
      character(len=*),                  intent(out) :: errmsg
      integer,                           intent(out) :: errflg

      character(kind=c_char) :: c_buf(128)
      character(len=64)      :: f_name
      integer                :: n_spec, i
      type(c_ptr)            :: state_mgr

      errmsg = ''
      errflg = 0

      if (.not. do_catchem) return

      ! 1. Fully initialize C++ Core manager dynamically with host grid and soil bounds
      call cc_model%initialize(catchem_configfile_in, im, 1, kte, &
                               nsoil=nsoil, nsoiltype=nsoilcat, nsurftype=nlndcat, &
                               rc=errflg)
      if (errflg /= CC_SUCCESS) then
         errmsg = 'CATChem Init: Failed to initialize C++ Core via cc_model'
         return
      end if

      state_mgr = cc_model%get_state_manager()
      n_spec = int(catchem_state_get_species_count(state_mgr))

      ! 2. Resolve index mappings from CCPP global constituent properties pointer
      if (allocated(catchem_indices_constituent_props)) deallocate(catchem_indices_constituent_props)
      allocate(catchem_indices_constituent_props(n_spec))

      do i = 1, n_spec
         call catchem_state_get_species_name_at(state_mgr, int(i, c_int), c_buf)
         call c_to_f_string(c_buf, f_name)

         call ccpp_const_get_idx(constituent_props_ptr, trim(f_name), &
                                 catchem_indices_constituent_props(i), errmsg, errflg)
         if (errflg /= 0) then
            errmsg = "CATChem Init: Missing required tracer: " // trim(f_name)
            return
         end if
      end do

   end subroutine ccpp_catchem_interface_init

   !> \brief Finalize the CATChem CCPP interface and release allocated memory
   !!
   !! \details Shuts down the modern C++ Core orchestrator, releases all associated memory allocations
   !! inside the unmanaged Kokkos and C++ Views, and deallocates the dynamic host-constituent index mapping arrays.
   !!
   !! \param[in]  do_catchem Logical coupling flag controlling activation of CATChem
   !! \param[out] errmsg Error message string set if finalization fails
   !! \param[out] errflg Error code flag (0=success, non-zero=failure)
   subroutine ccpp_catchem_interface_finalize(do_catchem, errmsg, errflg)
     implicit none

     logical, intent(in) :: do_catchem

     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     errmsg = ''
     errflg = 0

     if (.not. do_catchem) return

     call cc_model%finalize(errflg)
     if (errflg /= CC_SUCCESS) then
         errmsg = 'CATChem Finalize: Error finalising C++ Core via cc_model'
     end if

     if (allocated(catchem_indices_constituent_props)) deallocate(catchem_indices_constituent_props)

  end subroutine ccpp_catchem_interface_finalize

  !> \brief Execute CATChem chemistry, deposition, and emission calculations
  !!
  !! \details This is the primary execution entry point of the CCPP driver. It dynamically extracts
  !! the active chemical tracer mixing ratios from the non-contiguous global host array `constituents`
  !! into a local contiguous buffer, binds it alongside standard meteorological, surface, vegetative,
  !! soil, and climatological dust parameters directly to the C++ Core via zero-copy C-API pointers,
  !! runs the thread-safe parallel Kokkos chemistry solvers, and synchronizes the updated tracer state
  !! back to the host model in-place.
  !!
  !! \param[in]    im Horizontal loop extent / dimension of the chunk
  !! \param[in]    kte Vertical layer dimension / number of layers
  !! \param[in]    kme Vertical interface dimension / number of layer edges
  !! \param[in]    garea Grid cell area (m2)
  !! \param[in]    nsoil Vertical dimension of soil layers
  !! \param[in]    nlndcat Number of vegetation categories
  !! \param[in]    nsoilcat Number of soil categories
  !! \param[in]    lat Grid cell latitude coordinates (degrees)
  !! \param[in]    lon Grid cell longitude coordinates (degrees)
  !! \param[in]    do_catchem Logical coupling flag controlling activation of CATChem
  !! \param[in]    dt Physics timestep length (s)
  !! \param[in]    jdate Forecast start date and time array
  !! \param[in]    xcosz Cosine of solar zenith angle
  !! \param[in]    lwi Land-water-ice mask index (0=water, 1=land, 2=ice)
  !! \param[in]    frlanduse Fractional vegetation cover by category
  !! \param[in]    gvf Green vegetation fraction (0-1)
  !! \param[in]    seaicefrac Sea ice concentration fraction
  !! \param[in]    oceanfrac Fraction of ocean
  !! \param[in]    lakefrac Fraction of lake
  !! \param[in]    landfrac Fraction of land
  !! \param[in]    stype Soil type index
  !! \param[in]    vtype Vegetation type index
  !! \param[in]    snowdepth Physical snow depth (m)
  !! \param[in]    frsnow Snow cover fraction
  !! \param[in]    lai Leaf area index (m2 m-2)
  !! \param[in]    frsoil Fractional soil type by category
  !! \param[in]    pores Dry soil porosity constants
  !! \param[in]    resid Residual soil moisture constants
  !! \param[in]    ustar friction velocity (m s-1)
  !! \param[in]    u10m 10m zonal wind velocity (m s-1)
  !! \param[in]    v10m 10m meridional wind velocity (m s-1)
  !! \param[in]    tskin Ground skin temperature (K)
  !! \param[in]    ts Surface temperature for fluxes (K)
  !! \param[in]    hf2d Sensible heat flux (W m-2)
  !! \param[in]    lf2d Latent heat flux (W m-2)
  !! \param[in]    znt Roughness length for momentum (m)
  !! \param[in]    prsfc Surface air pressure (Pa)
  !! \param[in]    pblh Planetary boundary layer height (m)
  !! \param[in]    dswsfc Downward shortwave flux at surface (W m-2)
  !! \param[in]    nirbmdi Downward near-infrared direct flux (W m-2)
  !! \param[in]    nirdfdi Downward near-infrared diffuse flux (W m-2)
  !! \param[in]    visbmdi Downward visible direct flux (W m-2)
  !! \param[in]    visdfdi Downward visible diffuse flux (W m-2)
  !! \param[in]    sfc_alb_nir_dir Surface albedo due to NIR direct
  !! \param[in]    sfc_alb_nir_dif Surface albedo due to NIR diffuse
  !! \param[in]    sfc_alb_uvvis_dir Surface albedo due to UV/VIS direct
  !! \param[in]    sfc_alb_uvvis_dif Surface albedo due to UV/VIS diffuse
  !! \param[in]    soilmoist Soil moisture fraction (m3 m-3)
  !! \param[in]    pr3d Air pressure at interfaces (Pa)
  !! \param[in]    phl3d Geopotential height at interfaces (m)
  !! \param[in]    prl3d Mid-point air pressure (Pa)
  !! \param[in]    tk3d Mid-point air temperature (K)
  !! \param[in]    q3d Specific humidity / water vapor mixing ratio (kg kg-1)
  !! \param[in]    us3d Eastward wind velocity (m s-1)
  !! \param[in]    vs3d Northward wind velocity (m s-1)
  !! \param[in]    rh Relative humidity (percent)
  !! \param[in]    delp Pressure thickness (Pa)
  !! \param[in]    airden Dry air density (kg m-3)
  !! \param[in]    pfl_lsan Cloud liquid water mixing ratio (kg kg-1)
  !! \param[in]    pfl_isan Cloud ice water mixing ratio (kg kg-1)
  !! \param[in]    rain_cpl Precipitation rate (kg m-2 s-1)
  !! \param[in]    cldf Cloud area fraction
  !! \param[in]    clayfrac Soil clay fraction (0-1)
  !! \param[in]    rdrag Wind drag partitioning parameter (unitless)
  !! \param[in]    sandfrac Soil sand fraction (0-1)
  !! \param[in]    ssm Sediment supply map parameter (unitless)
  !! \param[in]    ustar_threshold Threshold friction velocity for erosion (m s-1)
  !! \param[in]    constituent_props CCPP constituent properties pointer array
  !! \param[inout] constituents Global 3D CCPP unified tracer array [cols, levels, species]
  !! \param[out]   errmsg Error message string set if execution fails
  !! \param[out]   errflg Error code flag (0=success, non-zero=failure)
  subroutine ccpp_catchem_interface_run(im, kte, kme, garea, nsoil, nlndcat, nsoilcat, &
     lat, lon, &
     do_catchem, &
     dt, jdate, &
     xcosz, &
     lwi, frlanduse, gvf, seaicefrac, oceanfrac, lakefrac, landfrac, &
     stype, vtype, snowdepth, frsnow, lai, frsoil, pores, resid, &
     ustar, u10m, v10m, tskin, ts, hf2d, lf2d, znt, prsfc, pblh, &
     dswsfc, nirbmdi, nirdfdi, visbmdi, visdfdi, &
     sfc_alb_nir_dir, sfc_alb_nir_dif, sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, &
     soilmoist, pr3d, phl3d, prl3d, tk3d, q3d, us3d, vs3d, rh, &
     delp, airden, pfl_lsan, pfl_isan, &
     rain_cpl, cldf, &
     clayfrac, rdrag, sandfrac, ssm, ustar_threshold, &
     constituent_props, constituents, &
     errmsg, errflg)

     implicit none

     integer, intent(in) :: im
     integer, intent(in) :: kte
     integer, intent(in) :: kme
     integer, intent(in) :: nsoil
     integer, intent(in) :: nlndcat
     integer, intent(in) :: nsoilcat
     real(kind_phys), dimension(im), intent(in), target :: garea
     real(kind_phys), dimension(im), intent(in), target :: lat
     real(kind_phys), dimension(im), intent(in), target :: lon
     real(kind_phys), dimension(im), intent(in), target :: xcosz

     real(kind_phys), intent(in) :: dt
     integer, intent(in) :: jdate(8)

     logical, intent(in) :: do_catchem

     integer, dimension(im), intent(in), target                :: lwi
     integer, dimension(im), intent(in), target                :: stype
     integer, dimension(im), intent(in), target                :: vtype

     real(kind_phys), dimension(im, nlndcat), intent(in), target :: frlanduse
     real(kind_phys), dimension(im, nsoilcat), intent(in), target :: frsoil
     real(kind_phys), dimension(30), intent(in), target        :: pores
     real(kind_phys), dimension(30), intent(in), target        :: resid
     real(kind_phys), dimension(im), intent(in), target        :: seaicefrac
     real(kind_phys), dimension(im), intent(in), target        :: oceanfrac
     real(kind_phys), dimension(im), intent(in), target        :: frsnow
     real(kind_phys), dimension(im), intent(in), target        :: lakefrac
     real(kind_phys), dimension(im), intent(in), target        :: landfrac
     real(kind_phys), dimension(im), intent(in), target        :: gvf
     real(kind_phys), dimension(im), intent(in), target        :: lai

     real(kind_phys), dimension(im, nsoil), intent(in), target :: soilmoist
     real(kind_phys), dimension(im, snowdepth) :: snowdepth_not_used ! just matching names
     real(kind_phys), dimension(im), intent(in), target        :: snowdepth
     real(kind_phys), dimension(im), intent(in), target        :: prsfc
     real(kind_phys), dimension(im), intent(in), target        :: pblh
     real(kind_phys), dimension(im), intent(in), target        :: u10m
     real(kind_phys), dimension(im), intent(in), target        :: v10m
     real(kind_phys), dimension(im), intent(in), target        :: ustar
     real(kind_phys), dimension(im), intent(in), target        :: tskin
     real(kind_phys), dimension(im), intent(in), target        :: ts
     real(kind_phys), dimension(im), intent(in), target        :: hf2d
     real(kind_phys), dimension(im), intent(in), target        :: lf2d
     real(kind_phys), dimension(im), intent(in), target        :: znt
     real(kind_phys), dimension(im), intent(in), target        :: dswsfc
     real(kind_phys), dimension(im), intent(in), target        :: sfc_alb_nir_dir
     real(kind_phys), dimension(im), intent(in), target        :: sfc_alb_nir_dif
     real(kind_phys), dimension(im), intent(in), target        :: sfc_alb_uvvis_dir
     real(kind_phys), dimension(im), intent(in), target        :: sfc_alb_uvvis_dif
     real(kind_phys), dimension(im), intent(in), target        :: nirbmdi
     real(kind_phys), dimension(im), intent(in), target        :: nirdfdi
     real(kind_phys), dimension(im), intent(in), target        :: visbmdi
     real(kind_phys), dimension(im), intent(in), target        :: visdfdi

     real(kind_phys), dimension(im, kme), intent(in), target :: pr3d
     real(kind_phys), dimension(im, kte), intent(in), target :: prl3d
     real(kind_phys), dimension(im, kte), intent(in), target :: delp
     real(kind_phys), dimension(im, kte), intent(in), target :: phl3d
     real(kind_phys), dimension(im, kte), intent(in), target :: tk3d
     real(kind_phys), dimension(im, kte), intent(in), target :: us3d
     real(kind_phys), dimension(im, kte), intent(in), target :: vs3d
     real(kind_phys), dimension(im, kte), intent(in), target :: q3d
     real(kind_phys), dimension(im, kte), intent(in), target :: airden
     real(kind_phys), dimension(im, kte), intent(in), target :: rh
     real(kind_phys), dimension(im, kte), intent(in), target :: pfl_lsan
     real(kind_phys), dimension(im, kte), intent(in), target :: pfl_isan

     real(kind_phys), dimension(im), intent(in), target        :: rain_cpl
     real(kind_phys), dimension(im), intent(in), target        :: cldf
     real(kind_phys), dimension(im), intent(in), target        :: clayfrac
     real(kind_phys), dimension(im), intent(in), target        :: rdrag
     real(kind_phys), dimension(im), intent(in), target        :: sandfrac
     real(kind_phys), dimension(im), intent(in), target        :: ssm
     real(kind_phys), dimension(im), intent(in), target        :: ustar_threshold

     type(ccpp_constituent_prop_ptr_t), intent(in)    :: constituent_props(:)
     real(kind_phys), target,           intent(inout) :: constituents(:,:,:)

     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     ! Local variables
     real(kind_phys), target, allocatable :: subset_constituents(:,:,:)
     integer :: n_spec, i
     type(c_ptr) :: state_mgr

     errmsg = ''
     errflg = 0

     if (.not. do_catchem) return

     ! 1. Extract non-contiguous indices into a contiguous local subset array
     n_spec = size(catchem_indices_constituent_props)
     allocate(subset_constituents(im, kte, n_spec))

     do i = 1, n_spec
        subset_constituents(:,:,i) = constituents(:,:,catchem_indices_constituent_props(i))
     end do

     ! 2. Meterological pointer mapping (Direct Flat C-API bindings avoiding shape restrictions)
     state_mgr = cc_model%get_state_manager()

     ! 3D volumetric met fields
     call catchem_state_bind_met_3d(state_mgr, "T"//c_null_char, c_loc(tk3d(1,1)))
     call catchem_state_bind_met_3d(state_mgr, "QV"//c_null_char, c_loc(q3d(1,1)))
     call catchem_state_bind_met_3d(state_mgr, "RH"//c_null_char, c_loc(rh(1,1)))
     call catchem_state_bind_met_3d(state_mgr, "PMID"//c_null_char, c_loc(prl3d(1,1)))
     call catchem_state_bind_met_3d(state_mgr, "PEDGE"//c_null_char, c_loc(pr3d(1,1)))
     call catchem_state_bind_met_3d(state_mgr, "DELP"//c_null_char, c_loc(delp(1,1)))
     call catchem_state_bind_met_3d(state_mgr, "AIRDEN"//c_null_char, c_loc(airden(1,1)))
     call catchem_state_bind_met_3d(state_mgr, "AIRDEN_DRY"//c_null_char, c_loc(airden(1,1)))
     call catchem_state_bind_met_3d(state_mgr, "PFILSAN"//c_null_char, c_loc(pfl_isan(1,1)))
     call catchem_state_bind_met_3d(state_mgr, "PFLLSAN"//c_null_char, c_loc(pfl_lsan(1,1)))
     call catchem_state_bind_met_3d(state_mgr, "soil_moisture"//c_null_char, c_loc(soilmoist(1,1)))

     ! 2D surface met fields
     call catchem_state_bind_met_2d(state_mgr, "PS"//c_null_char, c_loc(prsfc(1)))
     call catchem_state_bind_met_2d(state_mgr, "TS"//c_null_char, c_loc(ts(1)))
     call catchem_state_bind_met_2d(state_mgr, "LAT"//c_null_char, c_loc(lat(1)))
     call catchem_state_bind_met_2d(state_mgr, "LON"//c_null_char, c_loc(lon(1)))
     call catchem_state_bind_met_2d(state_mgr, "PBLH"//c_null_char, c_loc(pblh(1)))
     call catchem_state_bind_met_2d(state_mgr, "USTAR"//c_null_char, c_loc(ustar(1)))
     call catchem_state_bind_met_2d(state_mgr, "HFLUX"//c_null_char, c_loc(hf2d(1)))
     call catchem_state_bind_met_2d(state_mgr, "AREA_M2"//c_null_char, c_loc(garea(1)))

     ! Auxiliary physics surface variables
     call catchem_state_bind_met_2d(state_mgr, "u_10m"//c_null_char, c_loc(u10m(1)))
     call catchem_state_bind_met_2d(state_mgr, "v_10m"//c_null_char, c_loc(v10m(1)))
     call catchem_state_bind_met_2d(state_mgr, "skin_temperature"//c_null_char, c_loc(tskin(1)))
     call catchem_state_bind_met_2d(state_mgr, "roughness_length"//c_null_char, c_loc(znt(1)))
     call catchem_state_bind_met_2d(state_mgr, "vegetation_fraction"//c_null_char, c_loc(gvf(1)))
     call catchem_state_bind_met_2d(state_mgr, "leaf_area_index"//c_null_char, c_loc(lai(1)))
     call catchem_state_bind_met_2d(state_mgr, "lake_fraction"//c_null_char, c_loc(lakefrac(1)))
     call catchem_state_bind_met_2d(state_mgr, "snow_fraction"//c_null_char, c_loc(frsnow(1)))

     ! Dynamic climatological dust variables mapped individually
     call catchem_state_bind_met_2d(state_mgr, "clay_fraction"//c_null_char, c_loc(clayfrac(1)))
     call catchem_state_bind_met_2d(state_mgr, "drag_coefficient"//c_null_char, c_loc(rdrag(1)))
     call catchem_state_bind_met_2d(state_mgr, "sand_fraction"//c_null_char, c_loc(sandfrac(1)))
     call catchem_state_bind_met_2d(state_mgr, "surface_soil_moisture"//c_null_char, c_loc(ssm(1)))
     call catchem_state_bind_met_2d(state_mgr, "threshold_friction_velocity"//c_null_char, c_loc(ustar_threshold(1)))

     ! 3. Bind local unified chemistry concentrations subset
     call catchem_state_bind_unified_chemistry(state_mgr, c_loc(subset_constituents(1,1,1)))

     ! 4. Central C++ Core timestep solver execution
     call cc_model%run_timestep(real(dt, fp), errflg)
     if (errflg /= CC_SUCCESS) then
         errmsg = 'CATChem Run Error: Scheduled process execution failed inside C++ Core.'
         deallocate(subset_constituents)
         return
     end if

     ! 5. Dynamic synchronisation back to the global host tracer array
     do i = 1, n_spec
        constituents(:,:,catchem_indices_constituent_props(i)) = subset_constituents(:,:,i)
     end do

     deallocate(subset_constituents)

   end subroutine ccpp_catchem_interface_run

end module ccpp_catchem_interface
