!>
!! \file species_mod.F90
!! \brief This file contains the module for catchem species
!!
!! \ingroup core_modules
!!
!! This file contains the module for catchem species
!!
!!!>

module species_mod

   use precision_mod
   implicit none

   !> \brief Module for catchem species
   !!
   !! This module contains subroutines and functions related to the catchem species
   !!
   !! \ingroup core_modules
   !!
   !! \param Config The input config object.
   !! \param Species The Species object to be initialized.
   !! \param RC The return code.
   !!
   !!!>
   type, public :: SpeciesType

      ! Names
      character(len=30) :: long_name  !< long name for species used for netcdf attribute "long_name"
      character(len=30) :: short_name !< short name for species
      character(len=50) :: description !< description of species

      ! Logcial switches
      logical :: is_gas               !< if true, species is a gas and not an aerosol
      logical :: is_aerosol           !< if true, species is aerosol and not a gas
      logical :: is_tracer            !< if true, species is a tracer and not an aerosol or gas that undergoes chemistry or photolysis
      logical :: is_advected          !< if true, species is advected
      logical :: is_drydep            !< if true, species undergoes dry depotiion
      logical :: is_wetdep            !< if true, species undergoes wet depotiion
      logical :: is_photolysis        !< if true, species undergoes photolysis
      logical :: is_gocart_aero       !< if true, species is a GOCART aerosol species
      logical :: is_dust              !< if true, species is a dust
      logical :: is_seasalt           !< if true, species is a seasalt

      ! Numerical properties
      real(kind=fp) :: mw_g                 !< gaseous molecular weight
      real(kind=fp) :: density              !< particle density (kg/m3)
      real(kind=fp) :: radius               !< mean molecular diameter in meters
      real(kind=fp) :: lower_radius         !< lower radius in meters
      real(kind=fp) :: upper_radius         !< upper radius in meters
      real(kind=fp) :: viscosity            !< kinematic viscosity (m2/s)

      ! used for dry deposition
      real(kind=fp) :: dd_f0                !< reactivity factor for oxidation of biological substances
      real(kind=fp) :: dd_hstar             !< Henry’s law constant
      real(kind=fp) :: dd_DvzAerSnow        !< fix dry deposition velocity (cm/s) over ice and snow for certain aerosol species
      real(kind=fp) :: dd_DvzMinVal_snow    !< minimum dry deposition velocity (cm/s) over snow and ice
      real(kind=fp) :: dd_DvzMinVal_land    !< minimum dry deposition velocity (cm/s) over land

      ! used for wet deposition
      real(kind=fp) :: radius_wet           !< mean molecular diameter in meters for wet conditions
      real(kind=fp) :: henry_k0             !< Henry’s law solubility constant ( M / atm)
      real(kind=fp) :: henry_cr             !< Henry’s law volatility constant (K)
      real(kind=fp) :: henry_pKa            !< Henry’s Law pH correction factor (seems zeros for all species now)
      real(kind=fp) :: wd_retfactor         !< retention efficiency of species in the liquid cloud condensate as it is converted to precipitation
      logical       :: wd_LiqAndGas         !< whether the ice-to-gas ratio can be computed for this species by co-condensation
      real(kind=fp) :: wd_convfacI2G        !< conversion factor for computing the ice-to-gas ratio by co-condensation when wd_LiqAndGas = .true.
      real(kind=fp) :: wd_rainouteff1       !< temperature-dependent scale factor that is used to multiply the fraction of species scavenged by rainout.
      real(kind=fp) :: wd_rainouteff2       !<  T < 237k for (1);  237 <= T < 258k for (2);  T >= 258k for (3)
      real(kind=fp) :: wd_rainouteff3       !< TODO: maybe we can read it in as a list directly, instead of separating to 3 ??


      ! Default background concentration
      real(kind=fp) :: BackgroundVV        !< Background conc [v/v]

      ! Indices
      integer :: species_index        !< species index in species array
      integer :: drydep_index         !< drydep index in drydep array
      integer :: wetdep_index         !< drydep index in wetdep array
      integer :: photolysis_index     !< photolysis index in photolysis array
      integer :: gocart_aero_index    !< gocart_aero index in gocart_aero array

      ! Concentration
      real(kind=fp), ALLOCATABLE :: conc(:)             !< species concentration [v/v] or kg/kg

   end type SpeciesType

   !
   ! !DEFINED PARAMETERS:
   !
   !=========================================================================
   ! Missing species concentration value if not in restart file and special
   ! background value not defined
   !=========================================================================
   REAL(fp), PARAMETER, PUBLIC :: MISSING_VV  = 1.0e-20_fp ! Missing spc conc

contains

   subroutine init(Species_State, species_name, atomic_num)
      type(SpeciesType), intent(inout) :: Species_State
      character(len=*), intent(in) :: species_name
      integer, intent(in) :: atomic_num

      Species_State%short_name = species_name
      Species_State%mw_g = atomic_num
   end subroutine init

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
