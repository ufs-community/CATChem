!> \file GasChemScheme_MUSICA_Mod.F90
!! \brief MICM gas phase chemistry solver without photolysis
!!
!! Pure science kernel for no_phot scheme in GasChem process.
!! This module contains ONLY the computational algorithm with NO infrastructure dependencies.
!! Uses only basic Fortran types for maximum portability and reusability.
!!
!! SCIENCE CUSTOMIZATION GUIDE:
!! 1. Modify the algorithm in compute_no_phot (search for "TODO")
!! 2. Add scheme-specific helper subroutines as needed
!! 3. Update physical constants for your scheme
!! 4. Customize the environmental response functions
!!
!! INFRASTRUCTURE RESPONSIBILITIES (handled by host model):
!! - Parameter initialization and validation
!! - Input array validation and error handling
!! - Memory management and array allocation
!! - Integration with host model time stepping
!!
!! Generated on: 2026-08-24T17:51:00.650243
!! Author: Maggie Bruckner
!! Reference: MUSICA library
module GasChemScheme_MUSICA_Mod

   use precision_mod, only: fp
   use GasChemCommon_Mod, only: GasChemSchemeMUSICAConfig
   use musica_micm, only: get_micm_version, Rosenbrock, RosenbrockStandardOrder
   use musica_micm, only: micm_t, solver_stats_t
   use musica_state, only: state_t, conditions_t
   use musica_util, only: assert, error_t, string_t, mapping_t

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_no_phot

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter :: T_STANDARD = 303.15_fp    ! Standard reference temperature [K]
   real(fp), parameter :: DEFAULT_SCALING = 1.0e-9_fp ! Default emission scaling factor

contains

   !> Pure science computation for no_phot scheme
   !!
   !! This is a pure computational kernel implementing MICM gas phase chemistry solver without photolysis.
   !! NO error checking, validation, or infrastructure concerns.
   !! Host model must ensure all inputs are valid before calling.
   !!
   !! @param[in]  num_layers     Number of vertical layers
   !! @param[in]  num_species    Number of chemical species
   !! @param[in]  params         Scheme parameters (pre-validated by host)
   !! @param[in]  RSTARG    Required constant from Constants module
   !! @param[in]  AIRMW    Required constant from Constants module
   !! @param[in]  airden    AIRDEN field [appropriate units]
   !! @param[in]  pmid    PMID field [appropriate units]
   !! @param[in]  sphu    SPHU field [appropriate units]
   !! @param[in]  t    T field [appropriate units]
   !! @param[in]  tstep    Time step [s] - retrieved from process interface
   !! @param[in]  species_short_name    Species short_name property
   !! @param[in]  species_conc   Species concentrations [ppm or ug/kg] (num_layers, num_species)
   !! @param[inout] species_tendencies  Species tendency terms [mol/mol/s] (num_layers, num_species)
   !! @param[inout] total_rate_per_species_per_level    Net chemical change per species per level [ppmv/s] (num_layers, num_species)
   !! @param[inout] net_chemical_rate_per_species    Net chem rate [ppmv/s] (num_species)
   !! @param[in] diagnostic_species_id Indices mapping diagnostic species to species array (optional, for per-species diagnostics)
   impure subroutine compute_no_phot( &                                                                                           
      num_layers, &
      num_species, &
      params, &
      RSTARG, &
      AIRMW, &
      airden, &
      pmid, &
      sphu, &
      t, &
      tstep, &
      species_short_name, &
      species_conc, &
      species_tendencies, &
      rc, &
      micm_in, state_in, &
      total_rate_per_species_per_level, &
      net_chemical_rate_per_species, &
      diagnostic_species_id &
   )

      ! Arguments
      integer, intent(in) :: num_layers
      integer, intent(in) :: num_species
      type(GasChemSchemeMUSICAConfig), intent(in) :: params
      real(fp), intent(in) :: RSTARG  ! Required constant from Constants module
      real(fp), intent(in) :: AIRMW  ! Required constant from Constants module
      real(fp), intent(in) :: airden(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: pmid(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: sphu(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: t(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: tstep  ! Time step [s] - from process interface
      character(len=32), intent(in) :: species_short_name(:)  ! Species short_name property
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(inout) :: species_tendencies(num_layers, num_species)
      integer, intent(out) :: rc  ! Return code (0 for success, non-zero for error)
      type(micm_t), pointer, optional :: micm_in
      type(state_t), pointer, optional :: state_in
      real(fp), intent(inout), optional :: total_rate_per_species_per_level(:,:)
      real(fp), intent(inout), optional :: net_chemical_rate_per_species(:)
      integer, intent(in), optional :: diagnostic_species_id(:)  ! Indices mapping diagnostic species to species array

      ! Local variables
      integer :: k, species_idx, micm_sp_idx, rp_idx, nrp, idx
      integer :: diag_idx  ! For diagnostic species indexing
      integer :: target_rp_internal_idx
      integer :: species_idx_map(num_species) ! Indices mapping catchem species to micm species
      real(fp) :: rate_val
      real(fp) :: conc
      real(fp) :: elapsed, remaining, advanced
      type(string_t) :: solver_state
      type(solver_stats_t) :: solver_stats
      type(error_t) :: micm_error
      type(micm_t), pointer :: micm
      type(state_t), pointer :: state
      logical :: micm_local_created = .false.
      logical :: state_local_created = .false.
      integer :: solver_type
      character(len=:), allocatable :: rp_name


      rc = 0
      elapsed = 0.0_fp
      ! Note: species_tendencies and diagnostic arrays are already initialized                                                    
      ! by the host ProcessInterface before calling this subroutine.                                                              
      ! Do not re-initialize them here.

      if (present(micm_in)) then
         micm => micm_in
      else
         solver_type = RosenbrockStandardOrder
         micm => micm_t(params%mechanism, solver_type, micm_error)

         if (.not. micm_error%is_success()) then
            write(*,'(A)') "Error creating MICM: ", micm_error%message()
            rc = 1
            return
         end if
         micm_local_created = .true.
      end if

      if (present(state_in)) then
         state => state_in
      else
         state => micm%get_state(num_layers,micm_error)
         if (.not. micm_error%is_success()) then
            write(*,'(A)') "Error creating MICM state: ", micm_error%message()
            rc = 1
            return
         end if
         state_local_created = .true.
      end if
      nrp = state%rate_parameters_ordering%size()

      ! initialize MICM rate parameters
      ! Note: LOSS parameters stay at 0 unless explicitly configured.
      ! default rate to 1.0 so the
      ! YAML scaling_factor defines the effective rate constant.
      do rp_idx = 1, nrp
         rp_name = trim(state%rate_parameters_ordering%name(rp_idx))
         if (rp_name(1:min(5,len(rp_name))) == 'LOSS.') then
            rate_val = 1.0_fp
         else
            rate_val = 0.0_fp
         endif
         target_rp_internal_idx = state%rate_parameters_ordering%index(state%rate_parameters_ordering%name(rp_idx), micm_error)

         do k = 1, num_layers
            ! 1-Based Struct Match: 1 + Grid Layer Stride + (Zero-indexed Variable Position * Variable Stride)
            idx = 1 + (k - 1) * state%rate_parameters_strides%grid_cell + &
               (target_rp_internal_idx - 1) * state%rate_parameters_strides%variable
            state%rate_parameters(idx) = rate_val
         end do
      end do

      call map_catchem_micm_species(num_species,              &
                                    state,            &
                                    species_short_name,            &
                                    species_idx_map)

      ! set up initial conditions for MICM
      do k = 1, num_layers
         state%conditions(k)%temperature = t(k)
         state%conditions(k)%pressure = pmid(k)
         state%conditions(k)%air_density = pmid(k) / (RSTARG * t(k))
         do species_idx = 1, num_species

            micm_sp_idx = 1 + (k-1)*state%species_strides%grid_cell + &
                    (species_idx_map(species_idx)-1)*state%species_strides%variable

            if (trim(species_short_name(species_idx)) == 'N2' .or. &
                trim(species_short_name(species_idx)) == 'n2') then
                ! Add N2 specific processing here
                conc = 0.7808_fp * state%conditions(k)%air_density

            else if (trim(species_short_name(species_idx)) == 'O2' .or. &
                     trim(species_short_name(species_idx)) == 'o2') then
               ! Add O2 specific processing here
               conc = 0.2095_fp * state%conditions(k)%air_density

            else if (trim(species_short_name(species_idx)) == 'H2O' .or. &
                     trim(species_short_name(species_idx)) == 'h2o') then

               ! Add H2O specific processing here
               conc = (sphu(k) * state%conditions(k)%air_density) / 0.01801528_fp 
            else
               ! Default processing for all other chemical species
               ! convert species concentrations from ppmv to mol/m3
               conc = species_conc(k, species_idx) * 1.0e-6_fp * state%conditions(k)%air_density
            endif
            state%concentrations(micm_sp_idx) = conc
         enddo
      enddo

      do while (elapsed < tstep)
         remaining = tstep - elapsed
         call micm%solve(REAL(remaining, 8),state,solver_state,solver_stats,micm_error)
         advanced = solver_stats%final_time()
         elapsed = elapsed + advanced
      enddo
      if (.not. micm_error%is_success()) then
         write(*,'(A)') "Error solving MICM: ", micm_error%message()
         rc = 1
         return
      end if

      do k = 1, num_layers
         do species_idx = 1, num_species
            micm_sp_idx = 1 + (k-1)*state%species_strides%grid_cell + &
                    (species_idx_map(species_idx)-1)*state%species_strides%variable
            ! convert final concentration back to ppmv from mol/m3
            conc = state%concentrations(micm_sp_idx)  / state%conditions(k)%air_density * 1e6_fp

            species_tendencies(k, species_idx) = conc

            ! Per-species-per-level diagnostic: 2D array (levels, species)
            if (present(total_rate_per_species_per_level) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == species_idx) then
                     total_rate_per_species_per_level(k, diag_idx) = species_conc(k,species_idx) - species_tendencies(k, species_idx)
                     exit
                  end if
               end do
            end if

            ! Per-species diagnostic: only update for diagnostic species
            if (present(net_chemical_rate_per_species) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == species_idx) then
                     net_chemical_rate_per_species(diag_idx) = species_tendencies(k, species_idx) * 1.0_fp  ! TODO: Replace with actual calculation
                     exit
                  end if
               end do
            end if
         end do

      end do
      
      if (state_local_created .and. associated(state)) then
         deallocate(state)
         state => null()
      end if

      if (micm_local_created .and. associated(micm)) then
         deallocate(micm)
         micm => null()
      end if
      if (allocated(rp_name)) deallocate(rp_name)
   end subroutine compute_no_phot

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   subroutine map_catchem_micm_species(nspecies,              &
      micm_state,            &
      catchem_species_names,            &
      mapped_indicies)

      IMPLICIT NONE

      INTEGER, intent(in) :: nspecies ! number of species
      type(state_t), pointer :: micm_state
      character(len=32), intent(in) :: catchem_species_names(:)  ! Species short_name property

      integer, intent(out) :: mapped_indicies(nspecies)

      integer :: species_idx, micm_idx
      type(error_t) :: micm_error

      do species_idx = 1, nspecies
         micm_idx = micm_state%species_ordering%index(trim(catchem_species_names(species_idx)), micm_error)

         ! CRITICAL CHECK: Did MICM find the species name?
         if (micm_idx < 0) then
            print *, "FATAL ERROR: CATChem species '", trim(catchem_species_names(species_idx)), &
                     "' could not be found in the MICM solver configuration!"
            ! Optional: print micm_error message if your binding supports it
            stop "CATChem-MICM Mapping Failed: Aborting execution during process_init."
         end if

         mapped_indicies(species_idx) = micm_idx
      end do

   end subroutine map_catchem_micm_species

end module GasChemScheme_MUSICA_Mod
