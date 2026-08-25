!> \file test_settling_science.f90
!! \brief Unit tests for settling science schemes and physical solvers
!!
program test_settling_science
   use testing_mod, only: assert
   use precision_mod, only: fp
   use SettlingPhysics_Mod, only: settling_calc_vsettle, settling_swelling_gerber, settling_compute

   implicit none

   write(*,*) 'Testing Settling Science Schemes...'
   write(*,*) ''

   ! Test 1: Particle swelling calculation
   write(*,*) 'Test 1: Settling Gerber swelling model'
   block
      integer, parameter :: nlev = 1
      real(fp) :: rh(nlev), radius_dry, rhop_dry
      real(fp) :: radius_wet(nlev), rhop_wet(nlev)

      rh = 0.8_fp           ! 80% relative humidity
      radius_dry = 1.0e-6_fp ! 1 micron dry radius
      rhop_dry = 2200.0_fp   ! Sea salt density
      call settling_swelling_gerber(nlev, rh, radius_dry, rhop_dry, radius_wet, rhop_wet)
      call assert(radius_wet(1) > radius_dry, "Wet radius at RH=80% should be larger than dry radius")
   end block
   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! Test 2: Terminal settling velocity calculation
   write(*,*) 'Test 2: Terminal settling velocity calculation'
   block
      real(fp) :: vsettle, density, radius, temp, rhoa, grav
      radius = 2.0e-6_fp  ! 2 micron radius
      density = 2200.0_fp ! Sea salt particle density [kg/m3]
      rhoa = 1.225_fp     ! Air density
      temp = 288.15_fp    ! 288.15 K
      grav = 9.80665_fp   ! Gravity

      call settling_calc_vsettle(radius, density, rhoa, temp, grav, vsettle)
      call assert(vsettle > 0.0_fp, "Settling velocity must be positive")
      call assert(vsettle < 1.0_fp, "Settling velocity for 2-micron particle should be < 1 m/s")
   end block
   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: 1D Vertical Column Settling Solver
   write(*,*) 'Test 3: 1D Vertical Column Settling Solver'
   block
      integer, parameter :: nlev = 5
      real(fp) :: conc(nlev), temp(nlev), airden(nlev), rh(nlev)
      real(fp) :: z_edge(nlev + 1), delp(nlev)
      real(fp) :: radius_dry, rhop_dry
      real(fp) :: dt, grav
      integer :: k, rc

      radius_dry = 5.0e-6_fp
      rhop_dry = 2000.0_fp
      dt = 3600.0_fp ! 1 hour
      grav = 9.80665_fp
      conc = 10.0_fp ! Uniform initial concentration
      temp = 280.0_fp
      airden = 1.2_fp
      rh = 0.5_fp
      delp = 2000.0_fp

      do k = 1, nlev + 1
         z_edge(k) = real(k - 1, fp) * 200.0_fp
      end do

      call settling_compute(nlev, 1, dt, grav, &
         radius_dry, rhop_dry, 0, & ! 0 = no swelling
         conc, temp, airden, rh, z_edge, delp, &
         rc=rc)

      call assert(rc == 0, "Settling compute return code should be 0")
      call assert(conc(nlev) <= 10.0_fp, "Top level concentration should decrease or remain constant due to settling")
   end block
   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   write(*,*) 'All Settling Science tests passed successfully!'

end program test_settling_science
