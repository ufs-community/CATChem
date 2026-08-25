!>
!! \file SO4chemPhysics_Mod.F90
!! \brief Internalized SO4 oxidation chemistry physics routines from GOCART2G.
!!
!! All routines operate on 1D column arrays in CATChem's native
!! bottom-to-top vertical ordering. No vertical reversal is needed.
!!
!! Ported from src/external/GOCART/Process_Library/GOCART2G_Process.F90
!! Routines: SulfateUpdateOxidants, SulfateChemDriver_DMS,
!!           SulfateChemDriver_SO2, SulfateChemDriver_SO4,
!!           SulfateChemDriver_MSA, idaynum
!!
!! \ingroup so4chem_process
!<
module SO4chemPhysics_Mod

   use precision_mod, only: fp
   use Met_Utilities_Mod, only: solar_zenith_angle
   implicit none
   private

   ! Public chemistry routines
   public :: so4chem_driver
   public :: so4chem_update_oxidants
   public :: so4chem_dms_oxidation
   public :: so4chem_so2_oxidation
   public :: so4chem_so4_update
   public :: so4chem_msa_update

   ! --- DMS oxidation constants ---
   ! OH addition channel: k_add = (C_ADD_A * exp(C_ADD_EA/T) * [O2])
   !                             / (1 + C_ADD_B * exp(C_ADD_EB/T) * [O2])
   real(fp), parameter :: C_ADD_A  = 1.7e-42_fp   !< OH addition pre-exponential [cm6 molecule-2 s-1]
   real(fp), parameter :: C_ADD_EA = 7810.0_fp     !< OH addition activation energy / R [K]
   real(fp), parameter :: C_ADD_B  = 5.5e-31_fp    !< OH addition denominator pre-exponential [cm3 molecule-1]
   real(fp), parameter :: C_ADD_EB = 7460.0_fp     !< OH addition denominator activation energy / R [K]

   ! OH abstraction channel: k_abs = C_ABS * exp(C_ABS_EA/T)
   real(fp), parameter :: C_ABS    = 1.2e-11_fp    !< OH abstraction pre-exponential [cm3 molecule-1 s-1]
   real(fp), parameter :: C_ABS_EA = -260.0_fp     !< OH abstraction activation energy / R [K]

   ! NO3 channel: k_no3 = C_NO3 * exp(C_NO3_EA/T)
   real(fp), parameter :: C_NO3    = 1.9e-13_fp    !< NO3 pre-exponential [cm3 molecule-1 s-1]
   real(fp), parameter :: C_NO3_EA = 500.0_fp      !< NO3 activation energy / R [K]

   ! Branching ratio for MSA from OH addition channel
   real(fp), parameter :: B_MSA    = 0.25_fp       !< MSA yield from DMS+OH addition

   ! O2 volume fraction in dry air
   real(fp), parameter :: F_O2     = 0.21_fp       !< O2 volume mixing ratio in dry air

   ! --- SO2 oxidation constants (Troe formulation) ---
   ! k0 = C_TROE_K0 * (300/T)^C_TROE_N
   real(fp), parameter :: C_TROE_K0 = 3.0e-31_fp   !< Low-pressure limit [cm6 molecule-2 s-1]
   real(fp), parameter :: C_TROE_N  = 4.3_fp        !< Temperature exponent
   real(fp), parameter :: C_TROE_KI = 1.5e-12_fp    !< High-pressure limit [cm3 molecule-1 s-1]
   real(fp), parameter :: C_TROE_FC = 0.6_fp        !< Broadening factor

   ! Aqueous chemistry temperature threshold
   real(fp), parameter :: T_AQ_MIN = 258.0_fp       !< Minimum temperature for aqueous SO2 oxidation [K]

   ! Pressure lid threshold
   real(fp), parameter :: PLID_HPA = 0.01_fp        !< Pressure lid [hPa]

contains

   !---------------------------------------------------------------------------
   !> Compute day-of-year from YYYYMMDD integer.
   !! Ported from GOCART2G idaynum.
   !!
   !! @param[in] nymd  Date in YYYYMMDD format
   !! @return          Day of year (1–366)
   !---------------------------------------------------------------------------
   pure integer function idaynum(nymd)
      integer, intent(in) :: nymd

      integer :: yyyy, mm, dd, imon, isleapyr
      integer :: ndays(12)

      ndays = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

      yyyy = nymd / 10000
      mm   = mod(nymd, 10000) / 100
      dd   = mod(nymd, 100)

      ! Leap year determination
      isleapyr = 0
      if (mod(yyyy, 4) == 0) then
         isleapyr = 1
         if (mod(yyyy, 100) == 0) then
            isleapyr = 0
            if (mod(yyyy, 400) == 0) then
               isleapyr = 1
            end if
         end if
      end if

      ! Accumulate day number
      idaynum = 0
      if (mm == 1) then
         idaynum = dd
      else
         do imon = 1, mm - 1
            if (imon == 2 .and. isleapyr == 1) then
               idaynum = idaynum + 29
            else
               idaynum = idaynum + ndays(imon)
            end if
         end do
         idaynum = idaynum + dd
      end if

   end function idaynum

   !---------------------------------------------------------------------------
   !> Scale climatological OH, NO3, H2O2 fields based on solar zenith angle
   !! diurnal variation for a single column.
   !! Replaces GOCART2G SulfateUpdateOxidants.
   !!
   !! @param[in]    nlev          Number of vertical levels
   !! @param[in]    cdt           Chemistry timestep [s]
   !! @param[in]    nymd          Date in YYYYMMDD format
   !! @param[in]    nhms          Time in HHMMSS format
   !! @param[in]    lat_rad       Latitude [radians]
   !! @param[in]    lon_rad       Longitude [radians]
   !! @param[in]    airMolWght    Molecular weight of air [kg/kmol]
   !! @param[in]    nAvogadro     Avogadro's number [molecules/kmol]
   !! @param[in]    oh_clim       Climatological OH [VMR]
   !! @param[in]    no3_clim      Climatological NO3 [VMR]
   !! @param[in]    h2o2_clim     Climatological H2O2 [VMR]
   !! @param[inout] xoh           Scaled OH [molecules cm-3]
   !! @param[inout] xno3          Scaled NO3 [VMR]
   !! @param[inout] xh2o2         H2O2 [VMR]
   !! @param[inout] h2o2_init     Saved H2O2 initial field [VMR]
   !! @param[inout] recycle_h2o2  Flag to recycle H2O2 to climatology
   !! @param[in]    rhoa          Air density [kg/m3]
   !! @param[out]   rc            Return code (0=success)
   !---------------------------------------------------------------------------
   subroutine so4chem_update_oxidants(nlev, cdt, nymd, nhms, lat_rad, lon_rad, &
      airMolWght, nAvogadro, oh_clim, no3_clim, h2o2_clim, &
      xoh, xno3, xh2o2, h2o2_init, recycle_h2o2, rhoa, rc)

      integer,  intent(in)    :: nlev
      real(fp), intent(in)    :: cdt
      integer,  intent(in)    :: nymd
      integer,  intent(in)    :: nhms
      real(fp), intent(in)    :: lat_rad
      real(fp), intent(in)    :: lon_rad
      real(fp), intent(in)    :: airMolWght
      real(fp), intent(in)    :: nAvogadro
      real(fp), intent(in)    :: oh_clim(nlev)
      real(fp), intent(in)    :: no3_clim(nlev)
      real(fp), intent(in)    :: h2o2_clim(nlev)
      real(fp), intent(inout) :: xoh(nlev)
      real(fp), intent(inout) :: xno3(nlev)
      real(fp), intent(inout) :: xh2o2(nlev)
      real(fp), intent(inout) :: h2o2_init(nlev)
      logical,  intent(inout) :: recycle_h2o2
      real(fp), intent(in)    :: rhoa(nlev)
      integer,  intent(out), optional :: rc

      ! Local variables
      integer  :: k, n, jday, ndystep
      real(fp) :: xhour, xhouruse
      real(fp) :: sza_deg, cossza, cossza_now
      real(fp) :: tcosz, tday, tnight

      if (present(rc)) rc = 0

      ! Compute day of year and hour
      jday = idaynum(nymd)
      xhour = (real(nhms / 10000, fp) * 3600.0_fp &
         + real(mod(nhms, 10000) / 100, fp) * 60.0_fp &
         + real(mod(nhms, 100), fp)) / 3600.0_fp

      ! Recycle H2O2 to climatology if flag is set
      if (recycle_h2o2) then
         xh2o2 = h2o2_clim
         recycle_h2o2 = .false.
      end if

      ! Initialize OH from climatology
      xoh = oh_clim

      ! Integrate cos(SZA) over the full day to compute normalization
      ndystep = nint(86400.0_fp / cdt)
      tcosz = 0.0_fp
      tday  = 0.0_fp
      xhouruse = xhour

      do n = 1, ndystep
         call solar_zenith_angle(jday, xhouruse, lat_rad, lon_rad, sza_deg, cossza)
         tcosz = tcosz + cossza
         xhouruse = xhouruse + cdt / 3600.0_fp
         if (xhouruse > 24.0_fp) xhouruse = xhouruse - 24.0_fp
         if (cossza > 0.0_fp) tday = tday + cdt
      end do

      ! Get cos(SZA) at current time
      call solar_zenith_angle(jday, xhour, lat_rad, lon_rad, sza_deg, cossza_now)

      tnight = 86400.0_fp - tday

      ! Scale OH: proportional to cossza_now / tcosz during daytime
      do k = 1, nlev
         if (tcosz > 0.0_fp) then
            xoh(k) = oh_clim(k) * (86400.0_fp / cdt) * cossza_now / tcosz
         else
            xoh(k) = 0.0_fp
         end if
      end do

      ! Clamp negative OH to zero
      where (xoh < 0.0_fp) xoh = 0.0_fp

      ! Convert OH from VMR to number density [molecules cm-3]
      do k = 1, nlev
         xoh(k) = xoh(k) * 1000.0_fp * rhoa(k) / airMolWght * nAvogadro * 1.0e-6_fp
      end do

      ! Scale NO3: zero during daytime, scaled by 86400/tnight at night
      do k = 1, nlev
         if (cossza_now > 0.0_fp .or. tnight < tiny(1.0_fp)) then
            xno3(k) = 0.0_fp
         else
            xno3(k) = no3_clim(k) * 86400.0_fp / tnight
         end if
      end do

   end subroutine so4chem_update_oxidants

   !---------------------------------------------------------------------------
   !> Compute DMS oxidation by OH (addition + abstraction) and NO3.
   !! Replaces GOCART2G SulfateChemDriver_DMS.
   !!
   !! @param[in]    nlev        Number of vertical levels
   !! @param[in]    klid        Index for pressure lid (bottom of active range)
   !! @param[in]    cdt         Chemistry timestep [s]
   !! @param[in]    airMolWght  Molecular weight of air [kg/kmol]
   !! @param[in]    nAvogadro   Avogadro's number [molecules/kmol]
   !! @param[in]    fMassMSA    Molecular weight of MSA [g/mol]
   !! @param[in]    fMassDMS    Molecular weight of DMS [g/mol]
   !! @param[in]    fMassSO2    Molecular weight of SO2 [g/mol]
   !! @param[inout] dms         DMS concentration [kg/kg]
   !! @param[in]    xoh         OH concentration [molecules cm-3]
   !! @param[in]    xno3        NO3 concentration [VMR]
   !! @param[in]    cossza      Cosine of solar zenith angle (scalar)
   !! @param[in]    tmpu        Temperature [K]
   !! @param[in]    rhoa        Air density [kg/m3]
   !! @param[out]   pso2_dms    SO2 production from DMS oxidation [kg/kg/s]
   !! @param[out]   pmsa_dms    MSA production from DMS oxidation [kg/kg/s]
   !! @param[out]   rc          Return code (0=success)
   !---------------------------------------------------------------------------
   subroutine so4chem_dms_oxidation(nlev, klid, cdt, airMolWght, nAvogadro, &
      fMassMSA, fMassDMS, fMassSO2, dms, xoh, xno3, cossza, tmpu, rhoa, &
      pso2_dms, pmsa_dms, rc)

      integer,  intent(in)    :: nlev
      integer,  intent(in)    :: klid
      real(fp), intent(in)    :: cdt
      real(fp), intent(in)    :: airMolWght
      real(fp), intent(in)    :: nAvogadro
      real(fp), intent(in)    :: fMassMSA
      real(fp), intent(in)    :: fMassDMS
      real(fp), intent(in)    :: fMassSO2
      real(fp), intent(inout) :: dms(nlev)
      real(fp), intent(in)    :: xoh(nlev)
      real(fp), intent(in)    :: xno3(nlev)
      real(fp), intent(in)    :: cossza
      real(fp), intent(in)    :: tmpu(nlev)
      real(fp), intent(in)    :: rhoa(nlev)
      real(fp), intent(out)   :: pso2_dms(nlev)
      real(fp), intent(out)   :: pmsa_dms(nlev)
      integer,  intent(out), optional :: rc

      ! Local variables
      integer  :: k
      real(fp) :: tk, oh, air, o2, no3
      real(fp) :: rk1, rk2, rk3
      real(fp) :: dms0, dms_oh, dms_final

      if (present(rc)) rc = 0

      pso2_dms = 0.0_fp
      pmsa_dms = 0.0_fp

      do k = klid, nlev

         rk1 = 0.0_fp
         rk2 = 0.0_fp
         rk3 = 0.0_fp

         tk  = tmpu(k)
         oh  = xoh(k)

         ! Air molecules in # cm-3
         air = 1000.0_fp * rhoa(k) / airMolWght * nAvogadro * 1.0e-6_fp
         ! Oxygen molecules in # cm-3
         o2  = F_O2 * air
         ! NO3: go from volume mixing ratio to # cm-3
         no3 = xno3(k) * air

         ! Initial DMS concentration (kg/kg)
         dms0 = dms(k)
         dms0 = max(dms0, tiny(dms0))

         ! 1 & 2) DMS + OH: rk1 = addition, rk2 = abstraction
         if (oh > 0.0_fp) then
            rk1 = (C_ADD_A * exp(C_ADD_EA / tk) * o2) / &
               (1.0_fp + C_ADD_B * exp(C_ADD_EB / tk) * o2) * oh
            rk2 = (C_ABS * exp(C_ABS_EA / tk)) * oh
         end if

         ! 3) DMS + NO3: only happens at night
         if (cossza <= 0.0_fp) then
            rk3 = (C_NO3 * exp(C_NO3_EA / tk)) * no3
         end if

         ! DMS loss via OH then NO3
         dms_oh    = dms0 * exp(-(rk1 + rk2) * cdt)
         dms_final = dms_oh * exp(-rk3 * cdt)

         ! MSA production from OH addition channel
         if (abs(rk1 + rk2) <= 1.0e-30_fp) then
            pmsa_dms(k) = 0.0_fp
         else
            pmsa_dms(k) = (dms0 - dms_oh) * B_MSA * rk1 / (rk1 + rk2) &
               * (fMassMSA / fMassDMS) / cdt
         end if

         ! SO2 production: everything else
         pso2_dms(k) = (dms0 - dms_final &
            - pmsa_dms(k) * cdt * (fMassDMS / fMassMSA)) &
            * (fMassSO2 / fMassDMS) / cdt

         ! Update DMS
         dms_final = max(dms_final, tiny(dms_final))
         dms(k) = dms_final

      end do

   end subroutine so4chem_dms_oxidation

   !---------------------------------------------------------------------------
   !> Compute SO2 gas-phase and aqueous-phase oxidation to SO4.
   !! Replaces GOCART2G SulfateChemDriver_SO2.
   !! No dry deposition in this version (rk2=0 always).
   !!
   !! @param[in]    nlev        Number of vertical levels
   !! @param[in]    klid        Index for pressure lid
   !! @param[in]    cdt         Chemistry timestep [s]
   !! @param[in]    airMolWght  Molecular weight of air [kg/kmol]
   !! @param[in]    nAvogadro   Avogadro's number [molecules/kmol]
   !! @param[in]    grav        Gravitational acceleration [m/s2]
   !! @param[in]    fMassSO4    Molecular weight of SO4 [g/mol]
   !! @param[in]    fMassSO2    Molecular weight of SO2 [g/mol]
   !! @param[inout] so2         SO2 concentration [kg/kg]
   !! @param[in]    xoh         OH concentration [molecules cm-3]
   !! @param[inout] xh2o2       H2O2 concentration [VMR]
   !! @param[in]    tmpu        Temperature [K]
   !! @param[in]    rhoa        Air density [kg/m3]
   !! @param[in]    delp        Pressure thickness [Pa]
   !! @param[in]    cloud       Cloud fraction [0-1]
   !! @param[in]    pso2_dms    SO2 production from DMS oxidation [kg/kg/s]
   !! @param[out]   pso4g_so2   SO4 gas-phase production [kg/kg/s]
   !! @param[out]   pso4aq_so2  SO4 aqueous-phase production [kg/kg/s]
   !! @param[out]   rc          Return code (0=success)
   !---------------------------------------------------------------------------
   subroutine so4chem_so2_oxidation(nlev, klid, cdt, airMolWght, nAvogadro, grav, &
      fMassSO4, fMassSO2, so2, xoh, xh2o2, tmpu, rhoa, delp, cloud, &
      pso2_dms, pso4g_so2, pso4aq_so2, rc)

      integer,  intent(in)    :: nlev
      integer,  intent(in)    :: klid
      real(fp), intent(in)    :: cdt
      real(fp), intent(in)    :: airMolWght
      real(fp), intent(in)    :: nAvogadro
      real(fp), intent(in)    :: grav
      real(fp), intent(in)    :: fMassSO4
      real(fp), intent(in)    :: fMassSO2
      real(fp), intent(inout) :: so2(nlev)
      real(fp), intent(in)    :: xoh(nlev)
      real(fp), intent(inout) :: xh2o2(nlev)
      real(fp), intent(in)    :: tmpu(nlev)
      real(fp), intent(in)    :: rhoa(nlev)
      real(fp), intent(in)    :: delp(nlev)
      real(fp), intent(in)    :: cloud(nlev)
      real(fp), intent(in)    :: pso2_dms(nlev)
      real(fp), intent(out)   :: pso4g_so2(nlev)
      real(fp), intent(out)   :: pso4aq_so2(nlev)
      integer,  intent(out), optional :: rc

      ! Local variables
      integer  :: k
      real(fp) :: tk, oh, h2o2, air
      real(fp) :: k0, kk, f1, rk1, rk, rkt
      real(fp) :: so20, so2_cd, so2_final
      real(fp) :: L1, L2, fc, fMR

      if (present(rc)) rc = 0

      pso4g_so2  = 0.0_fp
      pso4aq_so2 = 0.0_fp

      ! Conversion factor: SO2 mmr to SO2 vmr
      fMR = airMolWght / fMassSO2

      do k = klid, nlev

         rk1 = 0.0_fp
         L1  = 0.0_fp
         L2  = 0.0_fp

         tk   = tmpu(k)
         oh   = xoh(k)
         h2o2 = max(xh2o2(k), tiny(xh2o2(k)))

         ! Air molecules in # cm-3
         air = 1000.0_fp * rhoa(k) / airMolWght * nAvogadro * 1.0e-6_fp

         ! 1) SO2 + OH(g): Troe formulation
         k0  = C_TROE_K0 * (300.0_fp / tk)**C_TROE_N
         kk  = k0 * air / C_TROE_KI
         f1  = (1.0_fp + (log10(kk))**2)**(-1.0_fp)
         rk1 = (k0 * air / (1.0_fp + kk)) * C_TROE_FC**f1 * oh

         ! No dry deposition: rk2 = 0
         rk  = rk1
         rkt = rk * cdt

         ! Initial SO2 concentration after adding DMS source
         so20 = so2(k) + pso2_dms(k) * cdt
         so20 = max(so20, tiny(so20))

         ! Gas-phase SO2 loss
         if (rk > 0.0_fp) then
            so2_cd = so20 * exp(-rkt)
            L1     = (so20 - so2_cd) * rk1 / rk
         else
            so2_cd = so20
            L1     = 0.0_fp
         end if

         ! Aqueous-phase cloud chemistry
         fc = cloud(k)
         if (fc > 0.0_fp .and. so2_cd > 0.0_fp .and. tk > T_AQ_MIN) then
            ! Check if H2O2 vmr limits SO2 oxidation
            if (fMR * so2_cd > h2o2) then
               fc   = fc * (h2o2 / (fMR * so2_cd))
               h2o2 = h2o2 * (1.0_fp - cloud(k))
            else
               h2o2 = h2o2 * (1.0_fp - cloud(k) * (fMR * so2_cd) / h2o2)
            end if
            so2_final = so2_cd * (1.0_fp - fc)
            L2 = so2_cd * fc
         else
            so2_final = so2_cd
            L2 = 0.0_fp
         end if

         ! Update H2O2
         xh2o2(k) = max(h2o2, tiny(h2o2))

         ! Update SO2
         so2_final = max(so2_final, tiny(so2_final))
         so2(k) = so2_final

         ! SO4 production rates
         pso4g_so2(k)  = L1 * (fMassSO4 / fMassSO2) / cdt
         pso4aq_so2(k) = L2 * (fMassSO4 / fMassSO2) / cdt

      end do

   end subroutine so4chem_so2_oxidation

   !---------------------------------------------------------------------------
   !> Update SO4 concentration from gas-phase and aqueous-phase production.
   !! Replaces GOCART2G SulfateChemDriver_SO4.
   !! No dry deposition in this version.
   !!
   !! @param[in]    nlev        Number of vertical levels
   !! @param[in]    klid        Index for pressure lid
   !! @param[in]    cdt         Chemistry timestep [s]
   !! @param[in]    grav        Gravitational acceleration [m/s2]
   !! @param[inout] so4         SO4 concentration [kg/kg]
   !! @param[in]    delp        Pressure thickness [Pa]
   !! @param[in]    pso4g_so2   SO4 gas-phase production [kg/kg/s]
   !! @param[in]    pso4aq_so2  SO4 aqueous-phase production [kg/kg/s]
   !! @param[out]   rc          Return code (0=success)
   !---------------------------------------------------------------------------
   subroutine so4chem_so4_update(nlev, klid, cdt, grav, so4, delp, &
      pso4g_so2, pso4aq_so2, rc)

      integer,  intent(in)    :: nlev
      integer,  intent(in)    :: klid
      real(fp), intent(in)    :: cdt
      real(fp), intent(in)    :: grav
      real(fp), intent(inout) :: so4(nlev)
      real(fp), intent(in)    :: delp(nlev)
      real(fp), intent(in)    :: pso4g_so2(nlev)
      real(fp), intent(in)    :: pso4aq_so2(nlev)
      integer,  intent(out), optional :: rc

      ! Local variables
      integer  :: k
      real(fp) :: so40, pso4

      if (present(rc)) rc = 0

      do k = klid, nlev
         so40 = so4(k)
         so40 = max(so40, tiny(so40))

         pso4 = pso4g_so2(k) + pso4aq_so2(k)

         ! Simple additive tendency (no dry deposition)
         so4(k) = so40 + pso4 * cdt
         so4(k) = max(so4(k), tiny(so4(k)))
      end do

   end subroutine so4chem_so4_update

   !---------------------------------------------------------------------------
   !> Update MSA concentration from DMS oxidation production.
   !! Replaces GOCART2G SulfateChemDriver_MSA.
   !! No dry deposition in this version.
   !!
   !! @param[in]    nlev        Number of vertical levels
   !! @param[in]    klid        Index for pressure lid
   !! @param[in]    cdt         Chemistry timestep [s]
   !! @param[in]    grav        Gravitational acceleration [m/s2]
   !! @param[inout] msa         MSA concentration [kg/kg]
   !! @param[in]    delp        Pressure thickness [Pa]
   !! @param[in]    pmsa_dms    MSA production from DMS oxidation [kg/kg/s]
   !! @param[out]   rc          Return code (0=success)
   !---------------------------------------------------------------------------
   subroutine so4chem_msa_update(nlev, klid, cdt, grav, msa, delp, &
      pmsa_dms, rc)

      integer,  intent(in)    :: nlev
      integer,  intent(in)    :: klid
      real(fp), intent(in)    :: cdt
      real(fp), intent(in)    :: grav
      real(fp), intent(inout) :: msa(nlev)
      real(fp), intent(in)    :: delp(nlev)
      real(fp), intent(in)    :: pmsa_dms(nlev)
      integer,  intent(out), optional :: rc

      ! Local variables
      integer  :: k
      real(fp) :: msa0

      if (present(rc)) rc = 0

      do k = klid, nlev
         msa0 = msa(k)
         msa0 = max(msa0, tiny(msa0))

         ! Simple additive tendency (no dry deposition)
         msa(k) = msa0 + pmsa_dms(k) * cdt
         msa(k) = max(msa(k), tiny(msa(k)))
      end do

   end subroutine so4chem_msa_update

   !---------------------------------------------------------------------------
   !> Top-level SO4 chemistry driver for a single column.
   !! Calls sub-routines in sequence: update_oxidants -> dms_oxidation ->
   !! so2_oxidation -> so4_update -> msa_update.
   !!
   !! @param[in]    nlev          Number of vertical levels
   !! @param[in]    klid          Index for pressure lid
   !! @param[in]    cdt           Chemistry timestep [s]
   !! @param[in]    nymd          Date in YYYYMMDD format
   !! @param[in]    nhms          Time in HHMMSS format
   !! @param[in]    lat_rad       Latitude [radians]
   !! @param[in]    lon_rad       Longitude [radians]
   !! @param[in]    airMolWght    Molecular weight of air [kg/kmol]
   !! @param[in]    nAvogadro     Avogadro's number [molecules/kmol]
   !! @param[in]    grav          Gravitational acceleration [m/s2]
   !! @param[in]    fMassDMS      Molecular weight of DMS [g/mol]
   !! @param[in]    fMassSO2      Molecular weight of SO2 [g/mol]
   !! @param[in]    fMassSO4      Molecular weight of SO4 [g/mol]
   !! @param[in]    fMassMSA      Molecular weight of MSA [g/mol]
   !! @param[inout] dms           DMS concentration [kg/kg]
   !! @param[inout] so2           SO2 concentration [kg/kg]
   !! @param[inout] so4           SO4 concentration [kg/kg]
   !! @param[inout] msa           MSA concentration [kg/kg]
   !! @param[in]    oh_clim       Climatological OH [VMR]
   !! @param[in]    no3_clim      Climatological NO3 [VMR]
   !! @param[in]    h2o2_clim     Climatological H2O2 [VMR]
   !! @param[inout] xoh           Working OH [molecules cm-3]
   !! @param[inout] xno3          Working NO3 [VMR]
   !! @param[inout] xh2o2         Working H2O2 [VMR]
   !! @param[inout] h2o2_init     Saved H2O2 initial field [VMR]
   !! @param[inout] recycle_h2o2  Flag to recycle H2O2 to climatology
   !! @param[in]    tmpu          Temperature [K]
   !! @param[in]    rhoa          Air density [kg/m3]
   !! @param[in]    delp          Pressure thickness [Pa]
   !! @param[in]    cloud         Cloud fraction [0-1]
   !! @param[in]    lwi           Land-water-ice flag (scalar)
   !! @param[out]   pso2_dms      SO2 production from DMS [kg/kg/s]
   !! @param[out]   pmsa_dms      MSA production from DMS [kg/kg/s]
   !! @param[out]   pso4g_so2     SO4 gas-phase production [kg/kg/s]
   !! @param[out]   pso4aq_so2    SO4 aqueous-phase production [kg/kg/s]
   !! @param[out]   rc            Return code (0=success)
   !---------------------------------------------------------------------------
   subroutine so4chem_driver(nlev, klid, cdt, nymd, nhms, lat_rad, lon_rad, &
      airMolWght, nAvogadro, grav, fMassDMS, fMassSO2, fMassSO4, fMassMSA, &
      dms, so2, so4, msa, oh_clim, no3_clim, h2o2_clim, &
      xoh, xno3, xh2o2, h2o2_init, recycle_h2o2, &
      tmpu, rhoa, delp, cloud, lwi, &
      pso2_dms, pmsa_dms, pso4g_so2, pso4aq_so2, rc)

      integer,  intent(in)    :: nlev
      integer,  intent(in)    :: klid
      real(fp), intent(in)    :: cdt
      integer,  intent(in)    :: nymd
      integer,  intent(in)    :: nhms
      real(fp), intent(in)    :: lat_rad
      real(fp), intent(in)    :: lon_rad
      real(fp), intent(in)    :: airMolWght
      real(fp), intent(in)    :: nAvogadro
      real(fp), intent(in)    :: grav
      real(fp), intent(in)    :: fMassDMS
      real(fp), intent(in)    :: fMassSO2
      real(fp), intent(in)    :: fMassSO4
      real(fp), intent(in)    :: fMassMSA
      real(fp), intent(inout) :: dms(nlev)
      real(fp), intent(inout) :: so2(nlev)
      real(fp), intent(inout) :: so4(nlev)
      real(fp), intent(inout) :: msa(nlev)
      real(fp), intent(in)    :: oh_clim(nlev)
      real(fp), intent(in)    :: no3_clim(nlev)
      real(fp), intent(in)    :: h2o2_clim(nlev)
      real(fp), intent(inout) :: xoh(nlev)
      real(fp), intent(inout) :: xno3(nlev)
      real(fp), intent(inout) :: xh2o2(nlev)
      real(fp), intent(inout) :: h2o2_init(nlev)
      logical,  intent(inout) :: recycle_h2o2
      real(fp), intent(in)    :: tmpu(nlev)
      real(fp), intent(in)    :: rhoa(nlev)
      real(fp), intent(in)    :: delp(nlev)
      real(fp), intent(in)    :: cloud(nlev)
      integer,  intent(in)    :: lwi
      real(fp), intent(out)   :: pso2_dms(nlev)
      real(fp), intent(out)   :: pmsa_dms(nlev)
      real(fp), intent(out)   :: pso4g_so2(nlev)
      real(fp), intent(out)   :: pso4aq_so2(nlev)
      integer,  intent(out), optional :: rc

      ! Local variables
      real(fp) :: sza_deg, cossza
      integer  :: jday

      if (present(rc)) rc = 0

      ! Initialize production rate arrays
      pso2_dms   = 0.0_fp
      pmsa_dms   = 0.0_fp
      pso4g_so2  = 0.0_fp
      pso4aq_so2 = 0.0_fp

      ! 1) Update oxidant fields (diurnal scaling)
      call so4chem_update_oxidants(nlev, cdt, nymd, nhms, lat_rad, lon_rad, &
         airMolWght, nAvogadro, oh_clim, no3_clim, h2o2_clim, &
         xoh, xno3, xh2o2, h2o2_init, recycle_h2o2, rhoa)

      ! Compute cos(SZA) for the column at current time
      jday = idaynum(nymd)
      call solar_zenith_angle(jday, &
         (real(nhms / 10000, fp) * 3600.0_fp &
         + real(mod(nhms, 10000) / 100, fp) * 60.0_fp &
         + real(mod(nhms, 100), fp)) / 3600.0_fp, &
         lat_rad, lon_rad, sza_deg, cossza)

      ! 2) DMS oxidation
      call so4chem_dms_oxidation(nlev, klid, cdt, airMolWght, nAvogadro, &
         fMassMSA, fMassDMS, fMassSO2, dms, xoh, xno3, cossza, tmpu, rhoa, &
         pso2_dms, pmsa_dms)

      ! 3) SO2 oxidation
      call so4chem_so2_oxidation(nlev, klid, cdt, airMolWght, nAvogadro, grav, &
         fMassSO4, fMassSO2, so2, xoh, xh2o2, tmpu, rhoa, delp, cloud, &
         pso2_dms, pso4g_so2, pso4aq_so2)

      ! 4) SO4 update
      call so4chem_so4_update(nlev, klid, cdt, grav, so4, delp, &
         pso4g_so2, pso4aq_so2)

      ! 5) MSA update
      call so4chem_msa_update(nlev, klid, cdt, grav, msa, delp, pmsa_dms)

      ! Save H2O2 after chemistry
      h2o2_init = xh2o2

   end subroutine so4chem_driver

end module SO4chemPhysics_Mod
