!>
!! \file
!! \brief CCPr Scheme for Volcanic Emissions
!!
!!
!! Reference: Benchmarking GOCART-2G in the Goddard Earth Observing System (GEOS)
!! Allison B. Collow, Peter R. Colarco, Arlindo M. da Silva, Virginie Buchard,
!! Huisheng Bian, M Chin, Sampa Das, Ravi Govindaraju, Dongchul Kim, and Valentina Aquila,
!! Geosci. Model Development, 17, 14431468, 2024
!! https://doi.org/10.5194/gmd-17-1443-2024
!!
!! \author Lacey Holland and Wei Li
!! \date 07/2024
!!!>
module CCPr_Scheme_Volcanic_GOCART_Mod

   implicit none

   private

   public :: CCPr_Scheme_Volcanic_GOCART
   public :: VolcanicEmisData
   public :: ReadASCIIPointEmis

   !> \brief VolcanicEmissionData
   !!
   !! VolcanicEmissionData is to hold volcanic emission data.
   !!
   !! \param vlat Volcano latitude
   !! \param vlon Volcano longitude
   !! \param VEmis Volcanic emissions [kg S/s]
   !! \param vbase Bottom elevation of emissions [m]
   !! \param vtop Top elevation of emissions [m]
   !! \param nPts Number of volcanic sources in the current file
   !! \param emissfile Emissions file name
   !! \param label Label for emissions
   !!
   !! \ingroup CCPr_Scheme_Volcanic_GOCART_Mod
   !!!>
   type :: VolcanicEmisData
      real :: vlat                        !volcano latitude
      real :: vlon                        !volcano longitude
      real :: VEmis                       !volcanic emissions [kg S/s]
      integer :: vbase                    !bottom elevation of emissions [m]
      integer :: vtop                     !top elevation of emissions [m]
      integer :: nPts                     !number of volcanic sources in the current file
      character(len=255) :: emissfile     !emissions file name
      character(len=255) :: label         !label for emissions
   end type VolcanicEmisData

contains

   !> \brief Brief description of the subroutine
   !!
   !! \param km          Number of vertical levels
   !! \param cdt         Model timestep [sec]
   !! \param VStart      Emissions Start time [sec]
   !! \param VEnd        Emissions end time [sec]
   !! \param nVolc       Number of volcanic sources
   !! \param iPoint      Grid cell index i of each volcanic source
   !! \param jPoint      Grid cell index j of each volcanic source
   !! \param hms         Current model time [sec]
   !! \param g0          Gravity [m/s^2]
   !! \param zbox        Geopotential Height difference [m] for layer
   !! \param delp        Pressure Thickness for layer [Pa]
   !! \param area        Area of grid cell [m^2]
   !! \param vSO2        Volcanic emissions  [kg S/s]
   !! \param nSO2        Index of SO2 relative to other sulfate tracers
   !! \param SO2         SO2 emissions [kg kg-1]
   !! \param SU_emis     SU emissions, kg/m2/s
   !! \param vCloud      Top elevation of emissions [m]
   !! \param vElev       Bottom elevation of emissions [m]
   !! \param vLat        Latitude specified in file [degree]
   !! \param VLon        Longitude specified in file [degree]
   !! \param RC           Success or Failure
   !!
   !!!>
   subroutine CCPr_Scheme_Volcanic_GOCART(km, &
      cdt, &
      VStart, &
      VEnd, &
      nVolc, &
      iPoint, &
      jPoint, &
      hms, &
      g0, &
      zbox, &
      delp, &
      area, &
      vSO2, &
      nSO2, &
      SO2, &
      vCloud, &
      vElev, &
      vLat, &
      VLon, &
      RC )

      USE GOCART2G_Process, only: SUVolcanicEmissions

      IMPLICIT NONE

      ! Arguments
      INTEGER, intent(in)                   :: km          ! number of vertical levels
      REAL, intent(in)                      :: cdt         ! model timestep [sec]
      INTEGER, intent(inout),dimension(:)   :: vStart      ! Emissions Start time [sec]
      INTEGER, intent(inout),dimension(:)   :: vEnd        ! Emissions end time [sec]
      INTEGER, intent(inout)                :: nVolc       ! number of volcanic sources
      INTEGER, intent(inout),dimension(:)   :: iPoint, jPoint ! grid cell index of each volcanic source
      INTEGER, intent(in)                   :: hms    ! current model time [sec]
      REAL, intent(in)                      :: g0
      REAL, intent(in), dimension(:) :: zbox  ! geopotential Height difference [m] for layer
      REAL, intent(in), dimension(:) :: delp   ! Pressure Thickness for layer [Pa]
      REAL, intent(inout),dimension(:,:)    :: area     ! area of grid cell [m^2]
      REAL, intent(inout),dimension(:)      :: vSO2   ! volcanic emissions  [kg S/s]
      INTEGER, intent(in)                   :: nSO2     ! index of SO2 relative to other sulfate tracers
      REAL, intent(inout),dimension(:,:,:),pointer  :: SO2       ! SO2 emissions [kg kg-1]
      REAL, intent(inout),dimension(:)        :: vCloud    ! top elevation of emissions [m]
      REAL, intent(inout),dimension(:)        :: vElev     ! bottom elevation of emissions [m]
      REAL, intent(inout),dimension(:)        :: vLat     ! latitude specified in file [degree]
      REAL, intent(inout),dimension(:)        :: VLon     ! longitude specified in file [degree]
      INTEGER, intent(inout)                  :: rc          ! error code

      !local variables
      REAL, dimension(:,:,:),pointer  :: SU_emis   ! SU emissions [kg/m2/s]
      REAL, dimension(:,:),pointer    :: SO2EMVol  ! volcanic emissions [kg m-2 s-1]
      REAL, parameter :: fMassSulfur = 32.  !  gram molecular weights of species
      REAL, parameter :: fMassSO2 = 64.     !  gram molecular weights of species
      real, pointer :: GOCART_ZBOX(:,:,:)
      real, pointer :: GOCART_DELP(:,:,:)
      character(len=256) :: errMsg
      character(len=256) :: thisLoc

      ! Initialize
      errMsg = ''
      thisLoc = ' -> at CCPr_Scheme_GOCART_SUVolcanicEmissions &
      & (in CCPr_Scheme_GOCART_SUVolcanicEmissions_mod.F90)'
      RC = 0

      ! transform data for GOCART SUVolcanicEmissions call
      call PrepMetVarsForGOCARTSUV(km,  &
         delp,                          &
         zbox,                          &
         GOCART_DELP,                   &
         GOCART_ZBOX)

      !convert SO2 unit from Kg S to Kg SO2
      vSO2 = vSO2 * fMassSO2 / fMassSulfur

      !call gocart emission function
      if (nVolc > 0) then

         !iPoint(1) = 0
         !jPoint(1) = 0

         allocate(SU_emis(1,1,nSO2)) !TODO: nSO2 =1 for now
         allocate(SO2EMVol, mold=area)
         !allocate(SO2EMVN, SO2EMVE, mold=area)

         call SUvolcanicEmissions (nVolc, vStart, vEnd, vSO2, &
            vElev, vCloud, &
            iPoint, jPoint, &
            hms, SO2EMVol, SO2, nSO2, &
            SU_emis, km, cdt, g0, gocart_ZBOX, gocart_DELP, area, &
            vLat, vLon, rc)

      end if

      if (associated(GOCART_DELP)) nullify(GOCART_DELP)
      if (associated(GOCART_zbox)) nullify(GOCART_zbox)
      if (associated(SU_emis)) nullify(SU_emis)
      if (associated(SO2EMVol)) nullify(SO2EMVol)


   end subroutine CCPr_Scheme_Volcanic_GOCART


   !> \brief Brief description of the subroutine
   !!
   !! \param filename           Emissions file name
   !! \param label              Label for emissions
   !! \param VolcanicEmissions  Volcanic emissions data
   !! \param rc                 Success or Failure
   !!
   !!!>
   subroutine ReadASCIIPointEmis (filename, label, VolcanicEmissions, rc )

      implicit none

      character(len=1055), intent(in) :: filename
      character(len=7), intent(in) :: label
      type(VolcanicEmisData), intent(inout), allocatable :: VolcanicEmissions(:)
      integer, intent(inout) :: rc
      !local variables
      integer :: num_emiss_sources=0
      integer :: num_lines=0
      integer :: num_skip=0
      integer :: i
      character(1056) :: line
      character(len=255) :: errmsg

      ! Open the file
      open(unit=10, file=filename, status='old', action='read', iostat=rc)

      if (rc /= 0) then
         print *, "Error opening file: ", filename, "  RC=", rc
         return
      end if

      ! Count the number of lines in the file
      readloop:  do while (rc >= 0)

         read(10, '(A)', iostat=rc) line
         num_lines = num_lines+1
         line = trim(line)

         if (rc /= 0) then
            print *, "Error reading file:", filename, "  RC=", rc
            return
         end if

         if (line(1:1)=="#") then
            num_skip = num_skip + 1
            continue
         else if (trim(line)==trim(label)//"::") then
            num_skip = num_skip + 1
            continue
         else if (line(1:2)=="::") then
            exit
         else
            num_emiss_sources = num_emiss_sources + 1
         end if

      end do readloop

      rewind(10)

      ! Allocate the array to hold all entries
      allocate( VolcanicEmissions(num_emiss_sources))

      do i = 1, num_skip
         read(10, '(A)', iostat=rc) line
         if (rc /= 0) return
      end do

      do i = 1, num_emiss_sources
         read(10, *, iostat=rc, iomsg=errmsg)  VolcanicEmissions(i)%vlat, &
            VolcanicEmissions(i)%vlon, &
            VolcanicEmissions(i)%vemis, &
            VolcanicEmissions(i)%vbase, &
            VolcanicEmissions(i)%vtop
         if (rc /= 0) then
            print *, "Error reading file:", trim(filename), "  RC=", rc
            print *, "Error message:", trim(errmsg)
            return
         end if
      end do

      ! Close the file and transfer data to output array
      close(10)

      VolcanicEmissions%nPts = num_emiss_sources
      VolcanicEmissions%emissfile = trim(filename)
      VolcanicEmissions%label = trim(label)

   end subroutine ReadASCIIPointEmis


   !> \brief Brief description of the subroutine
   !!
   !! \param km                Number of vertical levels
   !! \param delp              Pressure Thickness for layer [Pa]
   !! \param zbox              Geopotential Height difference [m] for layer
   !! \param GOCART_DELP       Pressure Thickness for layer in GOCART format [Pa]
   !! \param GOCART_ZBOX       Geopotential Height difference in GOCART format [m] for layer
   !!
   !!!>
   subroutine PrepMetVarsForGOCARTSUV(km, delp, zbox, GOCART_DELP, GOCART_ZBOX)

      IMPLICIT NONE

      ! INPUTS
      INTEGER, intent(in)                     :: km     ! number of vertical levels
      REAL,  intent(in), DIMENSION(:), target :: delp   ! Pressure Thickness for layer [Pa]
      REAL,  intent(in), DIMENSION(:), target :: zbox  ! Geopotential Height difference [m] for layer

      ! INPUT/OUTPUTS
      REAL, intent(inout), pointer :: GOCART_DELP(:,:,:)   !< pressure thickness for layer in GOCART format [Pa]
      REAL, intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_ZBOX  !< Geopotential Height difference in GOCART format [m] for layer

      allocate(GOCART_DELP(1, 1, km))
      allocate(GOCART_ZBOX(1, 1, km))

      GOCART_DELP(1,1,:) = delp !  pressure  in middle of layer
      GOCART_ZBOX(1,1,:) = zbox    ! mid layer geopotential height [m]

   end subroutine PrepMetVarsForGOCARTSUV


end module CCPr_Scheme_Volcanic_GOCART_Mod
