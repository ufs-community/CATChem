MODULE ReadEmissions

!
! Author: LDH  (13 Jan 2025)
!


   implicit none

   private

   public :: ReadASCIIPointEmissions
   public :: VolcanicEmissionData


   ! Define a derived type to hold the data for each emission entry
   ! This may need to go into the Emissions State type???
   type :: VolcanicEmissionData
      real :: vlat
      real :: vlon
      real :: VEmis
      integer :: vbase
      integer :: vtop
      integer :: nPts
      character(len=255) :: emissfile
      character(len=255) :: label
   end type VolcanicEmissionData

   !integer :: rc
   !character(len=1055) :: filename
   !character(len=1055) :: label

   !type(VolcanicEmissionData), allocatable :: VolcanicEmiss(:)

   !filename="./so2_volcanic_emissions_Carns.20220101.rc"
   !label="volcano"
   !call ReadASCIIPointEmissions(filename, label, VolcanicEmiss, rc )

contains



   subroutine ReadASCIIPointEmissions (filename, label, VolcanicEmissions, rc )

      implicit none

      integer, intent(out) :: rc
      integer :: num_emiss_sources=0
      integer :: num_lines=0
      integer :: num_skip=0
      integer :: i
      character(1056) :: line
      character(len=1055), intent(in) :: filename
      character(len=7), intent(in) :: label
      character(len=255) :: errmsg

      type(VolcanicEmissionData), allocatable :: VolcanicEmissions(:)
      type(VolcanicEmissionData), allocatable :: temp_emissions(:)
      !allocate(VolcanicEmiss(9))


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
      allocate(temp_emissions(num_emiss_sources))

      do i = 1, num_skip
         read(10, '(A)', iostat=rc) line
         if (rc /= 0) return
      end do

      do i = 1, num_emiss_sources
         read(10, *, iostat=rc, iomsg=errmsg) temp_emissions(i)%vlat, &
            temp_emissions(i)%vlon, &
            temp_emissions(i)%vemis, &
            temp_emissions(i)%vbase, &
            temp_emissions(i)%vtop
         if (rc /= 0) then
            print *, "Error reading file:", trim(filename), "  RC=", rc
            print *, "Error message:", trim(errmsg)
            return
         end if
      end do

      ! Close the file and transfer data to output array
      close(10)

      temp_emissions%nPts = num_emiss_sources
      temp_emissions%emissfile = trim(filename)
      temp_emissions%label = trim(label)

      VolcanicEmissions = temp_emissions

   end subroutine ReadASCIIPointEmissions


end MODULE ReadEmissions
