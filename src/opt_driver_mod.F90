MODULE opt_driver_mod
!**********************************************************************************  
! This computer software was prepared by Battelle Memorial Institute, hereinafter
! the Contractor, under Contract No. DE-AC05-76RL0 1830 with the Department of 
! Energy (DOE). NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!**********************************************************************************  
!
! WRF-chem V3.0 : Original version of optical_driver written by Jerome Fast (PNNL)
!                 and James Barnard (PNNL)
!
!WRF:MODEL_LAYER:CHEMISTRY
!
  use catchem_constants ,        only : kind_chem

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: optical_driver

CONTAINS
      SUBROUTINE optical_driver(curr_secs,dtstep,&
               chem,dz8w,alt,relhum,                                     &
              h2oai,h2oaj,                                              &
               tauaersw,gaersw,waersw,bscoefsw,tauaerlw,                 &
               l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,                      &
               num_chem,chem_opt,ids,ide, jds,jde, kds,kde,              &
               ims,ime, jms,jme, kms,kme,                                &
               its,ite, jts,jte, kts,kte                                 )

!------------------------------------------------------------------------

   USE opt_averaging_mod

   INTEGER,      INTENT(IN   ) :: chem_opt,num_chem,ids,ide,           &
                                           jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte
   REAL(KIND=kind_chem), INTENT(IN   ) :: curr_secs
   REAL(kind_chem),      INTENT(IN   ) :: dtstep
!
! array that holds all advected chemical species
!
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),             &
         INTENT(INOUT ) ::  chem
!
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme ),                       &  
         INTENT(IN ) ::  relhum, dz8w, alt, h2oai, h2oaj
!
! arrays that hold the aerosol optical properties
!
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme, 1:4 ),                       &  
         INTENT(INOUT ) ::                                             &
           tauaersw,gaersw,waersw,bscoefsw
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme, 1:16),                       &  
         INTENT(INOUT ) ::                                             &
           tauaerlw
   REAL(kind_chem), DIMENSION( ims:ime, kms:kme, jms:jme, 1:4 ),                  &  
         INTENT(INOUT ) ::                                             &
           l2aer, l3aer, l4aer, l5aer, l6aer, l7aer
!

!         
! local variables
!
      logical processingAerosols
      integer nbin_o
      integer option_method, option_mie

!-----------------------------------------------------------------
! compute only if simulating aerosols and aer_ra_feedback=1

!  IF (config_flags%aer_ra_feedback .eq. 0) THEN
!        call wrf_debug(15,'no feedback, return from optical driver')
!    return
!  ENDIF
!  select case (config_flags%chem_opt)
!  case ( RADM2SORG,           RADM2SORG_KPP,      RADM2SORG_AQ, &
!         GOCART_SIMPLE,       RACMSORG_KPP,       RACMSORG_AQ,  &
!         CBMZ_MOSAIC_4BIN,    CBMZ_MOSAIC_8BIN,   &
!         CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ )
      processingAerosols = .true.
!     call wrf_debug(15,'optical driver: process aerosols true')
!  case default
!     processingAerosols = .false.
!     call wrf_debug(15,'optical driver: process aerosols false')
!  end select

  if( processingAerosols ) then
!
! select aerosol optical property option
! VOLUME: volume averaging of refractive indicies
! * for MADE/SORGAM, assume same 8 size bins as MOSAIC by default
! SHELL: shell-core approach, placeholder
!
!  select case (config_flags%chem_opt)
!  case ( RADM2SORG,           RADM2SORG_KPP,      RADM2SORG_AQ, &
!         GOCART_SIMPLE,       RACMSORG_KPP,       RACMSORG_AQ   )
     nbin_o = 8
!  case (CBMZ_MOSAIC_4BIN,    CBMZ_MOSAIC_8BIN,   &
!        CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ )
!    nbin_o = nbin_a
!  end select
!
!    call wrf_debug(15,'optical averaging')
!    aer_op_opt_select: SELECT CASE(config_flags%aer_op_opt)
!    CASE (VOLUME_APPROX)
       option_method=1
       option_mie=1
!    CASE (MAXWELL_APPROX)
!      option_method=2
!      option_mie=1
!    CASE (VOLUME_EXACT)
!      option_method=1
!      option_mie=2
!    CASE (MAXWELL_EXACT)
!      option_method=2
!      option_mie=2
!    CASE (SHELL_EXACT)
!      option_method=3
!      option_mie=2
!    CASE DEFAULT
!       if( config_flags%aer_op_opt > 0 ) then
!          call wrf_message('WARNING: Invalid aer_op_opt. Defaulting to VOLUME_APPROX.')
!          option_method=1
!          option_mie=1
!       end if
!    END SELECT aer_op_opt_select

!    if( config_flags%aer_op_opt > 0 ) then
!       call wrf_debug(15,'optical driver: call optical averaging')
        call optical_averaging(curr_secs,dtstep,                     &
             nbin_o,option_method,option_mie,chem,dz8w,alt,  &
             relhum,h2oai,h2oaj,                                     &
             tauaersw,gaersw,waersw,bscoefsw,tauaerlw,               &
             l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,                    &
             num_chem,chem_opt,ids,ide, jds,jde, kds,kde,            &
             ims,ime, jms,jme, kms,kme,                              &
             its,ite, jts,jte, kts,kte                               )
!    else
        !If aer_op_opt==0 then the optical arrays are already set to
        !zero in chemics_init so there will not be a problem if the
        !user has selected aer_ra_feedback=1.
!    end if
!
   endif
   return

END SUBROUTINE optical_driver
END MODULE opt_driver_mod
