!
! Haiqin.Li@noaa.gov  
! 01/2020
! constant parameters and chemistry configurations and tracers
! (This will be splited into three subroutine for configuration, constant and tracers later)
! 08/2020 move configuration into chem nml
! Jian.He@noaa.gov, 05/2023
! Move to parameters folder for CATChem

module catchem_config

  use catchem_constants, only : kind_chem

  implicit none

  !-- constant paramters
  real(kind=kind_chem), parameter :: epsilc     = 1.e-30

  !-- chemistyr module configurations
  integer :: chem_opt = 300
  integer :: kemit = 1
  integer :: dust_opt = 5
  integer :: dmsemis_opt = 1
  integer :: seas_opt = 2
  integer :: biomass_burn_opt=1
  integer :: plumerise_flag = 2  ! 1=MODIS, 2=GBBEPx
  integer :: plumerisefire_frq=60
  integer :: chem_conv_tr  = 0
  integer :: aer_ra_feedback=1 !0
  integer :: aer_ra_frq=60
  integer :: wetdep_ls_opt = 1

  real(kind=kind_chem), parameter :: depo_fact=0.
  integer, parameter :: CHEM_OPT_GOCART= 300
  integer, parameter :: CHEM_OPT_GOCART_RACM  = 301
  integer, parameter :: CHEM_OPT_RACM_SOA_VBS = 108
  INTEGER, PARAMETER :: gocartracm_kpp = 301
  integer, parameter :: chem_tune_tracers = 20
  integer, parameter :: CHEM_OPT_MAX      = 500
  integer, parameter :: DUST_OPT_NONE = 0
  integer, parameter :: SEAS_OPT_NONE = 0
  ! -- DMS emissions
  integer, parameter :: DMSE_OPT_NONE   = 0
  integer, parameter :: DMSE_OPT_ENABLE = 1
  ! -- subgrid convective transport
  integer, parameter :: CTRA_OPT_NONE  = 0
  integer, parameter :: CTRA_OPT_GRELL = 2
  ! -- large scale wet deposition
  integer, parameter :: WDLS_OPT_NONE  = 0
  integer, parameter :: WDLS_OPT_GSD   = 1
  integer, parameter :: WDLS_OPT_NGAC  = 2

  ! --
  integer, parameter :: call_chemistry     = 1
  integer, parameter :: num_ebu            = 7
  integer, parameter :: num_ebu_in         = 7
  integer, parameter :: num_moist=3, num_chem=20, num_emis_seas=5, num_emis_dust=5
  integer, parameter :: num_emis_ant = 7, num_emis_vol =4

  integer, parameter :: SEAS_OPT_DEFAULT = 1

  integer, parameter :: DUST_OPT_GOCART  = 1
  integer, parameter :: DUST_OPT_AFWA    = 3
  integer, parameter :: DUST_OPT_FENGSHA = 5

  ! -- biomass burning emissions
  integer, parameter :: BURN_OPT_NONE   = 0
  integer, parameter :: BURN_OPT_ENABLE = 1
  integer, parameter :: FIRE_OPT_NONE   = 0
  integer, parameter :: FIRE_OPT_MODIS  = 1
  integer, parameter :: FIRE_OPT_GBBEPx = 2

  ! -- hydrometeors
  integer, parameter :: p_qv=1
  integer, parameter :: p_qc=2
  integer, parameter :: p_qi=3
  ! -- set pointers to predefined atmospheric tracers
  ! -- FV3 GFDL microphysics
  integer, parameter :: p_atm_shum = 1
  integer, parameter :: p_atm_cldq = 2
  integer, parameter :: p_atm_o3mr = 7

  integer :: numgas = 0


  !integer :: num_plume_data  


  real(kind=kind_chem) :: wetdep_ls_alpha(chem_tune_tracers)=-999.


  !-- tracers
  integer, parameter :: p_so2=1
  integer, parameter :: p_sulf=2
  integer, parameter :: p_dms=3
  integer, parameter :: p_msa=4
  integer, parameter :: p_p25=5
  integer, parameter :: p_bc1=6
  integer, parameter :: p_bc2=7
  integer, parameter :: p_oc1=8
  integer, parameter :: p_oc2=9
  integer, parameter :: p_dust_1=10
  integer, parameter :: p_dust_2=11
  integer, parameter :: p_dust_3=12
  integer, parameter :: p_dust_4=13
  integer, parameter :: p_dust_5=14
  integer, parameter :: p_seas_1=15
  integer, parameter :: p_seas_2=16
  integer, parameter :: p_seas_3=17
  integer, parameter :: p_seas_4=18
  integer, parameter :: p_seas_5=19
  integer, parameter :: p_p10   =20

  integer, parameter :: p_edust1=1,p_edust2=2,p_edust3=3,p_edust4=4,p_edust5=5
  integer, parameter :: p_eseas1=1,p_eseas2=2,p_eseas3=3,p_eseas4=4,p_eseas5=5
 
  integer :: p_ho=0,p_h2o2=0,p_no3=0

  ! constants
  real(kind=kind_chem), PARAMETER :: airmw      = 28.97
  real(kind=kind_chem), PARAMETER :: mw_so2_aer = 64.066
  real(kind=kind_chem), PARAMETER :: mw_so4_aer = 96.066
  real(kind=kind_chem), parameter :: smw        = 32.00
  real(kind=kind_chem), parameter :: mwdry      = 28.
!  <mw>d is the molecular weight of dry air (28.966), <mw>w/<mw>d = 0.62197, and
!  (<mw>d - <mw>w)/<mw>d = 0.37803
!  http://atmos.nmsu.edu/education_and_outreach/encyclopedia/humidity.htm
  !-----------------------------------------------------------------------  
  ! tracer info
  !-----------------------------------------------------------------------
  ! Tracer index:
  ! default initialization for all sulfur and carbon species is 0 (undefined)
  !     1. DMS       = Dimethyl sulfide            = CH3SCH3  
  !     2. SO2       = Sulfur dioxide              = SO2               
  !     3. SO4       = Sulfate                     = SO4            
  !     4. MSA       = Methane sulfonic acid       = CH3SO3H             
  integer            :: ndms=1, nso2=2, nso4=3, nmsa=4

  ! optical
  integer :: nbands   = 14
  integer :: nbandlw  = 16
  REAL,    PARAMETER :: oc_mfac = 1.8 
  REAL,    PARAMETER :: nh4_mfac = 1.375


  INTEGER, PARAMETER :: cbmz_mosaic_4bin_aq = 9
  INTEGER, PARAMETER :: cbmz_mosaic_8bin_aq = 10
  INTEGER, PARAMETER :: radm2sorg_aq = 11
  INTEGER, PARAMETER :: racmsorg_aq = 12
  INTEGER, PARAMETER :: mozart_kpp = 111
  INTEGER, PARAMETER :: mozart_mosaic_4bin_kpp = 201
  INTEGER, PARAMETER :: mozart_mosaic_4bin_aq_kpp = 202
  INTEGER, PARAMETER :: radm2 = 1
  INTEGER, PARAMETER :: radm2sorg = 2
  INTEGER, PARAMETER :: radm2sorg_aqchem = 41
  INTEGER, PARAMETER :: cbmz = 5
  INTEGER, PARAMETER :: cbmz_bb = 6
  INTEGER, PARAMETER :: cbmz_bb_kpp = 120
  INTEGER, PARAMETER :: cbmz_mosaic_kpp = 170
  INTEGER, PARAMETER :: cbmz_mosaic_4bin = 7
  INTEGER, PARAMETER :: cbmz_mosaic_8bin = 8
  INTEGER, PARAMETER :: cbmz_mosaic_dms_4bin_aq = 32
  INTEGER, PARAMETER :: cbmz_mosaic_dms_8bin_aq = 34
  INTEGER, PARAMETER :: cbmz_mosaic_dms_4bin = 31
  INTEGER, PARAMETER :: cbmz_mosaic_dms_8bin = 33
  INTEGER, PARAMETER :: cbmz_cam_mam3_noaq = 501
  INTEGER, PARAMETER :: cbmz_cam_mam3_aq = 503
  INTEGER, PARAMETER :: cbmz_cam_mam7_noaq = 502
  INTEGER, PARAMETER :: cbmz_cam_mam7_aq = 504
  INTEGER, PARAMETER :: mozcart_kpp = 112
  INTEGER, PARAMETER :: racmsoavbs_kpp = 108
  INTEGER, PARAMETER :: cbm4_kpp = 110
  INTEGER, PARAMETER :: gocartradm2 = 303
  INTEGER, PARAMETER :: crimech_kpp = 600
  INTEGER, PARAMETER :: cri_mosaic_8bin_aq_kpp = 601
  INTEGER, PARAMETER :: cri_mosaic_4bin_aq_kpp = 611
  INTEGER, PARAMETER :: cb05_sorg_aq_kpp = 131
  INTEGER, PARAMETER :: cb05_sorg_vbs_aq_kpp = 132

  !-- plumerise 
  integer, parameter :: p_e_bc  =1
  integer, parameter :: p_e_oc  =2
  integer, parameter :: p_e_sulf=3
  integer, parameter :: p_e_pm_25=4
  integer, parameter :: p_e_so2=5
  integer, parameter :: p_e_pm_10=6
  integer, parameter :: p_e_dms=7
  integer, parameter :: p_ebu_bc  =1
  integer, parameter :: p_ebu_oc  =2
  integer, parameter :: p_ebu_sulf=3
  integer, parameter :: p_ebu_pm25=4
  integer, parameter :: p_ebu_so2=5
  integer, parameter :: p_ebu_pm10=6
  integer, parameter :: p_ebu_dms=7
  integer, parameter :: p_ebu_in_bc  =1
  integer, parameter :: p_ebu_in_oc  =2
  integer, parameter :: p_ebu_in_sulf=3
  integer, parameter :: p_ebu_in_pm25=4
  integer, parameter :: p_ebu_in_so2=5
  integer, parameter :: p_ebu_in_pm10=6
  integer, parameter :: p_ebu_in_dms=7


  !-- gocartracm
  integer, parameter :: p_ebu_no = 3
  integer, parameter :: p_ebu_no2=4
  integer, parameter :: p_ebu_co =5
  integer, parameter :: p_ebu_eth=6
  integer, parameter :: p_ebu_hc3=7
  integer, parameter :: p_ebu_hc5=8
  integer, parameter :: p_ebu_hc8=9
  integer, parameter :: p_ebu_xyl=10
  integer, parameter :: p_ebu_olt=11
  integer, parameter :: p_ebu_oli=12
  integer, parameter :: p_ebu_tol=13
  integer, parameter :: p_ebu_csl=14
  integer, parameter :: p_ebu_hcho=15
  integer, parameter :: p_ebu_in_no = 3
  integer, parameter :: p_ebu_in_no2=4
  integer, parameter :: p_ebu_in_co =5
  integer, parameter :: p_ebu_in_eth=6
  integer, parameter :: p_ebu_in_hc3=7
  integer, parameter :: p_ebu_in_hc5=8
  integer, parameter :: p_ebu_in_hc8=9
  integer, parameter :: p_ebu_in_xyl=10
  integer, parameter :: p_ebu_in_olt=11
  integer, parameter :: p_ebu_in_oli=12
  integer, parameter :: p_ebu_in_tol=13
  integer, parameter :: p_ebu_in_csl=14
  integer, parameter :: p_ebu_in_hcho=15

  INTEGER , PARAMETER :: PARAM_FIRST_SCALAR = 2

  integer :: p_qr  = 1
  integer :: p_pan = 1
  integer :: p_o3  = 1
  integer :: p_nh3 = 1

  integer :: p_tr2 = 0

  integer :: p_vash_1      = 0
  integer :: p_vash_10     = 0
  integer :: p_vash_2      = 0
  integer :: p_vash_3      = 0
  integer :: p_vash_4      = 0
  integer :: p_vash_5      = 0
  integer :: p_vash_6      = 0
  integer :: p_vash_7      = 0
  integer :: p_vash_8      = 0
  integer :: p_vash_9      = 0

  integer :: p_e_vash1     = 1
  integer :: p_e_vash2     = 2
  integer :: p_e_vash3     = 3
  integer :: p_e_vash4     = 4
  !--- optical properties
  integer :: num_ext_coef    = 5
  integer :: num_bscat_coef  = 3
  integer :: num_asym_par    = 3
  integer :: p_extcof3    = 1
  integer :: p_extcof55   = 2
  integer :: p_extcof106  = 3
  integer :: p_extcof3_5  = 4
  integer :: p_extcof8_12 = 5
  integer :: p_bscof3     = 1
  integer :: p_bscof55    = 2
  integer :: p_bscof106   = 3
  integer :: p_asympar3   = 1
  integer :: p_asympar55  = 2
  integer :: p_asympar106 = 3

    ! -- fire options
    !select case (plumerise_flag)
    !  case (FIRE_OPT_NONE)
    !    ! -- valid option
    !  case (FIRE_OPT_MODIS)
   integer::     num_plume_data = 8
    !  case (FIRE_OPT_GBBEPx)
    !    ! -- valid option
    !    num_plume_data = 1
    !  case default
    !    return
    !end select


end module catchem_config
