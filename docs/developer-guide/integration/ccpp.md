# CCPP Integration

This guide covers integrating CATChem inline with the Common Community Physics Package (CCPP) framework used in the UFS directly in the atmosphere (ATM) component. This option is important for research as some processes require close interactions between chemistry and physics. Considering the interactive and interdependent nature of chemistry and physics, it is useful for certain applications to insert the aerosol or chemistry modules directly inside the physics suites.

## Overview

CCPP integration enables CATChem to work seamlessly with the CCPP framework, providing:

- **Standardized Interface** - CCPP-compliant subroutine signatures and metadata
- **Automatic Integration** - CCPP-generated calling sequences and type definitions
- **Host Model Compatibility** - Works with any CCPP-enabled host model
- **Performance Optimization** - Optimized data structures and memory layout

## CCPP Interface Architecture

### CCPP Wrapper Module

```fortran
module ccpp_catchem_interface
  use ccpp_types, only: ccpp_t
  use ccpp_fields, only: ccpp_field_t
  use catchem_mod
  implicit none
  private

  public :: ccpp_catchem_init
  public :: ccpp_catchem_run
  public :: ccpp_catchem_finalize

contains

  subroutine ccpp_catchem_init(ccpp_data, errmsg, errflg)
    type(ccpp_t), intent(inout) :: ccpp_data
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! Initialize CATChem using CCPP data
    call catchem_ccpp_init(ccpp_data, errmsg, errflg)
  end subroutine ccpp_catchem_init

  subroutine ccpp_catchem_run(ccpp_data, dt, errmsg, errflg)
    type(ccpp_t), intent(inout) :: ccpp_data
    real(kind_phys), intent(in) :: dt
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! Run CATChem processes
    call catchem_ccpp_run(ccpp_data, dt, errmsg, errflg)
  end subroutine ccpp_catchem_run

  subroutine ccpp_catchem_finalize(ccpp_data, errmsg, errflg)
    type(ccpp_t), intent(inout) :: ccpp_data
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! Finalize CATChem
    call catchem_ccpp_finalize(ccpp_data, errmsg, errflg)
  end subroutine ccpp_catchem_finalize

end module ccpp_catchem_interface
```

### CCPP Metadata File (`ccpp_catchem_interface.meta`)

```fortran
[ccpp-table-properties]
  name = ccpp_catchem_interface
  type = scheme
  dependencies = catchem_mod.F90,catchem_types.F90

########################################################################
[ccpp-arg-table]
  name = ccpp_catchem_init
  type = scheme
[im]
  standard_name = horizontal_loop_extent
  long_name = horizontal loop extent
  units = count
  dimensions = ()
  type = integer
  intent = in
[km]
  standard_name = vertical_layer_dimension
  long_name = vertical layer dimension
  units = count
  dimensions = ()
  type = integer
  intent = in
[dt]
  standard_name = timestep_for_chemistry
  long_name = chemistry timestep
  units = s
  dimensions = ()
  type = real
  kind = kind_phys
  intent = in
[temperature]
  standard_name = air_temperature
  long_name = layer mean air temperature
  units = K
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real
  kind = kind_phys
  intent = in
[pressure]
  standard_name = air_pressure
  long_name = layer mean air pressure
  units = Pa
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real
  kind = kind_phys
  intent = in
[specific_humidity]
  standard_name = specific_humidity
  long_name = water vapor specific humidity
  units = kg kg-1
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real
  kind = kind_phys
  intent = in
[o3_mixing_ratio]
  standard_name = ozone_mixing_ratio
  long_name = ozone mixing ratio
  units = kg kg-1
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real
  kind = kind_phys
  intent = inout
[no2_mixing_ratio]
  standard_name = nitrogen_dioxide_mixing_ratio
  long_name = nitrogen dioxide mixing ratio
  units = kg kg-1
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real
  kind = kind_phys
  intent = inout
[errmsg]
  standard_name = ccpp_error_message
  long_name = error message for error handling in CCPP
  units = none
  dimensions = ()
  type = character
  kind = len=*
  intent = out
[errflg]
  standard_name = ccpp_error_code
  long_name = error code for error handling in CCPP
  units = 1
  dimensions = ()
  type = integer
  intent = out
```

## Implementation

### CATChem CCPP Wrapper Implementation

```fortran
module catchem_ccpp_wrapper
  use ccpp_types
  use catchem_mod
  use catchem_types
  implicit none

contains

  subroutine catchem_ccpp_init(ccpp_data, errmsg, errflg)
    type(ccpp_t), intent(inout) :: ccpp_data
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    type(catchem_config_type) :: config
    integer :: im, km
    integer :: rc

    errflg = 0
    errmsg = ''

    ! Extract dimensions from CCPP data
    call ccpp_field_get(ccpp_data, 'horizontal_loop_extent', im, rc)
    if (rc /= 0) then
      errflg = 1
      errmsg = 'Failed to get horizontal_loop_extent from CCPP'
      return
    end if

    call ccpp_field_get(ccpp_data, 'vertical_layer_dimension', km, rc)
    if (rc /= 0) then
      errflg = 1
      errmsg = 'Failed to get vertical_layer_dimension from CCPP'
      return
    end if

    ! Initialize CATChem configuration
    call setup_catchem_config_from_ccpp(ccpp_data, config, rc)
    if (rc /= 0) then
      errflg = 1
      errmsg = 'Failed to setup CATChem configuration'
      return
    end if

    ! Initialize CATChem
    call catchem_init(config, im, km, rc)
    if (rc /= 0) then
      errflg = 1
      errmsg = 'CATChem initialization failed'
      return
    end if

  end subroutine catchem_ccpp_init

  subroutine catchem_ccpp_run(ccpp_data, dt, errmsg, errflg)
    type(ccpp_t), intent(inout) :: ccpp_data
    real(kind_phys), intent(in) :: dt
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    real(kind_phys), pointer :: temperature(:,:) => null()
    real(kind_phys), pointer :: pressure(:,:) => null()
    real(kind_phys), pointer :: specific_humidity(:,:) => null()
    real(kind_phys), pointer :: o3_mixing_ratio(:,:) => null()
    real(kind_phys), pointer :: no2_mixing_ratio(:,:) => null()

    type(catchem_state_type) :: catchem_state
    integer :: rc

    errflg = 0
    errmsg = ''

    ! Get meteorological fields from CCPP
    call ccpp_field_get(ccpp_data, 'air_temperature', temperature, rc)
    if (rc /= 0) then
      errflg = 1
      errmsg = 'Failed to get air_temperature from CCPP'
      return
    end if

    call ccpp_field_get(ccpp_data, 'air_pressure', pressure, rc)
    if (rc /= 0) then
      errflg = 1
      errmsg = 'Failed to get air_pressure from CCPP'
      return
    end if

    call ccpp_field_get(ccpp_data, 'specific_humidity', specific_humidity, rc)
    if (rc /= 0) then
      errflg = 1
      errmsg = 'Failed to get specific_humidity from CCPP'
      return
    end if

    ! Get chemical species from CCPP
    call ccpp_field_get(ccpp_data, 'ozone_mixing_ratio', o3_mixing_ratio, rc)
    if (rc /= 0) then
      errflg = 1
      errmsg = 'Failed to get ozone_mixing_ratio from CCPP'
      return
    end if

    call ccpp_field_get(ccpp_data, 'nitrogen_dioxide_mixing_ratio', no2_mixing_ratio, rc)
    if (rc /= 0) then
      errflg = 1
      errmsg = 'Failed to get nitrogen_dioxide_mixing_ratio from CCPP'
      return
    end if

    ! Setup CATChem state from CCPP data
    call setup_catchem_state_from_ccpp( &
      temperature, pressure, specific_humidity, &
      o3_mixing_ratio, no2_mixing_ratio, &
      catchem_state, rc)
    if (rc /= 0) then
      errflg = 1
      errmsg = 'Failed to setup CATChem state'
      return
    end if

    ! Run CATChem processes
    call catchem_run(catchem_state, dt, rc)
    if (rc /= 0) then
      errflg = 1
      errmsg = 'CATChem run failed'
      return
    end if

    ! Update CCPP fields with CATChem results
    call update_ccpp_from_catchem_state( &
      catchem_state, &
      o3_mixing_ratio, no2_mixing_ratio, &
      rc)
    if (rc /= 0) then
      errflg = 1
      errmsg = 'Failed to update CCPP fields'
      return
    end if

  end subroutine catchem_ccpp_run

end module catchem_ccpp_wrapper
```

## Configuration

### CCPP Suite Definition File (`suite_FV3_CATChem.xml`)

```xml
<?xml version="1.0" encoding="UTF-8"?>
<suite name="FV3_CATChem" version="1">
  <group name="time_vary">
    <subcycle loop="1">
      <scheme>GFS_time_vary_pre</scheme>
      <scheme>GFS_rrtmg_setup</scheme>
      <scheme>GFS_rad_time_vary</scheme>
      <scheme>GFS_phys_time_vary</scheme>
    </subcycle>
  </group>

  <group name="radiation">
    <subcycle loop="1">
      <scheme>GFS_suite_interstitial_rad_reset</scheme>
      <scheme>GFS_rrtmg_pre</scheme>
      <scheme>rrtmg_sw</scheme>
      <scheme>rrtmg_lw</scheme>
      <scheme>GFS_rrtmg_post</scheme>
    </subcycle>
  </group>

  <group name="physics">
    <subcycle loop="1">
      <scheme>GFS_suite_interstitial_phys_reset</scheme>
      <scheme>GFS_suite_stateout_reset</scheme>
      <scheme>get_prs_fv3</scheme>
      <scheme>GFS_suite_interstitial_1</scheme>
      <scheme>GFS_surface_generic_pre</scheme>
      <scheme>GFS_surface_composites_pre</scheme>
      <scheme>dcyc2t3</scheme>
      <scheme>GFS_surface_composites_inter</scheme>
      <scheme>GFS_suite_interstitial_2</scheme>
    </subcycle>

    <!-- Chemistry Integration -->
    <subcycle loop="1">
      <scheme>ccpp_catchem_interface</scheme>
    </subcycle>

    <subcycle loop="1">
      <scheme>GFS_suite_interstitial_3</scheme>
      <scheme>GFS_suite_stateout_update</scheme>
      <scheme>GFS_suite_interstitial_4</scheme>
      <scheme>GFS_MP_generic_pre</scheme>
      <scheme>mp_thompson_pre</scheme>
      <scheme>mp_thompson</scheme>
      <scheme>mp_thompson_post</scheme>
      <scheme>GFS_MP_generic_post</scheme>
    </subcycle>
  </group>

  <group name="stochastics">
    <subcycle loop="1">
      <scheme>GFS_stochastics</scheme>
    </subcycle>
  </group>
</suite>
```

### CCPP Configuration

```yaml
# CCPP-specific CATChem configuration
ccpp_integration:
  enabled: true

  # Field mapping between CCPP standard names and CATChem species
  field_mapping:
    air_temperature: temperature
    air_pressure: pressure
    specific_humidity: humidity
    ozone_mixing_ratio: O3
    nitrogen_dioxide_mixing_ratio: NO2
    sulfur_dioxide_mixing_ratio: SO2

  # Process configuration for CCPP integration
  processes:
    - name: chemistry
      enabled: true
      ccpp_compliant: true
      timestep_coupling: synchronous

    - name: settling
      enabled: true
      ccpp_compliant: true

    - name: dry_deposition
      enabled: false  # Handled by host model

  # CCPP-specific optimizations
  optimizations:
    column_processing: true
    memory_layout: ccpp_optimized
    data_locality: true
```

## Advanced Integration Features

### Field Mapping and Unit Conversion

```fortran
subroutine setup_ccpp_field_mapping(field_map)
  type(field_mapping_type), intent(out) :: field_map

  ! Meteorological fields
  call add_field_mapping(field_map, &
    ccpp_name='air_temperature', &
    catchem_name='temperature', &
    units_conversion=1.0_fp, &  ! K to K (no conversion)
    required=.true.)

  call add_field_mapping(field_map, &
    ccpp_name='air_pressure', &
    catchem_name='pressure', &
    units_conversion=1.0_fp, &  ! Pa to Pa (no conversion)
    required=.true.)

  ! Chemical species with unit conversion
  call add_field_mapping(field_map, &
    ccpp_name='ozone_mixing_ratio', &
    catchem_name='O3', &
    units_conversion=1.0_fp / molar_mass_o3, &  ! kg/kg to mol/mol
    required=.true.)

  call add_field_mapping(field_map, &
    ccpp_name='nitrogen_dioxide_mixing_ratio', &
    catchem_name='NO2', &
    units_conversion=1.0_fp / molar_mass_no2, &  ! kg/kg to mol/mol
    required=.true.)
end subroutine setup_ccpp_field_mapping
```

### Performance Optimization

```fortran
! CCPP-optimized column processing
subroutine catchem_ccpp_column_processing(ccpp_data, dt)
  type(ccpp_t), intent(inout) :: ccpp_data
  real(kind_phys), intent(in) :: dt

  integer :: im, km, i
  real(kind_phys), pointer :: temperature(:,:), pressure(:,:)
  real(kind_phys), pointer :: o3(:,:), no2(:,:)

  ! Get field pointers (more efficient than copying)
  call ccpp_field_get_ptr(ccpp_data, 'air_temperature', temperature)
  call ccpp_field_get_ptr(ccpp_data, 'air_pressure', pressure)
  call ccpp_field_get_ptr(ccpp_data, 'ozone_mixing_ratio', o3)
  call ccpp_field_get_ptr(ccpp_data, 'nitrogen_dioxide_mixing_ratio', no2)

  im = size(temperature, 1)
  km = size(temperature, 2)

  ! Process columns in parallel
  !$OMP PARALLEL DO PRIVATE(i)
  do i = 1, im
    call catchem_process_column( &
      temperature(i,:), pressure(i,:), &
      o3(i,:), no2(i,:), &
      dt)
  end do
  !$OMP END PARALLEL DO

end subroutine catchem_ccpp_column_processing
```

## Testing and Validation

### CCPP Integration Tests

```fortran
program test_ccpp_integration
  use ccpp_types
  use ccpp_catchem_interface
  implicit none

  type(ccpp_t) :: ccpp_data
  character(len=512) :: errmsg
  integer :: errflg
  real(kind_phys) :: dt = 60.0_fp

  ! Setup test CCPP data
  call setup_test_ccpp_data(ccpp_data)

  ! Test initialization
  call ccpp_catchem_init(ccpp_data, errmsg, errflg)
  if (errflg /= 0) then
    print *, 'CCPP init failed: ', trim(errmsg)
    stop 1
  end if

  ! Test run
  call ccpp_catchem_run(ccpp_data, dt, errmsg, errflg)
  if (errflg /= 0) then
    print *, 'CCPP run failed: ', trim(errmsg)
    stop 1
  end if

  ! Test finalization
  call ccpp_catchem_finalize(ccpp_data, errmsg, errflg)
  if (errflg /= 0) then
    print *, 'CCPP finalize failed: ', trim(errmsg)
    stop 1
  end if

  print *, 'CCPP integration tests passed!'

end program test_ccpp_integration
```

## Best Practices

### 1. CCPP Compliance
- Follow CCPP metadata standards exactly
- Use CCPP standard names for all fields
- Implement proper error handling with errmsg/errflg

### 2. Performance
- Use field pointers instead of copying data
- Implement column-wise processing for parallelization
- Minimize memory allocations in run routine

### 3. Error Handling
```fortran
! Always check CCPP field operations
call ccpp_field_get(ccpp_data, 'field_name', field_ptr, rc)
if (rc /= 0) then
  errflg = 1
  errmsg = 'Failed to get field: field_name'
  return
end if
```

CATChem integrated inline with the CCPP framework provides close interactions between chemistry and physics while maintaining performance and standards compliance.
