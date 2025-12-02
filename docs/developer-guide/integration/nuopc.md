# NUOPC Integration

This guide covers coupling CATChem with the National Unified Operational Prediction Capability (NUOPC) software layer used in Earth system models. In this integration, CATChem is coupled to the atmosphere (ATM) component within the UFS and CATChem acts as any other modeling component within an Earth system model.

## Overview

NUOPC integration enables CATChem to function as a component in coupled Earth system models, providing:

- **ESMF/NUOPC Compliance** - Standard NUOPC component interfaces
- **Coupling Infrastructure** - Standard field exchanges and time management
- **Multi-Model Integration** - Coupling with atmosphere, ocean, land, and ice models
- **Scalable Parallelization** - ESMF-based parallel decomposition

## NUOPC Component Architecture

### CATChem NUOPC Cap

```fortran
module CATCHEM_NUOPC_Cap
  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS    => SetServices, &
    model_label_Advance => label_Advance
  use catchem_mod
  implicit none

  private

  public :: SetServices

contains

  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! Register the generic methods
    call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

    ! Set entry point for initialization phases
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

    ! Set entry point for run phase
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_RUN, &
      phaseLabelList=(/model_label_Advance/), userRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

    ! Set entry point for finalization
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_FINALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=ModelFinalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  end subroutine SetServices

end module CATCHEM_NUOPC_Cap
```

### Initialization Phases

```fortran
subroutine InitializeP1(model, importState, exportState, clock, rc)
  type(ESMF_GridComp)  :: model
  type(ESMF_State)     :: importState, exportState
  type(ESMF_Clock)     :: clock
  integer, intent(out) :: rc

  ! Advertise import and export fields
  call NUOPC_Advertise(importState, StandardName="air_temperature", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  call NUOPC_Advertise(importState, StandardName="air_pressure", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  call NUOPC_Advertise(importState, StandardName="specific_humidity", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  call NUOPC_Advertise(exportState, StandardName="mass_fraction_of_ozone_in_air", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  call NUOPC_Advertise(exportState, StandardName="mass_fraction_of_nitrogen_dioxide_in_air", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

end subroutine InitializeP1

subroutine InitializeP2(model, importState, exportState, clock, rc)
  type(ESMF_GridComp)  :: model
  type(ESMF_State)     :: importState, exportState
  type(ESMF_Clock)     :: clock
  integer, intent(out) :: rc

  type(ESMF_Grid) :: grid
  type(ESMF_Field) :: field
  type(ESMF_Config) :: config
  character(len=ESMF_MAXSTR) :: config_file
  integer :: localPet, petCount

  ! Get ESMF context
  call ESMF_GridCompGet(model, config=config, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Read configuration
  call ESMF_ConfigGetAttribute(config, config_file, &
    label='catchem_config_file:', default='catchem_config.yml', rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Create grid (or get from parent)
  call CreateCATChemGrid(model, grid, rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Realize import fields
  call RealizeConnectedFields(importState, grid, "import", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Realize export fields
  call RealizeConnectedFields(exportState, grid, "export", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Initialize CATChem
  call InitializeCATChem(config_file, grid, rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

end subroutine InitializeP2
```

### Model Advance

```fortran
subroutine ModelAdvance(model, importState, exportState, clock, rc)
  type(ESMF_GridComp)  :: model
  type(ESMF_State)     :: importState, exportState
  type(ESMF_Clock)     :: clock
  integer, intent(out) :: rc

  type(ESMF_Time) :: currTime, stopTime
  type(ESMF_TimeInterval) :: timeStep
  type(ESMF_Field) :: field_temp, field_pres, field_humid
  type(ESMF_Field) :: field_o3, field_no2
  real(ESMF_KIND_R8), pointer :: temp_ptr(:,:,:), pres_ptr(:,:,:)
  real(ESMF_KIND_R8), pointer :: humid_ptr(:,:,:)
  real(ESMF_KIND_R8), pointer :: o3_ptr(:,:,:), no2_ptr(:,:,:)
  real(ESMF_KIND_R8) :: dt_seconds

  ! Get current time and timestep
  call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  call ESMF_TimeIntervalGet(timeStep, s_r8=dt_seconds, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Get import fields
  call ESMF_StateGet(importState, "air_temperature", field_temp, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  call ESMF_StateGet(importState, "air_pressure", field_pres, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  call ESMF_StateGet(importState, "specific_humidity", field_humid, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Get export fields
  call ESMF_StateGet(exportState, "mass_fraction_of_ozone_in_air", field_o3, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  call ESMF_StateGet(exportState, "mass_fraction_of_nitrogen_dioxide_in_air", field_no2, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Get field data pointers
  call ESMF_FieldGet(field_temp, farrayPtr=temp_ptr, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  call ESMF_FieldGet(field_pres, farrayPtr=pres_ptr, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  call ESMF_FieldGet(field_humid, farrayPtr=humid_ptr, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  call ESMF_FieldGet(field_o3, farrayPtr=o3_ptr, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  call ESMF_FieldGet(field_no2, farrayPtr=no2_ptr, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Run CATChem
  call RunCATChemNUOPC(temp_ptr, pres_ptr, humid_ptr, &
                       o3_ptr, no2_ptr, dt_seconds, rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

end subroutine ModelAdvance
```

## Configuration

### NUOPC Configuration (`nuopc_config.yml`)

```yaml
# NUOPC-specific CATChem configuration
nuopc_integration:
  enabled: true
  component_name: "CATCHEM"

  # Grid configuration
  grid:
    type: "structured"
    coordinate_system: "spherical"
    decomposition: "2d_block"

  # Field definitions with NUOPC standard names
  import_fields:
    - standard_name: "air_temperature"
      long_name: "Air Temperature"
      units: "K"
      grid_location: "center"

    - standard_name: "air_pressure"
      long_name: "Air Pressure"
      units: "Pa"
      grid_location: "center"

    - standard_name: "specific_humidity"
      long_name: "Specific Humidity"
      units: "kg kg-1"
      grid_location: "center"

  export_fields:
    - standard_name: "mass_fraction_of_ozone_in_air"
      long_name: "Ozone Mass Fraction"
      units: "kg kg-1"
      grid_location: "center"

    - standard_name: "mass_fraction_of_nitrogen_dioxide_in_air"
      long_name: "NO2 Mass Fraction"
      units: "kg kg-1"
      grid_location: "center"

  # Coupling configuration
  coupling:
    coupling_interval: 3600  # seconds
    interpolation_method: "bilinear"
    extrapolation_method: "none"

  # Performance settings
  performance:
    threading: "esmf_openmp"
    decomposition_optimization: true
    memory_allocation: "esmf_managed"
```

### ESMF/NUOPC Application Configuration

```fortran
program NUOPC_CATChem_App
  use ESMF
  use NUOPC
  use NUOPC_Driver, only: driver_routine_SS => SetServices

  implicit none

  type(ESMF_GridComp) :: driver
  integer :: rc, finalrc = ESMF_SUCCESS

  ! Initialize ESMF
  call ESMF_Initialize(defaultlogfilename="NUOPC_CATChem_App.Log", &
    logkindflag=ESMF_LOGKIND_MULTI, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Create driver component
  driver = ESMF_GridCompCreate(name="NUOPC_CATChem_Driver", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Set driver services
  call ESMF_GridCompSetServices(driver, driver_routine_SS, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Initialize driver
  call ESMF_GridCompInitialize(driver, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Run driver
  call ESMF_GridCompRun(driver, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Finalize driver
  call ESMF_GridCompFinalize(driver, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) &
    finalrc = ESMF_FAILURE

  ! Destroy driver
  call ESMF_GridCompDestroy(driver, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) &
    finalrc = ESMF_FAILURE

  ! Finalize ESMF
  call ESMF_Finalize(endflag=ESMF_END_NORMAL, rc=rc)
  if (rc /= ESMF_SUCCESS) finalrc = ESMF_FAILURE

  if (finalrc == ESMF_SUCCESS) then
    print *, "NUOPC CATChem Application completed successfully"
  else
    print *, "NUOPC CATChem Application failed"
  end if

end program NUOPC_CATChem_App
```

## Advanced Features

### Regridding and Interpolation

```fortran
subroutine SetupRegridding(srcGrid, dstGrid, regridMethod, rc)
  type(ESMF_Grid), intent(in) :: srcGrid, dstGrid
  type(ESMF_RegridMethod_Flag), intent(in) :: regridMethod
  integer, intent(out) :: rc

  type(ESMF_RouteHandle) :: routeHandle
  type(ESMF_Field) :: srcField, dstField

  ! Create fields for regridding
  srcField = ESMF_FieldCreate(srcGrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  dstField = ESMF_FieldCreate(dstGrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Store regridding weights
  call ESMF_FieldRegridStore(srcField, dstField, &
    regridmethod=regridMethod, &
    routeHandle=routeHandle, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Store route handle for later use
  call StoreRouteHandle('chemistry_regrid', routeHandle, rc)

end subroutine SetupRegridding
```

### Multi-Component Coupling

```fortran
! Driver for coupled atmosphere-chemistry system
subroutine CoupledModelDriver(model, rc)
  type(ESMF_GridComp), intent(inout) :: model
  integer, intent(out) :: rc

  type(ESMF_GridComp) :: atmComp, chemComp
  type(ESMF_CplComp) :: connector
  type(ESMF_State) :: atmExportState, chemImportState
  type(ESMF_Clock) :: clock

  ! Create atmosphere component
  atmComp = ESMF_GridCompCreate(name="ATM", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Create chemistry component (CATChem)
  chemComp = ESMF_GridCompCreate(name="CATCHEM", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Create connector for coupling
  connector = ESMF_CplCompCreate(name="ATM-CATCHEM", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Setup coupling fields and dependencies
  call SetupCouplingFields(atmComp, chemComp, connector, rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

  ! Run coupled system
  call RunCoupledSystem(atmComp, chemComp, connector, clock, rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

end subroutine CoupledModelDriver
```

## Testing and Validation

### NUOPC Compliance Testing

```fortran
program test_nuopc_compliance
  use ESMF
  use NUOPC
  use CATCHEM_NUOPC_Cap
  implicit none

  type(ESMF_GridComp) :: comp
  type(ESMF_State) :: importState, exportState
  type(ESMF_Clock) :: clock
  integer :: rc
  logical :: compliance_passed = .true.

  ! Initialize ESMF
  call ESMF_Initialize(rc=rc)

  ! Create test component
  comp = ESMF_GridCompCreate(name="CATCHEM_TEST", rc=rc)
  if (rc /= ESMF_SUCCESS) compliance_passed = .false.

  ! Test SetServices
  call SetServices(comp, rc)
  if (rc /= ESMF_SUCCESS) compliance_passed = .false.

  ! Create states
  importState = ESMF_StateCreate(name="import", stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
  exportState = ESMF_StateCreate(name="export", stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)

  ! Create clock
  call CreateTestClock(clock, rc)
  if (rc /= ESMF_SUCCESS) compliance_passed = .false.

  ! Test initialization phases
  call ESMF_GridCompInitialize(comp, importState=importState, &
    exportState=exportState, clock=clock, phase=1, rc=rc)
  if (rc /= ESMF_SUCCESS) compliance_passed = .false.

  ! Cleanup
  call ESMF_GridCompDestroy(comp, rc=rc)
  call ESMF_StateDestroy(importState, rc=rc)
  call ESMF_StateDestroy(exportState, rc=rc)
  call ESMF_ClockDestroy(clock, rc=rc)
  call ESMF_Finalize(rc=rc)

  if (compliance_passed) then
    print *, "NUOPC compliance test PASSED"
  else
    print *, "NUOPC compliance test FAILED"
    stop 1
  end if

end program test_nuopc_compliance
```

## Best Practices

### 1. ESMF Resource Management
```fortran
! Always check return codes
call ESMF_FieldGet(field, farrayPtr=data_ptr, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) return

! Properly destroy ESMF objects
call ESMF_FieldDestroy(field, rc=rc)
call ESMF_GridDestroy(grid, rc=rc)
```

### 2. Standard Names
- Use CF-compliant standard names for all fields
- Follow NUOPC naming conventions
- Document all field mappings clearly

### 3. Performance
- Minimize field copies between components
- Use ESMF regridding for efficient interpolation
- Implement proper parallel decomposition

This NUOPC integration enables CATChem to work in coupled Earth system modeling frameworks with standard interfaces and efficient coupling capabilities.
