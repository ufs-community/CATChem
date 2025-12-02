# Integration Guide

This section covers integrating CATChem with external modeling systems.

!!! note "We are actively integrating CATChem into the UFS. Please check back frequently for updates to this Integration Guide. Further updates on this section will be coming soon!"


## Overview

CATChem is designed to integrate seamlessly with various atmospheric modeling frameworks:

- **[CCPP Integration](ccpp.md)**: Common Community Physics Package
- **[NUOPC Integration](nuopc.md)**: National Unified Operational Prediction Capability
- **Standalone Operation**: Independent execution mode

## Integration Patterns

### Host Model Interface

```fortran
! Standard integration pattern
use CATChemAPI_Mod

type(CATChemType) :: catchem
integer :: rc

! Initialize
call catchem%init(config_file='catchem.yml', rc=rc)

! Main time loop
do time_step = 1, num_steps
  ! Update meteorological fields
  call catchem%update_met_fields(met_data, rc)

  ! Run chemistry and transport
  call catchem%run(dt=time_step_size, rc=rc)

  ! Extract chemical fields
  call catchem%get_chem_fields(chem_data, rc)
end do

! Clean up
call catchem%finalize(rc)
```

## Supported Frameworks

### CCPP (Common Community Physics Package)
- **[CCPP Integration](ccpp.md)** - CCPP-compliant interface
- Physics scheme integration
- Metadata specification
- Build system integration

### NUOPC (National Unified Operational Prediction Capability)
- **[NUOPC Integration](nuopc.md)** - ESMF/NUOPC component
- ESMF field and grid handling
- Multi-component coupling
- Timekeeping and scheduling

## Data Exchange

### Field Mapping

```yaml
# Field mapping configuration
field_mapping:
  # Meteorological fields
  temperature:
    host_name: "T"
    units: "K"
    description: "Air temperature"

  pressure:
    host_name: "P"
    units: "Pa"
    description: "Air pressure"

  # Chemical species
  ozone:
    host_name: "O3"
    units: "kg/kg"
    description: "Ozone mixing ratio"
```

### Grid Compatibility

CATChem supports various grid types:

- **Latitude-Longitude**: Regular and Gaussian grids
- **Cubed-Sphere**: FV3 native grid
- **Unstructured**: MPAS and other unstructured grids
- **Nested Grids**: Regional high-resolution domains

## Configuration

### Host Model Settings

```yaml
# Integration configuration
integration:
  framework: "NUOPC"           # CCPP, NUOPC, FV3, standalone
  grid_type: "cubed_sphere"    # lat_lon, cubed_sphere, unstructured

  # Data exchange
  coupling_frequency: 300      # Seconds
  field_interpolation: true    # Enable field interpolation

  # Performance
  shared_memory: true          # Use shared memory when possible
  parallel_io: true           # Parallel I/O operations
```

## Best Practices

### Performance Optimization

1. **Minimize Data Copying**: Use shared memory when possible
2. **Efficient Field Mapping**: Avoid unnecessary interpolation
3. **Balanced Coupling**: Match time scales appropriately
4. **Memory Management**: Clean up temporary arrays

### Error Handling

```fortran
! Robust error handling in host model
call catchem%run(dt, rc)
if (rc /= CC_SUCCESS) then
  call error_handler%log_error("CATChem execution failed", rc)
  call emergency_cleanup()
  return
end if
```

### Debugging Integration

```yaml
# Debug integration issues
debug:
  integration: true
  field_exchange: true
  timing: true

logging:
  level: "debug"
  integration_events: true
```

## Testing Integration

### Unit Tests

```fortran
! Test host model interface
program test_integration
  use CATChemAPI_Mod
  use TestFramework_Mod

  call test_initialization()
  call test_field_exchange()
  call test_time_stepping()
  call test_finalization()

end program
```

### System Tests

```bash
# Integration test suite
cd tests/integration
./run_ccpp_tests.sh
./run_nuopc_tests.sh
```

## Troubleshooting

### Common Issues

- **Grid Mismatch**: Verify grid compatibility and interpolation
- **Unit Conversion**: Check field units and conversions
- **Memory Leaks**: Monitor memory usage during coupling
- **Timing Issues**: Ensure consistent time stepping

### Diagnostic Tools

```yaml
# Integration diagnostics
diagnostics:
  field_exchange:
    enabled: true
    frequency: 300

  coupling_metrics:
    enabled: true
    fields: ["temperature", "pressure", "ozone"]
```

## See Also

- [Process Development](../processes/index.md)
- [Core Systems](../core/index.md)
- [Configuration Guide](../core/configuration.md)
