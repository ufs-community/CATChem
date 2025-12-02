# Scheme State Architecture: When to Use Persistent State

## Overview

The updated process templates now support **optional scheme state** - scheme state types are only generated when actually needed, reducing memory overhead and complexity.

## Decision Matrix: When Do You Need Scheme State?

### ❌ **NO Scheme State Needed** (Most Common)

**Use Case**: Simple schemes that compute emissions/reactions using only local variables

**Examples**:
- Basic emission schemes (GONG97, GEOS12)
- Simple chemical reaction schemes
- Lookup table-based calculations
- Schemes that operate independently on each column/timestep

**Generated Code**:
```fortran
! Only configuration type generated - NO state type
type :: SeaSaltSchemeGONG97Config
   real(fp) :: scale_factor = 1.0_fp
   logical :: weibull_flag = .false.
   ! ... parameters only
end type

! No SeaSaltSchemeGONG97State type generated
```

**YAML Configuration**:
```yaml
schemes:
  - name: "gong97"
    # No has_persistent_state, working_arrays, or persistent_variables
    parameters:
      scale_factor: 1.0
```

### ✅ **Scheme State NEEDED** (Special Cases)

#### **Case 1: Working Arrays for Complex Calculations**

**Use Case**: Schemes need temporary arrays that persist between calls or are expensive to allocate

**Example**: Size distribution calculations, multi-step algorithms

**YAML Configuration**:
```yaml
schemes:
  - name: "complex_scheme"
    has_persistent_state: true
    working_arrays:
      - name: "size_distribution"
        dimensions: ["n_size_bins", "n_columns"]
        description: "Size distribution working array"
      - name: "intermediate_results"
        dimensions: ["n_species", "n_levels", "n_columns"]
        description: "Multi-step calculation results"
```

**Generated Code**:
```fortran
type :: ProcessComplexSchemeState
   real(fp), allocatable :: size_distribution(:,:)
   real(fp), allocatable :: intermediate_results(:,:,:)
   ! Allocation/deallocation methods included
end type
```

#### **Case 2: Persistent Variables Between Timesteps**

**Use Case**: Schemes need to remember state from previous computations

**Example**: Time-dependent processes, iterative solvers, accumulated diagnostics

**YAML Configuration**:
```yaml
schemes:
  - name: "iterative_scheme"
    has_persistent_state: true
    persistent_variables:
      - name: "last_computation_time"
        type: "real(fp)"
        default: "-1.0_fp"
      - name: "convergence_history"
        type: "real(fp), dimension(10)"
        default: "0.0_fp"
```

#### **Case 3: Scheme-Specific Diagnostics**

**Use Case**: Schemes generate internal diagnostics beyond process-level ones

**YAML Configuration**:
```yaml
schemes:
  - name: "diagnostic_scheme"
    diagnostics:
      - name: "internal_rate_constant"
        description: "Scheme-internal rate diagnostic"
      - name: "convergence_metric"
        description: "Iterative convergence diagnostic"
```

## Performance Impact

### **Memory Savings** (No Scheme State)
```fortran
! Old approach - always generated:
type(SeaSaltSchemeStateGONG97) :: gong97_state    ! ~100KB per scheme
type(SeaSaltSchemeStateGONG03) :: gong03_state    ! ~100KB per scheme  
type(SeaSaltSchemeStateGEOS12) :: geos12_state    ! ~100KB per scheme
! Total: ~300KB per process instance

! New approach - only when needed:
! No state types for simple schemes
! Total: 0KB for typical emission processes
```

### **Code Simplicity** (No Scheme State)
```fortran
! No allocation/deallocation
! No state management  
! No memory errors
! Simpler debugging
```

## Migration Guide

### **For Existing Schemes**

1. **Evaluate if persistent state is actually needed**
   - Do you use arrays between timesteps? → YES, need state
   - Do you accumulate values over time? → YES, need state  
   - Do you only compute with inputs? → NO, no state needed

2. **Update YAML configuration**
   ```yaml
   # Add ONLY if needed:
   has_persistent_state: true
   working_arrays: [...]
   persistent_variables: [...]
   ```

3. **Regenerate process files** - templates will automatically:
   - Skip state generation for simple schemes
   - Generate proper state for complex schemes
   - Handle allocation/deallocation correctly

### **For New Schemes**

1. **Start simple** - Don't add persistent state initially
2. **Add state only when needed** - Most schemes don't need it
3. **Use specific array names** - Better than generic `work_array_1`

## Examples

### **Simple Emission Scheme** (No State)
```fortran
subroutine compute_gong97(n_levels, n_species, config, &
                          frocean, u10m, v10m, sst, &
                          species_conc, species_tendencies)
   ! All variables are local - no persistent state needed
   real(fp) :: local_emission_rate
   real(fp) :: wind_speed

   wind_speed = sqrt(u10m**2 + v10m**2)
   local_emission_rate = config%scale_factor * wind_function(wind_speed)
   species_tendencies(:) = local_emission_rate * size_distribution(:)
end subroutine
```

### **Complex Scheme** (With State)  
```fortran
subroutine compute_complex(n_levels, n_species, config, state, &
                           met_fields, species_conc, species_tendencies)
   ! Uses persistent working arrays and variables
   if (state%last_computation_time < 0.0) then
      call initialize_size_distribution(state%size_distribution)
   end if

   call update_distribution(state%size_distribution, met_fields)
   state%last_computation_time = current_time
end subroutine
```

## Recommendation

**Default to NO persistent state** - only add it when you have a specific need for:
- Working arrays that are expensive to allocate repeatedly
- Variables that must persist between timesteps  
- Scheme-specific diagnostics beyond process-level ones

This approach results in **cleaner, more efficient, and easier to maintain** code for the majority of atmospheric process schemes.
