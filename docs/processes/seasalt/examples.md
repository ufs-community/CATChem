# SeaSalt Examples

This document provides usage examples for the SeaSalt process.

## Basic Configuration

```yaml
processes:
  - name: seasalt
    enabled: true
    scheme: Monahan
    species: [O3, NO2, SO2]
    diagnostics: [process_rate, tendency]
```

## Fortran Usage

```fortran
use SeaSaltProcess_Mod
type(SeaSaltProcessType) :: process
type(StateContainerType) :: container
integer :: rc

! Initialize
call process%init(container, rc)

! Run for timestep
call process%run(container, rc)

! Clean up
call process%finalize(rc)
```
