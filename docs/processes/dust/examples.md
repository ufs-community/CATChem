# Dust Examples

This document provides usage examples for the Dust process.

## Basic Configuration

```yaml
processes:
  - name: dust
    enabled: true
    scheme: Fengsha
    species: [O3, NO2, SO2]
    diagnostics: [process_rate, tendency]
```

## Fortran Usage

```fortran
use DustProcess_Mod
type(DustProcessType) :: process
type(StateContainerType) :: container
integer :: rc

! Initialize
call process%init(container, rc)

! Run for timestep
call process%run(container, rc)

! Clean up
call process%finalize(rc)
```
