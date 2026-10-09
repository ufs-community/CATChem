# Legacy/current diagnostic-name mapping

This page maps diagnostic names emitted by the legacy `upstream/develop`
process interfaces to names emitted by the current C++ diagnostic contracts
and NUOPC writer.

The legacy column describes the final per-species or per-bin output variable,
not only the internal packed/parent field returned by
`get_required_diagnostic_fields`. The current column describes the final
variable written to NetCDF by `write_process_diagnostics`.

## Summary

| Process | Legacy output names | Current output names | Compatibility |
|---|---|---|---|
| Dry deposition | `drydep_con_<species>`; `drydep_velocity_<species>` | `drydep_con_per_species_<species>`; `drydep_velocity_per_species_<species>` | Renamed because the current packed parent name is retained in the unpacked child name |
| Sea-salt emission | `seasalt_mass_emission_total`; `seasalt_number_emission_total`; `seasalt_mass_emission_<bin>`; `seasalt_number_emission_<bin>` | Totals unchanged; bins are `seasalt_mass_emission_bins_<bin>` and `seasalt_number_emission_bins_<bin>` | Per-bin names changed; totals unchanged |
| Settling | `settling_velocity_<species>`; `settling_flux_<species>` | `settling_velocity_per_species_per_level_<species>`; `settling_flux_per_species_<species>` | Renamed because the current packed parent name is retained in the unpacked child name |
| Sulfate chemistry | `Production_rate_<species>`; `PSO4_from_gaseous_SO2_per_level`; `PSO4_from_aqueous_SO2_per_level`; `DMS_emission_flux` | Same names | Direct match |
| Wet deposition | `wetdep_mass_<species>`; `wetdep_flux_<species>` | Same names | Direct match |
| Dust emission | No implementation in `upstream/develop` | See [Dust](#dust-emission) | Current-only |
| Carbon chemistry | No implementation in `upstream/develop` | See [Carbon chemistry](#carbon-chemistry) | Current-only |
| Photolysis | No corresponding legacy process diagnostic | `photolysis_rate_<reaction>` | Current-only |
| Gas chemistry | No registered process diagnostics | None registered; gas chemistry consumes photolysis-rate diagnostics | No output diagnostic |

`<species>` is the configured short name, normally lowercase (`so2`, `so4`,
`bc1`, etc.). `<bin>` is the configured short name for a size bin, such as
`seas1` or `dust3`.

## Dry deposition

Legacy `ProcessDryDepInterface_Mod.F90` declares the parent names
`drydep_con_per_species` and `drydep_velocity_per_species`, but registers and
writes one child per selected species:

| Legacy | Current | Units |
|---|---|---|
| `drydep_con_<species>` | `drydep_con_per_species_<species>` | `ug/kg` or `ppm` |
| `drydep_velocity_<species>` | `drydep_velocity_per_species_<species>` | `m/s` |

Current internal packed fields are `drydep_con_per_species` and
`drydep_velocity_per_species`, with a `Species` axis. The NUOPC writer unpacks
them into the current child names. See the current registration in
[catchem_process_drydep.cpp](../../../src/process/drydep/catchem_process_drydep.cpp:125).

## Sea-salt emission

| Legacy | Current | Units |
|---|---|---|
| `seasalt_mass_emission_total` | `seasalt_mass_emission_total` | `kg/m2/s` |
| `seasalt_number_emission_total` | `seasalt_number_emission_total` | Legacy registers `kg/m2/s`; current registers `#/m2/s` |
| `seasalt_mass_emission_<bin>` | `seasalt_mass_emission_bins_<bin>` | `kg/m2/s` |
| `seasalt_number_emission_<bin>` | `seasalt_number_emission_bins_<bin>` | Legacy registers `kg/m2/s`; current registers `#/m2/s` |

The current packed parent fields are `seasalt_mass_emission_bins` and
`seasalt_number_emission_bins`; those parent names are not written as final
NetCDF variables. The current totals and packed fields are registered in
[catchem_process_seasalt.cpp](../../../src/process/seasalt/catchem_process_seasalt.cpp:146).

The canonical default bins are `seas1`, `seas2`, `seas3`, `seas4`, and
`seas5`.

## Settling

| Legacy | Current | Units |
|---|---|---|
| `settling_velocity_<species>` | `settling_velocity_per_species_per_level_<species>` | `m/s` |
| `settling_flux_<species>` | `settling_flux_per_species_<species>` | `kg/m2/s` |

Current internal packed fields are `settling_velocity_per_species_per_level`
and `settling_flux_per_species`, with a `Species` axis. See
[catchem_process_settling.cpp](../../../src/process/settling/catchem_process_settling.cpp:228).

## Sulfate chemistry

These names are unchanged between legacy and current:

| Name | Units |
|---|---|
| `Production_rate_<species>` | `kg/kg/s` |
| `PSO4_from_gaseous_SO2_per_level` | `kg/kg/s` |
| `PSO4_from_aqueous_SO2_per_level` | `kg/kg/s` |
| `DMS_emission_flux` | `kg/m2/s` |

The current per-species `Production_rate_<species>` fields are registered
individually rather than exposed through a packed parent. See
[catchem_process_so4chem.cpp](../../../src/process/so4chem/catchem_process_so4chem.cpp:137).

## Wet deposition

| Legacy | Current | Units |
|---|---|---|
| `wetdep_mass_<species>` | `wetdep_mass_<species>` | `kg/m2` |
| `wetdep_flux_<species>` | `wetdep_flux_<species>` | `kg/m2/s` |

Legacy also declares the internal parent names
`wetdep_mass_per_species_per_level` and
`wetdep_flux_per_species_per_level`; current uses the same concepts internally
through the Fortran bridge but exposes only the per-species fields.

The current registration is in
[catchem_process_wetdep.cpp](../../../src/process/wetdep/catchem_process_wetdep.cpp:92),
and the scheme fills the mass and flux diagnostics in
[WetDepScheme_JACOB_Mod.F90](../../../src/process/wetdep/schemes/WetDepScheme_JACOB_Mod.F90:439).

With the default species metadata, the possible suffixes are:

```text
so2, h2o2, so4, msa, bc1, bc2, oc1, oc2,
dust1, dust2, dust3, dust4, dust5,
seas1, seas2, seas3, seas4, seas5
```

## Dust emission

Dust diagnostics are current-only relative to `upstream/develop`:

| Current name | Units | Notes |
|---|---|---|
| `dust_emission_total` | `kg/m2/s` | Total dust emission |
| `dust_emission_bin_<bin>` | `kg/m2/s` | Unpacked from current parent `dust_emission_bin` |
| `dust_horizontal_flux` | `kg/m/s` | Column diagnostic |
| `dust_moisture_correction` | unitless | Column diagnostic |
| `dust_effective_threshold` | `m/s` | Column diagnostic |
| `dust_utar_threshold_<bin>` | `m/s` | Unpacked from current parent `dust_utar_threshold` |

The canonical default bins are `dust1` through `dust5`. See
[catchem_process_dust.cpp](../../../src/process/dust/catchem_process_dust.cpp:216).

## Carbon chemistry

Carbon chemistry is current-only relative to `upstream/develop`:

| Current name | Units | Notes |
|---|---|---|
| `carbchem_prod_mass_<species>` | `kg/kg` | Unpacked from `carbchem_prod_mass` |
| `carbchem_loss_flux_<species>` | `kg/m2/s` | Unpacked from `carbchem_loss_flux` |
| `carbchem_phobic_mass_<species>` | `kg/kg` | Unpacked from `carbchem_phobic_mass` |
| `carbchem_phobic_flux_<species>` | `kg/m2/s` | Unpacked from `carbchem_phobic_flux` |

See [catchem_process_carbchem.cpp](../../../src/process/carbchem/catchem_process_carbchem.cpp:109).

## Photolysis and gas chemistry

The current photolysis process dynamically registers:

```text
photolysis_rate_<reaction>
```

with units `s-1`, where `<reaction>` comes from the TUV-x photolysis-rate
ordering. Gas chemistry does not register an output diagnostic; it reads
`photolysis_rate_<reaction>` fields when a mechanism rate parameter has the
`PHOTO.<reaction>` form.

## Diagnostic-selection differences

The output name mapping is not the only difference between the two paths:

- Legacy drydep and settling emit child names without the packed parent
  qualifier. Current output retains the packed parent qualifier.
- Legacy sea-salt emits `seasalt_*_emission_<bin>`; current output retains the
  `_bins` parent qualifier.
- Current packed fields are unpacked by the generic NUOPC writer into one
  variable per species or bin. The packed parent itself is not written.
- For wetdep, current `diag_species: []` means all species marked
  `is_wetdep`. Legacy uses `All` when the setting is missing, but an explicit
  empty list is not converted to `All`; use an explicit species list when
  comparing runs.
- A safe cross-path selector list for the default wetdep species is the list
  in the [Wet deposition](#wet-deposition) section, with the process prefix
  appropriate to the diagnostic being selected.

## Source of truth

The legacy names are defined by the `get_required_diagnostic_fields` and
registration/update routines in the legacy Fortran process interfaces on the
`upstream/develop` ref. Current names are defined by the C++ process
registrations and the axes-driven writer in
`drivers/nuopc/catchem_nuopc_interface.F90`.
