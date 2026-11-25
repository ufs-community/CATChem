# Photolysis (TUV-x)

CATChem's `PhotolysisProcess` uses NCAR's TUV-x solver from MUSICA to calculate
photolysis rates in `s-1`. The process runs on the CPU, one atmospheric column
at a time. Build CATChem with `CATCHEM_BUILD_MUSICA=ON` and
`CATCHEM_ENABLE_KOKKOS=ON` to enable it.

## Configure photolysis and gas-phase chemistry

Schedule photolysis before gas chemistry so MICM reads the rates from the
current timestep. This example shows the schedule and file settings for a
coupled run. Replace both paths with your matching TUV-x and MICM configuration
files, and set `simulation.species_filename` to the CATChem species file.

```yaml
run_phases:
  chemistry:
    processes: [photolysis, gaschem]

processes:
  photolysis:
    activate: true
    config_file: "/path/to/tuvx/config.yml"
  gaschem:
    activate: true
    config_file: "/path/to/micm/config.yaml"
```

If you omit the photolysis `config_file`, CATChem looks for
`src/external/musica/configs/tuvx/tuv_5_4.yml` relative to the working directory.
Make the data files referenced by your TUV-x configuration available at their
configured paths. The shipped TUV 5.4 and Chapman files use different reaction
labels and require label alignment before they can exchange those rates.

## Feed gas species into TUV-x

In the CATChem species file, tag any gas species whose concentration should feed a
TUV-x profile. For example, this entry supplies nitrogen dioxide to a TUV-x
radiator that references the profile `NO2`.

```yaml
- name: no2
  __is_gas: true
  molecular weight [kg mol-1]: 0.046
  roles: [photolysis.NO2]
```

The suffix in `photolysis.<profile>` must match the TUV-x profile name, including
case. CATChem creates the profile if the TUV-x file omits it. For an explicit
species binding, any profile already defined under that name must be
host-updatable and use `molecule cm-3`. The TUV-x configuration still defines
how the profile is used, including its radiator and cross section where
applicable.

CATChem converts each tagged species from dry-air ppmv to molecules per cubic
centimeter using the local dry-air density. It multiplies that value by layer
thickness in centimeters to supply the layer column density in molecules per
square centimeter. Negative mixing ratios contribute zero absorption without
changing the host concentrations. A nonfinite tagged mixing ratio raises an error.

The legacy `photolysis.ozone` role maps to `O3` and preserves any fixed ozone
profile in an existing configuration. For untagged `O2` and `O3`, CATChem uses
configured profiles or background mixing ratios of 20.95% and 0.3 ppmv,
respectively.

## Supply temperatures

CATChem uses layer-center `T` in kelvin and layer thickness `BXHEIGHT` in meters.
Temperature-dependent quantum yields also use interface temperatures. A host
can supply `TEDGE` in kelvin using the axis-aware meteorology binding with the
interface axis and shape `(columns, levels + 1, 1)`. When `TEDGE` is absent or
stale, CATChem linearly interpolates `T` in height at interior interfaces and
uses the nearest layer-center temperature at the bottom and top. A fixed
temperature profile defined in the TUV-x configuration retains its values.
CATChem rejects nonfinite or nonpositive interface temperatures supplied to a
host-updatable profile.

## Cloud feedback

Set `processes.photolysis.cloud_feedback: true` to enable the CAM-SIMA
approximation. Supply current layer fields `CLDF` (cloud fraction, 0–1),
`CLDLIQ` (liquid water mixing ratio, kg/kg), and `DELP` (pressure thickness,
Pa), each shaped `(columns, levels, 1)`. When `DELP` is absent, CATChem
derives it from current `PEDGE` pressure interfaces. CATChem computes
`155 * CLDLIQ * DELP / g * sqrt(CLDF)` for each layer and applies that optical
depth at every TUV-x wavelength. The radiator's single-scattering albedo and
asymmetry remain zero, matching the current CAM-SIMA implementation. Fields
must be current, finite and physically valid. With feedback disabled, the
existing TUV-x configuration is used.

The current CCPP interface binds cloud liquid water and pressure thickness but
only has a column cloud fraction. It must be extended to provide layer `CLDF`
before enabling this option through that interface. The NUOPC mapping imports
`inst_cloud_frac_levels` as `CLDF` and `inst_pres_interface` as `PEDGE`,
which can provide `DELP`. It still needs a current cloud liquid-water mixing
ratio as `CLDLIQ`. Aerosol feedback remains future work.

## Read rates and verify coupling

CATChem computes the solar zenith angle from the column latitude, longitude,
and midpoint time. It calls `musica::RunTuvx` and averages the two interface
rates that bound each layer. It publishes the layer-center rates as
`photolysis_rate_<label>` diagnostics in `s-1`, with shape `(columns, levels)`.
If the solver reports an error for a column, CATChem logs it and clears that
column's rates.

MICM parameters named `PHOTO.<label>` read the diagnostic with the same
case-sensitive reaction label. A missing diagnostic gives a zero rate. Profile
names and reaction labels are independent. See the
[gas-phase chemistry guide](../gaschem/index.md) for mechanism settings.

In a build with MUSICA enabled, run both integration tests.

```bash
ctest --test-dir build/tests -R '^test_catchem_(photolysis|gaschem)$' --output-on-failure
```

`test_catchem_photolysis` compares a temperature-dependent nitrate quantum yield
with a constant-yield reference. It checks both interpolated and supplied
interface temperatures. `test_catchem_gaschem` checks daylight and nighttime
rates, attenuation by tagged absorbers, and the resulting change in MICM
concentrations for Chapman v0 and v1 mechanisms.
