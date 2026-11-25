# Input Files

CATChem requires various input files to perform atmospheric chemistry simulations. This guide describes the required and optional input files, their formats, and how to prepare them.

## Overview

CATChem input files can be organized into several categories:

- **Configuration files** - YAML files controlling model behavior
- **Meteorological data** - Weather and atmospheric state variables
- **Emission data** - Anthropogenic and biogenic emission sources
- **Chemical data** - Species properties and reaction mechanisms
- **Geographic data** - Land use, topography, and surface properties

## Configuration Files

### Main Configuration File

The main CATChem YAML file defines the simulation inputs, process settings,
and execution order. Put paths to the species and emission mapping files under
`simulation`. List scheduled process names under `run_phases` and activate
each one under its key in the `processes` map.

See the [configuration reference](configuration.md) for the full layout and
the [coupled TUV-x and MICM example](../processes/photolysis/index.md#configure-photolysis-and-gas-phase-chemistry)
for a chemistry schedule.

### Process-Specific Configuration

Each process can have detailed configuration in separate files:

**Chemistry configuration**

CATChem reads species properties from `simulation.species_filename`. Set
`processes.gaschem.config_file` for a MUSICA v1 mechanism or
`processes.gaschem.config_dir` for a v0 CAMP mechanism. TUV-x reads its own
file under `processes.photolysis.config_file` and calculates photolysis rates
during the timestep. See the [coupled configuration example](../processes/photolysis/index.md#configure-photolysis-and-gas-phase-chemistry).

**Emissions Configuration (`emissions_config.yml`):**
```yaml
emissions:
  anthropogenic:
    inventory: EDGAR_v6
    sectors:
      - power_plants
      - transportation
      - industrial
      - residential
    temporal_profiles:
      diurnal: "diurnal_profiles.nc"
      weekly: "weekly_profiles.nc"
      seasonal: "seasonal_profiles.nc"

  biogenic:
    model: MEGAN
    land_use_data: "land_use.nc"
    plant_functional_types: "pft_data.nc"
    emission_factors: "megan_ef.nc"

  point_sources:
    file: "point_sources.csv"
    format: csv
    plume_rise: briggs
```

## Meteorological Data

### Required Variables

CATChem requires the following meteorological variables:

**3D Variables:**

- Temperature (K)
- Pressure (Pa)
- Specific humidity (kg/kg)
- Zonal wind (m/s)
- Meridional wind (m/s)
- Vertical velocity (m/s)
- Turbulent kinetic energy (m²/s²)

**2D Variables:**

- Surface pressure (Pa)
- 2m temperature (K)
- 2m specific humidity (kg/kg)
- 10m wind speed (m/s)
- Precipitation rate (kg/m²/s)
- Boundary layer height (m)
- Surface heat flux (W/m²)

### File Format

Meteorological data should be in NetCDF format with CF conventions:

```bash
ncdump -h met_data_20240101_00.nc
```

```
netcdf met_data_20240101_00 {
dimensions:
    time = 1 ;
    lon = 100 ;
    lat = 80 ;
    lev = 50 ;

variables:
    double time(time) ;
        time:units = "hours since 2024-01-01 00:00:00" ;
        time:calendar = "gregorian" ;

    float lon(lon) ;
        lon:units = "degrees_east" ;
        lon:standard_name = "longitude" ;

    float lat(lat) ;
        lat:units = "degrees_north" ;
        lat:standard_name = "latitude" ;

    float lev(lev) ;
        lev:units = "Pa" ;
        lev:standard_name = "air_pressure" ;
        lev:positive = "down" ;

    float temperature(time, lev, lat, lon) ;
        temperature:units = "K" ;
        temperature:standard_name = "air_temperature" ;

    float pressure(time, lev, lat, lon) ;
        pressure:units = "Pa" ;
        pressure:standard_name = "air_pressure" ;
}
```

### Data Sources

**Common Meteorological Data Sources:**

- **NCEP GFS** - Global Forecast System
- **ECMWF ERA5** - Reanalysis data
- **NCAR WRF** - Weather Research and Forecasting model
- **NOAA NAM** - North American Mesoscale model

**Data Preparation:**
```bash
# Convert GRIB to NetCDF
cdo -f nc copy gfs_data.grib met_data.nc

# Regrid to target domain
cdo remapbil,target_grid.txt input_met.nc output_met.nc

# Extract variables
cdo select,name=t,u,v,q,ps input_met.nc extracted_met.nc
```

## Emission Data

### Anthropogenic Emissions

**File Format:**
```
netcdf emis_anthro_202401.nc {
dimensions:
    time = 744 ;  // hourly for January
    lon = 100 ;
    lat = 80 ;
    sector = 10 ;

variables:
    float NOx_emissions(time, sector, lat, lon) ;
        NOx_emissions:units = "mol/m2/s" ;
        NOx_emissions:long_name = "NOx emissions" ;

    float CO_emissions(time, sector, lat, lon) ;
        CO_emissions:units = "mol/m2/s" ;

    float SO2_emissions(time, sector, lat, lon) ;
        SO2_emissions:units = "mol/m2/s" ;

    char sector_names(sector, 50) ;
        sector_names:long_name = "Emission sector names" ;
}
```

**Recommended Minimal Species:**

- NOx (as NO2 equivalent)
- CO
- SO2
- NH3
- VOCs (individual species or lumped)
- PM2.5 and PM10
- Black carbon (BC)
- Organic carbon (OC)

### Biogenic Emissions

**MEGAN Format:**
```
netcdf emis_bio_20240101.nc {
dimensions:
    time = 24 ;  // hourly
    lon = 100 ;
    lat = 80 ;
    pft = 16 ;   // plant functional types

variables:
    float isoprene_ef(pft, lat, lon) ;
        isoprene_ef:units = "mg/m2/h" ;
        isoprene_ef:long_name = "Isoprene emission factor" ;

    float monoterpene_ef(pft, lat, lon) ;
        monoterpene_ef:units = "mg/m2/h" ;

    float lai(time, lat, lon) ;
        lai:units = "m2/m2" ;
        lai:long_name = "Leaf area index" ;

    float temperature_2m(time, lat, lon) ;
        temperature_2m:units = "K" ;
}
```

### Point Source Emissions

**CSV Format:**
```csv
source_id,longitude,latitude,stack_height,stack_diameter,exit_velocity,exit_temperature,NOx_rate,SO2_rate,CO_rate
PLANT001,-74.5,40.2,100.0,3.5,15.2,423.0,0.025,0.015,0.008
PLANT002,-74.3,40.1,80.0,2.8,12.1,398.0,0.018,0.012,0.006
```

**NetCDF Format:**
```
netcdf point_sources.nc {
dimensions:
    source = 1000 ;
    time = 8760 ;  // hourly for year

variables:
    double longitude(source) ;
    double latitude(source) ;
    float stack_height(source) ;
    float NOx_emissions(time, source) ;
        NOx_emissions:units = "kg/s" ;
}
```

## Chemical Data

### Species Properties

**Species Definition File (`species.yml`):**
```yaml
species:
  - name: O3
    molecular_weight: 48.0
    henry_constant: 1.0e-2  # M/atm
    diffusivity: 1.8e-5     # m2/s
    reactivity: moderate

  - name: NO2
    molecular_weight: 46.0
    henry_constant: 7.0e-3
    diffusivity: 1.6e-5
    reactivity: high

  - name: SO2
    molecular_weight: 64.1
    henry_constant: 1.4
    diffusivity: 1.2e-5
    reactivity: high
```

### Chemical Mechanisms

The CATChem species file describes host species and their process roles. The
MICM mechanism defines the gas-phase reactions. To supply a gas concentration
to TUV-x, give its CATChem species entry a `photolysis.<profile>` role that
matches a TUV-x profile name. The TUV-x file defines radiators, cross sections,
quantum yields, and reaction labels. See
[photolysis species bindings](../processes/photolysis/index.md#feed-gas-species-into-tuv-x)
and the [MICM mechanism settings](../processes/gaschem/index.md#configuration).

## Geographic Data

### Land Use Data

**File Format:**
```
netcdf land_use.nc {
dimensions:
    lon = 100 ;
    lat = 80 ;
    category = 24 ;  // USGS land use categories

variables:
    float land_use_fraction(category, lat, lon) ;
        land_use_fraction:units = "1" ;
        land_use_fraction:long_name = "Land use category fraction" ;

    int dominant_category(lat, lon) ;
        dominant_category:long_name = "Dominant land use category" ;

    char category_names(category, 100) ;
}
```

**Land Use Categories:**

1. Urban and Built-up
2. Dryland Cropland and Pasture
3. Irrigated Cropland and Pasture
4. Mixed Dryland/Irrigated Cropland
5. Cropland/Grassland Mosaic
6. Cropland/Woodland Mosaic
7. Grassland
8. Shrubland
9. Mixed Shrubland/Grassland
10. Savanna
11. Deciduous Broadleaf Forest
12. Deciduous Needleleaf Forest
13. Evergreen Broadleaf Forest
14. Evergreen Needleleaf Forest
15. Mixed Forest
16. Water Bodies
17. Herbaceous Wetland
18. Wooded Wetland
19. Barren or Sparsely Vegetated
20. Herbaceous Tundra
21. Wooded Tundra
22. Mixed Tundra
23. Bare Ground Tundra
24. Snow or Ice

### Topography Data

```
netcdf topography.nc {
dimensions:
    lon = 100 ;
    lat = 80 ;

variables:
    float elevation(lat, lon) ;
        elevation:units = "m" ;
        elevation:standard_name = "surface_altitude" ;

    float roughness_length(lat, lon) ;
        roughness_length:units = "m" ;
        roughness_length:long_name = "Surface roughness length" ;
}
```

## Data Preparation Tools

### Python Scripts

**Meteorological Data Processing:**
```python
import xarray as xr
import numpy as np

def prepare_met_data(input_file, output_file):
    """Convert meteorological data to CATChem format"""

    # Read input data
    ds = xr.open_dataset(input_file)

    # Rename variables to CATChem conventions
    ds = ds.rename({
        't': 'temperature',
        'u': 'u_wind',
        'v': 'v_wind',
        'q': 'specific_humidity'
    })

    # Add required attributes
    ds.temperature.attrs['units'] = 'K'
    ds.temperature.attrs['standard_name'] = 'air_temperature'

    # Save in CATChem format
    ds.to_netcdf(output_file, unlimited_dims=['time'])

# Usage
prepare_met_data('gfs_data.nc', 'catchem_met.nc')
```

**Emission Data Processing:**
```python
def prepare_emission_data(edgar_file, output_file):
    """Convert EDGAR emissions to CATChem format"""

    ds = xr.open_dataset(edgar_file)

    # Convert units from kg/m2/s to mol/m2/s
    mw_nox = 46.0  # g/mol
    ds['NOx_emissions'] = ds['NOx'] / mw_nox * 1000

    # Add temporal variation
    time_factors = create_temporal_profiles()
    ds = apply_temporal_profiles(ds, time_factors)

    ds.to_netcdf(output_file)
```

### Command Line Tools

**CDO (Climate Data Operators):**
```bash
# Regrid meteorological data
cdo remapbil,target_grid.nc input_met.nc output_met.nc

# Extract time period
cdo seldate,2024-01-01,2024-01-31 input.nc output.nc

# Calculate daily means
cdo daymean input.nc daily_output.nc

# Merge multiple files
cdo mergetime input_*.nc merged_output.nc
```

**NCO (NetCDF Operators):**
```bash
# Extract variables
ncks -v temperature,pressure,humidity input.nc output.nc

# Rename variables
ncrename -v t,temperature input.nc

# Add attributes
ncatted -a units,temperature,c,c,"K" input.nc

# Concatenate along time dimension
ncrcat input_*.nc output.nc
```

## Data Quality Control

### Validation Checks

**Meteorological Data:**
```python
def validate_met_data(filename):
    """Validate meteorological input data"""

    ds = xr.open_dataset(filename)

    # Check for required variables
    required_vars = ['temperature', 'pressure', 'humidity']
    missing_vars = [v for v in required_vars if v not in ds.variables]
    if missing_vars:
        raise ValueError(f"Missing variables: {missing_vars}")

    # Check for reasonable ranges
    if ds.temperature.min() < 200 or ds.temperature.max() > 350:
        raise ValueError("Temperature values outside reasonable range")

    # Check for missing data
    if ds.temperature.isnull().any():
        raise ValueError("Missing temperature data detected")

    print("Meteorological data validation passed")
```

**Emission Data:**
```python
def validate_emission_data(filename):
    """Validate emission input data"""

    ds = xr.open_dataset(filename)

    # Check for negative emissions
    for var in ds.data_vars:
        if (ds[var] < 0).any():
            print(f"Warning: Negative emissions found in {var}")

    # Check total emissions
    total_nox = ds.NOx_emissions.sum()
    print(f"Total NOx emissions: {total_nox.values:.2e} mol/s")
```

## File Organization

### Recommended Directory Structure

```
/data/catchem_run/
├── config/
│   ├── catchem_config.yml
│   ├── chemistry_config.yml
│   └── emissions_config.yml
├── meteorology/
│   ├── met_20240101_00.nc
│   ├── met_20240101_06.nc
│   └── ...
├── emissions/
│   ├── anthropogenic/
│   │   ├── emis_anthro_202401.nc
│   │   └── emis_anthro_202402.nc
│   ├── biogenic/
│   │   ├── emis_bio_20240101.nc
│   │   └── emis_bio_20240102.nc
│   └── point_sources.csv
├── geographic/
│   ├── land_use.nc
│   ├── topography.nc
│   └── roughness.nc
└── chemistry/
    ├── species.yml
    ├── micm/
    │   └── config.yaml
    └── tuvx/
        └── config.yml
```

### File Naming Conventions

**Time-based Files:**

- `%Y` - 4-digit year (2024)
- `%m` - 2-digit month (01-12)
- `%d` - 2-digit day (01-31)
- `%H` - 2-digit hour (00-23)

## Troubleshooting

### Common Issues

**File Not Found Errors:**
```bash
# Check file paths in configuration
grep -r "file_template\|directory" config/

# Verify file existence
ls -la /data/meteorology/met_20240101_00.nc

# Check file permissions
stat /data/meteorology/met_20240101_00.nc
```

**Format Errors:**
```bash
# Check NetCDF file structure
ncdump -h problematic_file.nc

# Validate CF conventions
cfchecker problematic_file.nc

# Check for corrupted files
ncdump -t problematic_file.nc > /dev/null
```

**Unit Conversion Issues:**
```python
# Common unit conversions
# Temperature: Celsius to Kelvin
temp_k = temp_c + 273.15

# Pressure: hPa to Pa
pressure_pa = pressure_hpa * 100

# Mixing ratio: kg/kg to mol/mol
mixing_ratio_mol = mixing_ratio_kg * mw_air / mw_species
```

## Best Practices

1. **Validate all input data** before running simulations
2. **Use consistent units** throughout all input files
3. **Follow CF conventions** for NetCDF metadata
4. **Document data sources** and processing steps
5. **Check for missing data** and handle appropriately
6. **Use appropriate temporal resolution** for your application
7. **Verify spatial alignment** between different datasets
8. **Test with small domains** before full simulations

## References

- [CF Conventions](http://cfconventions.org/)
- [NetCDF User Guide](https://docs.unidata.ucar.edu/netcdf-c/current/)
- [EDGAR Emission Database](https://edgar.jrc.ec.europa.eu/)
- [MEGAN Biogenic Emissions](https://bai.ess.uci.edu/megan/)
- [Configuration Guide](configuration.md)
