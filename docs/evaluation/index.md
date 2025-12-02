# Evaluation Studies

This section contains comprehensive evaluation studies and benchmarking results for CATChem model components and configurations.

## Overview

CATChem evaluation encompasses multiple aspects:

- **Scientific Accuracy**: Comparison with observations and reference models
- **Performance Benchmarks**: Computational efficiency assessments
- **Code Quality**: Software engineering metrics and testing coverage
- **Reproducibility**: Verification of consistent results across platforms

More information coming soon!

## Evaluation Plans Against Observations

In order to identify and correct model biases, models must be frequently evaluated against a variety of observations. Evaluation is key to ensure the success of UFS-Chem and the CATChem library and modeling component. This necessitates further investment of model evaluation tools like MELODIES MONET that are fully incorporated into the workflow and go beyond traditional meteorological variables and surface measurements to easily and consistently compare model simulations against a variety of surface, aircraft, and satellite observations important for air quality and atmospheric composition. MELODIES MONET combines the Model EvaLuation using Observations, DIagnostics and Experiments Software (MELODIES) project at NSF NCAR with the Model and ObservatioN Evaluation Toolkit (MONET) project at NOAA to develop a Python diagnostic package that is open source, generic, portable, and model-agnostic.

MELODIES MONET already evaluates a variety of research and operational models against surface, aircraft, and satellite observations. MELODIES MONET source code is publicly available on [**GitHub**](https://github.com/NOAA-CSL/MELODIES-MONET) and the user guide is on [**ReadTheDocs**](https://melodies-monet.readthedocs.io/). MELODIES MONET builds from existing tools like [**MONETIO**](https://monetio.readthedocs.io), which reads in various model and observational datasets and converts them to standard formats and [**MONET**](https://monet-arl.readthedocs.io), which has general functions for pairing, analysis, and plotting.

MELODIES MONET has already been used to evaluate current and soon to be operational NOAA air quality forecasts against surface observations leading to a better understanding of model biases and identification of possible solutions.

As part of this effort we will evaluate UFS-Chem against NOAA collected field campaign data, which includes the following field campaigns to better understand urban air quality:

- [**SUNVEx**](https://csl.noaa.gov/projects/sunvex/) - Southwest Urban NOx and VOC Experiment
- [**AEROMMA**](https://csl.noaa.gov/projects/aeromma/) - Atmospheric Emissions and Reactions Observed from Megacities to Marine Areas
- [**CUPiDS**](https://csl.noaa.gov/projects/aeromma/cupids/) - Coastal Urban Plume Dynamics Study
- [**AGES+**](https://csl.noaa.gov/projects/ages/), which is the conglomerate of all field campaigns that occurred in the summer of 2023 in various cities across North America to study urban air quality from ground, mobile, and aircraft platforms
- [**AiRMAPS**](https://csl.noaa.gov/projects/airmaps/) - Airborne and Remote sensing Multi Air Pollutant Surveys

And the following field campaigns to better understand emissions, processing, and dynamics of fires:

- [**FIREX-AQ**](https://csl.noaa.gov/projects/firex-aq/) - Fire Influence on Regional to Global Environments and Air Quality
- [**CALFIDE**](https://csl.noaa.gov/groups/csl7/measurements/2022calfide/) - California Fire Dynamics Experiment

## Evaluation Schedule

### Regular Activities

| Activity | Frequency | Responsibility |
|----------|-----------|----------------|
| Regression and unit tests | Every commit | Automated CI |
| Performance benchmarks | Every new version | Dev team |
| Scientific evaluation | Every new version | Science team |

### Release Evaluation

#### Pre-release Checklist
- [ ] All regression and unit tests pass
- [ ] Performance benchmarks meet criteria
- [ ] Scientific evaluation complete
- [ ] Documentation updated
- [ ] External review completed

## Related Documentation

- [Scientific References](../references.md)
- [Testing Guide](../developer-guide/processes/testing.md)
- [Performance Guide](../developer-guide/performance.md)
- [Contributing Guide](../developer-guide/contributing.md)

## Support and Issues

### Reporting Problems
If evaluation tests fail or produce unexpected results:

1. **Check known issues**: Review current evaluation status
2. **Reproduce locally**: Verify issue on your system
3. **Document findings**: Include configuration, data, results
4. **Submit issue**: Submit [**GitHub Issue**](https://github.com/UFS-Community/CATChem/issues) with "evaluation" label

### Getting Help
- **Discussion forum**: Ask questions about evaluation methods on [**GitHub Discussions**](https://github.com/UFS-Community/CATChem/discussions)

---

*This evaluation framework is continuously evolving. Contributions and suggestions are welcome through the standard development process.*
