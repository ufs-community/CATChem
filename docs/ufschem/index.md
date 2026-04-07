# UFS-Chem

![UFS-Chem_logo](../diagrams/UFS-Chem_logo.png){: width="35%", align="right"}
The Unified Forecast System with Chemistry (UFS-Chem) is an effort to unify atmospheric chemistry and composition in NOAA's UFS and is a large collaboration between NOAA OAR's Chemical Sciences Laboratory, Air Resources Laboratory, Global Systems Laboratory, Geophysical Fluid Dynamics Laboratory, Global Monitoring Laboratory and external collaborators including NSF NCAR's Atmospheric Chemistry Observations & Modeling Laboratory, U.S. EPA's Office of Research and Development, NASA's Global Modeling and Assimilation Office, and university collaborators at the University of Colorado Boulder Cooperative Institute for Research in Environmental Sciences (CIRES), George Mason University Cooperative Institute for Satellite Earth System Studies (CISESS), and University of Wisconsin.

NOAA's **[Unified Forecast System (UFS)](https://www.ufs.epic.noaa.gov/)** is a community-based Earth modeling system that plans to provide a framework to efficiently incorporate research advances on atmospheric composition and chemistry into NOAA's operational forecasts and application areas. Currently, chemistry-related code for different applications including weather, air quality, and smoke & dust forecasting is incorporated into the UFS through different methods. This non-unified framework is inefficient, difficult for developers to maintain, and not conducive for adding capabilities within the UFS for research applications. Through this work, we plan to unify atmospheric chemistry and composition within the UFS by creating CATChem, the Configurable ATmospheric Chemistry library and modeling component. CATChem will be flexible such that users can select the correct level of chemical complexity for their research or operational application. CATChem will include the following processes: passive tracers, chemical kinetics, aerosols, photolysis, wet deposition, dry deposition, connections to emissions, and connection to physics schemes. We will link CATChem to the UFS to create UFS-Chem. When possible, we will use tools already developed or being developed by the research community like the **[Model Independent Chemistry Module (MICM)](https://www2.acom.ucar.edu/modeling/model-independent-chemistry-module-micm)**, which is a component of the **[MUlti-Scale Infrastructure for Chemistry and Aerosols (MUSICA)](https://www2.acom.ucar.edu/sections/multi-scale-infrastructure-chemistry-modeling-musica)**, led by NSF NCAR and the **[Community Regional Atmospheric Chemistry Multiphase Mechanism (CRACMM)](https://www.epa.gov/cmaq/cracmm)**, led by U.S. EPA in collaboration with NOAA CSL and university partners.

![UFS-Chem Overview](../diagrams/UFS-Chem_concept_figure.png)
*Figure 1: Conceptual Framework of UFS-Chem:* ***[He et al., JAMES, 2026](https://doi.org/10.1029/2025MS005299)***

We will also add enhanced research capabilities into UFS-Chem, which will include:

  1. Options to use gas and aerosol chemical mechanisms of varying complexity
  2. Options for passive tracers, which will also allow benchmark verification of mass conservation across UFS-Chem
  3. Novel parameterizations for fire weather
  4. Ability to easily couple chemistry in-line (i.e., fully coupled in the physics routines) to different physics options
  5. Development of a more flexible emissions processing system
  6. Interfacing with state-of-science atmospheric composition data assimilation capabilities
  7. Further investment of model evaluation tools like **[MELODIES-MONET](https://melodies-monet.readthedocs.io/)** that enable model comparison against a variety of observations

UFS-Chem will increase efficiency in code development, reduce costs for code maintenance, reduce time and effort for research transitions to operations, and enhance collaborations with the research community. Continued engagement with the atmospheric chemistry research community is critical to ensure that research advances are efficiently and promptly included within the UFS, so that NOAA continues to provide state-of-the-art forecasts and monitoring of atmospheric composition to inform key societal challenges and optimal decision-making for current and future generations.

## UFS-Chem version 1.0

As a first step towards the innovative and ambitious goal described above, the gas-phase tropospheric and stratospheric chemistry from the Atmosphere Model version 4.1 (AM4.1) is incorporated into CATChem version 1.0 and linked to the UFS to create UFS-Chem version 1.0. This is the first UFS-Chem configuration for global air quality applications and is a major step in our goal of unifying atmospheric chemistry and composition within the UFS for research and operational applications. UFS-Chem version 1.0 is fully developed, evaluated against surface, aircraft and satellite evaluations, and documented in **[He et al., JAMES, 2026](https://doi.org/10.1029/2025MS005299)**.

## UFS-Chem version 2.0

Development of CATChem version 2.0 and UFS-Chem version 2.0 as our fully configurable system is now underway. Please check back frequently for updates as this is an area of active development.
