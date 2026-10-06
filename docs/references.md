# Scientific References

This document contains scientific references and citations relevant to CATChem development and usage.

## Core Scientific References

### Atmospheric Chemistry Modeling

1. **Stockwell, W. R., et al.** (1997). The second generation regional acid deposition model chemical mechanism for regional air quality modeling. *Journal of Geophysical Research*, 102(D22), 25847-25879.

2. **Seinfeld, J. H., & Pandis, S. N.** (2016). *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change* (3rd ed.). John Wiley & Sons.

3. **Zaveri, R. A., & Peters, L. K.** (1999). A new lumped structure photochemical mechanism for large‐scale applications. *Journal of Geophysical Research*, 104(D23), 30387-30415.

### Numerical Methods

4. **Verwer, J. G., & Simpson, D.** (1995). Explicit methods for stiff ODEs from atmospheric chemistry. *Applied Numerical Mathematics*, 18(4), 413-430.

5. **Jacobson, M. Z.** (2005). *Fundamentals of Atmospheric Modeling* (2nd ed.). Cambridge University Press.

### Earth System Modeling

6. **Grell, G. A., et al.** (2005). Fully coupled "online" chemistry within the WRF model. *Atmospheric Environment*, 39(37), 6957-6975.

7. **Lamarque, J. F., et al.** (2012). CAM-chem: description and evaluation of interactive atmospheric chemistry in the Community Earth System Model. *Geoscientific Model Development*, 5(2), 369-411.

## Process-Specific References

### Aerosol Processes

8. **Ackermann, I. J., et al.** (1998). Modal aerosol dynamics model for Europe: development and first applications. *Atmospheric Environment*, 32(17), 2981-2999.

9. **Binkowski, F. S., & Shankar, U.** (1995). The regional particulate matter model: 1. Model description and preliminary results. *Journal of Geophysical Research*, 100(D12), 26191-26209.

### Gas-Phase Chemistry

10. **Carter, W. P.** (2010). Development of the SAPRC-07 chemical mechanism. *Atmospheric Environment*, 44(40), 5324-5335.

11. **Yarwood, G., et al.** (2005). Updates to the Carbon Bond chemical mechanism: CB05. Final Report prepared for US EPA.

### Dry Deposition

12. **Zhang, L., et al.** (2003). A size-segregated particle dry deposition scheme for an atmospheric aerosol module. *Atmospheric Environment*, 37(4), 549-560.

13. **Wesely, M. L.** (1989). Parameterization of surface resistances to gaseous dry deposition in regional-scale numerical models. *Atmospheric Environment*, 23(6), 1293-1304.

### Emissions Processing

14. **Guenther, A., et al.** (2006). Estimates of global terrestrial isoprene emissions using MEGAN (Model of Emissions of Gases and Aerosols from Nature). *Atmospheric Chemistry and Physics*, 6(11), 3181-3210.

15. **Wiedinmyer, C., et al.** (2011). The Fire INventory from NCAR (FINN): a high resolution global model to estimate the emissions from open burning. *Geoscientific Model Development*, 4(3), 625-641.

## Computational Methods

### High-Performance Computing

16. **Balaji, V., et al.** (2017). ESMF v6.3.0rp1 User's Guide. Earth System Modeling Framework.

17. **Theurich, G., et al.** (2016). The Earth System Modeling Framework. In *Encyclopedia of Computational Science and Engineering* (pp. 1-6). Springer.

### Parallel Computing

18. **Gropp, W., et al.** (1999). *Using MPI: portable parallel programming with the message-passing interface* (Vol. 1). MIT press.

19. **OpenMP Architecture Review Board** (2018). OpenMP Application Programming Interface Version 5.0.

## Model Validation Studies

### Regional Applications

20. **Dennis, R., et al.** (2010). A framework for evaluating regional‐scale numerical photochemical modeling systems. *Environmental Fluid Mechanics*, 10(4), 471-489.

21. **Emery, C., et al.** (2017). Recommendations on statistics and benchmarks to assess photochemical model performance. *Journal of the Air & Waste Management Association*, 67(5), 582-598.

### Global Applications

22. **Young, P. J., et al.** (2018). Tropospheric Ozone Assessment Report: Assessment of global-scale model performance for global and regional ozone distributions, variability, and trends. *Elementa*, 6(1), 10.

23. **Stevenson, D. S., et al.** (2006). Multimodel ensemble simulations of present‐day and near‐future tropospheric ozone. *Journal of Geophysical Research*, 111(D8).

## Software Engineering

### Code Architecture

24. **Gamma, E., et al.** (1995). *Design patterns: elements of reusable object-oriented software*. Addison-Wesley.

25. **Martin, R. C.** (2017). *Clean Architecture: A Craftsman's Guide to Software Structure and Design*. Prentice Hall.

### Scientific Computing

26. **Wilson, G., et al.** (2014). Best practices for scientific computing. *PLoS biology*, 12(1), e1001745.

27. **Sandve, G. K., et al.** (2013). Ten simple rules for reproducible computational research. *PLoS computational biology*, 9(10), e1003285.

## CATChem Process Scheme References

References for the schemes implemented in the [CATChem processes](processes/index.md).

### Dust

28. **Ginoux, P., et al.** (2001). Sources and distributions of dust aerosols simulated with the GOCART model. *Journal of Geophysical Research*, 106(D17), 20255-20273. <https://doi.org/10.1029/2000JD000053>
29. **Zhang, L., et al.** (2022). Development and evaluation of the Aerosol Forecast Member in the National Center for Environment Prediction (NCEP)'s Global Ensemble Forecast System (GEFS-Aerosols v1). *Geoscientific Model Development*, 15, 5337-5369. <https://doi.org/10.5194/gmd-15-5337-2022>

### Sea Salt

30. **Gong, S. L., Barrie, L. A., & Blanchet, J.-P.** (1997). Modeling sea-salt aerosols in the atmosphere: 1. Model development. *Journal of Geophysical Research*, 102(D3), 3805-3818. <https://doi.org/10.1029/96JD02953>
31. **Gong, S. L.** (2003). A parameterization of sea-salt aerosol source function for sub- and super-micron particles. *Global Biogeochemical Cycles*, 17(4), 1097. <https://doi.org/10.1029/2003GB002079>
32. **Jaeglé, L., et al.** (2011). Global distribution of sea salt aerosols: new constraints from in situ and remote sensing observations. *Atmospheric Chemistry and Physics*, 11(7), 3137-3157. <https://doi.org/10.5194/acp-11-3137-2011>

### Dry Deposition

33. **Wesely, M. L.** (1989). Parameterization of surface resistances to gaseous dry deposition in regional-scale numerical models. *Atmospheric Environment*, 23(6), 1293-1304. <https://doi.org/10.1016/0004-6981(89)90153-4>
34. **Zhang, L., Gong, S., Padro, J., & Barrie, L.** (2001). A size-segregated particle dry deposition scheme for an atmospheric aerosol module. *Atmospheric Environment*, 35(3), 549-560. <https://doi.org/10.1016/S1352-2310(00)00326-5>
35. **Emerson, E. W., et al.** (2020). Revisiting particle dry deposition and its role in radiative effect estimates. *Proceedings of the National Academy of Sciences*, 117(42), 26076-26082. <https://doi.org/10.1073/pnas.2014761117>

### Wet Deposition

36. **Liu, H., Jacob, D. J., Bey, I., & Yantosca, R. M.** (2001). Constraints from 210Pb and 7Be on wet deposition and transport in a global three-dimensional chemical tracer model driven by assimilated meteorological fields. *Journal of Geophysical Research*, 106(D11), 12109-12128. <https://doi.org/10.1029/2000JD900839>

### Chemistry and Settling (GOCART-2G)

37. **Collow, A. B., et al.** (2024). Benchmarking GOCART-2G in the Goddard Earth Observing System (GEOS). *Geoscientific Model Development*, 17, 1443-1468. <https://doi.org/10.5194/gmd-17-1443-2024>

## UFS-Chem References

The two papers documenting UFS-Chem version 1.0 (see [UFS-Chem](ufschem/index.md)).

38. **He, J., Zhang, L., Schwantes, R. H., Baker, B., et al.** (2026). Incorporating gas-phase chemistry into the Unified Forecast System (UFS) for global air quality applications. *Journal of Advances in Modeling Earth Systems*, 18(3). <https://doi.org/10.1029/2025MS005299>

39. **Zhang, L., Li, H., Grell, G. A., Bhattacharjee, P. S., et al.** (2026). Development of the CCPP-based GEFS-aerosols component in the Unified Forecast System for subseasonal prediction (UFS-Chem v1.0). *Geoscientific Model Development*, 19, 8597–8626. <https://doi.org/10.5194/gmd-19-8597-2026>

## Contributing References

If you have scientific papers, reports, or other references that should be included in this list, please:

1. **Follow the format**: Use consistent citation style (APA format preferred)
2. **Provide context**: Include a brief description of relevance to CATChem
3. **Check accessibility**: Ensure references are publicly accessible when possible
4. **Submit via pull request**: Add new references through the standard contribution process

### Reference Categories

- **Core Science**: Fundamental atmospheric chemistry and physics
- **Numerical Methods**: Mathematical and computational techniques
- **Validation**: Model evaluation and comparison studies
- **Applications**: Specific use cases and scientific applications
- **Software Engineering**: Code design and development practices

## Related Documentation

- [Contributing Guide](developer-guide/contributing.md)
- [Validation Studies](evaluation/index.md)
- [Process Documentation](processes/index.md)
- [Developer Guide](developer-guide/index.md)

## Citation Guidelines

When citing CATChem in scientific publications:

```
CATChem Development Team (2024). CATChem: Configurable ATmospheric CHEmistry model.
Available at: https://github.com/ufs-community/CATChem
```

For specific process implementations, please also cite the relevant scientific references listed above.

---

*This reference list is maintained by the CATChem development team and is updated regularly. Last updated: 2024*
