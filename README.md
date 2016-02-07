[![DOI](https://zenodo.org/badge/20322/tomhull/airsea.svg)](https://zenodo.org/badge/latestdoi/20322/tomhull/airsea)
# The airsea R package #

### What is this repository for? ###

* Tools for marine air-sea gas exchange studies.
* Version 0.1

### What can the package do? ###

* Calculate water-side and air-side gas transfer velocity using a range of modern parametisations.
* Calculate Schmit numbers for (almost) any gas
* Calculate bubble supersaturation effects (bubble induced equilibrium concentration)
* Calculate net community production from in-situ oxygen observations (also known as open water metabolism).
* Calculate various seawater properies, density, viscosity...

### airsea and marelac ###

marelac and airsea have some similar features with differing implementation:

* Schmit number calculation, however marelac assumes 35 salinity, airsea has more gasses.
* Gas saturation concentration, airsea only includes O2 based on Benson and Krause. marelac uses Weiss for several gases and includes pressure compensation.

#### Gasses in marelac `gas_schmit` but not in airsea ####
* Helium (He)
* Sulfur hexafluoride (SF6)
* Trichlorofluoromethane (CCl3F)
* Dichlorodifluoromethane (CCl2f2)

### TODO ###
* explore marelac vs liang solubility
* fully implement liang2013 bubble and kw

### Installation ###

```r
library(devtools)
install_github("tomhull/airsea")
```

### References ###

* Hull T, Greenwood N, Kaiser J, Johnson M. Uncertainty and sensitivity in optode-based shelf-sea net community production estimates. Biogeosciences Discuss. Copernicus GmbH; 2015;12: 15611–15654. doi:10.5194/bgd-12-15611-2015

* Core gas transfer work taken from:
Johnson MT. A numerical scheme to calculate temperature and salinity dependent air-water transfer velocities for any gas. Ocean Sci Discuss. 2010;7: 251–290. doi:10.5194/osd-7-251-2010
