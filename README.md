# The airsea R package #

### What is this repository for? ###

* Tools for marine air-sea gas exchange studies.
* Version 0.1pre
* 
### What can the package do? ###
* Calculate water-side gas transfer velocity ($K_w$) using a range of modern parametisations.
* Calculate Schmit numbers for any gas
* Calculate net community production from in-situ oxygen observations.

### Installation ###

```{r install_airsea}
library(devtools)
install_github("tomhull/airsea")
```

### References ###

* Core gas transfer work taken from:
Johnson MT. A numerical scheme to calculate temperature and salinity dependent air-water transfer velocities for any gas. Ocean Sci Discuss. 2010;7: 251–290. doi:10.5194/osd-7-251-2010
