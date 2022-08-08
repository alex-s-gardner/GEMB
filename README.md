[![DOI](https://zenodo.org/badge/433556103.svg)](https://zenodo.org/badge/latestdoi/433556103)


# GEMB
A 1D column that simulates snow/firn/ice processes and surface-atmosphere mass and energy exchanges

The Glacier Energy and Mass Balance model (GEMB, the “B” is silent) is a vertical 1-D column model, i.e. no horizontal communication between model nodes, that simulates atmosphere-surface mass and energy exchanges and the internal evolution of snow, firn and ice. The model shares many characteristics with earlier published firn models that also simulate atmosphere-surface exchanges (e.g. Bassford, 2002; Bougamont & Bamber, 2005; Greuell & Konzelmann, 1994). The model is a finite-difference model with tens to hundreds of layers, the thickness of which are managed dynamically. It is forced at its surface with near-surface (2-10 m) estimates of precipitation, air temperature, wind speed, vapor pressure, surface pressure, and downwelling longwave and shortwave radiation fluxes and optional inputs of solar zenith angle, cloud optical thickness and bare ice albedo. At its bottom boundary, the model applies a constant thermal flux. Internally, the model simulates thermal diffusion, shortwave sub-surface penetration, meltwater retention, percolation and refreeze, effective snow grain size, dentricity, and sphericity, and compaction. In this section we detail specific implementation of various processes, and their options, within the model.

