[![DOI](https://zenodo.org/badge/433556103.svg)](https://zenodo.org/badge/latestdoi/433556103)


# GEMB
A 1D column that simulates snow/firn/ice processes and surface-atmosphere mass and energy exchanges

The Glacier Energy and Mass Balance model (GEMB, the “B” is silent) is a vertical 1-D column model, i.e. no horizontal communication between model nodes, that simulates atmosphere-surface mass and energy exchanges and the internal evolution of snow, firn and ice. The model shares many characteristics with earlier published firn models that also simulate atmosphere-surface exchanges (e.g. Bassford, 2002; Bougamont & Bamber, 2005; Greuell & Konzelmann, 1994). The model is a finite-difference model with tens to hundreds of layers, the thickness of which are managed dynamically. It is forced at its surface with near-surface (2-10 m) estimates of precipitation, air temperature, wind speed, vapor pressure, surface pressure, and downwelling longwave and shortwave radiation fluxes and optional inputs of solar zenith angle, cloud optical thickness and bare ice albedo. At its bottom boundary, the model applies a constant thermal flux. Internally, the model simulates thermal diffusion, shortwave sub-surface penetration, meltwater retention, percolation and refreeze, effective snow grain size, dentricity, and sphericity, and compaction. In this section we detail specific implementation of various processes, and their options, within the model.

# Citation Information 
If you use the GEMB software, please cite the following: 

* Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, *Geosci. Model Dev.*, 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 

If you use GEMB model outputs, please use the citation below:

* Schlegel, N.-J., & Gardner, A. (2024). Output from the Glacier Energy and Mass Balance (GEMB v1.0) forced with 3-hourly ERA5 fields and gridded to 10km, Greenland and Antarctica 1979-2023 (1.3-5day) [Data set]. *Zenodo*. [https://doi.org/10.5281/zenodo.10806250](https://doi.org/10.5281/zenodo.10806250).