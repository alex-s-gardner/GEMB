# GEMB: Glacier Energy and Mass Balance Model

[![DOI](https://zenodo.org/badge/433556103.svg)](https://zenodo.org/badge/latestdoi/433556103)

## Overview

The Glacier Energy and Mass Balance (GEMB, the "B" is silent) model is a comprehensive one-dimensional physical model designed to simulate the surface energy balance and vertical firn evolution of glaciers and ice sheets. It couples atmospheric forcing with subsurface thermodynamics and densification physics to resolve the evolution of temperature, density, water content, and grain properties over time.

GEMB is a column model (no horizontal communication) of intermediate complexity, prioritizing computational efficiency to accommodate the multi-millennial spin-ups required for initializing deep firn columns. It is used for interpreting satellite altimetry data, firn studies, surface mass balance inversion from satellite data, ice core studies, uncertainty quantification and model exploration in cryosphere research. A complete description of GEMB can be found in [*Gardner et al*., 2023](https://doi.org/10.5194/gmd-16-2277-2023).

## Key Capabilities

GEMB simulates a wide range of physical processes critical to glacier health:

* **Surface Energy Balance (SEB):** Resolves radiative fluxes (shortwave/longwave) and turbulent heat fluxes (sensible/latent) using Monin-Obukhov similarity theory.
* **Subsurface Thermodynamics:** Solves the heat equation with phase change, meltwater percolation, and refreezing.
* **Firn Densification:** Simulates the compaction of snow into firn and ice using empirical or semi-empirical schemes.
* **Hydrology:** Tracks meltwater retention (irreducible water content), percolation, and refreezing using a "bucket" scheme.
* **Dynamic Albedo:** Models albedo evolution with long-term memory, accounting for grain growth and specific surface area.
* **Grid Management:** Utilizes a dynamic Lagrangian-style vertical grid that evolves with accumulation and ablation, automatically merging and splitting layers to maintain numerical stability.

## Prerequisites
 This version of GEMB requires MATLAB. Some examples in the documentation generate sample data using the `simulate_precipitation` function, which depends on the Statistics Toolbox, but otherwise GEMB can be run using a base MATLAB license without any additional toolboxes. If you do not have a MATLAB license, GEMB may be compatible with [Octave](https://octave.org/), but has only been tested with MATLAB. 

## Documentation

### GEMB Contents

[List of functions](docs/GEMB_overview.md) shows the contents of GEMB and describes function dependencies. 

[List of variables](docs/GEMB_variables.md) describes variable names and units used throughout GEMB. 

### Basic Workflow

Using GEMB requires four basic steps: 

1. Define a Climate Forcing structure that contains time series of surface forcing using these structure field names: 

    ```matlab
    CF.time                           = ; % datetime
    CF.temperature_air                = ; % K 
    CF.pressure_air                   = ; % Pa
    CF.precipitation                  = ; % kg m^-2
    CF.wind_speed                     = ; % m s^-1
    CF.shortwave_downward             = ; % W m^-2
    CF.longwave_downward              = ; % W m^-2
    CF.vapor_pressure                 = ; % Pa
    CF.temperature_air_mean           = ; % K
    CF.wind_speed_mean                = ; % m s^-1
    CF.precipitation_mean             = ; % kg m^-2
    CF.temperature_observation_height = ; % m
    CF.wind_observation_height        = ; % m
    ```
2. Use [`model_initialize_parameters`](docs/model_initialize_parameters_documentation.md) to define model parameters.
3. Use [`model_initialize_column`](docs/model_initialize_column_documentation.md) to create an initial profile of temperature, density, grid spacing, and other column parameters. 
3. Run the [`gemb`](docs/gemb_documentation.md) function.

### Tutorials

[Getting ERA5 reanalysis data](docs/ERA5_time_series_data.md) describes how to download ERA5 time series data for a single location. 

[Analyzing ERA5 reanalysis data](docs/ERA5_analysis.md) describes how to read ERA5 hourly data into MATLAB and analyse it with GEMB.

## Citation
If you use the GEMB software in your research, please cite the following paper:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023.

If you use GEMB model outputs, please cite:

Schlegel, N.-J., & Gardner, A. (2025). Output from the Glacier Energy and Mass Balance (GEMB v1.0) forced with 3-hourly ERA5 fields and gridded to 10km, Greenland and Antarctica 1979-2024 (1.4) [Data set]. Zenodo. [https://doi.org/10.5281/zenodo.14714746](https://doi.org/10.5281/zenodo.14714746)
