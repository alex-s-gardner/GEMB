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

1. **Define Climate Forcing:** Use [`model_initialize_forcing`](docs/model_initialize_forcing_documentation.md) to create a timetable of surface forcing from modeled or observed time series, or use [`simulate_climate_forcing`](simulate_climate_forcing_documentation.md) to create synthetic data for testing. 
2. **Define Model Parameters:** Use [`model_initialize_parameters`](docs/model_initialize_parameters_documentation.md) to set model parameters such as the number of spinup cycles or which densification model is used.
3. **Initialize a Column:** Use [`model_initialize_column`](docs/model_initialize_column_documentation.md) to create an initial profile of temperature, density, grid spacing, and other column properties. 
4. **Run GEMB:** Enter the Climate Forcing, Model Parameters, and initial Profile into the [`gemb`](docs/gemb_documentation.md) function.

### Tutorials

* [Getting ERA5 reanalysis data](docs/ERA5_time_series_data.md) describes how to download ERA5 time series data for a single location. 
* [Analyzing ERA5 reanalysis data](docs/ERA5_analysis.md) describes how to read ERA5 hourly data into MATLAB and analyse it with GEMB.

### Example scripts

The tutorials linked above and documentation for each function provide step-by-step instructions for how to use GEMB. If you prefer to run bare-bones scripts without verbose explanations, check out these examples: 

* `GEMB_example_synthetic.m` runs a simple example using synthetic climate forcing data.
* `GEMB_example_ERA5.m` examines ERA5 data near Summit Station, Greenland.

## Citation
If you use the GEMB software in your research, please cite the following paper:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023.

If you use GEMB model outputs, please cite:

Schlegel, N.-J., & Gardner, A. (2025). Output from the Glacier Energy and Mass Balance (GEMB v1.0) forced with 3-hourly ERA5 fields and gridded to 10km, Greenland and Antarctica 1979-2024 (1.4) [Data set]. Zenodo. [https://doi.org/10.5281/zenodo.14714746](https://doi.org/10.5281/zenodo.14714746)
