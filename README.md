# GEMB: Glacier Energy and Mass Balance Model

[![DOI](https://zenodo.org/badge/433556103.svg)](https://zenodo.org/badge/latestdoi/433556103)

## Overview

The Glacier Energy and Mass Balance (GEMB, the "B" is silent) model is a comprehensive one-dimensional physical model designed to simulate the surface energy balance and vertical firn evolution of glaciers and ice sheets. It couples atmospheric forcing with subsurface thermodynamics and densification physics to resolve the evolution of temperature, density, water content, and grain properties over time.

GEMB is a column model (no horizontal communication) of intermediate complexity, prioritizing computational efficiency to accommodate the multi-millennial spin-ups required for initializing deep firn columns. It is widely used for uncertainty quantification and model exploration in cryosphere research. A complete description of GEMB can be found in [*Gardner et al*., 2023](https://doi.org/10.5194/gmd-16-2277-2023).

## Key Capabilities

GEMB simulates a wide range of physical processes critical to glacier health:

* **Surface Energy Balance (SEB):** Resolves radiative fluxes (shortwave/longwave) and turbulent heat fluxes (sensible/latent) using Monin-Obukhov similarity theory.
* **Subsurface Thermodynamics:** Solves the heat equation with phase change, meltwater percolation, and refreezing.
* **Firn Densification:** Simulates the compaction of snow into firn and ice using empirical or semi-empirical schemes.
* **Hydrology:** Tracks meltwater retention (irreducible water content), percolation, and refreezing using a "bucket" scheme.
* **Dynamic Albedo:** Models albedo evolution with long-term memory, accounting for grain growth and specific surface area.
* **Grid Management:** Utilizes a dynamic Lagrangian-style vertical grid that evolves with accumulation and ablation, automatically merging and splitting layers to maintain numerical stability.


## Model Configuration

The model is highly configurable via the `model_initialize_parameters` function. Key configuration categories include:

1.  **Densification Physics:**
    * Options include **"HerronLangway"** (empirical), **"Arthern"** (semi-empirical), and **"Ligtenberg"** (semi-empirical).
    * Calibration coefficients for the Ligtenberg model can be specified for Antarctica or Greenland (e.g., `"Gre_RACMO_GS_SW0"`).

2.  **Grid Geometry:**
    * Vertical discretization is controlled by `column_ztop` (depth of high-resolution surface zone) and `column_zy` (stretching factor for deep layers).
    * Users define minimum (`column_dzmin`) and maximum (`column_dzmax`) layer thicknesses to ensure stability.

3.  **Energy Balance & Optical Properties:**
    * **Albedo Schemes:** Options include `"GardnerSharp"`, `"GreuellKonzelmann"`, `"BruneLeFebre"`, or `"BougamontBamber"`.
    * **Solar Penetration:** The `sw_absorption_method` toggles between surface-only (0) or subsurface extinction (1).
    * **Thermal Conductivity:** Select between `"Sturm"` or `"Calonne"` parameterizations.
    * **Emissivity:** Configurable methods (0, 1, 2) and thresholds based on grain size.

4.  **Initialization:**
    * Fresh snow density can be set via models like `"350kgm2"`, `"Fausto"`, or `"Kaspers"`.
    * Spin-up cycles (`n_spinup_cycles`) allow the model to reach equilibrium before saving output.

## Climate Forcing Inputs

GEMB requires high-frequency meteorological forcing data. The `simulate_climate_forcing` function generates synthetic data containing the following variables:

* **Radiation:**
    * Downward Shortwave (`dsw0`): Calculated based on solar geometry and location.
    * Downward Longwave (`dlw0`): Derived from air temperature and vapor pressure, with adjustments for cloud cover.
* **Thermodynamics:**
    * Air Temperature (`T_air0`): Simulated time series.
    * Air Pressure (`p_air0`): Calculated from elevation and temperature.
    * Humidity: Relative humidity (`rh0`) and vapor pressure (`e_air0`).
* **Dynamics & Mass:**
    * Wind Speed (`V0`): Stochastic generation with seasonal noise.
    * Precipitation (`P0`): Accumulated mass input.
* **Metadata:**
    * Includes location (lat, lon, elev) and measurement heights (Vz, Tz).

## Getting Started

### Prerequisites
* MATLAB (The core physics are written in MATLAB/C++).

### Basic Usage
The following MATLAB snippet demonstrates the standard workflow to configure, initialize, and run a GEMB simulation.

```matlab
% 1. Specify Model Parameters
% Initialize the default parameter structure. 
ModelParam = model_initialize_parameters();

% 2. Load Climate Forcing
% Generate or load meteorological forcing data (Temperature, Wind, Precip, etc.)
% 'test_1' loads a pre-defined test case.
ClimateForcing = simulate_climate_forcing("test_1");

% 3. Initialize Grid
% Set up the initial vertical column state (Temperature, Density, Water content, etc.)
[T, dz, d, W, re, gdn, gsp, a, a_diffuse] = ...
    model_initialize_column(ModelParam, ClimateForcing);

% 4. Run GEMB
% Execute the model. 'verbose' enables detailed logging.
verbose = true;
OutData = gemb(T, dz, d, W, re, gdn, gsp, a, a_diffuse, ...
               ClimateForcing, ModelParam, verbose);
```
### Visualizing Results
Once the model run is complete, the results are stored in the `OutData` structure. You can visualize the evolution of the firn column using standard MATLAB plotting commands.

```matlab
% Example: Plotting Surface Mass Balance Components
figure;
plot(OutData.time, OutData.M, 'r', 'DisplayName', 'Melt');
hold on;
plot(OutData.time, OutData.R, 'b', 'DisplayName', 'Runoff');
legend;
title('Modeled Surface Mass Balance Components');
ylabel('Mass Flux (kg m^{-2})');
xlabel('Time');

% Example: Plotting Temperature Profile Evolution
% OutData.T contains temperature profiles [depth x time]
% Note: Using the initial depth profile for visualization purposes.
figure;
imagesc(OutData.time, cumsum(OutData.dz(:,1)), OutData.T); 
colorbar;
title('Firn Temperature Evolution');
ylabel('Depth (m)');
xlabel('Time');
```

## Citation
If you use the GEMB software in your research, please cite the following paper:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023.

If you use GEMB model outputs, please cite:

Schlegel, N.-J., & Gardner, A. (2025). Output from the Glacier Energy and Mass Balance (GEMB v1.0) forced with 3-hourly ERA5 fields and gridded to 10km, Greenland and Antarctica 1979-2024 (1.4) [Data set]. Zenodo. [https://doi.org/10.5281/zenodo.14714746](https://doi.org/10.5281/zenodo.14714746)