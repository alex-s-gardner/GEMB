# gemb documentation
`gemb` returns a center coordinates from GEMB column spacings. 

# Syntax 
```matlab
OutData = gemb(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, ClimateForcing, ModelParam)
OutData = gemb(..., verbose=true)
OutData = gemb(..., display_waitbar=false)
```

# Description
`OutData = gemb(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, ClimateForcing, ModelParam)` produces time series of snow, firn, and ice properties OutData from input vectors of the initial column T, dz, d, W, re, gdn, gsp, and a, a_diffuse. Input `ClimateForcing` is a structure containing time series of surface forcing parameters descibed below. Input `ModelParam` is from the function `model_initialize_parameters.m`. 

| Input | Units | Description
|-------|-------|-------------|
|`temperature`      | K     | Vector of initial layer temperatures |
|`dz`     | m     | Vector of initial layer thicknesses | 
|`density`      | kg m^-2 |        Vector of initial layer densities |
|`water`      | m     | Vector of initial layer water content |
|`grain_radius`     | m     | Vector of initial grain sizes (effective radius) |
|`grain_dendricity`    |       | Vector of initial grain dendricity (0-1) |
|`grain_sphericity`    |       | Vector of initial grain sphericity (0-1)|
|`albedo`      |       | Initial surface albedo (0-1) |
|`albedo_diffuse`|     | Initial diffuse albedo accumulator |
|**`ClimateForcing`** | | Structure containing time-series meteorological data|
| `ClimateForcing.dates`      | datenum   |  Time vector|
| `ClimateForcing.shortwave_downward`       | W m^-2     |  Downward shortwave radiation|
| `ClimateForcing.longwave_downward`       | W m^-2     |  Downward longwave radiation|
| `ClimateForcing.temperature_air`    | K          |  Air temperature|
| `ClimateForcing.pressure_air`    | Pa         | Air pressure|
| `ClimateForcing.relative_humidity`        | %          |  Relative humidity|
| `ClimateForcing.vapor_pressure`     | Pa         |  Vapor pressure|
| `ClimateForcing.wind_speed`         | m s^-1     |  Wind speed|
| `ClimateForcing.precipitation`         | kg m^-2    |  Precipitation |
| `ClimateForcing.wind_observation_height`         |            |     | 
| `ClimateForcing.temperature_observation_height` |            |     |  
| `ClimateForcing.temperature_air_mean` |            |     |  
| `ClimateForcing.wind_speed_mean`     |            |     |  
| `ClimateForcing.precipitation_mean`     |            |     |  
|**`ModelParam`** |     | Structure containing model configuration parameters|
| `ModelParam.run_prefix`      | datenum   |  Time vector|
| `ModelParam.spinpup_cycles`      | datenum   |  Time vector|
| `ModelParam.daten`      | datenum   |  Time vector|

`OutData = gemb(..., verbose=true)` turns on additional checks to ensure the model conserves mass and energy for every timestep and logs the results. Note: the verbose=true option may add ~10% to total processing time. 

`OutData = gemb(..., display_waitbar=false)` disables the graphical waitbar. 

# Example
Run a basic example. 

```matlab
% Initialize model parameters:
ModelParam = model_initialize_parameters();

ModelParam.output_frequency = "daily"; 

% Generate sample data: 
time_step_hours = 3;
ClimateForcing = simulate_climate_forcing("test_1", time_step_hours);

% Initialize grid:
[temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = model_initialize_column(ModelParam, ClimateForcing);

% Run GEMB: 
OutData = gemb(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, ClimateForcing, ModelParam);
```

# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 
