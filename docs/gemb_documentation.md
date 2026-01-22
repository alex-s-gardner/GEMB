# gemb documentation
`gemb` returns a center coordinates from GEMB column spacings. 

# Syntax 
```matlab
OutData = gemb(T, dz, d, W, re, gdn, gsp, a, a_diffuse, ClimateForcing, ModelParam)
OutData = gemb(..., verbose=true)
OutData = gemb(..., display_waitbar=false)
```

# Description
`OutData = gemb(T, dz, d, W, re, gdn, gsp, a, a_diffuse, ClimateForcing, ModelParam)` produces time series of snow, firn, and ice properties OutData from input vectors of the initial column T, dz, d, W, re, gdn, gsp, and a, a_diffuse. Input `ClimateForcing` is a structure containing time series of surface forcing parameters descibed below. Input `ModelParam` is from the function `model_initialize_parameters.m`. 

| Input | Units | Description
|-------|-------|-------------|
|`T`      | K     | Vector of initial layer temperatures |
|`dz`     | m     | Vector of initial layer thicknesses | 
|`d`      | kg m^-2 |        Vector of initial layer densities |
|`W`      | m     | Vector of initial layer water content |
|`re`     | m     | Vector of initial grain sizes (effective radius) |
|`gdn`    |       | Vector of initial grain dendricity (0-1) |
|`gsp`    |       | Vector of initial grain sphericity (0-1)|
|`a`      |       | Initial surface albedo (0-1) |
|`a_diffuse`|     | Initial diffuse albedo accumulator |
|**`ClimateForcing`** | | Structure containing time-series meteorological data|
| `ClimateForcing.daten`      | datenum   |  Time vector|
| `ClimateForcing.dsw0`       | W m^-2     |  Downward shortwave radiation|
| `ClimateForcing.dlw0`       | W m^-2     |  Downward longwave radiation|
| `ClimateForcing.T_air0`    | K          |  Air temperature|
| `ClimateForcing.p_air0`    | Pa         | Air pressure|
| `ClimateForcing.rh0`        | %          |  Relative humidity|
| `ClimateForcing.e_air0`     | Pa         |  Vapor pressure|
| `ClimateForcing.V0`         | m s^-1     |  Wind speed|
| `ClimateForcing.P0`         | kg m^-2    |  Precipitation |
| `ClimateForcing.Vz`         |            |     | 
| `ClimateForcing.Tz`         |            |     |  
| `ClimateForcing.T_air_mean` |            |     |  
| `ClimateForcing.V_mean`     |            |     |  
| `ClimateForcing.P_mean`     |            |     |  
|**`ModelParam`** |     | Structure containing model configuration parameters|
| `ModelParam.run_prefix`      | datenum   |  Time vector|
| `ModelParam.n_spinpup_cycles`      | datenum   |  Time vector|
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
[T, dz, d, W, re, gdn, gsp, a, a_diffuse] = model_initialize_column(ModelParam, ClimateForcing);

% Run GEMB: 
OutData = gemb(T, dz, d, W, re, gdn, gsp, a, a_diffuse, ClimateForcing, ModelParam);
```


# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 
