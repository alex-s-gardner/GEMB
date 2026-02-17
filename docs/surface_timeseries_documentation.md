# `surface_timeseries` documentation
`surface_timeseries` returns the top row of finite data in a matrix.

# Syntax 
```matlab
A_surface = surface_timeseries(A)
```

# Description

`A_surface = surface_timeseries(A)` for input matrix `A` of dimensions MxN returns a 1xN vector `A_surface` of the uppermost finite values of the matrix `A`.

# Example  
Solve a firn column time series using synthetic data, then plot the skin temperature. 

```matlab
% Generate sample data:
time_step_hours = 3;
ClimateForcing = simulate_climate_forcing("test_1", time_step_hours);

% Initialize model parameters:
ModelParam = model_initialize_parameters(output_frequency="daily");

% Initialize grid:
Profile = model_initialize_column(ModelParam, ClimateForcing);

% Run GEMB: 
OutData = gemb(Profile, ClimateForcing, ModelParam);

% Get a time series of skin temperature: 
temperature_skin = surface_timeseries(OutData.temperature); 

% Create a datetime array to make it easy: 
dates_datetime = datetime(OutData.dates,'convertfrom','datenum'); 

% Plot the timeseries: 
plot(dates_datetime, temperature_skin)
ylabel 'Skin temperature (K)'

exportgraphics(gcf,'surface_timeseries_documentation_01.png',Resolution=300)
```

![](https://chadagreene.com/GEMB_figures/surface_timeseries_documentation_01.png)

In the figure above, the skin temperature never exceeds 273.15 K because if any additional energy is added it will be used to melt the ice! 

# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
[https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 
