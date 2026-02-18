# Overview of GEMB

![](https://chadagreene.com/GEMB_figures/gemb_dependency_graph.jpg)

The figure above shows all of the functions that are called by the `gemb_run` example script. It shows that before running `gemb`, you must first: 

1. Use the `model_initialize_parameters` function to initialize a `ModelParam` structure that contains default values for physics modules, grid geometry, and output controls.
2. Use the `simulate_climate_forcing` function to initialize a `ClimateForcing` structure that contains time series of synthetic air temperature, downward shortwave and longwave radiation, air pressure, relative humidty, vapor pressure, wind speed, and precipitation. 
3. Use the `model_initialize_column` function to create vertical columns that contain initial values of temperature, density, water content, grain properties, albedo, and column spacing. 

With model parameters, climate forcing, and the initial state of the column defined, the `gemb` function calls `gemb_core` for each time step of the climate forcing. At each time step, `gemb_core` calls a series of functions that update the column grain size, albedo, shortwave radiation, temperature, accumulation, meltwater, and density. In this process, the `manage_layers` function adjusts the depth and number of vertical layers in the model to ensure that the thickness of any single layer does not exceed thresholds set for the minimum and maximum allowable layer thickness. 

# Function List 

### Primary functions:  

[**`dz2z`**](dz2z_documentation.md) returns a center coordinates from GEMB column spacings. 

[**`gemb`**](gemb_documentation.md) runs the GEMB model.

[**`model_initialize_column`**](model_initialize_column_documentation.md) initializes a GEMB column based on specified model and climate forcing parameters.

**[`model_initialize_parameters`](model_initialize_parameters_documentation.md)** initializes and validates the model configuration options, setting default values for physics modules, grid geometry, and output controls.

**[`surface_timeseries`](surface_timeseries_documentation.md)** returns the top row of finite data in a column timeseries of GEMB output data.

**[`dewpoint_to_vapor_pressure`](dewpoint_to_vapor_pressure_documentation.md)** converts dewpoint temperature to actual vapor pressure. 

### Climate simulation and fitting functions:

**`fit_air_temperature`** estimates simulation coefficients from observed data.

**`fit_longwave_irradiance_delta`**

**`fit_precipitation`** fits a seasonal Markov-Gamma model to hourly precipitation.

**`fit_seasonal_daily_noise`** fits yearly/daily sinusoids and noise stats.

**`vapor_pressure_to_relative_humidity`** calculates relative humidity from vapor pressure and temperature.

**`simulate_air_pressure`** simulates screen-level atmospheric pressure.

**`simulate_air_temperature`** simulates air temp using fitted coefficients.

**`simulate_climate_forcing`** generates synthetic climate forcing data for GEMB simulations based on predefined parameter sets.

**`simulate_coeffs_disp`** displays structure fields as executable MATLAB code.

**`simulate_longwave_irradiance`** estimates downward longwave radiation.

**`simulate_longwave_irradiance_delta`**

**`simulate_precipitation`** generates synthetic rainfall using fitted coefficients.

**`simulate_seasonal_daily_noise`** generates data from fitted coefficients

**`simulate_shortwave_irradiance`** simulates clear sky shortwave irradiance.

**`relative_humidity_to_vapor_pressure`** estimates actual vapor pressure from temperature and relative humidity.

**`simulation_parameter_sets`** retrieves predefined parameter sets for climate forcing simulations.

**`varname2longname`** returns a short description of given variable names.

### Internal functions: 

**`calculate_accumulation`** adds precipitation and deposition to the model column.

**`calculate_albedo`** calculates snow, firn and ice albedo as a function of effective grain radius, density and cloud amount, and decay time & wetness.

**`calculate_density`** computes the densification of snow/firn.

**`calculate_grain_size`** models the evolution of effective snow grain size, dendricity, and sphericity.

**`calculate_melt`** computes the quantity of meltwater due to snow temperature in excess of 0 degrees C, determines pore water content, and adjusts grid spacing.

**`calculate_shortwave_radiation`** distributes absorbed shortwave radiation vertically within snow/ice.

**`calculate_temperature`** computes new temperature profile accounting for energy absorption and thermal diffusion.

**`decyear2datenum`** converts decimal year to MATLAB serial date number. [Currently only called by `simulate_climate_forcing`.]

[**`densification_lookup_M01`**](densification_lookup_M01_documentation.md) returns calibrated coefficients of a densification model. [Currently only called by `calculate_density`.]

**`fast_divisors`** is a faster implementation of the `divisors` function and does NOT require the Symbolic Math Toolbox. 

**`gemb_core`**  performs a single time-step of the GEMB model. Calculates grain growth, albedo, radiative transfer, thermodynamics, accumulation, melt, layer management, and densification.

**`manage_layers`** adjusts the depth and number of vertical layers in the model to ensure that the thickness of any single layer does not exceed thresholds set for the minimum and maximum allowable layer thickness.

**`thermal_conductivity`** computes the thermal conductivity profile for snow, firn, and ice based on density and temperature. [Currently only called by `calculate_temperature`.]

**`turbulent_heat_flux`** computes sensible and latent heat fluxes using Monin-Obukhov similarity theory [Currently only called by `calculate_temperature`.]

# List of variables 
For a complete list of variables and descriptions, see [GEMB variables](GEMB_variables.md)

# Climate Forcing Data 

### Test Data
We keep some test data that we use for examples in the documentation [here](https://chadagreene.com/GEMB_data/).

### ERA5 time series
If you would like to know where we got our example data and learn how you can some of your own, check out the [ERA5 time series data page](ERA5_time_series_data). 

# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
[https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 