# `forcing_climatology` documentation
`forcing_climatology` creates a climatological average forcing structure that can be used to spinup a gemb run.

See also [`model_initialize_forcing`](model_initialize_forcing_documentation.md) and [`simulate_climate_forcing`](simulate_climate_forcing_documentation.md).

# Syntax 
```matlab
ClimateForcingSpinup = forcing_climatology(ClimateForcing)
ClimateForcingSpinup = forcing_climatology(ClimateForcing, datetime_range)
```

# Description

`ClimateForcingSpinup = forcing_climatology(ClimateForcing)` returns the climatological mean forcing of ClimateForcing structure for GEMB spinup runs. 

`ClimateForcingSpinup = forcing_climatology(ClimateForcing, datetime_range)` specifies a range of dates to include in the climatology. The `datetime_range` must be 1x2 datetime in the form `[datetime_start datetime_end]`. 

# Example 
Simulate a multi-decade time series of three-hourly climate forcing and convert it to climatological average. Then plot the raw data and climatological average temperature as a function of day of year. 

```matlab
% Generate sample data: 
time_step_hours = 3;
ClimateForcing = simulate_climate_forcing("test_1", time_step_hours);

% Convert full forcing time series to climatology:
ClimateForcingSpinup = forcing_climatology(ClimateForcing); 

figure
scatter(day(ClimateForcing.time, "dayofyear"), ClimateForcing.temperature_air, 10, year(ClimateForcing.time), "filled")
cb = colorbar; 
ylabel(cb, "Year")
hold on
plot(day(ClimateForcingSpinup.time, "dayofyear"), ClimateForcingSpinup.temperature_air,'k',linewidth=2)
axis tight
xlabel("Day of year")
ylabel("Air temperature (K)")
legend("ClimateForcing","Climatology","location","southwest")
legend boxoff
```

![](https://chadagreene.com/GEMB_figures/forcing_climatology_documentation_01.jpg)



# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
[https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 
