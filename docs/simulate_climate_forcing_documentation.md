# `simulate_climate_forcing` documentation
`simulate_climate_forcing` reproducibly generates synthetic climate forcing data for GEMB simulations based on predefined parameter sets.

# Syntax
```matlab
ClimateForcing = simulate_climate_forcing(set_id)
ClimateForcing = simulate_climate_forcing(set_id, time_step_hours)
```

# Description

**`ClimateForcing = simulate_climate_forcing(set_id)`** reproducibly simulates climate forcing based on a set of `set_id` parameters. The `set_id` parameters are defined in the `simulation_parameter_sets` function.
 
**`ClimateForcing = simulate_climate_forcing(set_id, time_step_hours)`** defines the temporal resolution of the forcing time series in hours. The default time step is defined by the `set_id. 

# Examples 
Create a climate forcing timetable at the default hourly resolution:

```matlab
ClimateForcing = simulate_climate_forcing("test_1");
```

Create a climate forcing timetable at three-hourly resoluion:

```matlab
ClimateForcing = simulate_climate_forcing("test_1", 3);
```
View the forcing variables: 

```matlab
>> head(ClimateForcing)
              time            temperature_air    pressure_air    precipitation    wind_speed    shortwave_downward    longwave_downward    vapor_pressure
      ____________________    _______________    ____________    _____________    __________    __________________    _________________    ______________
      01-Jan-1994 00:00:00        265.75            92638              0            2.2673            72.661               281.67              247.35    
      01-Jan-1994 02:59:52        266.33            92720              0            2.7777            152.38               265.65              261.84    
      01-Jan-1994 05:59:45        266.91            92802              0            3.5251            353.08               205.53               271.8    
      01-Jan-1994 08:59:37        267.49            92884              0            3.1225            556.57                219.4              296.43    
      01-Jan-1994 11:59:30        268.07            92966              0             2.507            641.43               217.26              312.81    
      01-Jan-1994 14:59:23        268.56            93043              0            2.1936            557.44               247.37              323.19    
      01-Jan-1994 17:59:15        269.04            93120              0            2.5887            354.31               225.75              333.43    
      01-Jan-1994 20:59:08        269.53            93197              0            2.1563            153.23               221.89               336.9    
```

View additional climate forcing properties: 

```matlab
>> ClimateForcing.Properties.CustomProperties
   ans = 
   CustomProperties with properties:
                 temperature_air_mean: 259.4000
                      wind_speed_mean: 5.1951
                   precipitation_mean: 1.1773e+03
       temperature_observation_height: 2
              wind_observation_height: 10
```

# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 
