# `model_initialize_forcing` documentation
`model_initialize_forcing` creates a timetable of climate forcing data for GEMB.

# Syntax
```matlab
ClimateForcing = model_initialize_forcing(time, temperature_air, pressure_air, precipitation, wind_speed, shortwave_downward, longwave_downward, vapor_pressure)
ClimateForcing = model_initialize_forcing(..., temperature_air_mean=value)
ClimateForcing = model_initialize_forcing(..., wind_speed_mean=value)
ClimateForcing = model_initialize_forcing(..., precipitation_mean=value)
ClimateForcing = model_initialize_forcing(..., temperature_observation_height=value)
ClimateForcing = model_initialize_forcing(..., wind_observation_height=value)
```

# Description

**`ClimateForcing = model_initialize_forcing(time, temperature_air, pressure_air, precipitation, wind_speed, shortwave_downward, longwave_downward, vapor_pressure)`** builds a timetable from an Nx1 time vector in datetime format and corresponding surface forcing vectors that must be of equal length. Inputs must be: 

| input | units | dimensions |
|-----|-----|-----|
|`time`             |  datenum  |(Nx1)|
|`temperature_air`  |  K    |    (Nx1)|
|`pressure_air`   |    Pa     |  (Nx1)|
|`precipitation`    |  kg m<sup>-2</sup>  |(Nx1) |
|`wind_speed`       |  m s<sup>-1</sup>   |(Nx1)|
|`shortwave_downward` |W m<sup>-2</sup>  | (Nx1)|
|`longwave_downward`  |W m<sup>-2</sup>  | (Nx1)|
|`vapor_pressure`   |  Pa   |   (Nx1)|
 
**`ClimateForcing = model_initialize_forcing(..., temperature_air_mean=value)`** specifies a climatological mean air temperature value. If a value is not specified, a warning message will indicate that the mean of `temperature_air` vector is assumed. 

**`ClimateForcing = model_initialize_forcing(..., wind_speed_mean=value)`** specifies a climatological mean wind speed value. If a value is not specified, a warning message will indicate that the mean of `wind_speed` vector is assumed. 

**`ClimateForcing = model_initialize_forcing(..., precipitation_mean=value)`** specifies a climatological mean precipitation value. If a value is not specified, a warning message will indicate that the mean of `precipitation` vector is assumed. 

**`ClimateForcing = model_initialize_forcing(..., temperature_observation_height=value)`** specifies the height of the observed or modeled temperature_air above the surface, in meters. If a value is not specified, a warning message will indicate that `2` m is assumed. 

**`ClimateForcing = model_initialize_forcing(..., wind_observation_height=value)`** specifies the height of the observed or modeled wind_speed above the surface, in meters. If a value is not specified, a warning message will indicate that `10` m is assumed. 

# Example
This example uses the ERA5 hourly land surface time series data described [here](ERA5_time_series_data.md). 

```matlab
% Define climate data filenames: 
filename_2m_temperature  = 'reanalysis-era5-land-timeseries-sfc-2m-temperature_summit.nc';
filename_pressure_precip = 'reanalysis-era5-land-timeseries-sfc-pressure-precipitation_summit.nc';
filename_radiation       = 'reanalysis-era5-land-timeseries-sfc-radiation-heat_summit.nc';
filename_wind            = 'reanalysis-era5-land-timeseries-sfc-wind_summit.nc';

% Reference date plus hours since the reference date: 
time_vector = datetime(1970,1,1) + duration(ncread(filename_2m_temperature,'valid_time'),0,0);

% Read temperature and pressure: 
temperature_air = ncread(filename_2m_temperature,'t2m'); 
pressure_air    = ncread(filename_pressure_precip,'sp');

% Multiply hourly precipitation by 1000 to convert from m to kg/m^2:  
precipitation   = ncread(filename_pressure_precip,'tp') * 1000;
```

At this point we run into a minor issue, because some precipitation values are less than zero: 

```matlab
>> min(precipitation)
ans =
  -4.4035e-05
```

The small negative precipitation does not represent evaporation. Instead, it is simply numerical noise we can set to zero and move on:

```matlab
% Fix numerical noise: 
precipitation(precipitation<0) = 0; 

% Wind speed is the hypotenuse of the vector components: 
wind_speed = hypot(ncread(filename_wind,'u10'),ncread(filename_wind,'v10')); 

% ERA5's ssrd is "surface solar radiation downwards" (Surface downwelling shortwave flux in air)
% Divide accumulated hourly flux by 3600 to convert the total number of joules per square meter to average watts per square meter: 
shortwave_downward   = ncread(filename_radiation,'ssrd') / 3600;

% ERA5's strd is "Surface thermal radiation downwards". Divide it by 3600 to convert accumulated energy to average flux:  
longwave_downward    = ncread(filename_radiation,'strd') / 3600;

% Need to convert dewpoint temperature to vapor pressure:
temperature_dewpoint = ncread(filename_2m_temperature,'d2m');
vapor_pressure = dewpoint_to_vapor_pressure(temperature_dewpoint); 

% Build a Climate Forcing timetable:
ClimateForcing = model_initialize_forcing(time_vector,... % time
    temperature_air,...
    pressure_air,...
    precipitation,...      
    wind_speed,...
    shortwave_downward,...
    longwave_downward,...
    vapor_pressure,...
    temperature_air_mean = mean(temperature_air),...
    wind_speed_mean = mean(wind_speed),...
    precipitation_mean = mean(precipitation),...
    temperature_observation_height = 2,...
    wind_observation_height = 10); 
```

This is what the complete climate forcing timetable looks like: 

```matlab
>> head(ClimateForcing)
            time            temperature_air    pressure_air    precipitation    wind_speed    shortwave_downward    longwave_downward    vapor_pressure
    ____________________    _______________    ____________    _____________    __________    __________________    _________________    ______________
    01-Jan-2000 00:00:00        245.19            64100          0.046432         13.465              0                  174.47              44.843    
    01-Jan-2000 01:00:00        245.07            64059          0.047861         12.809              0                  176.84              44.312    
    01-Jan-2000 02:00:00        244.91            64045          0.052294         11.945              0                   177.7              43.752    
    01-Jan-2000 03:00:00        244.75            64006          0.051484         10.962              0                  178.29               43.27    
    01-Jan-2000 04:00:00        244.63            63976          0.041652         10.096              0                  180.72              42.704    
    01-Jan-2000 05:00:00         244.3            63926          0.029244         9.1997              0                  173.98               41.65    
    01-Jan-2000 06:00:00        243.91            63891          0.021435         8.3843              0                  170.73              40.317    
    01-Jan-2000 07:00:00        243.06            63997          0.012085         7.1716              0                  162.37              36.928    
```

And here are the scalar values we specified: 

```matlab
>> ClimateForcing.Properties.CustomProperties
ans = 
CustomProperties with properties:

              temperature_air_mean: 245.9898
                   wind_speed_mean: 5.6444
                precipitation_mean: 0.0208
    temperature_observation_height: 2
           wind_observation_height: 10
```

# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 
