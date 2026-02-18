# ERA5 time series analysis

This tutorial describes how to use ERA5 reanalysis data to model the firn column in an accumulation zone  near Greenland's [Summit Station](https://en.wikipedia.org/wiki/Summit_Camp).  

## 0. Download data

To follow along, download the same data files that we explore below [here](https://chadagreene.com/GEMB_data/), or check out [this tutorial](ERA5_time_series_data.md) to see how to get the same type of data for any other location on Earth. 

## 1. Define Climate Forcing
The [`gemb`](gemb_documentation.md) function requires a Climate Forcing structure that contains the following input variables: 

| `gemb` variable        | ERA5 source     | Conversion notes                                                                                                                                  |
|------------------------|-----------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| `time`                 | `valid_time`    | Convert "hours since 1970-01-01" to datetime or datenum format.   |
| `temperature_air`      | `t2m`           | No conversion necessary. |
| `pressure_air`         | `sp`            | No conversion necessary.|
| `precipitation`        | `tp`            | Divide by 1000 to convert from m to kg m<sup>-2</sup>   |
| `wind_speed`           | `u10` and `v10` | Take the hypotenuse of the u and v components.  |
| `shortwave_downward`   | `ssrd`          | Divide by 3600 to convert from joules per hour to watts.  |
| `longwave_downward`    | `strd`          | Divide by 3600 to convert from joules per hour to watts.   |
| `vapor_pressure`       | `d2m`           | Requires conversion from dewpoint temperature to vapor pressure. Use [`dewpoint_to_vapor_pressure`](dewpoint_to_vapor_pressure_documentation.md). |
| `temperature_air_mean` | `t2m`           | Can use mean of `temperature_air` if the time series is long enough to represent climatological average.  |
| `wind_speed_mean`      | `u10` and `v10` | Can use mean of `wind_speed` if the time series is long enough to represent climatological average.  |
| `precipitation_mean`   | `tp`            | Can use mean of `precipitation` if the time series is long enough to represent climatological average.   |
| `temperature_observation_height`      |  | The modeled 2 m air temperature is equivalent to a thermometer placed 2 m above the ground. |
| `wind_observation_height`   |  | The modeled 10 m wind speed is equivalent to an anemometer placed 10 m above the ground. |

Defining a Climate Forcing structure `CF` for the [`gemb`](gemb_documentation.md) function will require us to use the table above to fill in the 13 variables in the following blank template: 

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

Assuming you've downloaded the example data and placed it in your current directory or another location where MATLAB can find it, define each of the filenames:

```matlab
filename_2m_temperature  = 'reanalysis-era5-land-timeseries-sfc-2m-temperature_summit.nc';
filename_pressure_precip = 'reanalysis-era5-land-timeseries-sfc-pressure-precipitation_summit.nc';
filename_radiation       = 'reanalysis-era5-land-timeseries-sfc-radiation-heat_summit.nc';
filename_wind            = 'reanalysis-era5-land-timeseries-sfc-wind_summit.nc';
```

MATLAB doesn't offer an elegant way to read time from a NetCDF, so you have to build it yourself. Start by looking at how the `valid_time` variable is defined in the file: 

```matlab
>> ncdisp(filename_2m_temperature, 'valid_time')
Source:
           /Users/cgreene/Documents/GEMB_data/reanalysis-era5-land-timeseries-sfc-2m-temperature_summit.nc
Format:
           netcdf4
Dimensions:
           valid_time = 227928
Variables:
    valid_time
           Size:       227928x1
           Dimensions: valid_time
           Datatype:   int64
           Attributes:
                       units    = 'hours since 1970-01-01'
                       calendar = 'proleptic_gregorian'
```

The units of time in the ERA5 NetCDF are `'hours since 1970-01-01'`, so create a time variable by defining a reference date in datetime format and add the number of hours since the reference date. Here, we'll let `time` be the the first variable defined in a Climate Forcing structure we'll call `CF`:

```matlab
% Reference date plus number of hours since the reference date: 
CF.time = datetime(1970,1,1) + duration(ncread(filename_2m_temperature,'valid_time'),0,0);
```

The remaining variables in the Climate Forcing structure `CF` are defined individually by reading the ERA5 variables listed in the table at the top of this page. Note that some units are converted from units used by ERA5 to units required by GEMB. Also note that while ERA5 does not explicitly contain vapor pressure data, GEMB provides a [`dewpoint_to_vapor_pressure`](dewpoint_to_vapor_pressure_documentation.md) function that converts dewpoint temperature to corresponding vapor pressure.  

```matlab
CF.temperature_air      = ncread(filename_2m_temperature,'t2m'); 
CF.pressure_air         = ncread(filename_pressure_precip,'sp');

% Multiply hourly precipitation by 1000 to convert from m to kg/m^2:  
CF.precipitation        = ncread(filename_pressure_precip,'tp') * 1000;

% Wind speed is the hypotenuse of the vector components: 
CF.wind_speed           = hypot(ncread(filename_wind,'u10'),ncread(filename_wind,'v10')); 

% ERA5's ssrd is "surface solar radiation downwards" (Surface downlwelling shortwave flux in air)
% Divide accumulated hourly flux by 3600 to convert the total number of joules per square meter to average watts per square meter: 
CF.shortwave_downward   = ncread(filename_radiation,'ssrd') / 3600;

% ERA5's strd is "Surface thermal radiation downwards". Divide it by 3600 to convert accumulated energy to average flux:  
CF.longwave_downward    = ncread(filename_radiation,'strd') / 3600;

temperature_dewpoint    = ncread(filename_2m_temperature,'d2m');
CF.vapor_pressure       = dewpoint_to_vapor_pressure(temperature_dewpoint); 

CF.temperature_air_mean = mean(CF.temperature_air); 

CF.wind_speed_mean      = mean(CF.wind_speed); 
CF.precipitation_mean   = mean(CF.precipitation); 

CF.temperature_observation_height = 2;  % 2 m air temperature
CF.wind_observation_height        = 10; % 10 m wind speed
```
## 3. Define Model Parameters

Use the [`model_initialize_parameters`](model_initialize_parameters_documentation.md) function to intialize model parameters. Although the climate forcing data are hourly, we will output to daily frequency to minimize the amount of data we will need to plot later: 

```matlab
ModelParam = model_initialize_parameters(output_frequency="daily");
```

## 4. Initialize a GEMB Column

Use the [`model_initialize_column`](model_initialize_column_documentation.md) function to create an initial GEMB column based on our climate forcing structure `CF` and model parameters `ModelParam`: 

```matlab
Profile = model_initialize_column(ModelParam, CF);
```

## 5. Run GEMB

With the climate forcing, model parameters, and column all set, running GEMB is very easy, although it may take ~30 seconds to run this example data: 

```matlab
OutData = gemb(Profile, CF, ModelParam);
```

## 6. Visualize albedo evolution

It's always a good idea to check the data after running GEMB. Here we'll look at the time series of surface albedo and zoom in to the first 8 years of the time series so we can what's going on: 

```
figure
plot(OutData.time, OutData.albedo_surface, LineWidth=2)
box off
xlim(datetime([2000 2008], 1, 1))
ylabel 'Surface albedo'

exportgraphics(gcf,'ERA5_analysis_albedo.jpg',Resolution=300)
```

![](https://chadagreene.com/GEMB_figures/ERA5_analysis_albedo.jpg)

Above, we see that the surface albedo generally hovers around 0.846, and drops a little on some warm summer days, as expected. Also notice the _much_ lower albedo values in the first few weeks of the time series, indicating that the model needs some spinup time for surface albedo to stabilize. 

## 7. Update Model Parameters 

To include a model spinup cycle, update the model parameters and re-run GEMB like this: 

```matlab
% Update the model parameters: 
ModelParam.spinup_cycles = 1; 

% Re-run GEMB: 
OutData = gemb(Profile, CF, ModelParam);

figure
plot(OutData.time, OutData.albedo_surface, LineWidth=2)
box off
xlim(datetime([2000 2008], 1, 1))
ylabel 'Surface albedo'

exportgraphics(gcf,'ERA5_analysis_albedo2.jpg',Resolution=300)
```

![](https://chadagreene.com/GEMB_figures/ERA5_analysis_albedo2.jpg)

Above we see that by including a spinup cycle, the final results are now stable from the beginning of the time series. Comparing to the non-spun-up version, the results may look very different, but note the differeny _y_ axis limits. 

## 8. Interpret Firn Air Content 

When snow falls and compacts into firn and ice, air gets trapped in the column. GEMB keeps track of the total `firn_air_content` of the column. Plotting the firn air content time series shows that the total air content in the column tends to increase over time, as expected, as air is only added at the top and never removed from the bottom or sides of the column. In the real world, however, the surface elevation does not show much of a trend at Summit Station, where ice and air slowly exit the sides of the imaginary column due to ice flow divergence. 

The easiest way to interpret firn air content anomalies is to assume that the total firn air content tends to be stable over long time horizons if the climate and ice dynamics remain stable. With this assumption, we can simply detrend the `firn_air_content` time series to analyze shorter-term anomalies in the residuals.  

```matlab
% Detrend FAC for an "anomaly" time series: 
firn_air_content_anomaly = detrend(OutData.firn_air_content); 

figure
subplot(2,1,1)
plot(OutData.time, OutData.firn_air_content)
box off
ylabel 'Firn air content (m)'

subplot(2,1,2)
plot(OutData.time, firn_air_content_anomaly)
box off
ylabel 'Firn air content anomaly (m)'

exportgraphics(gcf,'ERA5_analysis_firn_air_content.jpg',Resolution=300)
```

![](https://chadagreene.com/GEMB_figures/ERA5_analysis_firn_air_content.jpg)

## 9. Surface temperature vs air temperature


## 10. Column temperature 


## 11. Column melt and density 


# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
[https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 