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

% ERA5's ssrd is "surface solar radiation downwards" (Surface downwelling shortwave flux in air)
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

To include a specified number of model spinup cycles, update the model parameters and re-run GEMB like this: 

```matlab
% Update the model parameters: 
ModelParam.spinup_cycles = 3; 

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

Above we see that by including a spinup cycle, the final results are now stable from the beginning of the time series. Comparing to the non-spun-up version, the results may look very different, but note the different _y_ axis limits. 

## 8. Interpret Firn Air Content 

When snow falls and compacts into firn and ice, air gets trapped in the column. GEMB keeps track of the total `firn_air_content` of the column. Plotting the firn air content time series shows that the total air content in the column tends to increase over time, as expected, as air is only added at the top and never removed from the bottom or sides of the column. In the real world, however, the surface elevation does not exhibit much of a trend at Summit Station, where ice and air slowly exit the sides of the imaginary column due to ice flow divergence. 

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

Near-surface air temperature does not always match the skin temperature of Earth's surface. (This is why ERA5 provides separate 2 m air temperature and skin temperature products.) GEMB uses near-surface air temperature as forcing and accounts for the height of the modeled or observed air temperature above the surface when solving energy balance equations. 

Here we explore the relationship between near-surface air temperature and the surface temperature modeled by GEMB. To get a time series of surface temperature, use the [`surface_timeseries`](surface_timeseries_documentation.md) function to extract the surface values from `OutData.temperature`. 

Below, the optional [`cmocean phase`](https://github.com/chadagreene/cmocean) colormap is used to display day-of-year using a perceptually uniform, cyclic colormap. 

```matlab
% Get a surface temperature time series: 
temperature_surface = surface_timeseries(OutData.temperature); 

figure
scatter(OutData.temperature_air, temperature_surface, 2, ...
    day(OutData.time,"dayofyear"), "filled")

axis image 
hold on
plot(xlim,xlim,'k',LineWidth=1) % 1:1 line
xlabel '2 m air temperature (K) forcing from ERA5'
ylabel 'Surface temperature (K) modeled by GEMB '

cb = colorbar; 
ylabel(cb,'Day of year')
clim([1 365])
cmocean phase % optional colormap 

exportgraphics(gcf,'ERA5_analysis_surface_temperature.jpg',Resolution=300)
```

![](https://chadagreene.com/GEMB_figures/ERA5_analysis_surface_temperature.jpg)

Above, we see that GEMB's modeled surface temperatures are consistently lower than the near-surface air temperature in winter. Lower winter surface temperatures correspond to increased evaporative cooling due to drier, windier conditions in winter. Evaporation and condensation are tracked by GEMB, which we plot below as a function of day of year:

```matlab
figure
scatter(day(OutData.time,"dayofyear"), OutData.evaporation_condensation, 2,...
    year(OutData.time), "filled")
axis([0 366 -0.5 0.5])
xlabel 'Day of year'
ylabel('Daily evaporation\_condensation (kg m^{-2})')
colorbar

exportgraphics(gcf,'ERA5_analysis_evaporation_condensation.jpg',Resolution=300)
```

![](https://chadagreene.com/GEMB_figures/ERA5_analysis_evaporation_condensation.jpg)

Above, positive values indicate net evaporation and negative values indicate net condensation. High wind and low humidity are primarily responsible for increased wintertime evaporation and the reduced surface temperatures. Use the [`vapor_pressure_to_relative_humidity`](vapor_pressure_to_relative_humidity_documentation.md) function to convert vapor pressure to relative humidity and plot it along with hourly wind speed as a function of day of year:

```matlab
relative_humidity = vapor_pressure_to_relative_humidity(CF.vapor_pressure, CF.temperature_air);

figure
subplot(2,1,1)
binscatter(day(CF.time,"dayofyear"), CF.wind_speed, 150)
axis([0 366 0 25])
ylabel 'Wind speed (m s^{-1})'

subplot(2,1,2)
binscatter(day(CF.time,"dayofyear"), relative_humidity, 150)
axis([0 366 50 100])
ylabel 'Relative humidity (%)'
xlabel 'Day of year'

exportgraphics(gcf,'ERA5_analysis_wind_humidity.jpg',Resolution=300)
```

![](https://chadagreene.com/GEMB_figures/ERA5_analysis_wind_humidity.jpg)

## 10. Column temperature 

[`gemb`](gemb_documentation.md)'s output structure contains several MxN fields that represent the evolution of M column grid cells over N output time intervals. Use `imagesc` to plot `OutData.temperature` to see how the grid works: 

```matlab
figure
imagesc(OutData.temperature)
xlabel 'N time steps'
ylabel 'M column grid cells'

exportgraphics(gcf,'ERA5_analysis_temperature.jpg',Resolution=300)
```

![](https://chadagreene.com/GEMB_figures/ERA5_analysis_temperature.jpg)

The dark blue grid cells at the top of the figure above correspond to NaNs that are present as padding in case the column needs to grow beyond its initial size. Indeed, the column did grow during the spinup cycle, which is evident by the fact that the initital padding was set to 1000 grid cells. Checking the dimensions of the data, we can see how the `output_padding` plus the number of grid cells in the initial column define the output height of the grid, and more than 900 unused grid cells of padding in the figure above suggests that the `output_padding` parameter could safely be set to 100 in a future run of this test.

```matlab
>> ModelParam.output_padding
ans =
       1000.00
>> height(Profile)
ans =
        264.00
>> size(OutData.temperature)
ans =
       1264.00       9497.00
```

The temperature data in the figure above are placed vertically in grid-cell space, but the grid cells near the surface are much thinner than the deeper grid cells, and their dimensions evolve with each timestep. Use [`dz2z`](dz2z_documentation.md) to convert from grid-cell space to true vertical dimensions and plot with `pcolor` instead of `imagesc` because `pcolor` allows variable spacing wheras `imagesc` equal spacing. 

```matlab
% Get a 2D matrix of grid cell centers: 
z_center = dz2z(OutData.dz);

% Convert time to 2D so pcolor can plot it: 
time_2D = repmat(OutData.time,size(OutData.temperature,1),1);

figure
pcolor(time_2D,z_center,OutData.temperature)
shading flat
clim([230 273])
ylabel 'Column height (m)'
ylim([-10 1])
cb = colorbar;
ylabel(cb,'Temperature (K)')
cmocean thermal % optional colormap

exportgraphics(gcf,'ERA5_analysis_temperature2.jpg',Resolution=300)
```
![](https://chadagreene.com/GEMB_figures/ERA5_analysis_temperature2.jpg)

Above, we see warming near the surface each summer that decays with depth and travels slowly down the column. 

## 11. Column density and surface melt

Snow and ice layers keep a history of surface forcing. Below we plot density of the top 20 m of the column under a time series of surface melt. The optional [`cmocean dense`](https://github.com/chadagreene/cmocean) colormap is used here. 

```matlab
figure
pcolor(time_2D, z_center, OutData.density); 
shading interp
ylim([-20 8])
box off
ylabel('Column depth (m) or daily surface melt (kg m^{-2})')
cb = colorbar; 
ylabel(cb,'Column density (kg m^{-3})')
clim([350 400])
cmocean dense % optional colormap

% Add a time series of surface melt: 
hold on
plot(OutData.time, OutData.melt, 'linewidth',2)

exportgraphics(gcf,'ERA5_analysis_density.jpg',Resolution=300)
```

![](https://chadagreene.com/GEMB_figures/ERA5_analysis_density.jpg)

Above, notice that surface melt tends to refreeze and create dense layers that get buried by subsequent snowfall events. Also notice that some dense layers are present below the surface throughout the time series. Any layers present at the beginning of the Climate Forcing time series are in fact generated by melt events that were essentially imported from the spinup cycles. 

# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 