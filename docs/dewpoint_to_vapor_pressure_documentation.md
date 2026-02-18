# `dewpoint_to_vapor_pressure` documentation
`dewpoint_to_vapor_pressure ` returns the top row of finite data in a matrix A.

# Syntax 
```matlab
vapor_pressure = dewpoint_to_vapor_pressure(temperature_dewpoint)
```

# Description

`vapor_pressure = dewpoint_to_vapor_pressure(temperature_dewpoint)` converts `temperature_dewpoint` (K) to actual vapor pressure in Pa. 

# Example  
This example uses ERA5 forcing data that you can download [here](https://chadagreene.com/GEMB_data/). Learn how to get your own ERA5 data [here](ERA5_time_series_data). 

After downloading the `reanalysis-era5-land-timeseries-sfc-2m-temperature_summit.nc` file, use `ncread` to read the `d2m` variable, which is the dewpoint temperature 2 m above the surface. 

Reading dates of a NetCDF requires using `ncdisp` to check how dates are defined. In this file, the `valid_time` vector represents the number of hours since Jan 1, 1970, so we read in that vector and divide it by 24 to convert hours to days, add it to the datenum referenced to the start of 1970, and then convert from datenum to datetime format so the plot function will automatically recognize the horizontal axis as dates. 

```matlab
filename  = 'reanalysis-era5-land-timeseries-sfc-2m-temperature_summit.nc';

% Read the dewpoint temperature data:
temperature_dewpoint = ncread(filename,'d2m');

% Convert dewpoint temperature to actual vapor pressure: 
vapor_pressure = dewpoint_to_vapor_pressure(temperature_dewpoint); 

% Reference date plus hours since the referene date: 
time = datetime(1970,1,1) + duration(ncread(filename,'valid_time'),0,0);
    
figure
subplot(2,1,1)
plot(time,temperature_dewpoint)
box off
axis tight
ylabel 'Dewpoint temperature (K)'

subplot(2,1,2)
plot(time,vapor_pressure)
box off
axis tight
ylabel 'Vapor pressure (Pa)'

exportgraphics(gcf,'dewpoint_to_vapor_pressure_documentation_01.jpg', ...
    Resolution=300)
```

![](https://chadagreene.com/GEMB_figures/dewpoint_to_vapor_pressure_documentation_01.jpg) 

# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
[https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 
