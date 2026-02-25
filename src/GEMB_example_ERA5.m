% This script contains the code that's described in greater detail here:
% https://github.com/alex-s-gardner/GEMB/blob/main/docs/ERA5_analysis.md

%% Generate Climate Forcing table: 

save_figures = false; 

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

%% Run GEMB

% Initialize model parameter structure:
ModelParam = model_initialize_parameters(output_frequency="daily");

% Initialize grid:
Profile = model_initialize_column(ModelParam, ClimateForcing);

% Run GEMB: 
OutData = gemb(Profile, ClimateForcing, ModelParam);

%%

figure
plot(OutData.time, OutData.albedo_surface, LineWidth=2)
box off
xlim(datetime([2000 2008], 1, 1))
ylabel 'Surface albedo'

set(gcf,'renderer','painters')

if save_figures
    exportgraphics(gcf,'ERA5_analysis_albedo.jpg',Resolution=300)
end

%%

% Update the model parameters: 
ModelParam.spinup_cycles = 3; 

% Re-run GEMB: 
OutData = gemb(Profile, ClimateForcing, ModelParam);

figure
plot(OutData.time, OutData.albedo_surface, LineWidth=2)
box off
xlim(datetime([2000 2008], 1, 1))
ylabel 'Surface albedo'

set(gcf,'renderer','painters')

if save_figures
    exportgraphics(gcf,'ERA5_analysis_albedo2.jpg',Resolution=300)
end

%%

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

if save_figures
    exportgraphics(gcf,'ERA5_analysis_firn_air_content.jpg',Resolution=300)
end

%%

% Get a surface temperature time series: 
temperature_surface = surface_timeseries(OutData.temperature); 

figure
scatter(OutData.temperature_air, temperature_surface, 2, ...
    day(OutData.time,"dayofyear"), "filled")

axis image
hold on
plot(xlim,xlim,'k',LineWidth=1)
xlabel '2 m air temperature (K) forcing from ERA5'
ylabel 'Surface temperature (K) modeled by GEMB '

cb = colorbar; 
ylabel(cb,'Day of year')
clim([1 365])
cmocean phase % optional colormap 

if save_figures
    exportgraphics(gcf,'ERA5_analysis_surface_temperature.jpg',Resolution=300)
end
%%

figure
scatter(day(OutData.time,"dayofyear"), OutData.evaporation_condensation, 2,...
    year(OutData.time), "filled")
axis([0 366 -0.5 0.5])
xlabel 'Day of year'
ylabel('Daily evaporation\_condensation (kg m^{-2})')
colorbar

if save_figures
    exportgraphics(gcf,'ERA5_analysis_evaporation_condensation.jpg',Resolution=300)
end

%%

relative_humidity = vapor_pressure_to_relative_humidity(ClimateForcing.vapor_pressure, ClimateForcing.temperature_air);

figure
subplot(2,1,1)
binscatter(day(ClimateForcing.time,"dayofyear"), ClimateForcing.wind_speed, 150)
axis([0 366 0 25])
ylabel 'Wind speed (m s^{-1})'

subplot(2,1,2)
binscatter(day(ClimateForcing.time,"dayofyear"), relative_humidity, 150)
axis([0 366 50 100])
ylabel 'Relative humidity (%)'
xlabel 'Day of year'

if save_figures
    exportgraphics(gcf,'ERA5_analysis_wind_humidity.jpg',Resolution=300)
end

%%

figure
imagesc(OutData.temperature)
xlabel 'N time steps'
ylabel 'M column grid cells'

if save_figures
    exportgraphics(gcf,'ERA5_analysis_temperature.jpg',Resolution=300)
end

%%

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

if save_figures
    exportgraphics(gcf,'ERA5_analysis_temperature2.jpg',Resolution=300)
end

%%

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

hold on
plot(OutData.time, OutData.melt, 'linewidth',2)

if save_figures
    exportgraphics(gcf,'ERA5_analysis_density.jpg',Resolution=300)
end
