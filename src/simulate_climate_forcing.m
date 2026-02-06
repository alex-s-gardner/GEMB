function ClimateForcing = simulate_climate_forcing(set_id, time_step_hours)
% simulate_climate_forcing generates synthetic climate forcing data for GEMB
% simulations based on predefined parameter sets.
%
%% Syntax
%
%  ClimateForcing = simulate_climate_forcing(set_id)
%  ClimateForcing = simulate_climate_forcing(set_id, time_step_hours)
%
%% Description
%
% This function orchestrates the generation of synthetic meteorological forcing 
% time series required to drive the GEMB model. It uses a specified parameter 
% set ID to load location-specific constants and coefficients, then calls a 
% suite of sub-functions to simulate individual climate variables. 
%
% The simulation process includes:
% 1. Initialization: Sets up the time vector and initializes the random number 
%    generator for reproducibility.
% 2. Radiation: Simulates downward shortwave irradiance based on solar geometry 
%    and downward longwave irradiance based on air temperature and vapor pressure.
% 3. Thermodynamics: Generates time series for air temperature, air pressure, 
%    relative humidity, and vapor pressure.
% 4. Dynamics & Mass: Simulates wind speed and precipitation events.
% 5. Structuring: Aggregates all generated variables into a standardized 
%    ClimateForcing structure.
%
%% Inputs
%
%  set_id                : integer      Identifier for the simulation parameter set to use (defined in simulation_parameter_sets.m).
%
%% Outputs
%
%  ClimateForcing        : struct       Structure containing generated climate data:
%    .dates              : datenum      Time vector.
%    .shortwave_downward : W m^-2       Downward shortwave radiation.
%    .longwave_downward  : W m^-2       Downward longwave radiation.
%    .temperature_air    : K            Air temperature.
%    .pressure_air       : Pa           Air pressure.
%    .relative_humidity  : %            Relative humidity.
%    .vapor_pressure     : Pa           Vapor pressure.
%    .wind_speed         : m s^-1       Wind speed.
%    .precipitation      : kg m^-2      Precipitation.
%    .latitude           : degrees      Decimal degrees
%    .longitude          : degrees      Decimal degrees
%    .elevation          : meters              
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

% load climate simulation parameter set
[location_parameters, coeffs] = simulation_parameter_sets(set_id);

if nargin == 1
    time_step = location_parameters.time_step;
else
    time_step = time_step_hours / (365.25*24);
end
% initialize times and random seed
dec_year = location_parameters.start_date:time_step:location_parameters.end_date+1;
dec_year = dec_year(:);
dates = decyear2datenum(dec_year);

rng(location_parameters.rand_seed);

% simulate downward shortave radiation
shortwave_downward  = simulate_shortwave_irradiance(dec_year, location_parameters.latitude);

% simulate downward longwave radiation
temperature_air     = simulate_air_temperature(dec_year, location_parameters.latitude, location_parameters.elevation, ...
mean_offset         = coeffs.temperature_air.mean_offset,...
lat_scale           = coeffs.temperature_air.lat_scale ,...
daily_amp_scale     = coeffs.temperature_air.daily_amp_scale,...
weather_sigma_scale = coeffs.temperature_air.weather_sigma_scale,...
weather_corr        = coeffs.temperature_air.weather_corr);

% screen level air temperature [K]
pressure_air = simulate_air_pressure(dec_year, temperature_air, location_parameters.latitude, location_parameters.elevation);

% screen level relative humidity [%]
varname = "relative_humidity";
relative_humidity = simulate_seasonal_daily_noise(dec_year, coeffs.(varname));
relative_humidity(relative_humidity<coeffs.(varname).min_max(1)) = coeffs.(varname).min_max(1);
relative_humidity(relative_humidity>coeffs.(varname).min_max(2)) = coeffs.(varname).min_max(2);

% screen level vapor pressure [Pa]
vapor_pressure = relative_humidity_to_vapor_pressure(temperature_air, relative_humidity);

% downward logwave radiation [W m⁻²]
varname = "longwave_downward";
longwave_downward = simulate_longwave_irradiance(temperature_air, vapor_pressure);
longwave_downward = longwave_downward  + simulate_longwave_irradiance_delta(dec_year, coeffs.(varname));
longwave_downward(longwave_downward<coeffs.(varname).min_max(1)) = coeffs.(varname).min_max(1);
longwave_downward(longwave_downward>coeffs.(varname).min_max(2)) = coeffs.(varname).min_max(2);

% screen level wind speed [m s⁻¹]
varname = "wind_speed";
wind_speed = simulate_seasonal_daily_noise(dec_year, coeffs.(varname));
wind_speed(wind_speed<coeffs.(varname).min_max(1)) = coeffs.(varname).min_max(1);
wind_speed(wind_speed>coeffs.(varname).min_max(2)) = coeffs.(varname).min_max(2);

% precipitation [kg m⁻²]
varname       = "precipitation";
precipitation = simulate_precipitation(dec_year, coeffs.(varname));

% populate structure
ClimateForcing.dates              = dates;
ClimateForcing.shortwave_downward = shortwave_downward;
ClimateForcing.temperature_air    = temperature_air;
ClimateForcing.pressure_air       = pressure_air;
ClimateForcing.relative_humidity  = relative_humidity;
ClimateForcing.vapor_pressure     = vapor_pressure;
ClimateForcing.longwave_downward  = longwave_downward;
ClimateForcing.wind_speed         = wind_speed;
ClimateForcing.precipitation      = precipitation;

% Location specifc parameters
ClimateForcing.wind_observation_height        = location_parameters.wind_observation_height;
ClimateForcing.temperature_observation_height = location_parameters.temperature_observation_height;
ClimateForcing.temperature_air_mean           = location_parameters.temperature_air_mean;
ClimateForcing.wind_speed_mean                = mean(ClimateForcing.wind_speed);
ClimateForcing.precipitation_mean             = location_parameters.precipitation_mean;
ClimateForcing.elevation                      = location_parameters.elevation;
ClimateForcing.latitude                       = location_parameters.latitude;
ClimateForcing.longitude                      = location_parameters.longitude;

end