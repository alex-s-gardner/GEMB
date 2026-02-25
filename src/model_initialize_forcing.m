function ClimateForcing = model_initialize_forcing(time, temperature_air, pressure_air, precipitation, ...
    wind_speed, shortwave_downward, longwave_downward, vapor_pressure, options)
% model_initialize_forcing creates a timetable of climate forcing data for GEMB.
%
%% Syntax
%
%  ClimateForcing = model_initialize_forcing(time, temperature_air, pressure_air, precipitation, wind_speed, shortwave_downward, longwave_downward, vapor_pressure)
%  ClimateForcing = model_initialize_forcing(..., temperature_air_mean=value)
%  ClimateForcing = model_initialize_forcing(..., wind_speed_mean=value)
%  ClimateForcing = model_initialize_forcing(..., precipitation_mean=value)
%  ClimateForcing = model_initialize_forcing(..., temperature_observation_height=value)
%  ClimateForcing = model_initialize_forcing(..., wind_observation_height=value)
%
%% Description
%
% ClimateForcing = model_initialize_forcing(time, temperature_air, pressure_air, precipitation, wind_speed, shortwave_downward, longwave_downward, vapor_pressure)
% builds a timetable from an Nx1 time vector in datetime format and corresponding 
% surface forcing vectors that must be of equal length. Inputs must be: 
%    time               datenum  (Nx1)
%    temperature_air    K        (Nx1)
%    pressure_air       Pa       (Nx1)
%    precipitation      kg m^-2  (Nx1) 
%    wind_speed         m s^-1   (Nx1)
%    shortwave_downward W m^-2   (Nx1)
%    longwave_downward  W m^-2   (Nx1)
%    vapor_pressure     Pa       (Nx1)
%
% ClimateForcing = model_initialize_forcing(..., temperature_air_mean=value) 
% specifies a climatological mean air temperature value. If a value is not specified, 
% a warning message will indicate that the mean of temperature_air vector
% is assumed. 
%
% ClimateForcing = model_initialize_forcing(..., wind_speed_mean=value)
% specifies a climatological mean wind speed value. If a value is not specified, 
% a warning message will indicate that the mean of wind_speed vector
% is assumed. 
%
% ClimateForcing = model_initialize_forcing(..., precipitation_mean=value)
% specifies a climatological mean precipitation value. If a value is not specified, 
% a warning message will indicate that the mean of precipitation vector
% is assumed. 
% 
% ClimateForcing = model_initialize_forcing(..., temperature_observation_height=value)
% specifies the height of the observed or modeled temperature_air above the
% surface, in meters. If a value is not specified, a warning message will indicate 
% that 2 m is assumed. 
% 
% ClimateForcing = model_initialize_forcing(..., wind_observation_height=value)
% specifies the height of the observed or modeled wind_speed above the
% surface, in meters. If a value is not specified, a warning message will indicate 
% that 10 m is assumed. 
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 
% 
% See also model_initialize_column and model_initialize_parameters.


%% Check Inputs:

arguments
    time                                   (:,1) datetime
    temperature_air                        (:,1) {mustBeNumeric,mustBeFinite} 
    pressure_air                           (:,1) {mustBeNumeric,mustBeFinite} 
    precipitation                          (:,1) {mustBeNumeric,mustBeFinite} 
    wind_speed                             (:,1) {mustBeNumeric,mustBeFinite} 
    shortwave_downward                     (:,1) {mustBeNumeric,mustBeFinite} 
    longwave_downward                      (:,1) {mustBeNumeric,mustBeFinite} 
    vapor_pressure                         (:,1) {mustBeNumeric,mustBeFinite} 
    options.temperature_air_mean           (1,1) {mustBeNumeric} = NaN % Defaults are redefined below, with a warning message.
    options.wind_speed_mean                (1,1) {mustBeNumeric} = NaN % Defaults are redefined below, with a warning message.
    options.precipitation_mean             (1,1) {mustBeNumeric} = NaN % Defaults are redefined below, with a warning message.
    options.temperature_observation_height (1,1) {mustBeNumeric} = NaN % Defaults are redefined below, with a warning message.
    options.wind_observation_height        (1,1) {mustBeNumeric} = NaN % Defaults are redefined below, with a warning message.

end

assert(isequal(size(time), size(temperature_air), size(pressure_air), size(precipitation), ...
    size(wind_speed), size(shortwave_downward), size(longwave_downward), size(vapor_pressure)), ...
    'All input variables must have the same size.');

assert(all(temperature_air>100),                            "The temperature_air values provided are unrealistic. Ensure the units are kelvin.")
assert(all(pressure_air>=0)   & all(pressure_air<150000),   "The pressure_air values provided are unrealistic for Earth. Ensure the units are pascals.")
assert(all(precipitation>=0)  & all(precipitation<20000),   "The precipitation values provided are unrealistic for Earth. Ensure the units are kg/m^2 per timestep.")
assert(all(wind_speed>=0)     & all(wind_speed<1000),       "The wind_speed values provided are unrealistic for Earth. Ensure the units are m/s.")
assert(all(shortwave_downward<10000),                       "The shortwave_downward values provided are unrealistic for Earth. Ensure the units are W/m^2.")
assert(all(longwave_downward<10000),                        "The longwave_downward values provided are unrealistic for Earth. Ensure the units are W/m^2.")
assert(all(vapor_pressure>=0) & all(vapor_pressure<150000), "The vapor_pressure values provided are unrealistic for Earth. Ensure the units are W/m^2.")

if isnan(options.temperature_air_mean)
    warning("Undeclared temperature_air_mean. Assuming mean(temperature_air) represents the climatological mean temperature.")
    options.temperature_air_mean = mean(temperature_air); 
end

if isnan(options.wind_speed_mean)
    warning("Undeclared wind_speed_mean. Assuming mean(wind_speed_mean) represents the climatological mean wind speed.")
    options.wind_speed_mean = mean(wind_speed);
end

if isnan(options.precipitation_mean)
    warning("Undeclared precipitation_mean. Assuming mean(precipitation_mean) represents the climatological mean precipitation per timestep.")
    options.precipitation_mean = mean(precipitation_mean);
end

if isnan(options.temperature_observation_height)
    warning("Undeclared temperature_observation_height. Assuming air temperature is measured or modeled 2 m above the surface.")
    options.temperature_observation_height = 2;
end

if isnan(options.wind_observation_height)
    warning("Undeclared wind_observation_height. Assuming wind speed is measured or modeled 10 m above the surface.")
    options.wind_observation_height = 10;
end

%% Create Table:

ClimateForcing = timetable(time, temperature_air, pressure_air, precipitation, ...
    wind_speed, shortwave_downward, longwave_downward, vapor_pressure); 

ClimateForcing.Properties.VariableUnits = ["K", "Pa", "kg m^-2", "m s^-1", "W m^-2", "W m^-2", "Pa"]; 

ClimateForcing = addprop(ClimateForcing,...
    ["temperature_air_mean", "wind_speed_mean", "precipitation_mean", "temperature_observation_height", "wind_observation_height"],...
    [            "table",              "table",              "table",                          "table",                   "table"]);

ClimateForcing.Properties.CustomProperties.temperature_air_mean           = options.temperature_air_mean; 
ClimateForcing.Properties.CustomProperties.wind_speed_mean                = options.wind_speed_mean; 
ClimateForcing.Properties.CustomProperties.precipitation_mean             = options.precipitation_mean; 
ClimateForcing.Properties.CustomProperties.temperature_observation_height = options.temperature_observation_height; 
ClimateForcing.Properties.CustomProperties.wind_observation_height        = options.wind_observation_height; 

end

