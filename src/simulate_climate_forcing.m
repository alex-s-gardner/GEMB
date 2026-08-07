function ClimateForcing = simulate_climate_forcing(set_id, time_step_hours)
% simulate_climate_forcing reproducibly generates synthetic climate forcing data
% for GEMB simulations based on predefined parameter sets.
%
%% Syntax
%
%  ClimateForcing = simulate_climate_forcing(set_id)
%  ClimateForcing = simulate_climate_forcing(set_id, time_step_hours)
%
%% Description
%
% ClimateForcing = simulate_climate_forcing(set_id) reproducibly simulates
% climate forcing based on a set of set_id parameters. The set_id parameters
% are defined in the simulation_parameter_sets function.
%
% ClimateForcing = simulate_climate_forcing(set_id, time_step_hours)
% defines the temporal resolution of the forcing time series in hours. The
% default time step is defined by the set_id.
%
%% Example 1
% Create a climate forcing timetable at the default hourly resolution:
%
%  ClimateForcing = simulate_climate_forcing("test_1");
%
%% Example 2
% Create a climate forcing timetable at three-hourly resoluion:
%
%   ClimateForcing = simulate_climate_forcing("test_1", 3);
%
% View the forcing variables:
%
%   head(ClimateForcing)
%               time            temperature_air    pressure_air    precipitation    wind_speed    shortwave_downward    longwave_downward    vapor_pressure
%       ____________________    _______________    ____________    _____________    __________    __________________    _________________    ______________
%       01-Jan-1994 00:00:00        265.75            92638              0            2.2673            72.661               281.67              247.35
%       01-Jan-1994 03:00:00        266.33            92720              0            2.7777            152.38               265.65              261.84
%
%% Inputs
%
%  set_id                : string       Identifier for the simulation parameter set to use (defined in simulation_parameter_sets.m).
%  time_step_hours       : double       Temporal resolution in hours (optional).
%
%% Outputs
%
%  ClimateForcing        : timetable    Time table containing generated climate data.
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
% See also simulation_parameter_sets and model_initialize_forcing.

% load climate simulation parameter set
[location_parameters, coeffs] = simulation_parameter_sets(set_id);

% determine time step
if nargin == 1
    time_step_hours_val = location_parameters.time_step * (365.25*24); % convert fraction of year to hours
else
    time_step_hours_val = time_step_hours;
end

% create datetime time vector directly. The interval is half-open:
% [start_date-01-01, (end_date+1)-01-01), i.e. the trailing boundary instant
% (end_date+1)-01-01 00:00 is excluded. Including it would append a lone orphan
% step in a new day/year that pollutes daily output binning; dropping it yields
% whole days only (e.g. 32 years * 365.25 d * 8 steps).
time_vec = datetime(location_parameters.start_date, 1, 1) : hours(time_step_hours_val) : ...
           datetime(location_parameters.end_date + 1, 1, 1) - hours(time_step_hours_val);
time_vec = time_vec(:);

% convert datetime to decimal year for internal simulation functions
dec_year = year(time_vec) + (datenum(time_vec) - datenum(datetime(year(time_vec),1,1))) ./ ...
    (datenum(datetime(year(time_vec)+1,1,1)) - datenum(datetime(year(time_vec),1,1)));

rng(location_parameters.rand_seed);

% simulate downward shortave radiation
shortwave_downward  = simulate_shortwave_irradiance(dec_year, location_parameters.latitude);

% simulate air temperature [K]
temperature_air     = simulate_air_temperature(dec_year, location_parameters.latitude, location_parameters.elevation, ...
mean_offset         = coeffs.temperature_air.mean_offset,...
lat_scale           = coeffs.temperature_air.lat_scale ,...
daily_amp_scale     = coeffs.temperature_air.daily_amp_scale,...
weather_sigma_scale = coeffs.temperature_air.weather_sigma_scale,...
weather_corr        = coeffs.temperature_air.weather_corr);

% screen level air pressure [Pa]
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

% create climate forcing timetable using datetime directly
ClimateForcing = model_initialize_forcing(time_vec,...
    temperature_air, pressure_air, precipitation, ...
    wind_speed, shortwave_downward, longwave_downward, vapor_pressure,...
    temperature_air_mean = location_parameters.temperature_air_mean, ...
    precipitation_mean = location_parameters.precipitation_mean, ...
    wind_speed_mean = mean(wind_speed), ...
    temperature_observation_height = location_parameters.temperature_observation_height, ...
    wind_observation_height = location_parameters.wind_observation_height);

end
