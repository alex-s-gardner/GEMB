function ClimateForcing = simulate_climate_forcing(set_id, time_step_hours)
% simulate_climate_forcing generates synthetic climate forcing data for GEMB
% simulations based on predefined parameter sets.
%
%% Syntax
%
% ClimateForcing = simulate_climate_forcing(set_id)
% ClimateForcing = simulate_climate_forcing(set_id, time_step_hours)
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
%  set_id              : integer      Identifier for the simulation parameter set to use (defined in simulation_parameter_sets.m).
%
%% Outputs
%
%  ClimateForcing      : struct       Structure containing generated climate data:
%    .daten            : datenum      Time vector.
%    .dsw0             : W m^-2       Downward shortwave radiation.
%    .dlw0             : W m^-2       Downward longwave radiation.
%    .T_air0           : K            Air temperature.
%    .p_air0           : Pa           Air pressure.
%    .rh0              : %            Relative humidity.
%    .e_air0           : Pa           Vapor pressure.
%    .V0               : m s^-1       Wind speed.
%    .P0               : kg m^-2      Precipitation.
%    .lat, .lon, .elev : double       Location metadata (latitude, longitude, elevation).
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
daten = decyear2datenum(dec_year);
%% 
rng(location_parameters.rand_seed);

% simulate downward shortave radiation
dsw0 = simulate_shortwave_irradiance(dec_year, location_parameters.lat);

% simulate downward longwave radiation
T_air0 = simulate_air_temperature(dec_year, location_parameters.lat, location_parameters.elev, ...
    mean_offset = coeffs.T_air.mean_offset,...
    lat_scale=coeffs.T_air.lat_scale ,...
    daily_amp_scale=coeffs.T_air.daily_amp_scale,...
    weather_sigma_scale=coeffs.T_air.weather_sigma_scale,...
    weather_corr=coeffs.T_air.weather_corr);

% screen level air temperature [K]
p_air0 = simulate_air_pressure(dec_year, T_air0, location_parameters.lat, location_parameters.elev);

% screen level relative humidity [%]
varname = "rh";
rh0 = simulate_seasonal_daily_noise(dec_year, coeffs.(varname));
rh0(rh0<coeffs.(varname).min_max(1)) = coeffs.(varname).min_max(1);
rh0(rh0>coeffs.(varname).min_max(2)) = coeffs.(varname).min_max(2);

% screen level vapor pressure [Pa]
e_air0 = simulate_vapor_pressure(T_air0, rh0);

% downward logwave radiation [W m⁻²]
varname = "dlw";
dlw0 = simulate_longwave_irradiance(T_air0, e_air0);
dlw0 = dlw0  + simulate_longwave_irradiance_delta(dec_year, coeffs.(varname));
dlw0(dlw0<coeffs.(varname).min_max(1)) = coeffs.(varname).min_max(1);
dlw0(dlw0>coeffs.(varname).min_max(2)) = coeffs.(varname).min_max(2);

% screen level wind speed [m s⁻¹]
varname = "V";
V0 = simulate_seasonal_daily_noise(dec_year, coeffs.(varname));
V0(V0<coeffs.(varname).min_max(1)) = coeffs.(varname).min_max(1);
V0(V0>coeffs.(varname).min_max(2)) = coeffs.(varname).min_max(2);

% precipitation [kg m⁻²]
varname = "P";
P0 = simulate_precipitation(dec_year, coeffs.(varname));

% populate structure
ClimateForcing.daten = daten;
ClimateForcing.dsw0 = dsw0;
ClimateForcing.T_air0 = T_air0;
ClimateForcing.p_air0 = p_air0;
ClimateForcing.rh0 = rh0;
ClimateForcing.e_air0 = e_air0;
ClimateForcing.dlw0 = dlw0;
ClimateForcing.V0 = V0;
ClimateForcing.P0 = P0;

% Location specifc parameters
ClimateForcing.Vz         = location_parameters.Vz;
ClimateForcing.Tz         = location_parameters.Tz;
ClimateForcing.T_air_mean = location_parameters.T_air_mean;
ClimateForcing.V_mean     = mean(ClimateForcing.V0);
ClimateForcing.P_mean     = location_parameters.P_mean;
ClimateForcing.elev       = location_parameters.elev;
ClimateForcing.lat        = location_parameters.lat;
ClimateForcing.lon        = location_parameters.lon;

end