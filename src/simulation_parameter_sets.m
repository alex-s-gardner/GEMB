function [location_parameters, coeffs] = simulation_parameter_sets(set_id)
% simulation_parameter_sets retrieves predefined parameter sets for climate
% forcing simulations.
%
%% Syntax
%
%  [location_parameters, coeffs] = simulation_parameter_sets(set_id)
%
%% Description
%
% This function acts as a library of calibrated parameter sets for the GEMB 
% climate forcing generator. Based on the requested identification string 
% (set_id), it returns site-specific geographic constants and statistical 
% coefficients required to simulate synthetic weather data.
%
% The function defines:
% 1. Location Parameters: Geographic metadata (latitude, longitude, elevation), 
%    temporal bounds (start/end dates), and long-term climatological means 
%    (temperature, accumulation).
% 2. Statistical Coefficients: Empirical coefficients used to generate stochastic 
%    time series for specific climate variables:
%    * Air Temperature (temperature_air): Mean offsets and scaling factors.
%    * Relative Humidity (relative_humidity): Beta coefficients and noise characteristics.
%    * Longwave Radiation (longwave_downward): Gaussian mixture model parameters.
%    * Wind Speed (wind_speed): Regression coefficients and noise terms.
%    * Precipitation (precipitation): Harmonic coefficients for occurrence and magnitude.
%
%% Inputs
%
%  set_id              : string       Identifier for the desired parameter set (e.g., "test_1").
%
%% Outputs
%
%  location_parameters     : struct     Structure containing site metadata:
%    .latitude             : double     decimal degrees
%    .longitude            : double     decimal degrees
%    .elevation            : double     elevation (m).
%    .start_date           : double     Simulation start year (decimal).
%    .end_date             : double     Simulation end year (decimal).
%    .temperature_air_mean : K          Mean annual temperature.
%    .precipitation_mean   : kg m^-2    Mean annual accumulation.
%  coeffs                  : struct     Structure containing statistical generation coefficients:
%    .temperature_air      : K          Temperature coefficients.
%    .relative_humidity    : %          Relative humidity coefficients.
%    .longwave_downward    :            Longwave radiation coefficients.
%    .wind_speed           : m/s        Wind speed coefficients.
%    .precipitation        :            Precipitation coefficients.
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

% list valid sets to give user a hint
valid_sets = ["test_1", ];

if set_id == "test_1"

    % Location and time parameters:
    location_parameters.description                    = "parameters estimated using simulation_parameters_estimate_from_data.m as fit to original TEST_INPUT_1.mat data";
    location_parameters.latitude                       = -73.3307 ;     % [º]
    location_parameters.longitude                      = 290.6250 ;     % [º]
    location_parameters.elevation                      = 700 ;          % [m]
    location_parameters.start_date                     = 1994.00 ;      % [decimal year]
    location_parameters.end_date                       = 2025.00 ;      % [decimal year]
    location_parameters.wind_observation_height        = 10.0 ;         % wind observation height above surface [m]
    location_parameters.temperature_observation_height = 2.0 ;          % temperature observation height above surface [m]
    location_parameters.temperature_air_mean           = 259.4 ;        % average annual temerature [K]
    location_parameters.precipitation_mean             = 1177.3 ;       % average annual accumulation rate of snow or ice [kg m⁻² yr⁻¹]
    location_parameters.time_step                      = 1/(365.25*24); % [fraction of a year]
    location_parameters.rand_seed                      = 42 ;           % [seed for random number generator]
     
    % Screen level air temperature [K]:
    coeffs.temperature_air.mean_offset         = 9.8847 ;
    coeffs.temperature_air.lat_scale           = 0.1635 ;
    coeffs.temperature_air.daily_amp_scale     = 0.0000 ;
    coeffs.temperature_air.weather_sigma_scale = 0.8418 ;
    coeffs.temperature_air.weather_corr        = 0.7315 ;
     
    % Screen level relative humidity [%]
    coeffs.relative_humidity.beta       = [79.9128 -0.2619 1.3055 5.0727 0.5081 ]';
    coeffs.relative_humidity.noise_std  = 11.9734 ;
    coeffs.relative_humidity.noise_lag1 = 0.9896 ;
    coeffs.relative_humidity.min_max    = [0.0000 100.0000 ]';
     
    % Downward logwave radiation [W m⁻²]
    coeffs.longwave_downward.mu             = [ 113.2225  23.6664 ]';
    coeffs.longwave_downward.sigma          = [  30.9535  38.1762 ]';
    coeffs.longwave_downward.P              = [   0.3590   0.6410 ];
    coeffs.longwave_downward.limits         = [-142.9820 222.0189 ];
    coeffs.longwave_downward.prctile_bounds = [ -67.5562 147.9489 ];
    coeffs.longwave_downward.min_max        = [   0.0000 Inf      ]';
             
    % Screen level wind speed [m s⁻¹]
    coeffs.wind_speed.beta       = [5.2495 0.1296 -0.1237 -0.3030 -0.2289 ]';
    coeffs.wind_speed.noise_std  = 2.9324 ;
    coeffs.wind_speed.noise_lag1 = 0.9800 ;
    coeffs.wind_speed.min_max    = [0.0000 Inf ]';

    % Precipitation [kg m⁻²]
    coeffs.precipitation.P01_harmonics   = [0.0213 -0.0002 -0.0019 ]';
    coeffs.precipitation.P11_harmonics   = [0.9443 -0.0013 -0.0036 ]';
    coeffs.precipitation.Alpha_harmonics = [1.1357  0.0610  0.1906 ]';
    coeffs.precipitation.Beta_harmonics  = [0.4060 -0.0376 -0.1312 ]';
    coeffs.precipitation.wet_threshold   = 0.1000 ;
else
    msg = strjoin(valid_sets, ', '); 
    error(set_id + " is not a valid simulation parameter set. Valid sets include: [" + msg + "]")
end
    
end