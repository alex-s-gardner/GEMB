function [location_parameters, coeffs] = simulation_parameter_sets(set_id)
% simulation_parameter_sets Retrieves predefined parameter sets for climate
% forcing simulations.
%
%% Syntax
%
% [location_parameters, coeffs] = simulation_parameter_sets(set_id)
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
%    * Air Temperature (T_air): Mean offsets and scaling factors.
%    * Relative Humidity (rh): Beta coefficients and noise characteristics.
%    * Longwave Radiation (dlw): Gaussian mixture model parameters.
%    * Wind Speed (V): Regression coefficients and noise terms.
%    * Precipitation (P): Harmonic coefficients for occurrence and magnitude.
%
%% Inputs
%
%  set_id              : string       Identifier for the desired parameter set (e.g., "test_1").
%
%% Outputs
%
%  location_parameters : struct       Structure containing site metadata:
%    .lat, .lon, .elev : double       Coordinates and elevation.
%    .start_date       : double       Simulation start year (decimal).
%    .end_date         : double       Simulation end year (decimal).
%    .T_air_mean       : K            Mean annual temperature.
%    .P_mean           : kg m^-2      Mean annual accumulation.
%  coeffs              : struct       Structure containing statistical generation coefficients:
%    .T_air            : struct       Temperature coefficients.
%    .rh               : struct       Relative humidity coefficients.
%    .dlw              : struct       Longwave radiation coefficients.
%    .V                : struct       Wind speed coefficients.
%    .P                : struct       Precipitation coefficients.
%
%% Documentation
%
% For complete documentation, see: https://github.com/alex-s-gardner/GEMB
%
%% References
% If you use GEMB, please cite the following:
%
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass
% Balance (GEMB): a model of firn processes for cryosphere research, Geosci.
% Model Dev., 16, 2277–2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023.

    % list valid sets to give user a hint
    valid_sets = ["test_1", ];

    if set_id == "test_1"
        %%  location and time parameters
        location_parameters.description = "parameters estimated using simulation_parameters_estimate_from_data.m as fit to original TEST_INPUT_1.mat data";
        location_parameters.lat = -73.3307 ; % [º]
        location_parameters.lon = 290.6250 ; % [º]
        location_parameters.elev = 700 ; % [m]
        location_parameters.start_date = 1994.00 ; % [decimal year]
        location_parameters.end_date = 2025.00 ; % [decimal year]
        location_parameters.Vz = 10.0 ; % wind observation height above surface [m]
        location_parameters.Tz = 2.0 ; % temperature observation height above surface [m]
        location_parameters.T_air_mean = 259.4 ; % average annual temerature [K]
        location_parameters.P_mean = 1177.3 ; % average annual accumulation rate of snow or ice [kg m⁻² yr⁻¹]
        location_parameters.time_step = 0.0001 ; % [fraction of a year]
        location_parameters.rand_seed = 42 ; % [seed for random number generator]
         
        %% screen level air temperature [K]
        coeffs.T_air.mean_offset = 9.8847 ;
        coeffs.T_air.lat_scale = 0.1635 ;
        coeffs.T_air.daily_amp_scale = 0.0000 ;
        coeffs.T_air.weather_sigma_scale = 0.8418 ;
        coeffs.T_air.weather_corr = 0.7315 ;
         
        %% screen level relative humidity [%]
        coeffs.rh.beta = [79.9128 -0.2619 1.3055 5.0727 0.5081 ]';
        coeffs.rh.noise_std = 11.9734 ;
        coeffs.rh.noise_lag1 = 0.9896 ;
        coeffs.rh.min_max = [0.0000 100.0000 ]';
         
        coeffs.rh.beta = [79.9128 -0.2619 1.3055 5.0727 0.5081 ]';
        coeffs.rh.noise_std = 11.9734 ;
        coeffs.rh.noise_lag1 = 0.9896 ;
        coeffs.rh.min_max = [0.0000 100.0000 ]';
         
        %% downward logwave radiation [W m⁻²]
        coeffs.dlw.mu = [113.2225 23.6664 ]';
        coeffs.dlw.sigma = [30.9535 38.1762 ]';
        coeffs.dlw.P = [0.3590 0.6410 ];
        coeffs.dlw.limits = [-142.9820 222.0189 ];
        coeffs.dlw.prctile_bounds = [-67.5562 147.9489 ];
        coeffs.dlw.min_max = [0.0000 Inf ]';
                 
        %% screen level wind speed [m s⁻¹]
        coeffs.V.beta = [5.2495 0.1296 -0.1237 -0.3030 -0.2289 ]';
        coeffs.V.noise_std = 2.9324 ;
        coeffs.V.noise_lag1 = 0.9800 ;
        coeffs.V.min_max = [0.0000 Inf ]';

        %% precipitation [kg m⁻²]
        coeffs.P.P01_harmonics = [0.0213 -0.0002 -0.0019 ]';
        coeffs.P.P11_harmonics = [0.9443 -0.0013 -0.0036 ]';
        coeffs.P.Alpha_harmonics = [1.1357 0.0610 0.1906 ]';
        coeffs.P.Beta_harmonics = [0.4060 -0.0376 -0.1312 ]';
        coeffs.P.wet_threshold = 0.1000 ;
    else
        msg = strjoin(valid_sets, ', '); 
        error(set_id + " is not a valid simulation parameter set. Valid sets include: [" + msg + "]")
    end
    
end