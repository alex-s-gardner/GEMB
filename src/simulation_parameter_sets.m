function [location_parameters, coeffs] = simulation_parameter_sets(set_id)
% SIMULATION_PARAMETER_SETS Returns site-specific location and coefficients for simulations.
%
%   [location_parameters, coeffs] = SIMULATION_PARAMETER_SETS(set_id)
%   returns structured data containing geographic parameters and regression
%   coefficients used for generating synthetic environmental time series.
%
%   INPUT:
%       set_id - String or char array identifying the parameter set (e.g., "test_1").
%
%   OUTPUT:
%       location_parameters - Struct containing lat, elev, start_date, end_date, 
%                             time_step, and rand_seed.
%       coeffs - Struct containing sub-structs (Ta, rh, dlw, V) for temperature, 
%                humidity, longwave radiation, and wind speed coefficients.
%
%   EXAMPLE:
%       [loc, cf] = simulation_parameter_sets("test_1");
%
%   See also: simulation_parameters_estimate_from_data

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
        location_parameters.Tmean = 259.4 ; % average annual temerature [K]
        location_parameters.C = 1177.3 ; % average annual accumulation rate of snow or ice [kg m⁻² yr⁻¹]
        location_parameters.time_step = 0.0001 ; % [fraction of a year]
        location_parameters.rand_seed = 42 ; % [seed for random number generator]
         
        %% screen level air temperature [K]
        coeffs.Ta.mean_offset = 9.8847 ;
        coeffs.Ta.lat_scale = 0.1635 ;
        coeffs.Ta.daily_amp_scale = 0.0000 ;
        coeffs.Ta.weather_sigma_scale = 0.8418 ;
        coeffs.Ta.weather_corr = 0.7315 ;
         
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
        coeffs.dlw.beta = [55.8829 0.0948 0.4670 -10.9419 5.5415 ]';
        coeffs.dlw.noise_std = 54.8108 ;
        coeffs.dlw.noise_lag1 = 0.9885 ;
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