function [location_parameters, coeffs] = simulation_parameters_test1()
    %%  location and time parameters
    location_parameters.lat = -73.3307 ; % [º]
    location_parameters.elev = 700 ; % [m]
    location_parameters.start_date = 1994.00 ; % [decimal year]
    location_parameters.end_date = 2025.00 ; % [decimal year]
    location_parameters.time_step = 0.0001 ; % [fraction of a year]
    location_parameters.rand_seed = 42 ; % [seed for random number generator]
     
    %% screen level air temperature [K]
    coeffs.(varname).Ta.mean_offset = 9.8847 ;
    coeffs.(varname).Ta.lat_scale = 0.1635 ;
    coeffs.(varname).Ta.daily_amp_scale = 0.0000 ;
    coeffs.(varname).Ta.weather_sigma_scale = 0.8418 ;
    coeffs.(varname).Ta.weather_corr = 0.7315 ;
     
    %% screen level relative humidity [%]
    coeffs.(varname).rh.beta = [79.9128 -0.2619 1.3055 5.0727 0.5081 ];
    coeffs.(varname).rh.noise_std = 11.9734 ;
    coeffs.(varname).rh.noise_lag1 = 0.9896 ;
    coeffs.(varname).rh.min_max = [0.0000 100.0000 ];
     
    coeffs.(varname).rh.beta = [79.9128 -0.2619 1.3055 5.0727 0.5081 ];
    coeffs.(varname).rh.noise_std = 11.9734 ;
    coeffs.(varname).rh.noise_lag1 = 0.9896 ;
    coeffs.(varname).rh.min_max = [0.0000 100.0000 ];
     
    %% downward logwave radiation [W m⁻²]
    coeffs.(varname).dlw.beta = [56.5485 0.0995 0.4739 -12.7450 7.7032 ];
    coeffs.(varname).dlw.noise_std = 55.1521 ;
    coeffs.(varname).dlw.noise_lag1 = 0.9886 ;
     
    %% screen level wind speed [m s⁻¹]
    coeffs.(varname).V.beta = [5.2495 0.1296 -0.1237 -0.3030 -0.2289 ];
    coeffs.(varname).V.noise_std = 2.9324 ;
    coeffs.(varname).V.noise_lag1 = 0.9800 ;
    coeffs.(varname).V.min_max = [0.0000 Inf ];
    
end