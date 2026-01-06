function p_air = simulate_air_pressure(dec_year, T_air, lat, elev)
% SIMULATE_AIR_PRESSURE Simulates screen-level atmospheric pressure [Pa].
%
%   p_air = simulate_air_pressure(dec_year, T_air, lat, elev)
%
%   INPUTS:
%   dec_year - Decimal year (e.g., 2024.5). Vector.
%   T_air       - Near-surface air temperature in Kelvin [K]. Vector.
%   lat      - Latitude in degrees.
%   elev     - Elevation in meters.
%
%   OUTPUT:
%   p_air     - Screen-level air pressure in Pascals [Pa].
%
%   PHYSICS MODEL:
%   1. Simulates Mean Sea Level Pressure (MSLP) as a stochastic process.
%      - Base: 101,325 Pa.
%      - Weather: AR(1) "Red Noise" simulating High/Low pressure systems.
%      - Variance: Scales with latitude (Storm tracks have higher variance).
%   2. Calculates Local Pressure using the Hypsometric Equation.
%      - Uses inputs 'elev' and 'T_air' to determine air column density.
%      - P_local = P_mslp * exp( -g*z / (R * T_column) )

    %% 1. Input Sanitization
    if ~isscalar(dec_year), dec_year = dec_year(:); end
    if ~isscalar(T_air), T_air = T_air(:); end
    
    % Ensure time and temp vectors match
    if length(T_air) ~= length(dec_year)
        error('Input vectors dec_year and T_air must be the same length.');
    end

    %% 2. Simulate Mean Sea Level Pressure (MSLP) Weather Patterns
    % MSLP fluctuates around 1013.25 hPa due to synoptic weather (storms/highs).
    
    % A. Define Weather Variance based on Latitude
    %    Tropics have stable pressure; Mid-latitudes (storm tracks) vary wildly.
    %    Low Lat (0): ~300 Pa sigma. High Lat (60): ~1200 Pa sigma.
    phi = deg2rad(lat);
    sigma_pressure = 300 + (900 .* sin(phi).^2); % Pa
    
    % B. Generate Red Noise (AR1 Process)
    %    Pressure systems move slower than daily temp fluctuations.
    %    Correlation ~0.8 to 0.9 for daily data.
    alpha = 0.85; 
    
    n = length(dec_year);
    mslp_noise = zeros(n, 1);
    
    % Interpolate noise to avoid "static" if user inputs sub-daily steps
    % (Similar logic to the temperature function)
    start_day = floor(min(dec_year) * 365.25);
    end_day   = ceil(max(dec_year) * 365.25);
    num_days  = end_day - start_day + 1;
    
    daily_noise = zeros(num_days, 1);
    white_noise = randn(num_days, 1) * sigma_pressure * sqrt(1 - alpha^2);
    
    daily_noise(1) = white_noise(1);
    for t = 2:num_days
        daily_noise(t) = (alpha * daily_noise(t-1)) + white_noise(t);
    end
    
    % Map daily noise to user time steps
    user_day_indices = (dec_year * 365.25) - start_day + 1;
    weather_anomaly = interp1(1:num_days, daily_noise, user_day_indices, 'linear', 'extrap');
    
    % C. Final Simulated MSLP
    P_msl_std = 101325; % Standard mean sea level pressure (Pa)
    P_msl = P_msl_std + weather_anomaly;

    %% 3. Calculate Local Station Pressure (Hypsometric Reduction)
    % We need to estimate the temperature of the air column between sea level
    % and the station elevation.
    % Approximation: T_column_avg = T_surface + (0.5 * Lapse_Rate * Elevation)
    
    LapseRate = 0.0065; % K/m
    T_column_avg = T_air + (0.5 * LapseRate * elev);
    
    % Physical Constants
    g = 9.80665; % Gravity (m/s^2)
    R = 287.05;  % Specific gas constant for dry air (J/kg*K)
    
    % Hypsometric Equation
    % P_local = P_msl * exp( -g*z / (R*T) )
    exponent = (-g * elev) ./ (R .* T_column_avg);
    
    p_air = P_msl .* exp(exponent);

end