function T_air = simulate_air_temperature(dec_year, lat, elev, coeffs)
% SIMULATE_AIR_TEMPERATURE Simulates air temp using fitted coefficients.
%
%   T_air = simulate_air_temperature(dec_year, lat, elev, coeffs)
%
%   INPUTS:
%   dec_year - Decimal year (e.g., 2024.5).
%   lat      - Latitude in degrees.
%   elev     - Elevation in meters.
%   coeffs   - Structure containing fitted parameters (from fit_air_temperature):
%              .lat_scale           (Scaling for annual seasonal amplitude)
%              .daily_amp_scale     (Scaling for diurnal/daily amplitude)
%              .weather_sigma_scale (Scaling for weather noise magnitude)
%              .weather_corr        (Day-to-day weather persistence, 0 to 1)
%              .mean_offset (mean offset adjustment inKelvin [K])
%
%   OUTPUT:
%   T_air       - Near-surface air temperature in Kelvin [K].

    %% 1. Input Setup & Extraction
    if ~isscalar(dec_year), dec_year = dec_year(:); end
    if ~isscalar(lat) && length(lat)==length(dec_year), lat = lat(:); end

    % Extract parameters from struct with defaults if missing
    % (This allows manual struct creation if needed)
    if ~isfield(coeffs, 'lat_scale'), coeffs.lat_scale = 1.0; end
    if ~isfield(coeffs, 'daily_amp_scale'), coeffs.daily_amp_scale = 1.0; end
    if ~isfield(coeffs, 'weather_sigma_scale'), coeffs.weather_sigma_scale = 1.0; end
    if ~isfield(coeffs, 'weather_corr'), coeffs.weather_corr = 0.7; end

    lat_scale = coeffs.lat_scale;
    daily_amp_scale = coeffs.daily_amp_scale;
    weather_sigma_scale = coeffs.weather_sigma_scale;
    weather_corr = coeffs.weather_corr;
    mean_offset = coeffs.mean_offset;

    % Validations
    if weather_corr > 1 || weather_corr < 0
        error("coeffs.weather_corr must be between 0 and 1");
    end
    if lat_scale < 0 || daily_amp_scale < 0 || weather_sigma_scale < 0
        error("Scale coefficients must be positive");
    end

    %% 2. Climatology (Mean Temp)
    phi = deg2rad(lat); 
    T_sea_level = 300 - 50 .* sin(phi).^2;
    LapseRate = 0.0065; 
    T_mean = T_sea_level - (elev .* LapseRate) + mean_offset;

    %% 3. Seasonal Cycle (Annual Wave)
    % Base amplitude is 3K minimum + 22K scaled by latitude
    T_amp_annual = 3 + (22 .* abs(sin(phi)) * lat_scale);
    
    % Hemisphere Logic
    phase_shift_annual = zeros(size(lat));
    if isscalar(lat)
        if lat > 0, phase_shift_annual = 0.5; end
    else
        phase_shift_annual(lat > 0) = 0.5;
    end
    
    year_frac = dec_year - floor(dec_year);
    seasonal_signal = cos(2*pi * (year_frac - phase_shift_annual));

    %% 4. Diurnal Cycle (Daily Wave)
    % A. Calculate time of day (0.0 to 1.0)
    day_fraction = mod(dec_year * 365.25, 1);
    
    % B. Amplitude (Diurnal Temperature Range - DTR)
    DTR = 10 + (elev / 1000); 
    T_amp_daily = (DTR / 2) * daily_amp_scale;
    
    % C. Phase Shift (Peak at ~15:00 or 3 PM = 0.625)
    diurnal_signal = cos(2*pi * (day_fraction - 0.625));
    
    %% 5. Synoptic Weather (Correlated Noise)
    % A. Determine total number of days involved
    start_day = floor(min(dec_year) * 365.25);
    end_day   = ceil(max(dec_year) * 365.25);
    num_days  = end_day - start_day + 1;
    
    % B. Generate Daily Red Noise (AR1)
    base_sigma = 8.0; 
    sigma_weather = base_sigma * weather_sigma_scale; 
    
    daily_noise = zeros(num_days, 1);
    
    % Scale white noise by sqrt(1-alpha^2) to maintain variance over time
    white_noise = randn(num_days, 1) * sigma_weather * sqrt(1 - weather_corr^2);
    
    daily_noise(1) = white_noise(1);
    for t = 2:num_days
        daily_noise(t) = (weather_corr * daily_noise(t-1)) + white_noise(t);
    end
    
    % C. Interpolate Daily Noise to User Time Steps
    user_day_indices = (dec_year * 365.25) - start_day + 1;
    weather_signal = interp1(1:num_days, daily_noise, user_day_indices, 'linear', 'extrap');

    %% 6. Final Combination
    T_air = T_mean ...
       + (T_amp_annual .* seasonal_signal) ... % Annual Season
       + (T_amp_daily  .* diurnal_signal) ...  % Daily Cycle
       + weather_signal;                       % Random Weather

end