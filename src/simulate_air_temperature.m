function temperature_air = simulate_air_temperature(dec_year, latitude, elevation, coeffs)
% simulate_air_temperature simulates air temp using fitted coefficients.
%
%% Syntax
%
%  temperature_air = simulate_air_temperature(dec_year, latitude, elevation)
%  temperature_air = simulate_air_temperature(..., lat_scale=value)
%  temperature_air = simulate_air_temperature(..., daily_amp_scale=value)
%  temperature_air = simulate_air_temperature(..., weather_sigma_scale=value)
%  temperature_air = simulate_air_temperature(..., weather_corr=value)
%  temperature_air = simulate_air_temperature(..., mean_offset=value)
%
%% Description
%
% temperature_air = simulate_air_temperature(dec_year, latitude, elevation) simulates near-surface
% air temperatures (K) as a function of decimal year dec_year, latitude latitude,
% and elevation elevation (m). 
%  
% temperature_air = simulate_air_temperature(..., lat_scale=value) specifies scaling
% for annual seasonal amplitude. Must be positive. By default, lat_scale=1. 
%  
% temperature_air = simulate_air_temperature(..., daily_amp_scale=value) specifies
% scaling for diurnal amplitude. Must be positive. By default,
% daily_amp_scale=1.
%  
% temperature_air = simulate_air_temperature(..., weather_sigma_scale=value)
% specifies scaling for weather noise magnitude. Must be positive. By 
% default, weather_sigma_scale=1. 
%  
% temperature_air = simulate_air_temperature(..., weather_corr=value) specifies
% day-to-day weather persistence in a range of 0 to 1. By default,
% weather_corr=0.7.
%  
% temperature_air = simulate_air_temperature(..., mean_offset=value) specifies a mean
% temperature offset. By default, mean_offset=0. 
%
%% Example
% Consider Summit Station in Greenland (72.579583°N, 38.459186°W), 
% whose surface elevation is 3207 m. Simulate hourly surface temperatures
% at Summit for the first half of the 2023: 
% 
%   dec_year = 2023:1/(365*24):2023.5; % hourly time array
%  
%   temperature_air = simulate_air_temperature(dec_year, 72.579583, 3207);
%  
%   figure
%   plot(dec_year,temperature_air,'DisplayName','Default options')
%   ylabel 'Air temperature (K)'
%
% In the figure above, air temperature increases from winter to summer,
% with some randomness to simulate weather noise. Daily cycles are represented 
% as a pure sinusoid with a 24 hour period that remains constant throughout
% the time series. 
% 
% Now assume a mean offset of 9.8842 K, meaning XXXX. Also assume
% weather_corr = 0.7315, meaning YYYY, and weather_sigma_scale = 0.8418
% because ZZZZ. 
% 
%   T_air2 = simulate_air_temperature(dec_year, 72.579583, 3207,...
%       mean_offset = 9.8842,...
%       weather_corr = 0.7315,...
%       weather_sigma_scale = 0.8418);
%  
%   hold on
%   plot(dec_year,T_air2,'DisplayName','Fancy options')
%   legend('location','northwest')
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
% See also: fit_air_temperature.

%% 1. Input Setup & Extraction

arguments
    dec_year                   (:,1) {mustBeNumeric} 
    latitude                        (:,1) {mustBeGreaterThanOrEqual(latitude,-90),mustBeLessThanOrEqual(latitude,90)} 
    elevation                       (:,1) {mustBeGreaterThanOrEqual(elevation,0),mustBeLessThanOrEqual(elevation,10e3)} 
    coeffs.lat_scale           (1,1) {mustBePositive} = 1;
    coeffs.daily_amp_scale     (1,1) {mustBeGreaterThanOrEqual(coeffs.daily_amp_scale,0)} = 1;
    coeffs.weather_sigma_scale (1,1) {mustBePositive} = 1;
    coeffs.weather_corr        (1,1) {mustBeGreaterThanOrEqual(coeffs.weather_corr,0),mustBeLessThanOrEqual(coeffs.weather_corr,1)} = 0.7;
    coeffs.mean_offset         (1,1)                  = 0;
end

lat_scale           = coeffs.lat_scale;
daily_amp_scale     = coeffs.daily_amp_scale;
weather_sigma_scale = coeffs.weather_sigma_scale;
weather_corr        = coeffs.weather_corr;
mean_offset         = coeffs.mean_offset;

%% 2. Climatology (Mean Temp)

phi = deg2rad(latitude); 
T_sea_level = 300 - 50 .* sin(phi).^2;
LapseRate = 0.0065; 
temperature_air_mean = T_sea_level - (elevation .* LapseRate) + mean_offset;

%% 3. Seasonal Cycle (Annual Wave)

% Base amplitude is 3K minimum + 22K scaled by latitude
T_amp_annual = 3 + (22 .* abs(sin(phi)) * lat_scale);

% Hemisphere Logic
phase_shift_annual = zeros(size(latitude));
if isscalar(latitude)
    if latitude > 0, phase_shift_annual = 0.5; end
else
    phase_shift_annual(latitude > 0) = 0.5;
end

year_frac = dec_year - floor(dec_year);
seasonal_signal = cos(2*pi * (year_frac - phase_shift_annual));

%% 4. Diurnal Cycle (Daily Wave)

% A. Calculate time of day (0.0 to 1.0)
day_fraction = mod(dec_year * 365.25, 1);

% B. Amplitude (Diurnal Temperature Range - DTR)
DTR = 10 + (elevation / 1000); 
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
previous_state = rng;
rng(42);
white_noise = randn(num_days, 1) * sigma_weather * sqrt(1 - weather_corr^2);
rng(previous_state);

daily_noise(1) = white_noise(1);
for t = 2:num_days
    daily_noise(t) = (weather_corr * daily_noise(t-1)) + white_noise(t);
end

% C. Interpolate Daily Noise to User Time Steps
user_day_indices = (dec_year * 365.25) - start_day + 1;
weather_signal = interp1(1:num_days, daily_noise, user_day_indices, 'linear', 'extrap');

%% 6. Final Combination

temperature_air = temperature_air_mean...
   + (T_amp_annual .* seasonal_signal) ... % Annual Season
   + (T_amp_daily  .* diurnal_signal) ...  % Daily Cycle
   + weather_signal;                       % Random Weather

end