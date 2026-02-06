function pressure_air = simulate_air_pressure(dec_year, temperature_air, latitude, elevation)
% simulate_air_pressure simulates screen-level atmospheric pressure [Pa].
%
%   pressure_air = simulate_air_pressure(dec_year, temperature_air, latitude, elevation)
%
%   INPUTS:
%   dec_year - Decimal year (e.g., 2024.5). Vector.
%   temperature_air       - Near-surface air temperature in Kelvin [K]. Vector.
%   latitude      - Latitude in degrees.
%   elevation     - Elevation in meters.
%
%   OUTPUT:
%   pressure_air     - Screen-level air pressure in Pascals [Pa].
%
%   PHYSICS MODEL:
%   1. Simulates Mean Sea Level Pressure (MSLP) as a stochastic process.
%      - Base: 101,325 Pa.
%      - Weather: AR(1) "Red Noise" simulating High/Low pressure systems.
%      - Variance: Scales with latitude (Storm tracks have higher variance).
%   2. Calculates Local Pressure using the Hypsometric Equation.
%      - Uses inputs 'elevation' and 'temperature_air' to determine air column density.
%      - P_local = P_mslp * exp( -g*z / (R * T_column) )
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

%% 1. Input Sanitization

% Columnate: 
if ~isscalar(dec_year)
    dec_year = dec_year(:); 
end

if ~isscalar(temperature_air)
    temperature_air = temperature_air(:); 
end

% Ensure time and temp vectors match
if length(temperature_air) ~= length(dec_year)
    error('Input vectors dec_year and temperature_air must be the same length.');
end

%% 2. Simulate Mean Sea Level Pressure (MSLP) Weather Patterns
% MSLP fluctuates around 1013.25 hPa due to synoptic weather (storms/highs).

% A. Define Weather Variance based on Latitude
%    Tropics have stable pressure; Mid-latitudes (storm tracks) vary wildly.
%    Low Lat (0): ~300 Pa sigma. High Lat (60): ~1200 Pa sigma.
phi = deg2rad(latitude);
sigma_pressure = 300 + (900 .* sin(phi).^2); % Pa

% B. Generate Red Noise (AR1 Process)
%    Pressure systems move slower than daily temp fluctuations.
%    Correlation ~0.8 to 0.9 for daily data.
alpha = 0.85; 

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
T_column_avg = temperature_air + (0.5 * LapseRate * elevation);

% Physical Constants
g = 9.81; % Gravity (m/s^2)
R = 287.05;  % Specific gas constant for dry air (J/kg*K)

% Hypsometric Equation
% P_local = P_msl * exp( -g*z / (R*T) )
exponent = (-g * elevation) ./ (R .* T_column_avg);

pressure_air = P_msl .* exp(exponent);

end