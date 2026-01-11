function coeffs = fit_air_temperature(dec_year, T_air, lat, elev)
% fit_air_temperature estimates simulation coefficients from observed data.
%
%% Syntax
% 
%  coeffs = fit_air_temperature(dec_year, T_air, lat, elev)
% 
%% Description
%
% coeffs = fit_air_temperature(dec_year, T_air, lat, elev)
%
%%
%   INPUTS:
%   dec_year - Vector of decimal years (e.g. 2024.0, 2024.002).
%   T_air       - Vector of observed air temperatures in Kelvin [K].
%   lat      - Latitude in degrees (scalar).
%   elev     - Elevation in meters (scalar).
%
%   OUTPUT:
%   coeffs   - Structure compatible with simulate_air_temperature.
%              Contains: .mean_offset, .lat_scale, .daily_amp_scale,
%                        .weather_sigma_scale, .weather_corr
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

dec_year = dec_year(:);
T_air    = T_air(:);

% Constants from the simulation model
BASE_SIGMA = 8.0; 
LAPSE_RATE = 0.0065;

%% 2. Calculate Mean Offset
% The model calculates a theoretical mean based on lat/elev.
% We compare the observed mean to this theoretical mean.

phi = deg2rad(lat); 
T_sea_level_theoretical = 300 - 50 .* sin(phi).^2;
T_mean_theoretical = T_sea_level_theoretical - (elev .* LAPSE_RATE);

T_obs_mean = mean(T_air, 'omitnan');

% RESULT 1: mean_offset
coeffs.mean_offset = T_obs_mean - T_mean_theoretical;

%% 3. Harmonic Analysis (Least Squares Fit)
% We remove the mean and fit the Annual and Diurnal cosines.

T_anomaly = T_air - T_obs_mean;

% A. Construct Annual Basis Function
% Match simulation logic: cos(2*pi * (year_frac - phase))
% Phase is 0.5 for North, 0 for South.
phase_annual = 0;
if lat > 0, phase_annual = 0.5; end

year_frac = dec_year - floor(dec_year);
basis_annual = cos(2*pi * (year_frac - phase_annual));

% B. Construct Diurnal Basis Function
% Match simulation logic: cos(2*pi * (day_fraction - 0.625))
day_fraction = mod(dec_year * 365.25, 1);
basis_diurnal = cos(2*pi * (day_fraction - 0.625));

% C. Solve Linear Regression: Anomaly = c1*Annual + c2*Diurnal
X = [basis_annual, basis_diurnal];
betas = X \ T_anomaly;

fitted_amp_annual = betas(1);
fitted_amp_daily  = betas(2);

%% 4. Derive Scale Coefficients
% We now solve the simulation equations for the scalar variables.

% --- lat_scale ---
% Model: Amp = 3 + (22 * |sin(lat)| * lat_scale)
% Inverse: lat_scale = (fitted_amp - 3) / (22 * |sin(lat)|)
denom_annual = 22 * abs(sin(phi));

if denom_annual < 1e-4
    % Handle Equator case (sin(0)=0) to avoid Inf. 
    % At equator, scale doesn't matter mathematically in this model.
    coeffs.lat_scale = 1.0; 
else
    coeffs.lat_scale = (fitted_amp_annual - 3) / denom_annual;
end
% Enforce bounds (scale cannot be negative)
coeffs.lat_scale = max(0, coeffs.lat_scale);

% --- daily_amp_scale ---
% Model: Amp = (DTR_base / 2) * daily_amp_scale
dtr_base = 10 + (elev / 1000);
base_amp_daily = dtr_base / 2;

coeffs.daily_amp_scale = fitted_amp_daily / base_amp_daily;
coeffs.daily_amp_scale = max(0, coeffs.daily_amp_scale);

%% 5. Estimate Weather Parameters (Residual Analysis)
% Remove the fitted deterministic cycles to get pure weather noise

T_deterministic = (fitted_amp_annual * basis_annual) + (fitted_amp_daily * basis_diurnal);
residuals = T_anomaly - T_deterministic;

% The simulation assumes weather is a DAILY process (AR1).
% We must average residuals by day to estimate these parameters accurately.
% Otherwise, high-freq sampling (e.g. hourly) would skew the correlation.

day_indices = floor(dec_year * 365.25);
[~, ~, idx] = unique(day_indices);

% Compute daily mean of residuals
daily_res = accumarray(idx, residuals, [], @mean);

% --- weather_sigma_scale ---
% Model: sigma = 8.0 * scale
std_daily = std(daily_res, 'omitnan');
coeffs.weather_sigma_scale = std_daily / BASE_SIGMA;

% --- weather_corr ---
% Calculate lag-1 autocorrelation
if length(daily_res) > 2
    % 'coeff' normalizes the covariance so lag 0 is 1.0
    c = xcov(daily_res, 1, 'coeff'); 
    % xcov output structure is [lag-1, lag0, lag+1]. We want lag+1 (index 3).
    raw_corr = c(3);
    
    % Clamp between 0 and 1
    coeffs.weather_corr = max(0, min(0.99, raw_corr));
else
    % Fallback for insufficient data
    coeffs.weather_corr = 0.7; 
end

end