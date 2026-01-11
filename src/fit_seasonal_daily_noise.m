function coeffs = fit_seasonal_daily_noise(dec_year, y_data)
% fit_seasonal_daily_noise fits yearly/daily sinusoids and noise stats.
%
%% Syntax
% 
%  coeffs = fit_seasonal_daily_noise(dec_year, y_data)
%
%% Description
%
% coeffs = fit_seasonal_daily_noise(dec_year, y_data) determines
% coefficients of a model that represents y_data as the linear sum of
% daily and annual sinusoids, plus noise and a constant offset.
%
%% Details 
%
%   Inputs:
%       frac_year : Vector of time points in decimal years.
%       y_data    : Vector of observed data.
%
%   Output:
%       coeffs    : Structure containing fitted parameters:
%                   .beta (5x1 vector of linear weights)
%                   .noise_std (Standard Deviation of residuals)
%                   .noise_lag1 (Lag-1 Autocorrelation of residuals)
%
%   Model:
%       y = c1 + c2*cos(Yr) + c3*sin(Yr) + c4*cos(Day) + c5*sin(Day) + noise
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
% See also: simulate_seasonal_daily_noise. 

% Ensure inputs are column vectors
t = dec_year(:);
y = y_data(:);

% --- 1. Construct the Design Matrix (X) ---
% Frequency Constants
omega_yr  = 2 * pi;           % Once per year
omega_day = 2 * pi * 365.25;  % Approx 365.25 times per year

% Basis Functions
% Col 1: Mean (Intercept)
% Col 2: Daily Cosine
% Col 3: Daily Sine
% Col 4: Yearly Cosine
% Col 5: Yearly Sine
X = [ones(size(t)), ...
     cos(omega_day * t), sin(omega_day * t),...
     cos(omega_yr * t), sin(omega_yr * t)];

% --- 2. Linear Regression (Least Squares) ---
% Solve y = X*beta
beta = X \ y;

% Store deterministic coefficients
coeffs.beta = beta;

% --- 3. Residual Analysis (Noise Fitting) ---
% Calculate the deterministic model prediction
y_model = X * beta;

% Calculate residuals (What's left over)
residuals = y - y_model;

% Calculate Noise Statistics
coeffs.noise_std = std(residuals);

% Calculate Lag-1 Autocorrelation
% Note: This assumes reasonably uniform sampling.
if length(residuals) > 1
    % Correlation between x[t] and x[t-1]
    c_matrix = corrcoef(residuals(1:end-1), residuals(2:end));
    coeffs.noise_lag1 = c_matrix(1, 2);
else
    coeffs.noise_lag1 = 0;
end

end