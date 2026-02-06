function y_sim = simulate_seasonal_daily_noise(dec_year, coeffs)
% simulate_seasonal_daily_noise generates data from fitted coefficients.
%
%   y_sim = seasonal_daily_noise(frac_year, coeffs)
%
%   Inputs:
%       frac_year : Vector of time points to simulate.
%       coeffs    : Structure from seasonal_daily_noise_fit containing:
%                   .beta, .noise_std, .noise_lag1
%
%   Output:
%       y_sim     : Simulated y values (Deterministic + Random Noise).
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
% See also fit_seasonal_daily_noise.

% Ensure input is column vector
t = dec_year(:);
n = length(t);

% --- 1. Reconstruct Deterministic Component ---
omega_yr  = 2 * pi;
omega_day = 2 * pi * 365.25;

% Rebuild Design Matrix
X = [ones(size(t)), ...=
     cos(omega_day * t), sin(omega_day * t), ...
     cos(omega_yr  * t), sin(omega_yr  * t)];

% Calculate pure signal
y_deterministic = X * coeffs.beta;

% --- 2. Generate Synthetic Correlated Noise ---
% Uses AR(1) process: x(t) = phi*x(t-1) + epsilon

sigma = coeffs.noise_std;
phi   = coeffs.noise_lag1;

% Calculate scaling for the driving white noise to match target std
% Var(process) = Var(noise) / (1 - phi^2) -> Var(noise) = Var(process)*(1-phi^2)
white_noise_scale = sigma * sqrt(1 - phi^2);

% Generate white noise
u = randn(n, 1) * white_noise_scale;

% Correct initialization transient if using 'filter' 
% (Filter assumes 0 initial state, but we manually set noise(1))
% Alternatively, simple loop for clarity (slower but robust for simulation):
noise_loop = zeros(n, 1);
noise_loop(1) = randn * sigma;
for i = 2:n
    noise_loop(i) = phi * noise_loop(i-1) + u(i);
end

% --- 3. Combine ---
y_sim = y_deterministic + noise_loop;

end