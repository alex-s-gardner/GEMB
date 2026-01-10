function dlw = longwave_irradiance(temp_c, rh_percent)
% longwave_irradiance Simulates downwelling longwave irradiance.
%
%   lw_down = calc_longwave_t_rh(temp_c, rh_percent)
%
%   Inputs:
%       temp_c     : 2-meter Atmospheric Temperature in degrees Celsius.
%                    (Scalar or Vector)
%       rh_percent : Relative Humidity in percent (0-100).
%                    (Scalar or Vector)
%
%   Output:
%       lw_down    : Downwelling Longwave Irradiance in W/m^2.
%
%   Model:
%       Uses the Brutsaert (1975) clear sky emissivity model.
%       1. Calculates Saturation Vapor Pressure (Tetens Formula).
%       2. Calculates Actual Vapor Pressure from RH.
%       3. Estimates Emissivity = 1.24 * (e_a / T_K)^(1/7).
%       4. applies Stefan-Boltzmann Law.
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

% --- Constants ---
sigma = 5.67e-8;  % Stefan-Boltzmann constant (W m^-2 K^-4)

% --- Unit Conversions ---
temp_k = temp_c + 273.15; % Celsius to Kelvin

% --- 1. Calculate Vapor Pressure ---
% Saturation Vapor Pressure (e_s) in hPa using Tetens Formula
e_s = 6.112 .* exp((17.67 .* temp_c) ./ (temp_c + 243.5));

% Actual Vapor Pressure (e_a) in hPa
e_a = e_s .* (rh_percent ./ 100);

% --- 2. Calculate Emissivity ---
% Brutsaert (1975) Equation: epsilon = 1.24 * (e_a / T_K)^(1/7)
% Note: e_a is in hPa (mbar) and temp_k is in Kelvin.
epsilon_sky = 1.24 .* (e_a ./ temp_k).^(1/7);

% --- 3. Calculate Longwave Irradiance ---
% Stefan-Boltzmann Law: L = epsilon * sigma * T^4
dlw = epsilon_sky .* sigma .* (temp_k.^4);

end