function rh = relative_humidity(e_air, T_air)
% relative_humidity calculates relative humidity from vapor pressure and temperature.
% using the formula: 
%
%  rh = (e_air / e_saturation(T_air)) * 100
% 
%% Syntax 
% 
%  rh = relative_humidity(e_air, T_air)
%
%% Description
% 
% rh = relative_humidity(e_air, T_air) returns percent relative humidity rh 
% (from 0 to 100) as a function of vapor pressure e_air (Pa) and
% near-surface temperature T_air (K). 
%
%% Example 
% If the vapor pressure 313.9 Pa and temperature is 265.3 K, 
% then the relative humidity is: 
% 
%   rh = relative_humidity(313.9, 265.3)
%      = 92.7913
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

arguments
    e_air (:,:) {mustBeNumeric,mustBeGreaterThanOrEqual(e_air,0)}
    T_air (:,:) {mustBeNumeric,mustBeGreaterThan(T_air,0)}
end

% 1. Calculate Saturated Vapor Pressure (es) in Pa
%    Using Tetens' Formula: es = 610.78 * exp(17.27 * Tc / (Tc + 237.3))

Tc = T_air - 273.15; % Convert Kelvin to Celsius

A = 610.78; % Pa
B = 17.27;
C = 237.3;

es = A .* exp((B .* Tc) ./ (Tc + C));

% 2. Calculate Relative Humidity
rh = (e_air ./ es) * 100;

% 3. Clamp logical bounds (0% to 100%)
%    Physically, RH > 100% is supersaturation (rain/fog), but for
%    general variables, we typically cap it at 100.
rh = max(0, min(100, rh));

end