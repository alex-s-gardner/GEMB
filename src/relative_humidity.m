function rh = relative_humidity(e_air, T_air)
% RELATIVE_HUMIDITY Calculates relative humidity from vapor pressure and temp.
%
%   rh = relative_humidity(e_air, T_air)
%
%   INPUTS:
%   e_air - Actual Vapor Pressure in Pascals [Pa]
%   T_air   - Near-surface air temperature in Kelvin [K]
%
%   OUTPUT:
%   rh   - Relative Humidity in Percent [%] (clamped 0 to 100)
%
%   FORMULA:
%   RH = (e_air / e_saturation(T_air)) * 100
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

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