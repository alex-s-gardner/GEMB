function e_air = simulate_vapor_pressure(T_air, rh)
% simulate_vapor_pressure estimates actual vapor pressure from temperature and relative humidity.
%
%   e_air = simulate_vapor_pressure(T_air, rh)
%
%   INPUTS:
%   T_air   - Near-surface air temperature in Kelvin [K] (Scalar or Vector)
%   rh   - Relative Humidity in Percent [%] (0-100) (Scalar or Vector)
%
%   OUTPUT:
%   e_air - Actual Vapor Pressure in Pascals [Pa]
%
%   METHOD:
%   1. Calculate saturation vapor pressure (es) using Tetens' formula.
%   2. Apply humidity: e_air = es * (rh / 100)
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

% 1. Convert Kelvin to Celsius (required for Tetens' coefficients)
Tc = T_air - 273.15;

% 2. Calculate Saturation Vapor Pressure (es)
%    Constants for saturation over liquid water (Buck/Tetens)
A = 610.78;  % Pascals (Pressure at 0 deg C)
B = 17.27;   % Dimensionless
C = 237.3;   % Degrees Celsius

es = A .* exp((B .* Tc) ./ (Tc + C));

% 3. Calculate Actual Vapor Pressure
e_air = es .* (rh ./ 100);

end