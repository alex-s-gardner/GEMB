function vapor_pressure = relative_humidity_to_vapor_pressure(temperature_air, relative_humidity)
% relative_humidity_to_vapor_pressure estimates actual vapor pressure from temperature and relative humidity.
%
%% Syntax
% 
%  vapor_pressure = relative_humidity_to_vapor_pressure(temperature_air, relative_humidity)
%
%% Description
% 
% vapor_pressure = relative_humidity_to_vapor_pressure(temperature_air, relative_humidity)
% 
%% Example 
% 
%
%% 
%   INPUTS:
%   temperature_air   - Near-surface air temperature in Kelvin [K] (Scalar or Vector)
%   relative_humidity - Relative Humidity in Percent [%] (0-100) (Scalar or Vector)
%
%   OUTPUT:
%   vapor_pressure    - Actual Vapor Pressure in Pascals [Pa]
%
%   METHOD:
%   1. Calculate saturation vapor pressure (es) using Tetens' formula.
%   2. Apply humidity: vapor_pressure = es * (relative_humidity / 100)
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
Tc = temperature_air - 273.15;

% 2. Calculate Saturation Vapor Pressure (es)
%    Constants for saturation over liquid water (Buck/Tetens)
A = 610.78;  % pascals (Pressure at 0 deg C)
B = 17.27;   % dimensionless
C = 237.3;   % degrees celsius

es = A .* exp((B .* Tc) ./ (Tc + C));

% 3. Calculate Actual Vapor Pressure
vapor_pressure = es .* (relative_humidity ./ 100);

end