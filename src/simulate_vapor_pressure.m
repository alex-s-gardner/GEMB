function e_air = simulate_vapor_pressure(T_air, rh)
% SIMULATE_VAPOR_PRESSURE Estimates actual vapor pressure from Temp and RH.
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