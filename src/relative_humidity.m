function rh = relative_humidity(eAir, Ta)
% RELATIVE_HUMIDITY Calculates relative humidity from vapor pressure and temp.
%
%   rh = relative_humidity(eAir, Ta)
%
%   INPUTS:
%   eAir - Actual Vapor Pressure in Pascals [Pa]
%   Ta   - Near-surface air temperature in Kelvin [K]
%
%   OUTPUT:
%   rh   - Relative Humidity in Percent [%] (clamped 0 to 100)
%
%   FORMULA:
%   RH = (eAir / e_saturation(Ta)) * 100

    % 1. Calculate Saturated Vapor Pressure (es) in Pa
    %    Using Tetens' Formula: es = 610.78 * exp(17.27 * Tc / (Tc + 237.3))
    
    Tc = Ta - 273.15; % Convert Kelvin to Celsius
    
    A = 610.78; % Pa
    B = 17.27;
    C = 237.3;
    
    es = A .* exp((B .* Tc) ./ (Tc + C));

    % 2. Calculate Relative Humidity
    rh = (eAir ./ es) * 100;

    % 3. Clamp logical bounds (0% to 100%)
    %    Physically, RH > 100% is supersaturation (rain/fog), but for
    %    general variables, we typically cap it at 100.
    rh = max(0, min(100, rh));

end