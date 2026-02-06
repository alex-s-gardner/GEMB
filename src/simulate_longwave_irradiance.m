function longwave_downward = simulate_longwave_irradiance(temperature_air, vapor_pressure)
% simulate_longwave_irradiance estimates downward longwave radiation (Vectorized).
%
%   longwave_downward = simulate_longwave_irradiance(temperature_air, vapor_pressure)
%
%   INPUTS:
%   temperature_air   - Near-surface air temperature in Kelvin [K] (Scalar or Vector)
%   vapor_pressure - Near-surface vapor pressure in Pascals [Pa] (Scalar or Vector)
%
%   OUTPUT:
%   longwave_downward  - Downward Longwave Radiation in [W/m^2] (Column Vector)
%
%   METHOD:
%   Uses Brutsaert's (1975) parameterization for effective clear-sky emissivity.
%   Formula: L_down = epsilon * sigma * temperature_air^4
%   Where:   epsilon = 1.24 * (e_hPa / temperature_air)^(1/7)
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

% Force inputs to column vectors if they are arrays
if ~isscalar(temperature_air)
    temperature_air = temperature_air(:);
end

if ~isscalar(vapor_pressure)
    vapor_pressure = vapor_pressure(:);
end

% Check for dimension mismatch (unless one is a scalar for broadcasting)
if ~isscalar(temperature_air) && ~isscalar(vapor_pressure) && (length(temperature_air) ~= length(vapor_pressure))
    error('Input vectors temperature_air and vapor_pressure must have the same length.');
end

%% 2. Physical Constants & Conversions

sigma = 5.670374419e-8;  % Stefan-Boltzmann constant [W m^-2 K^-4]

% Brutsaert's coefficient (1.24) requires vapor pressure in hPa (millibars).
% Convert Pa -> hPa
e_hPa = vapor_pressure ./ 100;

%% 3. Calculate Effective Emissivity

% Use element-wise power (.^) and division (./)
% Brutsaert (1975) for clear sky:
epsilon_clear = 1.24 .* (e_hPa ./ temperature_air).^(1/7);

%% 4. Calculate Downward Longwave Radiation

% Stefan-Boltzmann Law: L = epsilon * sigma * T^4
longwave_downward = epsilon_clear .* sigma .* (temperature_air.^4);

end