function dlw = simulate_longwave_irradiance(T_air, e_air)
% SIMULATE_LONGWAVE_IRRADIANCE Estimates downward longwave radiation (Vectorized).
%
%   dlw = simulate_longwave_irradiance(T_air, e_air)
%
%   INPUTS:
%   T_air   - Near-surface air temperature in Kelvin [K] (Scalar or Vector)
%   e_air - Near-surface vapor pressure in Pascals [Pa] (Scalar or Vector)
%
%   OUTPUT:
%   dlw  - Downward Longwave Radiation in [W/m^2] (Column Vector)
%
%   METHOD:
%   Uses Brutsaert's (1975) parameterization for effective clear-sky emissivity.
%   Formula: L_down = epsilon * sigma * T_air^4
%   Where:   epsilon = 1.24 * (e_hPa / T_air)^(1/7)

    %% 1. Input Sanitization
    % Force inputs to column vectors if they are arrays
    if ~isscalar(T_air), T_air = T_air(:); end
    if ~isscalar(e_air), e_air = e_air(:); end
    
    % Check for dimension mismatch (unless one is a scalar for broadcasting)
    if ~isscalar(T_air) && ~isscalar(e_air) && (length(T_air) ~= length(e_air))
        error('Input vectors T_air and e_air must have the same length.');
    end

    %% 2. Physical Constants & Conversions
    sigma = 5.670374419e-8;  % Stefan-Boltzmann constant [W m^-2 K^-4]

    % Brutsaert's coefficient (1.24) requires vapor pressure in hPa (millibars).
    % Convert Pa -> hPa
    e_hPa = e_air ./ 100;

    %% 3. Calculate Effective Emissivity
    % Use element-wise power (.^) and division (./)
    % Brutsaert (1975) for clear sky:
    epsilon_clear = 1.24 .* (e_hPa ./ T_air).^(1/7);

    %% 4. Calculate Downward Longwave Radiation
    % Stefan-Boltzmann Law: L = epsilon * sigma * T^4
    dlw = epsilon_clear .* sigma .* (T_air.^4);

end