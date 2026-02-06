function shortwave_downward = simulate_shortwave_irradiance(decimal_year, latitude)
% simulate_shortwave_irradiance simulates clear sky shortwave irradiance.
%
%   shortwave_downward = shortwave_irradiance(decimal_year, latitude)
%
%   Inputs:
%       decimal_year : Decimal fractional year (e.g., 2024.5).
%                      Can be a scalar or a vector.
%                      Assumes Local Solar Time.
%       latitude     : Latitude in degrees (Positive for North).
%
%   Output:
%       shortwave_downward          : Downwelling Shortwave Irradiance (W/m^2).
%                      (Value is 0 when the sun is below the horizon).
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

% --- 1. Time Conversion (Decimal Year -> Day & Hour) ---
% Separate year integer to check for leap year
year_val = floor(decimal_year);
frac_year = decimal_year - year_val;

% Check for leap year (366 days)
is_leap = (mod(year_val, 4) == 0 & mod(year_val, 100) ~= 0) | ...
          (mod(year_val, 400) == 0);
days_in_year = 365 + is_leap;

% Continuous Day of Year (e.g., 1.5 = Jan 1st at Noon)
doy_continuous = frac_year .* days_in_year + 1;

% Integer Day Number (n) for Declination
n = floor(doy_continuous);

% Solar Hour (0 to 24)
solar_hour = (doy_continuous - n) * 24;

% --- 2. Solar Geometry ---
% Convert Latitude to Radians
phi = deg2rad(latitude);

% Solar Declination (delta)
% Angle of the sun relative to the equatorial plane
% Cooper's Equation: 23.45 * sin(360/365 * (284 + n))
delta = deg2rad(23.45 * sind((360 ./ 365.25) .* (284 + n)));

% Hour Angle (omega)
% 0 at Solar Noon, -ve morning, +ve afternoon. 15 deg per hour.
omega = deg2rad(15 * (solar_hour - 12));

% Cosine of Solar Zenith Angle (theta_z)
% cos(z) = sin(latitude)sin(delta) + cos(latitude)cos(delta)cos(omega)
cos_theta_z = (sin(phi) .* sin(delta)) + ...
              (cos(phi) .* cos(delta) .* cos(omega));

% --- 3. Calculate Irradiance (Haurwitz Model) ---
% DSW = 1098 * cos(z) * exp(-0.057 / cos(z))

% Initialize output vector
shortwave_downward = zeros(size(decimal_year));

% Find indices where sun is above horizon
daylight_mask = cos_theta_z > 0;

if any(daylight_mask)
    ctz = cos_theta_z(daylight_mask);
    % Apply model only to daylight hours
    shortwave_downward(daylight_mask) = 1098 .* ctz .* exp(-0.057 ./ ctz);
end

end