function longname = varname2longname(varname)
% varname2longname returns a short description of given variable names.
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

% Define Superscript Characters
sup_minus = char(8315); % ⁻
sup_1     = char(185);  % ¹
sup_2     = char(178);  % ²

% Construct Unit Strings
unit_wm2  = " [W m" + sup_minus + sup_2 + "]";  % [W m⁻²]
unit_kgm2 = " [kg m" + sup_minus + sup_2 + "]"; % [kg m⁻²]
unit_ms1  = " [m s" + sup_minus + sup_1 + "]";  % [m s⁻¹]

keys = [
    "longwave_downward", ...
    "shortwave_downward", ...
    "relative_humidity", ...
    "temperature_air", ...
    "precipitation", ...
    "wind_speed", ...
    "vapor_pressure", ...
    "pressure_air"
];

values = [
    "downward logwave radiation" + unit_wm2, ...
    "downward shortwave radiation" + unit_wm2, ...
    "screen level relative humidity [%]", ...
    "screen level air temperature [K]", ...
    "precipitation" + unit_kgm2, ...
    "screen level wind speed" + unit_ms1, ...
    "vapor pressure [Pa]", ...
    "screen level air pressure [Pa]"
];

% Create the dictionary
var_defs = dictionary(keys, values);

longname = var_defs(varname);
end