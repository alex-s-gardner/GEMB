function vapor_pressure = dewpoint_to_vapor_pressure(temperature_dewpoint)
% dewpoint_to_vapor_pressure converts dewpoint temperature to actual vapor pressure. 
% 
%% Syntax 
%
% vapor_pressure = dewpoint_to_vapor_pressure(temperature_dewpoint)
%
%% Description 
%
% vapor_pressure = dewpoint_to_vapor_pressure(temperature_dewpoint)
% converts temperature_dewpoint (K) to actual vapor pressure in Pa. 
%
%% Reference
% Formula by Tim Brice and Todd Hall of NOAA. 
% https://www.weather.gov/epz/wxcalc_vaporpressure
% 
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

%% Input checks 

arguments
    temperature_dewpoint  (:,:) {mustBeNumeric}
end

if temperature_dewpoint<100
    warning('Unrealistic dewpoint temperature. Ensure its units are kelvin.')
end

%% Mathematics

Td = temperature_dewpoint - 273.15;

vapor_pressure = 611.*10.^( 7.5*Td./ (237.3 + Td) );

end

