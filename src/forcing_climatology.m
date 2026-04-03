function ClimateForcingSpinup = forcing_climatology(ClimateForcing, datetime_range)
% forcing_climatology creates a climatological average forcing structure that can 
% be used to spinup a gemb run.
%
%% Syntax
% 
%  ClimateForcingSpinup = forcing_climatology(ClimateForcing)
%  ClimateForcingSpinup = forcing_climatology(ClimateForcing, datetime_range)
% 
%% Description 
% 
% ClimateForcingSpinup = forcing_climatology(ClimateForcing) returns the climatological
% mean forcing of ClimateForcing structure for GEMB spinup runs. 
%
% ClimateForcingSpinup = forcing_climatology(ClimateForcing, datetime_range) specifies 
% a range of dates to include in the climatology. The datetime_range must be 
% 1x2 datetime in the form [datetime_start datetime_end]. 
% 
%% Examples
% For examples, go to https://github.com/alex-s-gardner/GEMB
% 
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 
% 
% See also model_initialize_forcing and simulate_climate_forcing.

%% Input checks

assert(istimetable(ClimateForcing), "First input to forcing_climatology function must be an Mx7 timetable containing climate forcing.")

if nargin>1
    % Only use the specified date range:
    keep = ClimateForcing.time>=datetime_range(1) & ClimateForcing.time<=datetime_range(2);
    ClimateForcing = ClimateForcing(keep,:);
end

%% Calculate climatology

% Eliminate the 366th day of the year: 
leap_day = day(ClimateForcing.time, "dayofyear")==366;
ClimateForcing(leap_day,:) = [];

% Tabulate a histogram of how many timesteps correspond to each unique year:
tbl = tabulate(year(ClimateForcing.time)); 

% Find out which years have the maximum number of entries (so we can exclude partial years): 
table_index = tbl(:,2) == max(tbl(:,2));

% These are the years that will contribute to the climatology: 
unique_years = tbl(table_index,1);

% Get the ClimateForcing indices that will contribute to the climatology: 
forcing_index = ismember(year(ClimateForcing.time), unique_years);

% Reshape everything into 2D matrices and take the mean: 
time_2D            =      reshape(ClimateForcing.time(forcing_index),               [], numel(unique_years));
temperature_air    = mean(reshape(ClimateForcing.temperature_air(forcing_index),    [], numel(unique_years)),2);
pressure_air       = mean(reshape(ClimateForcing.pressure_air(forcing_index),       [], numel(unique_years)),2);
precipitation      = mean(reshape(ClimateForcing.precipitation(forcing_index),      [], numel(unique_years)),2);
wind_speed         = mean(reshape(ClimateForcing.wind_speed(forcing_index),         [], numel(unique_years)),2);
shortwave_downward = mean(reshape(ClimateForcing.shortwave_downward(forcing_index), [], numel(unique_years)),2);
longwave_downward  = mean(reshape(ClimateForcing.longwave_downward(forcing_index),  [], numel(unique_years)),2);
vapor_pressure     = mean(reshape(ClimateForcing.vapor_pressure(forcing_index),     [], numel(unique_years)),2);
time = time_2D(:,1);

%% Create output timetable

ClimateForcingSpinup = timetable(time, temperature_air, pressure_air, precipitation, ...
    wind_speed, shortwave_downward, longwave_downward, vapor_pressure); 

ClimateForcingSpinup.Properties.VariableUnits = ["K", "Pa", "kg m^-2", "m s^-1", "W m^-2", "W m^-2", "Pa"]; 

ClimateForcingSpinup = addprop(ClimateForcingSpinup,...
    ["temperature_air_mean", "wind_speed_mean", "precipitation_mean", "temperature_observation_height", "wind_observation_height"],...
    [            "table",              "table",              "table",                          "table",                   "table"]);

ClimateForcingSpinup.Properties.CustomProperties.temperature_air_mean           = ClimateForcing.Properties.CustomProperties.temperature_air_mean; 
ClimateForcingSpinup.Properties.CustomProperties.wind_speed_mean                = ClimateForcing.Properties.CustomProperties.wind_speed_mean; 
ClimateForcingSpinup.Properties.CustomProperties.precipitation_mean             = ClimateForcing.Properties.CustomProperties.precipitation_mean; 
ClimateForcingSpinup.Properties.CustomProperties.temperature_observation_height = ClimateForcing.Properties.CustomProperties.temperature_observation_height; 
ClimateForcingSpinup.Properties.CustomProperties.wind_observation_height        = ClimateForcing.Properties.CustomProperties.wind_observation_height; 

end