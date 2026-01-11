function ClimateForcingStep = model_inputs_single_timestep(index, dt, ClimateForcing, ModelParam)
% model_inputs_single_timestep
%
%% Syntax
%
%
%
%% Description
%
%
%
%% Example
% 
%
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

ClimateForcingStep.dt      = dt;                              % time step in seconds
ClimateForcingStep.T_air   = ClimateForcing.T_air0(index);    % screen level air temperature [K]   
ClimateForcingStep.V       = ClimateForcing.V0(index);        % wind speed [m s-1]
ClimateForcingStep.dlw     = ClimateForcing.dlw0(index);      % downward longwave radiation flux [W m-2]
ClimateForcingStep.dsw     = ClimateForcing.dsw0(index);      % downward shortwave radiation flux [W m-2]
ClimateForcingStep.e_air   = ClimateForcing.e_air0(index);    % screen level vapor pressure [Pa]
ClimateForcingStep.p_air   = ClimateForcing.p_air0(index);    % screen level air pressure [Pa]
ClimateForcingStep.P       = ClimateForcing.P0(index);        % precipitation [kg m-2]


% Location specifc parameters
ClimateForcingStep.Vz = ClimateForcing.Vz;
ClimateForcingStep.Tz = ClimateForcing.Tz;
ClimateForcingStep.T_air_mean = ClimateForcing.T_air_mean;
ClimateForcingStep.V_mean = ClimateForcing.V_mean;
ClimateForcingStep.P_mean = ClimateForcing.P_mean;
ClimateForcingStep.elev = ClimateForcing.elev;
ClimateForcingStep.lat = ClimateForcing.lat;
ClimateForcingStep.lon = ClimateForcing.lon;

% if we are provided with cc and cot values, extract for the timestep
if numel(ModelParam.black_carbon_snow)>1
    ClimateForcingStep.black_carbon_snow = ModelParam.black_carbon_snow(index);
else
    ClimateForcingStep.black_carbon_snow = ModelParam.black_carbon_snow;
end

if numel(ModelParam.black_carbon_ice)>1
    ClimateForcingStep.black_carbon_ice = ModelParam.black_carbon_ice(index);
else
    ClimateForcingStep.black_carbon_ice = ModelParam.black_carbon_ice;
end

if numel(ModelParam.cloud_optical_thickness)>1
    ClimateForcingStep.cloud_optical_thickness = ModelParam.cloud_optical_thickness(index);
else
    ClimateForcingStep.cloud_optical_thickness = ModelParam.cloud_optical_thickness;
end

if numel(ModelParam.solar_zenith_angle)>1
    ClimateForcingStep.solar_zenith_angle = ModelParam.solar_zenith_angle(index);
else
    ClimateForcingStep.solar_zenith_angle = ModelParam.solar_zenith_angle;
end

if numel(ModelParam.dsw_diffuse)>1
    ClimateForcingStep.dsw_diffuse = ModelParam.dsw_diffuse(index);
else
    ClimateForcingStep.dsw_diffuse = ModelParam.dsw_diffuse;
end

if numel(ModelParam.cloud_fraction)>1
    ClimateForcingStep.cloud_fraction = ModelParam.cloud_fraction(index);
else
    ClimateForcingStep.cloud_fraction = ModelParam.cloud_fraction;
end