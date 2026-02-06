function [temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse] = model_initialize_column(ModelParam, ClimateForcing)
% model_initialize_column initializes a GEMB column based on specified
% model and climate forcing parameters.
%
%% Syntax 
%
%  [temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse] = model_initialize_column(ModelParam, ClimateForcing)
%
%% Description
%
% [temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse] = model_initialize_column(ModelParam, ClimateForcing)
% uses inputs ModelParam from model_initialize_parameters and input structure
% ClimateForcing containing the field temperature_air_mean. Outputs are as follows: 
% 
%    Variable          Units     Description                Initial Value
%    temperature       K         temperature                mean annual surface temperature 
%    dz                m         thickness                  array generated from model parameters
%    density           kg m^-3   ice density                ModelParam.density_ice
%    water             kg m^-2   old-snow water content     0
%    grain_radius      mm        old-snow grain size        2.5
%    grain_dendricity  fraction  old-snow grain dendricity  0
%    grain_sphericity  fraction  old-snow grain sphericity  0
%    albedo            fraction  surface albedo             ModelParam.albedo_snow
%    albedo_diffuse    fraction  diffuse surface albedo     ModelParam.albedo_snow
%   
%% Example
% Initialize a GEMB column based on default model parameters and a mean
% surface temperature of -20 C. 
% 
%   % Initialize parameters: 
%   ModelParam = model_initialize_parameters;
%   ClimateForcing.temperature_air_mean = 253.15; % -20 C
%  
%   % Initialize Column: 
%   [temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse] = model_initialize_column(ModelParam, ClimateForcing);
%
%   % Inspect outputs: 
%   whos
%   Name                    Size            Bytes  Class     Attributes
% 
%   ClimateForcing          1x1               176  struct              
%   ModelParam              1x1              7972  struct              
%   albedo                264x1              2112  double              
%   albedo_diffuse        264x1              2112  double              
%   density               264x1              2112  double              
%   dz                    264x1              2112  double              
%   grain_dendricity      264x1              2112  double              
%   grain_radius          264x1              2112  double              
%   grain_sphericity      264x1              2112  double              
%   temperature           264x1              2112  double              
%   water                 264x1              2112  double              
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

assert(isscalar(ClimateForcing.temperature_air_mean),'Temperature ClimateForcing.temperature_air_mean must be a scalar.')
assert(ClimateForcing.temperature_air_mean>0,'Temperature ClimateForcing.temperature_air_mean must exceed 0 K.')
if ClimateForcing.temperature_air_mean<100
    warning('Temperature ClimateForcing.temperature_air_mean should be in kelvin, but is below 100, suggesting an error. Confirm that the units are kelvin.')
end

% initialze column variables 
dz               = model_initialize_grid(ModelParam);
m                = length(dz);
temperature      = zeros(m,1) + ClimateForcing.temperature_air_mean; % initial grid cell temperature to the annual mean temperature [K]
density          = zeros(m,1) + ModelParam.density_ice;    % density to that of ice [kg m-3]
water            = zeros(m,1);                             % water content of zero [kg m-2]
grain_radius     = zeros(m,1) + 2.5;                       % grain size of old snow [mm]
grain_dendricity = zeros(m,1);                             % grain dentricity of old snow
grain_sphericity = zeros(m,1);                             % grain sphericity of old snow
albedo           = zeros(m,1) + ModelParam.albedo_snow;    % albedo equal to fresh snow [fraction]
albedo_diffuse   = zeros(m,1) + ModelParam.albedo_snow;    % albedo equal to fresh snow [fraction]  

end

function dz = model_initialize_grid(ModelParam)
% model_initialize_grid sets up the initial layer thickness and total grid depth.  
% 
%% Syntax 
% 
%  dz = model_initialize_grid(ModelParam)
%
%% Description
% 
% dz = model_initialize_grid(ModelParam) creates a 1D grid structure
% containing the thickness dz of each cell in the column, where inputs match 
% Fig. 1 of Gardner et al., 2023 (https://doi.org/10.5194/gmd-16-2277-2023)
% and all inputs are scalars as follows: 
% 
%  * ModelParam.column_ztop (m): Thickness of the upper portion of the model grid, in which grid spacing is constant.
%  * ModelParam.column_ztop (m): Spacing of the upper portion of the model grid. 
%  * ModelParam.column_zmax (m): Maximum thickness of the total column. 
%  * ModelParam.column_zy (unitless): Grid cell stretching parameter for the lower portion of the model grid, in which grid length increases linearly with depth. 
% 
%% Example 
% 
% % Define inputs: 
% ModelParam.column_ztop  = 10; 
% ModelParam.column_dztop = 0.05; 
% ModelParam.column_zmax  = 250; 
% ModelParam.column_zy    = 1.025; 
% 
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

%% Error checks: 

% Calculate number of top grid points:
n_top = ModelParam.column_ztop/ModelParam.column_dztop;

d_tolerance  = 1e-11; % Depth tolerance 

% Check to see if the top grid cell structure length (ModelParam.column_dztop) goes evenly 
% into specified top structure depth (ModelParam.column_ztop)
assert(mod(n_top,1)==0,['Top grid cell structure length does not go evenly into ' ...
        'specified top structure depth, adjust ModelParam.column_dztop or ModelParam.column_ztop.'])

% Make sure top grid cell structure length (ModelParam.column_dztop) is greater than 5 cm
if ModelParam.column_dztop < 0.05-d_tolerance 
    warning('Initial top grid cell length (ModelParam.column_dztop) is < 0.05 m.')
end

%% Generate grid: 

% Initialize top grid depth vector:
dzT = ones(n_top,1)*ModelParam.column_dztop;

% Bottom grid cell depth = x*ModelParam.column_zy^(cells from to structure)

% Initialize bottom vectors
dzB = zeros(((ModelParam.column_zmax - ModelParam.column_ztop)/ModelParam.column_dztop),1);
gp0 = ModelParam.column_dztop;
z0  = ModelParam.column_ztop;
k   = 1;

while ModelParam.column_zmax > (z0 + d_tolerance) 
    dzB(k,1) = gp0 * ModelParam.column_zy;
    gp0 = dzB(k,1);
    z0 = z0 + gp0;
    k = k + 1;
end

% Delete excess cells from bottom vector: 
dzB(dzB == 0) = [];

% Combine top and bottom dz vectors
dz = [dzT ; dzB];

end