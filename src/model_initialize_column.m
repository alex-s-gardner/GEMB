function [T, dz, d, W, re, gdn, gsp, a, a_diffuse] = model_initialize_column(ModelParam, ClimateForcing)
% model_initialize_column
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

% initialze column variables 
dz        = model_initialize_grid(ModelParam);
m         = length(dz);
T         = zeros(m,1) + ClimateForcing.T_air_mean; % initial grid cell temperature to the annual mean temperature [K]
d         = zeros(m,1) + ModelParam.density_ice;    % density to that of ice [kg m-3]
W         = zeros(m,1);                             % water content to zero [kg m-2]
re        = zeros(m,1) + 2.5;                       % grain size to old snow [mm]
gdn       = zeros(m,1);                             % grain dentricity to old snow
gsp       = zeros(m,1);                             % grain sphericity to old snow
a         = zeros(m,1) + ModelParam.albedo_snow;    % albedo equal to fresh snow [fraction]
a_diffuse = zeros(m,1) + ModelParam.albedo_snow;    % albedo equal to fresh snow [fraction]  
end