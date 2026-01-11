function [dz, z_center] = model_initialize_grid(ModelParam)
% model_initialize_grid sets up the initial layer thickness and total grid depth.  
% 
%% Syntax 
% 
%  dz = model_initialize_grid(ModelParam)
%  [dz,z_center] = model_initialize_grid(ModelParam)
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
% [dz,z_center] =model_initialize_grid(ModelParam) also returns a
% 1D array z_center containing the center depths of each grid cell in
% meters.
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

% Optional output: 
if nargout>1
    z_center = -cumsum(dz) + dz/2; 
end

%% ---------NEED TO IMPLEMENT A PROPER GRID STRECHING ALGORITHM------------
% See https://github.com/alex-s-gardner/GEMB/issues/8. 

end