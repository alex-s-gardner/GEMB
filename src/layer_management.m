function [T, dz, d, W, re, gdn, gsp, a, a_diffuse, M_added, E_added] = ...
    layer_management(T, dz, d, W, re, gdn, gsp, a, a_diffuse, ModelParam, verbose)
% layer_management adjusts the depth and number of vertical layers in the model
% to ensure that the thickness of any single layer does not exceed thresholds
% set for the minimum and maximum allowable layer thickness.
%
%% Syntax
%
% [T, dz, d, W, re, gdn, gsp, a, a_diffuse, M_added, E_added] = ...
%    layer_management(T, dz, d, W, re, gdn, gsp, a, a_diffuse, ModelParam, verbose)
%
%% Description
%
% This function manages the vertical grid discretization of the firn column. 
% It performs three main operations to maintain numerical stability and 
% physical realism:
%
% 1. Merging: Scans the column for grid cells thinner than the minimum 
%    threshold (ModelParam.column_dzmin). These cells are merged with 
%    neighbors, conserving mass and energy via weighted averaging.
% 2. Splitting: Scans for grid cells thicker than the maximum threshold 
%    (ModelParam.column_dzmax). These cells are split in half, duplicating 
%    intensive properties and halving extensive properties.
% 3. Depth Adjustment: Ensures the total column depth stays within defined 
%    limits (ModelParam.column_zmax). It adds or removes layers at the 
%    bottom boundary as necessary and enforces the Dirichlet temperature 
%    boundary condition.
%
% The function tracks any mass (M_added) or energy (E_added) introduced or 
% removed during the bottom boundary adjustment to ensure closure of the 
% mass and energy balance budgets.
%
%% Inputs
%
%  T                       : K            Grid cell temperature.
%  dz                      : m            Grid cell thickness.
%  d                       : kg m^-3      Grid cell density.
%  W                       : kg m^-2      Water content.
%  re                      : mm           Grain size (effective radius).
%  gdn                     : unitless     Grain dendricity.
%  gsp                     : unitless     Grain sphericity.
%  a                       : fraction     Albedo.
%  a_diffuse               : fraction     Diffuse albedo.
%  ModelParam              : struct       Structure containing model parameters:
%    .column_dzmin         : m            Minimum allowable grid spacing.
%    .column_dzmax         : m            Maximum allowable grid spacing.
%    .column_zmax          : m            Maximum depth of the total column.
%    .column_zmin          : m            Minimum depth of the total column.
%    .column_ztop          : m            Thickness of the upper portion of the grid with constant spacing.
%    .column_zy            : unitless     Grid stretching parameter for the lower portion.
%  verbose                 : logical      Flag to enable mass/energy conservation checks.
%
%% Outputs
%
%  T            : K            Updated grid cell temperature.
%  dz           : m            Updated grid cell thickness.
%  d            : kg m^-3      Updated grid cell density.
%  W            : kg m^-2      Updated water content.
%  re           : mm           Updated grain size.
%  gdn          : unitless     Updated grain dendricity.
%  gsp          : unitless     Updated grain sphericity.
%  a            : fraction     Updated albedo.
%  a_diffuse    : fraction     Updated diffuse albedo.
%  M_added      : kg m^-2      Mass added to (positive) or removed from (negative) the column bottom.
%  E_added      : J m^-2       Energy added to (positive) or removed from (negative) the column bottom.
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

d_tolerance  = 1e-11; % tolerance for numerical comparison.

% Specify constants:
CtoK = 273.15;   % Celsius to Kelvin conversion
CI   = 2102;     % specific heat capacity of snow/ice (J kg-1 K-1)
LF   = 0.3345E6; % latent heat of fusion (J kg-1)

% store initial mass [kg] and energy [J]
M        = dz .* d;               % grid cell mass [kg]
M_total_initial = sum(W) + sum(M);       % total mass [kg]
E_total_initial = sum(M .* T * CI) + ...
    sum(W .* (LF + CtoK * CI));   % total energy [J] = initial enegy of snow/ice + initial enegy of water


T_bottom = T(end);
m        = length(T);

z_cumulative = cumsum(dz);

% A logical "mask" that indicates which cells are in the top layers:
top_layers = z_cumulative <= (ModelParam.column_ztop + d_tolerance);

% Define column_dzmin2 array using the top-layers' ModelParam.column_dzmin value for the entire column:
column_dzmin2 = ModelParam.column_dzmin * ones(m,1);

% Overwrite the bottom layers as the cumulative product times the stretching factor:
column_dzmin2(~top_layers) = cumprod(ModelParam.column_zy * ones(sum(~top_layers),1)) * ModelParam.column_dzmin;

% Define column_dzmax2 array using the top-layers' ModelParam.column_dzmin value for the entire column:
column_dzmax2 = ModelParam.column_dzmax * ones(m,1);

% Overwrite the bottom layers as the cumulative product times the stretching factor:
column_dzmax2(~top_layers) = cumprod(ModelParam.column_zy * ones(sum(~top_layers),1)) * ModelParam.column_dzmax;

% Preallocate a logical array that will be true for any cell to be deleted:
delete_cell = false(m,1);

% Check to see if any cells are too small and need to be merged
for i=1:m

    if dz(i) < (column_dzmin2(i) - d_tolerance)

        % dz has not met minimum thickness requirements, so we will delete it
        % and merge its contents into another cell:
        delete_cell(i) = true;

        % Detemine the target location for the cell contents to go:
        if i == m
            % If the very bottom cell (i==m) is too small, find the lowermost
            % cell that isn't going to be deleted:
            i_target = find(~delete_cell,1,'last');
        else
            i_target = i + 1;
        end

        % Move the quantities to the target location. Quantities are
        % calculated as linearly weighted functions of mass:
        m_new               = M(i) + M(i_target);
        T(i_target)         = (T(i)         * M(i) + T(i_target)         * M(i_target)) / m_new;
        a(i_target)         = (a(i)         * M(i) + a(i_target)         * M(i_target)) / m_new;
        a_diffuse(i_target) = (a_diffuse(i) * M(i) + a_diffuse(i_target) * M(i_target)) / m_new;

        % Use grain properties from lower cell:
        re(i_target)  = re(i);
        gdn(i_target) = gdn(i);
        gsp(i_target) = gsp(i);

        % Merge with underlying grid cell and delete old cell:
        dz(i_target) = dz(i) + dz(i_target);         % combine cell depths
        d(i_target)  = m_new / dz(i_target);         % combine top densities
        W(i_target)  = W(i) + W(i_target);           % combine liquid water
        M(i_target)  = m_new;                        % combine top masses
    end
end

% Delete combined cells:
W(delete_cell)                 = [];
dz(delete_cell)                = [];
d(delete_cell)                 = [];
T(delete_cell)                 = [];
a(delete_cell)                 = [];
re(delete_cell)                = [];
gdn(delete_cell)               = [];
gsp(delete_cell)               = [];
a_diffuse(delete_cell)         = [];
column_dzmax2(delete_cell)     = [];

% Calculate *new* length of cells:
m = length(T);

%% Split cells
% * An early implementation of this code used a loop which is included in comments at the bottom of this function for posterity.

% Find the cells that exceed tolerances:
f = find(dz > (column_dzmax2 + d_tolerance));

% Conserve quantities among the cells that will be split:
dz(f) = dz(f)/2;
W(f)  = W(f) /2;

% Sort the indices of all the cells including the ones that will be duplicated:
fs = sort([(1:m)';f]);

% Recreate the variables with split cells:
dz        = dz(fs);
W         = W(fs);
T         = T(fs);
d         = d(fs);
a         = a(fs);
a_diffuse = a_diffuse(fs);
re        = re(fs);
gdn       = gdn(fs);
gsp       = gsp(fs);

%% CORRECT FOR TOTAL MODEL DEPTH
% WORKS FINE BUT HAS BEEN DISABLED FOR CONVIENCE OF MODEL OUTPUT
% INTERPRETATION

% Calculate total model depth:
z_total = sum(dz);

if z_total < (ModelParam.column_zmax - d_tolerance)

    % Mass and energy to be added:
    M_added   = (dz(end) * d(end)) + W(end);
    E_added   = T(end) * (dz(end) * d(end)) * CI + W(end) * (LF + CtoK * CI);

    % Add a grid cell of the same size and temperature to the bottom:
    dz        = [   dz;    dz(end)];
    T         = [    T;     T(end)];
    W         = [    W;     W(end)];
    d         = [    d;     d(end)];
    a         = [    a;     a(end)];
    a_diffuse = [a_diffuse; a_diffuse(end)];
    re        = [   re;    re(end)];
    gdn       = [  gdn;   gdn(end)];
    gsp       = [  gsp;   gsp(end)];

elseif z_total > ModelParam.column_zmax+d_tolerance 

    % Mass and energy loss:
    M_added   = -((dz(end)*d(end)) + W(end));
    E_added   = -(T(end) * (dz(end)*d(end)) * CI) - W(end) * (LF+CtoK*CI);

    % Remove a grid cell from the bottom:
    dz(end)        = [];
    T(end)         = [];
    W(end)         = [];
    d(end)         = [];
    a(end)         = [];
    re(end)        = [];
    gdn(end)       = [];
    gsp(end)       = [];
    a_diffuse(end) = [];
else
    % No mass or energy is added or removed:
    M_added   = 0;
    E_added   = 0;
end

% The temperature of the bottom grid cell may have been modified when
% layers were merged. After checking for energy conservation, set:
%       T(end) = T_bottom
% This is to satisfy the Constant Temperature (Dirichlet) boundary
% condition. If this is not done then then thermal diffusion will blow up
E_added   = E_added + ((T_bottom - T(end)) * (dz(end)*d(end)) * CI);
T(end)    = T_bottom;

%% CHECK FOR MASS AND ENERGY CONSERVATION
if verbose
    % Calculate final mass [kg] and energy [J]
    M             = dz .* d;               % grid cell mass [kg]
    M_total_final = sum(W) + sum(M);       % total mass [kg]
    E_total_final = sum(M .* T * CI) + ...
        sum(W .* (LF + CtoK * CI));        % total energy [J] = initial enegy of snow/ice + initial enegy of water

    M_delta = M_total_initial - M_total_final + M_added;
    E_delta = E_total_initial - E_total_final + E_added;

    if (abs(M_delta) > 1E-3) || (abs(E_delta) > 1E-3)
        error(['Mass and/or energy are not conserved in melt equations:' newline ' M_delta: ' ...
            num2str(M_delta) ' E_delta: ' num2str(E_delta) newline])
    end

end

end