function [T, dz, d, W, re, gdn, gsp, a, a_diffuse, M_added, E_added] = ...
    layer_management(T, dz, d, W, re, gdn, gsp, a, a_diffuse, column_dzmin, ...
    column_zmax, column_zmin, column_ztop, column_zy, verbose)
% managelayers adjusts the depth and number of vertical layers in the model
% to ensure that the thickness of any single layer does not exceed thresholds
% set for the minimum and maximum allowable layer thickness.
%
%% Syntax
%
%
%
%% Description
%
%
%
%% Inputs
%
%  T            : K            Grid cell temperature.
%  d            : kg m^-3      Grid cell density.
%  dz           : m            Grid cell thickness.
%  W            : kg m^-2      Water content.
%  a            : fraction     Albedo.
%  a_diffuse    : fraction     Diffuse albedo.
%  M            : kg m^-2      Grid cell mass.
%  EI           : J m^-2       Initial energy of snow/ice.
%  EW           : J m^-2       Initial energy of water.
%  column_dzmin : m            Minimum allowable grid spacing.
%  column_zmax  : m            Maximum depth of the total column.
%  column_zmin  : m            Minimum depth of the total column.
%  re           : mm           Grain size
%  gdn          : unitless     Grain dendricity
%  gsp          : unitless     Grain sphericity
%  column_ztop  : m            Thickness of the upper portion of the model grid, in which grid spacing is constant.
%  column_zy    : unitless     Grid cell stretching parameter for the lower portion of the model grid, in which grid length increases linearly with depth.
%  CI           : J kg^-1 K^-1 Specific heat capacity of snow/ice.
%  LF           : J kg^-1      Latent heat of fusion.
%  CtoK         : K            273.15 conversion from C to K.
%
%% Outputs
%
%  d            : kg m^-3      Grid cell density.
%  T            : K            Grid cell temperature.
%  W            : kg m^-2      Water content.
%  M_added      : kg m^-2      Mass added to the column.
%  E_added      : J m^-2       Energy added to the column.
%  a            : fraction     Albedo.
%  a_diffuse    : fraction     Diffuse albedo.
%  m            : kg m^-2      Grid cell mass.
%  EI           : J m^-2       Initial energy of snow/ice.
%  EW           : J m^-2       Initial energy of water.
%  re           : mm           Grain size
%  gdn          : unitless     Grain dendricity
%  gsp          : unitless     Grain sphericity
%
%% Documentation
%
% For complete documentation, see: https://github.com/alex-s-gardner/GEMB
%
%% References
% If you use GEMB, please cite the following:
%
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass
% Balance (GEMB): a model of firn processes for cryosphere research, Geosci.
% Model Dev., 16, 2277-2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023

d_tolerance  = 1e-11; % tolerance for numerical comparison.

% Specify constants:
CtoK = 273.15;   % Celsius to Kelvin conversion
CI   = 2102;     % specific heat capacity of snow/ice (J kg-1 K-1)
LF   = 0.3345E6; % latent heat of fusion (J kg-1)

% store initial mass [kg] and energy [J]
M  = dz .* d;                  % grid cell mass [kg]
EI = M .* T * CI;              % initial enegy of snow/ice
EW = W .* (LF + CtoK * CI);    % initial enegy of water

M0_total = sum(W) + sum(M);       % total mass [kg]
E0_total = sum(EI) + sum(EW);     % total energy [J]

T_bottom = T(end);
m        = length(T);

% store initial mass [kg] and energy [J]
M  = dz .* d;                  % grid cell mass [kg]
EI = M .* T * CI;              % initial enegy of snow/ice
EW = W .* (LF + CtoK * CI);    % initial enegy of water

Z_cumulative = cumsum(dz);

% A logical "mask" that indicates which cells are in the top layers:
top_layers = Z_cumulative <= (column_ztop + d_tolerance );

% Define column_dzmin2 array using the top-layers' column_dzmin value for the entire column:
column_dzmin2 = column_dzmin * ones(m,1);

% Overwrite the bottom layers as the cumulative product times the stretching factor:
column_dzmin2(~top_layers) = cumprod(column_zy * ones(sum(~top_layers),1)) * column_dzmin;

% Define dzMax2 array using the top-layers' column_dzmin value for the entire column:
dzMax2 = 2 * column_dzmin * ones(m,1);

% In the bottom layers, dzMax2 is the larger of (column_zy * column_dzmin2) or (2 * column_dzmin)
dzMax2(~top_layers) = max(column_zy * column_dzmin2(~top_layers), 2 * column_dzmin);

%%

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
        T(i_target)         = (T(i) * M(i) + T(i_target) * M(i_target)) / m_new;
        a(i_target)         = (a(i) * M(i) + a(i_target) * M(i_target)) / m_new;
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
M(delete_cell)          = [];
W(delete_cell)          = [];
dz(delete_cell)         = [];
d(delete_cell)          = [];
T(delete_cell)          = [];
a(delete_cell)          = [];
re(delete_cell)         = [];
gdn(delete_cell)        = [];
gsp(delete_cell)        = [];
a_diffuse(delete_cell)  = [];
EI(delete_cell)         = [];
EW(delete_cell)         = [];
dzMax2(delete_cell)     = [];

% Calculate *new* length of cells:
m = length(T);

%% Split cells
% * An early implementation of this code used a loop which is included in comments at the bottom of this function for posterity.

% Find the cells that exceed tolerances:
f = find(dz > (dzMax2 + d_tolerance));

% Conserve quantities among the cells that will be split:
dz(f) = dz(f)/2;
W(f)  = W(f)/2;
M(f)  = M(f)/2;
EI(f) = EI(f)/2;
EW(f) = EW(f)/2;

% Sort the indices of all the cells including the ones that will be duplicated:
fs = sort([(1:m)';f]);

% Recreate the variables with split cells:
dz        = dz(fs);
W         = W(fs);
M         = M(fs);
T         = T(fs);
d         = d(fs);
a         = a(fs);
a_diffuse = a_diffuse(fs);
EI        = EI(fs);
EW        = EW(fs);
re        = re(fs);
gdn       = gdn(fs);
gsp       = gsp(fs);

%% CORRECT FOR TOTAL MODEL DEPTH
% WORKS FINE BUT HAS BEEN DISABLED FOR CONVIENCE OF MODEL OUTPUT
% INTERPRETATION

% Calculate total model depth:
Z_total = sum(dz);

if Z_total < (column_zmin - d_tolerance)

    % Mass and energy to be added:
    M_added   = M(end) + W(end);
    E_added   = T(end) * M(end) * CI + W(end) * (LF + CtoK * CI);

    % Add a grid cell of the same size and temperature to the bottom:
    dz        = [   dz;    dz(end)];
    T         = [    T;     T(end)];
    W         = [    W;     W(end)];
    M         = [    M;     M(end)];
    d         = [    d;     d(end)];
    a         = [    a;     a(end)];
    a_diffuse = [a_diffuse; a_diffuse(end)];
    re        = [   re;    re(end)];
    gdn       = [  gdn;   gdn(end)];
    gsp       = [  gsp;   gsp(end)];

elseif Z_total > column_zmax+d_tolerance 

    % Mass and energy loss:
    M_added   = -(M(end) + W(end));
    E_added   = -(T(end) * M(end) * CI) - W(end) * (LF+CtoK*CI);

    % Remove a grid cell from the bottom:
    dz(end)        = [];
    T(end)         = [];
    W(end)         = [];
    M(end)         = [];
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
E_added   = E_added + ((T_bottom - T(end)) * M(end) * CI);
T(end)    = T_bottom;

%% CHECK FOR MASS AND ENERGY CONSERVATION
if verbose
    % Calculate final mass [kg] and energy [J]

    EI    = M .* T * CI;
    EW    = W .* (LF + CtoK * CI);

    M1_total = sum(W) + sum(M);
    E1_total = sum(EI) + sum(EW);

    M_delta = round((M0_total - M1_total + M_added)*100)/100.;
    E_delta = round(E0_total - E1_total + E_added);

    if M_delta ~= 0 || E_delta ~= 0
        error(['Mass and energy are not conserved in melt equations:' newline ' M_delta: ' ...
            num2str(M_delta) ' E_delta: ' num2str(E_delta) newline])
    end

end

end