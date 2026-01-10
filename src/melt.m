function [T, dz, d, W, re, gdn, gsp, a, a_diffuse, M_total, M_surf, R_total, F_total] = ...
    melt(T, dz, d, W, re, gdn, gsp, a, a_diffuse, Ra, density_ice, verbose)

% melt computes the quantity of meltwater due to snow temperature in excess
% of 0 deg C, determines pore water content and adjusts grid spacing.
%
%% Syntax
%
% [T, dz, d, W, re, gdn, gsp, a, a_diffuse, M_total, M_surf, R_total, F_total] = ...
%    melt(T, dz, d, W, re, gdn, gsp, a, a_diffuse, Ra, density_ice, verbose)
%
%% Description
%
% This function calculates the thermodynamic and hydraulic evolution of the 
% firn column during surface melt events . 
% It employs a "tipping bucket" approach to simulate the percolation, refreezing, 
% and retention of liquid water. The specific processes include:
%
% 1. Initial Refreeze: Existing pore water in layers with $T < 0^\circ C$ is 
%    refrozen, releasing latent heat and warming the layer.
% 2. Melt Generation: If layer temperatures exceed $0^\circ C$, the excess 
%    energy is converted into liquid meltwater. If the energy excess is greater 
%    than that required to melt the entire layer, the surplus energy is 
%    transferred to the layer below.
% 3. Percolation: Liquid water (melt + rain) percolates downward. It flows 
%    through the snowpack until it either:
%    * Refreezes: In cold underlying layers, increasing density and temperature.
%    * Retains: As pore water, up to the irreducible water content saturation 
%        ($S_{wi} = 0.07$).
%    * Runs off: If it encounters an impermeable ice lens (density > 
%        pore close-off) or reaches the bottom of the column.
%
% The function maintains strict conservation of mass and energy, checking budgets 
% if the verbose flag is enabled.
%
%% Inputs
%
%  T            : K            Grid cell temperature.
%  dz           : m            Grid cell thickness.
%  d            : kg m^-3      Grid cell density.
%  W            : kg m^-2      Water content.
%  re           : mm           Grain size (effective radius).
%  gdn          : unitless     Grain dendricity.
%  gsp          : unitless     Grain sphericity.
%  a            : fraction     Albedo.
%  a_diffuse    : fraction     Diffuse albedo.
%  Ra           : kg m^-2      Rainfall amount added to the column.
%  density_ice  : kg m^-3      Density of ice (constant).
%  verbose      : logical      Flag to enable mass/energy conservation checks.
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
%  M_total      : kg m^-2      Total column melt mass.
%  M_surf       : kg m^-2      Melt mass generated in the surface layer only.
%  R_total      : kg m^-2      Total runoff leaving the column.
%  F_total      : kg m^-2      Total mass refrozen within the column.
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

%% INITIALIZATION

T_tolerance  = 1e-10;
d_tolerance  = 1e-11;
W_tolerance  = 1e-13;

% Specify constants:
CtoK        = 273.15;   % Celsius to Kelvin conversion
CI          = 2102;     % specific heat capacity of snow/ice (J kg-1 K-1)
LF          = 0.3345E6; % latent heat of fusion (J kg-1)
density_phc = 830.0;    % pore hole close off density [kg m-3]
ice_layer_dzmin = 0.1; % if the density is greater than density_phc and the layer thickness exceeds ice_layer_dzmin [m], then all meltwater percolating down is counted as runoff

m         = length(T);
W_delta   = zeros(m,1);

% store initial mass [kg] and energy [J]
M0  = dz .* d;                   % grid cell mass [kg]
EI  = M0 .* T * CI;              % initial enegy of snow/ice
EW  = W .* (LF + CtoK * CI);     % initial enegy of water

M0_total = sum(W) + sum(M0);       % total mass [kg]
E0_total = sum(EI) + sum(EW);      % total energy [J]

% initialize melt and runoff scalars
R_total   = 0;   % sum runoff [kg m^-2]
M_total   = 0;   % total melt [kg m^-2]
M_surf    = 0;   % surface layer melt

% output
E_surplus = 0;

% calculate temperature excess above 0 degC
T_excess = max(0, T - CtoK);        % [K] to [degC]

% new grid point center temperature, T [K]
T = min(T,CtoK);

% specify irreducible water content saturation [fraction]
Swi = 0.07;  % assumed constant after Colbeck, 1974

%% REFREEZE PORE WATER
% check if any pore water
if sum(W) > W_tolerance 
    % disp('PORE WATER REFREEZE')
    % calculate maximum freeze amount, F_max [kg]
    F_max = max(0, -((T - CtoK) .* M0 * CI) / LF);

    % freeze pore water and change snow/ice properties
    W_delta  = min(F_max, W);                                                        % freeze mass [kg]
    W        = W - W_delta;                                                          % pore water mass [kg]
    M0       = M0 + W_delta;                                                         % new mass [kg]
    d        = M0 ./ dz;                                                             % density [kg m-3]
    T        = T + double(M0>W_tolerance ).*(W_delta.*(LF+(CtoK - T)*CI)./(M0.*CI)); % temperature [K]

    % if pore water froze in ice then adjust d and dz thickness
    d(d > density_ice-d_tolerance ) = density_ice;
    dz = M0 ./ d;

end

% squeeze water from snow pack
W_irreducible = (density_ice - d) .* Swi .* (M0 ./ d);          % irreducible water content [kg m^-2]
W_excess      = max(0, W - W_irreducible);                      % water "squeezed" from snow [kg m^-2]

%% MELT, PERCOLATION AND REFREEZE
F = zeros(m,1);

% Add previous refreeze to F and reset W_delta
F          = F + W_delta;
W_delta(:) = 0;

% run melt algorithm if there is melt water or excess pore water
if (sum(T_excess) > T_tolerance) || (sum(W_excess) > W_tolerance)
    % disp ('MELT OCCURS')
    % check to see if thermal energy exceeds energy to melt entire cell
    % if so redistribute temperature to lower cells (temperature surplus)
    % (maximum T of snow before entire grid cell melts is a constant
    % LF/CI = 159.1342)
    T_surplus = max(0, T_excess - LF/CI);

    if sum(T_surplus) > T_tolerance

        % calculate surplus energy
        E_surplus = T_surplus .* CI .* M0;
        i         = 1;

        while (sum(E_surplus) > T_tolerance)  && (i < (m+1))

            if i<m
                % use surplus energy to increase the temperature of lower cell
                T(i+1)        = E_surplus(i) / M0(i+1) / CI + T(i+1);

                T_excess(i+1) = max(0, T(i+1) - CtoK) + T_excess(i+1);
                T(i+1)        = min(CtoK, T(i+1));

                T_surplus(i+1) = max(0, T_excess(i+1) - LF/CI);
                E_surplus(i+1) = T_surplus(i+1) * CI * M0(i+1);
            else
                E_surplus      = E_surplus(i);
                warning('surplus energy at the base of GEMB column' + newline)
            end

            % adjust current cell properties (again 159.1342 is the max T)
            T_excess(i)  = LF/CI;
            E_surplus(i) = 0;
            i            = i + 1;

        end
    end

    % convert temperature excess to melt [kg]
    M_max    = T_excess .* d .* dz * CI / LF;
    M        = min(M_max, M0);                 % melt
    M_surf   = M(1);
    M_total  = max(0,sum(M)-Ra);              % total melt [kg] minus the liquid rain that had been added

    % calculate maximum refreeze amount, F_max [kg]
    F_max = max(0, -((T - CtoK) .* d .* dz * CI)/ LF);

    % initialize refreeze, runoff, flux_dn and W_delta vectors [kg]
    R       = zeros(m,1);
    flux_dn = [R; 0];

    % determine the deepest grid cell where melt/pore water is generated
    X             = find((M > W_tolerance)  | (W_excess > W_tolerance), 1, 'last');
    X(isempty(X)) = 1;

    Xi = 1;
    m  = length(T);

    % meltwater percolation
    for i = 1:m
        % calculate total melt water entering cell
        M_input = M(i) + flux_dn(i);

        ice_depth = 0;
        % If this grid cell's density exceeds the pore closeoff density:
        if d(i) >= density_phc-d_tolerance 
            for l=i:m
                if d(l)>=density_phc-d_tolerance 
                    ice_depth = ice_depth+dz(l);
                    if ice_depth > ice_layer_dzmin + d_tolerance % OPTIMIZATION: Early break
                        break;
                    end
                else
                    break
                end
            end
        end

        % break loop if there is no meltwater and if depth is > mw_depth
        if abs(M_input) < W_tolerance  && i > X
            break

            % if reaches impermeable ice layer all liquid water runs off (R)
        elseif (d(i) >= (density_ice-d_tolerance))  || ((d(i) >= density_phc-d_tolerance)  && (ice_depth > ice_layer_dzmin+d_tolerance))  % density_phc = pore hole close off [kg m-3]
            % disp('ICE LAYER')
            % no water freezes in this cell
            % no water percolates to lower cell
            % cell ice temperature & density do not change

            M0(i)         = M0(i) - M(i);                                    % mass after melt
            W_irreducible = (density_ice-d(i)) * Swi * (M0(i)/d(i));         % irreducible water
            W_delta(i)    = max(min(M_input, W_irreducible - W(i)),-1*W(i)); % change in pore water
            R(i)          = max(0, M_input - W_delta(i));                    % runoff

            % check if no energy to refreeze meltwater
        elseif abs(F_max(i)) < d_tolerance 
            % disp('REFREEZE == 0')
            % no water freezes in this cell
            % cell ice temperature & density do not change

            M0(i)         = M0(i) - M(i);                                    % mass after melt
            W_irreducible = (density_ice-d(i)) * Swi * (M0(i)/d(i));         % irreducible water
            W_delta(i)    = max(min(M_input, W_irreducible - W(i)),-1*W(i)); % change in pore water
            flux_dn(i+1)  = max(0, M_input - W_delta(i));                    % meltwater out
            R(i)          = 0;

            % some or all meltwater refreezes
        else
            % change in density and temperature
            % disp('MELT REFREEZE')
            %-----------------------melt water-----------------------------
            M0(i) = M0(i) - M(i);
            dz_0  = M0(i)/d(i);
            d_max = (density_ice - d(i))*dz_0;         % d max = density_ice
            F1    = min(min(M_input,d_max),F_max(i));  % maximum refreeze
            M0(i) = M0(i) + F1;                        % mass after refreeze
            d(i)  = M0(i)/dz_0;

            %-----------------------pore water-----------------------------
            W_irreducible = (density_ice-d(i))* Swi * dz_0;                     % irreducible water
            W_delta(i)    = max(min(M_input - F1, W_irreducible-W(i)),-1*W(i)); % change in pore water
            F2            = 0;

            % ---------------- THIS HAS NOT BEEN CHECKED-----------------_
            if W_delta(i) < 0.0-W_tolerance            % excess pore water
                d_max  = (density_ice - d(i))*dz_0;    % maximum refreeze
                F2_max = min(d_max, F_max(i)-F1);      % maximum refreeze
                F2     = min(-1.0*W_delta(i), F2_max); % pore water refreeze
                M0(i)  = M0(i) + F2;                   % mass after refreeze
                d(i)   = M0(i)/dz_0;
            end
            % -------------------------------------------------------------

            F(i)         = F(i) + F1 + F2;
            flux_dn(i+1) = max(0.0, M_input - F1 - W_delta(i)); % meltwater out

            if M0(i) > W_tolerance 
                T(i) = T(i) + ...                               % change in temperature
                    ((F1+F2)*(LF+(CtoK - T(i))*CI)./(M0(i).*CI));
            end

            % check if an ice layer forms
            if abs(d(i) - density_ice) < d_tolerance 
                % disp('ICE LAYER FORMS')
                % excess water runs off
                R(i) = flux_dn(i+1);

                % no water percolates to lower cell
                flux_dn(i+1) = 0;
            end
        end

        Xi = Xi+1;
    end

    % GRID CELL SPACING AND MODEL DEPTH
    if verbose
        if any(W < 0.0-W_tolerance )
            error('Negative pore water generated in melt equations.')
        end
    end

    % delete all cells with zero mass
    % adjust pore water
    W = W + W_delta;

    % calculate R_total:
    R_total = sum(R) + flux_dn(Xi);

    % delete all cells with zero mass
    D            = (M0 <= 0+W_tolerance );
    M0(D)        = [];
    W(D)         = [];
    d(D)         = [];
    T(D)         = [];
    a(D)         = [];
    re(D)        = [];
    gdn(D)       = [];
    gsp(D)       = [];
    a_diffuse(D) = [];

    % calculate new grid lengths
    dz = M0 ./ d;
end

F_total = sum(F);

%% CHECK FOR MASS AND ENERGY CONSERVATION
if verbose
    % Calculate final mass [kg] and energy [J]
    ER_total = R_total * (LF + CtoK * CI);
    EI       = M0 .* T * CI;
    EW       = W .* (LF + CtoK * CI);

    M1_total = sum(W) + sum(M0) + R_total;
    E1_total = sum(EI) + sum(EW);

    M_delta = round((M0_total - M1_total)*100)/100.;
    E_delta = round(E0_total - E1_total - ER_total - E_surplus);

    if M_delta ~= 0 || E_delta ~= 0
        error(['Mass and energy are not conserved in melt equations:' newline ' M_delta: ' ...
            num2str(M_delta) ' E_delta: ' num2str(E_delta) newline])
    end

    if any(W < 0.0-W_tolerance )
        error('Negative pore water generated in melt equations.')
    end
end

