function [T, dz, d, water, re, gdn, gsp, a, a_diffuse, melt_total, melt_surface, runoff_total, freeze_total] = ...
    calculate_melt(T, dz, d, water, re, gdn, gsp, a, a_diffuse, rain, ModelParam, verbose)
% calculate_melt computes the quantity of meltwater due to snow temperature in excess
% of 0 deg C, determines pore water content and adjusts grid spacing.
%
%% Syntax
%
% [T, dz, d, water, re, gdn, gsp, a, a_diffuse, melt_total, melt_surface, runoff_total, freeze_total] = ...
%    calculate_melt(T, dz, d, water, re, gdn, gsp, a, a_diffuse, rain, ModelParam, verbose)
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
%  T                                : K            Grid cell temperature.
%  dz                               : m            Grid cell thickness.
%  d                                : kg m^-3      Grid cell density.
%  water                            : kg m^-2      Water content.
%  re                               : mm           Grain size (effective radius).
%  gdn                              : unitless     Grain dendricity.
%  gsp                              : unitless     Grain sphericity.
%  a                                : fraction     Albedo.
%  a_diffuse                        : fraction     Diffuse albedo.
%  rain                             : kg m^-2      Rainfall amount added to the column.
%  ModelParam                       : struct       Model parameters:
%     .density_ice                  : kg m^-3   Density of ice.
%     .water_irreducible_saturation :fractionirreducible water content saturation [fraction]
%     .density_ice          
%  verbose                          : logical      Flag to enable mass/energy conservation checks.
%
%% Outputs
%
%  T            : K            Updated grid cell temperature.
%  dz           : m            Updated grid cell thickness.
%  d            : kg m^-3      Updated grid cell density.
%  water        : kg m^-2      Updated water content.
%  re           : mm           Updated grain size.
%  gdn          : unitless     Updated grain dendricity.
%  gsp          : unitless     Updated grain sphericity.
%  a            : fraction     Updated albedo.
%  a_diffuse    : fraction     Updated diffuse albedo.
%  melt_total   : kg m^-2      Total column melt mass.
%  melt_surface : kg m^-2      Melt mass generated in the surface layer only.
%  runoff_total : kg m^-2      Total runoff leaving the column.
%  freeze_total : kg m^-2      Total mass refrozen within the column.
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

T_tolerance      = 1e-10;
d_tolerance      = 1e-11;
water_tolerance  = 1e-13;

% Specify constants:
CtoK            = 273.15;   % Celsius to Kelvin conversion
CI              = 2102;     % specific heat capacity of snow/ice (J kg-1 K-1)
LF              = 0.3345E6; % latent heat of fusion (J kg-1)
d_phc           = 830.0;    % pore hole close off density [kg m-3]
ice_layer_dzmin = 0.1;      % if the density is greater than d_phc and the layer thickness exceeds ice_layer_dzmin [m], then all meltwater percolating down is counted as runoff

m            = length(T);
water_delta  = zeros(m,1);

% store initial mass [kg] and energy [J]
M = dz .* d;                   % grid cell mass [kg]

if verbose
    M_total_initial = sum(water)  + sum(M);      % total mass [kg]

    E_total_initial = sum(M .* T * CI) + ...
        sum(water .* (LF + CtoK * CI));          % total energy [J] = initial enegy of snow/ice + initial enegy of water
end

% initialize melt and runoff scalars
runoff_total = 0;   % sum runoff [kg m^-2]
melt_total   = 0;   % total melt [kg m^-2]
melt_surface = 0;   % surface layer melt

% calculate temperature excess above 0 degC
T_excess = max(0, T - CtoK);        % [K] to [degC]

% new grid point center temperature, T [K]
T = min(T, CtoK);

% specify irreducible water content saturation [fraction]
water_irreducible_saturation = 0.07;  % assumed constant after Colbeck, 1974

%% REFREEZE PORE WATE_runoff

% check if any pore water
if sum(water) > water_tolerance 
    % disp('PORE WATE_runoff REFREEZE')
    % calculate maximum freeze amount, freeze_max [kg]
    freeze_max = max(0, -((T - CtoK) .* M * CI) / LF);

    % freeze pore water and change snow/ice properties
    water_delta  = min(freeze_max, water);                                                     % freeze mass [kg]
    water        = water - water_delta;                                                        % pore water mass [kg]
    M            = M + water_delta;                                                            % new mass [kg]
    d            = M ./ dz;                                                                    % density [kg m-3]
    T            = T + double(M>water_tolerance ).*(water_delta.*(LF+(CtoK - T)*CI)./(M.*CI)); % temperature [K]

    % if pore water froze in ice then adjust d and dz thickness
    d(d > ModelParam.density_ice - d_tolerance ) = ModelParam.density_ice;
    dz = M ./ d;

end

% squeeze water from snow pack
water_irreducible = (ModelParam.density_ice - d) .* ModelParam.water_irreducible_saturation .* (M ./ d);          % irreducible water content [kg m^-2]
water_excess      = max(0, water - water_irreducible);             % water "squeezed" from snow [kg m^-2]

%% MELT, PE_runoffCOLATION AND REFREEZE
freeze = zeros(m,1);

% Add previous freeze to freeze and reset water_delta
freeze          = freeze + water_delta;
water_delta(:)  = 0;

% run melt algorithm if there is melt water or excess pore water
if (sum(T_excess) > T_tolerance) || (sum(water_excess) > water_tolerance)
    % disp ('MELT OCCURS')
    
    % Check to see if thermal energy exceeds energy to melt entire cell
    % if so redistribute temperature to lower cells (temperature surplus)
    % (maximum T of snow before entire grid cell melts is a constant
    % LF/CI = 159.1342)
    T_surplus = max(0, T_excess - LF/CI);

    if sum(T_surplus) > T_tolerance

        % calculate surplus energy
        E_surplus = T_surplus .* CI .* M;
        i         = 1;

        while (sum(E_surplus) > T_tolerance)  && (i < (m+1))
            % entire cell melts

            if i<m
                % use surplus energy to increase the temperature of lower cell
                T(i+1)        = E_surplus(i) / M(i+1) / CI + T(i+1);

                T_excess(i+1) = max(0, T(i+1) - CtoK) + T_excess(i+1);
                T(i+1)        = min(CtoK, T(i+1));

                T_surplus(i+1) = max(0, T_excess(i+1) - LF/CI);
                E_surplus(i+1) = T_surplus(i+1) * CI * M(i+1);
            else
                error('surplus energy reached the base of gemb column (i.e. entire column melted out in a single time step)')
            end

            % adjust current cell properties (again 159.1342 is the max T)
            T_excess(i)  = LF/CI;
            E_surplus(i) = 0;
            i            = i + 1;
        end
    end

    % convert temperature excess to melt [kg]
    melt_maximum  = T_excess .* d .* dz * CI / LF;
    melt          = min(melt_maximum, M);               % melt
    melt_surface  = melt(1);
    melt_total    = max(0,sum(melt)-rain);              % total melt [kg] minus the liquid rain that had been added

    % calculate maximum refreeze amount, freeze_max [kg]
    freeze_max = max(0, -((T - CtoK) .* d .* dz * CI)/ LF);

    % initialize refreeze, runoff, flux_dn and water_delta vectors [kg]
    runoff  = zeros(m,1);
    flux_dn = [runoff; 0];

    % determine the deepest grid cell where melt/pore water is generated
    X             = find((melt > water_tolerance)  | (water_excess > water_tolerance), 1, 'last');
    X(isempty(X)) = 1;

    Xi = 1;
    m  = length(T);

    % meltwater percolation
    for i = 1:m
        % calculate total melt water entering cell
        melt_input = melt(i) + flux_dn(i);

        ice_depth = 0;
        % If this grid cell's density exceeds the pore closeoff density:
        if d(i) >= d_phc-d_tolerance 
            for l=i:m
                if d(l)>=d_phc-d_tolerance 
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
        if abs(melt_input) < water_tolerance  && i > X
            break

            % if reaches impermeable ice layer all liquid water runs off (runoff)
        elseif (d(i) >= (ModelParam.density_ice-d_tolerance))  || ((d(i) >= d_phc-d_tolerance)  && (ice_depth > ice_layer_dzmin+d_tolerance))  % d_phc = pore hole close off [kg m-3]
            % disp('ICE LAYE_runoff')
            % no water freezes in this cell
            % no water percolates to lower cell
            % cell ice temperature & density do not change

            M(i)              = M(i) - melt(i);                                                 % mass after melt
            water_irreducible = (ModelParam.density_ice-d(i)) * water_irreducible_saturation * (M(i)/d(i));                         % irreducible water
            water_delta(i)    = max(min(melt_input, water_irreducible - water(i)), -water(i)); % change in pore water
            runoff(i)         = max(0, melt_input - water_delta(i));                            % runoff

            % check if no energy to refreeze meltwater
        elseif abs(freeze_max(i)) < d_tolerance 
            % disp('REFREEZE == 0')
            % no water freezes in this cell
            % cell ice temperature & density do not change

            M(i)              = M(i) - melt(i);                                    % mass after melt
            water_irreducible = (ModelParam.density_ice-d(i)) * water_irreducible_saturation * (M(i)/d(i));         % irreducible water
            water_delta(i)    = max(min(melt_input, water_irreducible - water(i)),-1*water(i)); % change in pore water
            flux_dn(i+1)      = max(0, melt_input - water_delta(i));                    % meltwater out
            runoff(i)         = 0;

            % some or all meltwater refreezes
        else
            % change in density and temperature
            % disp('MELT REFREEZE')
            %-----------------------melt water-----------------------------
            M(i)    = M(i) - melt(i);
            dz_0    = M(i)/d(i);
            d_max   = (ModelParam.density_ice - d(i))*dz_0;                 % d max = ModelParam.density_ice
            freeze1 = min(min(melt_input,d_max),freeze_max(i));  % maximum refreeze
            M(i)    = M(i) + freeze1;                            % mass after refreeze
            d(i)    = M(i)/dz_0;

            %-----------------------pore water-----------------------------
            water_irreducible = (ModelParam.density_ice-d(i))* water_irreducible_saturation * dz_0;                     % irreducible water
            water_delta(i)    = max(min(melt_input - freeze1, water_irreducible-water(i)),-1 * water(i)); % change in pore water
            freeze2           = 0;

            % ---------------- THIS HAS NOT BEEN CHECKED-----------------_
            if water_delta(i) < 0.0-water_tolerance                  % excess pore water
                d_max       = (ModelParam.density_ice - d(i))*dz_0;  % maximum refreeze
                freeze2_max = min(d_max, freeze_max(i)-freeze1);     % maximum refreeze
                freeze2     = min(-1.0*water_delta(i), freeze2_max); % pore water refreeze
                M(i)        = M(i) + freeze2;                        % mass after refreeze
                d(i)        = M(i)/dz_0;
            end
            % -------------------------------------------------------------

            freeze(i)    = freeze(i) + freeze1 + freeze2;
            flux_dn(i+1) = max(0.0, melt_input - freeze1 - water_delta(i)); % meltwater out

            if M(i) > water_tolerance 
                T(i) = T(i) + ...                               % change in temperature
                    ((freeze1+freeze2)*(LF+(CtoK - T(i))*CI)./(M(i).*CI));
            end

            % check if an ice layer forms
            if abs(d(i) - ModelParam.density_ice) < d_tolerance 
                % disp('ICE LAYE_runoff FORMS')
                % excess water runs off
                runoff(i)    = flux_dn(i+1);

                % no water percolates to lower cell
                flux_dn(i+1) = 0;
            end
        end

        Xi = Xi+1;
    end

    % GRID CELL SPACING AND MODEL DEPTH
    if verbose
        if any(water < 0.0-water_tolerance )
            error('Negative pore water generated in melt equations.')
        end
    end

    % delete all cells with zero mass
    % adjust pore water
    water = water + water_delta;

    % calculate runoff_total:
    runoff_total = sum(runoff) + flux_dn(Xi);

    % delete all cells with zero mass
    D            = (M <= 0+water_tolerance );
    M(D)         = [];
    water(D)     = [];
    d(D)         = [];
    T(D)         = [];
    a(D)         = [];
    re(D)        = [];
    gdn(D)       = [];
    gsp(D)       = [];
    a_diffuse(D) = [];

    % calculate new grid lengths
    dz = M ./ d;
end

freeze_total = sum(freeze);

%% CHECK FOR MASS AND ENE_runoffGY CONSE_runoffVATION

if verbose
    % Calculate final mass [kg] and energy [J]
    E_total_runoff = runoff_total * (LF + CtoK * CI);

    M_total_final = sum(water) + sum(M) + runoff_total;
    E_total_final = sum(M .* T * CI) + sum(water .* (LF + CtoK * CI));

    M_delta = M_total_initial - M_total_final;
    E_delta = E_total_initial - E_total_final - E_total_runoff;

    if (abs(M_delta) > 1E-3) || (abs(E_delta) > 1E-3)
        error(['Mass and/or energy are not conserved in melt equations:' newline ' M_delta: ' ...
            num2str(M_delta) ' E_delta: ' num2str(E_delta) newline])
    end

    if any(water < 0.0-water_tolerance )
        error('Negative pore water generated in melt equations.')
    end
end

