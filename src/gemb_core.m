function [T, dz, d, W, re, gdn, gsp, a, a_diffuse, EC, melt_surface, sw_net, shf, ...
    lhf, ulw, Ra, melt, R, F, M_added, E_added, compaction_dens, compaction_melt] = ...
   gemb_core(T, dz, d, W, re, gdn, gsp, a, a_diffuse, EC, melt_surface, ...
    ClimateForcingStep, ModelParam, verbose)
% GEMB_STEP Performs a single time-step of the GEMB model.
%   Calculates grain growth, albedo, radiative transfer, thermodynamics, 
%   accumulation, melt, layer management, and densification.
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

if verbose
    % Specify constants:
    CtoK = 273.15;   % Celsius to Kelvin conversion
    CI   = 2102;     % specific heat capacity of snow/ice (J kg-1 K-1)
    LF   = 0.3345E6; % latent heat of fusion (J kg-1)

    M               = dz .* d;
    M_total_initial = sum(M) + sum(W);        % total mass [kg]
    E_total_initial = sum(M .* T * CI) + ...
        sum(W .* (LF + CtoK * CI));           % total energy [J] = initial enegy of snow/ice + initial enegy of water
    T_bottom        = T(end);

    % Determine initial energy [kg]:

end

% 1. Snow grain metamorphism
% [always calculate, used in thermo and albedo]
[re, gdn, gsp] = ...
    grain_growth(T, dz, d, W, re, gdn, gsp, ClimateForcingStep, ModelParam);

% 2. Calculate snow, firn, and ice albedo
% Uses EC and melt_surface from the previous time step
[a, a_diffuse] = ...
    albedo(T, dz, d, W, re, a, a_diffuse, EC, melt_surface, ClimateForcingStep, ModelParam);

% 3. Determine distribution of absorbed SW radiation with depth
swf = ...
    shortwave(dz, d, re, a(1), a_diffuse(1), ClimateForcingStep, ModelParam);

% 4. Calculate net shortwave [W m-2]
sw_net = ...
    sum(swf);

% 5. Calculate new temperature-depth profile and turbulent heat fluxes [W m-2]
[T, ulw, shf, lhf, EC] = ...
    thermo(T, dz, d, W(1), re, swf, ClimateForcingStep, ModelParam, verbose);

% 6. Change in thickness of top cell due to evaporation/condensation
% Assuming same density as top cell
% ## NEED TO FIX THIS IN CASE ALL OR MORE OF CELL EVAPORATES ##
dz(1) = ...
    dz(1) + EC / d(1);

% 7. Add snow/rain to top grid cell adjusting cell depth, temperature, and density
[T, dz, d, W, re, gdn, gsp, a, a_diffuse, Ra] = ...
    accumulation(T, dz, d, W, re, gdn, gsp, a, a_diffuse, ClimateForcingStep, ModelParam);

% 8. Melt and wet compaction
% Calculate water production melt [kg m-2], runoff R [kg m-2], and resulting changes
compaction_melt = ...
    sum(dz); % Track thickness before melt

[T, dz, d, W, re, gdn, gsp, a, a_diffuse, melt, melt_surface, R, F] = ...
    melting(T, dz, d, W, re, gdn, gsp, a, a_diffuse, Ra, ModelParam.density_ice, verbose);

compaction_melt = ...
    (compaction_melt - sum(dz)); % Calculate wet compaction

% 9. Manage the layering to match user defined requirements
[T, dz, d, W, re, gdn, gsp, a, a_diffuse, M_added, E_added] = ...
    layer_management(T, dz, d, W, re, gdn, gsp, a, a_diffuse, ModelParam, verbose);

% 10. Allow non-melt densification and determine compaction [m]
compaction_dens = ...
    sum(dz); % Track thickness before densification

[dz, d] = ...
    densification(T, dz, d, re, ClimateForcingStep, ModelParam);

compaction_dens = ...
    (compaction_dens - sum(dz)); % Calculate dry compaction


if verbose
    % calculate total system mass
    M = dz .* d;
    M_total_final = sum(M) + sum(W);        % total mass [kg]
    M_change   = M_total_final - M_total_initial + R - ClimateForcingStep.P - EC  - M_added;


    ER_total = R_total * (LF + CtoK * CI);
    EI       = M0 .* T * CI;
    EW       = W .* (LF + CtoK * CI);
    dt = ClimateForcingStep.dt;
    E_total_final = sum(M .* T * CI) + ...
        sum(W .* (LF + CtoK * CI));         % total energy [J] = initial enegy of snow/ice + initial enegy of water

    E_change   = E_total_final - E_total_initial - (sw_net*dt) - (ulw*dt) - (shf*dt) - (lhf*dt);

    % check mass conservation
    if abs(M_change) > 1E-3
        error('total system mass not conserved in MB function')
    end

    % check bottom grid cell T is unchanged
    if abs(T(end)-T_bottom) > 1E-3
        error('temperature of bottom grid cell changed: original = %0.10g J, updated = %0.10g J',T_bottom,T(end))
    end
end

end