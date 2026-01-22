function [T, dz, d, water, re, gdn, gsp, a, a_diffuse, EC, melt_surface, sw_net, shf, ...
    lhf, ulw, Ra, melt, R, F, M_added, E_added, compaction_dens, compaction_melt] = ...
   gemb_core(T, dz, d, water, re, gdn, gsp, a, a_diffuse, EC, melt_surface, ...
    ClimateForcingStep, ModelParam, verbose)
% gemb_core performs a single time-step of the GEMB model.
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
    M_total_initial = sum(M) + sum(water);        % total mass [kg]
    E_total_initial = sum(M .* T * CI) + ...
        sum(water .* (LF + CtoK * CI));           % total energy [J] = initial enegy of snow/ice + initial enegy of water
    T_bottom        = T(end);

    % Determine initial energy [kg]:

end

% 1. Snow grain metamorphism
% [always calculate, used in thermo and albedo]
[re, gdn, gsp] = ...
    calculate_grain_size(T, dz, d, water, re, gdn, gsp, ClimateForcingStep, ModelParam);

% 2. Calculate snow, firn, and ice albedo
% Uses EC and melt_surface from the previous time step
[a, a_diffuse] = ...
    calculate_albedo(T, dz, d, water, re, a, a_diffuse, EC, melt_surface, ClimateForcingStep, ModelParam);

% 3. Determine distribution of absorbed SW radiation with depth
swf = calculate_shortwave_radiation(dz, d, re, a(1), a_diffuse(1), ClimateForcingStep, ModelParam);

% 4. Calculate net shortwave [W m-2]
sw_net = sum(swf);

% 5. Calculate new temperature-depth profile and turbulent heat fluxes [W m-2]
[T, ulw, shf, lhf, ghf, EC] = ...
    calculate_temperature(T, dz, d, water(1), re, swf, ClimateForcingStep, ModelParam, verbose);

% 4. Calculate net longwave [W m-2]
lw_net = ClimateForcingStep.dlw - ulw;

% 6. Change in thickness of top cell due to evaporation/condensation
% Assuming same density as top cell
% ## NEED TO FIX THIS IN CASE ALL OR MORE OF CELL EVAPORATES ##
dz(1) = dz(1) + EC / d(1);

if verbose
    E_EC = EC * T(1) * CI;
end

% 7. Add snow/rain to top grid cell adjusting cell depth, temperature, and density
[T, dz, d, water, re, gdn, gsp, a, a_diffuse, Ra] = ...
    calculate_accumulation(T, dz, d, water, re, gdn, gsp, a, a_diffuse, ClimateForcingStep, ModelParam, verbose);

% 8. Melt and wet compaction
% Calculate water production melt [kg m-2], runoff R [kg m-2], and resulting changes
compaction_melt = sum(dz); % Track thickness before melt

[T, dz, d, water, re, gdn, gsp, a, a_diffuse, melt, melt_surface, R, F] = ...
    calculate_melt(T, dz, d, water, re, gdn, gsp, a, a_diffuse, Ra, ModelParam, verbose);

compaction_melt = (compaction_melt - sum(dz)); % Calculate wet compaction

% 9. Manage the layering to match user defined requirements
[T, dz, d, water, re, gdn, gsp, a, a_diffuse, M_added, E_added] = ...
    manage_layers(T, dz, d, water, re, gdn, gsp, a, a_diffuse, ModelParam, verbose);

% 10. Allow non-melt densification and determine compaction [m]
compaction_dens = sum(dz); % Track thickness before densification

[dz, d] = calculate_density(T, dz, d, re, ClimateForcingStep, ModelParam);

compaction_dens = ...
    (compaction_dens - sum(dz)); % Calculate dry compaction


if verbose
    dt = ClimateForcingStep.dt;

    % calculate total system mass
    M = dz .* d;
    M_total_final = sum(M) + sum(water);        % total mass [kg]
    M_delta       = M_total_final - M_total_initial + R - ClimateForcingStep.P - EC - M_added;

    % check mass conservation
    if abs(M_delta) > 1E-3
        error('total system mass not conserved: M_delta = %0.4f', M_delta)
    end

    % need to account for rain
    E_snow     = ((ClimateForcingStep.P - Ra) * (ClimateForcingStep.T_air) * CI);
    E_rain     = Ra * (ClimateForcingStep.T_air * CI + LF);
    E_R        = R * (LF + CtoK * CI);
    E_thermal  = sum((dz .* d) .* T * CI);
    E_water    = sum(water .* (LF + CtoK * CI));
    E_sw       = (sw_net*dt);
    E_lw       = (lw_net*dt);
    E_thf      = ((shf+lhf)*dt);
    E_ghf      = (ghf*dt);

    E_total_final = E_thermal + E_water + E_R; % total energy [J] = initial enegy of snow/ice + initial enegy of water
  
    E_used     = E_total_final - E_total_initial;
    E_supplied = E_sw + E_lw + E_thf + E_snow + E_rain + E_ghf + E_EC + E_added ;
    E_delta   = E_used - E_supplied;

    % check energy conservation
    if abs(E_delta) > 1E-3
        error('total system energy not conserved: E_delta = %0.4f \n E_sw = %0.4f, E_lw = %0.4f, E_thf = %0.4f, E_added = %0.4f', E_delta, E_sw, E_lw, E_thf, E_added)
    end

    % check bottom grid cell T is unchanged
    if abs(T(end)-T_bottom) > 1E-3
        error('temperature of bottom grid cell changed: original = %0.10g J, updated = %0.10g J',T_bottom,T(end))
    end
end

end