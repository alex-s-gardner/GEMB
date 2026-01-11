function [T, dz, d, W, re, gdn, gsp, a, a_diffuse, EC, M_surf, sw_net, shf, ...
    lhf, ulw, Ra, M, R, F, M_added, E_added, compaction_dens, compaction_melt] = ...
   gemb_core(T, dz, d, W, re, gdn, gsp, a, a_diffuse, EC, M_surf, ...
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
    % Determine initial mass [kg]:
    M_total_initial = sum (dz .* d) + sum(W);
    T_bottom        = T(end);
end

% 1. Snow grain metamorphism
% [always calculate, used in thermo and albedo]
[re, gdn, gsp] = ...
    grain_growth(T, dz, d, W, re, gdn, gsp, ClimateForcingStep, ModelParam);

% 2. Calculate snow, firn, and ice albedo
% Uses EC and M_surf from the previous time step
[a, a_diffuse] = ...
    albedo(T, dz, d, W, re, a, a_diffuse, EC, M_surf, ClimateForcingStep, ModelParam);

% 3. Determine distribution of absorbed SW radiation with depth
swf = ...
    shortwave(dz, d, re, a(1), a_diffuse(1), ClimateForcingStep, ModelParam);

% 4. Calculate net shortwave [W m-2]
sw_net = ...
    sum(swf);

% 5. Calculate new temperature-depth profile and turbulent heat fluxes [W m-2]
[T, shf, lhf, EC, ulw] = ...
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
% Calculate water production M [kg m-2], runoff R [kg m-2], and resulting changes
compaction_melt = ...
    sum(dz); % Track thickness before melt

[T, dz, d, W, re, gdn, gsp, a, a_diffuse, M, M_surf, R, F] = ...
    melt(T, dz, d, W, re, gdn, gsp, a, a_diffuse, Ra, ModelParam.density_ice, verbose);

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
    M_total_final = sum (dz .* d) + sum(W);
    M_change   = M_total_final + R - ClimateForcingStep.P - EC - M_total_initial - M_added;

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