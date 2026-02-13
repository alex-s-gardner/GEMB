function [temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, evaporation_condensation, melt_surface, shortwave_net, heat_flux_sensible, ...
    heat_flux_latent, longwave_upward, rain, melt, runoff, refreeze, mass_added, E_added, densification_from_compaction, densification_from_melt] = ...
   gemb_core(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, evaporation_condensation, melt_surface, ...
    ClimateForcingStep, ModelParam, verbose)
% gemb_core performs a single time-step of the GEMB model.
%   Calculates grain growth, albedo, radiative transfer, thermodynamics, 
%   accumulation, melt, layer management, and densification.
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
    C_ice   = 2102;     % specific heat capacity of snow/ice (J kg-1 K-1)
    LF   = 0.3345E6; % latent heat of fusion (J kg-1)

    M               = dz .* density;
    M_total_initial = sum(M) + sum(water);        % total mass [kg]
    E_total_initial = sum(M .* temperature * C_ice) + ...
        sum(water .* (LF + CtoK * C_ice));           % total energy [J] = initial enegy of snow/ice + initial enegy of water
    T_bottom        = temperature(end);

    % Determine initial energy [kg]:

end

% 1. Snow grain metamorphism
% [always calculate, used in thermo and albedo]
[grain_radius, grain_dendricity, grain_sphericity] = ...
    calculate_grain_size(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, ClimateForcingStep, ModelParam);

% 2. Calculate snow, firn, and ice albedo
% Uses evaporation_condensation and melt_surface from the previous time step
[albedo, albedo_diffuse] = ...
    calculate_albedo(temperature, dz, density, water, grain_radius, albedo, albedo_diffuse, evaporation_condensation, melt_surface, ClimateForcingStep, ModelParam);

% 3. Determine distribution of absorbed SW radiation with depth
shortwave_flux = calculate_shortwave_radiation(dz, density, grain_radius, albedo(1), albedo_diffuse(1), ClimateForcingStep, ModelParam);

% 4. Calculate net shortwave [W m-2]
shortwave_net = sum(shortwave_flux);

% 5. Calculate new temperature-depth profile and turbulent heat fluxes [W m-2]
[temperature, longwave_upward, heat_flux_sensible, heat_flux_latent, ghf, evaporation_condensation] = ...
    calculate_temperature(temperature, dz, density, water(1), grain_radius, shortwave_flux, ClimateForcingStep, ModelParam, verbose);

% 4. Calculate net longwave [W m-2]
longwave_net = ClimateForcingStep.longwave_downward - longwave_upward;

% 6. Change in thickness of top cell due to evaporation/condensation
% Assuming same density as top cell
% ## NEED TO FIX THIS IN CASE ALL OR MORE OF CELL EVAPORATES ##
dz(1) = dz(1) + evaporation_condensation / density(1);

if verbose
    E_evaporation_condensation = evaporation_condensation * temperature(1) * C_ice;
end

% 7. Add snow/rain to top grid cell adjusting cell depth, temperature, and density
[temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, rain] = ...
    calculate_accumulation(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, ClimateForcingStep, ModelParam, verbose);

% 8. Melt and wet compaction
% Calculate water production melt [kg m-2], runoff [kg m-2], and resulting changes
densification_from_melt = sum(dz); % Track thickness before melt

[temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, melt, melt_surface, runoff, refreeze] = ...
    calculate_melt(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, rain, ModelParam, verbose);

densification_from_melt = (densification_from_melt - sum(dz)); % Calculate wet compaction

% 9. Manage the layering to match user defined requirements
[temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, mass_added, E_added] = ...
    manage_layers(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, ModelParam, verbose);

% 10. Allow non-melt densification and determine compaction [m]
densification_from_compaction = sum(dz); % Track thickness before densification

[dz, density] = calculate_density(temperature, dz, density, grain_radius, ClimateForcingStep, ModelParam);

densification_from_compaction = ...
    (densification_from_compaction - sum(dz)); % Calculate dry compaction


if verbose
    dt = ClimateForcingStep.dt;

    % calculate total system mass
    M = dz .* density;
    M_total_final = sum(M) + sum(water);        % total mass [kg]
    M_delta       = M_total_final - M_total_initial + runoff - ClimateForcingStep.precipitation - evaporation_condensation - mass_added;

    % check mass conservation
    if abs(M_delta) > 1E-3
        error('total system mass not conserved: M_delta = %0.4f', M_delta)
    end

    % need to account for rain
    E_snow      = ((ClimateForcingStep.precipitation - rain) * (ClimateForcingStep.temperature_air) * C_ice);
    E_rain      = rain * (ClimateForcingStep.temperature_air * C_ice + LF);
    E_runoff    = runoff * (LF + CtoK * C_ice);
    E_thermal   = sum((dz .* density) .* temperature * C_ice);
    E_water     = sum(water .* (LF + CtoK * C_ice));
    E_shortwave = (shortwave_net*dt);
    E_longwave  = (longwave_net*dt);
    E_thf       = ((heat_flux_sensible+heat_flux_latent)*dt);
    E_ghf       = (ghf*dt);

    E_total_final = E_thermal + E_water + E_runoff; % total energy [J] = initial enegy of snow/ice + initial enegy of water
  
    E_used      = E_total_final - E_total_initial;
    E_supplied  = E_shortwave + E_longwave + E_thf + E_snow + E_rain + E_ghf + E_evaporation_condensation + E_added ;
    E_delta     = E_used - E_supplied;

    % check energy conservation
    if abs(E_delta) > 1E-3
        error('total system energy not conserved: E_delta = %0.4f \n E_shortwave = %0.4f, E_longwave = %0.4f, E_thf = %0.4f, E_added = %0.4f', E_delta, E_shortwave, E_longwave, E_thf, E_added)
    end

    % check bottom grid cell temperature is unchanged
    if abs(temperature(end)-T_bottom) > 1E-3
        error('temperature of bottom grid cell changed: original = %0.10g J, updated = %0.10g J',T_bottom,temperature(end))
    end
end

end