function [T, dz, d, W, re, gdn, gsp, a, adiff, EC, M_surf, sw_net, shf, ...
    lhf, ulw, Ra, M, R, F, mass_added, energy_added, compaction_dens, compaction_melt] = ...
    gemb_core(T, dz, d, W, re, gdn, gsp, a, adiff, dt, P, EC, M_surf, ...
    density_ice, black_carbon_snow, black_carbon_ice, solar_zenith_angle, cloud_optical_thickness, cloud_fraction, dsw, ...
    dsw_diffuse, dlw, Ta, V, eAir, pAir, S, verbose)

    % GEMB_STEP Performs a single time-step of the GEMB model.
    %   Calculates grain growth, albedo, radiative transfer, thermodynamics, 
    %   accumulation, melt, layer management, and densification.

    % 1. Snow grain metamorphism
    % [always calculate, used in thermo and albedo]
    [re, gdn, gsp] = grainGrowth(T, dz, d, W, re, gdn, gsp, dt, S.albedo_method);

    % 2. Calculate snow, firn, and ice albedo
    % Uses EC and M_surf from the previous time step
    [a, adiff] = albedo(T, dz, d, W, re, a, adiff, dt, P, EC, M_surf, density_ice, ...
        black_carbon_snow, black_carbon_ice, solar_zenith_angle, cloud_optical_thickness, cloud_fraction, S.albedo_method, ...
        S.albedo_ice, S.albedo_snow, S.albedo_fixed, S.albedo_desnity_threshold, S.albedo_wet_snow_t0, S.albedo_dry_snow_t0, S.albedo_K);

    % 3. Determine distribution of absorbed SW radiation with depth
    swf = shortwave(S.sw_absorption_method, S.albedo_method, dsw, dsw_diffuse, a(1), adiff(1), d, dz, re, density_ice);

    % 4. Calculate net shortwave [W m-2]
    sw_net = sum(swf);

    % 5. Calculate new temperature-depth profile and turbulent heat fluxes [W m-2]
    [T, shf, lhf, EC, ulw] = thermo(T, dz, d, W(1), re, dt, swf, dlw, Ta, V, eAir, pAir, density_ice, ...
        S.thermal_conductivity_method, S.emissivity_method, S.emissivity, S.ulw_delta, S.emissivity_re_threshold, S.Vz, S.Tz, verbose);

    % 6. Change in thickness of top cell due to evaporation/condensation
    % Assuming same density as top cell
    % ## NEED TO FIX THIS IN CASE ALL OR MORE OF CELL EVAPORATES ##
    dz(1) = dz(1) + EC / d(1);

    % 7. Add snow/rain to top grid cell adjusting cell depth, temperature, and density
    [T, dz, d, W, re, gdn, gsp, a, adiff, Ra] = ...
        accumulation(T, dz, d, W, re, gdn, gsp, a, adiff, Ta, P, V, density_ice, ...
        S.albedo_method, S.new_snow_method, S.T_mean, S.column_dzmin, S.C, S.V_mean, S.albedo_snow);

    % 8. Melt and wet compaction
    % Calculate water production M [kg m-2], runoff R [kg m-2], and resulting changes
    compaction_melt = sum(dz); % Track thickness before melt
    [T, dz, d, W, re, gdn, gsp, a, adiff, M, M_surf, R, F] = ...
        melt(T, dz, d, W, re, gdn, gsp, a, adiff, Ra, density_ice, verbose);
    compaction_melt = (compaction_melt - sum(dz)); % Calculate wet compaction

    % 9. Manage the layering to match user defined requirements
    [T, dz, d, W, re, gdn, gsp, a, adiff, mass_added, energy_added] = ...
        managelayers(T, dz, d, W, re, gdn, gsp, a, adiff, S.column_dzmin, S.column_zmax, S.column_zmin, S.column_ztop, S.column_zy, verbose);

    % 10. Allow non-melt densification and determine compaction [m]
    compaction_dens = sum(dz); % Track thickness before densification
    [dz, d] = densification(T, dz, d, re, dt, density_ice, S.albedo_method, S.densification_method, S.T_mean, S.C, S.sw_absorption_method, S.albedo_desnity_threshold);
    compaction_dens = (compaction_dens - sum(dz)); % Calculate dry compaction

end