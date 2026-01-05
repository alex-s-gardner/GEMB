function [T, dz, d, W, re, gdn, gsp, a, adiff, EC, Msurf, netSW, shf, ...
    lhf, ulw, Ra, M, R, F, mAdd, addE, comp1, comp2] = ...
    gemb_step(T, dz, d, W, re, gdn, gsp, a, adiff, dt, P, EC, Msurf, ...
    dIce, ccsnowValue, cciceValue, szaValue, cotValue, cldFrac, dsw, ...
    dswdiffrf, dlw, Ta, V, eAir, pAir, S, verbose)

    % GEMB_STEP Performs a single time-step of the GEMB model.
    %   Calculates grain growth, albedo, radiative transfer, thermodynamics, 
    %   accumulation, melt, layer management, and densification.

    % 1. Snow grain metamorphism
    % [always calculate, used in thermo and albedo]
    [re, gdn, gsp] = grainGrowth(T, dz, d, W, re, gdn, gsp, dt, S.aIdx);

    % 2. Calculate snow, firn, and ice albedo
    % Uses EC and Msurf from the previous time step
    [a, adiff] = albedo(T, dz, d, W, re, a, adiff, dt, P, EC, Msurf, dIce, ...
        ccsnowValue, cciceValue, szaValue, cotValue, cldFrac, S.aIdx, ...
        S.aIce, S.aSnow, S.aValue, S.adThresh, S.t0wet, S.t0dry, S.K);

    % 3. Determine distribution of absorbed SW radiation with depth
    swf = shortwave(S.swIdx, S.aIdx, dsw, dswdiffrf, a(1), adiff(1), d, dz, re, dIce);

    % 4. Calculate net shortwave [W m-2]
    netSW = sum(swf);

    % 5. Calculate new temperature-depth profile and turbulent heat fluxes [W m-2]
    [T, shf, lhf, EC, ulw] = thermo(T, dz, d, W(1), re, dt, swf, dlw, Ta, V, eAir, pAir, dIce, ...
        S.tcIdx, S.eIdx, S.teValue, S.dulwrfValue, S.teThresh, S.dzMin, S.Vz, S.Tz, ...
        S.ThermoDeltaTScaling, S.isdeltaLWup, verbose);

    % 6. Change in thickness of top cell due to evaporation/condensation
    % Assuming same density as top cell
    % ## NEED TO FIX THIS IN CASE ALL OR MORE OF CELL EVAPORATES ##
    dz(1) = dz(1) + EC / d(1);

    % 7. Add snow/rain to top grid cell adjusting cell depth, temperature, and density
    [T, dz, d, W, re, gdn, gsp, a, adiff, Ra] = ...
        accumulation(T, dz, d, W, re, gdn, gsp, a, adiff, Ta, P, V, dIce, ...
        S.aIdx, S.dsnowIdx, S.Tmean, S.dzMin, S.C, S.Vmean, S.aSnow);

    % 8. Melt and wet compaction
    % Calculate water production M [kg m-2], runoff R [kg m-2], and resulting changes
    comp2 = sum(dz); % Track thickness before melt
    [T, dz, d, W, re, gdn, gsp, a, adiff, M, Msurf, R, F] = ...
        melt(T, dz, d, W, re, gdn, gsp, a, adiff, Ra, dIce, verbose);
    comp2 = (comp2 - sum(dz)); % Calculate wet compaction

    % 9. Manage the layering to match user defined requirements
    [T, dz, d, W, re, gdn, gsp, a, adiff, mAdd, addE] = ...
        managelayers(T, dz, d, W, re, gdn, gsp, a, adiff, S.dzMin, S.zMax, S.zMin, S.zTop, S.zY, verbose);

    % 10. Allow non-melt densification and determine compaction [m]
    comp1 = sum(dz); % Track thickness before densification
    [dz, d] = densification(T, dz, d, re, dt, dIce, S.aIdx, S.denIdx, S.Tmean, S.C, S.swIdx, S.adThresh);
    comp1 = (comp1 - sum(dz)); % Calculate dry compaction

end