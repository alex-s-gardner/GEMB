function [temperature, longwave_upward, heat_flux_sensible, heat_flux_latent, ghf, evaporation_condensation] = ...
    calculate_temperature(temperature, dz, density, water_surface, grain_radius, shortwave_flux, ClimateForcingStep, ModelParam, verbose)
% calculate_temperature computes new temperature profile accounting for energy absorption
% and thermal diffusion.
%
%% Syntax
%
% [temperature, longwave_upward, heat_flux_sensible, heat_flux_latent, ghf, evaporation_condensation] = ...
%    calculate_temperature(temperature, dz, density, water_surface, grain_radius, shortwave_flux, ClimateForcingStep, ModelParam, verbose)
%
%% Description
%
% This function solves the 1D heat transfer equation to update the vertical 
% temperature profile of the snow/firn column. It accounts for:
%
% 1. Surface Energy Balance: Calculates the net energy flux at the surface 
%    boundary, including:
%    * Turbulent Fluxes: Sensible and latent heat fluxes derived from 
%      bulk aerodynamic formulas, using roughness lengths from Bougamont (2005) 
%      and Foken (2008).
%    * Radiative Fluxes: Incoming/outgoing longwave radiation and surface 
%      shortwave absorption.
% 2. Subsurface Physics: 
%    * Thermal Diffusion: Solves the discretized heat equation using a 
%      finite-volume scheme (Patankar, 1980).
%    * Shortwave Penetration: Adds absorbed shortwave energy (calculated 
%      externally) as a source term at depth.
%    * Thermal Conductivity: Updates conductivity based on density and 
%      temperature (Sturm, 1997).
%
% The function iterates over sub-time steps determined by the stability 
% criterion (Von Neumann stability analysis) to ensure numerical stability 
% of the explicit diffusion scheme.
%
%% Inputs
%
%  temperature              : K            Grid cell temperature (vector).
%  dz                       : m            Grid cell thickness (vector).
%  density                  : kg m^-3      Grid cell density (vector).
%  water_surface            : kg m^-2      Surface water content.
%  grain_radius             : mm           Grain radius (vector).
%  shortwave_flux           : W m^-2       Absorbed shortwave radiation flux per layer.
%  ClimateForcingStep       : struct       Forcing data for the current time step:
%    .longwave_downward     : W m^-2       Downward longwave radiation flux.
%    .temperature_air       : K            2m air temperature.
%    .wind_speed            : m s^-1       Wind velocity.
%    .vapor_pressure        : Pa           Vapor pressure.
%    .pressure_air          : Pa           Air pressure.
%    .dt                    : s            Time step duration.
%  ModelParam               : struct       Model parameters:
%    .density_ice           : kg m^-3      Ice density (threshold for roughness changes).
%    .surface_roughness_effective_ratio : - Ratio of scalar roughness (zT, zQ) to momentum roughness (z0).
%    .dt_divisors           : -            Factors used to find stable time steps.
%    .emissivity_method     : string       Method for calculating emissivity ("uniform", "re_threshold", etc.).
%    .emissivity            : -            Base longwave emissivity.
%    .emissivity_grain_radius_threshold : mm         Grain radius threshold for emissivity switching.
%    .emissivity_grain_radius_large   : -            Emissivity for large grain sizes.
%  verbose                  : logical      Flag to enable energy conservation checks.
%
%% Outputs
%
%  temperature              : K            Updated grid cell temperature (vector).
%  longwave_upward          : W m^-2       Mean upward longwave radiation flux.
%  heat_flux_sensible       : W m^-2       Mean sensible heat flux.
%  heat_flux_latent         : W m^-2       Mean latent heat flux.
%  ghf                      : W m^-2       Mean ground heat flux (flux at bottom boundary).
%  evaporation_condensation : kg m^-2      Cumulative evaporation (+) / condensation (-) mass.
%
%% References
%
% Physics implementations based on:
% Bougamont, M., et al. (2005). (Surface roughness).
% Foken, T. (2008). Micrometeorology. (Roughness lengths).
% Patankar, S. V. (1980). Numerical Heat Transfer and Fluid Flow. (Discretization).
% Sturm, M., et al. (1997). (Thermal conductivity).
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023.

%% INITIALIZE

C_ice = 2102;       % heat capacity of snow/ice (J kg-1 k-1)
RC  = 8.314;        % gas constant [mol-1 K-1]
SB = 5.67E-8;       % Stefan-Boltzmann constant [W m-2 K-4]

d_tolerance    = 1e-11;
T_tolerance    = 1e-4;
W_tolerance    = 1e-13;
ds             = density(1);      % density of top grid cell

% calculated air density [kg/m3]
density_air = 0.029 * ClimateForcingStep.pressure_air / (RC * ClimateForcingStep.temperature_air);

% thermal capacity of top grid cell [J/k]
TCs = density(1) * dz(1) * C_ice;

% determine grid point 'center' vector size
m = length(density);
if m == 0
    error('column has not no gridcells: length(density) = 0')
end

% initialize Evaporation - Condensation
longwave_upward_cumulative = 0.0;
EC_cumulative  = 0.0;
lhf_cumulative = 0.0;
shf_cumulative = 0.0;
ghf_cumulative = 0.0;
CtoK           = 273.15;

if verbose
    T_bottom = temperature(end);
end

%% SURFACE ROUGHNESS (Bougamont, 2005)

% wind/temperature surface roughness height [m]
if (ds < (ModelParam.density_ice - d_tolerance)) && (water_surface < W_tolerance)
    z0 = 0.00012;       % 0.12 mm for dry snow
elseif ds >= (ModelParam.density_ice - d_tolerance)
    z0 = 0.0032;        % 3.2 mm for ice
else
    z0 = 0.0013;        % 1.3 mm for wet snow
end

% determine emissivity
[emissivity, emissivity_melt_switch] = emissivity_initialize(grain_radius(1), ModelParam);

% zT and zQ are percentage of z0
zT     = z0 * ModelParam.surface_roughness_effective_ratio;
zQ     = z0 * ModelParam.surface_roughness_effective_ratio;

% if ClimateForcingStep.wind_speed = 0 goes to infinity therfore if ClimateForcingStep.wind_speed = 0 change
ClimateForcingStep.wind_speed(ClimateForcingStep.wind_speed < 0.01-d_tolerance) = 0.01;

%% THERMAL CONDUCTIVITY (Sturm, 1997: J. Glaciology)

% calculate new thermal conductivity (K) profile [W m-1 K-1]
K = thermal_conductivity(temperature, density, ModelParam);

%% THERMAL DIFFUSION COEFFICIENTS
% A discretization scheme which truncates the Taylor-Series expansion
% after the 3rd term is used. See Patankar 1980, Ch. 3&4
%
% discretized heat equation:
%
%                 Tp = (Au*Tuo+ Ad*Tdo+ (Ap-Au-Ad)Tpo+ S) / Ap
%
% where neighbor coefficients Au, Ap, & Ad are
%
%                   Au = [dz_u/2KU + dz/2KP]^-1
%                   Ad = [dz_d/2KD + dz/2KP]^-1
%                   Ap = density*C_ice*dz/Dt
%
% and u & density represent grid points up and down from the center grid point
% point p and o identifies previous time step values. S is a source term.

% u, density, and p conductivities
KU = [NaN   ; K(1:m-1)];
KD = [K(2:m); NaN];
KP =  K;

% determine u, density & p cell widths
dzU = [NaN    ; dz(1:m-1)];
dzD = [dz(2:m); NaN];

% find stable dt for thermodynamics loop 
dt = thermo_optimal_dt(dz, density, C_ice, K, ModelParam.dt_divisors);

% determine mean (harmonic mean) of K/dz for u, density, & p
Au = (dzU./(2*KU) + dz./(2*KP)).^(-1);
Ad = (dzD./(2*KD) + dz./(2*KP)).^(-1);
Ap = (density.*dz*C_ice)/dt;

% Create neighbor arrays for diffusion calculations instead of a
% tridiagonal matrix

% create "neighbor" coefficient matrix
Nu = Au ./ Ap;
Nd = Ad ./ Ap;
Np = 1 - Nu - Nd;

% specify boundary conditions
% Constant Temperature (Dirichlet) boundary condition
Nu(m) = 0;
Np(m) = 1;
Nd(m) = 0;

% zero flux at surface
Nu(1) = 0;         % Disconnect from the node above (Air/Ghost node)
Np(1) = 1 - Nd(1); % Balance the center node to conserve energy (Weights must sum to 1)

%% RADIATIVE FLUXES

% energy supplied by shortwave radiation [J]
sw        = shortwave_flux * dt;

% ensure no sw reaches bottom cell, add any flux to bottom cell to the cell above
sw(end-1) = sw(end-1) + sw(end);
sw(end)   = 0;

% temperature change due to SW
T_delta_sw      = sw ./ (C_ice * density .* dz);

% Upward longwave radiation flux is calculated from the snow surface
% temperature which is set equal to the average temperature of the
% top grid cells.

% energy supplied by downward longwave radiation to the top grid cell [J]
longwave_downward = ClimateForcingStep.longwave_downward * dt;

% temperature change due to longwave_downward_surf
T_delta_longwave_downward = longwave_downward / TCs;

% only update turbulent head flux every thf_trigger_threshold
thf_trigger_threshold = 1 * 60 * 60; % this will update thf every simulation hour
thf_trigger           = thf_trigger_threshold;

%% PREALLOCATE ARRAYS BEFORE LOOP FOR IMPROVED PERFORMANCE

Tu = zeros(m,1);
Td = zeros(m,1);

%% CALCULATE ENERGY SOURCES AND DIFFUSION FOR EVERY TIME STEP [dt]

for i = 1:dt:ClimateForcingStep.dt
    %      PART OF ENERGY CONSERVATION CHECK
    %      store initial temperature
    if verbose
        % total initial heat energy
        E_initial = sum(temperature .* (C_ice * density .* dz));
    end

    % calculate temperature of snow surface (T_surface)
    % when incoming SW radition is allowed to penetrate the surface,
    % the modeled energy balance becomes very sensitive to how T_surface is
    % calculated.  Here, we take the surface temperature to be temperature(1), but
    % note that the estimated enegy balance & melt are significanly
    % less when T_surface is taken as the mean of the x top grid cells (temperature(1) + temperature(2))/2.0.
    T_surface = temperature(1);
    T_surface = min(273.15, T_surface);    % don't allow T_surface to exceed 273.15 K (0 deg C)

    % TURBULENT HEAT FLUX
    if thf_trigger >= thf_trigger_threshold
        [heat_flux_sensible, heat_flux_latent, latent_heat] = turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, ClimateForcingStep);
        thf_trigger = 0;
    end

    lhf_cumulative = lhf_cumulative + heat_flux_latent * dt;
    shf_cumulative = shf_cumulative + heat_flux_sensible * dt;

    % mass loss (-)/accretion(+) due to evaporation/condensation [kg]
    evaporation_condensation = heat_flux_latent / latent_heat * dt;

    % temperature change due turbulent fluxes
    thf = (heat_flux_sensible + heat_flux_latent) * dt;
    T_delta_thf = thf  / TCs;

    % upward longwave radiation
    longwave_upward            = -(SB * T_surface.^4.0 * emissivity) * dt;
    longwave_upward_cumulative = longwave_upward_cumulative - longwave_upward; % accumulated for output
    T_delta_longwave_upward    = longwave_upward / TCs;

    % new grid point temperature

    % SW penetrates surface
    temperature    = temperature    + T_delta_sw;
    temperature(1) = temperature(1) + T_delta_longwave_downward + T_delta_longwave_upward + T_delta_thf;

    % energy flux across lower boundary (energy supplied by underling ice)
    ghf = Ad(end-1) * (temperature(end) - temperature(end-1)) * dt;
    ghf_cumulative = ghf_cumulative + ghf;


    % temperature diffusion

    % Tu: Shift temperature down one step.
    % The first element is the 'Ghost Node' above surface.
    % For Zero Flux, T_ghost = temperature(1).
    Tu(1)     = temperature(1);
    Tu(2:end) = temperature(1:end-1);

    % Td: Shift temperature up one step.
    % The last element is the 'Ghost Node' below the bottom.
    % Since temperature(m) is fixed (Nu(m)=0, Nd(m)=0), this value is unused,
    % but duplicating temperature(end) keeps the vector size correct.
    Td(1:end-1) = temperature(2:end);
    Td(end)     = temperature(end);

    temperature = (Np .* temperature) + (Nu .* Tu) + (Nd .* Td);

    % calculate cumulative evaporation (+)/condensation(-)
    EC_cumulative = EC_cumulative + evaporation_condensation;

    % if emissivity_method == "re_w_threshold" then check if the surface is melting
    if emissivity_melt_switch
        if temperature(1) < (CtoK - T_tolerance)
            emissivity = ModelParam.emissivity;
        else
            emissivity = ModelParam.emissivity_grain_radius_large;
        end
    end

    % CHECK FOR ENERGY (E) CONSERVATION [UNITS: J]

    if verbose
        E_used     = sum(temperature .* (C_ice * density .* dz)) - E_initial;
        E_supplied = sum(sw) + longwave_downward + longwave_upward + thf + ghf;
        E_delta    = E_used - E_supplied;

        E_tolerance = 1e-3;
        if (abs(E_delta) > E_tolerance) || isnan(E_delta)
           
            fprintf('inputs : T_surface = %0.4f K, water_surface = %0.4f kg m-2, re_surface = %0.04f mm, shortwave_flux = %0.4f W m-2, longwave_downwardf = %0.4f W m-2, ClimateForcingStep.temperature_air = %0.4f K, ClimateForcingStep.wind_speed = %0.4f m/s, ClimateForcingStep.vapor_pressure = %0.3f Pa, ClimateForcingStep.pressure_air = %0.4f Pa \n', ...
                              temperature(1)               , water_surface               , grain_radius(1)                 , sum(shortwave_flux)         , ClimateForcingStep.longwave_downward             , ClimateForcingStep.temperature_air          , ClimateForcingStep.wind_speed            , ClimateForcingStep.vapor_pressure           , ClimateForcingStep.pressure_air)

            fprintf('internals : sw = %0.10g J, longwave_downward = %0.10g J, longwave_upward = %0.10g J, thf = %0.10g J, ghf = %0.10g J \n', ...
                                 sum(sw)      , longwave_downward           , longwave_upward           , thf           , ghf)

            error('energy not conserved in thermodynamics equations: supplied = %0.10g J, used = %0.10g J',...
                                                                     E_supplied         , E_used)
        end

        if T_bottom ~= temperature(end)
            error('temperature of bottom grid cell changed inside of thermal function: original = %0.10g J, updated = %0.10g J',T_bottom,temperature(end))
        end
    end
    thf_trigger = thf_trigger + dt;
end

heat_flux_latent = lhf_cumulative / ClimateForcingStep.dt; % J -> W/m2
heat_flux_sensible = shf_cumulative / ClimateForcingStep.dt; % J -> W/m2
longwave_upward = longwave_upward_cumulative / ClimateForcingStep.dt; % J -> W/m2
ghf = ghf_cumulative / ClimateForcingStep.dt; % J -> W/m2
evaporation_condensation = EC_cumulative;

end

function [emissivity, emissivity_melt_switch] = emissivity_initialize(re_surface, ModelParam)
gdn_tolerance  = 1e-10;
switch ModelParam.emissivity_method
    case "uniform"
        emissivity = ModelParam.emissivity;
        emissivity_melt_switch = false;
    case "re_threshold"
        if re_surface <= (ModelParam.emissivity_grain_radius_threshold + gdn_tolerance)
            emissivity = ModelParam.emissivity;
        else
            emissivity = ModelParam.emissivity_grain_radius_large;
        end
        emissivity_melt_switch = false;
    case "re_w_threshold"
        
        if re_surface <= (ModelParam.emissivity_grain_radius_threshold + gdn_tolerance) 
            % populate emissivity for first thermo interation then update
            emissivity = ModelParam.emissivity;
            emissivity_melt_switch = true;
        else
            emissivity = ModelParam.emissivity_grain_radius_large;
            emissivity_melt_switch = false; % no need to check further as grain_radius(1) > emissivity_grain_radius_threshold
        end
end
end

%% SUBFUNCTIONS

function dt = thermo_optimal_dt(dz, density, C_ice, K, global_dt_or_dt_divisors)
% thermo_optimal_dt 
%
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

% determine minimum acceptable delta t (diffusion number > 1/2) [s]
% 1. Calculate the theoretical limit for every single grid cell
%    (Using 0.5 is the absolute limit; we usually aim lower for safety)
stability_limit_per_cell = 0.5 * (density .* C_ice .* dz.^2) ./ K;

% 2. Find the bottleneck: The smallest allowable dt in the entire column
max_safe_dt = min(stability_limit_per_cell);

% 3. Apply a Safety Factor (0.9 or 0.8 is standard)
%    This accounts for floating point errors or slight non-linearities in K
dt_target = max_safe_dt * 0.8;

% 4. (Optional) Sanity check to prevent extremely small steps if bad data enters
if dt_target < 1e-4
    warning('Timestep is extremely small (%e). Check for near-zero dz layers.', dt_target);
end

% 5. Fit this target into your input data frequency (global_dt)
%    Find the largest divisor of global_dt that is <= dt_target  
if isscalar(global_dt_or_dt_divisors)
    dt_divisors = fast_divisors(global_dt_or_dt_divisors * 10000)/10000; % ClimateForcingStep.dt is in seconds  
else
    dt_divisors = global_dt_or_dt_divisors;
end

dt = dt_divisors(find(dt_divisors <= dt_target, 1, 'last'));

if isempty(dt)
    warning("thermo dt_target < all dt_divisors, setting thermo == to the smallest dt_divisors... this may make termo diffusion unstable")
    dt = dt_divisors(1); % Fallback to smallest possible step
end

end
