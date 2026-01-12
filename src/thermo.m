function [T, ulw, shf, lhf, ghf, EC] = ...
    thermo(T, dz, d, water_surface, re, swf, ClimateForcingStep, ModelParam, verbose)
% thermo computes new temperature profile accounting for energy absorption
% and thermal diffusion.
%
%% Syntax
%
% [T, ulw, shf, lhf, ghf, EC] = ...
%    thermo(T, dz, d, water_surface, re, swf, ClimateForcingStep, ModelParam, verbose)
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
%  T                        : K            Grid cell temperature (vector).
%  dz                       : m            Grid cell thickness (vector).
%  d                        : kg m^-3      Grid cell density (vector).
%  water_surface            : kg m^-2      Surface water content.
%  re                       : mm           Grain radius (vector).
%  swf                      : W m^-2       Absorbed shortwave radiation flux per layer.
%  ClimateForcingStep       : struct       Forcing data for the current time step:
%    .dlw                   : W m^-2       Downward longwave radiation flux.
%    .T_air                 : K            2m air temperature.
%    .V                     : m s^-1       Wind velocity.
%    .e_air                 : Pa           Vapor pressure.
%    .p_air                 : Pa           Air pressure.
%    .dt                    : s            Time step duration.
%  ModelParam               : struct       Model parameters:
%    .density_ice           : kg m^-3      Ice density (threshold for roughness changes).
%    .surface_roughness_effective_ratio : - Ratio of scalar roughness (zT, zQ) to momentum roughness (z0).
%    .dt_divisors           : -            Factors used to find stable time steps.
%    .emissivity_method     : string       Method for calculating emissivity ("uniform", "re_threshold", etc.).
%    .emissivity            : -            Base longwave emissivity.
%    .emissivity_re_threshold : mm         Grain radius threshold for emissivity switching.
%    .emissivity_re_large   : -            Emissivity for large grain sizes.
%  verbose                  : logical      Flag to enable energy conservation checks.
%
%% Outputs
%
%  T                        : K            Updated grid cell temperature (vector).
%  ulw                      : W m^-2       Mean upward longwave radiation flux.
%  shf                      : W m^-2       Mean sensible heat flux.
%  lhf                      : W m^-2       Mean latent heat flux.
%  ghf                      : W m^-2       Mean ground heat flux (flux at bottom boundary).
%  EC                       : kg m^-2      Cumulative evaporation (+) / condensation (-) mass.
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
CI = 2102;          % heat capacity of snow/ice (J kg-1 k-1)
RC  = 8.314;         % gas constant [mol-1 K-1]
SB = 5.67E-8;       % Stefan-Boltzmann constant [W m-2 K-4]

d_tolerance    = 1e-11;
T_tolerance    = 1e-4;
W_tolerance    = 1e-13;
ds             = d(1);      % density of top grid cell

% calculated air density [kg/m3]
density_air = 0.029 * ClimateForcingStep.p_air / (RC * ClimateForcingStep.T_air);

% thermal capacity of top grid cell [J/k]
TCs = d(1) * dz(1) * CI;

% determine grid point 'center' vector size
m = length(d);
if m == 0
    error('column has not no gridcells: length(d) = 0')
end

% initialize Evaporation - Condensation
ulw_cumulative = 0.0;
EC_cumulative  = 0.0;
lhf_cumulative = 0.0;
shf_cumulative = 0.0;
ghf_cumulative = 0.0;
CtoK           = 273.15;

if verbose
    T_bottom = T(end);
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
[emissivity, emissivity_melt_switch] = emissivity_initialize(re(1), ModelParam);

% zT and zQ are percentage of z0
zT     = z0 * ModelParam.surface_roughness_effective_ratio;
zQ     = z0 * ModelParam.surface_roughness_effective_ratio;

% if ClimateForcingStep.V = 0 goes to infinity therfore if ClimateForcingStep.V = 0 change
ClimateForcingStep.V(ClimateForcingStep.V < 0.01-d_tolerance) = 0.01;

%% THERMAL CONDUCTIVITY (Sturm, 1997: J. Glaciology)
% calculate new thermal conductivity (K) profile [W m-1 K-1]
K = thermal_conductivity(T, d, ModelParam);

%% THERMAL DIFFUSION COEFFICIENTS

% A discretization scheme which truncates the Taylor-Series expansion
% after the 3rd term is used. See Patankar 1980, Ch. 3&4

% discretized heat equation:

%                 Tp = (Au*Tuo+ Ad*Tdo+ (Ap-Au-Ad)Tpo+ S) / Ap

% where neighbor coefficients Au, Ap, & Ad are

%                   Au = [dz_u/2KU + dz/2KP]^-1
%                   Ad = [dz_d/2KD + dz/2KP]^-1
%                   Ap = d*CI*dz/Dt

% and u & d represent grid points up and down from the center grid point
% point p and o identifies previous time step values. S is a source term.

% u, d, and p conductivities
KU = [NaN   ; K(1:m-1)];
KD = [K(2:m); NaN];
KP =  K;

% determine u, d & p cell widths
dzU = [NaN    ; dz(1:m-1)];
dzD = [dz(2:m); NaN];

% find stable dt for thermodynamics loop 
dt = thermo_optimal_dt(dz, d, CI, K, ModelParam.dt_divisors);

% determine mean (harmonic mean) of K/dz for u, d, & p
Au = (dzU./(2*KU) + dz./(2*KP)).^(-1);
Ad = (dzD./(2*KD) + dz./(2*KP)).^(-1);
Ap = (d.*dz*CI)/dt;

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
sw        = swf * dt;

% ensure no sw reaches bottom cell, add any flux to bottom cell to the cell above
sw(end-1) = sw(end-1) + sw(end);
sw(end)   = 0;

% temperature change due to SW
T_delta_sw      = sw ./ (CI * d .* dz);

% Upward longwave radiation flux is calculated from the snow surface
% temperature which is set equal to the average temperature of the
% top grid cells.

% energy supplied by downward longwave radiation to the top grid cell [J]
dlw = ClimateForcingStep.dlw * dt;

% temperature change due to dlw_surf
T_delta_dlw = dlw / TCs;

% only update turbulent head flux every thf_trigger_threshold
thf_trigger_threshold = 1 * 60 * 60; % this will update thf every simulation hour
thf_trigger           = thf_trigger_threshold;

%% PREALLOCATE ARRAYS BEFORE LOOP FOR IMPROVED PERFORMANCE
Tu = zeros(m,1);
Td = zeros(m,1);

%% CALCULATE ENERGY SOURCES AND DIFFUSION FOR EVERY TIME STEP [dt]
for i = 1:dt:ClimateForcingStep.dt
    %     % PART OF ENERGY CONSERVATION CHECK
    %     % store initial temperature
    if verbose
        % total initial heat energy
        E_initial = sum(T .* (CI * d .* dz));
    end

    % calculate temperature of snow surface (T_surface)
    % when incoming SW radition is allowed to penetrate the surface,
    % the modeled energy balance becomes very sensitive to how T_surface is
    % calculated.  Here, we take the surface temperature to be T(1), but
    % note that the estimated enegy balance & melt are significanly
    % less when T_surface is taken as the mean of the x top grid cells (T(1) + T(2))/2.0.
    T_surface = T(1);
    T_surface = min(273.15, T_surface);    % don't allow T_surface to exceed 273.15 K (0 deg C)

    % TURBULENT HEAT FLUX
    if thf_trigger >= thf_trigger_threshold
        [shf, lhf, latent_heat] = turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, ClimateForcingStep);
        thf_trigger = 0;
    end
    lhf_cumulative = lhf_cumulative + lhf * dt;
    shf_cumulative = shf_cumulative + shf * dt;

    % mass loss (-)/accretion(+) due to evaporation/condensation [kg]
    EC = lhf / latent_heat * dt;

    % temperature change due turbulent fluxes
    thf = (shf + lhf) * dt;
    T_delta_thf = thf  / TCs;

    % upward longwave radiation
    ulw            = -(SB * T_surface.^4.0 * emissivity) * dt;
    ulw_cumulative = ulw_cumulative - ulw; % accumulated for output
    T_delta_ulw    = ulw / TCs;

    % new grid point temperature

    % SW penetrates surface
    T    = T    + T_delta_sw;
    T(1) = T(1) + T_delta_dlw + T_delta_ulw + T_delta_thf;

    % energy flux across lower boundary (energy supplied by underling ice)
    ghf = Ad(end-1) * (T(end) - T(end-1)) * dt;
    ghf_cumulative = ghf_cumulative + ghf;


    % temperature diffusion

    % Tu: Shift T down one step.
    % The first element is the 'Ghost Node' above surface.
    % For Zero Flux, T_ghost = T(1).
    Tu(1)     = T(1);
    Tu(2:end) = T(1:end-1);

    % Td: Shift T up one step.
    % The last element is the 'Ghost Node' below the bottom.
    % Since T(m) is fixed (Nu(m)=0, Nd(m)=0), this value is unused,
    % but duplicating T(end) keeps the vector size correct.
    Td(1:end-1) = T(2:end);
    Td(end)     = T(end);

    T = (Np .* T) + (Nu .* Tu) + (Nd .* Td);

    % calculate cumulative evaporation (+)/condensation(-)
    EC_cumulative = EC_cumulative + EC;

    % if emissivity_method == "re_w_threshold" then check if the surface is melting
    if emissivity_melt_switch
        if T(1) < (CtoK - T_tolerance)
            emissivity = ModelParam.emissivity;
        else
            emissivity = ModelParam.emissivity_re_large;
        end
    end

    %% CHECK FOR ENERGY (E) CONSERVATION [UNITS: J]
    if verbose
        E_used     = sum(T .* (CI * d .* dz)) - E_initial;
        E_supplied = sum(sw) + dlw + ulw + thf + ghf;
        E_delta    = E_used - E_supplied;

        E_tolerance = 1e-3;
        if (abs(E_delta) > E_tolerance) || isnan(E_delta)
           
            fprintf('inputs : T_surface = %0.4f K, water_surface = %0.4f kg m-2, re_surface = %0.04f mm, swf = %0.4f W m-2, dlwf = %0.4f W m-2, ClimateForcingStep.T_air = %0.4f K, ClimateForcingStep.V = %0.4f m/s, ClimateForcingStep.e_air = %0.3f Pa, ClimateForcingStep.p_air = %0.4f Pa \n', ...
                              T(1)               , water_surface               , re(1)                 , sum(swf)         , ClimateForcingStep.dlw             , ClimateForcingStep.T_air          , ClimateForcingStep.V            , ClimateForcingStep.e_air           , ClimateForcingStep.p_air)

            fprintf('internals : sw = %0.10g J, dlw = %0.10g J, ulw = %0.10g J, thf = %0.10g J, ghf = %0.10g J \n', ...
                                 sum(sw)      , dlw           , ulw           , thf           , ghf_flux)

            error('energy not conserved in thermodynamics equations: supplied = %0.10g J, used = %0.10g J',...
                                                                     E_supplied         , E_used)
        end

        if T_bottom ~= T(end)
            error('temperature of bottom grid cell changed inside of thermal function: original = %0.10g J, updated = %0.10g J',T_bottom,T(end))
        end
    end
    thf_trigger = thf_trigger + dt;
end

lhf = lhf_cumulative / ClimateForcingStep.dt; % J -> W/m2
shf = shf_cumulative / ClimateForcingStep.dt; % J -> W/m2
ulw = ulw_cumulative / ClimateForcingStep.dt; % J -> W/m2
ghf = ghf_cumulative / ClimateForcingStep.dt; % J -> W/m2
EC = EC_cumulative;
end

function [emissivity, emissivity_melt_switch] = emissivity_initialize(re_surface, ModelParam)
gdn_tolerance  = 1e-10;
switch ModelParam.emissivity_method
    case "uniform"
        emissivity = ModelParam.emissivity;
        emissivity_melt_switch = false;
    case "re_threshold"
        if re_surface <= (ModelParam.emissivity_re_threshold + gdn_tolerance)
            emissivity = ModelParam.emissivity;
        else
            emissivity = ModelParam.emissivity_re_large;
        end
        emissivity_melt_switch = false;
    case "re_w_threshold"
        
        if re_surface <= (ModelParam.emissivity_re_threshold + gdn_tolerance) 
            % populate emissivity for first thermo interation then update
            emissivity = ModelParam.emissivity;
            emissivity_melt_switch = true;
        else
            emissivity = ModelParam.emissivity_re_large;
            emissivity_melt_switch = false; % no need to check further as re(1) > emissivity_re_threshold
        end
end
end
