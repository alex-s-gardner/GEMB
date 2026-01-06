function [T, shf_cum, lhf_cum, EC, ulwrf] = ...
    thermo(T, dz, d, W_surface, re, swf, dlwrf, T_air, V, e_air, p_air, dt0, ...
    density_ice, Vz, Tz, emissivity, ulw_delta, emissivity_re_threshold, ...
    emissivity_method, thermal_conductivity_method, verbose)

% thermo computes new temperature profile accounting for energy absorption
% and thermal diffusion.
%
%% Syntax
%
%
%
%% Description
%
%
%
%% Inputs
%
% * T:         grid cell temperature [k]
% * dz:        grid cell depth [m]
% * d:         grid cell density [kg m-3]
% * swf:       shortwave radiation fluxes [W m-2]
% * dlwrf:     downward longwave radiation fluxes [W m-2]
% * T_air:        2 m air temperature
% * V:         wind velocity [m s-1]
% * e_air:      screen level vapor pressure [Pa]
% * W_surface: surface water content [kg]
% * dt0:       time step of input data [s]
% * Vz:        air temperature height above surface [m]
% * Tz:        wind height above surface [m]
%
%% Outputs
%
% * T:     grid cell temperature [k]
% * EC:    evaporation/condensation [kg]
% * ulwrf: upward longwave radiation flux [W m-2]
%
%% Documentation
%
% For complete documentation, see: https://github.com/alex-s-gardner/GEMB
%
%% References
% If you use GEMB, please cite the following:
%
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass
% Balance (GEMB): a model of firn processes for cryosphere research, Geosci.
% Model Dev., 16, 2277–2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023.

%% INITIALIZE
CI = 2102;          % heat capacity of snow/ice (J kg-1 k-1)
R  = 8.314;         % gas constant [mol-1 K-1]
SB = 5.67E-8;       % Stefan-Boltzmann constant [W m-2 K-4]

d_tolerance    = 1e-11;
gdn_tolerance  = 1e-10;
W_tolerance    = 1e-13;

ds = d(1);      % density of top grid cell

% calculated air density [kg/m3]
density_air = 0.029 * p_air / (R * T_air);

% thermal capacity of top grid cell [J/k]
TCs = d(1) * dz(1) * CI;

% determine grid point 'center' vector size
m = length(d);
if m == 0
    error('column has not no gridcells: length(d) = 0')
end

% initialize Evaporation - Condensation
EC      = 0.0;
ulwrf   = 0.0;
lhf_cum = 0.0;
shf_cum = 0.0;

if verbose
    T_bottom = T(end);
end

%% SURFACE ROUGHNESS (Bougamont, 2005)
% wind/temperature surface roughness height [m]
if (ds < (density_ice - d_tolerance))  && (W_surface < W_tolerance)
    z0 = 0.00012;       % 0.12 mm for dry snow
elseif ds >= (density_ice - d_tolerance)
    z0 = 0.0032;        % 3.2 mm for ice
else
    z0 = 0.0013;        % 1.3 mm for wet snow
end

% zT and zQ are percentage of z0 (Foken 2008)
zratio = 10.0;
zT     = z0 / zratio;
zQ     = z0 / zratio;

% if V = 0 goes to infinity therfore if V = 0 change
V(V < 0.01-d_tolerance) = 0.01;

%% THERMAL CONDUCTIVITY (Sturm, 1997: J. Glaciology)
% calculate new thermal conductivity (K) profile [W m-1 K-1]
K = thermal_conductivity(T, d, density_ice, thermal_conductivity_method);

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
KU = [NaN ; K(1:m-1)];
KD = [K(2:m) ; NaN];
KP = K;

% determine u, d & p cell widths
dzU = [NaN; dz(1:m-1)];
dzD = [dz(2:m) ; NaN];

% determine minimum acceptable delta t (diffusion number > 1/2) [s]
% 1. Calculate the theoretical limit for every single grid cell
%    (Using 0.5 is the absolute limit; we usually aim lower for safety)
stability_limit_per_cell = 0.5 * (d .* CI .* dz.^2) ./ K;

% 2. Find the bottleneck: The smallest allowable dt in the entire column
max_safe_dt = min(stability_limit_per_cell);

% 3. Apply a Safety Factor (0.9 or 0.8 is standard)
%    This accounts for floating point errors or slight non-linearities in K
dt_target = max_safe_dt * 0.8;

% 4. (Optional) Sanity check to prevent extremely small steps if bad data enters
if dt_target < 1e-4
    warning('Timestep is extremely small (%e). Check for near-zero dz layers.', dt_target);
end

% 5. Fit this target into your input data frequency (dt0)
%    Find the largest divisor of dt0 that is <= dt_target
%    (Your existing divisor logic works well here)

if rem(dt0,1) ~= 0
    warning('rounding dt0 as it is not an exact integer: dt0 = %0.4f', dt0)
    dt0 = round(dt0);
end

f = divisors(dt0 * 10000) / 10000; % assuming dt0 is in seconds
dt = f(find(f <= dt_target, 1, 'last'));

if isempty(dt)
    dt = f(1); % Fallback to smallest possible step
end

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
sw = swf * dt;

% temperature change due to SW
dT_sw = sw ./ (CI * d .* dz);

% Upward longwave radiation flux is calculated from the snow surface
% temperature which is set equal to the average temperature of the
% top grid cells.

% energy supplied by downward longwave radiation to the top grid cell [J]
dlw = dlwrf * dt;

% temperature change due to dlw_surf
dT_dlw = dlw / TCs;

%% PREALLOCATE ARRAYS BEFORE LOOP FOR IMPROVED PERFORMANCE
Tu = zeros(m,1);
Td = zeros(m,1);

%% CALCULATE ENERGY SOURCES AND DIFFUSION FOR EVERY TIME STEP [dt]
for i = 1:dt:dt0
    %     % PART OF ENERGY CONSERVATION CHECK
    %     % store initial temperature
    if verbose
        % total initial heat energy
        E_init = sum(T .* (CI * d .* dz));

        % energy flux across lower boundary (energy supplied by underling ice)
        base_flux = Ad(end-1) * (T(end) - T(end-1)) * dt;
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
    [shf, lhf, L] = turbulent_heat_flux(T_air, T_surface, p_air, e_air, ...
        V, density_air, Vz, Tz, z0, zT, zQ);

    % mass loss (-)/accretion(+) due to evaporation/condensation [kg]
    EC_day = lhf * 86400 / L;

    % temperature change due turbulent fluxes
    turb = (shf + lhf) * dt;
    dT_turb = turb  / TCs;

    % If user wants to directly set emissivity, or grain radius is larger than the
    % threshold, or emissivity_method is 2 and we have wet snow or ice, use prescribed emissivity
    if (emissivity_method == 0) || ((emissivity_re_threshold - re(1)) <= gdn_tolerance)  || ((emissivity_method == 2) && (z0 > 0.001+gdn_tolerance))
        emissivity = emissivity;
    else
        emissivity = 1.0;
    end

    ulw    = -(SB * T_surface.^4.0 * emissivity + ulw_delta) * dt;

    ulwrf  = ulwrf - ulw/dt0;
    dT_ulw = ulw / TCs;

    % new grid point temperature

    % SW penetrates surface
    T    = T    + dT_sw;
    T(1) = T(1) + dT_dlw + dT_ulw + dT_turb;

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
    EC = EC + (EC_day / 86400) * dt;

    lhf_cum = lhf_cum + lhf * dt / dt0;
    shf_cum = shf_cum + shf * dt / dt0;

    %% CHECK FOR ENERGY (E) CONSERVATION [UNITS: J]
    if verbose
        E_used = sum(T .* (CI * d .* dz)) - E_init;
        E_sup = sum(sw) + dlw + ulw + turb + base_flux;
        E_delta = E_used - E_sup;

        if abs(E_delta) > 1E-4 || isnan(E_delta)
            fprintf('sw = %0.10g J, dlw = %0.10g J, ulw = %0.10g J, turb = %0.10g J, base_flux = %0.10g J \n', sum(sw) , dlw , ulw , turb , base_flux)
            error('energy not conserved in thermodynamics equations: supplied = %0.10g J, used = %0.10g J', E_sup, E_used)
        end

        if T_bottom ~= T(end)
            error('temperature of bottom grid cell changed inside of thermal function: original = %0.10g J, updated = %0.10g J',T_bottom,T(end))
        end
    end
end