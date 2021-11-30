function [T, EC] = thermo(T, dz, d, swf, dlwrf, Ta, V, eAir, pAir, Ws, dt0, ...
    Vz, Tz)

%% ENGLACIAL THERMODYNAMICS
 
% Description: 
% computes new temperature profile accounting for energy absorption and 
% thermal diffusion.

%% INPUTS
% * T: grid cell temperature [k]
% * dz: grid cell depth [m]
% * d: grid cell density [kg m-3]
% * swf: shortwave radiation fluxes [W m-2]
% * dlwrf: downward longwave radiation fluxes [W m-2]
% * Ta: 2 m air temperature
% * V:  wind velocity [m s-1]
% * eAir: screen level vapor pressure [Pa]
% * Ws: surface water content [kg]
% * dt0: time step of input data [s]
% * elev: surface elevation [m a.s.l.] 
% * Vz: air temperature height above surface [m]
% * Tz: wind height above surface [m]

%% OUTPUTS
% * T: grid cell temperature [k]
% * EC: evaporation/condensation [kg]

%% INITIALIZE

CI = 2102;          % heat capacity of snow/ice (J kg-1 k-1)
dtScaling = 1/11;   % This is a necessary fudge factor to avoid instability errors in diffusion [Nicole discovery]
% CA = 1005;        % heat capacity of air (J kg-1 k-1)
% LF = 0.3345E6;    % latent heat of fusion(J kg-1)
% LV = 2.495E6;     % latent heat of vaporization(J kg-1)
% dIce = 910;       % density of ice [kg m-3]
% dSnow = 300;      % density of snow [kg m-3]
SB = 5.67E-8;       % Stefan-Boltzmann constant [W m-2 K-4]

ds = d(1);      % density of top grid cell

% calculated air density [kg/m3]
dAir = 0.029 * pAir /(8.314 * Ta);

% theremal capacity of top grid cell [J/k]
TCs = d(1)*dz(1)*CI; 

% determine grid point 'center' vector size
m = length(d);

% initialize Evaporation - Condenstation 
EC = 0;

% check if all SW applied to surface or distributed throught subsurface
% swIdx = length(swf) > 1

%% SURFACE ROUGHNESS (Bougamont, 2005)
% wind/temperature surface roughness height [m]
if ds < 910 && Ws == 0
    z0 = 0.00012;       % 0.12 mm for dry snow
elseif ds >= 910
    z0 = 0.0032;        % 3.2 mm for ice
else
    z0 = 0.0013;        % 1.3 mm for wet snow
end

% if V = 0 goes to infinity therfore if V = 0 change
V(V < 0.01) = 0.01;

% Bulk-transfer coefficient for turbulent fluxes
An =  0.4^2 / log(Tz/z0)^2;     % Bulk-transfer coefficient
C = An * dAir * V;              % shf & lhf common coefficient

%% THERMAL CONDUCTIVITY (Sturm, 1997: J. Glaciology)
% calculate new K profile [W m-1 K-1]

% logical index of snow and firn
sfIdx = d < 910;

% initialize conductivity
K = zeros(m,1);

% for snow and firn (density < 910 kg m-3) (Sturn et al, 1997)
K(sfIdx) = 0.138 - 1.01E-3 * d(sfIdx) + 3.233E-6 * (d(sfIdx).^2);

% for ice (density >= 910 kg m-3)
K(~sfIdx) = 9.828 * exp(-5.7E-3*T(~sfIdx));

%% THERMAL DIFFUSION COEFFICIENTS
 
% A discretization scheme which truncates the Taylor-Series expansion
% after the 3rd term is used. See Patankar 1980, Ch. 3&4
 
% discretized heat equation:
 
%                 Tp = (Au*Tu° + Ad*Td° + (Ap-Au-Ad)Tp° + S) / Ap
 
% where neighbor coefficients Au, Ap, & Ad are
 
%                   Au = [dz_u/2KP + dz_p/2KE]^-1
%                   Ad = [dz_d/2KP + dz_d/2KD]^-1 
%                   Ap = d*CI*dz/Dt 
 
% and u & d represent grid points up and down from the center grid point 
% p and % u & d represent grid points up and down from the center grid 
% point p and ° identifies previous time step values. S is a source term.
 
% u, d, and p conductivities
KU = [NaN ; K(1:m-1)];
KD = [K(2:m) ; NaN];
KP = K;

% determine u, d & p cell widths
dzU = [NaN; dz(1:m-1)];
dzD = [dz(2:m) ; NaN];
 
% determine minimum acceptable delta t (diffusion number > 1/2) [s]
dt = min(CI * dz.^2 .* d  ./ (3 * K)) .* dtScaling;

% smallest possible even integer of 60 min where diffusion number > 1/2
% must go evenly into one hour or the data frequency it it is smaller

% all integer factors of the number of second in a day (8600 [s])
f = [1 2 3 4 5 6 8 9 10 12 15 16 18 20 24 25 30 36 40 45 48 50 60 ...
    72 75 80 90 100 120 144 150 180 200 225 240 300 360 400 450 600 ...
    720 900 1200 1800 3600];

% return the min integer factor that is < dt
dt = max(f(f<dt));

% determine mean (harmonic mean) of K/dz for u, d, & p
Au = (dzU./(2*KP) + dz./(2*KU)).^(-1);
Ad = (dzD./(2*KP) + dz./(2*KD)).^(-1);
Ap = (d.*dz*CI)/dt;

% create "neighbor" coefficient matrix
Nu = Au ./ Ap;
Nd = Ad ./ Ap;
Np = 1 - Nu - Nd;

% specify boundary conditions
% constant flux at bottom
Nu(m) = 0;
Np(m) = 1;

% zero flux at surface
Np(1) = 1 - Nd(1);

% Create neighbor arrays for diffusion calculations instead of a
% tridiagonal matrix
Nu(1) = 0;
Nd(m) = 0;


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
T0 = zeros(m+2,1);

%% CALCULATE ENERGY SOURCES AND DIFFUSION FOR EVERY TIME STEP [dt]
for i = 1:dt:dt0
%     % PART OF ENERGY CONSERVATION CHECK
%     % store initial temperature
%     T_init = T;
    
    % calculate temperature of snow surface (Ts)
    % when incoming SW radition is allowed to penetrate the surface,
    % the modeled energy balance becomes very sensitive to how Ts is
    % calculated.  The estimated enegy balance & melt are significanly
    % less when Ts is taken as the mean of the x top grid cells.
    Ts = (T(1) + T(2))/2;
    Ts = min(273.15,Ts);    % don't allow Ts to exceed 273.15 K (0°C)
    
    % TURBULENT HEAT FLUX
    
    % Monin–Obukhov Stability Correction
    % Reference:
    % Ohmura, A., 1982: Climate and Energy-Balance on the Arctic Tundra.
    % Journal of Climatology, 2, 65-84.
    
    % calculate the Bulk Richardson Number (Ri)
    Ri = (2*9.81* (Vz - z0) * (Ta - Ts)) / ((Ta + Ts)* V^2);
    
    % calculate Monin-–Obukhov stability factors 'coef_M' and 'coef_H'
    
    % do not allow Ri to exceed 0.19
    Ri = min(Ri, 0.19);
    
    % calculate momentum 'coef_M' stability factor
    if Ri > 0
        % if stable
        coefM = 1/(1-5.2*Ri);
    else
        coefM = (1-18*Ri)^-0.25;
    end
    
    % calculate heat/wind 'coef_H' stability factor
    if Ri < -0.03
        coefH = coefM / 1.3; % well done to Nicole for spotting /1.3 instead of *1.3 13/03/2018
    else
        coefH = coefM;
    end
    
    %% Sensible Heat
    % calculate the sensible heat flux [W m-2](Patterson, 1998)
    shf = C * 1005 * (Ta - Ts);
    
    % adjust using Monin–Obukhov stability theory
    shf = shf / (coefM * coefH);
    
    %% Latent Heat
    % determine if snow pack is melting & calcualte surface vapour pressure
    % over ice or liquid water
    if Ts >= 273.15
        L = 2.495E6; % for liquid water at 273.15 k to vapor
        
        % for liquid surface (assume liquid on surface when Ts == 0 deg C)
        % Wright (1997), US Meteorological Handbook from Murphy and Koop,
        % 2005 Apendix A
        eS = 611.21 * exp(17.502 * (Ts - 273.15) / (240.97 + Ts - 273.15));
        
    else
        L = 2.8295E6; % latent heat of sublimation
        
        % for an ice surface Murphy and Koop, 2005 [Equation 7]
        eS = exp(9.550426 - 5723.265/Ts + 3.53068 * log(Ts)...
            - 0.00728332 * Ts);
    end

    % Latent heat flux [W m-2]
    lhf = C * L * (eAir - eS) * 0.622 / pAir;
    
    % adjust using Monin–Obukhov stability theory (if lhf '+' then there is
    % energy and mass gained at the surface, if '-' then there is mass and
    % energy loss at the surface.
    lhf = lhf / (coefM * coefH);
    
    % mass loss (-)/acreation(+) due to evaporation/condensation [kg]
    EC_day = lhf * 86400 / L;
    
    % temperature change due turbulent fluxes
    turb = (shf + lhf)* dt;
    dT_turb = turb  / TCs;
    
    % upward longwave contribution
    ulw = - SB * Ts^4 * dt;
    dT_ulw = ulw / TCs;
    
    % new grid point temperature
    
    % SW penetrates surface
    T = T + dT_sw;
    T(1) = T(1) + dT_dlw + dT_ulw + dT_turb;
    
    % temperature diffusion
    T0(2:m+1) = T;
    Tu = T0(1:m);
    Td = T0(3:m+2);
    
    T = (Np .* T) + (Nu .* Tu) + (Nd .* Td);
    
    % calculate cumulative evaporation (+)/condensation(-)
    EC = EC + (EC_day/86400)*dt;
    
     %% CHECK FOR ENERGY (E) CONSERVATION [UNITS: J]
%     % energy flux across lower boundary (energy supplied by underling ice)
%     base_flux = Ad(end-1)*(T_init(end)-T_init(end-1)) * dt;
%     
%     E_used = sum((T - T_init) .* (d.*dz*CI));
%     E_sup = ((sum(swf)  * dt) + dlw + ulw + turb + base_flux);
%     
%     E_diff = E_used - E_sup;
%     
%     if abs(E_diff) > 1E-6 || isnan(E_diff)
%         disp(T(1))
%         error('energy not conserved in thermodynamics equations')
%     end
end