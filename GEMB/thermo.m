function [shf_cum, lhf_cum, T, EC, ulwrf] = thermo(T, re, dz, d, swf, dlwrf, Ta, V, eAir, pAir, tcIdx, eIdx, ...
        teValue, dulwrfValue, teThresh, Ws, dt0, dzMin, Vz, Tz, dtScaling, dIce, isdeltaLWup)
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
% * Vz: air temperature height above surface [m]
% * Tz: wind height above surface [m]
% 
%% Outputs
% 
% * T: grid cell temperature [k]
% * EC: evaporation/condensation [kg]
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
CA = 1005.0;        % heat capacity of air (J kg-1 K-1)
CtoK = 273.15;      % Kelvin to Celcius conversion/ice melt. point T in K
LS = 2.8295e6;      % latent heat of sublimation (J kg-1)
R  = 8.314;         % gas constant [mol-1 K-1]
% dtScaling = 1/11; % This is a necessary fudge factor to avoid instability errors in diffusion [Nicole discovery]
% CA = 1005;        % heat capacity of air (J kg-1 k-1)
% LF = 0.3345E6;    % latent heat of fusion(J kg-1)
LV = 2.495E6;       % latent heat of vaporization(J kg-1)
% dSnow = 300;      % density of snow [kg m-3]
SB = 5.67E-8;       % Stefan-Boltzmann constant [W m-2 K-4]

Ttol = 1e-10;
Dtol = 1e-11;
Gdntol = 1e-10;
Wtol = 1e-13;

ds = d(1);      % density of top grid cell

% calculated air density [kg/m3]
dAir = 0.029 * pAir /(R * Ta);

% thermal capacity of top grid cell [J/k]
TCs = d(1)*dz(1)*CI; 

% determine grid point 'center' vector size
m = length(d);

% initialize Evaporation - Condensation 
EC      = 0.0;
ulwrf   = 0.0;
lhf_cum = 0.0;
shf_cum = 0.0;

% check if all SW applied to surface or distributed throught subsurface
% swIdx = length(swf) > 1

%% SURFACE ROUGHNESS (Bougamont, 2005)
% wind/temperature surface roughness height [m]
if ds < dIce-Dtol && Ws < Wtol
    z0 = 0.00012;       % 0.12 mm for dry snow
elseif ds >= dIce-Dtol
    z0 = 0.0032;        % 3.2 mm for ice
else
    z0 = 0.0013;        % 1.3 mm for wet snow
end

% zT and zQ are percentage of z0 (Foken 2008)
zratio = 10.0;
zT = z0/zratio;
zQ = z0/zratio;

% if V = 0 goes to infinity therfore if V = 0 change
V(V < 0.01-Dtol) = 0.01;

% Bulk-transfer coefficient for turbulent fluxes
An =  0.4^2; % Bulk-transfer coefficient
C = An*V;  % shf & lhf common coefficient

%% THERMAL CONDUCTIVITY (Sturm, 1997: J. Glaciology)
% calculate new K profile [W m-1 K-1]

% logical index of snow and firn
sfIdx = d < dIce-Dtol;

% initialize conductivity
K = zeros(m,1);

% for snow and firn (density < 910 kg m-3) (Sturm et al, 1997) or (Calonne et al., 2011)
if tcIdx == 2
    K(sfIdx) = 0.024 - 1.23E-4 * d(sfIdx) + 2.5e-6 * (d(sfIdx).^2);
else %default (Sturm et al, 1997)
    K(sfIdx) = 0.138 - 1.01E-3 * d(sfIdx) + 3.233E-6 * (d(sfIdx).^2);
end

% for ice (density >= 910 kg m-3)
K(~sfIdx) = 9.828 * exp(-5.7E-3*T(~sfIdx));

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
if m>0
    dt = min(CI * dz.^2 .* d  ./ (3 * K) .* dtScaling);
else
    dt=1e12;
end

% smallest possible even integer of 60 min where diffusion number > 1/2
% must go evenly into one hour or the data frequency it it is smaller

% all integer factors of the number of second in a day (8600 [s])
f = [1 2 3 4 5 6 8 9 10 12 15 16 18 20 24 25 30 36 40 45 48 50 60 ...
    72 75 80 90 100 120 144 150 180 200 225 240 300 360 400 450 600 ...
    720 900 1200 1800 3600];
fi=1./f(end-12:-1:2);
f = [fi((fi*1e12-floor(fi*1e12)) == 0) f];

% return the min integer factor that is < dt
I=f<dt-Dtol;
if sum(I)
	dt = max(f(I));
else
	dt = f(1);
	display([' WARNING: calculated timestep for thermal loop is < ' num2str(f(1)) ' second. (' num2str(dt) ' sec) ' sprintf('\n')])
end

if m>0
    % determine mean (harmonic mean) of K/dz for u, d, & p
    Au = (dzU./(2*KU) + dz./(2*KP)).^(-1);
    Ad = (dzD./(2*KD) + dz./(2*KP)).^(-1);
    Ap = (d.*dz*CI)/dt;

    % create "neighbor" coefficient matrix
    Nu = Au ./ Ap;
    Nd = Ad ./ Ap;
    Np = 1 - Nu - Nd;
else
    Nu = 0; 
    Nd = 0; 
    Np = 0;
end

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
    % calculated.  Here, we take the surface temperature to be T(1), but
     % note that the estimated enegy balance & melt are significanly
    % less when Ts is taken as the mean of the x top grid cells (T(1) + T(2))/2.0.
    Ts = T(1);
    Ts = min(273.15,Ts);    % don't allow Ts to exceed 273.15 K (0 deg C)
    
    % TURBULENT HEAT FLUX
    
    % Monin-Obukhov Stability Correction
    % Reference:
    % Ohmura, A., 1982: Climate and Energy-Balance on the Arctic Tundra.
    % Journal of Climatology, 2, 65-84.
    
    % calculate the Bulk Richardson Number (Ri)
    Ri = ((100000./pAir).^0.286).*(2.0*9.81*(Ta - Ts)) ./ (Tz.*(Ta + Ts).*(((V/Vz).^2.0)));
    
    a1     = 1.0;
    b1     = 2.0/3.0;
    c1     = 5.0;
    d1     = 0.35;
    PhiMz0 = 0.0;
    PhiHzT = 0.0;
    PhiHzQ = 0.0;

    if (Ri > 0.0+Ttol)
        % if stable
        if (Ri < 0.2-Ttol)
            zL = Ri./(1.0-5.0*Ri);
        else
            zL = Ri;
        end
        
        %zL = min(zL, 0.5); %Sjoblom, 2014
        zLM = max(zL./Vz.*z0,1e-3);
        zLT = max(zL./Tz.*zT,1e-3);
        
        % Ding et al. 2020, from Beljaars and Holtslag (1991)
        PhiMz  = -1.*(a1*zL + b1*(zL-c1/d1)*exp(-1.*d1*zL) + b1*c1/d1);
        PhiHz  = -1.*((1.+2.*a1*zL/3.).^1.5 + b1*(zL-c1/d1)*exp(-1.*d1*zL) + b1*c1/d1 - 1.0);
        PhiMz0 = -1.*(a1*zLM + b1*(zLM-c1/d1)*exp(-1.*d1*zLM) + b1*c1/d1);
        PhiHzT = -1.*((1.+2.*a1*zLT/3.).^1.5 + b1*(zLT-c1/d1)*exp(-1.*d1*zLT) + b1*c1/d1 - 1.0);
        
        PhiHzQ=PhiHzT;
    
    else 
    
        zL  = Ri/1.5; %max(Ri, -0.5+Ttol)/1.5; % Hogstrom (1996)
            
        %Sjoblom, 2014
        xm=(1.0-19.0*zL).^-0.25;
        PhiMz=2.0*log((1.+xm)/2.0) + log((1.+xm.^2)/2.0) - 2.*atan(xm) + pi/2.;
        
        xh=0.95*(1.0-11.6*zL).^(-0.5);
        PhiHz=2.0*log((1.0+xh.^2)/2.0);

    end
    
    coefM  = log(Vz./z0) - PhiMz + PhiMz0; % Ding et al., 2019
    coefHT = log(Tz./zT) - PhiHz + PhiHzT; % Sjoblom, 2014, after Foken 2008
    coefHQ = log(Tz./zQ) - PhiHz + PhiHzQ; % Sjoblom, 2014, after Foken 2008

    %% Sensible Heat
    % calculate the sensible heat flux [W m-2](Patterson, 1998)
    shf = dAir .* C .* CA .* (Ta - Ts) .* (100000./pAir).^0.286;
    
    % adjust using Monin-Obukhov stability theory
    shf = shf./(coefM.*coefHT);
    
    %% Latent Heat
    %   determine if snow pack is melting & calcualte surface vapour pressure over ice or liquid water
    if (Ts >= CtoK-Ttol)
        L = LV; %for liquid water at 273.15 k to vapor
        
        %for liquid surface (assume liquid on surface when Ts == 0 deg C)
        % Wright (1997), US Meteorological Handbook from Murphy and Koop, 2005 Appendix A
        %eS = 611.21 * exp(17.502 * (Ts - CtoK) / (240.97 + Ts - CtoK));
        % Murray 1967, https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
        eS = 610.78 * exp(17.2693882 .* (Ts - CtoK - 0.01) ./ (Ts - 35.86));
    else
        L = LS; % latent heat of sublimation
        
        % for an ice surface Murphy and Koop, 2005 [Equation 7]
        %eS = exp(9.550426 - 5723.265/Ts + 3.53068 * log(Ts) - 0.00728332 * Ts);
        % for an ice surface Ding et al., 2019 after Bolton, 1980
        eS = 610.78 * exp(21.8745584 .* (Ts - CtoK - 0.01) ./ (Ts - 7.66));
    end

    %Latent heat flux [W m-2]
    lhf = C .* L .* (eAir - eS) / (461.9*(Ta+Ts)/2.0);

    % adjust using Monin-Obukhov stability theory (if lhf '+' then there is energy and mass gained at the surface,
    % if '-' then there is mass and energy loss at the surface.
    lhf = lhf./(coefM.*coefHQ);

    % mass loss (-)/accretion(+) due to evaporation/condensation [kg]
    EC_day = lhf * 86400 / L;
    
    % temperature change due turbulent fluxes
    turb = (shf + lhf)* dt;
    dT_turb = turb  / TCs;

    % upward longwave contribution
    deltaULW   = 0.0;
    emissivity = 1.0;

    %If user wants to set a upward long wave bias
    if isdeltaLWup
        deltaULW = dulwrfValue; 
    end

    %If user wants to directly set emissivity, or grain radius is larger than the
    %threshold, or eIdx is 2 and we have wet snow or ice, use prescribed emissivity
    if (eIdx==0 | (teThresh - re(1))<=Gdntol | (eIdx==2 & z0>0.001+Gdntol)) 
        emissivity = teValue; 
    end

    ulw    = - (SB * Ts.^4.0 * emissivity + deltaULW) * dt;
    ulwrf  = ulwrf - ulw/dt0;
    dT_ulw = ulw / TCs;
    
    % new grid point temperature
    
    % SW penetrates surface
    T    = T    + dT_sw;
    T(1) = T(1) + dT_dlw + dT_ulw + dT_turb;
    
    % temperature diffusion
    T0(2:m+1) = T;
    T0(1)     = Ta;
    T0(m+2)   = T(m);
    Tu        = T0(1:m);
    Td        = T0(3:m+2);
    
    T = (Np .* T) + (Nu .* Tu) + (Nd .* Td);

    % calculate cumulative evaporation (+)/condensation(-)
    EC = EC + (EC_day/86400)*dt;

    lhf_cum = lhf_cum+lhf*dt/dt0;
    shf_cum = shf_cum+shf*dt/dt0;
    
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