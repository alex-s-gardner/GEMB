function [shf, lhf, EC] = turbulentFlux(Ta, Ts, V, eAir, pAir, ds, Ws, Vz, Tz)
% turbulentFlux computes the surface sensible and latent heat fluxes and
% calculates mass loss and accretion due to condensation/evaporation.
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
% 
% 
%% Outputs
% 
% 
%% Documentation
% 
% For complete documentation, see: https://github.com/alex-s-gardner/GEMB 
% 
%% References 
% This function uses formulations from the following references: 
% 
% Murphy, D. M. and Koop, T.: Review of the Vapour Pressures of Ice and 
% Supercooled Water for Atmospheric Applications, Q. J. Roy. Meteor. Soc., 
% 131, 1539–1565, https://doi.org/10.1256/qj.04.94, 2005. 
% 
% Ohmura, A., 1982: Climate and Energy-Balance on the Arctic Tundra.
% Journal of Climatology, 2, 65-84.
% 
% If you use GEMB, please cite the following: 
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass 
% Balance (GEMB): a model of firn processes for cryosphere research, Geosci. 
% Model Dev., 16, 2277–2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023.

% Reference: 
% Dingman, 2002.

%% INPUTS:
%   Ta: 2 m air temperature [K]
%   Ts: snow/firn/ice surface temperature [K]
%   V: wind speed [m s^-^1]
%   eAir: screen level vapor pressure [Pa]
%   pAir: surface pressure [Pa]
%   ds: surface density [kg/m^3]
%   Ws: surface liquid water content [kg/m^2]
%   Vz: height above ground at which wind (V) eas sampled [m]
%   Tz: height above ground at which temperature (T) was sampled [m]

%% FUNCTION INITILIZATION 

% CA = 1005;                    % heat capacity of air (J kg-1 k-1)
% LF = 0.3345E6;                % latent heat of fusion(J kg-1)
% LV = 2.495E6;                 % latent heat of vaporization(J kg-1)
% dIce = 910;                   % density of ice [kg m-3]

% calculated air density [kg/m3]
d_air = 0.029 * pAir /(8.314 * Ta);

%% Determine Surface Roughness
% Bougamont, 2006

% wind/temperature surface roughness height [m]
if ds < 910 && Ws == 0
    z0 = 0.00012;               % 0.12 mm for dry snow
elseif ds >= 910
    z0 = 0.0032;                % 3.2 mm for ice 
else
    z0 = 0.0013;                % 1.3 mm for wet snow
end

%% Monin Obukhov Stability Correction
% Reference:
% Ohmura, A., 1982: Climate and Energy-Balance on the Arctic Tundra.
% Journal of Climatology, 2, 65-84.

% if V = 0 goes to infinity therfore if V = 0 change
V(V < 0.01) = 0.01;

% calculate the Bulk Richardson Number (Ri)
Ri = (2*9.81* (Vz - z0) * (Ta - Ts)) / ((Ta + Ts)* V^2);
    
% calculate MoninObukhov stability factors 'coef_M' and 'coef_H'

% do not allow Ri to exceed 0.19
Ri(Ri > 0.19) = 0.19;

% calculate momentum 'coef_M' stability factor
if Ri > 0
    % if stable
    coef_M = (1-5.2*Ri)^-1;
else
    coef_M = (1-18*Ri)^-0.25;
end

% calculate heat/wind 'coef_H' stability factor
if Ri < -0.03
    coef_H = 1.3 * coef_M;
else
    coef_H = coef_M;
end
%% Bulk-transfer coefficient

An =  0.4^2 / log(Tz/z0)^2;     % Bulk-transfer coefficient
C = An * d_air * V;             % shf & lhf common coefficient

%% Sensible Heat
% calculate the sensible heat flux [W m-2](Patterson, 1998)

shf = C * 1005 * (Ta - Ts);

% adjust using Monin-Obukhov stability theory
shf = shf / (coef_M * coef_H);

%% Latent Heat
% determine if snow pack is melting & calcualte surface vapour pressure
% over ice or liquid water

if Ts >= 273.15
    L = 2.495E6;
    
    % for an ice surface Murphy and Koop, 2005 [Equation 7]
    eS = exp(9.550426 - 5723.265/Ts + 3.53068 * log(Ts)...
        - 0.00728332 * Ts);
else
    L = 2.8295E6; % latent heat of sublimation
    % for liquid surface (assume liquid on surface when Ts == 0 deg C)
    % Wright (1997), US Meteorological Handbook from Murphy and Koop,
    % 2005 Apendix A
    eS = 611.21 * exp(17.502 * (Ts - 273.15) / (240.97 + Ts - 273.15));
end

% Latent heat flux [W m-2]
lhf = C * L * (eAir - eS) * 0.622 / pAir;

% adjust using MoninObukhov stability theory (if lhf '+' then there is
% energy and mass gained at the surface, if '-' then there is mass and 
% energy loss at the surface. 
lhf = lhf / (coef_M * coef_H);

% mass loss (-)/acreation(+) due to evaporation/condensation [kg]
EC = lhf * 86400 / L;

end