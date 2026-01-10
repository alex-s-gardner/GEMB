function [shf, lhf, L] = turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, ClimateForcingStep)
% turbulent_heat_flux computes sensible and latent heat fluxes using 
% Monin-Obukhov similarity theory.
%
%% Syntax
%
% [shf, lhf, L] = turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, ClimateForcingStep)
%
%% Description
%
% This function calculates the turbulent exchange of energy (sensible heat) 
% and mass (moisture/latent heat) between the snow/ice surface and the 
% atmosphere . It employs the bulk aerodynamic method corrected for 
% atmospheric stability using Monin-Obukhov similarity theory.
%
% The calculation proceeds in three stages:
% 1. Stability Analysis: Calculates the Bulk Richardson Number (Ri) to 
%    characterize the stability of the surface layer based on the vertical 
%    temperature gradient and wind shear.
% 2. Stability Corrections: Adjusts the bulk transfer coefficients using 
%    integrated stability functions (Psi) derived from Beljaars and Holtslag 
%    (1991) for stable conditions and standard logarithmic profiles for 
%    unstable conditions.
% 3. Flux Computation:
%    * Sensible Heat Flux (SHF): Driven by the temperature difference between 
%      the surface and the air.
%    * Latent Heat Flux (LHF): Driven by the vapor pressure gradient. The 
%      function automatically detects the surface phase (ice vs. water) to 
%      apply the appropriate latent heat constant (sublimation vs. vaporization) 
%      and saturation vapor pressure formula.
%
%% Inputs
%
%  T_surface          : K            Surface temperature.
%  density_air        : kg m^-3      Air density.
%  z0                 : m            Aerodynamic roughness length for momentum.
%  zT                 : m            Roughness length for heat.
%  zQ                 : m            Roughness length for moisture.
%  ClimateForcingStep : struct       Forcing data for the current time step:
%    .V               : m s^-1       Wind speed.
%    .p_air           : Pa           Air pressure.
%    .T_air           : K            Air temperature.
%    .e_air           : Pa           Air vapor pressure.
%    .Tz              : m            Measurement height for temperature.
%    .Vz              : m            Measurement height for wind.
%
%% Outputs
%
%  shf                : W m^-2       Sensible heat flux (positive toward surface).
%  lhf                : W m^-2       Latent heat flux (positive toward surface).
%  L                  : J kg^-1      Latent heat of vaporization or sublimation used.
%
%% References
%
% Physics implementations based on:
% Beljaars, A. C. M., & Holtslag, A. A. M. (1991). Flux parameterization over 
%   land surfaces for atmospheric models. Journal of Applied Meteorology.
% Ohmura, A. (1982). Climate and Energy-Balance on the Arctic Tundra.
% Murray, F. W. (1967). (Saturation vapor pressure over water).
% Bolton, D. (1980). (Saturation vapor pressure over ice).
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

%% CONSTANTS & INITIALIZATION
T_tolerance  = 1e-10;       % Tolerance
CA           = 1005.0;      % Specific heat capacity of air [J kg-1 K-1]
g            = 9.81;        % Gravity [m s-2]

% Constants for Latent Heat
CtoK = 273.15;      % Kelvin to Celsius conversion
LV   = 2.495E6;     % Latent heat of vaporization [J kg-1]
LS   = 2.8295e6;    % Latent heat of sublimation [J kg-1]

% Bulk-transfer coefficient (Neutral)
An = 0.4^2;         
C  = An * ClimateForcingStep.V;        

%% STABILITY CORRECTION (Monin-Obukhov)
% Ohmura, A., 1982: Climate and Energy-Balance on the Arctic Tundra.

% Bulk Richardson Number (Ri)
Ri = ((100000 / ClimateForcingStep.p_air)^0.286) * ...
    (2.0 * g * (ClimateForcingStep.T_air - T_surface)) / ...
    (ClimateForcingStep.Tz * (ClimateForcingStep.T_air + T_surface) * ...
    (((ClimateForcingStep.V/ClimateForcingStep.Vz)^2.0)));

% Constants for Beljaars and Holtslag (1991)
a1 = 1.0; b1 = 2.0 / 3.0; c1 = 5.0; d1 = 0.35;
PhiMz0 = 0.0; PhiHzT = 0.0; PhiHzQ = 0.0;

if (Ri > 0.0 + T_tolerance ) % --- STABLE ---
    if (Ri < 0.2 - T_tolerance )
        zL = Ri / (1.0 - (5.0 * Ri));
    else
        zL = Ri;
    end
    
    zLM = max(zL / ClimateForcingStep.Vz * z0, 1e-3);
    zLT = max(zL / ClimateForcingStep.Tz * zT, 1e-3);
    
    % Integrated Stability Functions (Psi)
    PhiMz  = -1*(a1*zL + b1*(zL-c1/d1)*exp(-1*d1*zL) + b1*c1/d1);
    PhiHz  = -1*((1+2*a1*zL/3)^1.5 + b1*(zL-c1/d1)*exp(-1*d1*zL) + b1*c1/d1 - 1.0);
    
    PhiMz0 = -1*(a1*zLM + b1*(zLM-c1/d1)*exp(-1*d1*zLM) + b1*c1/d1);
    PhiHzT = -1*((1+2*a1*zLT/3)^1.5 + b1*(zLT-c1/d1)*exp(-1*d1*zLT) + b1*c1/d1 - 1.0);
    
    PhiHzQ = PhiHzT;

else % --- UNSTABLE ---
    zL  = Ri/1.5; 
    xm = (1.0 - 19.0*zL)^-0.25;
    PhiMz = 2.0*log((1+xm)/2.0) + log((1+xm^2)/2.0) - 2*atan(xm) + pi/2;
    
    xh = 0.95*(1.0 - 11.6*zL)^(-0.5);
    PhiHz = 2.0*log((1.0+xh^2)/2.0);
end

% Final Transfer Coefficients
coefM  = log(ClimateForcingStep.Vz/z0) - PhiMz + PhiMz0; 
coefHT = log(ClimateForcingStep.Tz/zT) - PhiHz + PhiHzT; 
coefHQ = log(ClimateForcingStep.Tz/zQ) - PhiHz + PhiHzQ; 

%% SENSIBLE HEAT FLUX [W m-2]
shf = density_air * C * CA * (ClimateForcingStep.T_air - T_surface) * (100000/ClimateForcingStep.p_air)^0.286;
shf = shf/(coefM*coefHT);

%% LATENT HEAT FLUX [W m-2]

% Determine Phase (Melting or Freezing) for Latent Heat Constant and Saturation Pressure
if (T_surface >= CtoK - T_tolerance )
    % Liquid water surface
    L = LV; 
    % Saturation Vapor Pressure (Murray 1967)
    eS = 610.78 * exp(17.2693882 * (T_surface - CtoK - 0.01) / (T_surface - 35.86));
else
    % Ice surface
    L = LS; 
    % Saturation Vapor Pressure (Ding et al., 2019 / Bolton 1980)
    eS = 610.78 * exp(21.8745584 * (T_surface - CtoK - 0.01) / (T_surface - 7.66));
end

% Calculate Latent Heat Flux
% 461.9 is the specific gas constant for water vapor [J kg-1 K-1]
lhf = C * L * (ClimateForcingStep.e_air - eS) / (461.9 * (ClimateForcingStep.T_air + T_surface) / 2.0);

% Adjust using Monin-Obukhov stability theory
lhf = lhf / (coefM * coefHQ);
end