function [shf, lhf, L] = turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, ClimateForcingStep)

% CALC_TURBULENT_HEAT_FLUX Calculates sensible and latent heat fluxes.
%
%   [shf, lhf, coefM, coefHT, coefHQ] = calc_turbulent_heat_flux(...)
%
%   Inputs:
%       ClimateForcingStep.T_air       : Air Temperature [K] (at height ClimateForcingStep.Tz)
%       T_surface   : Surface Temperature [K]
%       ClimateForcingStep.p_air       : Atmospheric Pressure [Pa]
%       ClimateForcingStep.e_air       : Screen level vapor pressure [Pa]
%       ClimateForcingStep.V           : Wind Speed [m s-1] (at height ClimateForcingStep.Vz)
%       density_air : Air Density [kg m-3]
%       ClimateForcingStep.Vz          : Measurement height for wind [m]
%       ClimateForcingStep.Tz          : Measurement height for temperature [m]
%       z0          : Aerodynamic roughness length [m]
%       zT          : Thermal roughness length [m]
%       zQ          : Humidity roughness length [m]
%
%   Outputs:
%       shf    : Sensible Heat Flux [W m-2]
%       lhf    : Latent Heat Flux [W m-2]
%       coefM  : Drag coefficient for momentum
%       coefHT : Transfer coefficient for heat
%       coefHQ : Transfer coefficient for humidity

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
        (2.0 * g * (ClimateForcingStep.T_air - T_surface)) / (ClimateForcingStep.Tz*(ClimateForcingStep.T_air + T_surface)*(((ClimateForcingStep.V/ClimateForcingStep.Vz)^2.0)));
    
    % Constants for Beljaars and Holtslag (1991)
    a1 = 1.0; b1 = 2.0 / 3.0; c1 = 5.0; d1 = 0.35;
    PhiMz0 = 0.0; PhiHzT = 0.0; PhiHzQ = 0.0;
    
    if (Ri > 0.0 + T_tolerance ) % --- STABLE ---
        if (Ri < 0.2 - T_tolerance )
            zL = Ri/(1.0 - 5.0*Ri);
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