function swf = calculate_shortwave_radiation(dz, d, re, a_surface, a_diffuse_surface, ...
    ClimateForcingStep, ModelParam)
% calculate_shortwave_radiation distributes absorbed shortwave radiation vertically within snow/ice.
%
%% Syntax 
% 
% swf = calculate_shortwave_radiation(dz, d, re, a_surface, a_diffuse_surface, ...
%     ClimateForcingStep, ModelParam)
%
%% Description
% 
% This function calculates the vertical distribution of absorbed shortwave (solar) 
% radiation within the snow and firn column . 
% Depending on the selected model configuration, the radiation is either:
%
% 1. Surface Absorption: All net shortwave energy is absorbed entirely by the 
%    top grid cell (sw_absorption_method = 0).
% 2. Subsurface Penetration: Shortwave energy penetrates the surface and is 
%    absorbed by deeper layers (sw_absorption_method = 1).
%
% When subsurface penetration is enabled, the extinction of radiation is modeled 
% using one of two methods:
% * Density-Dependent Extinction: The extinction coefficient scales linearly 
%   between snow and ice values based on density (Bassford, 2002 formulation).
%   % * Spectral-Dependent Extinction (BruneLeFebre): Radiation is split into three 
%   spectral bands (UV/visible, near-IR, IR) with band-specific albedos and 
%   grain-size dependent extinction coefficients (Lefebre et al., 2003).
%
%% Inputs
% 
%   dz                   : m            Grid cell thickness (vector).
%   d                    : kg m^-3      Grid cell density (vector).
%   re                   : mm           Grid cell effective grain radius (vector).
%   a_surface            : fraction     Broadband surface albedo.
%   a_diffuse_surface    : fraction     Diffuse surface albedo.
%   ClimateForcingStep   : struct       Current time-step forcing data:
%     .dsw               : W m^-2       Downward shortwave radiative flux.
%     .dsw_diffuse       : W m^-2       Downward diffuse shortwave flux.
%   ModelParam           : struct       Model parameters:
%     .sw_absorption_method : integer   0 (surface only) or 1 (subsurface penetration).
%     .albedo_method        : string    Albedo scheme selection (e.g., "GardnerSharp", "BruneLeFebre").
%     .density_ice          : kg m^-3   Density of ice.
% 
%% Outputs
% 
%   swf                  : W m^-2       Vector of absorbed shortwave radiation flux for each grid cell.
%
%% Example 
% 
% % Example call assuming initialized variables:
% sw_flux = calculate_shortwave_radiation(dz, density, grain_radius, 0.85, 0.85, forcing, params);
% 
%% References 
% This function uses formulations from the following references: 
% 
% Lefebre, F., Gallée, H., van Ypersele, J.-P., and Greuell, W.: Modeling of 
% snow and ice melt at ETH Camp (West Greenland): A study of surface albedo, 
% J. Geophys. Res., 108, 4231, https://doi.org/10.1029/2001JD001160, 2003. 
% 
% Greuell, W. and Konzelmann, T.: Numerical modelling of the energy balance 
% and the englacial temperature of the Greenland Ice Sheet, Calculations for 
% the ETH-Camp location (West Greenland, 1155 m a.s.l.), Global Planet. 
% Change, 9, 91–114, 1994. 
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

%% SHORTWAVE FUNCTION
d_tolerance  = 1e-11;

% Initialize variables:
m = length(d);
swf = zeros(m,1);

if (ModelParam.sw_absorption_method == 0) || ...
        ((ModelParam.density_ice - d(1))<d_tolerance)  % all sw radation is absorbed by the top grid cell

    % calculate surface shortwave radiation fluxes [W m-2]
    if (ModelParam.albedo_method == "GardnerSharp") % ModelParam.albedo_method = "gardner_2009"
        swf(1) = (1.0 - a_surface) * max(0.0, (ClimateForcingStep.dsw - ClimateForcingStep.dsw_diffuse)) ...
            +  (1.0 - a_diffuse_surface) * ClimateForcingStep.dsw_diffuse;
    else
        swf(1) = (1 - a_surface) * ClimateForcingStep.dsw;
    end
    
else % sw radation is absorbed at depth within the glacier
    
    if ModelParam.albedo_method == "BruneLeFebre"    % ModelParam.albedo_method = "brun_1992" function of effective radius (3 spectral bands)
        
        % convert effective radius [mm] to grain size [m]
        gsz = (re * 2) / 1000;
        
        % Spectral fractions [0.3-0.8um 0.8-1.5um 1.5-2.8um]
        % (Lefebre et al., 2003)
        sF = [0.606; 0.301; 0.093];
        
        % initialize variables
        B1_cum = ones(m+1,1);
        B2_cum = B1_cum;
        
        % spectral albedos:
        % 0.3 - 0.8um
        a1 = min(0.98, 0.95 - 1.58 * gsz(1)^0.5);

        % 0.8 - 1.5um
        a2 = max(0, 0.95 - 15.4 * gsz(1)^0.5);
        
        % 1.5 - 2.8um
        a3 = max(0.127, 0.88 + 346.3 * gsz(1) - 32.31 * gsz(1)^0.5);
        
        % seperate net shortwave radiative flux into spectral ranges
        swfS = (sF * ClimateForcingStep.dsw) .* (1 - [a1; a2; a3]);
        
        % absorption coefficient for spectral range:
        h = d ./(gsz.^0.5);
        B1 = .0192 * h;                 % 0.3 - 0.8um
        B2 = .1098 * h;                 % 0.8 - 1.5um
        % B3 = +inf                     % 1.5 - 2.8um
        
        % cumulative extinction factors
        B1_cum(2:end) =  cumprod(exp(-B1.*dz));
        B2_cum(2:end) =  cumprod(exp(-B2.*dz));
        
        % flux across grid cell boundaries
        Qs1 = swfS(1) * B1_cum;
        Qs2 = swfS(2) * B2_cum;
        
        % net energy flux to each grid cell
        swf = (Qs1(1:m)-Qs1(2:m+1)) + (Qs2(1:m)-Qs2(2:m+1));
        
        % add flux absorbed at surface
        swf(1) = swf(1) + swfS(3);
        
    else % function of grid cell density
        
        % fraction of sw radiation absorbed in top grid cell
        % (wavelength > 0.8um)
        SWs = 0.36;
        
        % SWs and SWss coefficients need to be better constranted. Greuell
        % and Konzelmann 1994 used SWs = 0.36 and SWss = 0.64 a_surface this the
        % the % of SW radiation with wavelengths > and < 800 nm
        % respectively.  This, however, may not account for the fact that
        % the albedo of wavelengths > 800 nm has a much lower albedo.
        
        % calculate surface shortwave radiation fluxes [W m-2]
        swf_s = SWs * (1 - a_surface) * ClimateForcingStep.dsw;
        
        % calculate surface shortwave radiation fluxes [W m-2]
        swf_ss = (1-SWs) * (1 - a_surface) * ClimateForcingStep.dsw;
        
        % SW allowed to penetrate into snowpack
        Bs = 10;    % snow SW extinction coefficient [m-1] (Bassford,2006)
        Bi = 1.3;   % ice SW extinction coefficient [m-1] (Bassford,2006)
        
        % calculate extinction coefficient B [m-1] vector
        B = Bs + (300 - d) .* ((Bs - Bi)/(ModelParam.density_ice - 300));
        
        % cumulative extinction factor
        B_cum =  [1; cumprod(exp(-B.*dz))];
        
        % flux across grid cell boundaries
        Qs = swf_ss * B_cum;
        
        % net energy flux to each grid cell
        swf = (Qs(1:m)-Qs(2:m+1));
        
        % add flux absorbed at surface
        swf(1) = swf(1) + swf_s;
    end
end
