function [dz, d] = densification(T, dz, d, re, ClimateForcingStep, ModelParam)
% densification computes the densification of snow/firn.
%
%% Syntax
%
%  [dz, d] = densification(T, dz, d, re, ClimateForcingStep, ModelParam)
%
%% Description
%
% [dz, d] = densification(T, dz, d, re, ClimateForcingStep, ModelParam)
%
%% Inputs:
%   ModelParam.densification_method = densification model to use:
%       1-"HerronLangway": emperical model of Herron and Langway (1980)
%       2-"Anthern": semi-emperical model of Anthern et al. (2010)
%       3-"AnthernB": DO NOT USE: physical model from Appendix B of Anthern et al. (2010)
%       4-"Li and Zwally": DO NOT USE: emperical model of Li and Zwally (2004)
%       5-"Helsen": DO NOT USE: modified emperical model (4) by Helsen et al. (2008)
%       6-"Ligtenberg": Antarctica semi-emperical model of Ligtenberg et al. (2011)
%
%   d                           = initial snow/firn density [kg m-3]
%   T                           = temperature [K]
%   dz                          = grid cell size [m]
%   re                          = effective grain radius [mm];
%   T_air                       = mean annual temperature
%   ClimateForcingStep.P_mean   = average accumulation rate [kg m-2 yr-1]
%   ClimateForcingStep.dt       = time lapsed [s]
%
%% FOR TESTING
% ModelParam.densification_method = "Anthern";
% d                               = 800;
% T                               = 270;
% dz                              = 0.005;
% ClimateForcingStep.P_mean       = 200;
% ClimateForcingStep.dt           = 60*60;
% re                              = 0.7;
% ClimateForcingStep.T_air_mean   = 273.15-18;
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

%% MAIN FUNCTION

d_tolerance  = 1e-11;

% specify constants
dt      = ClimateForcingStep.dt / 86400;  % convert from [s] to [d]
R       = 8.314;                          % gas constant [mol-1 K-1]
CtoK    = 273.15;                         % Kelvin to Celcius conversion/ice melt. point T in K
% Ec    = 60                              % activation energy for self-diffusion of water
%                                         % molecules through the ice tattice [kJ mol-1]
% Eg    = 42.4                            % activation energy for grain growth [kJ mol-1]

% initial mass
mass_init = d .* dz;

% calculate new snow/firn density for:
%   snow with densities <= 550 [kg m-3]
%   snow with densities > 550 [kg m-3]

idx = d <= (550.0 + d_tolerance);
switch ModelParam.densification_method
    case "HerronLangway" % Herron and Langway (1980)
        c0 = ( 11 * exp(-10160 ./ (T(idx)  * R))) .* ClimateForcingStep.P_mean/1000;
        c1 = (575 * exp(-21400 ./ (T(~idx) * R))) .* (ClimateForcingStep.P_mean/1000)^0.5;

    case "Anthern" % Arthern et al. (2010) [semi-emperical]
        % common variable
        % NOTE: Ec=60000, Eg=42400 (i.e. should be in J not kJ)
        H = exp((-60000./(T * R)) + (42400./(ClimateForcingStep.T_air_mean.* R))) ...
            .* (ClimateForcingStep.P_mean * 9.81);

        c0 = 0.07 * H(idx);
        c1 = 0.03 * H(~idx);

    case "AnthernB" % Arthern et al. (2010) [physical model eqn. B1]
        % calcualte overburden pressure
        obp =[0;  (cumsum(dz(1:end-1)) .* d(1:end-1))];

        % common variable
        H = exp((-60000./(T * R))) .* obp ./ (re/1000).^2;
        c0 = 9.2E-9 * H(idx);
        c1 = 3.7E-9 * H(~idx);

    case "LiZwally" % Li and Zwally (2004)
        c0 = (ClimateForcingStep.P_mean./ModelParam.density_ice) .* ...
            max(139.21 - 0.542*ClimateForcingStep.T_air_mean,1) .* 8.36 .* max(CtoK - T,1.0).^-2.061;
        c1 = c0;
        c0 = c0(idx);
        c1 = c1(~idx);

    case "Helsen" % Helsen et al. (2008)
        % common variable
        c0 = (ClimateForcingStep.P_mean./ModelParam.density_ice) .* ...
            max(76.138 - 0.28965 * ClimateForcingStep.T_air_mean, 1) .* 8.36 .* max(CtoK - T,1.0).^-2.061;
        c1 = c0;
        c0 = c0(idx);
        c1 = c1(~idx);

    case "Ligtenberg" % Ligtenberg and others (2011) [semi-emperical], Antarctica
        H = exp((-60000.0 ./ (T * R)) + (42400.0 ./ (ClimateForcingStep.T_air_mean.* R))) .* ...
            (ClimateForcingStep.P_mean .* 9.81);
      
        c0arth = 0.07 * H;
        c1arth = 0.03 * H;

        M01 = densification_lookup_M01(ModelParam.densification_coeffs_M01);

        if numel(M01) == 4
            M0 = max(M01(1) - (M01(2) * log(ClimateForcingStep.P_mean)),0.25);
            M1 = max(M01(3) - (M01(4) * log(ClimateForcingStep.P_mean)),0.25);
        else
            if abs(ModelParam.density_ice - 820.0) < d_tolerance
                M0 = max(M01(1,1) - (M01(1,2) * log(ClimateForcingStep.P_mean)),0.25);
                M1 = max(M01(1,3) - (M01(1,4) * log(ClimateForcingStep.P_mean)),0.25);
            else
                M0 = max(M01(2,1) - (M01(2,2) * log(ClimateForcingStep.P_mean)),0.25);
                M1 = max(M01(2,3) - (M01(2,4) * log(ClimateForcingStep.P_mean)),0.25);
            end
        end

        c0 = M0 * c0arth(idx);
        c1 = M1 * c1arth(~idx);
    otherwise
        error("unrecognized densification method")

end

% new snow density
d(idx)  = d(idx)  + (c0 .* (ModelParam.density_ice - d(idx)) / 365 * dt);
d(~idx) = d(~idx) + (c1 .* (ModelParam.density_ice - d(~idx)) / 365 * dt);

%disp((num2str(nanmean(c0 .* (ModelParam.density_ice - d(idx)) / 365 * dt))))

% do not allow densities to exceed the density of ice
d(d > (ModelParam.density_ice - d_tolerance)) = ModelParam.density_ice;

% calculate new grid cell length
dz = mass_init ./ d;

end