function [dz, density] = calculate_density(temperature, dz, density, grain_radius, ClimateForcingStep, ModelParam)
% calculate_density computes the densification of snow/firn.
%
%% Syntax
%
%  [dz, density] = calculate_density(temperature, dz, density, grain_radius, ClimateForcingStep, ModelParam)
%
%% Description
%
% [dz, density] = calculate_density(temperature, dz, density, grain_radius, ClimateForcingStep, ModelParam)
%
%% Inputs:
%   ModelParam.densification_method = densification model to use:
%      "Arthern": (default) semi-empirical model of Arthern et al. (2010)
%      "HerronLangway": empirical model of Herron and Langway (1980)
%      "Ligtenberg": Antarctica semi-empirical model of Ligtenberg et al. (2011)
%      "ArthernB": DO NOT USE: physical model from Appendix B of Arthern et al. (2010)
%      "Li and Zwally": DO NOT USE: empirical model of Li and Zwally (2004)
%      "Helsen": DO NOT USE: modified empirical model (4) by Helsen et al. (2008)
%
%   density                               = initial snow/firn density [kg m-3]
%   temperature                           = temperature [K]
%   dz                                    = grid cell height [m]
%   grain_radius                          = effective grain radius [mm];
%   temperature_air                       = mean annual temperature
%   ClimateForcingStep.precipitation_mean = average accumulation rate [kg m-2 yr-1]
%   ClimateForcingStep.dt                 = time lapsed [s]
%
%% FOR TESTING
% ModelParam.densification_method         = "Arthern";
% density                                 = 800;
% temperature                             = 270;
% dz                                      = 0.005;
% ClimateForcingStep.precipitation_mean   = 200;
% ClimateForcingStep.dt                   = 60*60;
% grain_radius                            = 0.7;
% ClimateForcingStep.temperature_air_mean = 273.15-18;
%
%% References 
% 
% Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and Thomas, E. R.: In situ measurements of  
% Antarctic snow compaction compared with predictions of models, J. Geophys. Res., 115, F03011,
% https://doi.org/10.1029/2009JF001306, 2010. 
%
% Helsen, M. M., van den Broeke, M. R., van de Wal, R. S. W., van de Berg, W. J., van Meijgaard, E., Davis, 
% C. H., Li, Y., and Goodwin, I.: Elevation changes in Antarctica mainly determined by accumulation 
% variability, Science, 320, 1626–1629, https://doi.org/10.1126/science.1153894, 2008. 
%
% Herron, M. and Langway, C.: Firn Densification: An Empirical Model, J. Glaciol., 25, 373–385, 
% https://doi.org/10.3189/S0022143000015239, 1980. 
%
% Li, J. and Zwally, H.: Modeling the density variation in the shallow firn layer, Ann. Glaciol., 38, 
% 309–313, https://doi.org/10.3189/172756404781814988, 2004. 
%
% Ligtenberg, S. R. M., Helsen, M. M., and van den Broeke, M. R.: An improved semi-empirical model for 
% the densification of Antarctic firn, The Cryosphere, 5, 809–819, https://doi.org/10.5194/tc-5-809-2011, 
% 2011. 
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
CtoK    = 273.15;                         % kelvin to celsius conversion/ice melt. point temperature in K
% Ec    = 60                              % activation energy for self-diffusion of water
%                                         % molecules through the ice tattice [kJ mol-1]
% Eg    = 42.4                            % activation energy for grain growth [kJ mol-1]

% initial mass
mass_init = density .* dz;

% calculate new snow/firn density for:
%   snow with densities <= 550 [kg m-3]
%   snow with densities > 550 [kg m-3]

idx = density <= (550.0 + d_tolerance);
switch ModelParam.densification_method
    case "HerronLangway" % Herron and Langway (1980)
        c0 = ( 11 * exp(-10160 ./ (temperature(idx)  * R))) .* ClimateForcingStep.precipitation_mean/1000;
        c1 = (575 * exp(-21400 ./ (temperature(~idx) * R))) .* (ClimateForcingStep.precipitation_mean/1000)^0.5;

    case "Arthern" % Arthern et al. (2010) [semi-empirical]
        % common variable
        % NOTE: Ec=60000, Eg=42400 (i.e. should be in J not kJ)
        H = exp((-60000./(temperature * R)) + (42400./(ClimateForcingStep.temperature_air_mean.* R))) ...
            .* (ClimateForcingStep.precipitation_mean * 9.81);

        c0 = 0.07 * H(idx);
        c1 = 0.03 * H(~idx);

    case "ArthernB" % Arthern et al. (2010) [physical model eqn. B1]
        % calcualte overburden pressure
        obp =[0;  (cumsum(dz(1:end-1)) .* density(1:end-1))];

        % common variable
        H = exp((-60000./(temperature * R))) .* obp ./ (grain_radius/1000).^2;
        c0 = 9.2E-9 * H(idx);
        c1 = 3.7E-9 * H(~idx);

    case "LiZwally" % Li and Zwally (2004)
        c0 = (ClimateForcingStep.precipitation_mean./ModelParam.density_ice) .* ...
            max(139.21 - 0.542*ClimateForcingStep.temperature_air_mean,1) .* 8.36 .* max(CtoK - temperature,1.0).^-2.061;
        c1 = c0;
        c0 = c0(idx);
        c1 = c1(~idx);

    case "Helsen" % Helsen et al. (2008)
        % common variable
        c0 = (ClimateForcingStep.precipitation_mean./ModelParam.density_ice) .* ...
            max(76.138 - 0.28965 * ClimateForcingStep.temperature_air_mean, 1) .* 8.36 .* max(CtoK - temperature,1.0).^-2.061;
        c1 = c0;
        c0 = c0(idx);
        c1 = c1(~idx);

    case "Ligtenberg" % Ligtenberg and others (2011) [semi-empirical], Antarctica
        H = exp((-60000.0 ./ (temperature * R)) + (42400.0 ./ (ClimateForcingStep.temperature_air_mean.* R))) .* ...
            (ClimateForcingStep.precipitation_mean .* 9.81);
      
        c0arth = 0.07 * H;
        c1arth = 0.03 * H;

        M01 = densification_lookup_M01(ModelParam.densification_coeffs_M01);

        if numel(M01) == 4
            M0 = max(M01(1) - (M01(2) * log(ClimateForcingStep.precipitation_mean)),0.25);
            M1 = max(M01(3) - (M01(4) * log(ClimateForcingStep.precipitation_mean)),0.25);
        else
            if abs(ModelParam.density_ice - 820.0) < d_tolerance
                M0 = max(M01(1,1) - (M01(1,2) * log(ClimateForcingStep.precipitation_mean)),0.25);
                M1 = max(M01(1,3) - (M01(1,4) * log(ClimateForcingStep.precipitation_mean)),0.25);
            else
                M0 = max(M01(2,1) - (M01(2,2) * log(ClimateForcingStep.precipitation_mean)),0.25);
                M1 = max(M01(2,3) - (M01(2,4) * log(ClimateForcingStep.precipitation_mean)),0.25);
            end
        end

        c0 = M0 * c0arth(idx);
        c1 = M1 * c1arth(~idx);
    otherwise
        error("unrecognized densification method")

end

% new snow density
density(idx)  = density(idx)  + (c0 .* (ModelParam.density_ice - density(idx)) / 365 * dt);
density(~idx) = density(~idx) + (c1 .* (ModelParam.density_ice - density(~idx)) / 365 * dt);

%disp((num2str(nanmean(c0 .* (ModelParam.density_ice - density(idx)) / 365 * dt))))

% do not allow densities to exceed the density of ice
density(density > (ModelParam.density_ice - d_tolerance)) = ModelParam.density_ice;

% calculate new grid cell length
dz = mass_init ./ density;

end