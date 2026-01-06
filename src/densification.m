function [dz, d] = densification(T, dz, d, re, dt, density_ice, T_mean, P_mean, albedo_desnity_threshold, albedo_method, sw_absorption_method, densification_method)

% densification computes the densification of snow/firn using the emperical model of
% Herron and Langway (1980) or the semi-emperical model of Anthern et al. (2010).
%
% Inputs:
%   densification_method = densification model to use:
%       1 = emperical model of Herron and Langway (1980)
%       2 = semi-imerical model of Anthern et al. (2010)
%       3 = physical model from Appendix B of Anthern et al. (2010)
%   d   = initial snow/firn density [kg m-3]
%   T   = temperature [K]
%   dz  = grid cell size [m]
%   P_mean   = average accumulation rate [kg m-2 yr-1]
%   dt  = time lapsed [s]
%   re  = effective grain radius [mm];
%   T_air  = mean annual temperature
%
%% FOR TESTING
% densification_method = 2;
% d = 800;
% T = 270;
% dz = 0.005;
% P_mean = 200;
% dt = 60*60;
% re = 0.7;
% T_mean = 273.15-18;
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
% Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and Thomas, E. R.:
% In situ measurements of Antarctic snow compaction compared with predictions of
% models, J. Geophys. Res., 115, F03011, https://doi.org/10.1029/2009JF001306, 2010. 
%
% Herron, M. and Langway, C.: Firn Densification: An Empirical Model, J. Glaciol.,
% 25, 373–385, https://doi.org/10.3189/S0022143000015239, 1980. 
%
% If you use GEMB, please cite the following:
%
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass
% Balance (GEMB): a model of firn processes for cryosphere research, Geosci.
% Model Dev., 16, 2277–2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023.

%% MAIN FUNCTION

d_tolerance  = 1e-11;

% specify constants
dt      = dt / 86400;  % convert from [s] to [d]
R       = 8.314;       % gas constant [mol-1 K-1]
CtoK    = 273.15;      % Kelvin to Celcius conversion/ice melt. point T in K
% Ec    = 60           % activation energy for self-diffusion of water
%                      % molecules through the ice tattice [kJ mol-1]
% Eg    = 42.4         % activation energy for grain growth [kJ mol-1]

% initial mass
mass_init = d .* dz;

% calculate new snow/firn density for:
%   snow with densities <= 550 [kg m-3]
%   snow with densities > 550 [kg m-3]
idx = d <= (550.0 + d_tolerance);
switch densification_method
    case 1 % Herron and Langway (1980)
        c0 = ( 11 * exp(-10160 ./ (T(idx)  * R))) .* P_mean/1000;
        c1 = (575 * exp(-21400 ./ (T(~idx) * R))) .* (P_mean/1000)^0.5;

    case 2 % Arthern et al. (2010) [semi-emperical]
        % common variable
        % NOTE: Ec=60000, Eg=42400 (i.e. should be in J not kJ)
        H = exp((-60000./(T * R)) + (42400./(T_mean .* R))) ...
            .* (P_mean * 9.81);

        c0 = 0.07 * H(idx);
        c1 = 0.03 * H(~idx);

    case 3 % Arthern et al. (2010) [physical model eqn. B1]
        % calcualte overburden pressure
        obp =[0;  (cumsum(dz(1:end-1)) .* d(1:end-1))];

        % common variable
        H = exp((-60000./(T * R))) .* obp ./ (re/1000).^2;
        c0 = 9.2E-9 * H(idx);
        c1 = 3.7E-9 * H(~idx);

    case 4 % Li and Zwally (2004)
        c0 = (P_mean./density_ice) .* max(139.21 - 0.542*T_mean,1) .* 8.36 .* max(CtoK - T,1.0).^-2.061;
        c1 = c0;
        c0 = c0(idx);
        c1 = c1(~idx);

    case 5 % Helsen et al. (2008)
        % common variable
        c0 = (P_mean./density_ice) .* max(76.138 - 0.28965 * T_mean, 1) .* 8.36 .* max(CtoK - T,1.0).^-2.061;
        c1 = c0;
        c0 = c0(idx);
        c1 = c1(~idx);

    case 6 % Ligtenberg and others (2011) [semi-emperical], Antarctica
        % common variable
        % From literature: H = exp((-60000.0/(T_mean * R)) + (42400.0/(T_mean * R))) * (P_mean * 9.81);
        H      = exp((-60000.0 ./ (T * R)) + (42400.0 ./ (T_mean .* R))) .* (P_mean .* 9.81);
        c0arth = 0.07 * H;
        c1arth = 0.03 * H;

        % ERA5 new albedo_method=1, sw_absorption_method=0
        if (albedo_method == 1) && (sw_absorption_method == 0)
            if abs(albedo_desnity_threshold - 820.0) < d_tolerance 
                % ERA5 v4 (Paolo et al., 2023)
                M0 = max(1.5131 - (0.1317 * log(P_mean)),0.25);
                M1 = max(1.8819 - (0.2158 * log(P_mean)),0.25);
            else
                % ERA5 new albedo_method=1, sw_absorption_method=0
                %M0 = max(1.8785 - (0.1811 * log(P_mean)),0.25);
                %M1 = max(2.0005 - (0.2346 * log(P_mean)),0.25);
                % ERA5 new albedo_method=1, sw_absorption_method=0, bare ice
                M0 = max(1.8422 - (0.1688 * log(P_mean)),0.25);
                M1 = max(2.4979 - (0.3225 * log(P_mean)),0.25);
            end
    		% ERA5 new albedo_method=2, sw_absorption_method=1
        elseif (albedo_method < 3) && (sw_absorption_method > 0)
            M0 = max(2.2191 - (0.2301 * log(P_mean)),0.25);
            M1 = max(2.2917 - (0.2710 * log(P_mean)),0.25);
        else
            % ERA5 new albedo_method=2, sw_absorption_method=0
            %elseif (albedo_method==2)
            %From Ligtenberg
            %H = exp((-60000.0/(T_mean * R)) + (42400.0/(T_mean * R))) * (P_mean * 9.81);
            %M0 = max(1.435 - (0.151 * log(P_mean)),0.25);
            %M1 = max(2.366 - (0.293 * log(P_mean)),0.25);
            %RACMO callibration, default (Gardner et al., 2023)
            M0 = max(1.6383 - (0.1691 * log(P_mean)),0.25);
            M1 = max(1.9991 - (0.2414 * log(P_mean)),0.25);
        end

        c0 = M0 * c0arth(idx);
        c1 = M1 * c1arth(~idx);

    case 7 % Kuipers Munneke and others (2015) [semi-emperical], Greenland
        % common variable
        % From literature: H = exp((-60000.0/(T[i] * R)) + (42400.0/(T[i] * R))) * (P_mean * 9.81);
        H = exp((-60000.0 ./ (T * R)) + (42400.0 ./ (T_mean .* R))) .* (P_mean .* 9.81);

        c0arth = 0.07 * H;
        c1arth = 0.03 * H;
        % ERA5 new albedo_method=1, sw_absorption_method=0
        if (albedo_method == 1) && (sw_absorption_method == 0)
            if abs(albedo_desnity_threshold - 820.0) < d_tolerance 
                % ERA5 v4
                M0 = max(1.3566 - (0.1350 * log(P_mean)),0.25);
                M1 = max(1.8705 - (0.2290 * log(P_mean)),0.25);
            else
                % ERA5 new albedo_method=1, sw_absorption_method=0
                %M0 = max(1.4574 - (0.1123 * log(P_mean)),0.25);
                %M1 = max(2.0238 - (0.2070 * log(P_mean)),0.25);
                % ERA5 new albedo_method=1, sw_absorption_method=0, bare ice
                M0 = max(1.4318 - (0.1055 * log(P_mean)),0.25);
                M1 = max(2.0453 - (0.2137 * log(P_mean)),0.25);
            end
    		% ERA5 new albedo_method=2, sw_absorption_method=1
        elseif albedo_method<3 && sw_absorption_method>0
            M0 = max(1.7834 - (0.1409 * log(P_mean)),0.25);
            M1 = max(1.9260 - (0.1527 * log(P_mean)),0.25);
        else
            % ERA5 new albedo_method=2, sw_absorption_method=0
            %elseif (albedo_method==2)
            % From Kuipers Munneke
            %M0 = max(1.042 - (0.0916 * log(P_mean)),0.25);
            %M1 = max(1.734 - (0.2039 * log(P_mean)),0.25);
            %RACMO callibration, default (Gardner et al., 2023)
            M0 = max(1.2691 - (0.1184 * log(P_mean)),0.25);
            M1 = max(1.9983 - (0.2511 * log(P_mean)),0.25);
        end

        c0 = M0*c0arth(idx);
        c1 = M1*c1arth(~idx);

end

% new snow density
d(idx)  = d(idx)  + (c0 .* (density_ice - d(idx)) / 365 * dt);
d(~idx) = d(~idx) + (c1 .* (density_ice - d(~idx)) / 365 * dt);

%disp((num2str(nanmean(c0 .* (density_ice - d(idx)) / 365 * dt))))

% do not allow densities to exceed the density of ice
d(d > (density_ice - d_tolerance)) = density_ice;

% calculate new grid cell length
dz = mass_init ./ d;

end