function [d, dz] = densification(denIdx, d, T, dz, C, dt, re, Tmean)

%% THIS NEEDS TO BE DOUBLE CHECKED AS THERE SEAMS TO BE LITTLE DENSIFICATION IN THE MODEL OUTOUT [MAYBE COMPATION IS COMPNSATED FOR BY TRACES OF SNOW???]

%% FUNCTION INFO

% Author: Alex Gardner, University of Alberta
% Date last modified: FEB, 2008 

% Description: 
%   computes the densification of snow/firn using the emperical model of
%   Herron and Langway (1980) or the semi-emperical model of Anthern et al.
%   (2010)

% Inputs:
%   denIdx = densification model to use:
%       1 = emperical model of Herron and Langway (1980)
%       2 = semi-imerical model of Anthern et al. (2010)
%       3 = physical model from Appendix B of Anthern et al. (2010)
%   d   = initial snow/firn density [kg m-3]
%   T   = temperature [K]
%   dz  = grid cell size [m]
%   C   = average accumulation rate [kg m-2 yr-1]
%   dt  = time lapsed [s]
%   re  = effective grain radius [mm];
%   Ta  = mean annual temperature                                          

% Reference: 
% Herron and Langway (1980), Anthern et al. (2010)

%% FOR TESTING
% denIdx = 2;
% d = 800;
% T = 270;
% dz = 0.005;
% C = 200;
% dt = 60*60;
% re = 0.7;
% Tmean = 273.15-18;

%% MAIN FUNCTION
% specify constants
dIce    = 910;         % density of ice [kg m-3]
dt      = dt / 86400;  % convert from [s] to [d]
% R     = 8.314        % gas constant [mol-1 K-1]
% Ec    = 60           % activation energy for self-diffusion of water
%                      % molecules through the ice tattice [kJ mol-1]
% Eg    = 42.4         % activation energy for grain growth [kJ mol-1]


% initial mass
mass_init = d .* dz;

% calculate new snow/firn density for:
%   snow with densities <= 550 [kg m-3]
%   snow with densities > 550 [kg m-3]
idx = d <= 550;
switch denIdx
    case 1 % Herron and Langway (1980)
        c0 = (11 * exp(-10160 ./ (T(idx) * 8.314))) .* C/1000;
        c1 = (575 * exp(-21400 ./ (T(~idx)* 8.314))) .* (C/1000)^0.5;
        
    case 2 % Arthern et al. (2010) [semi-emperical]
        % common variable
        % NOTE: Ec=60000, Eg=42400 (i.e. should be in J not kJ)
        H = exp((-60000./(T * 8.314)) + (42400./(Tmean .* 8.314))) ...     
            .* (C * 9.81);
        
        c0 = 0.07 * H(idx);
        c1 = 0.03 * H(~idx);
        
    case 3 % Arthern et al. (2010) [physical model eqn. B1]
        % calcualte overburden pressure
        obp =[0;  (cumsum(dz(1:end-1)) .* d(1:end-1))];
        
        % common variable
        H = exp((-60./(T * 8.314))) .* obp ./ (re/1000).^2;
        c0 = 9.2E-9 * H(idx);
        c1 = 3.7E-9 * H(~idx);
        
    case 4 % Li and Zwally (2004)
        c0 = (C/dIce) * (139.21 - 0.542*Tmean)*8.36*(273.15 - T) ^ -2.061;
        c1 = c0;
        
    case 5 % Helsen et al. (2008)
        % common variable
        c0 = (C/dIce) * (76.138 - 0.28965*Tmean)*8.36*(273.15 - T) ^ -2.061;
        c1 = c0;
end

% new snow density
d(idx) = d(idx) + (c0 .* (dIce - d(idx)) / 365 * dt);
d(~idx) = d(~idx) + (c1 .* (dIce - d(~idx)) / 365 * dt);

%disp((num2str(nanmean(c0 .* (dIce - d(idx)) / 365 * dt))))

% do not allow densities to exceed the density of ice
d(d > dIce) = dIce;

% calculate new grid cell length
dz = mass_init ./ d;


