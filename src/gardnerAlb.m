function [a] = gardnerAlb(re, dz, d, c1, c2, SZA, t)
% gardnerAlb is a Matlab implementation of the snow and ice broadband albedo 
% parameterization developed by Alex Gardner.
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
% ONE LAYER
%   - re    : effective radius [mm] to calculate S1 and S2 (specific surface area of the snow or ice [cm^2 g-1])
%   - c1    : concentration of light absorbing carbon  [ppmw]
%   - c2    : concentration of light absorbing carbon of bottom ice layer [ppmw]
%   - SZA   : solar zenith angle of the incident radiation [deg]
%   - t     : cloud optical thickness
%
% TWO LAYER
%   - z1    : depth of snow suface layer [mm w.e.]
%   - c2    : concentration of light absorbing carbon of bottom ice
%             layer [ppmw]
% 
%% Outputs
% 
% 
%% Documentation
% 
% For complete documentation, see: https://github.com/alex-s-gardner/GEMB 
% 
%% References 
% The formulations in this function are described in: 
% 
% Gardner, A. S., and Sharp, M. J.: A review of snow and ice albedo and the 
% development of a new physically based broadband albedo parameterization, 
% J. Geophys. Res., 115, F01009, 10.1029/2009jf001444, 2010. 
% 
% If you use GEMB, please cite the following: 
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass 
% Balance (GEMB): a model of firn processes for cryosphere research, Geosci. 
% Model Dev., 16, 2277–2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023.

%% Single layer albedo parameterization

Dtol = 1e-11;

%convert effective radius to specific surface area [cm2 g-1]
S1 = 3.0 / (0.091 * re(1));

% effective solar zenith angle
x = min((t./(3*cos(pi*SZA/180))).^0.5, 1);
u = 0.64*x + (1-x).* cos(pi*SZA/180);

% pure snow albedo
as = (1.48 - S1.^-0.07);
 
% change in pure snow albedo due to soot loading
dac = max(0.04 - as, ...
    ((-c1).^0.55)./(0.16 + 0.6*S1.^0.5 + (1.8*c1.^0.6).*(S1.^-0.25)));
 
% Two layer albedo parameterization
%   do two layer calculation if there is more than 1 layer
lice = find([d; 999]>=830-Dtol);
z1   = sum(dz(1:(lice(1)-1)).*d(1:(lice(1)-1)));

m=length(d);
if (m>0 & lice(1)<=m & z1 > Dtol)
    
    % determine albedo values for bottom layer
    S2 = 3.0 / (0.091 * re(lice(1)));

    % pure snow albedo
    as2 = (1.48 - S2.^-0.07);
    
    % change in pure snow albedo due to soot loading
    dac2 = max(0.04 - as2, ...
        ((-c2).^0.55)./(0.16 + 0.6*S2.^0.5 + (1.8*c2.^0.6).*(S2.^-0.25)));
    
    % determine the effective change due to finite depth and soot loading
    A = min(1, (2.1 * z1 .^(1.35*(1-as) - 0.1*c1 - 0.13)));
    
    dac =  (as2 + dac2 - as) + A .* ((as + dac) - (as2 + dac2));
end
 
% change in albedo due to solar zenith angle
dasz = 0.53 * as .* (1 - (as + dac)) .* (1 - u) .^ 1.2;
 
% change in albedo due to cloud (apart from change in diffuse fraction)
dat = (0.1 * t .* (as + dac).^ 1.3) ./ ((1 + 1.5*t).^as);
 
%% Broadband albedo 

a = as + dac + dasz + dat;

end