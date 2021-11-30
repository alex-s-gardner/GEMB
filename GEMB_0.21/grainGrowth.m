function [re, gdn, gsp]  = grainGrowth(T, dz, d, W, re, gdn, gsp, dt)

%% Function Documentation
% *Created by*: Alex S. Gardner, University of Alberta
%
% *Description*: models the effective snow grain size
%
% *Reference:*
% DENDRITIC SNOW METAMORPHISM:
% Brun, E., P. David, M. Sudul, and G. Brunot, 1992: A numerical model to
% simulate snow-cover stratigraphy for operational avalanche forecasting.
% Journal of Glaciology, 38, 13-22.
%
% NONDENDRITIC SNOW METAMORPHISM:
% Dry snow metamorphism:
% Marbouty, D., 1980: An experimental study of temperature-gradient
% metamorphism. Journal of Glaciology, 26, 303-312.
%
% WET SNOW METAMORPHISM:
% Brun, E., 1989: Investigation on wet-snow metamorphism in respect of
% liquid-water content. Annals of Glaciology, 13, 22-26.

%% INPUTS
% * T: grid cell temperature [k]
% * dz: grid cell depth [m]
% * d: grid cell density [kg m-3]
% * W: water content [kg]
% * re: effective grain size [mm]
% * gdn: grain dentricity
% * gsp: grain sphericity
% * dt: time step of input data [s]

%% OUTPUTS
% * re: effective grain size [mm]
% * gdn: grain dentricity
% * gsp: grain sphericity

%% Function
gsz = re * 2;

% convert dt from seconds to days
dt = dt/86400;

% convert T from k to deg C

% determine liquied-water content in terms
lwc = W ./ (d .* dz) * 100;

% set maximum water content by mass to 9 percent (Brun, 1980)
lwc(lwc > 9) = 9;

%% Calculate temperature gradiant across grid cells
% returns the avereage gradinet across the upper and lower grid cell

% initialize
dT = zeros(size(T));

% depth of grid point center from surface
zGPC = (cumsum(dz) - dz/2);

% Take forward differences on left and right edges
dT(1,1) = (T(2,1) - T(1,1))/(zGPC(2)-zGPC(1));
dT(end,1) = (T(end,1) - T(end-1,1))/(zGPC(end)-zGPC(end-1));

% Take centered differences on interior points
zGPC = zGPC(3:end) - zGPC(1:end-2);
dT(2:end-1,1) = (T(3:end,1)-T(1:end-2,1))./zGPC(:,1);

% take absolute value of temperature gradient
dT = abs(dT);

% index for dentricity > 0 & == 0
G  = gdn > 0;
J = ~G;

%% DENDRITIC SNOW METAMORPHISM
% FOR SNOW DENTRICITY > 0

% if there is snow dentricity > 0
if sum(G) ~= 0
%    disp ('DENDRITIC DRY SNOW METAMORPHISM')
    % index for dentricity > 0 and T gradients < 5 ºC m-1 and >= 5 ºC m-1
    H = dT < 5 & G & W == 0; % asg not wet accounted for on 19/08/29
    I = dT >= 5 & G & W == 0; % asg not wet accounted for on 19/08/29
    
    % determine coefficients
    A = - 2E8 * exp(-6E3 ./ T(H)) * dt;
    B = 1E9 * exp(-6E3 ./ T(H)) * dt;
    C = (-2E8 * exp(-6E3 ./ T(I)) * dt) .* (dT(I) .^ 0.4);
    
    % new dendricity and sphericity for dT < 5 ºC m-1
    gdn(H) = gdn(H) + A;
    gsp(H) = gsp(H) + B;
    
    % new dendricity and sphericity for dT >= 5 ºC m-1
    gdn(I) = gdn(I) + C;
    gsp(I) = gsp(I) + C;
    
    % WET SNOW METAMORPHISM
    
    % index for dendritic wet snow
    L = W > 0 & G;
    
    % check if snowpack is wet
    if sum(L) ~= 0
%        disp('DENDRITIC WET SNOW METAMORPHISM')
        % determine coefficient
        D = (1/16) * (lwc(L) .^ 3) * dt;
        
        % new dendricity and sphericity for wet snow
        gdn(L) = gdn(L) - D;
        gsp(L) = gsp(L) + D;
    end
    
    % dendricity and sphericity can not be > 1 or < 0
    gdn(gdn < 0) = 0;
    gsp(gsp < 0) = 0;
    gdn(gdn > 1) = 1;
    gsp(gsp > 1) = 1;
    
    % determine new grain size (mm)
    gsz(G) = 0.1 + (1-gdn(G))*0.25 + (0.5-gsp(G))*0.1;
end

% if there is snow dentricity == 0
if sum(J) ~= 0
%    disp('NONDENDRITIC SNOW METAMORPHISM')

    % DRY SNOW METAMORPHISM (Marbouty, 1980)
    % grouped model coefficinets from Marbouty, 1980: Figure 9
    P = J & W == 0; % asg not wet accounted for on 19/08/29
    Q = Marbouty(T(P), d(P), dT(P));
 
    % calculate grain growth
    gsz(P) = gsz(P) + Q * dt;
    
    % WET SNOW METAMORPHISM (Brun, 1989)
    
    % index for nondendritic wet snow
    K = W > 0 & J;
    
    % check if snowpack is wet
    if sum(K) ~= 0
%        disp('NONDENDRITIC WET SNOW METAMORPHISM')
        % wet rate of change coefficient
        E = 1.28E-8 + (4.22E-10 * (lwc(K).^3))* (dt *86400);   % [mm^3 s^-1]
        
        % calculate change in grain volume and convert to grain size
        gsz(K) = 2 * (3/(pi * 4)*((4 / 3)*pi*(gsz(K)/2).^3 + E)).^(1/3);
    end
    
    % grains with sphericity == 1 can not have grain sizes > 2 mm (Brun, 1992)
    gsz(gsp == 1 & gsz > 2) = 2;
    
    % grains with sphericity == 0 can not have grain sizes > 5 mm (Brun, 1992)
    gsz(gsp ~= 1 & gsz > 5) = 5;
end

% convert grain size back to effective grain radius
re = gsz ./2;

end

function Q = Marbouty(T, d, dT)
%% calculates grain growth according to Fig. 9 of Marbouty, 1980
% ------NO GRAIN GROWTH FOR d > 400 kg m-3 because (H is set to zero)------
% ---------------this is a major limitation of the model-------------------
% initialize
F = zeros( size(T));
H = F;
G = F;

% convert T from K to ºC
T = T - 273.15;

% temperature coefficient F
I = T > -6;
F(I) =  0.7 + ((T(I)/-6) * 0.3);
I = T <= -6 & T > -22;
F(I) =  1 - ((T(I)+6)/-16 * 0.8);
I = T <= -22 & T > -40;
F(I) =  0.2 - ((T(I)+22)/-18 * 0.2);

% density coefficient F
H(d < 150) = 1;
I = d >= 150 & d < 400;
H(I) = 1 - ((d(I)-150)/250);

% temperature gradient coefficient G
I = dT >= 0.16 & dT < 0.25;
G(I) = ((dT(I) - 0.16)/0.09) * 0.1;
I = dT >= 0.25 & dT < 0.4;
G(I) = 0.1 + (((dT(I) - 0.25)/0.15) * 0.57);
I = dT >= 0.4 & dT < 0.5;
G(I) = 0.67 + (((dT(I) - 0.4)/0.1) * 0.23);
I = dT >= 0.5 & dT < 0.7;
G(I) = 0.9 + (((dT(I) - 0.5)/0.2) * 0.1);
G(dT >= 0.7) = 1;

% model time growth constat E
E = 0.09;        %[mm d-1]

% grouped coefficinet Q
Q = F.*H.*G.*E;
end
