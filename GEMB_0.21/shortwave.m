function swf = shortwave(swIdx, aIdx, dsw, as, d, dz, re)

%% DISTRIBUTES ABSORBED SHORTWAVE RADIATION WITHIN SNOW/ICE

% swIdx = 0 : all absorbed SW energy is assigned to the top grid cell

% swIdx = 1 : absorbed SW is distributed with depth as a function of:
%   1 : snow density (taken from Bassford, 2004)
%   2 : grain size in 3 spectral bands (Brun et al., 1992)

%% Inputs
%   swIdx   = shortwave allowed to penetrate surface (0 = No, 1 = Yes)
%   aIdx    = method for calculating albedo (1-4)
%   dsw     = downward shortwave radiative flux [w m-2]
%   as      = surface albedo
%   d       = grid cell density [kg m-3]
%   dz      = grid cell depth [m]
%   re      = grid cell effective grain radius [mm]

%% Outputs
%   swf     = absorbed shortwave radiation [W m-2]

%% SHORTWAVE FUNCTION
% initialize variables
m = length(d);
swf = zeros(m,1);

if swIdx == 0 % all sw radation is absorbed in by the top grid cell
        
    % calculate surface shortwave radiation fluxes [W m-2]
    swf(1) = (1 - as) * dsw;
    
else % sw radation is absorbed at depth within the glacier
    
    if aIdx == 2    % function of effective radius (3 spectral bands)
        
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
        a1 = min(0.98, 1 - 1.58 *gsz(1)^0.5);
        % 0.8 - 1.5um
        a2 = max(0, 0.95 - 15.4 *gsz(1)^0.5);
        % 1.5 - 2.8um
        a3 = max(0.127, 0.88 + 346.3*gsz(1) - 32.31*gsz(1)^0.5);
        
        % seperate net shortwave radiative flux into spectral ranges
        swfS = (sF * dsw) .* (1 - [a1; a2; a3]);
        
        % absorption coeficient for spectral range:
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
        % and Konzelmann 1994 used SWs = 0.36 and SWss = 0.64 as this the
        % the % of SW radiation with wavelengths > and < 800 nm
        % respectively.  This, however, may not account for the fact that
        % the albedo of wavelengths > 800 nm has a much lower albedo.
        
        % calculate surface shortwave radiation fluxes [W m-2]
        swf_s = SWs * (1 - as) * dsw;
        
        % calculate surface shortwave radiation fluxes [W m-2]
        swf_ss = (1-SWs) * (1 - as) * dsw;
        
        % SW allowed to penetrate into snowpack
        Bs = 10;    % snow SW extinction coefficient [m-1] (Bassford,2006)
        Bi = 1.3;   % ice SW extinction coefficient [m-1] (Bassford,2006)
        
        % calculate extinction coefficient B [m-1] vector
        B = Bs + (300 - d) .* ((Bs - Bi)/(910 - 300));
        
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