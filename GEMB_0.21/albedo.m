function a = albedo(aIdx, re, d, n, aIce, aSnow, a, T, W, P, EC,...
    t0wet, t0dry, K, dt)

%% Calculates Snow, firn and ice albedo as a function of:
%   1 : effective grain radius (Gardner & Sharp, 2009)
%   2 : effective grain radius (Brun et al., 2009)
%   3 : density and cloud amount (Greuell & Konzelmann, 1994)
%   4 : exponential time decay & wetness (Bougamont & Bamber, 2005)

%% Inputs
% aIdx      = albedo method to use

% Methods 1 & 2
%   re      = surface effective grain radius [mm]

% Methods 3
%   d       = snow surface density [kg m-3]
%   n       = cloud amount
%   aIce    = albedo of ice
%   aSnow   = albedo of fresh snow

% Methods 4
%   aIce    = albedo of ice
%   aSnow   = albedo of fresh snow
%   a       = grid cell albedo from prevous time step;
%   T       = grid cell temperature [k]
%   W       = pore water [kg]
%   P       = precipitation [mm w.e.] or [kg m-3]
%   EC      = surface evaporation (-) condensation (+) [kg m-2]
%   t0wet   = time scale for wet snow (15-21.9) [d]
%   t0dry   = warm snow timescale [15] [d]
%   K       = time scale temperature coef. (7) [d]
%   dt      = time step of input data [s]

%% Output
%   a       = grid cell albedo 

%% Usage 
% Method 1
% a = albedo(1, 0.1); 

% Method 4
% a = albedo(4, [], [], [], 0.48, 0.85, [0.8 0.5 ... 0.48], ...
%   [273 272.5 ... 265], [0 0.001 ... 0], 0, 0.01, 15, 15, 7, 3600)

%% Function
switch aIdx
    case 1 % function of effective grain radius
        
        % convert effective radius to specific surface area [cm2 g-1]
        S = 3 ./ (.091 * re);
        
        % determine broadband albedo
        a(1) = 1.48 - S^-0.07;
        
    case 2 % function of effective grain radius
        % Spectral fractions  (Lefebre et al., 2003)
        % [0.3-0.8um 0.8-1.5um 1.5-2.8um]
        sF = [0.606 0.301 0.093];
        
        % convert effective radius to grain size in meters
        gsz = (re * 2) / 1000;
        
        % spectral range:
        % 0.3 - 0.8um
        a1 = min(0.98, 1 - 1.58 *gsz^0.5);
        % 0.8 - 1.5um
        a2 = max(0, 0.95 - 15.4 *gsz^0.5);
        % 1.5 - 2.8um
        a3 = max(0.127, 0.88 + 346.3*gsz - 32.31*gsz^0.5);
        
        % broadband surface albedo
        a(1) = sF * [a1; a2; a3];
        
    case 3 % a as a function of density
        
        % define constants
        dIce = 910;    % density of ice [kg m-3]
        dSnow = 300;   % density of fresh snow [kg m-3]       
        
        % calculate albedo
        a(1) = aIce + (d - dIce)*(aSnow - aIce) ...
            / (dSnow - dIce) + (0.05 * (n - 0.5));
        
    case 4 % exponential time decay & wetness
        
        % set variables
        dSnow = 300;       %[kg m-3]
        
        % change in albedo with time:
        %   (d_a) = (a - a_old)/(t0)
        % where: t0 = timescale for albedo decay
        
        dt = dt / 86400;    % convert from [s] to [d]
        
        % initialize variables
        t0 = zeros(size(a));
        
        % specify constants
        % a_wet = 0.15;         % water albedo (0.15)
        % a_new = aSnow        % new snow albedo (0.64 - 0.89)
        % a_old = aIce;        % old snow/ice albedo (0.27-0.53)
        % t0_wet = t0wet;      % time scale for wet snow (15-21.9) [d]
        % t0_dry = t0dry;      % warm snow timescale [15] [d]
        % K = 7                 % time scale temperature coef. (7) [d]
        % W0 = 300;             % 200 - 600 [mm]
        z_snow = 15;            % 16 - 32 [mm]
        
        % determine timescale for albedo decay
        t0(W > 0) = t0wet;                     % wet snow timescale
        T = T - 273.15;                         % change T from K to °C
        t0warm = abs(T) * K + t0dry;          % 'warm' snow timescale
        t0(W == 0 & T >= -10) = ...
            t0warm(W == 0 & T >= -10);
        t0(T < -10,1) =  10 * K + t0dry;       % 'cold' snow timescale
        
        % calculate new albedo
        d_a = (a - aIce) ./ t0 * dt;           % change in albedo
        a = a - d_a;                            % new albedo
        
        % modification of albedo due to thin layer of snow or solid
        % condensation (deposition) at the surface surface
        
        % check if condensation occurs & if it is deposited in solid phase
        if ( EC > 0 && T(1) < 0)
            P = P + (EC/dSnow) * 1000;  % add cond to precip [mm]
        end
        
        a(1) = aSnow - (aSnow - a(1)) * exp(-P/z_snow);
        
        %----------THIS NEEDS TO BE IMPLEMENTED AT A LATER DATE------------
        % modification of albedo due to thin layer of water on the surface
        % a_surf = a_wet - (a_wet - a_surf) * exp(-W_surf/W0);
end

%% Check for erroneous values
if a(1) > 1
    warning ('albedo > 1.0')
elseif a(1) < 0
    warning ('albedo is negative')
elseif isnan(a(1))
    error ('albedo == NAN')
end

end