function [a, adiff] = albedo(aIdx, re, dz, d, n, aIce, aSnow, aValue, adThresh, a, adiff, TK, W, P, EC,...
    Ms, clabSnow, clabIce, SZA, COT, t0wet, t0dry, K, dt, dIce)
% albedo calculates snow, firn and ice albedo as a function of:
%   1 : effective grain radius (Gardner & Sharp, 2009)
%   2 : effective grain radius (Brun et al., 2009)
%   3 : density and cloud amount (Greuell & Konzelmann, 1994)
%   4 : exponential time decay & wetness (Bougamont & Bamber, 2005)
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
% aIdx      = albedo method to use
%
% Method 0
%  aValue   = direct input value for albedo, override all changes to albedo
%
% adThresh
%  Apply below method to all areas with densities below this value,
%  or else apply direct input value, allowing albedo to be altered.
%
% Methods 1 & 2
%   re      = surface effective grain radius [mm]
% Method 1, optional
%  clabSnow = concentration of light absorbing carbon  [ppm1], default 0
%  SZA      = solar zenith angle of the incident radiation [deg], default 0
%  COT      = cloud optical thickness, default 0
%  For TWO LAYER
%  clabIce  = concentration of light absorbing carbon of first ice layer [ppm1], default 0
%
% Method 3
%   d       = snow surface density [kg m-3]
%   n       = cloud amount
%   aIce    = albedo of ice
%   aSnow   = albedo of fresh snow
%
% Method 4
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
%% Outputs
%
%  asdiff  = surface albedo for diffuse radiation
%
%% Documentation
%
% For complete documentation, see: https://github.com/alex-s-gardner/GEMB
%
%% References
% If you use GEMB, please cite the following:
%
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass
% Balance (GEMB): a model of firn processes for cryosphere research, Geosci.
% Model Dev., 16, 2277–2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023.

%% Usage
% Method 1
% a = albedo(1, 0.1);

% Method 4
% a = albedo(4, [], [], [], 0.48, 0.85, [0.8 0.5 ... 0.48], ...
%   [273 272.5 ... 265], [0 0.001 ... 0], 0, 0.01, 15, 15, 7, 3600)

Ttol = 1e-10;
Dtol = 1e-11;
Wtol = 1e-13;

dSnow  = 300.0; % density of fresh snow [kg m-3]
dPHC   = 830.0; % Pore closeoff density
ai_max = 0.58;  % maximum ice albedo, from Lefebre,2003
ai_min = aIce;  % minimum ice albedo
as_min = 0.65;  % minimum snow albedo, from Alexander 2014

%% Function
if(aIdx==0 || (adThresh - d(1))<Dtol)
    a(1) = aValue;
else
    switch aIdx
        case 1 % function of effective grain radius

            % clabSnow, IssmDouble clabIce, IssmDouble SZA, IssmDouble COT, int m
            a(1)    =gardnerAlb(re, dz, d, clabSnow, clabIce,  SZA, COT);
            adiff(1)=gardnerAlb(re, dz, d, clabSnow, clabIce, 50.0, COT);

        case 2 % function of effective grain radius
            % Spectral fractions  (Lefebre et al., 2003)
            % [0.3-0.8um 0.8-1.5um 1.5-2.8um]
            sF = [0.606 0.301 0.093];

            % convert effective radius to grain size in meters
            gsz = (re(1) * 2.) / 1000.;

            % spectral range:
            % 0.3 - 0.8um
            a1 = min(0.98, 0.95 - 1.58 *gsz^0.5);
            % 0.8 - 1.5um
            a2 = max(0., 0.95 - 15.4 *gsz^0.5);
            % 1.5 - 2.8um
            a3 = max(0.127, 0.88 + 346.3*gsz - 32.31*gsz^0.5);

            % broadband surface albedo
            a(1) = sF * [a1; a2; a3];

        case 3 % a as a function of density

            % calculate albedo
            a(1) = aIce + (d(1) - dIce)*(aSnow - aIce) ...
                / (dSnow - dIce) + (0.05 * (n - 0.5));

        case 4 % exponential time decay & wetness

            % change in albedo with time:
            %   (d_a) = (a - a_old)/(t0)
            % where: t0 = timescale for albedo decay

            dt = dt / 86400;    % convert from [s] to [d]

            % initialize variables
            t0 = zeros(size(a));

            % specify constants
            % a_wet = 0.15;        % water albedo (0.15)
            % a_new = aSnow        % new snow albedo (0.64 - 0.89)
            % a_old = aIce;        % old snow/ice albedo (0.27-0.53)
            % t0_wet = t0wet;      % time scale for wet snow (15-21.9) [d]
            % t0_dry = t0dry;      % warm snow timescale [15] [d]
            % K = 7                % time scale temperature coef. (7) [d]
            % W0 = 300;            % 200 - 600 [mm]
            z_snow = 15;           % 16 - 32 [mm]

            % determine timescale for albedo decay
            t0(W > 0+Wtol) = t0wet;               % wet snow timescale
            T = TK - 273.15;                       % change T from K to °C
            t0warm = abs(T) * K + t0dry;          % 'warm' snow timescale
            t0(abs(W)<Wtol & T >= -10-Ttol) = t0warm(abs(W)<Wtol & T >= -10-Ttol);
            t0(T < -10-Ttol,1) =  10 * K + t0dry; % 'cold' snow timescale

            % calculate new albedo
            d_a = (a - aIce) ./ t0 * dt;           % change in albedo
            a = a - d_a;                            % new albedo

            % modification of albedo due to thin layer of snow or solid
            % condensation (deposition) at the surface surface

            % check if condensation occurs & if it is deposited in solid phase
            if ( EC > 0+Dtol && T(1) < 0-Ttol)
                P = P + (EC/dSnow) * 1000;  % add cond to precip [mm]
            end

            a(1) = aSnow - (aSnow - a(1)) * exp(-P./z_snow);

    end

    %If we do not have fresh snow
    if (aIdx<3 && aIdx>0 && (adThresh - d(1))>=Dtol)
        % In a snow layer < 10cm, account for mix of ice and snow,
        % after P. Alexander et al., 2014
        lice = find([d; 999]>=dPHC-Dtol);
        depthsnow = sum(dz(1:(lice(1)-1)));

        if (depthsnow<=0.1+Dtol & lice(1)<=length(d) & d(lice(1))>=dPHC-Dtol)
            aice = ai_max + (as_min - ai_max)*(d(lice(1))-dIce)/(dPHC-dIce);
            a(1) = aice + max(a(1)-aice,0.0)*(depthsnow/0.1);
        end

        if (d(1)>=dPHC-Dtol)
            if (d(1)<dIce-Dtol) %For continuity of albedo in firn i.e. P. Alexander et al., 2014

                %ai=ai_max + (as_min - ai_max)*(dI-dIce)/(dPHC-dIce);
                %dPHC is pore close off (830 kg m^-3)
                %dI is density of the upper firn layer
                a(1) = ai_max + (as_min - ai_max)*(d(1)-dIce)/(dPHC-dIce);

            else %surface layer is density of ice

                %When density is > dIce (typically 910 kg m^-3, 920 is used by Alexander in MAR),
                %ai=ai_min + (ai_max - ai_min)*e^(-1*(Msw(t)/K))
                %K is a scale factor (set to 200 kg m^-2)
                %Msw(t) is the time-dependent accumulated amount of excessive surface meltwater
                %  before run-off in kg m^-2 (melt per GEMB timestep, i.e. 3 hourly)
                M = Ms+W(1);
                a(1)=max(ai_min + (ai_max - ai_min)*exp(-1.0*(M/200.0)), ai_min);
            end
        end
    end
end

%% Check for erroneous values

if a(1) > 1+Ttol
    warning ('albedo > 1.0')
elseif a(1) < 0-Dtol
    warning ('albedo is negative')
elseif isnan(a(1))
    error ('albedo == NAN')
end

end