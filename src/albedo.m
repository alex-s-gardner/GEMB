function [a, a_diffuse] = albedo(T, dz, d, W, re, a, a_diffuse, dt, P, EC, ...
    M_surf, density_ice, black_carbon_snow, black_carbon_ice, ...
    solar_zenith_angle, cloud_optical_thickness, cloud_fraction, ...
    albedo_method, albedo_ice, albedo_snow, albedo_fixed, ...
    albedo_desnity_threshold, albedo_wet_snow_t0, albedo_dry_snow_t0, ...
    albedo_K)

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
% albedo_method             = albedo method to use ["None", "GardnerSharp", "BruneLeFebre", "GreuellKonzelmann", "BougamontBamber"]
%
% Method 0
%  albedo_fixed             = direct input value for albedo, override all changes to albedo
%
% albedo_desnity_threshold
%  Apply below method to all areas with densities below this value,
%  or else apply direct input value, allowing albedo to be altered.
%
% Methods 1 & 2
%  re                      = surface effective grain radius [mm]
% Method 1, optional
%  black_carbon_snow        = concentration of light absorbing carbon  [ppm1], default 0
%  solar_zenith_angle       = solar zenith angle of the incident radiation [deg], default 0
%  cloud_optical_thickness  = cloud optical thickness, default 0
%  For TWO LAYER
%  black_carbon_ice         = concentration of light absorbing carbon of first ice layer [ppm1], default 0
%
% Method 3
%   d                    = snow surface density [kg m-3]
%   cloud_fraction       = cloud amount
%   albedo_ice           = albedo of ice
%   albedo_snow          = albedo of fresh snow
%
% Method 4
%   albedo_ice           = albedo of ice
%   albedo_snow          = albedo of fresh snow
%   a                    = grid cell albedo from prevous time step;
%   T                    = grid cell temperature [k]
%   W                    = pore water [kg]
%   P                    = precipitation [mm w.e.] or [kg m-3]
%   EC                   = surface evaporation (-) condensation (+) [kg m-2]
%   albedo_wet_snow_t0   = time scale for wet snow (15-21.9) [d]
%   albedo_dry_snow_t0   = warm snow timescale [15] [d]
%   albedo_K             = time scale temperature coef. (7) [d]
%   dt                   = time step of input data [s]
%
%% Outputs
%
%  a_diffuse  = surface albedo for diffuse radiation
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

T_tolerance  = 1e-10;
d_tolerance  = 1e-11;
W_tolerance  = 1e-13;

% constants
CtoK                = 273.15;      % Celsius to Kelvin conversion
density_fresh_snow  = 300.0;       % density of fresh snow [kg m-3]
density_phc         = 830.0;       % Pore closeoff density
albedo_ice_max      = 0.58;        % maximum ice albedo, from Lefebre,2003
albedo_ice_min      = albedo_ice;  % minimum ice albedo
albedo_snow_min     = 0.65;        % minimum snow albedo, from Alexander 2014

%% Function
if (albedo_method == "None") || ((albedo_desnity_threshold - d(1)) < d_tolerance)
    a(1) = albedo_fixed;
else
    switch albedo_method
        case "GardnerSharp" % function of effective grain radius

            % black_carbon_snow, IssmDouble black_carbon_ice, IssmDouble solar_zenith_angle, IssmDouble cloud_optical_thickness, int m
            a(1)         = albedo_gardner(re, dz, d, black_carbon_snow, black_carbon_ice,  solar_zenith_angle, cloud_optical_thickness);
            a_diffuse(1) = albedo_gardner(re, dz, d, black_carbon_snow, black_carbon_ice, 50.0, cloud_optical_thickness);

        case "BruneLeFebre" % function of effective grain radius
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

        case "GreuellKonzelmann" % a as a function of density

            % calculate albedo
            a(1) = albedo_ice + (d(1) - density_ice)*(albedo_snow - albedo_ice) ...
                / (density_fresh_snow - density_ice) + (0.05 * (cloud_fraction - 0.5));

        case "BougamontBamber" % exponential time decay & wetness

            % change in albedo with time:
            %   (d_a) = (a - a_old)/(t0)
            % where: t0 = timescale for albedo decay

            dt = dt / 86400;    % convert from [s] to [d]

            % initialize variables
            t0 = zeros(size(a));

            % specify constants
            % a_wet = 0.15;        % water albedo (0.15)
            % a_new = albedo_snow        % new snow albedo (0.64 - 0.89)
            % a_old = albedo_ice;        % old snow/ice albedo (0.27-0.53)
            % t0_wet = albedo_wet_snow_t0;      % time scale for wet snow (15-21.9) [d]
            % t0_dry = albedo_dry_snow_t0;      % warm snow timescale [15] [d]
            % albedo_K = 7                % time scale temperature coef. (7) [d]
            % W0 = 300;            % 200 - 600 [mm]
            z_snow = 15;           % 16 - 32 [mm]

            % determine timescale for albedo decay
            t0(W > 0+W_tolerance ) = albedo_wet_snow_t0;               % wet snow timescale
            TC = T - CtoK;                       % change T from K to °C
            t0warm = abs(TC) * albedo_K + albedo_dry_snow_t0;          % 'warm' snow timescale
            t0(abs(W)<W_tolerance  & TC >= -10-T_tolerance ) = t0warm(abs(W)<W_tolerance  & TC >= -10-T_tolerance );
            t0(TC < -10-T_tolerance ,1) =  10 * albedo_K + albedo_dry_snow_t0; % 'cold' snow timescale

            % calculate new albedo
            d_a = (a - albedo_ice) ./ t0 * dt;      % change in albedo
            a = a - d_a;                            % new albedo

            % modification of albedo due to thin layer of snow or solid
            % condensation (deposition) at the surface surface

            % check if condensation occurs & if it is deposited in solid phase
            if ( EC > 0+d_tolerance  && TC(1) < 0-T_tolerance )
                P = P + (EC/density_fresh_snow) * 1000;  % add cond to precip [mm]
            end

            a(1) = albedo_snow - (albedo_snow - a(1)) * exp(-P./z_snow);

    end

    %If we do not have fresh snow
    if ismember(albedo_method,["GreuellKonzelmann","BougamontBamber"]) && ismember(albedo_method,["BruneLeFebre", "GreuellKonzelmann", "BougamontBamber"]) && ((albedo_desnity_threshold - d(1)) >= d_tolerance)
        % In a snow layer < 10cm, account for mix of ice and snow,
        % after P. Alexander et al., 2014
        lice = find([d; 999]>=density_phc-d_tolerance );
        depthsnow = sum(dz(1:(lice(1)-1)));

        if (depthsnow <= (0.1+d_tolerance))  && (lice(1) <= length(d)) && (d(lice(1)) >= (density_phc-d_tolerance))
            aice = albedo_ice_max + (albedo_snow_min - albedo_ice_max) * (d(lice(1)) - density_ice) / (density_phc - density_ice);
            a(1) = aice + max(a(1) - aice,0.0) * (depthsnow/0.1);
        end

        if (d(1) >= density_phc-d_tolerance)
            if (d(1) < density_ice-d_tolerance) %For continuity of albedo in firn i.e. P. Alexander et al., 2014

                %ai=albedo_ice_max + (albedo_snow_min - albedo_ice_max)*(dI-density_ice)/(density_phc-density_ice);
                %density_phc is pore close off (830 kg m^-3)
                %dI is density of the upper firn layer
                a(1) = albedo_ice_max + (albedo_snow_min - albedo_ice_max)*(d(1)-density_ice)/(density_phc-density_ice);

            else %surface layer is density of ice

                %When density is > density_ice (typically 910 kg m^-3, 920 is used by Alexander in MAR),
                %ai=albedo_ice_min + (albedo_ice_max - albedo_ice_min)*e^(-1*(Msw(t)/albedo_K))
                %albedo_K is a scale factor (set to 200 kg m^-2)
                %Msw(t) is the time-dependent accumulated amount of excessive surface meltwater
                %  before run-off in kg m^-2 (melt per GEMB timestep, i.e. 3 hourly)
                M = M_surf+W(1);
                a(1)=max(albedo_ice_min + (albedo_ice_max - albedo_ice_min) * exp(-1.0 * (M / 200.0)), albedo_ice_min);
            end
        end
    end
end

%% Check for erroneous values

if a(1) > (1 + T_tolerance)
    warning ('albedo > 1.0')
elseif a(1) < (0 - d_tolerance)
    warning ('albedo is negative')
elseif isnan(a(1))
    error ('albedo == NAN')
end

end

function a = albedo_gardner(re, dz, d, c1, c2, SZA, t)
% albedo_gardner is a Matlab implementation of the snow and ice broadband albedo 
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

d_tolerance  = 1e-11;

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
lice = find([d; 999]>=830-d_tolerance );
z1   = sum(dz(1:(lice(1)-1)).*d(1:(lice(1)-1)));

m=length(d);
if (m > 0) && (lice(1) <= m) && (z1 > d_tolerance)
    
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