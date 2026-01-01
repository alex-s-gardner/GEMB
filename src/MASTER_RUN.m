% number of processors to use
mod.NW = 48;
mod.PC = 'off'; %'flux'; %'flux'; % 'on';
TEST = true;

%% GEMB INITIALIZATION
% unique model run ID to save output as
S.runPfx = 'S2A1D2';

% spin-up
S.spinUp = 2;   % number of cycles of met data run before output
% calcualted, set spinUp = 0 for no spin up

% select method of calculating albedo and subsurface absorption (default is 1)
%   0 : direct input from aValue parameter, no use of adThresh
%   1 : effective grain radius (Gardner & Sharp, 2009)
%   2 : effective grain radius (Brun et al., 1992; LeFebre et al., 2003), with swIdx=1, SW penetration follows grain size in 3 spectral bands (Brun et al., 1992) 
%   3 : density and cloud amount (Greuell & Konzelmann, 1994)
%   4 : exponential time decay & wetness (Bougamont & Bamber, 2005)
S.aIdx = 1;
% Apply aIdx method to all areas with densities below this value, or else apply direct input value from aValue, allowing albedo to be altered.
% Default value is rho water (1023 kg m-3).
S.adThresh = 1023;

% apply all SW to top grid cell (0) or allow SW to penetrate surface (1) 
% (default 0: if swIdx=1 and aIdx=2, function of effective radius (Brun et al., 1992) or else dependent on snow density (taken from Bassford, 2002))
S.swIdx = 1;

% select densification model to use (default is 2):
%   1 = emperical model of Herron and Langway (1980)
%   2 = semi-emperical model of Anthern et al. (2010)
%   3 = DO NOT USE: physical model from Appendix B of Anthern et al. (2010)
%   4 = DO NOT USE: emperical model of Li and Zwally (2004)
%   5 = DO NOT USE: modified emperical model (4) by Helsen et al. (2008)
%   6 = Antarctica semi-emperical model of Ligtenberg et al. (2011)
%   7 = Greenland semi-emperical model of Kuipers Munneke et al. (2015)
S.denIdx = 2;

% select model for fresh snow accumulation density (default is 1):
%   0 = Original GEMB value, 150 kg/m^3
%   1 = Antarctica value of fresh snow density, 350 kg/m^3
%   2 = Greenland value of fresh snow density, 315 kg/m^3, Fausto et al. (2018)
%   3 = Antarctica model of Kaspers et al. (2004)
%   4 = Greenland model of Kuipers Munneke et al. (2015)
S.dsnowIdx = 1;

% select method for calculating emissivity (default is 1)
%   0: direct input from teValue parameter, no use of teThresh
%   1: default value of 1, in areas with grain radius below teThresh
%   2: default value of 1, in areas with grain radius below teThresh and areas of dry snow (not bare ice or wet) at the surface
S.eIdx = 1;
% Apply eIdx method to all areas with grain radii above this value (mm), or else apply direct input value from teValue, allowing emissivity to be altered.
% Default value is a effective grain radius of 10 mm.
S.teThresh = 10;

% select method for calculating thermal conductivity (default is 1)
% 1: after Sturm et al, 1997
% 2: after Calonne et al., 2011
S.tcIdx = 1;

% GRID INITIALIZATION
% set depth of top grid cell structure (constant grid length) [m]
S.zTop = 10;

% set initial top vertical and min allowable grid spacings [m]
S.dzTop = 0.05; %0.05;
S.dzMin = S.dzTop / 2;

% set initial/max and min model depth [m]
S.zMax = 250;
S.zMin = ceil(S.zMax/2 /10)*10;

% strech grid cells bellow top_z by a [top_dz * y ^ (cells bellow top_z)]
S.zY = 1.10;

%thermal: scaling factor to multiply the thermal diffusion timestep (delta t)
S.ThermoDeltaTScaling = 1/11;

% mean annual wind velocity [m s-1], climatology
S.Vmean=10.0;

%optional inputs:
S.teValue = 1.0;     % Outward longwave radiation thermal emissivity forcing at every element (default in code is 1).
                     % Used only if eIdx==0, or effective grain radius exceeds teThresh

S.isdeltaLWup=false; % True to perturb the long wave radiation upwards.
S.dulwrfValue = 0.0; % Delta with which to perturb the long wave radiation upwards. Use if isdeltaLWup is true.

S.isrestart=false;   % True if we want to restart from *ini parameters set in S struct

% OTHER
% specify frequency to output data (density, grid length, and temperature)
%   - 'monthly'
%   - 'daily'
S.outputFreq = 'monthly';

% input data directory
if TEST
    S.outDIR = fullfile('..','TEST_DATA/');
else
    % output directory
    S.outDIR = '/Volumes/MasterBrain/data/GEMB/Output/';
end

% FIXED ALBEDO VARIABLES
switch S.aIdx
    case {1 2}
        % albedo tuning parameters
        % for methods of calculating albedo see albedo function
        % METHOD 1 & 2
        S.aSnow = 0.85;         % new snow albedo (0.64 - 0.89)
        S.aIce = 0.48;          % range 0.27-0.58 for old snow
	S.aValue = S.aSnow;  % Albedo forcing at every element.  Used only if aIdx == 0, or density exceeds adThresh

	%Defaut values, but these can also be set as time series forcing
	S.dswdiffrf=0.0;        % downward diffusive shortwave radiation flux [W/m^2]
	S.szaValue=0.0;         % Solar Zenith Angle [degree]
	S.cotValue=0.0;         % Cloud Optical Thickness
	S.ccsnowValue=0.0;      % concentration of light absorbing carbon for snow [ppm1]
	S.cciceValue=0.0;       % concentration of light absorbing carbon for ice [ppm1]
        
    case 3
        % RADIATION CORRECTION FACTORS
        % -> only used for met station data and Greuell & Konzelmann, 1994 albedo
        S.cldFrac = 0.1;        % average cloud amount
        
    case 4
        % additonal tuning parameters albedo as a funtion of age and water content
        % (Bougamont et al., 2005)
        S.t0wet = 15;           % time scale for wet snow (15-21.9) [d]
        S.t0dry = 30;           % warm snow timescale (30) [d]
        S.K = 7;                % time scale temperature coef. (7) [d]
end

%% RUN GEMB
runPfx = S.runPfx;

% open matlab pool for parallel processing
if strcmp(mod.PC, 'on')
    parpool(mod.NW)
elseif strcmp(mod.PC, 'flux')
    parpool('mpiexecflux', mod.NW)
end

if TEST
    [dateN, P0, Ta0, V0, dlw0, dsw0, eAir0, pAir0, LP] = simulate_climate_forcing(set_id);
    S0 = combineStrucData_GEMB(S,LP,1);

    GEMB(P0, Ta0, V0, dateN, dlw0, dsw0, eAir0, pAir0, S0, S.isrestart)
else
    error("input case not defined")
end
