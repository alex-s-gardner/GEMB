% number of processors to use
S.n_workers = 1;
S.run_id    = "test_1";

%% GEMB INITIALIZATION
% unique model run ID to save output as
S.run_prefix = 'S2A1D2';

% spin-up
S.n_spinup_cycles = 2;   % number of cycles of met data run before output
% calcualted, set spinUp = 0 for no spin up

% select method of calculating albedo and subsurface absorption (default is 1)
%   0 : direct input from albedo_fixed parameter, no use of albedo_desnity_threshold
%   1 : effective grain radius (Gardner & Sharp, 2009)
%   2 : effective grain radius (Brun et al., 1992; LeFebre et al., 2003), with sw_absorption_method=1, SW penetration follows grain size in 3 spectral bands (Brun et al., 1992)
%   3 : density and cloud amount (Greuell & Konzelmann, 1994)
%   4 : exponential time decay & wetness (Bougamont & Bamber, 2005)
S.albedo_method = 1;

% Apply albedo_method method to all areas with densities below this value, or else apply direct input value from albedo_fixed, allowing albedo to be altered.
% Default value is rho water (1023 kg m-3).
S.albedo_desnity_threshold = 1023;

% apply all SW to top grid cell (0) or allow SW to penetrate surface (1)
% (default 0: if sw_absorption_method=1 and albedo_method=2, function of effective radius (Brun et al., 1992) or else dependent on snow density (taken from Bassford, 2002))
S.sw_absorption_method = 0;

% select densification model to use (default is 2):
%   1 = emperical model of Herron and Langway (1980)
%   2 = semi-emperical model of Anthern et al. (2010)
%   3 = DO NOT USE: physical model from Appendix B of Anthern et al. (2010)
%   4 = DO NOT USE: emperical model of Li and Zwally (2004)
%   5 = DO NOT USE: modified emperical model (4) by Helsen et al. (2008)
%   6 = Antarctica semi-emperical model of Ligtenberg et al. (2011)
%   7 = Greenland semi-emperical model of Kuipers Munneke et al. (2015)
S.densification_method = 2;

% select model for fresh snow accumulation density (default is 1):
%   0 = Original GEMB value, 150 kg/m^3
%   1 = Antarctica value of fresh snow density, 350 kg/m^3
%   2 = Greenland value of fresh snow density, 315 kg/m^3, Fausto et al. (2018)
%   3 = Antarctica model of Kaspers et al. (2004)
%   4 = Greenland model of Kuipers Munneke et al. (2015)
S.new_snow_method = 1;

% select method for calculating emissivity (default is 1)
%   0: direct input from emissivity parameter, no use of emissivity_re_threshold
%   1: default value of 1, in areas with grain radius below emissivity_re_threshold
%   2: default value of 1, in areas with grain radius below emissivity_re_threshold and areas of dry snow (not bare ice or wet) at the surface
S.emissivity_method = 1;

% Apply emissivity_method method to all areas with effective grain radii (re) above this value (mm), or else apply direct input value from emissivity, allowing emissivity to be altered.
% Default value is a effective grain radius of 10 mm.
S.emissivity_re_threshold = 10;

% select method for calculating thermal conductivity (default is 1)
% 1: after Sturm et al, 1997
% 2: after Calonne et al., 2011
S.thermal_conductivity_method = 1;

% GRID INITIALIZATION
% set depth of top grid cell structure (constant grid length) [m]
S.column_ztop = 10;

% set initial top vertical and min allowable grid spacings [m]
S.column_dztop = 0.05;
S.column_dzmin = S.column_dztop / 2;

% set initial/max and min model depth [m]
S.column_zmax = 250;
S.column_zmin = ceil(S.column_zmax/2 /10)*10;

% strech grid cells bellow top_z by a [top_dz * y ^ (cells bellow top_z)]
S.column_zy = 1.10;

% mean annual wind velocity [m s-1], climatology
S.V_mean = 10.0;

% optional inputs:
S.emissivity = 1.0;     % Outward longwave radiation thermal emissivity forcing at every element (default in code is 1).
% Used only if emissivity_method==0, or effective grain radius exceeds emissivity_re_threshold

S.ulw_delta = 0.0; % Delta [W/m²] with which to perturb the long wave radiation upwards. ulw_delta = 0.0 unless you have very good reason

S.is_restart = false;   % True if we want to restart from *ini parameters set in S struct

% OTHER
% specify frequency to output data (density, grid length, and temperature)
%   - 'monthly'
%   - 'daily'
S.output_frequency = 'monthly';


% FIXED ALBEDO VARIABLES
% albedo tuning parameters 
% for methods of calculating albedo see albedo function
% ----------------------------- METHODS 1 & 2 -----------------------------
S.albedo_snow  = 0.85;              % new snow albedo (0.64 - 0.89)
S.albedo_ice   = 0.48;              % range 0.27-0.58 for old snow
S.albedo_fixed = S.albedo_snow;     % Albedo forcing at every element.  Used only if albedo_method == 0, or density exceeds albedo_desnity_threshold

%Defaut values, but these can also be set as time series forcing
S.dsw_diffuse             = 0.0;    % downward diffusive shortwave radiation flux [W/m^2]
S.solar_zenith_angle      = 0.0;    % Solar Zenith Angle [degree]
S.cloud_optical_thickness = 0.0;    % Cloud Optical Thickness
S.black_carbon_snow       = 0.0;    % concentration of light absorbing carbon for snow [ppm1]
S.black_carbon_ice        = 0.0;    % concentration of light absorbing carbon for ice [ppm1]

% ------------------------------- METHOD 3 --------------------------------
% RADIATION CORRECTION FACTORS
% -> only used for met station data and Greuell & Konzelmann, 1994 albedo
S.cloud_fraction = 0.1;        % average cloud amount

% ------------------------------- METHOD 4 --------------------------------
% additonal tuning parameters albedo as a funtion of age and water content
% (Bougamont et al., 2005)
S.albedo_wet_snow_t0 = 15;     % time scale for wet snow (15-21.9) [d]
S.albedo_dry_snow_t0 = 30;     % warm snow timescale (30) [d]
S.albedo_K           = 7;      % time scale temperature coef. (7) [d]

% Constants
S.density_ice = 910;     % density of ice [kg m-3]

%% RUN GEMB

% open matlab pool for parallel processing
if S.n_workers > 1
    parpool(S.n_workers)
end

switch S.run_id
    case "test_1"
        verbose = true;
        % output directory
        S.output_dir = '../test_1';

        [daten, P, Ta, V, dlw, dsw, eAir, pAir, LP] = simulate_climate_forcing(run_id);
        S = combineStrucData_GEMB(S,LP,1);
    
        GEMB(daten, Ta, V, dlw, dsw, eAir, pAir, P, S, S.is_restart, verbose)
    otherwise
        error("input case not defined")
end