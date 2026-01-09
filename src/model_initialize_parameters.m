function options = model_initialize_parameters(options)
% model_initialize_parameters initializes and validates the model configuration
% options, setting default values for physics modules, grid geometry, and
% output controls.
%
%% Syntax
%
% options = model_initialize_parameters(options)
%
%% Description
%
% This function utilizes MATLAB's argument validation framework to parse,
% validate, and set defaults for the GEMB model parameters. It allows the user
% to configure specific physical parameterizations and numerical settings.
%
% Key Configuration Areas:
% 1. Initialization: Sets spin-up cycles and run identifiers.
% 2. Densification & Mass: Selects firn compaction models, fresh snow density, 
%    and ice density.
% 3. Energy Balance: Configures longwave emissivity, thermal conductivity, 
%    and detailed albedo/shortwave penetration schemes.
% 4. Grid Geometry: Defines the vertical discretization, including the 
%    high-resolution surface capture zone and deep layer stretching.
% 5. Output: Controls temporal resolution and array padding.
%
%% Inputs
%
%  options                          : struct       Structure containing model configuration fields (Name-Value pairs).
%
%    --- GENERAL & INITIALIZATION ---
%    .run_prefix                    : string       Unique run identifier (default: "default").
%    .n_spinup_cycles               : integer      Number of spin-up cycles (default: 0).
%
%    --- DENSITY & DENSIFICATION ---
%    .densification_method          : string       Model: "HerronLangway", "Anthern", "Ligtenberg".
%    .densification_coeffs_M01      : string       Coeffs for Ligtenberg model (e.g., "Gre_RACMO_GS_SW0").
%    .new_snow_method               : string       Fresh snow density model (e.g., "350kgm2", "Fausto").
%    .density_ice                   : double       Density of glacier ice [kg m^-3].
%
%    --- LONGWAVE EMISSIVITY ---
%    .emissivity_method             : string       Method: "uniform", "re_threshold", "re_w_threshold".
%    .emissivity                    : double       Base longwave emissivity (0-1).
%    .emissivity_re_large           : double       Emissivity for large grain sizes (0-1).
%    .emissivity_re_threshold       : double       Grain radius threshold [mm] for emissivity switching.
%
%    --- THERMAL CONDUCTIVITY ---
%    .thermal_conductivity_method   : string       Model: "Sturm" or "Calonne".
%
%    --- ALBEDO & RADIATION ---
%    .albedo_method                 : string       Scheme: "GardnerSharp", "GreuellKonzelmann", etc.
%    .albedo_desnity_threshold      : double       Density threshold [kg m^-3] below which albedo_method is applied.
%    .sw_absorption_method          : double       0 (surface only) or 1 (subsurface penetration).
%    .albedo_snow                   : double       Albedo for fresh snow (0.5-0.95).
%    .albedo_ice                    : double       Albedo for bare ice (0.2-0.6).
%    .albedo_fixed                  : double       Fixed albedo used if albedo_method="None" or density > threshold.
%    .dsw_diffuse                   : double       Downward diffusive shortwave flux [W m^-2].
%    .solar_zenith_angle            : double       Solar Zenith Angle [degrees].
%    .cloud_optical_thickness       : double       Cloud Optical Thickness.
%    .black_carbon_snow             : double       BC concentration in snow [ppm].
%    .black_carbon_ice              : double       BC concentration in ice [ppm].
%    .cloud_fraction                : double       Cloud fraction (0-1) for GreuellKonzelmann method.
%    .albedo_wet_snow_t0            : double       Time scale [d] for wet snow (BougamontBamber).
%    .albedo_dry_snow_t0            : double       Time scale [d] for dry snow (BougamontBamber).
%    .albedo_K                      : double       Temperature coef. time scale [d] (BougamontBamber).
%
%    --- OUTPUT CONTROLS ---
%    .output_frequency              : string       Output resolution: "daily", "monthly", or "all".
%    .output_padding                : integer      Extra vertical levels for grid resizing (default: 1000).
%
%    --- GRID GEOMETRY ---
%    .column_ztop                   : double       Depth of constant grid spacing at the surface [m].
%    .column_dztop                  : double       Initial surface grid spacing [m].
%    .column_dzmin                  : double       Minimum allowable grid spacing [m].
%    .column_dzmax                  : double       Maximum allowable grid spacing [m].
%    .column_zmax                   : double       Maximum total column depth [m].
%    .column_zmin                   : double       Minimum total column depth [m].
%    .column_zy                     : double       Grid stretching factor for lower layers.
%
%    * See arguments block below for full list of valid defaults and ranges *
%
%% Outputs
%
%  options                          : struct       Validated structure with all necessary model parameters populated.
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

    arguments
        %% GEMB INITIALIZATION
        % unique model run ID to save output as
        options.run_prefix (1,1) string = "default";
        
        % spin-up
        % number of cycles of met data run before output is calculated.
        % set spinUp = 0 for no spin up
        options.n_spinup_cycles (1,1) double {mustBeInteger, mustBeInRange(options.n_spinup_cycles, 0, 10000)} = 0;  
         
        %% DENSITY AND DENSIFICATION 
        % select densification model to use (default is "Anthern"):
        %   1-"HerronLangway" : empirical model of Herron and Langway (1980)
        %   2-"Anthern"       : semi-empirical model of Anthern et al. (2010)
        %   3-"AnthernB"      : DO NOT USE: physical model from Appendix B of Anthern et al. (2010)
        %   4-"LiZwally"      : DO NOT USE: empirical model of Li and Zwally (2004)
        %   5-"Helsen"        : DO NOT USE: modified empirical model (4) by Helsen et al. (2008)
        %   6-"Ligtenberg"    : semi-empirical model of Ligtenberg et al. (2011)
        options.densification_method (1,1) string {mustBeMember(options.densification_method, ...
            ["HerronLangway", "Anthern", "Ligtenberg"])} = "Anthern";
        
        % Specify densification coefficients for Ligtenberg model. These 
        % coefficients have been calibrated to match observations (default is "Gre_RACMO_GS_SW0"):
        % ------------------------ Antarctic -------------------------
        %   "Ant_ERA5_GS_SW0"    : ERA5 new albedo_method="GardnerSharp", sw_absorption_method=0
        %   "Ant_ERA5v4_Paolo23" : ERA5 v4 (Paolo et al., 2023)
        %   "Ant_ERA5_BF_SW1"    : ERA5 new albedo_method="BruneLeFebre", sw_absorption_method=1
        %   "Ant_RACMO_GS_SW0"   : RACMO callibration, default (Gardner et al., 2023)
        %   "Ant_Ligtenberg"     : Ligtenberg and others (2011), Antarctica
        % ------------------------- Greenland ------------------------
        %   "Gre_ERA5_GS_SW0"    : ERA5 new albedo_method="GardnerSharp", sw_absorption_method=0, firn & bare ice
        %   "Gre_RACMO_GS_SW0"   : RACMO calibration, default (Gardner et al., 2023)
        %   "Gre_RACMO_GB_SW1"   : ismember(albedo_method,["GreuellKonzelmann","BougamontBamber"]) && sw_absorption_method>0
        %   "Gre_KuipersMunneke" : Kuipers Munneke and others (2015) [semi-empirical], Greenland
        options.densification_coeffs_M01 (1,1) string {mustBeMember(options.densification_coeffs_M01, ...
            ["Ant_ERA5v4_Paolo23", "Ant_ERA5_BF_SW1", "Ant_RACMO_GS_SW0", "Gre_ERA5_GS_SW0" ...
            "Gre_RACMO_GS_SW0", "Gre_RACMO_GB_SW1", "Gre_KuipersMunneke"])} = "Gre_RACMO_GS_SW0";

        % select model for fresh snow accumulation density (default is "350kgm2"):
        %   0-"150kgm2"       : Original GEMB value, 150 kg/m^3
        %   1-"350kgm2"       : Antarctica value of fresh snow density, 350 kg/m^3
        %   2-"Fausto"        : Greenland value of fresh snow density, 315 kg/m^3, Fausto et al. (2018)
        %   3-"Kaspers"       : Antarctica model of Kaspers et al. (2004)
        %   4-"KuipersMunneke": Greenland model of Kuipers Munneke et al. (2015)
        options.new_snow_method (1,1) string {mustBeMember(options.new_snow_method, ...
            ["150kgm2", "350kgm2", "Fausto", "Kaspers", "KuipersMunneke"])} = "350kgm2";
        
        % density of glacier ice
        options.density_ice (1,1) double {mustBeInRange(options.density_ice, 800, 950)} = 910; % density of ice [kg m-3]

        %% LONGWAVE EMISSIVITY
        % Select method for calculating emissivity (default is "uniform")
        %   0-"uniform"       : uses "emissivity" for all snow/firn/ice surfaces
        %   1-"re_threshold"  : uses "emissivity" for re <= emissivity_re_threshold 
        %                       & "emissivity_re_large" for re > emissivity_re_threshold
        %   2:"re_w_threshold": uses "emissivity" for (re <= emissivity_re_threshold & there is no liquid water at the surface) 
        %                       & "emissivity_re_large" for (re > emissivity_re_threshold or when there is liquid water at the surface)
        options.emissivity_method (1,1) string {mustBeMember(options.emissivity_method, ...
            ["uniform", "re_threshold", "re_w_threshold"])} = "uniform";

        % Specify longwave emissivity (emissivity_re_large only used for
        % "emissivity_method" == "re_threshold" or "re_w_threshold"
        % Default updated to 0.98 to satisfy range check [0, 1]
        options.emissivity (1,1) double {mustBeInRange(options.emissivity, 0, 1)} = 0.98;
        options.emissivity_re_large (1,1) double {mustBeInRange(options.emissivity_re_large, 0, 1)} = 0.98;
        
        % Specify the effective grain radii (re) used when "emissivity_method" == "re_threshold" or "re_w_threshold"
        % Default value is a effective grain radius of 10 mm.
        options.emissivity_re_threshold (1,1) double {mustBeInRange(options.emissivity_re_threshold, 0, 100)} = 10;
        
        %% THEMAL CONDUCTIVITY
        % select method for calculating thermal conductivity (default is "Sturm")
        % 1-"Sturm"  : after Sturm et al, 1997
        % 2-"Calonne": after Calonne et al., 2011
        options.thermal_conductivity_method (1,1) string {mustBeMember(options.thermal_conductivity_method, ...
            ["Sturm", "Calonne"])} = "Sturm";
         
        %% ALBEDO

        % Select method of calculating albedo and subsurface absorption (default is "GardnerSharp")
        %   0-"None"             : direct input from albedo_fixed parameter, no use of albedo_desnity_threshold
        %   1-"GardnerSharp"     : effective grain radius (Gardner & Sharp, 2009)
        %   2-"BruneLeFebre"     : effective grain radius (Brun et al., 1992; LeFebre et al., 2003), with sw_absorption_method=1, SW penetration follows grain size in 3 spectral bands (Brun et al., 1992)
        %   3-"GreuellKonzelmann": density and cloud amount (Greuell & Konzelmann, 1994)
        %   4-"BougamontBamber"  : exponential time decay & wetness (Bougamont & Bamber, 2005)
        options.albedo_method (1,1) string {mustBeMember(options.albedo_method, ...
            ["None", "GardnerSharp", "BruneLeFebre", "GreuellKonzelmann", "BougamontBamber"])} = "GardnerSharp";
        
        % Apply albedo_method method to all areas with densities below this value, or else apply direct input value from albedo_fixed, allowing albedo to be altered.
        % Default value is rho water (1023 kg m-3).
        options.albedo_desnity_threshold (1,1) double {mustBeInRange(options.albedo_desnity_threshold, 0, 1023)} = 1023;
        
        % apply all SW to top grid cell (0) or allow SW to penetrate surface (1)
        % (default 0: if sw_absorption_method=1 and albedo_method=2, function of effective radius (Brun et al., 1992) or else dependent on snow density (taken from Bassford, 2002))
        options.sw_absorption_method (1,1) double {mustBeMember(options.sw_absorption_method, [0,1])} = 0;

        % for methods of calculating albedo see albedo function
        % --------------- "GardnerSharp" & "BruneLeFebre" ----------------------
        options.albedo_snow (1,1) double {mustBeInRange(options.albedo_snow, 0.5, .95)}   = 0.85; % new snow albedo (0.64 - 0.89)
        options.albedo_ice (1,1) double {mustBeInRange(options.albedo_ice, 0.2, .6)}      = 0.48; % range 0.27-0.58 for old snow
        options.albedo_fixed (1,1) double {mustBeInRange(options.albedo_fixed, 0.2, .95)} = 0.85; % Albedo forcing at every element.  Used only if albedo_method == 0, or density exceeds albedo_desnity_threshold
        
        % Default values, but these can also be set as time series forcing
        options.dsw_diffuse (1,1) double {mustBeInRange(options.dsw_diffuse, 0, 1000)}                       = 0.0;  % downward diffusive shortwave radiation flux [W/m^2]
        options.solar_zenith_angle (1,1) double {mustBeInRange(options.solar_zenith_angle, 0, 90)}           = 0.0;  % Solar Zenith Angle [degree]
        options.cloud_optical_thickness (1,1) double {mustBeInRange(options.cloud_optical_thickness, 0, 30)} = 0.0;  % Cloud Optical Thickness
        options.black_carbon_snow (1,1) double {mustBeInRange(options.black_carbon_snow, 0, 2)}              = 0.0;  % concentration of light absorbing carbon for snow [ppm]
        options.black_carbon_ice (1,1) double {mustBeInRange(options.black_carbon_ice, 0, 2)}                = 0.0;  % concentration of light absorbing carbon for ice [ppm]
        
        % ------------------------------- "GreuellKonzelmann" --------------------------------
        % radiation correction factor
        % -> only used for met station data and Greuell & Konzelmann, 1994 albedo
        options.cloud_fraction (1,1) double {mustBeInRange(options.cloud_fraction, 0, 1)} = 0.1; % average cloud amount
        
        % ------------------------------- "BougamontBamber" --------------------------------
        % additional tuning parameters albedo as a funtion of age and water content
        % (Bougamont et al., 2005)
        options.albedo_wet_snow_t0 (1,1) double {mustBeInRange(options.albedo_wet_snow_t0, 5, 25)}  = 15; % time scale for wet snow (15-21.9) [d]
        options.albedo_dry_snow_t0 (1,1) double {mustBeInRange(options.albedo_dry_snow_t0, 20, 40)} = 30; % warm snow timescale (30) [d]
        options.albedo_K (1,1) double {mustBeInRange(options.albedo_K, 2, 12)}                      = 7;  % time scale temperature coef. (7) [d]
       
        %% MODEL OUTPUT
        % specify frequency to output data (density, grid length, and temperature)
        %   - "monthly"
        %   - "daily"
        %   - "all"
        options.output_frequency (1,1) string {mustBeMember(options.output_frequency, ...
            ["all", "monthly", "daily"])} = 'monthly';

        % number of additional vertical levels in output initialization to accommodate changing grid size
        options.output_padding  (1,1) double {mustBeInteger, mustBeInRange(options.output_padding, 0, 10000)} = 1000;   
    
        %% GRID INITIALIZATION
        % set depth of top grid cell structure (constant grid length) [m]
        options.column_ztop (1,1) double {mustBeInRange(options.column_ztop, 0, 100)} = 10;
        
        % set initial top vertical and min allowable grid spacings [m]
        options.column_dztop (1,1) double {mustBeInRange(options.column_dztop, 0, 0.2)} = 0.05;
        options.column_dzmin (1,1) double {mustBeInRange(options.column_dzmin, 0, 0.2)} = 0.05 / 2;
        options.column_dzmax (1,1) double {mustBeInRange(options.column_dzmax, 0, 0.2)} = 0.05 + (0.05 / 2);

        % set initial/max and min model depth [m]
        options.column_zmax (1,1) double {mustBeInRange(options.column_zmax, 0, 1000)} = 250;
        options.column_zmin (1,1) double {mustBeInRange(options.column_zmin, 0, 1000)} = ceil(250/2 /10)*10;
        
        % stretch grid cells below top_z by a [top_dz * y ^ (cells below top_z)]
        options.column_zy (1,1) double {mustBeInRange(options.column_zy, 1, 2)} = 1.10;
    end
end