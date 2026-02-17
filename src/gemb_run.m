%==========================================================================
% GEMB: Glacier Energy and Mass Balance Model (the B in GEMB is silent)
%==========================================================================
%
%  Version: 2.0
%  Author : Alex S. Gardner, NASA Jet Propulsion Laboratory
%  Author : Nicole Schlegel, NOAA
%  Author : Chad A. Greene , NASA Jet Propulsion Laboratory
%  Date   : 2026
%
%--------------------------------------------------------------------------
% DESCRIPTION
%--------------------------------------------------------------------------
% The Glacier Energy and Mass Balance (GEMB) model is a comprehensive 
% one-dimensional physical model simulating the surface energy balance and 
% vertical firn evolution of glaciers and ice sheets.
%
% It couples atmospheric forcing with subsurface thermodynamics and 
% densification physics to resolve the evolution of temperature, density, 
% water content, and grain properties over time.
%
%--------------------------------------------------------------------------
% KEY CAPABILITIES
%--------------------------------------------------------------------------
% 1. Surface Energy Balance (SEB):
%    - Resolves radiative fluxes (shortwave/longwave) and turbulent heat
%      fluxes (sensible/latent) using Monin-Obukhov similarity theory.
%    - Calculates surface temperature and melt generation iteratively.
%
% 2. Subsurface Physics:
%    - Solves the heat equation with phase change and meltwater percolation.
%    - Simulates densification using empirical or semi-empirical schemes.
%    - Tracks meltwater retention (irreducible water content) and refreezing.
%
% 3. Grid Management:
%    - Utilizes a Lagrangian-style vertical grid that evolves with accumulation.
%    - Automatically merges thin layers and splits thick layers to maintain
%      numerical stability.
%
%--------------------------------------------------------------------------
% MODEL PARAMETERS (Configuration)
%--------------------------------------------------------------------------
% The model is highly configurable via the model_initialize_parameters
% function. Key configuration categories include:
%
% 1. Densification Physics:
%    - Options include "HerronLangway" (empirical), "Arthern" (semi-empirical), 
%      and "Ligtenberg" (semi-empirical).
%    - Calibration coefficients for the Ligtenberg model can be specified for
%      Antarctica or Greenland (e.g., "Gre_RACMO_GS_SW0").
%      
% 2. Grid Geometry:
%    - Vertical discretization is controlled by column_ztop (depth of high-res
%      surface zone) and column_zy (stretching factor for deep layers).
%    - User defines minimum (column_dzmin) and maximum (column_dzmax) layer
%      thicknesses to ensure stability.
%
% 3. Energy Balance & Optical Properties:
%    - Albedo schemes: "GardnerSharp", "GreuellKonzelmann", "BruneLeFebre",
%      or "Bougamont2005".
%    - Solar Penetration: shortwave_absorption_method toggles between surface-only (0)
%      or subsurface extinction (1).
%    - Thermal Conductivity: Options for "Sturm" or "Calonne" parameterizations.
%    - Emissivity: Configurable methods (0, 1, 2) and thresholds based on 
%      grain size.
%
% 4. Initialization:
%    - Fresh snow density can be set via models like "350kgm2", "Fausto", 
%      or "Kaspers".
%    - Spin-up cycles (spinup_cycles) allow the model to reach equilibrium
%      before saving output.
%
%--------------------------------------------------------------------------
% CLIMATE FORCING (Inputs)
%--------------------------------------------------------------------------
% GEMB requires high-frequency meteorological forcing data. The 
% simulate_climate_forcing function generates synthetic data containing:
%
% 1. Radiation:
%    - Downward Shortwave (shortwave_downward): Calculated based on solar geometry and 
%      location.
%    - Downward Longwave (longwave_downward): Derived from air temperature and vapor 
%      pressure, with adjustments for cloud cover.
%      
% 2. Thermodynamics:
%    - Air Temperature (temperature_air): Simulated time series.
%    - Air Pressure (pressure_air): Calculated from elevation and temperature.
%    - Humidity: Relative humidity (relative_humidity) and vapor pressure (vapor_pressure).
%
% 3. Dynamics & Mass:
%    - Wind Speed (wind_speed): Stochastic generation with seasonal noise.
%    - Precipitation (precipitation): Accumulated mass input.
%
% 4. Metadata:
%    - Includes location (latitude, longitude, elevation) and measurement heights (wind_observation_height, temperature_observation_height).
%
%--------------------------------------------------------------------------
% WORKFLOW EXECUTION
%--------------------------------------------------------------------------
% A typical simulation follows this sequence (see gemb_run.m):
%
%   1. CONFIGURATION : Define physics options (ModelParam).
%                      -> model_initialize_parameters()
%
%   2. FORCING       : Load/Generate meteorological drivers (ClimateForcing).
%                      -> simulate_climate_forcing()
%
%   3. INITIALIZATION: Establish initial column state (temperature, density, water, etc.).
%                      -> model_initialize_column()
%
%   4. SIMULATION    : Execute the core time-stepping loop.
%                      -> gemb()
%
%--------------------------------------------------------------------------
% REFERENCES
%--------------------------------------------------------------------------
% If you use GEMB in your research, please cite:
%
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass
% Balance (GEMB): a model of firn processes for cryosphere research, Geosci.
% Model Dev., 16, 2277–2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023.
%
% Code Repository: https://github.com/alex-s-gardner/GEMB
%==========================================================================

% Specify run case
run_id = "test_1";

%% RUN GEMB
switch run_id
    case "test_1"
        verbose = true;

        % [1] specify model parameters
        % [if data is not modified then it can be passed as an stucture]
        ModelParam = model_initialize_parameters;

        time_step_hours = 3;

        ClimateForcing = simulate_climate_forcing(run_id, time_step_hours);

        % [3] Initialize grid -or- load in data to restart a simulation
        [temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse] = ...
            model_initialize_column(ModelParam, ClimateForcing);
        
        % [4] Rum GEMB
        OutData = ...
            gemb(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, ClimateForcing, ModelParam, verbose);
        
        % [5] Save model output and model run settings

        %save(fullfile(S.output_dir, S.run_id), '-struct', 'O', '-v7.3')
        %save(fullfile(S.output_dir, S.run_id), 'S', '-append')

    otherwise
        error("input case not defined")
end