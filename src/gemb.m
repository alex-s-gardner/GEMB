function OutData = gemb(Profile, ClimateForcing, ModelParam, display_options)
% gemb runs the Glacier Energy and Mass Balance (GEMB) model by Gardner et al., 2023.
%
% GEMB calculates a 1-D surface glacier mass balance, includes detailed
% representation of subsurface processes, and key features include:
%
% * melt water percolation and refreeze
% * pore water retention
% * dynamic albedo with long-term memory
% * subsurface temperature diffusion
% * subsurface penetration of shortwave radiation
%
% GEMB is the primary driver function for the Glacier Energy and Mass Balance
% model. It integrates a 1-D column of snow/firn/ice forward in time given
% initial state vectors and meteorological forcing.
%
% The function performs the following steps:
% 1. Initializes outputs and calculates time steps from the forcing data.
% 2. Runs a spin-up loop (defined by ModelParam.spinup_cycles) to equate
%    the model state with the climate forcing.
% 3. Iterates through each time step in the climate forcing data.
% 4. Calls the physics engine (gemb_core) to calculate energy fluxes,
%    melt, percolation, and layer updates.
% 5. Checks for mass conservation and boundary condition stability.
% 6. Aggregates and returns the simulation results.
%
%% Syntax
%
%  OutData = gemb(Profile, ClimateForcing, ModelParam)
%  OutData = gemb(..., verbose=true)
%  OutData = gemb(..., display_waitbar=false)
%
%% Description
%
% OutData = gemb(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, ClimateForcing, ModelParam)
% produces time series of snow, firn, and ice properties OutData from input
% vectors of the initial column temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, and albedo, albedo_diffuse.
% Input ClimateForcing is a structure containing time series of surface
% forcing parameters descibed below. Input ModelParam is from the function
% model_initialize_parameters.m. 
%
%   Profile                  : Table containing initial column state variables. Must include:
%     temperature            : Vector of initial layer temperatures [K].
%     dz                     : Vector of initial layer thicknesses [m].
%     density                : Vector of initial layer densities [kg m^-3].
%     water                  : Vector of initial layer water content [kg m^-2].
%     grain_radius           : Vector of initial grain sizes (effective radius) [mm].
%     grain_dendricity       : Vector of initial grain dendricity (0-1).
%     grain_sphericity       : Vector of initial grain sphericity (0-1).
%     albedo                 : Initial surface albedo (0-1).
%     albedo_diffuse         : Initial diffuse albedo accumulator.
%   ClimateForcing           : Structure containing time-series meteorological data.
%     .dates                 : datenum      Time vector.
%     .shortwave_downward    : W m^-2       Downward shortwave radiation.
%     .longwave_downward     : W m^-2       Downward longwave radiation.
%     .temperature_air       : K            Air temperature.
%     .pressure_air          : Pa           Air pressure.
%     .relative_humidity     : %            Relative humidity.
%     .vapor_pressure        : Pa           Vapor pressure.
%     .wind_speed            : m s^-1       Wind speed.
%     .precipitation         : kg m^-2      Precipitation.
%     .wind_observation_height         : 
%     .temperature_observation_height         : 
%     .temperature_air_mean  : 
%     .wind_speed_mean       : 
%     .precipitation_mean    : 
%   ModelParam     : Structure containing model configuration parameters.
%                    Must include 'run_prefix' and 'spinup_cycles'. See
%                    model_initialize_parameters for more information.
% 
% OutData = gemb(..., verbose=true) turns on additional checks to ensure
% the model conserves mass and energy for every timestep and logs the
% results. Note: the verbose=true option may add ~10% to total processing
% time. 
%
% OutData = gemb(..., display_waitbar=false) disables the graphical
% waitbar. 
%
%% Example
% Run a basic example: 
%   
%   Initialize model parameters:
%   ModelParam = model_initialize_parameters(output_frequency="daily");
%    
%   % Generate sample data: 
%   time_step_hours = 3;
%   ClimateForcing = simulate_climate_forcing("test_1", time_step_hours);
%   
%   % Initialize grid:
%   Profile = model_initialize_column(ModelParam, ClimateForcing);
%   
%   % Run GEMB: 
%   OutData = gemb(Profile, ClimateForcing, ModelParam);
% 
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 


%% Check Inputs 

arguments 
    Profile           (:,9) table {mustContainVariables(Profile, ["temperature", "dz", "density", "water", "grain_radius", "grain_dendricity", "grain_sphericity", "albedo", "albedo_diffuse"])}
    ClimateForcing    (1,1) struct {mustHaveFields(ClimateForcing, ["dates", "shortwave_downward", "longwave_downward", "temperature_air", "pressure_air", "vapor_pressure", "wind_speed", "precipitation","wind_observation_height","temperature_observation_height","temperature_air_mean","wind_speed_mean","precipitation_mean"])}
    ModelParam                      (1,1) struct {mustHaveFields(ModelParam, ["run_prefix", "spinup_cycles","output_frequency","output_padding","black_carbon_snow","black_carbon_ice","cloud_optical_thickness","solar_zenith_angle","shortwave_downward_diffuse","cloud_fraction","density_ice"])}
    display_options.verbose         (1,1) logical = false
    display_options.display_waitbar (1,1) logical = true
end

assert(ModelParam.rain_temperature_threshold>=270.15 & ModelParam.rain_temperature_threshold<=276.15,'ModelParam.rain_temperature_threshold should be within three degrees of 273.15. Ensure you are using kelvin.')

verbose = display_options.verbose;

% Convert table input to individual variables for backward compatibility with gemb_core and subfunctions. 
temperature      = Profile.temperature;
dz               = Profile.dz;
density          = Profile.density;
water            = Profile.water;
grain_radius     = Profile.grain_radius;
grain_dendricity = Profile.grain_dendricity;
grain_sphericity = Profile.grain_sphericity;
albedo           = Profile.albedo;
albedo_diffuse   = Profile.albedo_diffuse;

assert(all(temperature>0),'Temperature Profile must be in kelvin and exceed 0 K.')
if any(temperature<100)
    warning('Temperature Profile should be in kelvin, but some values are below 100, suggesting an error. Confirm that the units are kelvin.')
end
assert(all(dz>0),'Profile variable dz must be positive.')
assert(all(water>=0),'Profile variable water must be greater than or equal to zero.')
assert(all(grain_radius>=0),'Profile variable grain_radius must be greater than or equal to zero.')
assert(all(grain_dendricity>=0) & all(grain_dendricity<=1),'Profile variable grain_dendricity must be in the range of 0 to 1.')
assert(all(grain_sphericity>=0) & all(grain_dendricity<=1),'Profile variable grain_dendricity must be in the range of 0 to 1.')
assert(all(albedo>=0) & all(albedo<=1),'Profile variable albedo must be in the range of 0 to 1.')
assert(all(albedo_diffuse>=0) & all(albedo_diffuse<=1),'Profile variable albedo_diffuse must be in the range of 0 to 1.')

%% Begin GEMB

if verbose
    disp(['------------------ STARTING RUN # ' num2str(ModelParam.run_prefix) ' --------------------' ])
    tic                                        % start timer
end
ClimateForcing.dates = datenum(ClimateForcing.dates); %  
dates = ClimateForcing.dates;              % extract dates for convenience
dt    = (dates(2)-dates(1)) * (60*60*24);  % input time step in seconds

if rem(dt,1) ~= 0
    if verbose
        warning('Rounding dt as it is not an exact integer: dt = %0.4f', dt)
    end
    dt = round(dt);
end

% initialize monolevel variables
evaporation_condensation = 0;              % surface evaporation (-) condensation (+) [kg m-2]    
melt_surface             = 0;              % initialize surface melt for albedo parameterization

% pre calculate (this is a speed optimization for thermal)
ModelParam.dt_divisors = fast_divisors(dt * 10000)/10000;

%% Initialize output structure

column_length = length(dz);
[output_index, OutData, OutCum] = model_initialize_output(column_length, ClimateForcing, ModelParam);

%% Initialize Progress Bar Variables

total_cycles = ModelParam.spinup_cycles + 1;

if display_options.display_waitbar
    steps_per_cycle   = length(dates);
    total_steps       = total_cycles * steps_per_cycle;
    waitbar_step_mod  = max(round(total_steps/100),1);
    global_step_count = 0;
    
    % Create waitbar if running in a graphical environment
    if usejava('desktop')
        h_bar = waitbar(0, 'Initializing GEMB Simulation...', 'Name', 'GEMB Progress');
    else
        h_bar = [];
    end
end

%% Start spinup loop

for simulation_iteration = 1:total_cycles

    % Initialize cumulative variables:
    runoff_cumulative                   = 0;
    refreeze_cumulative                 = 0;
    melt_cumulative                     = 0;
    evaporation_condensation_cumulative = 0;
    precipitation_cumulative            = 0;
    mass_added_cumulative               = 0;
    melt_surface_cumulative             = 0;
    rain_cumulative                     = 0;

    % Start loop for data frequency

    % Specify the time range over which the mass balance is to be calculated:
    for date_ind = 1:length(dates)

        % Extract daily data:
        ClimateForcingStep = model_inputs_single_timestep(date_ind, dt, ClimateForcing, ModelParam);
        
        % run GEMB for a single time interval
        [temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, evaporation_condensation, melt_surface, shortwave_net, heat_flux_sensible, ...
            heat_flux_latent, longwave_upward, rain, melt, runoff, refreeze, mass_added, ~, ...
            densification_from_compaction, densification_from_melt] = ...
           gemb_core(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, evaporation_condensation, melt_surface, ...
            ClimateForcingStep, ModelParam, verbose);

        % calculate net longwave [water m-2]
        longwave_net = ClimateForcingStep.longwave_downward - longwave_upward;

        % sum component mass changes [kg m-2]
        mass_added_cumulative    = mass_added + mass_added_cumulative;
        melt_cumulative          = melt + melt_cumulative;
        melt_surface_cumulative  = melt_surface + melt_surface_cumulative;
        runoff_cumulative        = runoff + runoff_cumulative;
        precipitation_cumulative = ClimateForcingStep.precipitation +  precipitation_cumulative;
        evaporation_condensation_cumulative = evaporation_condensation + evaporation_condensation_cumulative;   % evap(-) / cond(+)
        rain_cumulative          = rain + rain_cumulative;
        refreeze_cumulative      = refreeze + refreeze_cumulative;
      
        % Check if not in spinup_cycle
        if simulation_iteration == ModelParam.spinup_cycles + 1

            % grow cumulative output values
            OutCum = model_cumulative_add(melt, runoff, refreeze, evaporation_condensation, rain, mass_added, ...
                shortwave_net, longwave_net, heat_flux_sensible, heat_flux_latent, densification_from_compaction, densification_from_melt, ...
                density, albedo, dz, ModelParam, OutCum);

            if output_index(date_ind)

                [OutData, OutCum] = ...
                    model_output_populate(density, temperature, water, dz, grain_radius, grain_dendricity, grain_sphericity, ...
                    output_index, date_ind, ModelParam, OutData, OutCum);
    
            end
        end

        if display_options.display_waitbar
            % Update Progress Bar
            global_step_count = global_step_count + 1;
            if ~isempty(h_bar) && (mod(global_step_count, waitbar_step_mod) == 0 || global_step_count == total_steps)
                 % Calculate percentage
                 pct_complete = global_step_count / total_steps;
                 
                 % Create message: Date | Cycle X of Y
                 msg = sprintf('Simulating: %s | Cycle: %d / %d', ...
                     datetime(dates(date_ind),'ConvertFrom', 'datenum'), ...
                     simulation_iteration, ...
                     total_cycles);
                 
                 % Update the bar
                 waitbar(pct_complete, h_bar, msg);
            end
        end
    end

    if verbose
        % Display cycle completed and time to screen:
        disp([num2str(ModelParam.run_prefix) ': cycle: ' num2str(simulation_iteration) ' of '  ...
            num2str(ModelParam.spinup_cycles + 1) ', cpu time: ' num2str(round(toc)) ' sec,'...
            ' avg melt: ' num2str(round(melt_cumulative/(dates(end)-dates(1))*365.25)) ...
            ' kg/m2/yr']);
    end
end

if display_options.display_waitbar
    % Close the progress bar
    if ~isempty(h_bar) && ishandle(h_bar)
        close(h_bar);
    end
end

end

%% SUBFUNCTIONS
function [output_index, OutData, OutCum] = model_initialize_output(column_length, ClimateForcing, ModelParam)
    % model_initialize_output initializes the data structures used to store model 
    % results and determines the timing of output generation.
    %
    %% Syntax
    %
    % [output_index, OutData, OutCum] = model_initialize_output(column_length, ClimateForcing, ModelParam)
    %
    %% Description
    %
    % This function sets up the output environment for the GEMB model. It performs
    % the following initialization tasks:
    %
    % 1. **Time Indexing:** Determines the logical indices (`output_index`) corresponding 
    %    to the requested output frequency (e.g., 'daily', 'monthly', or 'all').
    % 2. **Data Allocation:** Pre-allocates the `OutData` structure with NaNs for both 
    %    single-level variables (scalars per time step) and multi-level profile variables 
    %    (vectors per time step) to optimize memory usage.
    % 3. **Forcing Aggregation:** Pre-calculates time-averaged air temperature and 
    %    total precipitation for each output interval.
    % 4. **Accumulator Reset:** Initializes the `OutCum` structure, which tracks 
    %    cumulative fluxes (e.g., melt, runoff, energy) between output write steps, 
    %    setting all values to zero.
    %
    %% Inputs
    %
    %  column_length       : integer   Number of vertical grid cells in the model column.
    %  ClimateForcing      : struct    Forcing data structure containing:
    %    .dates            : datenum   Vector of time steps.
    %    .temperature_air  : K         Air temperature.
    %    .precipitation    : kg m^-2.  Precipitation.
    %  ModelParam          : struct    Model parameters structure containing:
    %    .output_frequency : string    Frequency of output ('daily', 'monthly', 'all').
    %    .output_padding   : integer   Extra buffer size for profile arrays.
    %
    %% Outputs
    %
    %  output_index   : logical      Boolean mask indicating which time steps are saved.
    %  OutData        : struct       Structure for storing time-series outputs, initialized with NaNs.
    %                                Contains fields for monolevel (e.g., 'melt', 'runoff') and profile (e.g., 'temperature', 'density') variables.
    %  OutCum         : struct       Structure for tracking cumulative values between outputs.
    %                                Initialized to zero for fields like 'mass_added', 'shortwave_net'.
    %
    %% Author Information
    % The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
    % from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
    % at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
    % 
    % Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
    % a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
    % https://doi.org/10.5194/gmd-16-2277-2023, 2023. 
    
    % Determine save time steps:
    date_vector = datevec([ClimateForcing.dates; (ClimateForcing.dates(end) ...
        + ClimateForcing.dates(end)-ClimateForcing.dates(end-1))]);
    switch ModelParam.output_frequency
        case "monthly"
            output_index = (date_vector(1:end-1,2) - date_vector(2:end,2)) ~= 0;
        case "daily"
            output_index = (date_vector(1:end-1,3) - date_vector(2:end,3)) ~= 0;
        case "all"
            output_index = true(numel(ClimateForcing.dates),1);
        otherwise
            error('ModelParam.output_frequency can only be "monthly", "daily", or "all".')
    end
    
    % single level time series
    varname.monolevel = {'dates', 'melt', 'runoff', 'refreeze', 'evaporation_condensation', 'shortwave_net', ...
        'longwave_net', 'heat_flux_sensible', 'heat_flux_latent', 'albedo_surface', 'valid_profile_length', ...
        'densification_from_compaction', 'densification_from_melt', 'ps'};
    
    n = sum(output_index);
    
    % Set single level time series to NaNs:
    for v = 1:length(varname.monolevel)
        OutData.(varname.monolevel{v}) = nan(1,n);
    end
    
    OutData.dates = ClimateForcing.dates(output_index)';

    % Time averages/totals:
    I = find(output_index);                      % save index
    for i = 1:n
        if i == 1
            OutData.temperature_air(i) = mean(ClimateForcing.temperature_air(1:I(i))) - 273.15;  % convert from K to deg C
            OutData.precipitation(i)   = sum(ClimateForcing.precipitation(1:I(i)));
        else
            OutData.temperature_air(i) = mean(ClimateForcing.temperature_air((I(i-1)+1):I(i))) - 273.15;
            OutData.precipitation(i)   = sum(ClimateForcing.precipitation((I(i-1)+1):I(i)));
        end
    end
    
    % Multi level time series:
    varname.profile = {'temperature', 'dz', 'density', 'water', 'grain_radius', 'grain_dendricity', 'grain_sphericity', 'albedo', 'albedo_diffuse', 'ps'};
    
    for v = 1:length(varname.profile)
        OutData.(varname.profile{v}) = nan(column_length + ModelParam.output_padding, n);
    end
    
    % Initialize cumulative output values:
    varname.cumulative = {'runoff', 'melt', 'refreeze', 'evaporation_condensation', 'rain', 'mass_added', 'shortwave_net', ...
        'longwave_net', 'heat_flux_sensible', 'heat_flux_latent', 'albedo_surface', 'densification_from_compaction', ...
        'densification_from_melt', 'firn_air_content'};
    
    % Set cumulative values zero:
    for v = 1:length(varname.cumulative)
        OutCum.(varname.cumulative{v}) = 0;
    end
    
    OutCum.count = 0;
end

function ClimateForcingStep = model_inputs_single_timestep(index, dt, ClimateForcing, ModelParam)
    % model_inputs_single_timestep
    %
    %% Syntax
    %
    %
    %
    %% Description
    %
    %
    %
    %% Example
    % 
    %
    %
    %% Author Information
    % The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
    % from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
    % at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
    % 
    % Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
    % a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
    % https://doi.org/10.5194/gmd-16-2277-2023, 2023. 
    
    ClimateForcingStep.dt                 = dt;                                       % time step in seconds
    ClimateForcingStep.temperature_air    = ClimateForcing.temperature_air(index);    % screen level air temperature [K]   
    ClimateForcingStep.wind_speed         = ClimateForcing.wind_speed(index);         % wind speed [m s-1]
    ClimateForcingStep.longwave_downward  = ClimateForcing.longwave_downward(index);  % downward longwave radiation flux [water m-2]
    ClimateForcingStep.shortwave_downward = ClimateForcing.shortwave_downward(index); % downward shortwave radiation flux [water m-2]
    ClimateForcingStep.vapor_pressure     = ClimateForcing.vapor_pressure(index);     % screen level vapor pressure [Pa]
    ClimateForcingStep.pressure_air       = ClimateForcing.pressure_air(index);       % screen level air pressure [Pa]
    ClimateForcingStep.precipitation      = ClimateForcing.precipitation(index);      % precipitation [kg m-2]
    
    % Location specifc parameters
    ClimateForcingStep.wind_observation_height        = ClimateForcing.wind_observation_height;
    ClimateForcingStep.temperature_observation_height = ClimateForcing.temperature_observation_height;
    ClimateForcingStep.temperature_air_mean           = ClimateForcing.temperature_air_mean;
    ClimateForcingStep.wind_speed_mean                = ClimateForcing.wind_speed_mean;
    ClimateForcingStep.precipitation_mean             = ClimateForcing.precipitation_mean;

    % if we are provided with cc and cot values, extract for the timestep
    if numel(ModelParam.black_carbon_snow)>1
        ClimateForcingStep.black_carbon_snow = ModelParam.black_carbon_snow(index);
    else
        ClimateForcingStep.black_carbon_snow = ModelParam.black_carbon_snow;
    end
    
    if numel(ModelParam.black_carbon_ice)>1
        ClimateForcingStep.black_carbon_ice = ModelParam.black_carbon_ice(index);
    else
        ClimateForcingStep.black_carbon_ice = ModelParam.black_carbon_ice;
    end
    
    if numel(ModelParam.cloud_optical_thickness)>1
        ClimateForcingStep.cloud_optical_thickness = ModelParam.cloud_optical_thickness(index);
    else
        ClimateForcingStep.cloud_optical_thickness = ModelParam.cloud_optical_thickness;
    end
    
    if numel(ModelParam.solar_zenith_angle)>1
        ClimateForcingStep.solar_zenith_angle = ModelParam.solar_zenith_angle(index);
    else
        ClimateForcingStep.solar_zenith_angle = ModelParam.solar_zenith_angle;
    end
    
    if numel(ModelParam.shortwave_downward_diffuse)>1
        ClimateForcingStep.shortwave_downward_diffuse = ModelParam.shortwave_downward_diffuse(index);
    else
        ClimateForcingStep.shortwave_downward_diffuse = ModelParam.shortwave_downward_diffuse;
    end
    
    if numel(ModelParam.cloud_fraction)>1
        ClimateForcingStep.cloud_fraction = ModelParam.cloud_fraction(index);
    else
        ClimateForcingStep.cloud_fraction = ModelParam.cloud_fraction;
    end
end

function OutCum = model_cumulative_add(melt, runoff, refreeze, evaporation_condensation, rain, mass_added, ...
    shortwave_net, longwave_net, heat_flux_sensible, heat_flux_latent, densification_from_compaction, densification_from_melt, ...
    density, albedo, dz, ModelParam, OutCum)
    % model_cumulative_add updates cumulative variables for model output.
    %
    %% Syntax
    %
    % OutCum = model_cumulative_add(OutCum, melt, runoff, refreeze, evaporation_condensation, rain, mass_added, ...
    %    shortwave_net, longwave_net, heat_flux_sensible, heat_flux_latent, longwave_upward, densification_from_compaction, densification_from_melt, ...
    %    density, albedo, dz, ModelParam)
    %
    %% Description
    %
    % This function updates the tracking structure `OutCum` by adding the current
    % time-step's values to the running totals. It replaces the use of `eval()`
    % in earlier versions of GEMB to significantly improve runtime performance.
    %
    % It performs two main tasks:
    % 1. Calculates derived variables for the current state (e.g., firn_air_content, 
    %    surface properties albedo_surface).
    % 2. Explicitly sums these values into the `OutCum` structure fields.
    %
    %% Inputs
    %
    %  OutCum                   : struct       Structure containing cumulative variables from previous steps.
    %  melt                     : kg m^-2      Melt mass.
    %  runoff                   : kg m^-2      Runoff mass.
    %  refreeze                 : kg m^-2      Refrozen mass.
    %  evaporation_condensation : kg m^-2      Evaporation/Condensation mass.
    %  rain                     : kg m^-2      Rain mass.
    %  mass_added               : kg m^-2      Mass added/removed by manage_layers.
    %  shortwave_net            : W m^-2       Net shortwave radiation.
    %  longwave_net             : W m^-2       Net longwave radiation.
    %  heat_flux_sensible       : W m^-2       Sensible heat flux.
    %  heat_flux_latent         : W m^-2       Latent heat flux.
    %  longwave_upward          : W m^-2       Upward longwave radiation.
    %  densification_from_compaction  : m      Compaction due to densification.
    %  densification_from_melt  : m            Compaction due to melt.
    %  density                  : kg m^-3      Density profile.
    %  albedo                   : fraction     Albedo profile.
    %  dz                       : m            Layer thickness profile.
    %  ModelParam               : struct       Model parameters (needs .density_ice).
    %
    %% Outputs
    %
    %  OutCum           : struct       Updated cumulative structure with incremented .count.
    %
    %% Author Information
    % The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
    % from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
    % at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
    % 
    % Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
    % a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
    % https://doi.org/10.5194/gmd-16-2277-2023, 2023. 
    
    % 1. Calculate derived variables for output
    albedo_surface       = albedo(1);
    
    % Firn Air Content (firn_air_content) [m]
    % Defined as the integrated column thickness of air equivalent.
    % firn_air_content = sum(dz * (rho_ice - rho) / rho_ice) for rho < rho_ice
    % Note: The original implementation divided by 1000 instead of ModelParam.density_ice?
    % Preserving original logic: sum(dz.*(ModelParam.density_ice - min(density,ModelParam.density_ice)))/1000;
    firn_air_content = sum(dz .* (ModelParam.density_ice - min(density, ModelParam.density_ice))) / 1000;
    
    % 2. Explicitly accumulate values
    % Using explicit assignment is significantly faster than dynamic field access
    OutCum.runoff                   = OutCum.runoff + runoff;
    OutCum.melt                     = OutCum.melt + melt;
    OutCum.refreeze                 = OutCum.refreeze + refreeze;
    OutCum.evaporation_condensation = OutCum.evaporation_condensation + evaporation_condensation;
    OutCum.rain                     = OutCum.rain + rain;
    OutCum.mass_added               = OutCum.mass_added + mass_added;
    OutCum.shortwave_net            = OutCum.shortwave_net + shortwave_net;
    OutCum.longwave_net             = OutCum.longwave_net + longwave_net;
    OutCum.heat_flux_sensible       = OutCum.heat_flux_sensible + heat_flux_sensible;
    OutCum.heat_flux_latent         = OutCum.heat_flux_latent + heat_flux_latent;
    OutCum.albedo_surface           = OutCum.albedo_surface + albedo_surface;
    
    OutCum.densification_from_compaction = OutCum.densification_from_compaction + densification_from_compaction;
    OutCum.densification_from_melt  = OutCum.densification_from_melt + densification_from_melt;
    OutCum.firn_air_content         = OutCum.firn_air_content + firn_air_content;
    
    % Increment the counter
    OutCum.count = OutCum.count + 1;

end


function [OutData, OutCum] = ...
    model_output_populate(density, temperature, water, dz, grain_radius, grain_dendricity, grain_sphericity, ...
     output_index, date_ind, ModelParam, OutData, OutCum)
    % model_output_populate stores model state and fluxes into the output structure.
    %
    %% Syntax
    %
    % [OutData, OutCum] = model_output_populate(OutData, OutCum, ...
    %    density, temperature, water, dz, grain_radius, grain_dendricity, grain_sphericity, ...
    %    ModelParam, output_index, date_ind)
    %
    %% Description
    %
    % This function checks if the current simulation step (`date_ind`) corresponds
    % to a requested output interval. If so, it performs the following:
    % 1. **Time-Averaging:** Calculates averages for state variables (e.g., Albedo, Temperature)
    %    stored in `OutCum` by dividing by the step count.
    % 2. **Profile Storage:** Saves instantaneous vertical profiles (Density, Temperature, 
    %    Grain properties) into the pre-allocated `OutData` arrays.
    % 3. **Reset:** Resets all cumulative trackers in `OutCum` to zero for the next interval.
    %
    %% Inputs
    %
    %  OutData          : struct       Main output structure containing time series.
    %  OutCum           : struct       Accumulator structure for the current interval.
    %  density, temperature, water, ... : vectors      Current vertical profiles of density, temperature, etc.
    %  ModelParam       : struct       Model parameters (needs .density_ice, .output_padding).
    %  output_index     : logical      Vector indicating which time steps are output steps.
    %  date_ind         : integer      Current time step index.
    %
    %% Outputs
    %
    %  OutData          : struct       Updated output structure.
    %  OutCum           : struct       Reset accumulator structure (if output occurred).
    %
    %% Author Information
    % The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
    % from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
    % at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
    % 
    % Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
    % a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
    % https://doi.org/10.5194/gmd-16-2277-2023, 2023. 
    
    if output_index(date_ind)
        
        % Store model output in structure format
        varname = fieldnames(OutCum);

        % Remove "count" because users don't need it: 
        varname(strcmpi(varname,'count')) = [];
        
        % Determine the index for the output array (where to save in OutData)
        r = sum(output_index(1:date_ind));
        
        % Define which variables are cumulative totals vs time-averages
        cumulative_vars = {'melt', 'runoff', 'refreeze', 'evaporation_condensation', 'precipitation', 'rain', 'mass_added', ...
                           'densification_from_compaction', 'densification_from_melt'};
        is_cumulative = ismember(varname, cumulative_vars);
        
        % Transfer data from OutCum to OutData
        for v = 1:length(varname)
            if is_cumulative(v)
                % Cumulative variables (e.g. Melt, Runoff) are just copied
                OutData.(varname{v})(r) = OutCum.(varname{v});
            else
                % Averaged variables (e.g. Temperature) are divided by step count
                OutData.(varname{v})(r) = OutCum.(varname{v}) / OutCum.count;
            end
        end
        
        % Calculate Total Mass for surface height change (ps) calculation
        mass_total = sum(dz .* density);
    
        % Store instantaneous level (profile) data
        o = (size(density,1) - 1);
        
        % Ensure the output array is large enough to hold the current column depth
        if (length(OutData.dz) - o) < 1
            error("the length of the simulation column [%0.0f] is larger than the lenght of the output array [%0.0f]\n    -> try increasing the value of ModelParam.output_padding", (o+1), size(OutData.grain_radius,1))
        end
        
        % Save vertical profiles
        OutData.grain_radius(end-o:end,r)     = grain_radius;
        OutData.density(end-o:end,r)          = density;
        OutData.temperature(end-o:end,r)      = temperature;
        OutData.water(end-o:end,r)            = water;
        OutData.dz(end-o:end,r)               = dz;
        OutData.grain_dendricity(end-o:end,r) = grain_dendricity;
        OutData.grain_sphericity(end-o:end,r) = grain_sphericity;
        
        % Calculate surface height change relative to ice equivalent
        OutData.ps(end-o:end,r)               = sum(dz) - mass_total/ModelParam.density_ice;
        OutData.valid_profile_length(r)       = o+1;
        
        % Reset cumulative values back to zero for the next interval
        for v = 1:length(varname)
            OutCum.(varname{v}) = 0;
        end
        
        OutCum.count = 0;
    end
end

% Custom Validation Function
function mustHaveFields(s, requiredFields)
    for f = requiredFields
        if ~isfield(s, f)
            error('InvalidInput:MissingField', ...
                'Structure is missing the required field: %s', f);
        end
    end
end

% Custom validation function
function mustContainVariables(tab, requiredVars)
    if ~all(ismember(requiredVars, tab.Properties.VariableNames))
        eid = 'Table:MissingVariables';
        msg = 'Input table must contain variables: ' + strjoin(requiredVars, ", ");
        throwAsCaller(MException(eid, msg));
    end
end