function OutData = gemb(T, dz, d, W, re, gdn, gsp, a, a_diffuse, ClimateForcing, ModelParam, verbose)
% GEMB runs the Glacier Energy and Mass Balance (GEMB) model by Gardner et al., 2023.
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
%% Syntax
%
% OutData = gemb(T, dz, d, W, re, gdn, gsp, a, a_diffuse, ClimateForcing, ModelParam, verbose)
%
%% Description
%
% GEMB is the primary driver function for the Glacier Energy and Mass Balance
% model. It integrates a 1-D column of snow/firn/ice forward in time given
% initial state vectors and meteorological forcing.
%
% The function performs the following steps:
% 1. Initializes outputs and calculates time steps from the forcing data.
% 2. Runs a spin-up loop (defined by ModelParam.n_spinup_cycles) to equate
%    the model state with the climate forcing.
% 3. Iterates through each time step in the climate forcing data.
% 4. Calls the physics engine (gemb_core) to calculate energy fluxes,
%    melt, percolation, and layer updates.
% 5. Checks for mass conservation and boundary condition stability.
% 6. Aggregates and returns the simulation results.
%
%% Inputs
%
% T              : Vector of initial layer temperatures [K].
% dz             : Vector of initial layer thicknesses [m].
% d              : Vector of initial layer densities [kg m^-3].
% W              : Vector of initial layer water content [kg m^-2].
% re             : Vector of initial grain sizes (effective radius) [m].
% gdn            : Vector of initial grain dendricity (0-1).
% gsp            : Vector of initial grain sphericity (0-1).
% a              : Initial surface albedo (0-1).
% a_diffuse      : Initial diffuse albedo accumulator.
% ClimateForcing : Structure containing time-series meteorological data.
%                  Must include 'daten' (dates), 'dlw' (downward longwave),
%                  'P' (precipitation), and other standard forcing fields.
% ModelParam     : Structure containing model configuration parameters.
%                  Must include 'run_prefix' and 'n_spinup_cycles'.
% verbose        : Logical flag or integer to control console output detail.
%
%% Outputs
%
% OutData        : Structure containing the model results. The specific
%                  fields are determined by 'model_initialize_output' and
%                  populated by 'model_output_populate'. Typically includes
%                  time series of mass balance, surface fluxes, and
%                  subsurface profiles.
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 


disp(['------------------ STARTING RUN # ' num2str(ModelParam.run_prefix) ' --------------------' ])
tic                                        % start timer
daten = ClimateForcing.daten;              % extract daten for convenience
dt    = (daten(2)-daten(1)) * (60*60*24);  % input time step in seconds

if rem(dt,1) ~= 0
    warning('rounding dt as it is not an exact integer: dt = %0.4f', dt)
    dt = round(dt);
end

% initialize monolevel variables
EC     = 0;                        % surface evaporation (-) condensation (+) [kg m-2]    
M_surf = 0;                        % initialize surface melt for albedo parameterization

% pre calculate (this is a speed optimization for thermal)
ModelParam.dt_divisors = fast_divisors(dt * 10000)/10000;

%% initialize output structure
column_length = length(dz);
[output_index, OutData, OutCum] = model_initialize_output(column_length, ClimateForcing, ModelParam);

%% Initialize Progress Bar Variables
total_cycles = ModelParam.n_spinup_cycles + 1;
steps_per_cycle = length(daten);
total_steps = total_cycles * steps_per_cycle;
waitbar_step_mod = max(round(total_steps/100),1);
global_step_count = 0;

% Create waitbar if running in a graphical environment
if usejava('desktop')
    h_bar = waitbar(0, 'Initializing GEMB Simulation...', 'Name', 'GEMB Progress');
else
    h_bar = [];
end


%% Start spinup loop
for simulation_iteration = 1:total_cycles

    % Initialize cumulative variables:
    R_cumulative          = 0;
    F_cumulative          = 0;
    M_cumulative          = 0;
    EC_cumulative         = 0;
    P_cumulative          = 0;
    M_added_cumulative    = 0;
    M_surf_cumulative     = 0;
    Ra_cumulative         = 0;

    %% Start loop for data frequency

    % Specify the time range over which the mass balance is to be calculated:
    for date_ind = 1:length(daten)

        % Extract daily data:
        [ClimateForcingStep] = model_inputs_single_timestep(date_ind, dt, ClimateForcing, ModelParam);
        
        % run GEMB for a single time interval
        [T, dz, d, W, re, gdn, gsp, a, a_diffuse, EC, M_surf, sw_net, shf, ...
            lhf, ulw, Ra, M, R, F, M_added, ~, ...
            compaction_dens, compaction_melt] = ...
           gemb_core(T, dz, d, W, re, gdn, gsp, a, a_diffuse, EC, M_surf, ...
            ClimateForcingStep, ModelParam, verbose);

        % calculate net longwave [W m-2]
        lw_net = ClimateForcingStep.dlw - ulw;

        % sum component mass changes [kg m-2]
        M_added_cumulative    = M_added + M_added_cumulative;
        M_cumulative          = M + M_cumulative;
        M_surf_cumulative     = M_surf + M_surf_cumulative;
        R_cumulative          = R + R_cumulative;
        P_cumulative          = ClimateForcingStep.P +  P_cumulative;
        EC_cumulative         = EC + EC_cumulative;   % evap(-) / cond(+)
        Ra_cumulative         = Ra + Ra_cumulative;
        F_cumulative          = F + F_cumulative;
      
        % Check if not in spinup_cycle
        if simulation_iteration == ModelParam.n_spinup_cycles + 1

            % grow cumulative output values
            OutCum = model_cumulative_add(M, R, F, EC, Ra, M_added, ...
                sw_net, lw_net, shf, lhf, ulw, compaction_dens, compaction_melt, ...
                d, a, re, dz, ModelParam, OutCum);

            if output_index(date_ind)

                [OutData, OutCum] = ...
                    model_output_populate(d, T, W, dz, re, gdn, gsp, ...
                    output_index, date_ind, ModelParam, OutData, OutCum);
    
            end
        end

        % Update Progress Bar
        global_step_count = global_step_count + 1;
        if ~isempty(h_bar) && (mod(global_step_count, waitbar_step_mod) == 0 || global_step_count == total_steps)
             % Calculate percentage
             pct_complete = global_step_count / total_steps;
             
             % Create message: Date | Cycle X of Y
             msg = sprintf('Simulating: %s | Cycle: %d / %d', ...
                 datetime(daten(date_ind),'ConvertFrom', 'datenum'), ...
                 simulation_iteration, ...
                 total_cycles);
             
             % Update the bar
             waitbar(pct_complete, h_bar, msg);
        end
    end

    % display cycle completed and time to screen
    disp([num2str(ModelParam.run_prefix) ': cycle: ' num2str(simulation_iteration) ' of '  ...
        num2str(ModelParam.n_spinup_cycles + 1) ', cpu time: ' num2str(round(toc)) ' sec,'...
        ' avg melt: ' num2str(round(M_cumulative/(daten(end)-daten(1))*365.25)) ...
        ' kg/m2/yr']);
end

% Close the progress bar
if ~isempty(h_bar) && ishandle(h_bar)
    close(h_bar);
end

end

%% SUBFUNCTIONS
function OutCum = model_cumulative_add(M, R, F, EC, Ra, M_added, ...
    sw_net, lw_net, shf, lhf, compaction_dens, compaction_melt, ...
    d, a, re, dz, ModelParam, OutCum)
% model_cumulative_add updates cumulative variables for model output.
%
%% Syntax
%
% OutCum = model_cumulative_add(OutCum, M, R, F, EC, Ra, M_added, ...
%    sw_net, lw_net, shf, lhf, ulw, compaction_dens, compaction_melt, ...
%    d, a, re, dz, ModelParam)
%
%% Description
%
% This function updates the tracking structure `OutCum` by adding the current
% time-step's values to the running totals. It replaces the use of `eval()`
% in earlier versions of GEMB to significantly improve runtime performance.
%
% It performs two main tasks:
% 1. Calculates derived variables for the current state (e.g., FAC, 
%    surface properties d1, a1, re1).
% 2. Explicitly sums these values into the `OutCum` structure fields.
%
%% Inputs
%
%  OutCum           : struct       Structure containing cumulative variables from previous steps.
%  M                : kg m^-2      Melt mass.
%  R                : kg m^-2      Runoff mass.
%  F                : kg m^-2      Refrozen mass.
%  EC               : kg m^-2      Evaporation/Condensation mass.
%  Ra               : kg m^-2      Rain mass.
%  M_added          : kg m^-2      Mass added/removed by layer management.
%  sw_net           : W m^-2       Net shortwave radiation.
%  lw_net           : W m^-2       Net longwave radiation.
%  shf              : W m^-2       Sensible heat flux.
%  lhf              : W m^-2       Latent heat flux.
%  ulw              : W m^-2       Upward longwave radiation.
%  compaction_dens  : m            Compaction due to densification.
%  compaction_melt  : m            Compaction due to melt.
%  d                : kg m^-3      Density profile.
%  a                : fraction     Albedo profile.
%  re               : mm           Grain radius profile.
%  dz               : m            Layer thickness profile.
%  ModelParam       : struct       Model parameters (needs .density_ice).
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
d1    = d(1);
a1    = a(1);
re1   = re(1);

% Firn Air Content (FAC) [m]
% Defined as the integrated column thickness of air equivalent.
% FAC = sum(dz * (rho_ice - rho) / rho_ice) for rho < rho_ice
% Note: The original implementation divided by 1000 instead of ModelParam.density_ice?
% Preserving original logic: sum(dz.*(ModelParam.density_ice - min(d,ModelParam.density_ice)))/1000;
FAC = sum(dz .* (ModelParam.density_ice - min(d, ModelParam.density_ice))) / 1000;

% 2. Explicitly accumulate values
% Using explicit assignment is significantly faster than dynamic field access
OutCum.R               = OutCum.R + R;
OutCum.M               = OutCum.M + M;
OutCum.F               = OutCum.F + F;
OutCum.EC              = OutCum.EC + EC;
OutCum.Ra              = OutCum.Ra + Ra;
OutCum.M_added         = OutCum.M_added + M_added;

OutCum.sw_net          = OutCum.sw_net + sw_net;
OutCum.lw_net          = OutCum.lw_net + lw_net;
OutCum.shf             = OutCum.shf + shf;
OutCum.lhf             = OutCum.lhf + lhf;

OutCum.a1              = OutCum.a1 + a1;
OutCum.re1             = OutCum.re1 + re1;
OutCum.d1              = OutCum.d1 + d1;

OutCum.compaction_dens = OutCum.compaction_dens + compaction_dens;
OutCum.compaction_melt = OutCum.compaction_melt + compaction_melt;
OutCum.FAC             = OutCum.FAC + FAC;

% Increment the counter
OutCum.count = OutCum.count + 1;

end


function [OutData, OutCum] = ...
    model_output_populate(d, T, water, dz, re, gdn, gsp, ...
     output_index, date_ind, ModelParam, OutData, OutCum)
% model_output_populate stores model state and fluxes into the output structure.
%
%% Syntax
%
% [OutData, OutCum] = model_output_populate(OutData, OutCum, ...
%    d, T, water, dz, re, gdn, gsp, ...
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
%  d, T, water, ... : vectors      Current vertical profiles of density, temperature, etc.
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
    
    % Determine the index for the output array (where to save in OutData)
    r = sum(output_index(1:date_ind));
    
    % Define which variables are cumulative totals vs time-averages
    cumulative_vars = {'M', 'R', 'F', 'EC', 'P', 'Ra', 'M_added', ...
                       'compaction_dens', 'compaction_melt'};
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
    M_total = sum(dz .* d);

    % Store instantaneous level (profile) data
    o = (size(d,1) - 1);
    
    % Ensure the output array is large enough to hold the current column depth
    if (length(OutData.dz) - o) < 1
        error("the length of the simulation column [%0.0f] is larger than the lenght of the output array [%0.0f]\n    -> try increasing the value of ModelParam.output_padding", (o+1), size(OutData.re,1))
    end
    
    % Save vertical profiles
    OutData.re(end-o:end,r)      = re;
    OutData.d(end-o:end,r)       = d;
    OutData.T(end-o:end,r)       = T;
    OutData.water(end-o:end,r)   = water;
    OutData.dz(end-o:end,r)      = dz;
    OutData.gdn(end-o:end,r)     = gdn;
    OutData.gsp(end-o:end,r)     = gsp;
    
    % Calculate surface height change relative to ice equivalent
    OutData.ps(end-o:end,r)  = sum(dz) - M_total/ModelParam.density_ice;
    OutData.m(r)             = o+1;
    
    % Reset cumulative values back to zero for the next interval
    for v = 1:length(varname)
        OutCum.(varname{v}) = 0;
    end
    
    OutCum.count = 0;
end

end