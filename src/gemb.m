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

% fixed lower temperature bounday condition - T is fixed
T_bottom = T(end);

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

    % Determine initial mass [kg]:
    M_initial = sum (dz .* d) + sum(W);

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
            lhf, ulw, Ra, M, R, F, M_added, E_added, ...
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
        W_total               = sum(W);
        P_cumulative          = ClimateForcingStep.P +  P_cumulative;
        EC_cumulative         = EC + EC_cumulative;   % evap(-) / cond(+)
        Ra_cumulative         = Ra + Ra_cumulative;
        F_cumulative          = F + F_cumulative;
      
        % calculate total system mass
        M_total    = sum(dz .* d);
        M_change   = M_total + R_cumulative + W_total- P_cumulative - EC_cumulative - M_initial - M_added_cumulative;
        M_change   = round(M_change * 100)/100;

        % check mass conservation
        if M_change ~= 0
            error('total system mass not conserved in MB function')
        end

        % check bottom grid cell T is unchanged
        if abs(T(end)-T_bottom) > 0.001
            error('temperature of bottom grid cell changed: original = %0.10g J, updated = %0.10g J',T_bottom,T(end))
        end

        % !! This needs to be made into a function call !!!
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
                 datestr(daten(date_ind), 'dd-mmm-yyyy'), ...
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