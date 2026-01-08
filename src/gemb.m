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
%
%
%% Description
%
%
%
%% Inputs
%
%
%
%% Outputs
%
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

% fixed lower temperature bounday condition - T is fixed
T_bottom = T(end);

%% initialize output structure
column_length = length(dz);
[output_index, OutData, OutCum] = model_initialize_output(column_length, ClimateForcing, ModelParam);

%% Start spinup loop
for simulation_iteration = 1:(ModelParam.n_spinup_cycles + 1)

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
    end

    % display cycle completed and time to screen
    disp([num2str(ModelParam.run_prefix) ': cycle: ' num2str(simulation_iteration) ' of '  ...
        num2str(ModelParam.n_spinup_cycles + 1) ', cpu time: ' num2str(round(toc)) ' sec,'...
        ' avg melt: ' num2str(round(M_cumulative/(daten(end)-daten(1))*365.25)) ...
        ' kg/m2/yr']);
end
end
