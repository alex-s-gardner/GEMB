function O = GEMB(T, dz, d, W, re, gdn, gsp, a, a_diffuse, daten, T_air0, V0, dlw0, dsw0, e_air0, p_air0, P0, S, LP, verbose)
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


disp(['------------------ STARTING RUN # ' num2str(S.run_prefix) ' --------------------' ])
tic                                     % start timer
dt = (daten(2)-daten(1)) * (60*60*24);  % input time step in seconds

if rem(dt,1) ~= 0
    warning('rounding dt as it is not an exact integer: dt = %0.4f', dt)
    dt = round(dt);
end

% initialize monolevel variables
EC     = 0;                        % surface evaporation (-) condensation (+) [kg m-2]    
M_surf = 0;                        % initialize surface melt for albedo parameterization

% fixed lower temperature bounday condition - T is fixed
T_bottom = T(end);

% deteremine save time steps
date_vector = datevec([daten; (daten(end) + daten(end)-daten(end-1))]);
switch S.output_frequency
    case "monthly"
        output_index = (date_vector(1:end-1,2) - date_vector(2:end,2)) ~= 0;
    case "daily"
        output_index = (date_vector(1:end-1,3) - date_vector(2:end,3)) ~= 0;
    case "3hourly"
        output_index = (date_vector(1:end-1,4) - date_vector(2:end,4)) ~= 0;
end

%% initialize output structure
column_length = length(dz);
[O, OC] = model_output_initialize(daten, T_air0, P0, column_length, output_index);

%% Start year loop for model spin up
for yIdx = 1:S.n_spinup_cycles + 1

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
    for dIdx = 1:length(daten)

        % Extract daily data:
        [T_air, V, dlw, dsw, e_air, p_air, P, black_carbon_snow, ...
            black_carbon_ice, cloud_optical_thickness, ...
            solar_zenith_angle, dsw_diffuse, cloud_fraction] = ...
            model_inputs_single_timestep(dIdx, T_air0, V0, dlw0, dsw0, ...
            e_air0, p_air0, P0, S);
        
        % run GEMB for a single time interval
        [T, dz, d, W, re, gdn, gsp, a, a_diffuse, EC, M_surf, sw_net, shf, ...
            lhf, ulw, Ra, M, R, F, M_added, E_added, ...
            compaction_dens, compaction_melt] = ...
           gemb_core(T, dz, d, W, re, gdn, gsp, a, a_diffuse, dt, P, EC, M_surf, ...
            black_carbon_snow, black_carbon_ice, solar_zenith_angle, ...
            cloud_optical_thickness, cloud_fraction, dsw, dsw_diffuse, dlw, T_air, ...
            V, e_air, p_air, S, LP, verbose);

        % calculate net longwave [W m-2]
        lw_net = dlw - ulw;

        % sum component mass changes [kg m-2]
        M_added_cumulative    = M_added + M_added_cumulative;
        M_cumulative          = M + M_cumulative;
        M_surf_cumulative     = M_surf + M_surf_cumulative;
        R_cumulative          = R + R_cumulative;
        W_total               = sum(W);
        P_cumulative          = P +  P_cumulative;
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
        if yIdx == S.n_spinup_cycles + 1
            % initialize cumulative and average variables for output
            d1    = d(1);
            a1    = a(1);
            re1   = re(1);
            Q_net = sw_net + lw_net + shf + lhf;
            FAC                   = sum(dz.*(S.density_ice - min(d,S.density_ice)))/1000;
            ov_varname = fieldnames(OC);
            for v = 1:length(ov_varname)
                if (ov_varname{v} == "count") 
                    continue
                else
                    OC.(ov_varname{v}) = eval(ov_varname{v}) + OC.(ov_varname{v});
                end
            end

            OC.count = OC.count + 1;

            if output_index(dIdx)
                % Store model output in structure format

                % time averaged monolevel values
                r = sum(output_index(1:dIdx));

                for v = 1:length(ov_varname)
                    % check if output is a cumulative value
                    if sum(strcmp(ov_varname{v}, {'M', 'R', 'F', 'EC', 'P', 'Ra','M_added', 'compaction_dens', 'compaction_melt'})) == 1
                        O.(ov_varname{v})(r) = OC.(ov_varname{v});
                    else
                        % if not cumulative then divide by time steps
                        O.(ov_varname{v})(r) = OC.(ov_varname{v}) / OC.count;
                    end
                end

                % instantaneous level data
                o = (size(d,1) - 1);
                O.re(end-o:end,r)  = re;
                O.d(end-o:end,r)   = d;
                O.T(end-o:end,r)   = T;
                O.W(end-o:end,r)   = W;
                O.dz(end-o:end,r)  = dz;
                O.gdn(end-o:end,r) = gdn;
                O.gsp(end-o:end,r) = gsp;
                O.ps(end-o:end,r)  = sum(dz) - M_total/S.density_ice;
                O.m(r)             = o+1;

                % set cumulative values back to zero
                for v = 1:length(ov_varname)
                    OC.(ov_varname{v}) = 0;
                end

                OC.count = 0;
            end
        end
    end

    % display cycle completed and time to screen
    disp([num2str(S.run_id) ': cycle: ' num2str(yIdx) ' of '  ...
        num2str(S.n_spinup_cycles + 1) ', cpu time: ' num2str(round(toc)) ' sec,'...
        ' avg melt: ' num2str(round(M_cumulative/(daten(end)-daten(1))*365.25)) ...
        ' kg/m2/yr']);
end
end
