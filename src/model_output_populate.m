function [OutData, OutCum] = ...
    model_output_populate(d, T, W, dz, re, gdn, gsp, ...
     output_index, date_ind, ModelParam, OutData, OutCum)
% model_output_populate stores model state and fluxes into the output structure.
%
%% Syntax
%
% [OutData, OutCum] = model_output_populate(OutData, OutCum, ...
%    d, T, W, dz, re, gdn, gsp, ...
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
%  d, T, W, ...     : vectors      Current vertical profiles of density, temperature, etc.
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
    OutData.re(end-o:end,r)  = re;
    OutData.d(end-o:end,r)   = d;
    OutData.T(end-o:end,r)   = T;
    OutData.W(end-o:end,r)   = W;
    OutData.dz(end-o:end,r)  = dz;
    OutData.gdn(end-o:end,r) = gdn;
    OutData.gsp(end-o:end,r) = gsp;
    
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