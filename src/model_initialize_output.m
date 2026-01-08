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
%  column_length  : integer      Number of vertical grid cells in the model column.
%  ClimateForcing : struct       Forcing data structure containing:
%    .daten       : datenum      Vector of time steps.
%    .T_air0      : K            Air temperature.
%    .P0          : m w.e.       Precipitation.
%  ModelParam     : struct       Model parameters structure containing:
%    .output_frequency : string  Frequency of output ('daily', 'monthly', 'all').
%    .output_padding   : integer Extra buffer size for profile arrays.
%
%% Outputs
%
%  output_index   : logical      Boolean mask indicating which time steps are saved.
%  OutData        : struct       Structure for storing time-series outputs, initialized with NaNs.
%                                Contains fields for monolevel (e.g., 'M', 'R') and profile (e.g., 'T', 'd') variables.
%  OutCum         : struct       Structure for tracking cumulative values between outputs.
%                                Initialized to zero for fields like 'M_added', 'sw_net', 'Q_net'.
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

    % deteremine save time steps
    date_vector = datevec([ClimateForcing.daten; (ClimateForcing.daten(end) ...
        + ClimateForcing.daten(end)-ClimateForcing.daten(end-1))]);
    switch ModelParam.output_frequency
        case "monthly"
            output_index = (date_vector(1:end-1,2) - date_vector(2:end,2)) ~= 0;
        case "daily"
            output_index = (date_vector(1:end-1,3) - date_vector(2:end,3)) ~= 0;
        case "all"
            output_index = true(length(date_vector)-1);
    end

    % single level time series
    varname.monolevel = {'time', 'M', 'R', 'F', 'EC', 'sw_net', ...
        'lw_net', 'shf', 'lhf', 'a1', 'Q_net', 're1', 'd1', 'm', ...
        'compaction_dens', 'compaction_melt', 'ps'};
    
    n = sum(output_index);
    
    % set single level time series to NaNs
    OutData.time = ClimateForcing.daten(output_index)';

    for v = 1:length(varname.monolevel)
        OutData.(varname.monolevel{v}) = nan(1,n);
    end

    % time averages/totals
    I = find(output_index);                      % save index
    for i = 1:n
        if i == 1
            OutData.T_air(i) = mean(ClimateForcing.T_air0(1:I(i))) - 273.15;  % convert from K to deg C
            OutData.P(i)  = sum(ClimateForcing.P0(1:I(i)));
        else
            OutData.T_air(i) = mean(ClimateForcing.T_air0((I(i-1)+1):I(i))) - 273.15;
            OutData.P(i)  = sum(ClimateForcing.P0((I(i-1)+1):I(i)));
        end
    end

    % multi level time series
    varname.profile = {'T', 'dz', 'd', 'W', 're', 'gdn', 'gsp', 'a', 'a_diffuse', 'ps'};

    for v = 1:length(varname.profile)
        OutData.(varname.profile{v}) = nan(column_length + ModelParam.output_padding, n);
    end
    
    % initialize cumulative output values
    varname.cumulative = {'R', 'M', 'F', 'EC', 'Ra', 'M_added', 'sw_net', ...
        'lw_net', 'shf', 'lhf', 'a1', 're1', 'ulw', 'd1', 'compaction_dens', ...
        'compaction_melt', 'Q_net', 'FAC'};
    
    % set cumulative values zero
    for v = 1:length(varname.cumulative)
        OutCum.(varname.cumulative{v}) = 0;
    end
    
    OutCum.count = 0;
end