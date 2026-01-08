function [output_index, OutData, OutCum] = model_initialize_output(column_length, ClimateForcing, ModelParam)

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