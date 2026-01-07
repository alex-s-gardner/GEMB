% Specify run case
run_id = "test_1";


%% RUN GEMB
switch run_id
    case "test_1"
        verbose = true;

        % [1] specify model parameters
        % [if data is not modified then it can be passed as an stucture]
        ModelParam = ...
            model_initialize_parameters();

        % [2] Load in climate data 
        %     CF = Climate Forcing
        %     LP = Location Specifc Parameters
        % [if data is not modified then it can be passed as an stucture]
        ClimateForcing = ...
            simulate_climate_forcing(run_id);

        % [3] Initialize grid -or- load in data to restart a simulation
        % [these variables are modified within the model and should not be
        % bundled in a stucture]
        [T, dz, d, W, re, gdn, gsp, a, a_diffuse] = ...
            model_initialize_column(ModelParam, ClimateForcing);
        
        % [4] Rum GEMB
        OutData = ...
            gemb(T, dz, d, W, re, gdn, gsp, a, a_diffuse, ClimateForcing, ModelParam, verbose);
        
        % [5] Save model output and model run settings

        %save(fullfile(S.output_dir, S.run_id), '-struct', 'O', '-v7.3')
        %save(fullfile(S.output_dir, S.run_id), 'S', '-append')

    otherwise
        error("input case not defined")
end