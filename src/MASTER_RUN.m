% Source of climate forcing
run_id    = "test_1";


%% RUN GEMB
switch run_id
    case "test_1"
        verbose = true;

        % [1] specify model parameters
        S = model_initialize_parameters()

        % [2] Load in climate data
        [daten, P, T_air, V, dlw, dsw, e_air, p_air, LP] = simulate_climate_forcing(run_id);

        % [3] Initialize grid -or- load in data to restart a simulation
        [T, dz, d, W, re, gdn, gsp, a, a_diffuse] = model_initialize_column(S, LP);
        
        % [4] Rum GEMB
        O = GEMB(T, dz, d, W, re, gdn, gsp, a, a_diffuse, daten, T_air, V, dlw, dsw, e_air, p_air, P, S, LP, verbose);
    
        % [5] Save model output and model run settings

        %save(fullfile(S.output_dir, S.run_id), '-struct', 'O', '-v7.3')
        %save(fullfile(S.output_dir, S.run_id), 'S', '-append')

    otherwise
        error("input case not defined")
end