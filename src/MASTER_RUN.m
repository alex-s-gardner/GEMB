% Source of climate forcing
run_id    = "test_1";
S = model_setup()

%% RUN GEMB
switch run_id
    case "test_1"
        verbose = true;
        % output directory
        S.output_dir = '../test_1';

        [daten, P, T_air, V, dlw, dsw, e_air, p_air, LP] = simulate_climate_forcing(run_id);
        S = combineStrucData_GEMB(S,LP,1);
        S.V_mean = mean(V);
    
        GEMB(daten, T_air, V, dlw, dsw, e_air, p_air, P, S, S.is_restart, verbose)
    otherwise
        error("input case not defined")
end