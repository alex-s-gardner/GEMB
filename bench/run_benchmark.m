% MATLAB benchmark for the GEMB hot path.
% Mirrors bench/opt_bench.jl: 1994-2025 3-hourly simulation with 75-year spinup.
function run_benchmark()
    here = fileparts(mfilename('fullpath'));
    repo = fileparts(here);
    addpath(fullfile(repo, 'src'));

    time_step_hours = 3;

    build = @() build_inputs(time_step_hours);
    [cf, mp, profile, cf_clim] = build();

    fprintf('Warmup...\n');
    run_hotpath(profile, cf_clim, mp, cf);

    nsamples = 5;
    times = zeros(nsamples, 1);
    for i = 1:nsamples
        [cf, mp, profile, cf_clim] = build();
        t0 = tic;
        run_hotpath(profile, cf_clim, mp, cf);
        times(i) = toc(t0);
        fprintf('  sample %d: %.3f s\n', i, times(i));
    end

    fprintf('\nMATLAB hot path (1994-2025 3-hourly, 75-year spinup):\n');
    fprintf('  min    = %.3f s\n', min(times));
    fprintf('  median = %.3f s\n', median(times));
    fprintf('  mean   = %.3f s\n', mean(times));
end

function [cf, mp, profile, cf_clim] = build_inputs(time_step_hours)
    cf = simulate_climate_forcing("test_1", time_step_hours);
    mp = model_initialize_parameters(output_frequency="daily");
    profile = model_initialize_profile(mp, cf);
    cf_clim = forcing_climatology(cf);
end

function out = run_hotpath(profile, cf_clim, mp, cf)
    profile_spunup = gemb_spinup(profile, cf_clim, mp, 75);
    out = gemb(profile_spunup, cf, mp);
end
