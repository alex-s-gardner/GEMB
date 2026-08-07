% Total-workflow benchmark: 75-yr climatological spinup + 32-yr transient run
% (~107 model-years, 3-hourly forcing). Times the physics (spinup + gemb),
% excluding forcing synthesis, post-warmup. Matches the Julia measurement.

addpath("../src")

time_step_hours = 3;

% Build inputs once (not timed)
cf   = simulate_climate_forcing("test_1", time_step_hours);
mp   = model_initialize_parameters(output_frequency="daily");
prof = model_initialize_profile(mp, cf);
cfc  = forcing_climatology(cf);
mps  = model_initialize_parameters(output_frequency="last");

% Warmup (JIT / caches)
ps = gemb_spinup(prof, cfc, mps, 75);
out = gemb(ps, cf, mp);

% Timed: full physics workflow, minimum of a few runs
nrun = 3;
tot = inf;
for k = 1:nrun
    tic;
    ps  = gemb_spinup(prof, cfc, mps, 75);   % 75-yr spinup
    out = gemb(ps, cf, mp);                   % 32-yr transient run
    tk = toc;
    fprintf('  run %d: %.3f s\n', k, tk);
    tot = min(tot, tk);
end

fprintf('\nMATLAB total workflow (~107 model-yr, 3-hourly): min = %.3f s over %d runs\n', tot, nrun);
