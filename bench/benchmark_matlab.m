% MATLAB Performance Benchmark
% Equivalent to Julia's benchmark_vs_matlab.jl

fprintf('======================================================================\n');
fprintf('GEMB MATLAB Performance Benchmark\n');
fprintf('======================================================================\n');

addpath("../src")

% Warmup run
fprintf('\n[1/3] Warmup run...\n');
tic;
time_step_hours_warmup = 3;
ClimateForcing_warmup = simulate_climate_forcing("test_1", time_step_hours_warmup);
ModelParam_warmup = model_initialize_parameters(output_frequency="daily");
Profile_warmup = model_initialize_profile(ModelParam_warmup, ClimateForcing_warmup);
ClimateForcingClimatology_warmup = forcing_climatology(ClimateForcing_warmup);
ModelParam_spinup_warmup = model_initialize_parameters(output_frequency="last");
Profile_spunup_warmup = gemb_spinup(Profile_warmup, ClimateForcingClimatology_warmup, ModelParam_spinup_warmup, 75);
OutData_warmup = gemb(Profile_spunup_warmup, ClimateForcing_warmup, ModelParam_warmup);
warmup_time = toc;
fprintf('   Warmup complete: %.2f seconds\n', warmup_time);

% Benchmark run #1
fprintf('\n[2/3] Benchmark run #1...\n');
tic;
time_step_hours = 3;
ClimateForcing1 = simulate_climate_forcing("test_1", time_step_hours);
ModelParam1 = model_initialize_parameters(output_frequency="daily");
Profile1 = model_initialize_profile(ModelParam1, ClimateForcing1);
ClimateForcingClimatology1 = forcing_climatology(ClimateForcing1);
ModelParam_spinup1 = model_initialize_parameters(output_frequency="last");
Profile_spunup1 = gemb_spinup(Profile1, ClimateForcingClimatology1, ModelParam_spinup1, 75);
OutData1 = gemb(Profile_spunup1, ClimateForcing1, ModelParam1);
run1_time = toc;
fprintf('   Run #1 complete: %.2f seconds\n', run1_time);

% Benchmark run #2
fprintf('\n[3/3] Benchmark run #2...\n');
tic;
ClimateForcing2 = simulate_climate_forcing("test_1", time_step_hours);
ModelParam2 = model_initialize_parameters(output_frequency="daily");
Profile2 = model_initialize_profile(ModelParam2, ClimateForcing2);
ClimateForcingClimatology2 = forcing_climatology(ClimateForcing2);
ModelParam_spinup2 = model_initialize_parameters(output_frequency="last");
Profile_spunup2 = gemb_spinup(Profile2, ClimateForcingClimatology2, ModelParam_spinup2, 75);
OutData2 = gemb(Profile_spunup2, ClimateForcing2, ModelParam2);
run2_time = toc;
fprintf('   Run #2 complete: %.2f seconds\n', run2_time);

% Verify consistency
melt1 = sum(OutData1.melt, 'omitnan');
melt2 = sum(OutData2.melt, 'omitnan');
fprintf('\n   Consistency check: %.6f == %.6f ', melt1, melt2);
if abs(melt1 - melt2) < 1e-10
    fprintf('✓ PASS\n');
else
    fprintf('✗ FAIL\n');
end

% Statistics
avg_time = (run1_time + run2_time) / 2;
fprintf('\n======================================================================\n');
fprintf('MATLAB PERFORMANCE SUMMARY\n');
fprintf('======================================================================\n');
fprintf('Warmup time: %.2f seconds\n', warmup_time);
fprintf('Run #1 time: %.2f seconds\n', run1_time);
fprintf('Run #2 time: %.2f seconds\n', run2_time);
fprintf('Average runtime (runs 1-2): %.2f seconds\n', avg_time);

% Output statistics
fprintf('\n======================================================================\n');
fprintf('OUTPUT STATISTICS (for cross-validation)\n');
fprintf('======================================================================\n');
n_timesteps = length(OutData1.time);
n_layers = size(OutData1.temperature, 1);
mean_albedo = mean(OutData1.albedo_surface, 'omitnan');
total_melt = sum(OutData1.melt, 'omitnan');
total_runoff = sum(OutData1.runoff, 'omitnan');
total_refreeze = sum(OutData1.refreeze, 'omitnan');
total_precipitation = sum(OutData1.precipitation, 'omitnan');
thickness_final = OutData1.thickness_cumulative(end);

fprintf('Time steps: %d\n', n_timesteps);
fprintf('Profile layers: %d\n', n_layers);
fprintf('Mean surface albedo: %.6f\n', mean_albedo);
fprintf('Total melt: %.6f kg/m²\n', total_melt);
fprintf('Total runoff: %.6f kg/m²\n', total_runoff);
fprintf('Total refreezing: %.6f kg/m²\n', total_refreeze);
fprintf('Total precipitation: %.6f kg/m²\n', total_precipitation);
fprintf('Final cumulative thickness: %.6f m\n', thickness_final);

% Timing breakdown
fprintf('\n======================================================================\n');
fprintf('TIMING BREAKDOWN (Run #2)\n');
fprintf('======================================================================\n');
fprintf('\nRunning detailed timing analysis...\n');

tic;
ClimateForcing_timing = simulate_climate_forcing("test_1", time_step_hours);
t_forcing = toc;
fprintf('  Climate forcing generation: %.3f s\n', t_forcing);

tic;
ModelParam_timing = model_initialize_parameters(output_frequency="daily");
Profile_timing = model_initialize_profile(ModelParam_timing, ClimateForcing_timing);
t_init = toc;
fprintf('  Profile initialization: %.3f s\n', t_init);

tic;
ClimateForcingClimatology_timing = forcing_climatology(ClimateForcing_timing);
t_clim = toc;
fprintf('  Climatology creation: %.3f s\n', t_clim);

tic;
ModelParam_spinup_timing = model_initialize_parameters(output_frequency="last");
Profile_spunup_timing = gemb_spinup(Profile_timing, ClimateForcingClimatology_timing, ModelParam_spinup_timing, 75);
t_spinup = toc;
fprintf('  Spinup (75 years): %.3f s\n', t_spinup);

tic;
OutData_timing = gemb(Profile_spunup_timing, ClimateForcing_timing, ModelParam_timing);
t_main = toc;
fprintf('  Main simulation: %.3f s\n', t_main);

total_timing = t_forcing + t_init + t_clim + t_spinup + t_main;
fprintf('\n  Total (sum of components): %.3f s\n', total_timing);

fprintf('\n======================================================================\n');
fprintf('Benchmark complete!\n');
fprintf('======================================================================\n');

% Save results
save('matlab_benchmark_results.mat', 'run1_time', 'run2_time', 'avg_time', ...
     'total_melt', 'total_runoff', 'mean_albedo', 'OutData1', '-v7.3');
fprintf('\nResults saved to matlab_benchmark_results.mat\n');
