% Run synthetic example and save outputs for Julia comparison
cd '/Users/gardnera/Documents/GitHub/GEMB'
addpath('src')

% Generate 3-hourly synthetic climate forcing data:
time_step_hours = 3;
ClimateForcing = simulate_climate_forcing("test_1", time_step_hours);

% Initialize model parameters:
ModelParam = model_initialize_parameters('output_frequency', "daily");

% Initialize a column:
Profile = model_initialize_profile(ModelParam, ClimateForcing);

% Create a climatological average time series:
ClimateForcingClimatology = forcing_climatology(ClimateForcing);

% Spinup a profile for 75 years of average forcing:
Profile_spunup = gemb_spinup(Profile, ClimateForcingClimatology, ModelParam, 75);

% Run model
OutData = gemb(Profile_spunup, ClimateForcing, ModelParam);

% Print summary
fprintf('MATLAB Simulation complete!\n')
fprintf('  Time steps: %d\n', length(OutData.melt))
fprintf('  Profile layers: %d\n', size(OutData.temperature, 1))
fprintf('  Mean surface albedo: %.3f\n', mean(OutData.albedo_surface))
fprintf('  Total melt: %.2f kg/m²\n', sum(OutData.melt))
fprintf('  Total runoff: %.2f kg/m²\n', sum(OutData.runoff))

% Save key outputs
save('/Users/gardnera/Documents/GitHub/GEMB.jl/examples/matlab_synthetic_output.mat', 'OutData', '-v7.3');
fprintf('Saved outputs to matlab_synthetic_output.mat\n')
